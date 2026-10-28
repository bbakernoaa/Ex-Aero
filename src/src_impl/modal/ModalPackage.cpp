/// @file ModalPackage.cpp
/// @brief YAML/config validation, table reuse, and host/device routing.
/// @author EX-aero contributors
/// @details History: 2026-10. Input: YAML/in-memory host arrays; no file
/// output.
#include <algorithm>
#include <CloudIceYaml.hpp>
#include <exaero/ModalPresets.hpp>
#include <exaero/PackageFactory.hpp>
#ifdef EXAERO_ENABLE_MAM4XX
#include <mam4xx/Mam4xxSolver.hpp>
#endif
#include <cmath>
#include <limits>
#include <loader/MieFileLoader.hpp>
#include <loader/MieTableStore.hpp>
#include <map>
#include <modal/ModalPackage.hpp>
#include <numeric>
#include <set>
#include <yaml-cpp/yaml.h>

namespace exaero {
namespace {
constexpr double pi = 3.14159265358979323846;
/// @brief Fail at the host boundary before changing any caller-owned memory.
void require(bool ok, const std::string &message) {
  if (!ok)
    throw std::invalid_argument("FATAL ERROR: modal: " + message);
}
/// @brief Positive finite scalar validation.
bool positive(double x) { return std::isfinite(x) && x > 0; }
/// @brief Reject typos and structurally invalid YAML instead of ignoring them.
void keys(const YAML::Node &n, std::initializer_list<const char *> names) {
  require(n.IsMap(), "expected a YAML mapping");
  for (const auto &e : n) {
    const auto key = e.first.as<std::string>();
    require(std::any_of(names.begin(), names.end(),
                        [&](const char *s) { return key == s; }),
            "unknown setting: " + key);
  }
}
/// @brief Extract a private, contiguous host-pointer environmental descriptor.
modal::Environment raw(const EnvironmentalStateView &e, bool full = true) {
  return {e.temperature.data_handle(),
          full ? e.pressure.data_handle() : nullptr,
          full ? e.air_density.data_handle() : nullptr,
          e.relative_humidity.data_handle(),
          full ? e.layer_thickness.data_handle() : nullptr};
}
/// @brief Validate four-dimensional output, including checked extent products.
void shape(const View4D<double> &v, std::size_t nc, std::size_t nl,
           std::size_t n2, std::size_t n3) {
  require(v.data_handle() && v.extent(0) == nc && v.extent(1) == nl &&
              v.extent(2) == n2 && v.extent(3) == n3,
          "output extent mismatch or null output");
  const auto limit = static_cast<std::size_t>(std::numeric_limits<int>::max());
  require(n2 > 0 && n3 > 0 && nc <= limit / nl / n2 / n3,
          "output slab exceeds execution-policy limit");
}
/// @brief Validate an increasing interpolation coordinate, allowing one node.
void axis(const std::vector<double> &v, bool strictly_positive) {
  require(!v.empty(), "optical table axis is empty");
  for (std::size_t i = 0; i < v.size(); ++i)
    require(std::isfinite(v[i]) && (strictly_positive ? v[i] > 0 : v[i] >= 0) &&
                (i == 0 || v[i] > v[i - 1]),
            "optical axes must be finite and strictly increasing");
}
/// @brief Copy a resolved source curve; its optical semantics remain explicit.
modal::Table table(const ModalTableBinding &binding,
                   const ModalSpeciesConfig &material) {
  const auto &name = material.name;
  require(binding.curve.interpolate == "linear",
          "unsupported optical interpolation strategy");
  // Use existing loader/merge/override behavior, then take an owned snapshot so
  // later package initializations cannot alter this package's tables.
  MieTableStore store;
  store.reset_to_baked_in();
  store.ensure_modal_loaded();
  if (!binding.data_file.empty())
    store.apply_runtime_file(binding.data_file);
  auto cc = binding.curve;
  cc.species_name = "modal_" + name;
  require(!cc.source_label.empty(), "optical source label is required");
  store.apply_curve_config({cc});
  const auto *curve = store.find_curve(cc.species_name);
  require(curve != nullptr, "optical source could not be resolved");
  const auto kind =
      binding.kind == ModalTableKind::Unknown && curve->modal_particle_resolved
          ? ModalTableKind::ParticleResolved
          : binding.kind;
  require(kind != ModalTableKind::Unknown,
          "optical table distribution semantics are unknown for " + name);
  if (curve->modal_particle_resolved) {
    require(kind == ModalTableKind::ParticleResolved,
            "generated single-particle curves cannot be declared distribution "
            "averages");
    require(std::abs(material.density / curve->modal_density - 1) < 1e-10 &&
                std::abs(material.hygroscopicity - curve->modal_kappa) < 1e-10,
            "material density/kappa must match the generated optical recipe");
  }
  modal::Table out;
  out.source_label = curve->source_label;
  out.source_version = curve->source_version;
  out.citation = curve->citation;
  out.kind = kind == ModalTableKind::ParticleResolved ? 1 : 2;
  out.sigma = binding.geometric_std_dev;
  out.radius = curve->radius;
  out.rh = curve->rh;
  const auto wavelengths =
      binding.wavelength_nodes.empty() && curve->physical_wavelengths
          ? curve->lambda
          : binding.wavelength_nodes;
  require(wavelengths.size() == curve->lambda.size(),
          "physical wavelength_nodes must map every source spectral node "
          "(RRTMG axes are band IDs)");
  out.wavelength = wavelengths;
  axis(out.radius, true);
  require(out.kind != 1 || out.radius.size() >= 2,
          "particle-resolved curves need at least two radius nodes");
  axis(out.rh, false);
  axis(out.wavelength, true);
  require(out.rh.back() <= 1, "table RH exceeds one");
  const auto count = out.radius.size() * out.rh.size() * out.wavelength.size();
  const auto field = [&](const std::string &key, const char *unit) {
    const auto it = curve->fields.find(key);
    require(it != curve->fields.end() && it->second.rank == 3 &&
                it->second.values.size() == count && it->second.unit == unit,
            "missing or incompatible optical field: " + key);
    return it->second.values;
  };
  out.extinction = field("bext", "m^2 (kg dry mass)^-1");
  out.scattering = field("bsca", "m^2 (kg dry mass)^-1");
  out.asymmetry = field("g", "1");
  out.backscatter = field("bbck", "m^2 (kg dry mass)^-1 sr^-1");
  for (std::size_t i = 0; i < count; ++i)
    require(std::isfinite(out.extinction[i]) && out.extinction[i] >= 0 &&
                std::isfinite(out.scattering[i]) && out.scattering[i] >= 0 &&
                out.scattering[i] <= out.extinction[i] * (1 + 1e-6) &&
                std::isfinite(out.asymmetry[i]) &&
                std::abs(out.asymmetry[i]) <= 1 &&
                std::isfinite(out.backscatter[i]) && out.backscatter[i] >= 0,
            "invalid optical coefficient");
  return out;
}
/// @brief Check and add one complete optical-provider result.
void addOptics(const ModalOpticalResult &a, View4D<double> &out, std::size_t c,
               std::size_t k, std::size_t b) {
  require(std::isfinite(a.asymmetry) && std::abs(a.asymmetry) <= 1,
          "provider returned invalid asymmetry");
  for (double x : {a.extinction, a.scattering, a.backscatter,
                   a.lidar_backscatter, a.fine_extinction, a.fine_scattering,
                   a.pm25_extinction, a.pm25_scattering})
    require(std::isfinite(x) && x >= 0,
            "provider returned invalid optical coefficient");
  require(a.scattering <= a.extinction * (1 + 1e-6) &&
              a.fine_extinction <= a.extinction * (1 + 1e-6) &&
              a.pm25_extinction <= a.extinction * (1 + 1e-6) &&
              a.fine_scattering <= a.scattering * (1 + 1e-6) &&
              a.pm25_scattering <= a.scattering * (1 + 1e-6),
          "provider violated optical bounds");
  using namespace optical_indices;
  out(c, k, b, EXTINCTION_COEFF) += a.extinction;
  out(c, k, b, SCATTERING_COEFF) += a.scattering;
  out(c, k, b, ASYMMETRY_FACTOR) += a.asymmetry * a.scattering;
  out(c, k, b, BACKSCATTER_COEFF) += a.backscatter;
  out(c, k, b, LIDAR_BACKSCATTER) += a.lidar_backscatter;
}
} // namespace

void ModalPackage::initialize(const std::string &yaml) {
  try {
    const auto root = YAML::Load(yaml);
    keys(root,
         {"backend", "preset", "processes", "constants", "species", "modes",
          "mode_overrides", "aging_routes", "transfer_owner", "coagulation",
          "coagulation_routes", "transfer_routes", "nucleation", "resuspension_routes"});
    require(!root["backend"] || root["backend"].as<std::string>() == "modal",
            "backend must be modal");
    ModalConfig config;
    if (root["preset"]) {
      const auto preset = root["preset"].as<std::string>();
      require(preset == "mam4" || preset == "mam4xx_reference",
              "unknown modal preset");
      config =
          preset == "mam4" ? mam4ModalConfig() : mam4xxReferenceModalConfig();
    }
    require((root["species"] && root["species"].IsSequence()) ||
                (!root["species"] && !config.species.empty()),
            "species must be a sequence or supplied by a preset");
    require((root["modes"] && root["modes"].IsSequence()) ||
                (!root["modes"] && !config.modes.empty()),
            "modes must be a sequence or supplied by a preset");
    if (root["species"])
      config.species.clear();
    if (root["modes"]) {
      config.modes.clear();
      // Routes name modes: a replacement mode list needs explicit routes.
      config.coagulation_routes.clear();
      config.transfer_routes.clear();
      config.resuspension_routes.clear();
    }
    if (const auto processes = root["processes"]) {
      keys(processes, {"coagulation", "water_uptake"});
      if (processes["coagulation"]) {
        require(!root["coagulation"],
                "coagulation and processes.coagulation conflict");
        const auto name = processes["coagulation"].as<std::string>();
        require(name == "none" || name == "brownian" || name == "mam4xx",
                "unknown coagulation process");
        config.coagulation = name != "none";
        config.coagulation_model = name == "mam4xx"
                                       ? ModalCoagulationModel::Mam4xx
                                       : ModalCoagulationModel::Brownian;
        if (name == "mam4xx")
          config.coagulation_routes.clear();
      }
      if (processes["water_uptake"]) {
        const auto name = processes["water_uptake"].as<std::string>();
        require(name == "kappa_equilibrium" || name == "mam4xx",
                "unknown water uptake process");
        config.water_uptake = name == "mam4xx"
                                  ? ModalWaterUptakeModel::Mam4xx
                                  : ModalWaterUptakeModel::KappaEquilibrium;
      }
    }
    if (const auto constants = root["constants"]) {
      keys(constants, {"gravity", "water_molecular_weight", "gas_constant",
                       "surface_tension"});
      if (constants["gravity"])
        config.constants.gravity = constants["gravity"].as<double>();
      if (constants["water_molecular_weight"])
        config.constants.water_molecular_weight =
            constants["water_molecular_weight"].as<double>();
      if (constants["gas_constant"])
        config.constants.gas_constant = constants["gas_constant"].as<double>();
      if (constants["surface_tension"])
        config.constants.surface_tension =
            constants["surface_tension"].as<double>();
    }
    if (root["coagulation"])
      config.coagulation = root["coagulation"].as<bool>();
    for (const auto &n : root["species"]) {
      keys(n,
           {"name", "density", "molecular_weight", "hygroscopicity", "optics",
            "ice_categories"});
      ModalSpeciesConfig s;
      s.name = n["name"].as<std::string>();
      s.density = n["density"].as<double>();
      s.molecular_weight = n["molecular_weight"].as<double>();
      if (n["hygroscopicity"])
        s.hygroscopicity = n["hygroscopicity"].as<double>();
      s.ice_categories = cloud_ice_yaml::categories(n["ice_categories"]);
      if (const auto o = n["optics"]) {
        keys(o, {"source", "kind", "geometric_std_dev", "data_file",
                 "radius_nodes", "wavelength_nodes", "overrides"});
        ModalTableBinding binding;
        binding.curve.source_label = o["source"].as<std::string>();
        const auto kind = o["kind"] ? o["kind"].as<std::string>() : "unknown";
        require(kind == "particle_resolved" || kind == "lognormal_average" ||
                    kind == "unknown",
                "unknown optical table kind");
        binding.kind = kind == "unknown" ? ModalTableKind::Unknown
                       : kind == "particle_resolved"
                           ? ModalTableKind::ParticleResolved
                           : ModalTableKind::LognormalAverage;
        if (o["geometric_std_dev"])
          binding.geometric_std_dev = o["geometric_std_dev"].as<double>();
        if (o["data_file"])
          binding.data_file = o["data_file"].as<std::string>();
        if (o["radius_nodes"])
          binding.curve.radius_nodes =
              o["radius_nodes"].as<std::vector<double>>();
        if (o["wavelength_nodes"])
          binding.wavelength_nodes =
              o["wavelength_nodes"].as<std::vector<double>>();
        if (o["overrides"]) {
          require(o["overrides"].IsSequence(), "overrides must be a sequence");
          for (const auto &ov : o["overrides"]) {
            keys(ov, {"category", "attribute_index", "values"});
            SpeciesCurveConfig::Override v;
            const auto cat = ov["category"].as<std::string>();
            require(cat == "microphysical" || cat == "spectral_optical",
                    "unknown override category");
            v.category = cat == "microphysical"
                             ? AttributeCategory::Microphysical
                             : AttributeCategory::SpectralOptical;
            v.attribute_index = ov["attribute_index"].as<int>();
            v.values = ov["values"].as<std::vector<double>>();
            binding.curve.overrides.push_back(std::move(v));
          }
        }
        s.optics = std::move(binding);
      }
      config.species.push_back(std::move(s));
    }
    for (const auto &n : root["modes"]) {
      keys(n, {"name", "geometric_radius", "geometric_std_dev", "species",
               "optical_treatment", "crystallization_rh", "deliquescence_rh",
               "dry_geometric_radius_min", "dry_geometric_radius_max",
               "ice_rules"});
      ModalModeConfig m;
      m.name = n["name"].as<std::string>();
      if (n["crystallization_rh"])
        m.crystallization_rh = n["crystallization_rh"].as<double>();
      if (n["deliquescence_rh"])
        m.deliquescence_rh = n["deliquescence_rh"].as<double>();
      m.geometric_radius = n["geometric_radius"].as<double>();
      m.geometric_std_dev = n["geometric_std_dev"].as<double>();
      if (n["dry_geometric_radius_min"])
        m.dry_geometric_radius_min = n["dry_geometric_radius_min"].as<double>();
      if (n["dry_geometric_radius_max"])
        m.dry_geometric_radius_max = n["dry_geometric_radius_max"].as<double>();
      if (n["species"])
        m.species = n["species"].as<std::vector<std::string>>();
      m.ice_rules = cloud_ice_yaml::rules(n["ice_rules"]);
      if (n["optical_treatment"]) {
        const auto treatment = n["optical_treatment"].as<std::string>();
        require(treatment == "component_tables" ||
                    treatment == "internal_mixture",
                "unknown optical treatment");
        if (treatment == "internal_mixture")
          m.optical_treatment = ModalOpticalTreatment::InternalMixture;
      }
      config.modes.push_back(std::move(m));
    }
    if (const auto overrides = root["mode_overrides"]) {
      require(overrides.IsMap(),
              "mode_overrides must be a name-to-settings mapping");
      std::set<std::string> overridden;
      for (const auto &entry : overrides) {
        const auto name = entry.first.as<std::string>();
        auto mode = std::find_if(config.modes.begin(), config.modes.end(),
                                 [&](const auto &m) { return m.name == name; });
        require(mode != config.modes.end() && overridden.insert(name).second,
                "unknown or duplicate mode override: " + name);
        const auto &n = entry.second;
        keys(n, {"geometric_radius", "geometric_std_dev", "crystallization_rh",
                 "deliquescence_rh", "dry_geometric_radius_min",
                 "dry_geometric_radius_max", "ice_rules"});
        if (n["geometric_radius"])
          mode->geometric_radius = n["geometric_radius"].as<double>();
        if (n["geometric_std_dev"])
          mode->geometric_std_dev = n["geometric_std_dev"].as<double>();
        if (n["crystallization_rh"])
          mode->crystallization_rh = n["crystallization_rh"].as<double>();
        if (n["deliquescence_rh"])
          mode->deliquescence_rh = n["deliquescence_rh"].as<double>();
        if (n["dry_geometric_radius_min"])
          mode->dry_geometric_radius_min = n["dry_geometric_radius_min"].as<double>();
        if (n["dry_geometric_radius_max"])
          mode->dry_geometric_radius_max = n["dry_geometric_radius_max"].as<double>();
        if (n["ice_rules"])
          mode->ice_rules = cloud_ice_yaml::rules(n["ice_rules"]);
      }
    }
    for (const char *key : {"coagulation_routes", "transfer_routes"})
      if (root[key])
        require(root[key].IsSequence(),
                std::string(key) + " must be a sequence");
    if (root["coagulation_routes"])
      config.coagulation_routes.clear();
    if (root["coagulation_routes"])
      for (const auto &n : root["coagulation_routes"]) {
        keys(n, {"first", "second", "destination"});
        config.coagulation_routes.push_back(
            {n["first"].as<std::string>(), n["second"].as<std::string>(),
             n["destination"].as<std::string>()});
      }
    if (root["transfer_owner"]) {
      const auto owner = root["transfer_owner"].as<std::string>();
      require(owner == "host" || owner == "microphysics",
              "unknown transfer owner");
      config.transfer_owner = owner == "host"
                                  ? ModalTransferOwner::Host
                                  : ModalTransferOwner::Microphysics;
    }
    if (root["aging_routes"]) {
      require(root["aging_routes"].IsSequence(),
              "aging_routes must be a sequence");
      config.aging_routes.clear();
      for (const auto &n : root["aging_routes"]) {
        keys(n,
             {"source", "destination", "coating_species", "coating_thickness"});
        config.aging_routes.push_back(
            {n["source"].as<std::string>(), n["destination"].as<std::string>(),
             n["coating_species"].as<std::vector<std::string>>(),
             n["coating_thickness"].as<double>()});
      }
    }
    if (root["transfer_routes"])
      config.transfer_routes.clear();
    if (root["transfer_routes"])
      for (const auto &n : root["transfer_routes"]) {
        keys(n, {"source", "destination", "radius_threshold"});
        config.transfer_routes.push_back({n["source"].as<std::string>(),
                                          n["destination"].as<std::string>(),
                                          n["radius_threshold"].as<double>()});
      }
    if (root["resuspension_routes"]) {
      require(root["resuspension_routes"].IsSequence(),
              "resuspension_routes must be a sequence");
      config.resuspension_routes.clear();
      for (const auto &n : root["resuspension_routes"]) {
        keys(n, {"species", "source", "destination"});
        config.resuspension_routes.push_back(
            {n["species"].as<std::string>(), n["source"].as<std::string>(),
             n["destination"].as<std::string>()});
      }
    }
    if (const auto n = root["nucleation"]) {
      keys(n, {"destination", "prefactor", "exponent", "seed_radius",
               "mass_fractions", "precursor_molecular_weight"});
      ModalNucleationConfig v;
      v.destination = n["destination"].as<std::string>();
      v.prefactor = n["prefactor"].as<double>();
      v.exponent = n["exponent"].as<double>();
      v.seed_radius = n["seed_radius"].as<double>();
      v.mass_fractions = n["mass_fractions"].as<std::vector<double>>();
      if (n["precursor_molecular_weight"])
        v.precursor_molecular_weight =
            n["precursor_molecular_weight"].as<double>();
      config.nucleation = std::move(v);
    }
    initializeModal(config);
  } catch (const YAML::Exception &e) {
    throw std::invalid_argument(std::string("FATAL ERROR: modal YAML: ") +
                                e.what());
  }
}

void ModalPackage::initializeModal(const ModalConfig &config) {
  require(!config.species.empty() && !config.modes.empty(),
          "at least one species and mode are required");
  const auto ns = config.species.size(), nm = config.modes.size();
  require(ns < 100000 && nm < 100000 &&
              nm <= static_cast<std::size_t>(std::numeric_limits<int>::max()) /
                        9 / (ns + 1) &&
              nm <= static_cast<std::size_t>(std::numeric_limits<int>::max()) /
                        nm,
          "configuration exceeds supported index range");
  modal::Parameters p;
  const auto &constants = config.constants;
  require(positive(constants.gravity) &&
              positive(constants.water_molecular_weight) &&
              positive(constants.gas_constant) &&
              positive(constants.surface_tension),
          "physical constants must be finite and positive");
  p.gravity = constants.gravity;
  p.water_molecular_weight = constants.water_molecular_weight;
  p.gas_constant = constants.gas_constant;
  p.surface_tension = constants.surface_tension;
  require(config.coagulation_model == ModalCoagulationModel::Brownian ||
              config.coagulation_model == ModalCoagulationModel::Mam4xx,
          "unknown coagulation model");
  require(config.water_uptake == ModalWaterUptakeModel::KappaEquilibrium ||
              config.water_uptake == ModalWaterUptakeModel::Mam4xx,
          "unknown water uptake model");
  const bool native_coagulation =
      config.coagulation &&
      config.coagulation_model == ModalCoagulationModel::Mam4xx;
  p.mam4xx_water = config.water_uptake == ModalWaterUptakeModel::Mam4xx;
  require((!native_coagulation && !p.mam4xx_water) ||
              packageAvailable("mam4xx"),
          "selected native process requires EXAERO_ENABLE_MAM4XX=ON");
  p.coagulation = config.coagulation && !native_coagulation;
  if (native_coagulation) {
    const auto reference = mam4xxReferenceModalConfig();
    require(ns == 7 && nm == 4 && p.mam4xx_water &&
                config.coagulation_routes.empty(),
            "native coagulation requires four MAM4 modes, seven materials, "
            "native water and fixed native routes");
    for (std::size_t s = 0; s < ns; ++s) {
      const auto &a = config.species[s], &b = reference.species[s];
      require(a.name == b.name && std::abs(a.density / b.density - 1) < 1e-12 &&
                  std::abs(a.hygroscopicity - b.hygroscopicity) < 1e-14,
              "native coagulation requires MAM4 material order/density/kappa");
    }
    for (std::size_t m = 0; m < nm; ++m) {
      const auto &a = config.modes[m], &b = reference.modes[m];
      require(
          a.name == b.name &&
              std::abs(a.geometric_radius / b.geometric_radius - 1) < 1e-12 &&
              std::abs(a.geometric_std_dev / b.geometric_std_dev - 1) < 1e-12 &&
              a.crystallization_rh == b.crystallization_rh &&
              a.deliquescence_rh == b.deliquescence_rh && a.species.empty(),
          "native coagulation requires native mode order/widths/water "
          "thresholds and coating slots");
    }
  }
  std::map<std::string, int> species, modes;
  for (std::size_t s = 0; s < ns; ++s) {
    const auto &v = config.species[s];
    require(!v.name.empty() && v.name != "num" &&
                species.emplace(v.name, s).second,
            "species names must be unique, nonempty, and not num");
    require(positive(v.density) && positive(v.molecular_weight) &&
                std::isfinite(v.hygroscopicity) && v.hygroscopicity >= 0,
            "invalid material properties");
    validateIceConfiguration(v.ice_categories, {});
    p.density.push_back(v.density);
    p.kappa.push_back(v.hygroscopicity);
  }
  p.allowed.resize(nm * ns);
  std::vector<StateField> fields;
  std::set<std::string> field_names;
  for (std::size_t m = 0; m < nm; ++m) {
    const auto &v = config.modes[m];
    require(v.dry_geometric_radius_min.has_value() ==
                v.dry_geometric_radius_max.has_value(),
            "modal dry-radius bounds must be paired");
    if (v.dry_geometric_radius_min)
      require(positive(*v.dry_geometric_radius_min) &&
                  positive(*v.dry_geometric_radius_max) &&
                  *v.dry_geometric_radius_min < *v.dry_geometric_radius_max,
              "invalid modal dry-radius bounds");
    require(!v.name.empty() && modes.emplace(v.name, m).second,
            "mode names must be unique and nonempty");
    require(positive(v.geometric_radius) &&
                std::isfinite(v.geometric_std_dev) && v.geometric_std_dev > 1 &&
                std::isfinite(std::exp(
                    18 * std::pow(std::log(v.geometric_std_dev), 2))) &&
                positive(std::pow(v.geometric_radius, 3)),
            "invalid modal radius or spread");
    require(v.optical_treatment == ModalOpticalTreatment::ComponentTables ||
                v.optical_treatment == ModalOpticalTreatment::InternalMixture,
            "invalid optical treatment");
    require(std::isfinite(v.crystallization_rh) &&
                std::isfinite(v.deliquescence_rh) &&
                v.crystallization_rh >= 0 &&
                v.crystallization_rh < v.deliquescence_rh &&
                v.deliquescence_rh <= 1,
            "water transition RH must satisfy 0 <= crystallization < "
            "deliquescence <= 1");
    p.crystallization_rh.push_back(v.crystallization_rh);
    p.deliquescence_rh.push_back(v.deliquescence_rh);
    p.radius.push_back(v.geometric_radius);
    p.sigma.push_back(v.geometric_std_dev);
    std::set<std::string> permitted;
    for (const auto &name : v.species)
      require(species.count(name) && permitted.insert(name).second,
              "unknown or duplicate permitted species");
    for (std::size_t s = 0; s < ns; ++s)
      p.allowed[m * ns + s] =
          v.species.empty() || permitted.count(config.species[s].name);
    std::vector<IceCategory> available_ice_categories;
    for (std::size_t s = 0; s < ns; ++s)
      if (p.allowed[m * ns + s])
        for (auto category : config.species[s].ice_categories)
          if (std::find(available_ice_categories.begin(),
                        available_ice_categories.end(),
                        category) == available_ice_categories.end())
            available_ice_categories.push_back(category);
    validateIceConfiguration(available_ice_categories, v.ice_rules);
    for (std::size_t s = 0; s <= ns; ++s) {
      const auto name =
          (s == ns ? "num" : config.species[s].name) + "_" + v.name + "_i";
      require(field_names.insert(name).second,
              "ambiguous packed field name: " + name);
      fields.push_back({name,
                        s == ns ? StateQuantity::Number : StateQuantity::Mass,
                        static_cast<int>(m), false});
    }
  }
  const auto mode = [&](const std::string &name) {
    require(modes.count(name), "unknown mode: " + name);
    return modes.at(name);
  };
  const auto compatible = [&](int from, int to) {
    for (std::size_t s = 0; s < ns; ++s)
      require(!p.allowed[from * ns + s] || p.allowed[to * ns + s],
              "destination cannot contain a source material");
  };
  std::set<std::pair<int, int>> resuspension_sources;
  std::vector<std::vector<int>> resuspension_next(ns,
                                                  std::vector<int>(nm, -1));
  for (const auto &route : config.resuspension_routes) {
    require(species.count(route.species), "unknown resuspension species");
    const int s = species.at(route.species), a = mode(route.source),
              b = mode(route.destination);
    require(a != b && resuspension_sources.insert({s, a}).second &&
                p.allowed[a * ns + s] && p.allowed[b * ns + s],
            "invalid, duplicate or disallowed resuspension route");
    resuspension_next[s][a] = b;
  }
  for (std::size_t s = 0; s < ns; ++s)
    for (std::size_t m = 0; m < nm; ++m) {
      std::set<int> path;
      for (int next = m; next >= 0; next = resuspension_next[s][next])
        require(path.insert(next).second,
                "resuspension routes contain a cycle");
    }
  require(config.transfer_owner == ModalTransferOwner::Microphysics ||
              config.transfer_owner == ModalTransferOwner::Host,
          "unknown transfer owner");
  std::set<int> aging_sources;
  std::vector<int> aging_next(nm, -1);
  for (const auto &route : config.aging_routes) {
    const int a = mode(route.source), b = mode(route.destination);
    require(a != b && aging_sources.insert(a).second &&
                positive(route.coating_thickness) &&
                positive(std::pow(route.coating_thickness, 3)) &&
                !route.coating_species.empty(),
            "invalid or duplicate aging route");
    compatible(a, b);
    aging_next[a] = b;
    std::set<std::string> coatings;
    for (const auto &name : route.coating_species)
      require(species.count(name) && p.allowed[a * ns + species.at(name)] &&
                  coatings.insert(name).second,
              "unknown, disallowed or duplicate coating material");
  }
  for (std::size_t m = 0; m < nm; ++m) {
    std::set<int> path;
    for (int next = m; next >= 0; next = aging_next[next])
      require(path.insert(next).second, "aging routes contain a cycle");
  }
  p.destinations.assign(nm * nm, -1);
  for (std::size_t m = 0; m < nm; ++m)
    p.destinations[m * nm + m] = m;
  for (const auto &route : config.coagulation_routes) {
    const int a = mode(route.first), b = mode(route.second),
              d = mode(route.destination);
    require(a != b && p.destinations[a * nm + b] == -1,
            "duplicate or self coagulation route");
    compatible(a, d);
    compatible(b, d);
    p.destinations[a * nm + b] = p.destinations[b * nm + a] = d;
  }
  if (p.coagulation)
    for (int dest : p.destinations)
      require(dest >= 0,
              "every unordered intermode pair needs a coagulation destination");
  std::vector<int> indegree(nm), outgoing(nm, -1), order;
  for (std::size_t r = 0; r < config.transfer_routes.size(); ++r) {
    const auto &v = config.transfer_routes[r];
    const int a = mode(v.source), b = mode(v.destination);
    require(a != b && outgoing[a] == -1 && positive(v.radius_threshold),
            "invalid or duplicate growth-transfer source");
    compatible(a, b);
    outgoing[a] = r;
    ++indegree[b];
  }
  for (std::size_t m = 0; m < nm; ++m)
    if (!indegree[m])
      order.push_back(m);
  for (std::size_t k = 0; k < order.size(); ++k)
    if (outgoing[order[k]] >= 0) {
      const auto &v = config.transfer_routes[outgoing[order[k]]];
      const int b = mode(v.destination);
      if (--indegree[b] == 0)
        order.push_back(b);
    }
  require(order.size() == nm, "growth-transfer routes contain a cycle");
  for (int m : order)
    if (outgoing[m] >= 0) {
      const auto &v = config.transfer_routes[outgoing[m]];
      p.transfer_source.push_back(m);
      p.transfer_destination.push_back(mode(v.destination));
      p.transfer_radius.push_back(v.radius_threshold);
    }
  if (config.nucleation) {
    const auto &n = *config.nucleation;
    const int m = mode(n.destination);
    require(std::isfinite(n.prefactor) && n.prefactor >= 0 &&
                positive(n.exponent) && positive(n.seed_radius) &&
                positive(n.precursor_molecular_weight) &&
                n.mass_fractions.size() == ns,
            "invalid nucleation configuration");
    double sum = 0;
    int populated = 0;
    for (std::size_t s = 0; s < ns; ++s) {
      require(std::isfinite(n.mass_fractions[s]) && n.mass_fractions[s] >= 0,
              "invalid seed mass fraction");
      sum += n.mass_fractions[s];
      if (n.mass_fractions[s] > 0) {
        ++populated;
        require(p.allowed[m * ns + s] &&
                    std::abs(config.species[s].molecular_weight /
                                 n.precursor_molecular_weight -
                             1) < 1e-10,
                "empirical H2SO4 seeds require a permitted H2SO4-equivalent "
                "species");
      }
    }
    require(std::abs(sum - 1) < 1e-12 && populated == 1,
            "empirical scheme supports pure dry H2SO4-equivalent seeds only; "
            "mixed chemistry belongs to MKPP");
  }
  std::vector<modal::Table> tables(nm * ns);
  for (std::size_t s = 0; s < ns; ++s)
    if (config.species[s].optics) {
      const auto &binding = *config.species[s].optics;
      const auto snapshot = table(binding, config.species[s]);
      for (std::size_t m = 0; m < nm; ++m)
        if (p.allowed[m * ns + s] &&
            config.modes[m].optical_treatment ==
                ModalOpticalTreatment::ComponentTables) {
          if (snapshot.kind == 2)
            require(std::isfinite(snapshot.sigma) &&
                        std::abs(snapshot.sigma /
                                     config.modes[m].geometric_std_dev -
                                 1) < 1e-10,
                    "table spread does not match mode spread");
          tables[m * ns + s] = snapshot;
        }
    }
  // Commit only after all shape, routing, semantic and source validation
  // passes.
  auto committed_config = config;
  config_ = std::move(committed_config);
  parameters_ = std::move(p);
  tables_ = std::move(tables);
  fields_ = std::move(fields);
  solver_.reset();
  initialized_ = true;
}

const modal::Solver &ModalPackage::solver() const {
  if (!solver_)
    solver_ = modal::create(parameters_, tables_);
  return *solver_;
}

std::vector<ModalOpticalTableInfo> ModalPackage::opticalTables() const {
  std::vector<ModalOpticalTableInfo> result;
  for (std::size_t m = 0; m < config_.modes.size(); ++m)
    for (std::size_t s = 0; s < config_.species.size(); ++s) {
      const auto &t = tables_[m * config_.species.size() + s];
      if (t.kind)
        result.push_back({config_.modes[m].name, config_.species[s].name,
                          t.source_label, t.source_version, t.citation,
                          t.kind == 1 ? ModalTableKind::ParticleResolved
                                      : ModalTableKind::LognormalAverage});
    }
  return result;
}

unsigned ModalPackage::capabilities() const {
  unsigned result = Microphysics | Diagnostics | CCN | Emissions;
  if (!initialized_)
    return result;
  const auto ns = config_.species.size();
  for (std::size_t m = 0; m < config_.modes.size(); ++m) {
    if (config_.modes[m].optical_treatment ==
        ModalOpticalTreatment::InternalMixture) {
      if (!optics_provider_)
        return result;
    } else
      for (std::size_t s = 0; s < ns; ++s)
        if (parameters_.allowed[m * ns + s] && !tables_[m * ns + s].kind)
          return result;
  }
  return result | Optics;
}

void ModalPackage::validate(const EnvironmentalStateView &env,
                            const View3D<const double> &state,
                            bool full) const {
  require(initialized_, "initialize before compute");
  const auto nc = state.extent(0), nl = state.extent(1), nf = fields_.size();
  require(state.data_handle() && nc > 0 && nl > 0 && state.extent(2) == nf,
          "invalid packed state shape or pointer");
  require(nc <= static_cast<std::size_t>(std::numeric_limits<int>::max()) / nl /
                    nf,
          "state exceeds execution-policy limit");
  for (const auto &v : {env.temperature, env.relative_humidity})
    require(v.data_handle() && v.extent(0) == nc && v.extent(1) == nl,
            "environment shape or pointer mismatch");
  if (full)
    for (const auto &v : {env.pressure, env.air_density, env.layer_thickness})
      require(v.data_handle() && v.extent(0) == nc && v.extent(1) == nl,
              "environment shape or pointer mismatch");
  const auto ns = config_.species.size();
  for (std::size_t k = 0; k < nl; ++k)
    for (std::size_t c = 0; c < nc; ++c) {
      require(positive(env.temperature(c, k)) &&
                  std::isfinite(env.relative_humidity(c, k)) &&
                  env.relative_humidity(c, k) >= 0 &&
                  env.relative_humidity(c, k) <= 1,
              "invalid temperature or RH");
      if (full)
        require(positive(env.pressure(c, k)) &&
                    positive(env.air_density(c, k)) &&
                    std::isfinite(env.layer_thickness(c, k)) &&
                    env.layer_thickness(c, k) >= 0,
                "invalid pressure/density/layer thickness");
      for (std::size_t f = 0; f < nf; ++f)
        require(std::isfinite(state(c, k, f)) && state(c, k, f) >= 0,
                "state must be finite and nonnegative");
      for (std::size_t m = 0; m < config_.modes.size(); ++m) {
        double mass = 0, volume = 0, hyg = 0;
        for (std::size_t s = 0; s < ns; ++s) {
          const double a = state(c, k, m * (ns + 1) + s);
          require(parameters_.allowed[m * ns + s] || a == 0,
                  "disallowed material populated in a mode");
          mass += a;
          volume += a / config_.species[s].density;
          hyg += a / config_.species[s].density *
                 config_.species[s].hygroscopicity;
        }
        const double number = state(c, k, m * (ns + 1) + ns);
        require(std::isfinite(mass) && std::isfinite(volume) &&
                    std::isfinite(hyg) && std::isfinite(volume + 49 * hyg),
                "material volume or mass overflow");
        require(
            (mass > 0) == (number > 0),
            "each mode requires both positive mass and number, or both zero");
        if (number > 0)
          require(positive(volume / number) &&
                      positive(std::pow(std::cbrt(volume / number), 2)),
                  "unrepresentable particle volume");
      }
    }
}

void ModalPackage::executeMicrophysics(const EnvironmentalStateView &e,
                                       View3D<double> &state, double dt) {
  advanceWithCoating(e, state, dt, nullptr);
}
void ModalPackage::advanceWithCoating(const EnvironmentalStateView &e,
                                      View3D<double> &state, double dt,
                                      View3D<double> *gain) {
  const View3D<const double> input(state.data_handle(), state.extent(0),
                                   state.extent(1), state.extent(2));
  validate(e, input);
  require(std::isfinite(dt) && dt >= 0, "invalid timestep");
  if (gain)
    require(gain->data_handle() && gain->extents() == state.extents(),
            "coating gain shape mismatch");
  std::vector<double> tmp(state.data_handle(),
                          state.data_handle() + state.size()),
      gained(gain ? state.size() : 0, 0);
  if (dt > 0) {
    if (config_.coagulation &&
        config_.coagulation_model == ModalCoagulationModel::Mam4xx) {
#ifdef EXAERO_ENABLE_MAM4XX
      std::vector<double> native(tmp.size() * 2, 0),
          native_gain(gain ? native.size() : 0);
      std::copy(tmp.begin(), tmp.end(), native.begin());
      const Mam4xxEnvironment env{
          e.temperature.data_handle(), e.pressure.data_handle(),
          e.air_density.data_handle(), e.relative_humidity.data_handle(),
          e.layer_thickness.data_handle()};
      mam4xxMicrophysics(state.extent(0), state.extent(1), env, native.data(),
                         dt, gain ? native_gain.data() : nullptr);
      std::copy_n(native.begin(), tmp.size(), tmp.begin());
      if (gain)
        std::copy_n(native_gain.begin(), gained.size(), gained.begin());
#endif
    }
    const bool transfer =
        config_.transfer_owner == ModalTransferOwner::Microphysics;
    if (parameters_.coagulation ||
        (transfer && !parameters_.transfer_radius.empty()))
      modal::advance(solver(), state.extent(0), state.extent(1), raw(e),
                     tmp.data(), dt,
                     gain && parameters_.coagulation ? gained.data() : nullptr,
                     true, transfer);
    validate(e, View3D<const double>(tmp.data(), state.extent(0),
                                     state.extent(1), state.extent(2)));
  }
  for (double value : gained)
    require(std::isfinite(value) && value >= 0, "invalid coagulation gain");
  std::copy(tmp.begin(), tmp.end(), state.data_handle());
  if (gain)
    std::copy(gained.begin(), gained.end(), gain->data_handle());
}
void ModalPackage::renameModes(const EnvironmentalStateView &e,
                               View3D<double> &state) {
  require(config_.transfer_owner == ModalTransferOwner::Host,
          "renaming requires transfer_owner: host");
  const View3D<const double> input(state.data_handle(), state.extent(0),
                                   state.extent(1), state.extent(2));
  validate(e, input);
  std::vector<double> tmp(state.data_handle(),
                          state.data_handle() + state.size());
  if (!parameters_.transfer_radius.empty())
    modal::advance(solver(), state.extent(0), state.extent(1), raw(e),
                   tmp.data(), 0, nullptr, false, true);
  validate(e, View3D<const double>(tmp.data(), state.extent(0), state.extent(1),
                                   state.extent(2)));
  std::copy(tmp.begin(), tmp.end(), state.data_handle());
}

void ModalPackage::modalProperties(const EnvironmentalStateView &e,
                                   const View3D<const double> &state,
                                   View4D<double> &out,
                                   View4D<double> *jac) const {
  validate(e, state, false);
  require(!jac || config_.water_uptake != ModalWaterUptakeModel::Mam4xx,
          "native water property Jacobians are unavailable; dry uptake "
          "Jacobians remain supported");
  shape(out, state.extent(0), state.extent(1), config_.modes.size(),
        modal_properties::COUNT);
  if (jac)
    shape(*jac, state.extent(0), state.extent(1), config_.modes.size(),
          modal_properties::COUNT * fields_.size());
  modal::properties(solver(), state.extent(0), state.extent(1), raw(e, false),
                    state.data_handle(), out.data_handle(),
                    jac ? jac->data_handle() : nullptr);
}

void ModalPackage::modalUptake(const EnvironmentalStateView &e,
                               const View3D<const double> &state,
                               std::span<const ModalVapor> gases,
                               View4D<double> &out, View4D<double> *jac) const {
  validate(e, state, false);
  shape(out, state.extent(0), state.extent(1), config_.modes.size(),
        gases.size());
  if (jac)
    shape(*jac, state.extent(0), state.extent(1), config_.modes.size(),
          gases.size() * fields_.size());
  std::vector<modal::Vapor> v;
  for (const auto &g : gases) {
    require(positive(g.diffusivity) && positive(g.molecular_weight) &&
                positive(g.accommodation) && g.accommodation <= 1,
            "invalid vapor transport properties");
    v.push_back({g.diffusivity, g.molecular_weight, g.accommodation});
  }
  modal::uptake(solver(), state.extent(0), state.extent(1), raw(e, false),
                state.data_handle(), v, out.data_handle(),
                jac ? jac->data_handle() : nullptr);
}

void ModalPackage::computeDerivedDiagnostics(const EnvironmentalStateView &e,
                                             const View3D<const double> &state,
                                             View3D<double> &out) {
  validate(e, state);
  require(out.data_handle() && out.extent(0) == state.extent(0) &&
              out.extent(1) == state.extent(1) &&
              out.extent(2) == diagnostic_indices::NUM_DIAGNOSTICS,
          "diagnostic shape mismatch");
  require(state.extent(0) <=
              static_cast<std::size_t>(std::numeric_limits<int>::max()) /
                  state.extent(1) / diagnostic_indices::NUM_DIAGNOSTICS,
          "diagnostic slab exceeds execution-policy limit");
  modal::diagnostics(solver(), state.extent(0), state.extent(1), raw(e),
                     state.data_handle(), out.data_handle());
}

void ModalPackage::computeCCN(const EnvironmentalStateView &e,
                              const View3D<const double> &state,
                              const View1D<const double> &ss,
                              View4D<double> &out) {
  validate(e, state, false);
  shape(out, state.extent(0), state.extent(1), ss.extent(0), 1);
  require(ss.data_handle(), "null supersaturation axis");
  for (std::size_t i = 0; i < ss.extent(0); ++i)
    require(std::isfinite(ss(i)) && ss(i) >= 0, "invalid supersaturation");
  modal::ccn(solver(), state.extent(0), state.extent(1), raw(e, false),
             state.data_handle(), ss.data_handle(), ss.extent(0),
             out.data_handle());
}

void ModalPackage::computeEmissions(const EnvironmentalStateView &e,
                                    const EmissionsInputView &input,
                                    View3D<double> &out) {
  require(initialized_, "initialize before emissions");
  const auto &flux = input.flux;
  require(flux.data_handle() && out.data_handle() &&
              flux.extents() == out.extents() && flux.extent(0) > 0 &&
              flux.extent(1) > 0 && flux.extent(2) == fields_.size(),
          "emissions shape mismatch");
  require(flux.extent(0) <=
              static_cast<std::size_t>(std::numeric_limits<int>::max()) /
                  flux.extent(1) / fields_.size(),
          "emissions slab exceeds execution-policy limit");
  const bool area = input.flux_type == FluxType::AREA_FLUX;
  require(area || input.flux_type == FluxType::MASS_CONCENTRATION_RATE,
          "invalid flux convention");
  if (area)
    require(e.layer_thickness.data_handle() &&
                e.layer_thickness.extent(0) == flux.extent(0) &&
                e.layer_thickness.extent(1) == flux.extent(1),
            "emissions layer shape mismatch");
  std::vector<double> tmp(out.size());
  View3D<double> target(tmp.data(), out.extent(0), out.extent(1),
                        out.extent(2));
  for (std::size_t k = 0; k < flux.extent(1); ++k)
    for (std::size_t c = 0; c < flux.extent(0); ++c) {
      const double dz = area ? e.layer_thickness(c, k) : 1;
      require(positive(dz), "area flux requires positive layer thickness");
      for (std::size_t f = 0; f < fields_.size(); ++f) {
        require(std::isfinite(flux(c, k, f)) && flux(c, k, f) >= 0,
                "invalid emission rate");
        if (fields_[f].quantity == StateQuantity::Mass)
          require(parameters_.allowed[fields_[f].mode * config_.species.size() +
                                      f % (config_.species.size() + 1)] ||
                      flux(c, k, f) == 0,
                  "emissions populate disallowed material");
        target(c, k, f) = flux(c, k, f) / dz;
        require(std::isfinite(target(c, k, f)), "emission conversion overflow");
      }
    }
  std::copy(tmp.begin(), tmp.end(), out.data_handle());
}

int ModalPackage::getSpeciesIndex(const std::string &name) const {
  for (std::size_t i = 0; i < fields_.size(); ++i)
    if (fields_[i].name == name)
      return i;
  return -1;
}

double ModalPackage::emittedNumber(std::size_t m,
                                   std::span<const double> masses) const {
  require(initialized_ && m < config_.modes.size() &&
              masses.size() == config_.species.size(),
          "invalid emitted mode or mass vector");
  double volume = 0;
  for (std::size_t s = 0; s < masses.size(); ++s) {
    require(std::isfinite(masses[s]) && masses[s] >= 0 &&
                (parameters_.allowed[m * masses.size() + s] || masses[s] == 0),
            "invalid emitted material");
    volume += masses[s] / config_.species[s].density;
  }
  const auto &mode = config_.modes[m];
  const double ln = std::log(mode.geometric_std_dev);
  const double number =
      volume / (4 * pi / 3 * std::pow(mode.geometric_radius, 3) *
                std::exp(4.5 * ln * ln));
  require(std::isfinite(number), "emitted number overflow");
  return number;
}

ModalNucleationTerms ModalPackage::nucleationTerms(double acid) const {
  require(initialized_ && config_.nucleation.has_value(),
          "nucleation scheme is not configured");
  require(std::isfinite(acid) && acid >= 0,
          "invalid sulfuric acid concentration");
  const auto &n = *config_.nucleation;
  const auto it =
      std::find_if(config_.modes.begin(), config_.modes.end(),
                   [&](const auto &m) { return m.name == n.destination; });
  const auto dest = static_cast<int>(it - config_.modes.begin());
  ModalNucleationTerms terms;
  terms.destination = dest;
  terms.number_rate =
      acid == 0 || n.prefactor == 0
          ? 0
          : std::exp(std::log(n.prefactor) +
                     n.exponent * (std::log(acid) -
                                   std::log(n.precursor_molecular_weight) +
                                   std::log(6.02214076e23)));
  const double ln = std::log(it->geometric_std_dev);
  const double particle_volume =
      4 * pi / 3 * std::pow(n.seed_radius, 3) * std::exp(4.5 * ln * ln);
  double specific_volume = 0;
  for (std::size_t s = 0; s < config_.species.size(); ++s)
    specific_volume += n.mass_fractions[s] / config_.species[s].density;
  const double particle_mass = particle_volume / specific_volume;
  require(positive(particle_mass) && std::isfinite(terms.number_rate),
          "nucleation rate/seed overflow");
  for (double fraction : n.mass_fractions)
    terms.mass_rates.push_back(terms.number_rate * particle_mass * fraction);
  terms.precursor_loss = terms.number_rate * particle_mass;
  require(std::isfinite(terms.precursor_loss), "nucleation mass rate overflow");
  return terms;
}

void ModalPackage::computeOptics(const EnvironmentalStateView &e,
                                 const View3D<const double> &state,
                                 const View1D<const double> &bands,
                                 View4D<double> &out) {
  validate(e, state);
  shape(out, state.extent(0), state.extent(1), bands.extent(0),
        optical_indices::NUM_OPTICS);
  require(bands.data_handle(), "null wavelength axis");
  for (std::size_t b = 0; b < bands.extent(0); ++b)
    require(positive(bands(b)), "invalid wavelength");
  require(capabilities() & Optics,
          "optics requires compatible tables for every permitted component and "
          "a provider for internally mixed modes");
  const auto nc = state.extent(0), nl = state.extent(1),
             nm = config_.modes.size(), ns = config_.species.size(),
             nb = bands.extent(0);
  require(nc <= static_cast<std::size_t>(std::numeric_limits<int>::max()) / nl /
                    nm / nb / 9,
          "modal optical workspace is too large");
  std::vector<double> tmp(out.size()), mode_values(nc * nl * nm * nb * 9),
      prop(nc * nl * nm * modal_properties::COUNT);
  View4D<double> target(tmp.data(), nc, nl, nb, optical_indices::NUM_OPTICS),
      properties_view(prop.data(), nc, nl, nm, modal_properties::COUNT);
  modalProperties(e, state, properties_view, nullptr);
  modal::optics(solver(), nc, nl, raw(e), state.data_handle(),
                bands.data_handle(), nb, mode_values.data());
  std::vector<double> masses(ns), properties(modal_properties::COUNT);
  using namespace optical_indices;
  for (std::size_t b = 0; b < nb; ++b)
    for (std::size_t k = 0; k < nl; ++k)
      for (std::size_t c = 0; c < nc; ++c) {
        for (std::size_t m = 0; m < nm; ++m) {
          for (std::size_t s = 0; s < ns; ++s)
            masses[s] = state(c, k, m * (ns + 1) + s);
          for (int q = 0; q < modal_properties::COUNT; ++q)
            properties[q] = properties_view(c, k, m, q);
          if (properties[modal_properties::NUMBER] == 0)
            continue;
          ModalOpticalResult a;
          if (config_.modes[m].optical_treatment ==
              ModalOpticalTreatment::InternalMixture)
            a = optics_provider_->evaluate(config_, m, e.temperature(c, k),
                                           e.relative_humidity(c, k), bands(b),
                                           masses, properties);
          else {
            const auto v = [&](int q) {
              return mode_values[c + nc * (k + nl * (m + nm * (b + nb * q)))];
            };
            a = {v(0), v(1), v(2), v(3), v(4), v(5), v(6), v(7), v(8)};
          }
          addOptics(a, target, c, k, b);
          target(c, 0, b, FINE_MODE_EXTINCTION_AOT) +=
              a.fine_extinction * e.layer_thickness(c, k);
          target(c, 0, b, FINE_MODE_SCATTERING_AOT) +=
              a.fine_scattering * e.layer_thickness(c, k);
          target(c, 0, b, PM2_5_EXTINCTION_AOT) +=
              a.pm25_extinction * e.layer_thickness(c, k);
          target(c, 0, b, PM2_5_SCATTERING_AOT) +=
              a.pm25_scattering * e.layer_thickness(c, k);
        }
        if (target(c, k, b, SCATTERING_COEFF) > 0)
          target(c, k, b, ASYMMETRY_FACTOR) /=
              target(c, k, b, SCATTERING_COEFF);
        target(c, 0, b, EXTINCTION_AOT) +=
            target(c, k, b, EXTINCTION_COEFF) * e.layer_thickness(c, k);
        target(c, 0, b, SCATTERING_AOT) +=
            target(c, k, b, SCATTERING_COEFF) * e.layer_thickness(c, k);
      }
  for (double x : tmp)
    require(std::isfinite(x), "optical result overflow");
  std::copy(tmp.begin(), tmp.end(), out.data_handle());
}
} // namespace exaero
