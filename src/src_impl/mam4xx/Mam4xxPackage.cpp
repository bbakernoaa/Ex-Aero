/// @file Mam4xxPackage.cpp
/// @brief Validated host interface to the isolated MAM4xx device solver.
#include <cmath>
#include <exaero/Mam4xxConfig.hpp>
#include <limits>
#include <mam4xx/Mam4xxPackage.hpp>
#include <mam4xx/Mam4xxSolver.hpp>
#include <yaml-cpp/yaml.h>

namespace exaero {
namespace {
/// @brief Reject invalid inputs before any output or state is modified.
void require(bool ok, const char *message) {
  if (!ok)
    throw std::invalid_argument(std::string("FATAL ERROR: MAM4xx: ") + message);
}
/// @brief Extract contiguous host storage without exposing mdspan to the
/// solver.
Mam4xxEnvironment raw(const EnvironmentalStateView &env) {
  return {env.temperature.data_handle(), env.pressure.data_handle(),
          env.air_density.data_handle(), env.relative_humidity.data_handle(),
          env.layer_thickness.data_handle()};
}
} // namespace

std::vector<StateField> Mam4xxPackage::stateFields() const {
  std::vector<StateField> result;
  for (int phase = 0; phase < 2; ++phase)
    for (int mode = 0; mode < mam4xx_layout::modes; ++mode)
      for (int species = 0; species <= mam4xx_layout::species; ++species) {
        const bool number = species == mam4xx_layout::species;
        result.push_back(
            {std::string(number ? "num"
                                : mam4xx_layout::species_names[species]) +
                 "_" + std::string(mam4xx_layout::mode_names[mode]) +
                 (phase ? "_c" : "_i"),
             number ? StateQuantity::Number : StateQuantity::Mass, mode,
             phase != 0});
      }
  return result;
}

void Mam4xxPackage::initialize(const std::string &config_yaml) {
  const auto config = YAML::Load(config_yaml);
  require(!config || config.IsNull() || config.IsMap(),
          "configuration must be a mapping");
  if (config && config.IsMap()) {
    for (const auto &entry : config) {
      const auto key = entry.first.as<std::string>();
      require(key == "backend", "unknown configuration key");
      require(entry.second.as<std::string>() == "mam4xx",
              "backend must be mam4xx");
    }
  }
  initialized_ = true;
}

void Mam4xxPackage::validate(const EnvironmentalStateView &env,
                             const View3D<const double> &state,
                             bool full_environment) const {
  require(initialized_, "initialize the package before computing");
  const auto cells = state.extent(0), levels = state.extent(1);
  require(cells > 0 && levels > 0 && state.data_handle(),
          "empty or null state");
  require(cells <= static_cast<std::size_t>(std::numeric_limits<int>::max()) /
                       levels / mam4xx_layout::fields,
          "grid is too large for the execution policy");
  require(state.extent(2) == mam4xx_layout::fields, "state requires 64 fields");
  for (const auto &field : {env.temperature, env.relative_humidity})
    require(field.data_handle() && field.extent(0) == cells &&
                field.extent(1) == levels,
            "environment extent mismatch or null data");
  if (full_environment)
    for (const auto &field :
         {env.pressure, env.air_density, env.layer_thickness})
      require(field.data_handle() && field.extent(0) == cells &&
                  field.extent(1) == levels,
              "environment extent mismatch or null data");
  for (std::size_t k = 0; k < levels; ++k)
    for (std::size_t c = 0; c < cells; ++c) {
      require(std::isfinite(env.temperature(c, k)) && env.temperature(c, k) > 0,
              "temperature must be finite and positive");
      if (full_environment) {
        for (const double value : {env.pressure(c, k), env.air_density(c, k)})
          require(
              std::isfinite(value) && value > 0,
              "temperature, pressure and density must be finite and positive");
        require(std::isfinite(env.layer_thickness(c, k)) &&
                    env.layer_thickness(c, k) >= 0,
                "layer thickness must be finite and nonnegative");
      }
      require(std::isfinite(env.relative_humidity(c, k)) &&
                  env.relative_humidity(c, k) >= 0 &&
                  env.relative_humidity(c, k) <= 1,
              "RH must be finite and in [0,1]");
      for (int f = 0; f < mam4xx_layout::fields; ++f)
        require(std::isfinite(state(c, k, f)) && state(c, k, f) >= 0,
                "state must be finite and nonnegative");
      for (int phase = 0; phase < 2; ++phase)
        for (int m = 0; m < mam4xx_layout::modes; ++m) {
          double mass = 0;
          for (int s = 0; s < mam4xx_layout::species; ++s)
            mass += state(c, k, mam4xx_layout::mass(m, s, phase));
          require(std::isfinite(mass), "modal mass overflows");
          require((mass > 0) ==
                      (state(c, k, mam4xx_layout::number(m, phase)) > 0),
                  "a populated mode requires both positive mass and number");
        }
    }
}

void Mam4xxPackage::executeMicrophysics(const EnvironmentalStateView &env,
                                        View3D<double> &state, double dt) {
  const View3D<const double> input(state.data_handle(), state.extent(0),
                                   state.extent(1), state.extent(2));
  validate(env, input);
  require(std::isfinite(dt) && dt >= 0,
          "timestep must be finite and nonnegative");
  if (dt == 0)
    return;
  mam4xxMicrophysics(state.extent(0), state.extent(1), raw(env),
                     state.data_handle(), dt);
}

void Mam4xxPackage::computeDerivedDiagnostics(const EnvironmentalStateView &env,
                                              const View3D<const double> &state,
                                              View3D<double> &output) {
  validate(env, state);
  require(output.data_handle() && output.extent(0) == state.extent(0) &&
              output.extent(1) == state.extent(1) &&
              output.extent(2) == diagnostic_indices::NUM_DIAGNOSTICS,
          "diagnostic extent mismatch or null output");
  mam4xxDiagnostics(state.extent(0), state.extent(1), raw(env),
                    state.data_handle(), output.data_handle());
}

void Mam4xxPackage::computeOptics(const EnvironmentalStateView &,
                                  const View3D<const double> &,
                                  const View1D<const double> &,
                                  View4D<double> &) {
  throw std::logic_error("FATAL ERROR: MAM4xx optics requires internally mixed "
                         "modal optical tables; unavailable in this backend");
}

void Mam4xxPackage::computeCCN(const EnvironmentalStateView &env,
                               const View3D<const double> &state,
                               const View1D<const double> &supersaturations,
                               View4D<double> &output) {
  validate(env, state, false);
  require(supersaturations.data_handle() && supersaturations.extent(0) > 0,
          "empty supersaturation axis");
  require(supersaturations.extent(0) <=
              static_cast<std::size_t>(std::numeric_limits<int>::max()) /
                  state.extent(0) / state.extent(1),
          "CCN slab is too large");
  require(output.data_handle() && output.extent(0) == state.extent(0) &&
              output.extent(1) == state.extent(1) &&
              output.extent(2) == supersaturations.extent(0) &&
              output.extent(3) == 1,
          "CCN extent mismatch or null output");
  for (std::size_t i = 0; i < supersaturations.extent(0); ++i)
    require(std::isfinite(supersaturations(i)) && supersaturations(i) >= 0,
            "supersaturations must be finite and nonnegative");
  mam4xxCCN(state.extent(0), state.extent(1), raw(env), state.data_handle(),
            supersaturations.data_handle(), supersaturations.extent(0),
            output.data_handle());
}

void Mam4xxPackage::computeEmissions(const EnvironmentalStateView &env,
                                     const EmissionsInputView &input,
                                     View3D<double> &output) {
  require(initialized_, "initialize the package before computing");
  const auto &flux = input.flux;
  require(flux.data_handle() && output.data_handle() && flux.extent(0) > 0 &&
              flux.extent(1) > 0 && flux.extent(2) == mam4xx_layout::fields &&
              flux.extents() == output.extents(),
          "emissions require matching 64-field arrays");
  const bool area = input.flux_type == FluxType::AREA_FLUX;
  require(area || input.flux_type == FluxType::MASS_CONCENTRATION_RATE,
          "invalid flux convention");
  if (area)
    require(env.layer_thickness.data_handle() &&
                env.layer_thickness.extent(0) == flux.extent(0) &&
                env.layer_thickness.extent(1) == flux.extent(1),
            "emissions layer extent mismatch");
  // Validate the whole input before writing any output (also permits in-place
  // mapping).
  for (std::size_t k = 0; k < flux.extent(1); ++k)
    for (std::size_t c = 0; c < flux.extent(0); ++c) {
      if (area)
        require(std::isfinite(env.layer_thickness(c, k)) &&
                    env.layer_thickness(c, k) > 0,
                "area flux requires positive finite layer thickness");
      for (int f = 0; f < mam4xx_layout::fields; ++f)
        require(std::isfinite(flux(c, k, f)) && flux(c, k, f) >= 0,
                "invalid emissions flux");
    }
  for (std::size_t k = 0; k < flux.extent(1); ++k)
    for (std::size_t c = 0; c < flux.extent(0); ++c)
      for (int f = 0; f < mam4xx_layout::fields; ++f)
        output(c, k, f) =
            flux(c, k, f) / (area ? env.layer_thickness(c, k) : 1.0);
}

int Mam4xxPackage::getSpeciesIndex(const std::string &name) const {
  const auto fields = stateFields();
  for (std::size_t i = 0; i < fields.size(); ++i)
    if (fields[i].name == name)
      return static_cast<int>(i);
  return -1;
}
std::string Mam4xxPackage::getSpeciesName(int index) const {
  return stateFields().at(index).name;
}
} // namespace exaero
