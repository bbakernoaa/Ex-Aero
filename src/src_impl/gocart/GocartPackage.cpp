#include <gocart/GocartPackage.hpp>
#include <loader/MieTableStore.hpp>
#include <yaml-cpp/yaml.h>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace exaero {

// Semantic view over a package configuration (YAML node or GocartConfig struct).
// initializeImpl() consumes ONLY this interface — it never touches YAML nodes or
// GocartConfig fields directly, so both entry points share one code path and the
// precedence order (reset -> file -> config) cannot drift between them.
// Impl-private: declared in GocartPackage.hpp (elaborated type-specifier), defined
// here at exaero:: scope so the member declaration binds to exactly this type.
struct PackageConfigReader {
    virtual ~PackageConfigReader() = default;
    virtual int numSpecies() const = 0;
    virtual std::string speciesName(int i) const = 0;
    virtual GocartSpeciesParams speciesScalars(int i) const = 0; // 8 doubles + has_optics_lookup; offsets filled by impl
    virtual bool hasOpticsLookup(int i) const = 0;
    virtual GocartLegacyOpticsLookup opticsLookup(int i) const = 0;      // validated: n>=2, equal lengths
    virtual bool hasMieTable(int i) const = 0;
    virtual SpeciesCurveConfig mieTable(int i) const = 0;                // species_name pre-filled
    virtual bool hasActivation() const = 0;
    virtual int activationExtraMask() const = 0;                          // category bits OR-ed on default
    virtual std::vector<std::string> activationSpecies() const = 0;
    virtual std::string activationDataFile() const = 0;
    virtual int numEmissionsMappings() const = 0;
    virtual GocartEmissionsMappingConfig emissionsMapping(int s) const = 0; // raw_name ignored downstream (parity)
};

namespace {
// Map a human-readable attribute name (YAML override key / index code) to the
// (category, index) pair the store resolves. Names follow data-model.md; the index
// codes mirror the *_indices constants in AttributeQuery.hpp. Returns false for an
// unknown name so the caller can fail fast rather than silently drop an override.
bool attribute_key(const std::string& name, AttributeCategory& cat, int& idx) {
    using namespace exaero::microphysical_indices;
    using namespace exaero::spectral_optical_indices;
    static const std::vector<std::pair<std::string, std::pair<AttributeCategory, int>>> table = {
        {"wet_particle_density", {AttributeCategory::Microphysical, WET_PARTICLE_DENSITY}},
        {"growth_factor",        {AttributeCategory::Microphysical, GROWTH_FACTOR}},
        {"effective_radius",     {AttributeCategory::Microphysical, EFFECTIVE_RADIUS}},
        {"mass_mean_radius",     {AttributeCategory::Microphysical, MASS_MEAN_RADIUS}},
        {"bin_lower_radius",     {AttributeCategory::Microphysical, BIN_LOWER_RADIUS}},
        {"bin_upper_radius",     {AttributeCategory::Microphysical, BIN_UPPER_RADIUS}},
        {"volume_per_mass",      {AttributeCategory::Microphysical, VOLUME_PER_MASS}},
        {"area_per_mass",        {AttributeCategory::Microphysical, AREA_PER_MASS}},
        {"particle_mass",        {AttributeCategory::Microphysical, PARTICLE_MASS}},
        {"qext", {AttributeCategory::SpectralOptical, EXTINCTION_EFFICIENCY}},
        {"qsca", {AttributeCategory::SpectralOptical, SCATTERING_EFFICIENCY}},
        {"qabs", {AttributeCategory::SpectralOptical, ABSORPTION_EFFICIENCY}},
        {"bext", {AttributeCategory::SpectralOptical, MASS_EXTINCTION}},
        {"bsca", {AttributeCategory::SpectralOptical, MASS_SCATTERING}},
        {"bbck", {AttributeCategory::SpectralOptical, MASS_BACKSCATTER}},
        {"lidar_ratio", {AttributeCategory::SpectralOptical, LIDAR_RATIO}},
        {"g", {AttributeCategory::SpectralOptical, ASYMMETRY_FACTOR}},
        {"ssa", {AttributeCategory::SpectralOptical, SINGLE_SCATTERING_ALBEDO}},
        {"refreal", {AttributeCategory::SpectralOptical, REFRACTIVE_INDEX_REAL}},
        {"refimag", {AttributeCategory::SpectralOptical, REFRACTIVE_INDEX_IMAG}},
    };
    for (const auto& [key, val] : table) {
        if (key == name) { cat = val.first; idx = val.second; return true; }
    }
    return false;
}

// Parse a top-level activation category token into a bitmask of (1 << AttributeCategory).
int category_mask_from_yaml(const YAML::Node& node) {
    int mask = 0;
    if (!node) return mask;
    if (node.IsSequence()) {
        for (const auto& c : node) {
            const std::string s = c.as<std::string>();
            if (s == "microphysical") mask |= attribute_category_bit(AttributeCategory::Microphysical);
            else if (s == "spectral") mask |= attribute_category_bit(AttributeCategory::SpectralOptical);
            else if (s == "polarized") mask |= attribute_category_bit(AttributeCategory::PolarizedMoment);
            else throw std::runtime_error("EX-aero Error: unknown activation category '" + s + "'");
        }
    }
    return mask;
}

// PackageConfigReader over a parsed YAML document. Each method is a direct port of
// the inline YAML access initialize() used before the seam refactor; error strings
// and fail-fast checks are preserved verbatim, and missing-scalar `.as<T>()`
// exceptions propagate unchanged (no catch, no softening).
class YamlReader final : public PackageConfigReader {
public:
    explicit YamlReader(YAML::Node config) : config_(std::move(config)) {
        if (!config_["species"]) {
            throw std::runtime_error("GOCART YAML config missing 'species' field");
        }
        species_ = config_["species"];
    }

    int numSpecies() const override { return static_cast<int>(species_.size()); }

    std::string speciesName(int i) const override {
        auto s = species_[i];
        return s["name"] ? s["name"].as<std::string>() : ("Species_" + std::to_string(i));
    }

    GocartSpeciesParams speciesScalars(int i) const override {
        auto s = species_[i];
        GocartSpeciesParams p{
            s["dry_density"].as<double>(),
            s["molecular_weight"].as<double>(),
            s["dry_particle_diameter"].as<double>(),
            s["hygroscopicity"].as<double>(),
            s["lognormal_sigma"].as<double>(),
            s["lognormal_dg"].as<double>(),
            s["refractive_index_real"].as<double>(),
            s["refractive_index_imag"].as<double>(),
        };
        p.has_optics_lookup = s["has_optics_lookup"] && s["has_optics_lookup"].as<bool>();
        return p;
    }

    bool hasOpticsLookup(int i) const override {
        auto s = species_[i];
        return s["has_optics_lookup"] && s["has_optics_lookup"].as<bool>();
    }

    GocartLegacyOpticsLookup opticsLookup(int i) const override {
        auto s = species_[i];
        auto rh_bins = s["rh_bins"];
        auto ext_lookup = s["ext_lookup"];
        auto ssa_lookup = s["ssa_lookup"];
        auto asm_lookup = s["asm_lookup"];
        if (!rh_bins || !ext_lookup || !ssa_lookup || !asm_lookup) {
            throw std::runtime_error("EX-aero Error: species '" + speciesName(i) +
                "' sets has_optics_lookup but omits a required lookup list");
        }
        const int n_rh = static_cast<int>(rh_bins.size());
        if (n_rh < 2) {
            throw std::runtime_error("EX-aero Error: species '" + speciesName(i) +
                "' rh_bins needs at least 2 points to interpolate");
        }
        if (static_cast<int>(ext_lookup.size()) != n_rh ||
            static_cast<int>(ssa_lookup.size()) != n_rh ||
            static_cast<int>(asm_lookup.size()) != n_rh) {
            throw std::runtime_error("EX-aero Error: species '" + speciesName(i) +
                "' lookup lists have mismatched lengths (fail fast, no silent fallback)");
        }
        GocartLegacyOpticsLookup lut;
        for (int j = 0; j < n_rh; ++j) lut.rh.push_back(rh_bins[j].as<double>());
        for (int j = 0; j < n_rh; ++j) lut.ext.push_back(ext_lookup[j].as<double>());
        for (int j = 0; j < n_rh; ++j) lut.ssa.push_back(ssa_lookup[j].as<double>());
        for (int j = 0; j < n_rh; ++j) lut.asm_.push_back(asm_lookup[j].as<double>());
        return lut;
    }

    bool hasMieTable(int i) const override {
        auto s = species_[i];
        return static_cast<bool>(s["mie_table"]);
    }

    SpeciesCurveConfig mieTable(int i) const override {
        auto s = species_[i];
        auto mt = s["mie_table"];
        SpeciesCurveConfig cc;
        cc.species_name = speciesName(i);
        if (mt["source"]) cc.source_label = mt["source"].as<std::string>();
        if (mt["interpolate"]) cc.interpolate = mt["interpolate"].as<std::string>();
        if (mt["solver_radius_node"]) cc.solver_radius_node = mt["solver_radius_node"].as<int>();
        if (mt["radius_nodes"]) {
            for (const auto& rn : mt["radius_nodes"]) cc.radius_nodes.push_back(rn.as<double>());
        }
        if (mt["overrides"]) {
            for (YAML::const_iterator it = mt["overrides"].begin(); it != mt["overrides"].end(); ++it) {
                const std::string attr = it->first.as<std::string>();
                AttributeCategory cat; int idx;
                if (!attribute_key(attr, cat, idx)) {
                    throw std::runtime_error("EX-aero Error: species '" + cc.species_name +
                        "' override names unknown attribute '" + attr + "' (fail fast)");
                }
                SpeciesCurveConfig::Override ov;
                ov.category = cat; ov.attribute_index = idx;
                for (const auto& v : it->second) ov.values.push_back(v.as<double>());
                cc.overrides.push_back(std::move(ov));
            }
        }
        return cc;
    }

    bool hasActivation() const override { return static_cast<bool>(config_["activation"]); }

    int activationExtraMask() const override {
        auto act = config_["activation"];
        return act["categories"] ? category_mask_from_yaml(act["categories"]) : 0;
    }

    std::vector<std::string> activationSpecies() const override {
        std::vector<std::string> activated_species;
        auto act = config_["activation"];
        if (act["species"]) {
            for (const auto& sp : act["species"]) activated_species.push_back(sp.as<std::string>());
        }
        return activated_species;
    }

    std::string activationDataFile() const override {
        auto act = config_["activation"];
        return act["data_file"] ? act["data_file"].as<std::string>() : std::string{};
    }

    int numEmissionsMappings() const override {
        auto mapping_node = config_["emissions_mapping"];
        return mapping_node ? static_cast<int>(mapping_node.size()) : 0;
    }

    GocartEmissionsMappingConfig emissionsMapping(int s) const override {
        auto mapping_node = config_["emissions_mapping"];
        GocartEmissionsMappingConfig em;
        em.raw_name = mapping_node[s]["raw_name"] ? mapping_node[s]["raw_name"].as<std::string>()
                                                  : "CECE_Raw_Species";
        auto mappings = mapping_node[s]["mappings"];
        if (mappings) {
            for (size_t m = 0; m < mappings.size(); ++m) {
                GocartEmissionsMappingConfig::Mapping mp;
                mp.target_species = mappings[m]["target_species"].as<std::string>();
                mp.mass_split_fraction = mappings[m]["mass_split_fraction"].as<double>();
                if (mappings[m]["is_modal_mode"] && mappings[m]["is_modal_mode"].as<bool>()) {
                    mp.is_modal_mode = true;
                    mp.emitted_particle_diameter = mappings[m]["emitted_particle_diameter"].as<double>();
                    mp.lognormal_sigma = mappings[m]["lognormal_sigma"].as<double>();
                }
                em.mappings.push_back(std::move(mp));
            }
        }
        return em;
    }

private:
    YAML::Node config_, species_;
};

// Reader over an in-memory GocartConfig (design 2026-09-10 §4). Near-trivial by
// design: the caller already supplied typed data, so there are no names to
// validate (curve overrides carry typed (category, attribute_index) pairs).
// Validation duplicates the YAML fail-fast guarantees with the same
// "EX-aero Error:" style (never silent): non-empty species, non-empty names,
// and legacy-lookup length consistency.
class StructReader final : public PackageConfigReader {
public:
    explicit StructReader(const GocartConfig& config) : config_(config) {
        if (config_.species.empty()) {
            throw std::runtime_error("EX-aero Error: structured config has no species");
        }
        for (std::size_t i = 0; i < config_.species.size(); ++i) {
            if (config_.species[i].name.empty()) {
                throw std::runtime_error("EX-aero Error: structured config species at index " +
                    std::to_string(i) + " has an empty name");
            }
        }
    }
    int numSpecies() const override { return static_cast<int>(config_.species.size()); }
    std::string speciesName(int i) const override { return config_.species[i].name; }
    GocartSpeciesParams speciesScalars(int i) const override {
        const auto& s = config_.species[i];
        GocartSpeciesParams p{s.dry_density, s.molecular_weight, s.dry_particle_diameter,
                              s.hygroscopicity, s.lognormal_sigma, s.lognormal_dg,
                              s.refractive_index_real, s.refractive_index_imag};
        p.has_optics_lookup = s.optics_lookup.has_value();
        return p;
    }
    bool hasOpticsLookup(int i) const override {
        return config_.species[i].optics_lookup.has_value();
    }
    GocartLegacyOpticsLookup opticsLookup(int i) const override {
        const auto& l = *config_.species[i].optics_lookup; // caller checked engaged
        const int n = static_cast<int>(l.rh.size());
        if (n < 2) {
            throw std::runtime_error("EX-aero Error: species '" + speciesName(i) +
                "' rh_bins needs at least 2 points to interpolate");
        }
        if (static_cast<int>(l.ext.size()) != n || static_cast<int>(l.ssa.size()) != n ||
            static_cast<int>(l.asm_.size()) != n) {
            throw std::runtime_error("EX-aero Error: species '" + speciesName(i) +
                "' lookup lists have mismatched lengths (fail fast, no silent fallback)");
        }
        return l;
    }
    bool hasMieTable(int i) const override { return config_.species[i].mie_table.has_value(); }
    SpeciesCurveConfig mieTable(int i) const override {
        SpeciesCurveConfig cc = *config_.species[i].mie_table; // caller checked engaged
        cc.species_name = config_.species[i].name;            // config name wins (R9)
        // Fail-fast parity with the YAML path: an override row must supply exactly one
        // value per declared radius node. When radius_nodes is empty the effective axis
        // is the SOURCE axis, which StructReader cannot see; the identical guarantee is
        // enforced downstream in MieTableStore::apply_curve_config against the effective
        // n_radius (same as YAML), so there is no divergence between the two paths.
        if (!cc.radius_nodes.empty()) {
            const std::size_t n_nodes = cc.radius_nodes.size();
            for (const auto& ov : cc.overrides) {
                if (ov.values.size() != n_nodes) {
                    throw std::runtime_error("EX-aero Error: species '" + cc.species_name +
                        "' override for attribute " + std::to_string(ov.attribute_index) +
                        " must supply one value per radius_node (fail fast)");
                }
            }
        }
        return cc;
    }
    bool hasActivation() const override { return config_.activation.has_value(); }
    int activationExtraMask() const override {
        int mask = 0;
        for (auto c : config_.activation->categories) mask |= attribute_category_bit(c);
        return mask;
    }
    std::vector<std::string> activationSpecies() const override {
        return config_.activation ? config_.activation->species : std::vector<std::string>{};
    }
    std::string activationDataFile() const override {
        return config_.activation ? config_.activation->data_file : std::string();
    }
    int numEmissionsMappings() const override {
        return static_cast<int>(config_.emissions_mapping.size());
    }
    GocartEmissionsMappingConfig emissionsMapping(int s) const override {
        return config_.emissions_mapping[s];
    }
private:
    const GocartConfig& config_;
};
} // namespace


    // Decoupled solver helper function declarations (implemented in GocartSolver.cpp)
    GocartSolverState* create_solver_state(int num_species, const GocartSpeciesParams* params,
                                           const double* pool, int pool_size);
    void free_solver_state(GocartSolverState* state);

    void run_gocart_diagnostics(
        GocartSolverState* state,
        int num_cells, int num_levels, int num_species,
        const double* rh_ptr, const double* thick_ptr, const double* state_ptr, double* diags_ptr
    );

    void run_gocart_optics(
        GocartSolverState* state,
        int num_cells, int num_levels, int num_bands, int num_species,
        const double* wavelengths_ptr, const double* rh_ptr, const double* thick_ptr,
        const double* state_ptr, double* optics_ptr
    );

    void run_gocart_ccn(
        GocartSolverState* state,
        int num_cells, int num_levels, int num_ss, int num_species,
        const double* ss_ptr, const double* temp_ptr, const double* rh_ptr,
        const double* state_ptr, double* ccn_ptr
    );

    void run_gocart_emissions(
        GocartSolverState* state,
        int num_cells, int num_levels, int num_raw_species, int num_target_species,
        int flux_type_code,
        const double* thick_ptr,
        const double* raw_emissions_ptr,
        double* target_emissions_out_ptr
    );

    GocartPackage::GocartPackage() : solver_state_(nullptr) {}

    GocartPackage::~GocartPackage() {
        if (solver_state_) {
            free_solver_state(solver_state_);
        }
    }

    // Shared orchestration for both configuration entry points (design 2026-09-10,
    // Approach A): species params -> store reset -> activation -> curve configs ->
    // pool wiring -> emissions -> solver-state upload. Reads config ONLY through the
    // PackageConfigReader interface, so the YAML path and the structured-config path
    // (Task 3) share one code path and the precedence order cannot drift between them.
    void GocartPackage::initializeImpl(PackageConfigReader& reader) {
        num_species_ = reader.numSpecies();

        h_species_params_.clear();
        h_species_params_.reserve(num_species_);
        species_names_.clear();
        species_names_.reserve(num_species_);
        h_curve_pool_.clear(); // rebuilt from scratch each initialize (single H2D upload)

        for (int i = 0; i < num_species_; ++i) {
            // Extract the dynamic species name
            std::string name = reader.speciesName(i);
            species_names_.push_back(name);

            GocartSpeciesParams p = reader.speciesScalars(i);

            if (p.has_optics_lookup) {
                // Legacy RH lookups (presence/length checks already validated by the
                // reader): arbitrary length, appended to the flat pool.
                // The block layout is the documented {rh,ext,ssa,asm} x n_rh order (R10).
                const GocartLegacyOpticsLookup lut = reader.opticsLookup(i);
                const int n_rh = static_cast<int>(lut.rh.size());
                p.n_radius = 1;   // band-integrated legacy lookup: single effective bin
                p.n_rh = n_rh;
                p.n_lambda = 0;
                p.curve_offset = static_cast<int>(h_curve_pool_.size());
                for (int j = 0; j < n_rh; ++j) h_curve_pool_.push_back(lut.rh[j]);
                for (int j = 0; j < n_rh; ++j) h_curve_pool_.push_back(lut.ext[j]);
                for (int j = 0; j < n_rh; ++j) h_curve_pool_.push_back(lut.ssa[j]);
                for (int j = 0; j < n_rh; ++j) h_curve_pool_.push_back(lut.asm_[j]);
            }
            h_species_params_.push_back(p);
        }

        // --- GEOSmie MIE curve mapping (ADR-003 / research R9): parse per-species
        //     `mie_table:` block into SpeciesCurveConfig, then resolve through the store.
        //     Order matters: reset -> optional runtime file (activation) -> config, so
        //     precedence is config > runtime-file > baked-in. ---
        MieTableStore& store = MieTableStore::instance();
        store.reset_to_baked_in();

        // Top-level activation (optional): species list, categories, runtime data file.
        // The default mask (microphysical + spectral-optical) is the impl's, not the
        // reader's: the reader only supplies extra category bits OR-ed onto it.
        int activation_mask = attribute_category_bit(AttributeCategory::Microphysical) |
                              attribute_category_bit(AttributeCategory::SpectralOptical);
        std::vector<std::string> activated_species;
        std::string runtime_file;
        if (reader.hasActivation()) {
            activation_mask |= reader.activationExtraMask();
            activated_species = reader.activationSpecies();
            runtime_file = reader.activationDataFile();
        }
        // Runtime file (if any) is applied through setAttributeActivation (FR-012); a bad
        // file aborts with FATAL ERROR and no fallback (FR-009).
        setAttributeActivation(activated_species, activation_mask, runtime_file);

        std::vector<SpeciesCurveConfig> curve_configs; // per-species mie_table: bindings (R9)
        for (int i = 0; i < num_species_; ++i) {
            if (!reader.hasMieTable(i)) continue;
            curve_configs.push_back(reader.mieTable(i));
        }
        if (!curve_configs.empty()) setSpeciesCurveConfig(curve_configs);

        // --- Device-resident microphysical curves (T021, ADR-003 R10): append each
        //     curve-bound species' {rh, growth_factor, wet_particle_density} block to the
        //     same flat pool the optics kernels use. Extents are data, never literals;
        //     the diagnostics kernel reads them on-device with zero H2D in the loop. ---
        for (int i = 0; i < num_species_; ++i) {
            const SpeciesCurve* c = store.find_curve(species_names_[i]);
            if (!c) continue;
            auto git = c->fields.find("growth_factor");
            auto dit = c->fields.find("wet_particle_density");
            if (git == c->fields.end() || dit == c->fields.end()) continue;
            const int nR = c->n_radius();
            const int nH = c->n_rh();
            if (static_cast<int>(git->second.values.size()) != nR * nH ||
                static_cast<int>(dit->second.values.size()) != nR * nH) {
                throw std::runtime_error("EX-aero Error: species '" + species_names_[i] +
                    "' microphysical field length does not match its axes (fail fast)");
            }
            auto& p = h_species_params_[i];
            p.micro_offset = static_cast<int>(h_curve_pool_.size());
            p.n_micro_radius = nR;
            p.n_micro_rh = nH;
            p.solver_radius_node = c->solver_radius_node;
            h_curve_pool_.insert(h_curve_pool_.end(), c->rh.begin(), c->rh.end());
            h_curve_pool_.insert(h_curve_pool_.end(), git->second.values.begin(), git->second.values.end());
            h_curve_pool_.insert(h_curve_pool_.end(), dit->second.values.begin(), dit->second.values.end());
        }

        // --- Device-resident spectral curves (T028, ADR-003 R10): append each curve-bound
        //     species' {rh, lambda, bext, ssa, g} block so the optics kernels read mass
        //     extinction / albedo / asymmetry from the table (RH-linear, log-band) instead
        //     of the ADT analytical solver. Extents are data; a band coordinate outside the
        //     curve lambda domain declines to the fallback path (never a silent wrong value). ---
        for (int i = 0; i < num_species_; ++i) {
            const SpeciesCurve* c = store.find_curve(species_names_[i]);
            if (!c) continue;
            auto bit = c->fields.find("bext");
            auto sit = c->fields.find("ssa");
            auto git = c->fields.find("g");
            if (bit == c->fields.end() || sit == c->fields.end() || git == c->fields.end()) continue;
            const int nR = c->n_radius();
            const int nH = c->n_rh();
            const int nL = c->n_lambda();
            const std::size_t need = static_cast<std::size_t>(nR) * nH * nL;
            if (c->lambda.size() != static_cast<std::size_t>(nL) ||
                bit->second.values.size() != need || sit->second.values.size() != need ||
                git->second.values.size() != need) {
                throw std::runtime_error("EX-aero Error: species '" + species_names_[i] +
                    "' spectral field length does not match its axes (fail fast)");
            }
            auto& p = h_species_params_[i];
            p.spec_offset = static_cast<int>(h_curve_pool_.size());
            p.n_spec_radius = nR;
            p.n_spec_rh = nH;
            p.n_spec_lambda = nL;
            h_curve_pool_.insert(h_curve_pool_.end(), c->rh.begin(), c->rh.end());
            h_curve_pool_.insert(h_curve_pool_.end(), c->lambda.begin(), c->lambda.end());
            h_curve_pool_.insert(h_curve_pool_.end(), bit->second.values.begin(), bit->second.values.end());
            h_curve_pool_.insert(h_curve_pool_.end(), sit->second.values.begin(), sit->second.values.end());
            h_curve_pool_.insert(h_curve_pool_.end(), git->second.values.begin(), git->second.values.end());
        }

        // --- Device-resident polarized moments (T033, FR-003/FR-015, ADR-003 R10): append
        //     each moment-bearing species' {rh, lambda, pmom} block to the same single-upload
        //     pool. Extents (pol, moment) are data from the curve (never literals); a species
        //     without a pmom field gets no block (pmom_offset stays -1 => FR-008 not-in-source).
        //     Hot-path consumers read stride-1 over the fastest moment axis. ---
        for (int i = 0; i < num_species_; ++i) {
            const SpeciesCurve* c = store.find_curve(species_names_[i]);
            if (!c) continue;
            auto pit = c->fields.find("pmom");
            if (pit == c->fields.end()) continue;
            const int nR = c->n_radius();
            const int nH = c->n_rh();
            const int nL = c->n_lambda();
            const int nP = c->n_pol();
            const int nM = c->n_moment();
            const std::size_t need =
                static_cast<std::size_t>(nR) * nH * nL * nP * nM;
            if (nP <= 0 || nM <= 0 ||
                pit->second.values.size() != need ||
                c->lambda.size() != static_cast<std::size_t>(nL)) {
                throw std::runtime_error("EX-aero Error: species '" + species_names_[i] +
                    "' polarized-moment field length does not match its axes (fail fast)");
            }
            auto& p = h_species_params_[i];
            p.pmom_offset = static_cast<int>(h_curve_pool_.size());
            p.n_pmom_radius = nR;
            p.n_pmom_rh = nH;
            p.n_pmom_lambda = nL;
            p.n_pmom_pol = nP;
            p.n_pmom_moment = nM;
            h_curve_pool_.insert(h_curve_pool_.end(), c->rh.begin(), c->rh.end());
            h_curve_pool_.insert(h_curve_pool_.end(), c->lambda.begin(), c->lambda.end());
            h_curve_pool_.insert(h_curve_pool_.end(), pit->second.values.begin(), pit->second.values.end());
        }

        // Parse emissions mapping schemas if present (SPEC-EMISSIONS-002).
        // raw_name is diagnostics-only downstream (parity); the sequential mapping
        // index s is the raw CECE index.
        const int num_raw_species = reader.numEmissionsMappings();
        for (int s = 0; s < num_raw_species; ++s) {
            const GocartEmissionsMappingConfig em = reader.emissionsMapping(s);
            const int raw_cece_idx = s; // sequentially map indices

            for (const auto& mp : em.mappings) {
                int target_idx = getSpeciesIndex(mp.target_species);
                if (target_idx >= 0) {
                    auto& params = h_species_params_[target_idx];
                    params.emissions_mapping.is_active = true;
                    params.emissions_mapping.raw_cece_index = raw_cece_idx;
                    params.emissions_mapping.mass_split_fraction = mp.mass_split_fraction;

                    if (mp.is_modal_mode) {
                        params.emissions_mapping.is_modal_mode = true;
                        params.emissions_mapping.emitted_particle_diameter = mp.emitted_particle_diameter;
                        params.emissions_mapping.lognormal_sigma = mp.lognormal_sigma;
                    }
                }
            }
        }

        // Initialize solver state (allocates and uploads parameters + flat curve pool to GPU)
        if (solver_state_) {
            free_solver_state(solver_state_);
        }
        solver_state_ = create_solver_state(num_species_, h_species_params_.data(),
                                            h_curve_pool_.data(), static_cast<int>(h_curve_pool_.size()));
    }

    void GocartPackage::initialize(const std::string& config_yaml) {
        YamlReader reader{YAML::Load(config_yaml)};
        initializeImpl(reader);
    }

    void GocartPackage::initialize(const GocartConfig& config) {
        StructReader reader{config};
        initializeImpl(reader);
    }

    void GocartPackage::executeMicrophysics(
        const EnvironmentalStateView& env,
        View3D<double>& state,
        double delta_time_sec) {
        // GOCART is passive bulk aerosol. No internal microphysics evolution.
    }

    void GocartPackage::computeDerivedDiagnostics(
        const EnvironmentalStateView& env,
        const View3D<const double>& state,
        View3D<double>& diagnostics_out) {

        if (!solver_state_) {
            throw std::runtime_error("GocartPackage not initialized");
        }

        // --- Defensive Dimension Verification Checks (Hole 2) ---
        int num_cells = state.extent(0);
        int num_levels = state.extent(1);
        int num_species = state.extent(2);

        if (num_species != num_species_) {
            throw std::runtime_error("EX-aero State Error: State array species dimension (" + 
                                     std::to_string(num_species) + 
                                     ") does not match GOCART initialized species count (" + 
                                     std::to_string(num_species_) + ")!");
        }
        if (env.relative_humidity.extent(0) != state.extent(0) || 
            env.relative_humidity.extent(1) != state.extent(1)) {
            throw std::runtime_error("EX-aero Grid Error: Environmental Relative Humidity dimensions do not match the state array spatial grid!");
        }

        const double* rh_ptr = env.relative_humidity.data_handle();
        const double* thick_ptr = env.layer_thickness.data_handle();
        const double* state_ptr = state.data_handle();
        double* diags_ptr = diagnostics_out.data_handle();

        // Delegate to decoupled Kokkos solver
        run_gocart_diagnostics(
            solver_state_,
            num_cells, num_levels, num_species_,
            rh_ptr, thick_ptr, state_ptr, diags_ptr
        );
    }

    void GocartPackage::computeEmissions(
        const EnvironmentalStateView& env,
        const EmissionsInputView& emissions_in,
        View3D<double>& emissions_out) {

        if (!solver_state_) {
            throw std::runtime_error("GocartPackage not initialized");
        }

        int num_cells = emissions_in.flux.extent(0);
        int num_levels = emissions_in.flux.extent(1);
        int num_raw_species = emissions_in.flux.extent(2);
        int num_target_species = emissions_out.extent(2);

        if (num_target_species != num_species_) {
            throw std::runtime_error("EX-aero State Error: Target emissions array species dimension (" + 
                                     std::to_string(num_target_species) + 
                                     ") does not match GOCART initialized species count (" + 
                                     std::to_string(num_species_) + ")!");
        }

        if (env.layer_thickness.extent(0) != num_cells || env.layer_thickness.extent(1) != num_levels) {
            throw std::runtime_error("EX-aero Grid Error: Environmental Layer Thickness dimensions do not match the emissions array spatial grid!");
        }

        if (emissions_out.extent(0) != num_cells || emissions_out.extent(1) != num_levels) {
            throw std::runtime_error("EX-aero Grid Error: Output emissions array spatial dimensions do not match the input emissions grid!");
        }

        const double* thick_ptr = env.layer_thickness.data_handle();
        const double* raw_emissions_ptr = emissions_in.flux.data_handle();
        double* target_emissions_out_ptr = emissions_out.data_handle();
        int flux_type_code = (emissions_in.flux_type == FluxType::AREA_FLUX) ? 1 : 0;

        run_gocart_emissions(
            solver_state_,
            num_cells, num_levels, num_raw_species, num_target_species,
            flux_type_code, thick_ptr, raw_emissions_ptr, target_emissions_out_ptr
        );
    }

    void GocartPackage::computeOptics(
        const EnvironmentalStateView& env,
        const View3D<const double>& state,
        const View1D<const double>& wavelengths,
        View4D<double>& optics_out) {

        if (!solver_state_) {
            throw std::runtime_error("GocartPackage not initialized");
        }

        // --- Defensive Dimension Verification Checks (Hole 2) ---
        int num_cells = state.extent(0);
        int num_levels = state.extent(1);
        int num_species = state.extent(2);
        int num_bands = wavelengths.extent(0);

        if (num_species != num_species_) {
            throw std::runtime_error("EX-aero State Error: State array species dimension (" + 
                                     std::to_string(num_species) + 
                                     ") does not match GOCART initialized species count (" + 
                                     std::to_string(num_species_) + ")!");
        }
        if (env.relative_humidity.extent(0) != state.extent(0) || 
            env.relative_humidity.extent(1) != state.extent(1)) {
            throw std::runtime_error("EX-aero Grid Error: Environmental Relative Humidity dimensions do not match the state array spatial grid!");
        }

        const double* wavelengths_ptr = wavelengths.data_handle();
        const double* rh_ptr = env.relative_humidity.data_handle();
        const double* thick_ptr = env.layer_thickness.data_handle();
        const double* state_ptr = state.data_handle();
        double* optics_ptr = optics_out.data_handle();

        // Delegate to decoupled Kokkos solver
        run_gocart_optics(
            solver_state_,
            num_cells, num_levels, num_bands, num_species_,
            wavelengths_ptr, rh_ptr, thick_ptr, state_ptr, optics_ptr
        );
    }

    void GocartPackage::computeCCN(
        const EnvironmentalStateView& env,
        const View3D<const double>& state,
        const View1D<const double>& supersaturations,
        View4D<double>& ccn_out) {

        if (!solver_state_) {
            throw std::runtime_error("GocartPackage not initialized");
        }

        // --- Defensive Dimension Verification Checks (Hole 2) ---
        int num_cells = state.extent(0);
        int num_levels = state.extent(1);
        int num_species = state.extent(2);
        int num_ss = supersaturations.extent(0);

        if (num_species != num_species_) {
            throw std::runtime_error("EX-aero State Error: State array species dimension (" + 
                                     std::to_string(num_species) + 
                                     ") does not match GOCART initialized species count (" + 
                                     std::to_string(num_species_) + ")!");
        }
        if (env.relative_humidity.extent(0) != state.extent(0) || 
            env.relative_humidity.extent(1) != state.extent(1)) {
            throw std::runtime_error("EX-aero Grid Error: Environmental Relative Humidity dimensions do not match the state array spatial grid!");
        }

        const double* ss_ptr = supersaturations.data_handle();
        const double* temp_ptr = env.temperature.data_handle();
        const double* rh_ptr = env.relative_humidity.data_handle();
        const double* state_ptr = state.data_handle();
        double* ccn_ptr = ccn_out.data_handle();

        // Delegate to decoupled Kokkos solver
        run_gocart_ccn(
            solver_state_,
            num_cells, num_levels, num_ss, num_species_,
            ss_ptr, temp_ptr, rh_ptr, state_ptr, ccn_ptr
        );
    }

    // --- Dynamic Species-to-Index Queries (Hole 4) ---
    int GocartPackage::getSpeciesIndex(const std::string& name) const {
        for (int i = 0; i < num_species_; ++i) {
            if (species_names_[i] == name) {
                return i;
            }
        }
        return -1; // Not found
    }

    std::string GocartPackage::getSpeciesName(int index) const {
        if (index < 0 || index >= num_species_) {
            throw std::out_of_range("EX-aero Error: Species index (" + std::to_string(index) + ") out of bounds!");
        }
        return species_names_[index];
    }

    // --- GEOSmie MIE attribute surface (ADR-003 / research R9) ---
    void GocartPackage::computeAttributes(
        const EnvironmentalStateView& env,
        const View3D<const double>& state,
        int species_index,
        AttributeCategory category,
        const View1D<const double>& wavelengths,
        View3D<double>& attributes_out,
        View3D<int>* status_out) {

        if (species_index < 0 || species_index >= num_species_) {
            throw std::runtime_error("EX-aero Error: computeAttributes species index out of range");
        }
        const std::string name = species_names_[species_index];
        auto& store = MieTableStore::instance();
        const SpeciesCurve* c = store.find_curve(name);
        const int radius_node = c ? c->solver_radius_node : 0;

        // Attribute count is derived from the category, never a fixed literal (R10). For
        // polarized moments the count is the number of addressable (element, moment) slots
        // carried by THIS species' curve (data, FR-003), gated by activation (FR-010) so it
        // agrees with momentCounts() used to size the caller's output.
        int num_attr = 0;
        int n_pol = 0, n_moment = 0;
        switch (category) {
            case AttributeCategory::Microphysical:   num_attr = microphysical_indices::NUM_ATTRIBUTES; break;
            case AttributeCategory::SpectralOptical: num_attr = spectral_optical_indices::NUM_ATTRIBUTES; break;
            case AttributeCategory::PolarizedMoment:
                if (c && store.is_activated(name, category)) {
                    n_pol = c->n_pol();
                    n_moment = c->n_moment();
                }
                num_attr = n_pol * n_moment; // 0 => no moments / not activated (FR-008/FR-010)
                break;
            default: num_attr = 0;
        }

        const int num_cells = static_cast<int>(state.extent(0));
        const int num_levels = static_cast<int>(state.extent(1));
        // Multi-band spectral path (T027, FR-002): one column per (band, attribute) in
        // band-major order matching computeOptics memory layout; Microphysical ignores
        // the wavelength axis entirely (single block, NaN wavelength).
        int num_bands = (category == AttributeCategory::SpectralOptical)
                            ? static_cast<int>(wavelengths.extent(0)) : 1;
        if (num_bands < 1) num_bands = 1;
        if (static_cast<int>(attributes_out.extent(2)) < num_bands * num_attr) {
            throw std::runtime_error("EX-aero Error: computeAttributes output extent (" +
                std::to_string(attributes_out.extent(2)) + ") < num_bands * num_attributes (" +
                std::to_string(num_bands * num_attr) + ")");
        }
        if (status_out && static_cast<int>(status_out->extent(2)) < num_bands * num_attr) {
            throw std::runtime_error("EX-aero Error: computeAttributes status extent too small");
        }

        for (int ic = 0; ic < num_cells; ++ic) {
            for (int il = 0; il < num_levels; ++il) {
                const double rh = env.relative_humidity(ic, il);
                for (int ib = 0; ib < num_bands; ++ib) {
                    const double wavelength =
                        (category == AttributeCategory::SpectralOptical)
                            ? wavelengths(static_cast<std::size_t>(ib)) : std::nan("");
                    for (int ia = 0; ia < num_attr; ++ia) {
                        double value = 0.0;
                        // For moments, slot ia enumerates (element, moment) -> the encoded
                        // attribute_index (idx = moment*ELEMENT_STRIDE + element); the
                        // stride is the documented encoding, ia itself is dense (no gaps).
                        const int attr_index =
                            (category == AttributeCategory::PolarizedMoment)
                                ? (ia / n_pol) * polarized_moment_indices::ELEMENT_STRIDE
                                      + (ia % n_pol)
                                : ia;
                        const AttributeStatus st = store.query(
                            name, category, attr_index, rh, wavelength, radius_node, &value, nullptr);
                        const int slot = ib * num_attr + ia; // band-major (T027)
                        attributes_out(ic, il, slot) = value; // 0 only when status says unavailable
                        if (status_out) (*status_out)(ic, il, slot) = static_cast<int>(st);
                    }
                }
            }
        }
    }

    void GocartPackage::setSpeciesCurveConfig(const std::vector<SpeciesCurveConfig>& curves) {
        // Resolve through the store: resample radius onto config nodes + apply overrides
        // (precedence config > file > baked). A bad strategy/source aborts with FATAL ERROR.
        MieTableStore::instance().apply_curve_config(curves);
    }

    void GocartPackage::setAttributeActivation(const std::vector<std::string>& species,
                                               int categories_mask,
                                               const std::string& runtime_file_path) {
        auto& store = MieTableStore::instance();
        if (!runtime_file_path.empty()) {
            store.apply_runtime_file(runtime_file_path); // FATAL on bad file, no fallback (FR-009)
        }
        store.set_activation(species, categories_mask);
    }

    AttributeStatus GocartPackage::queryAttribute(int species_index, AttributeCategory category,
                                                  int attribute_index, double rh, double wavelength_m,
                                                  double* value_out,
                                                  ProvenanceInfo* provenance_out) const {
        if (species_index < 0 || species_index >= num_species_)
            return AttributeStatus::NotInSource;
        const std::string name = species_names_[species_index];
        auto& store = MieTableStore::instance();
        const SpeciesCurve* c = store.find_curve(name);
        const int radius_node = c ? c->solver_radius_node : 0; // hot-path default bin (R9)

        double local_value = 0.0;
        Provenance prov;
        const AttributeStatus status = store.query(
            name, category, attribute_index, rh, wavelength_m, radius_node, &local_value, &prov);

        // Write value only on an available/interpolated result (FR-008: never a silent 0).
        const bool ok = status == AttributeStatus::Available ||
                        status == AttributeStatus::AvailableFile ||
                        status == AttributeStatus::AvailableConfig ||
                        status == AttributeStatus::Interpolated;
        if (value_out && ok) *value_out = local_value;
        if (provenance_out) {
            prov.status = status;
            *provenance_out = prov.to_public(c ? c->n_radius() : 0, c ? c->n_rh() : 0,
                                             c ? c->n_lambda() : 0,
                                             c ? c->n_pol() : 0, c ? c->n_moment() : 0);
        }
        return status;
    }

    void GocartPackage::momentCounts(int species_index, int* num_pol_out,
                                     int* num_moment_out) const {
        int n_pol = 0, n_moment = 0;
        if (species_index >= 0 && species_index < num_species_) {
            const std::string name = species_names_[species_index];
            const SpeciesCurve* c = MieTableStore::instance().find_curve(name);
            // Counts are curve data gated by activation (FR-003/FR-008/FR-010); they must
            // agree with the slot count computeAttributes() fills for the same species.
            if (c && MieTableStore::instance().is_activated(
                         name, AttributeCategory::PolarizedMoment)) {
                n_pol = c->n_pol();
                n_moment = c->n_moment();
            }
        }
        if (num_pol_out) *num_pol_out = n_pol;
        if (num_moment_out) *num_moment_out = n_moment;
    }

} // namespace exaero
