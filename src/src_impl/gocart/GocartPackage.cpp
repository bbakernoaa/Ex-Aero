#include <gocart/GocartPackage.hpp>
#include <yaml-cpp/yaml.h>
#include <stdexcept>

namespace exaero {

    // Decoupled solver helper function declarations (implemented in GocartSolver.cpp)
    GocartSolverState* create_solver_state(int num_species, const GocartSpeciesParams* params);
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

    void GocartPackage::initialize(const std::string& config_yaml) {
        YAML::Node config = YAML::Load(config_yaml);
        if (!config["species"]) {
            throw std::runtime_error("GOCART YAML config missing 'species' field");
        }

        auto species_node = config["species"];
        num_species_ = species_node.size();

        h_species_params_.clear();
        h_species_params_.reserve(num_species_);
        species_names_.clear();
        species_names_.reserve(num_species_);

        for (int i = 0; i < num_species_; ++i) {
            auto s = species_node[i];
            bool has_lookup = s["has_optics_lookup"] && s["has_optics_lookup"].as<bool>();
            
            // Extract the dynamic species name
            std::string name = s["name"] ? s["name"].as<std::string>() : ("Species_" + std::to_string(i));
            species_names_.push_back(name);

            GocartSpeciesParams p{
                s["dry_density"].as<double>(),
                s["molecular_weight"].as<double>(),
                s["dry_particle_diameter"].as<double>(),
                s["hygroscopicity"].as<double>(),
                s["lognormal_sigma"].as<double>(),
                s["lognormal_dg"].as<double>(),
                s["refractive_index_real"].as<double>(),
                s["refractive_index_imag"].as<double>(),
                has_lookup
            };

            if (has_lookup) {
                auto rh_bins = s["rh_bins"];
                auto ext_lookup = s["ext_lookup"];
                auto ssa_lookup = s["ssa_lookup"];
                auto asm_lookup = s["asm_lookup"];
                for (int j = 0; j < 8; ++j) {
                    p.rh_bins[j] = rh_bins[j].as<double>();
                    p.ext_lookup[j] = ext_lookup[j].as<double>();
                    p.ssa_lookup[j] = ssa_lookup[j].as<double>();
                    p.asm_lookup[j] = asm_lookup[j].as<double>();
                }
            }
            h_species_params_.push_back(p);
        }

        // Parse emissions mapping schemas if present (SPEC-EMISSIONS-002)
        if (config["emissions_mapping"]) {
            auto mapping_node = config["emissions_mapping"];
            for (size_t s = 0; s < mapping_node.size(); ++s) {
                std::string raw_name = mapping_node[s]["raw_name"] ? mapping_node[s]["raw_name"].as<std::string>() : "CECE_Raw_Species";
                int raw_cece_idx = static_cast<int>(s); // sequentially map indices
                
                auto mappings = mapping_node[s]["mappings"];
                if (mappings) {
                    for (size_t m = 0; m < mappings.size(); ++m) {
                        std::string target_spec = mappings[m]["target_species"].as<std::string>();
                        double frac = mappings[m]["mass_split_fraction"].as<double>();
                        
                        int target_idx = getSpeciesIndex(target_spec);
                        if (target_idx >= 0) {
                            auto& params = h_species_params_[target_idx];
                            params.emissions_mapping.is_active = true;
                            params.emissions_mapping.raw_cece_index = raw_cece_idx;
                            params.emissions_mapping.mass_split_fraction = frac;
                            
                            if (mappings[m]["is_modal_mode"] && mappings[m]["is_modal_mode"].as<bool>()) {
                                params.emissions_mapping.is_modal_mode = true;
                                params.emissions_mapping.emitted_particle_diameter = mappings[m]["emitted_particle_diameter"].as<double>();
                                params.emissions_mapping.lognormal_sigma = mappings[m]["lognormal_sigma"].as<double>();
                            }
                        }
                    }
                }
            }
        }

        // Initialize solver state (allocates and uploads parameters to GPU)
        if (solver_state_) {
            free_solver_state(solver_state_);
        }
        solver_state_ = create_solver_state(num_species_, h_species_params_.data());
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

} // namespace exaero
