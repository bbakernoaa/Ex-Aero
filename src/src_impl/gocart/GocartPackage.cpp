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

        for (int i = 0; i < num_species_; ++i) {
            auto s = species_node[i];
            bool has_lookup = s["has_optics_lookup"] && s["has_optics_lookup"].as<bool>();
            
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

        // Extract raw data pointers and sizes from the dynamic C++20 mdspan views
        int num_cells = state.extent(0);
        int num_levels = state.extent(1);

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

    void GocartPackage::computeOptics(
        const EnvironmentalStateView& env,
        const View3D<const double>& state,
        const View1D<const double>& wavelengths,
        View4D<double>& optics_out) {

        if (!solver_state_) {
            throw std::runtime_error("GocartPackage not initialized");
        }

        int num_cells = state.extent(0);
        int num_levels = state.extent(1);
        int num_bands = wavelengths.extent(0);

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

} // namespace exaero
