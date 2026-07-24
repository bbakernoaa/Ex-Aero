#pragma once
#include <exaero/IAerosolPackage.hpp>
#include <gocart/GocartSpeciesParams.hpp>
#include <vector>

namespace exaero {

    struct GocartSolverState;

    class GocartPackage : public IAerosolPackage {
    private:
        int num_species_ = 0;
        std::vector<GocartSpeciesParams> h_species_params_;
        GocartSolverState* solver_state_ = nullptr;

    public:
        GocartPackage();
        ~GocartPackage() override;

        void initialize(const std::string& config_yaml) override;

        void executeMicrophysics(
            const EnvironmentalStateView& env,
            View3D<double>& state,
            double delta_time_sec) override;

        void computeDerivedDiagnostics(
            const EnvironmentalStateView& env,
            const View3D<const double>& state,
            View3D<double>& diagnostics_out) override;

        void computeOptics(
            const EnvironmentalStateView& env,
            const View3D<const double>& state,
            View3D<double>& optics_out) override;

        // Accessors for testing
        int get_num_species() const { return num_species_; }
        GocartSpeciesParams get_species_params(int i) const { return h_species_params_[i]; }
    };

} // namespace exaero
