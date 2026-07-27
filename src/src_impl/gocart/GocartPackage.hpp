#pragma once
#include <exaero/IAerosolPackage.hpp>
#include <gocart/GocartSpeciesParams.hpp>
#include <vector>
#include <string>

namespace exaero {

    struct GocartSolverState;

    class GocartPackage : public IAerosolPackage {
    private:
        int num_species_ = 0;
        std::vector<GocartSpeciesParams> h_species_params_;
        std::vector<std::string> species_names_; // Stores parsed YAML species names dynamically
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
            const View1D<const double>& wavelengths,
            View4D<double>& optics_out) override;

        void computeCCN(
            const EnvironmentalStateView& env,
            const View3D<const double>& state,
            const View1D<const double>& supersaturations,
            View4D<double>& ccn_out) override;

        // Dynamic Species-to-Index mapping query APIs
        int getSpeciesIndex(const std::string& name) const override;
        std::string getSpeciesName(int index) const override;

        // Accessors for testing
        int get_num_species() const { return num_species_; }
        GocartSpeciesParams get_species_params(int i) const { return h_species_params_[i]; }
    };

} // namespace exaero
