#pragma once
#include <exaero/IAerosolPackage.hpp>
#include <string>
#include <vector>

namespace exaero {

    class LutGenerator {
    private:
        int num_rh_ = 0;
        int num_bands_ = 0;
        int num_species_ = 0;
        int num_moments_ = 16;

        std::vector<double> rh_bins_;
        std::vector<double> wavelengths_;
        
        // Underlying allocations owned by LutGenerator
        std::vector<double> ext_data_;
        std::vector<double> sca_data_;
        std::vector<double> moments_data_;

        // Multi-dimensional output views (using layout_left)
        View3D<double> ext_coeff_view_; // (rh, band, species)
        View3D<double> sca_coeff_view_; // (rh, band, species)
        View4D<double> moments_view_;   // (rh, band, species, moment)

    public:
        LutGenerator(const std::string& config_yaml);
        ~LutGenerator() = default;

        // Prevent dangerous dangling pointer aliasing due to mdspan view storage
        LutGenerator(const LutGenerator&) = delete;
        LutGenerator& operator=(const LutGenerator&) = delete;
        LutGenerator(LutGenerator&&) = delete;
        LutGenerator& operator=(LutGenerator&&) = delete;

        void generate(const std::string& custom_spheroid_json_path = "");

        const View3D<double>& ext_coeff_view() const { return ext_coeff_view_; }
        const View3D<double>& sca_coeff_view() const { return sca_coeff_view_; }
        const View4D<double>& moments_view() const { return moments_view_; }

        int num_rh() const { return num_rh_; }
        int num_bands() const { return num_bands_; }
        int num_species() const { return num_species_; }
        int num_moments() const { return num_moments_; }
        
        std::string to_yaml() const;
    };

} // namespace exaero
