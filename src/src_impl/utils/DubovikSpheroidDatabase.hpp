#pragma once
#include <vector>

namespace exaero {

    struct SpheroidKernelPoint {
        double ext_efficiency;
        double sca_efficiency;
        double moments[16];
    };

    class SpheroidDatabase {
    public:
        static SpheroidDatabase& instance();

        // Performs trilinear interpolation over refractive index and size parameter
        SpheroidKernelPoint interpolate(double n_real, double n_imag, double size_parameter) const;

        const std::vector<double>& n_real_grid() const { return n_real_grid_; }
        const std::vector<double>& n_imag_grid() const { return n_imag_grid_; }
        const std::vector<double>& size_param_grid() const { return size_param_grid_; }

    private:
        SpheroidDatabase();
        ~SpheroidDatabase() = default;

        // Embedded static Dubovik spheroid database grid
        std::vector<double> n_real_grid_;
        std::vector<double> n_imag_grid_;
        std::vector<double> size_param_grid_;
        std::vector<SpheroidKernelPoint> kernel_data_;
    };

} // namespace exaero
