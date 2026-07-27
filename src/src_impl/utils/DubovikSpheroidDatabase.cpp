#include <utils/DubovikSpheroidDatabase.hpp>
#include <cmath>
#include <algorithm>

namespace exaero {

    SpheroidDatabase& SpheroidDatabase::instance() {
        static SpheroidDatabase db;
        return db;
    }

    SpheroidDatabase::SpheroidDatabase() {
        n_real_grid_ = { 1.50, 1.53, 1.56 };
        n_imag_grid_ = { 0.001, 0.003, 0.008 };
        size_param_grid_ = { 0.1, 1.0, 5.0, 10.0 };

        size_t total_points = n_real_grid_.size() * n_imag_grid_.size() * size_param_grid_.size();
        kernel_data_.resize(total_points);

        // Fill static grid with mock, physically scaled Dubovik spheroid coefficients
        for (size_t r = 0; r < n_real_grid_.size(); ++r) {
            for (size_t i = 0; i < n_imag_grid_.size(); ++i) {
                for (size_t s = 0; s < size_param_grid_.size(); ++s) {
                    double n = n_real_grid_[r];
                    double k = n_imag_grid_[i];
                    double x = size_param_grid_[s];

                    size_t idx = r * (n_imag_grid_.size() * size_param_grid_.size()) +
                                 i * size_param_grid_.size() + s;

                    // Mock physical values: ext scales with size parameter and refractive index
                    kernel_data_[idx].ext_efficiency = 2.0 * (1.0 - std::exp(-x * 0.5)) * (n / 1.53);
                    kernel_data_[idx].sca_efficiency = kernel_data_[idx].ext_efficiency * (1.0 - k * 10.0);
                    
                    // Set mock phase function Legendre moments (asymmetry decay)
                    for (int m = 0; m < 16; ++m) {
                        kernel_data_[idx].moments[m] = std::pow(0.9, m) / (1.0 + x * 0.1);
                    }
                }
            }
        }
    }

    SpheroidKernelPoint SpheroidDatabase::interpolate(double n_real, double n_imag, double size_parameter) const {
        // 1. Clamp input coordinates to grid bounds defensively (HPC safety gate)
        double r_val = std::clamp(n_real, n_real_grid_.front(), n_real_grid_.back());
        double i_val = std::clamp(n_imag, n_imag_grid_.front(), n_imag_grid_.back());
        double s_val = std::clamp(size_parameter, size_param_grid_.front(), size_param_grid_.back());

        // Lambda helper to find lower bound coordinate grid indices and linear weights
        auto find_grid_indices = [](const std::vector<double>& grid, double val, int& lower_idx, double& weight) {
            if (val <= grid.front()) {
                lower_idx = 0;
                weight = 0.0;
                return;
            }
            if (val >= grid.back()) {
                lower_idx = static_cast<int>(grid.size() - 2);
                weight = 1.0;
                return;
            }
            for (size_t i = 0; i < grid.size() - 1; ++i) {
                if (val >= grid[i] && val <= grid[i+1]) {
                    lower_idx = static_cast<int>(i);
                    weight = (val - grid[i]) / (grid[i+1] - grid[i]);
                    return;
                }
            }
        };

        int r0 = 0; double wr = 0.0; find_grid_indices(n_real_grid_, r_val, r0, wr);
        int i0 = 0; double wi = 0.0; find_grid_indices(n_imag_grid_, i_val, i0, wi);
        int s0 = 0; double ws = 0.0; find_grid_indices(size_param_grid_, s_val, s0, ws);

        int r1 = r0 + 1;
        int i1 = i0 + 1;
        int s1 = s0 + 1;

        // Retrieve the 8 surrounding grid point references once (HPC optimization!)
        auto get_point = [&](int r, int i, int s) -> const SpheroidKernelPoint& {
            size_t idx = r * (n_imag_grid_.size() * size_param_grid_.size()) +
                         i * size_param_grid_.size() + s;
            return kernel_data_[idx];
        };

        const auto& p000 = get_point(r0, i0, s0);
        const auto& p100 = get_point(r1, i0, s0);
        const auto& p010 = get_point(r0, i1, s0);
        const auto& p110 = get_point(r1, i1, s0);
        const auto& p001 = get_point(r0, i0, s1);
        const auto& p101 = get_point(r1, i0, s1);
        const auto& p011 = get_point(r0, i1, s1);
        const auto& p111 = get_point(r1, i1, s1);

        auto lerp = [](double v0, double v1, double w) { return v0 + w * (v1 - v0); };

        SpheroidKernelPoint out;

        // Trilinear Interpolation of Extinction Efficiency
        double ext_00 = lerp(p000.ext_efficiency, p100.ext_efficiency, wr);
        double ext_10 = lerp(p010.ext_efficiency, p110.ext_efficiency, wr);
        double ext_01 = lerp(p001.ext_efficiency, p101.ext_efficiency, wr);
        double ext_11 = lerp(p011.ext_efficiency, p111.ext_efficiency, wr);

        double ext_0 = lerp(ext_00, ext_10, wi);
        double ext_1 = lerp(ext_01, ext_11, wi);

        out.ext_efficiency = lerp(ext_0, ext_1, ws);

        // Trilinear Interpolation of Scattering Efficiency
        double sca_00 = lerp(p000.sca_efficiency, p100.sca_efficiency, wr);
        double sca_10 = lerp(p010.sca_efficiency, p110.sca_efficiency, wr);
        double sca_01 = lerp(p001.sca_efficiency, p101.sca_efficiency, wr);
        double sca_11 = lerp(p011.sca_efficiency, p111.sca_efficiency, wr);

        double sca_0 = lerp(sca_00, sca_10, wi);
        double sca_1 = lerp(sca_01, sca_11, wi);

        out.sca_efficiency = lerp(sca_0, sca_1, ws);

        // Trilinear Interpolation of Legendre Moments Array
        for (int m = 0; m < 16; ++m) {
            double m_00 = lerp(p000.moments[m], p100.moments[m], wr);
            double m_10 = lerp(p010.moments[m], p110.moments[m], wr);
            double m_01 = lerp(p001.moments[m], p101.moments[m], wr);
            double m_11 = lerp(p011.moments[m], p111.moments[m], wr);

            double m_0 = lerp(m_00, m_10, wi);
            double m_1 = lerp(m_01, m_11, wi);

            out.moments[m] = lerp(m_0, m_1, ws);
        }

        return out;
    }

} // namespace exaero
