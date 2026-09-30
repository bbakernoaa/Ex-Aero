#pragma once
/// @file DubovikSpheroidDatabase.hpp
/// @brief Trilinear lookup kernel over (n, k, x) spheroid optics.
///
/// Provides non-spherical (spheroid) extinction/scattering efficiencies
/// and Legendre phase moments as a function of complex refractive index
/// and Maxwell-Garnett size parameter @f$ x = \pi D / \lambda @f$,
/// following the Dubovik spheroid database methodology.
///
/// @warning Current status: the embedded grid is a physically-scaled
/// MOCK placeholder (see DubovikSpheroidDatabase.cpp). It is wired for
/// compilation and interface validation only; production runs use the
/// GEOSmie Mie tables in MieTableStore, not this database.
#include <vector>

namespace exaero {

/// @brief One kernel grid point: efficiencies plus 16 Legendre moments.
///
/// POD layout (no std::vector) so a point is trivially copyable during
/// the 8-corner trilinear blend.
struct SpheroidKernelPoint {
  double ext_efficiency; ///< Q_ext, dimensionless.
  double sca_efficiency; ///< Q_sca, dimensionless, <= Q_ext.
  double moments[16];    ///< Legendre phase-function moments m_0..m_15
                         ///< (m_0 = 1 normalization convention).
};

/// @brief Singleton spheroid optics database with trilinear lookup.
///
/// @warning Mock data — see the file-level warning above.
class SpheroidDatabase {
public:
  /// @brief Process-wide singleton (built once on first use).
  static SpheroidDatabase &instance();

  /// @brief Trilinear interpolation over (n_real, n_imag, size parameter).
  ///
  /// Inputs are clamped to the grid bounds first (HPC safety gate: never
  /// extrapolates, never reads out of bounds), then the surrounding 8
  /// corners are blended with linear weights per axis.
  /// @param n_real Real part of the complex refractive index n.
  /// @param n_imag Imaginary part k (absorption).
  /// @param size_parameter Maxwell-Garnett size parameter x = pi D / lambda.
  /// @return Interpolated efficiencies and moment set.
  SpheroidKernelPoint interpolate(double n_real, double n_imag,
                                  double size_parameter) const;

  /// @brief Ascending real-index grid (read-only accessor).
  const std::vector<double> &n_real_grid() const { return n_real_grid_; }
  /// @brief Ascending imaginary-index grid (read-only accessor).
  const std::vector<double> &n_imag_grid() const { return n_imag_grid_; }
  /// @brief Ascending size-parameter grid (read-only accessor).
  const std::vector<double> &size_param_grid() const {
    return size_param_grid_;
  }

private:
  /// @brief Builds the (mock) grid; private to enforce the singleton.
  SpheroidDatabase();
  ~SpheroidDatabase() = default;

  // Embedded static Dubovik spheroid database grid
  std::vector<double> n_real_grid_;     ///< n samples [dimensionless].
  std::vector<double> n_imag_grid_;     ///< k samples [dimensionless].
  std::vector<double> size_param_grid_; ///< x samples [dimensionless].
  /// @brief Flattened (n, k, x) kernel lattice, x fastest-varying
  /// (row-major ravel).
  std::vector<SpheroidKernelPoint> kernel_data_;
};

} // namespace exaero
