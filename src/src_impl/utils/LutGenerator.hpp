#pragma once
/// @file LutGenerator.hpp
/// @brief Offline lookup-table (LUT) builder over the (rh, band, species)
/// grid.
///
/// Allocates the flat storage for generated extinction/scattering
/// coefficients and phase moments and exposes them as column-major
/// mdspan views for downstream serialization. Sits on the
/// table-generation side of the pipeline (Python drives the final LUT
/// bake; see tools/generate_mie_tables.py), not the runtime hot path.
///
/// @warning Partial scaffolding: generate() is a compilation stub; the
/// production tables are baked by the Python toolchain from the pinned
/// GEOSmie snapshot.
#include <exaero/IAerosolPackage.hpp>
#include <string>
#include <vector>

namespace exaero {

/// @brief Builds and holds one generated aerosol optics LUT slab.
///
/// Owns the backing std::vector storage; the mdspan members are
/// non-owning views over it, which is why copy/move are deleted — a
/// copied object would dangle its views at the moved-from buffers.
class LutGenerator {
private:
  int num_rh_ = 0;       ///< RH bin count (LUT axis length).
  int num_bands_ = 0;    ///< Wavelength band count.
  int num_species_ = 0;  ///< Species count (single-species generator today).
  int num_moments_ = 16; ///< Legendre moment count (configurable via YAML).

  std::vector<double> rh_bins_;     ///< RH axis [fraction 0-1], ascending.
  std::vector<double> wavelengths_; ///< Band wavelengths [m], ascending.

  // Underlying allocations owned by LutGenerator
  std::vector<double> ext_data_;     ///< Flat extinction storage [m^-1] backing
                                     ///< ext_coeff_view_.
  std::vector<double> sca_data_;     ///< Flat scattering storage [m^-1] backing
                                     ///< sca_coeff_view_.
  std::vector<double> moments_data_; ///< Flat phase-moment storage backing
                                     ///< moments_view_.

  // Multi-dimensional output views (using layout_left)
  View3D<double> ext_coeff_view_; ///< Non-owning (rh, band, species) view.
  View3D<double> sca_coeff_view_; ///< Non-owning (rh, band, species) view.
  View4D<double> moments_view_;   ///< Non-owning (rh, band, species, moment)
                                  ///< view.

public:
  /// @brief Parse the generation grid from YAML and allocate the slabs.
  /// @param config_yaml YAML with a `generation_grid:` block (rh_bins,
  ///             wavelengths, optional legendre_moments).
  /// @throws std::runtime_error on a missing or malformed grid block.
  LutGenerator(const std::string &config_yaml);
  ~LutGenerator() = default;

  // Prevent dangerous dangling pointer aliasing due to mdspan view storage
  LutGenerator(const LutGenerator &) = delete;
  LutGenerator &operator=(const LutGenerator &) = delete;
  LutGenerator(LutGenerator &&) = delete;
  LutGenerator &operator=(LutGenerator &&) = delete;

  /// @brief Populate the LUT slabs (stub: not yet wired to a kernel).
  /// @param custom_spheroid_json_path Optional override optics source; empty
  ///             => built-in database path.
  void generate(const std::string &custom_spheroid_json_path = "");

  /// @brief Read-only (rh, band, species) extinction view.
  const View3D<double> &ext_coeff_view() const { return ext_coeff_view_; }
  /// @brief Read-only (rh, band, species) scattering view.
  const View3D<double> &sca_coeff_view() const { return sca_coeff_view_; }
  /// @brief Read-only (rh, band, species, moment) phase view.
  const View4D<double> &moments_view() const { return moments_view_; }

  /// @brief RH bin count.
  int num_rh() const { return num_rh_; }
  /// @brief Wavelength band count.
  int num_bands() const { return num_bands_; }
  /// @brief Species count.
  int num_species() const { return num_species_; }
  /// @brief Legendre moment count.
  int num_moments() const { return num_moments_; }

  /// @brief Serialize the grid metadata as a YAML fragment.
  std::string to_yaml() const;
};

} // namespace exaero
