#pragma once
/// @file ModalCoupling.hpp
/// @brief Additive modal/MKPP interface, independent of IAerosolPackage's ABI.
/// @author EX-aero contributors
/// @details History: 2026-10. Usage: dynamic_cast<IModalCoupling*>(package).
/// Input/output arrays are host memory, layout_left, unless documented
/// otherwise.
#include <exaero/IAerosolPackage.hpp>
#include <exaero/ModalConfig.hpp>
#include <exaero/RadiationBands.hpp>
#include <memory>
#include <span>

namespace exaero {
/// @brief Property indices in (cell,level,mode,property) host output.
namespace modal_properties {
inline constexpr int DRY_RADIUS = 0, WET_RADIUS = 1, NUMBER = 2, MASS = 3,
                     VOLUME = 4, KAPPA = 5, SURFACE_AREA = 6, WET_DENSITY = 7,
                     WATER = 8, COUNT = 9;
}
/// @brief Geometry-only vapor input; MKPP supplies thermochemistry separately.
struct ModalVapor {
  double diffusivity = 0;      ///< Molecular diffusivity [m2/s].
  double molecular_weight = 0; ///< kg/mol.
  double accommodation = 1;    ///< (0,1].
};
/// @brief Fully integrated mode optics per volume of air, including water once.
struct ModalOpticalResult {
  double extinction = 0, scattering = 0; ///< m^-1.
  double asymmetry = 0;                  ///< Scattering-weighted g in [-1,1].
  double backscatter = 0, lidar_backscatter = 0;   ///< m^-1 sr^-1.
  double fine_extinction = 0, fine_scattering = 0; ///< Dry radius < 0.5 um.
  double pm25_extinction = 0,
         pm25_scattering = 0; ///< Wet aero diameter <2.5 um.
};
/// @brief Synchronous host provider for validated internally mixed optics.
/// Implementations must validate size/composition/RH/wavelength coverage and
/// throw on incompatibility. Return complete mode integrals, not coefficients
/// to be multiplied by mass again. Called only for populated modes.
class IModalOpticsProvider {
public:
  virtual ~IModalOpticsProvider() = default;
  /// @brief Evaluate one mode; properties use modal_properties indices.
  virtual ModalOpticalResult
  evaluate(const ModalConfig &, std::size_t mode, double temperature, double rh,
           double wavelength, std::span<const double> species_mass,
           std::span<const double> properties) const = 0;
};
/// @brief Synchronous provider for internally mixed native band optics.
/// Results are complete mode integrals, already including water exactly once.
class IModalBandOpticsProvider {
public:
  virtual ~IModalBandOpticsProvider() = default;
  virtual bool supports(const ModalConfig &, std::size_t mode,
                        const RadiationBandInfo &) const = 0;
  virtual ModalOpticalResult evaluateBand(
      const ModalConfig &, std::size_t mode, const RadiationBandInfo &,
      double temperature, double rh, std::span<const double> species_mass,
      std::span<const double> properties) const = 0;
};
/// @brief Rate contributions for MKPP; EX-aero never applies these to state.
struct ModalNucleationTerms {
  int destination = -1;
  double number_rate = 0;    ///< particles/m3/s.
  double precursor_loss = 0; ///< kg H2SO4/m3/s; valid for pure H2SO4 seeds.
  std::vector<double> mass_rates; ///< kg/m3/s, configuration species order.
};
/// @brief Provenance and declared semantics of one resolved component table.
struct ModalOpticalTableInfo {
  std::string mode, species, source_label, source_version, citation;
  ModalTableKind kind = ModalTableKind::Unknown;
};
/// @brief Separate extension interface; existing package vtables are unchanged.
class IModalCoupling {
public:
  virtual ~IModalCoupling() = default;
  /// @brief Initialize with the same validation as YAML.
  virtual void initializeModal(const ModalConfig &) = 0;
  /// @brief Read immutable validated configuration.
  virtual const ModalConfig &modalConfig() const = 0;
  /// @brief Source identifiers/citations for immutable per-package snapshots.
  virtual std::vector<ModalOpticalTableInfo> opticalTables() const = 0;
  /// @brief Query properties and optional d(property)/d(state_field).
  /// Derivatives shape: (cell,level,mode,property*field_count+field).
  /// Native MAM4xx wet-property Jacobians are unavailable and fail atomically.
  /// Dry uptake Jacobians remain supported.
  /// Dry radius derivatives at an empty mode are defined as zero; its reference
  /// radius is returned. Wet-radius/density derivatives hold RH fixed.
  virtual void modalProperties(const EnvironmentalStateView &,
                               const View3D<const double> &, View4D<double> &,
                               View4D<double> *derivatives = nullptr) const = 0;
  /// @brief Fuchs-Sutugin uptake coefficients [s^-1], (cell,level,mode,vapor).
  /// Dry lognormal geometry is used; equilibrium concentrations belong to MKPP.
  /// Optional derivatives: (cell,level,mode,vapor*field_count+field).
  virtual void modalUptake(const EnvironmentalStateView &,
                           const View3D<const double> &,
                           std::span<const ModalVapor>, View4D<double> &,
                           View4D<double> *derivatives = nullptr) const = 0;
  /// @brief Evaluate J=A*C^p and balanced seed yields from kg H2SO4/m3.
  virtual ModalNucleationTerms nucleationTerms(double sulfuric_acid) const = 0;
  /// @brief Infer emitted particle number from species masses and mode sizing.
  virtual double emittedNumber(std::size_t mode,
                               std::span<const double> masses) const = 0;
  /// @brief Install/remove the optional internally mixed provider.
  virtual void
      setModalOpticsProvider(std::shared_ptr<const IModalOpticsProvider>) = 0;
  virtual void setModalBandOpticsProvider(
      std::shared_ptr<const IModalBandOpticsProvider>) = 0;
};
/// @brief Checked access to the separate modal interface (throws if
/// unsupported).
inline IModalCoupling &modalCoupling(IAerosolPackage &package) {
  auto *modal = dynamic_cast<IModalCoupling *>(&package);
  if (!modal)
    throw std::logic_error(
        "FATAL ERROR: package has no generic modal coupling");
  return *modal;
}
/// @brief Structured initialization without altering the base class vtable.
inline void initializeModal(IAerosolPackage &package,
                            const ModalConfig &config) {
  modalCoupling(package).initializeModal(config);
}
} // namespace exaero
