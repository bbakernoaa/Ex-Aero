#pragma once
/// @file ModalChemistry.hpp
/// @brief Conservative modal source terms for host-owned chemical integration.
/// @author EX-aero contributors
/// @details History: 2026-10. Host layout_left buffers; no file I/O or
/// chemistry solver.
#include <exaero/ModalCoupling.hpp>

namespace exaero {
/// @brief Bind a transported vapor to one dry aerosol material.
/// Gas/equilibrium concentrations use kg of this dry material equivalent/m3.
/// Transport molecular weight is the physical vapor kg/mol. The host supplies
/// chemical stoichiometry and conversion from actual vapor mass to this basis.
struct ModalCondensable {
  std::string species;
  ModalVapor transport;
};
/// @brief Evaluate F=uptake*(gas-equilibrium), without advancing state.
/// @param gas (cell,level,nv), kg dry-material-equivalent/m3.
/// @param equilibrium (cell,level,nm,nv), same concentration basis.
/// @param aerosol_rates (cell,level,nf), kg/m3/s; number rows are zero.
/// @param gas_rates (cell,level,nv), equal/opposite exchange in the same basis.
/// @param jac Optional (cell,level,nf+nv,ni), ni>=nf+nv; aerosol fields precede
/// gases.
/// @param equilibrium_jac Optional (cell,level,nm*nv,ni), row=mode*nv+vapor.
/// Extra columns beyond nf+nv are host-defined chemical coordinates.
/// Omitting equilibrium derivatives explicitly freezes thermochemistry in jac.
/// Derivatives include dry-geometry feedback. Evaporation of absent material is
/// blocked, with zero active-branch derivative. MKPP owns positivity and empty-
/// mode number cleanup. Inputs/outputs must not overlap. All outputs are
/// atomic.
void modalPhaseExchange(IAerosolPackage &, const EnvironmentalStateView &,
                        const View3D<const double> &state,
                        std::span<const ModalCondensable>,
                        const View3D<const double> &gas,
                        const View4D<const double> &equilibrium,
                        View3D<double> &aerosol_rates,
                        View3D<double> &gas_rates,
                        View4D<double> *jac = nullptr,
                        const View4D<const double> *equilibrium_jac = nullptr);
/// @brief Conservative coating-driven aging RHS; never advances state.
/// @param condensation,coagulation (cell,level,nf), coating mass rates kg/m3/s.
/// Condensation is signed net phase exchange. Coagulation is nonnegative gross
/// incoming mass; number entries in both budgets must be zero. Sum budgets
/// once. The positive net coating volume rate divided by the population's exact
/// lognormal shell volume defines a loss frequency; all source mass AND number
/// transfer together. This generic shell closure is not calibrated MAM4 aging.
/// @param rates (cell,level,nf), signed mass/number transfer tendencies.
/// @param jac Optional (cell,level,nf,ni), ni>=nf, aerosol columns first.
/// @param condensation_jac,coagulation_jac Optional budget derivatives, nf..ni
/// rows and ni columns; a full phase-exchange Jacobian can be passed directly.
/// Missing budget derivatives freeze that budget. No hidden finite differences.
void modalAgingTerms(IAerosolPackage &, const EnvironmentalStateView &,
                     const View3D<const double> &state,
                     const View3D<const double> &condensation,
                     const View3D<const double> &coagulation,
                     View3D<double> &rates, View4D<double> *jac = nullptr,
                     const View4D<const double> *condensation_jac = nullptr,
                     const View4D<const double> *coagulation_jac = nullptr);
/// @brief Advance EX-aero microphysics and return accepted gross coagulation
/// gains (cell,level,nf), kg/m3 accumulated over dt; number entries are zero.
/// Overwrites the gain buffer, including dt=0. Divide by the chemical interval
/// to use as a frozen aging source. Transfer routes require host ownership so
/// gains remain attached to the post-coagulation population. Atomic outputs.
void advanceModalWithCoating(IAerosolPackage &, const EnvironmentalStateView &,
                             View3D<double> &state, double dt,
                             View3D<double> &coagulation_gain);
/// @brief Apply configured conservative growth-tail renaming after chemistry.
/// Requires transfer_owner: host. No coagulation or chemistry is repeated.
/// Call once per accepted host split interval, never on rejected MKPP stages.
void renameModalModes(IAerosolPackage &, const EnvironmentalStateView &,
                      View3D<double> &state);
} // namespace exaero
