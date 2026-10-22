# Modal exchange, coating aging, and MKPP integration

`ModalChemistry.hpp` supplies additive C++ interfaces for MKPP/host-owned
integration. Equivalent C entry points are in `ModalChemistry_C.h` and the
Fortran `exaero_interface` module. These functions do not run a chemical solver,
choose thermodynamics, or modify an external MKPP repository/generated mechanism.
The optional chemistry driver in this repository still requires external registry
inputs. A mechanism's RHS/Jacobian callback must explicitly call these interfaces.

## Phase exchange and mass basis

Bind each transported vapor to one distinct configured dry aerosol material:

```cpp
std::array<exaero::ModalCondensable, 1> bindings{{
    {"organic", {8e-6, 0.15, 1.0}} // diffusivity m2/s, actual vapor kg/mol, accommodation
}};
```

For mode m and vapor v, the exchange flux is

```
F(m,v) = uptake(m,v) * (gas(v) - equilibrium(m,v))
aerosol_rate(m,material(v)) += F(m,v)
gas_rate(v) -= F(m,v)
```

Uptake is the existing dry-lognormal Fuchs-Sutugin integral, in s^-1, including
prognostic number and geometry. It is recomputed at the current chemical stage.
Positive flux condenses; negative flux evaporates. A material with zero mass
cannot evaporate; the Jacobian of that blocked branch is zero. Empty modes have
zero flux and cannot be populated by condensation without nucleation/seeding.
Disallowed material slots do not exchange. Phase exchange leaves number rows zero.
MKPP must enforce nonnegative states and remove residual number if a mode fully
evaporates; there is no timestep-dependent clipping inside the RHS.

**Gas and equilibrium concentrations use kg of target dry-material equivalent
per m3 of air**, matching the bound aerosol field. For the same compound this
is ordinary mass concentration. If H2SO4 gas is bound to SO4 aerosol, convert
to SO4-equivalent mass first; the transport molecular weight is still physical
H2SO4. Balance of these returned rates is in the declared equivalent-material
basis. Water or other coproduct stoichiometry belongs to MKPP. This interface
does not claim conservation of total physical mass across omitted coproducts.
The host must never interchange complete-salt mass, ionic mass and vapor mass.

Equilibrium concentrations (including activity, volatility and Kelvin corrections
chosen by the host) are supplied by MKPP/thermodynamics. The geometry kernel does
not add another equilibrium or humidity correction.

## Dense Jacobian contract

All arrays are host column-major. Let nf be the configured aerosol field count,
nm the mode count, and nv the number of vapor bindings. Choose ni >= nf+nv
Jacobian columns; extra columns represent explicitly chosen host coordinates. Local chemical state is
ordered **all nf aerosol fields, then nv vapor concentrations**.

| Buffer | Shape | Meaning |
| --- | --- | --- |
| gas | (cell,level,nv) | Equivalent-material concentration |
| equilibrium | (cell,level,nm,nv) | Host equilibrium concentration |
| aerosol_rates | (cell,level,nf) | Signed mass/number tendencies |
| gas_rates | (cell,level,nv) | Signed gas tendencies |
| jac | (cell,level,nf+nv,ni) | Combined exchange RHS derivatives |
| equilibrium_jac | (cell,level,nm*nv,ni) | Host thermodynamic derivatives; row=m*nv+v |

`modalPhaseExchange` includes analytic uptake-geometry derivatives and any supplied
equilibrium derivatives. Omitting `equilibrium_jac` explicitly freezes equilibrium
in the returned Jacobian; it is not the complete thermodynamic Jacobian unless
equilibrium is independent of state. Environmental and transport inputs are held
fixed. Extra columns have the units chosen for the corresponding supplied
equilibrium derivatives; the host must map those columns into its mechanism
with the matching scale. They can include noncondensing chemical species. Every material's gas/particle RHS and Jacobian budgets cancel.

`modalChemicalMechanismMap` extends the existing aerosol mapping with ordered
vapor bindings, target-material molecular weights and an explicit
vapor-molecules/material-unit ratio. `modalChemicalJacobianScale` supplies the
row/column chain-rule factor.

For a mechanism variable z_i = scale_i*x_i, convert a returned rate by scale_i
and a Jacobian entry by scale_i/scale_j. Include gas equivalent-material factors
as well as the existing aerosol/particle-number mappings. MKPP concentration
conventions and stoichiometry must be explicit. Evaluate geometry, equilibrium
and their derivatives at the same stage state.

## Coating-driven aging

Configure routes as in [the example](../../examples/modal_chemistry.yaml):

```yaml
aging_routes:
  - source: primary_carbon
    destination: accumulation
    coating_species: [SO4, SOA]
    coating_thickness: 1e-9
```

The thickness is a host-calibrated dry coating shell thickness, not a supplied
MAM4 aging calibration. `modalAgingTerms` accepts two distinct mass-rate budgets:

- Signed net condensation/evaporation rates for every packed aerosol field.
- Nonnegative **gross incoming coagulation** mass rates, not net mode changes.

Both budgets have shape (cell,level,nf), and number entries must be zero. For the
selected coating species, combine their signed volume rates once and use the
positive result B. The exact lognormal population shell volume is

```
S = 4*pi*N*(rg^2*exp(2*ln(sigma)^2)*h
          +rg*exp(0.5*ln(sigma)^2)*h^2 + h^3/3)
aging_frequency = max(B,0)/S
```

The resulting transfer rate is frequency times **each source material mass and
particle number**. Destinations receive equal amounts. Routes use a common stage
snapshot, so incoming particles are not aged a second time through another route
in the same RHS evaluation. Duplicate sources, cycles and incompatible destination
species are rejected. Nonpositive coating budgets and empty modes have zero
aging rates, with zero active-branch Jacobian there.

This is a configurable, continuous shell-coating closure. It is **not** the
native finite-step MAM4 aging fraction or a calibrated replacement. Stiffness and
positivity are handled by MKPP. Supply both budget derivatives for a complete
aging Jacobian. The full exchange Jacobian can be passed directly as the
condensation derivative input: leading aerosol rows are used. Missing budget
derivatives explicitly freeze that budget. Native-water wet-property Jacobian
limitations do not affect these dry-geometry derivatives.

## Gross coagulation gains and process ownership

`advanceModalWithCoating` advances the configured coagulation and overwrites a
caller-owned gain buffer with accepted incoming mass **integrated over dt**.
Generic Brownian substeps accumulate gross cross-mode gains. Native MAM4xx
exports its primary-carbon coating transfer and accumulation gains. Self
collisions add no coating; number gain entries are zero. dt=0 returns zero gains.
No aging transfer is silently applied inside this step.

Set `transfer_owner: host` when growth transfer is configured. This keeps coating
gains attached to the post-coagulation populations. A typical split interval is:

1. Advance EX-aero coagulation with the gain output.
2. Divide gains by the chemical interval to form a frozen coagulation-coating
   rate. This uniform-in-time use of an accepted split-step budget is an explicit
   approximation; a fully coupled coagulation RHS needs a different host schedule.
3. Let MKPP evaluate exchange and aging together at every stage, supplying the
   current exchange rate/Jacobian and the frozen coagulation budget. Add aging
   aerosol rows to the exchange/chemistry rows; gas rows are unchanged by aging.
4. After MKPP accepts the interval, call `renameModalModes` once.

Do not accumulate gain buffers across rejected steps, replay them in multiple
accepted intervals, or also apply a separate aging transfer. The host owns their
lifetime and rollback. These routines are stateless; repeated RHS evaluations
read the same budget and do not consume or mutate it.

`renameModalModes` applies the existing conservative lognormal growth-tail
reclosure using `transfer_routes`; it repeats no coagulation. It requires host
ownership. Legacy `executeMicrophysics` retains automatic growth transfer unless
`transfer_owner: host` is selected. Renaming is a post-acceptance redistribution,
not an RHS contribution or native MAM4 growth-renaming implementation.

## Validation and limits

Inputs and outputs must not overlap. Dense shapes, finite values, rate conventions
and material bindings are validated before publishing results. Outputs are atomic
on failure. New APIs preserve the existing IAerosolPackage and IModalCoupling
vtables; package chemistry interfaces are additive free functions.

Tests cover independent shell-volume references, material/number conservation,
geometry and thermodynamic Jacobian finite differences, evaporation boundaries,
gross coagulation gains, transfer ownership, native MAM4xx budgets, and C/Fortran
column-major interoperability. Run the normal CTest suite with MAM4xx ON and OFF
using Kokkos 5.1.1. These checks validate the coupling and generic closure; native
MAM4 aging parity, generated-mechanism integration, thermodynamic calibration,
and GPU validation remain separate work.
