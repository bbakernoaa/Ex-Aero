# Modal exchange, coating aging, and MKPP integration

`ModalChemistry.hpp` supplies additive C++ interfaces for MKPP/host-owned
integration. Equivalent C entry points are in `ModalChemistry_C.h` and the
Fortran `exaero_interface` module. These functions do not run a chemical solver,
choose thermodynamics, or modify an external MKPP repository/generated mechanism.
The optional chemistry driver in this repository still requires external registry
inputs. A mechanism using these interfaces must explicitly call them in its
RHS/Jacobian callback. The symbolic YAML example below evaluates its own
illustrative closures instead and must not add duplicate callback contributions.

## MKPP reaction YAML examples

[modal_processes.mkpp.yaml](../../examples/modal_processes.mkpp.yaml) is a complete
OpenAtmos mechanism illustrating exchange, nucleation, and aging. It is separate
from [modal_chemistry.yaml](../../examples/modal_chemistry.yaml), which configures
EX-aero. The mechanism uses **frozen, illustrative coefficients**, not atmospheric
calibrations. It declares only the fields needed for these examples; a complete
four-mode host callback mechanism must bind every configured aerosol field, including unused
fields, through `modalChemicalMechanismMap`.

The example was checked against
[MKPP fb2a87b5e6ad](https://github.com/bbakernoaa/MKPP/tree/fb2a87b5e6adc46aa07bfb5a5130633cd3ba76c5).
That revision reads `A` at reaction level and multiplies it by reactant
concentrations raised to their stoichiometric coefficients. Its documentation's
nested `parameters.rate_expression` form is not implemented by its lowering
code. No `MODAL_AGING`, `NUCLEATION`, or EX-aero callback reaction type is wired
into that compiler. These examples therefore use supported `ARRHENIUS` entries
with B=C=0. Do not replace a numeric `A` with a callback name and assume the
generated solver will invoke EX-aero.

### State-dependent YAML with the MKPP extension

[modal_processes_dynamic.mkpp.yaml](../../examples/modal_processes_dynamic.mkpp.yaml)
uses the [companion MKPP extension](https://github.com/bbakernoaa/MKPP/pull/2)
for actual stage-dependent expressions,
`rate_definitions`, `rate_basis: flux`, and `reaction_orders`. The corresponding
[EX-aero material configuration](../../examples/modal_processes_dynamic.aerosol.yaml)
uses complete H2SO4, organic material, and black carbon in three modes. It uses
complete acid mass, **not** the MAM4 preset's SO4-equivalent material convention.

For each mode `m` in `aitken`, `primary_carbon`, and `accumulation`, bind:

| EX-aero field | MKPP variable | SI-to-mechanism scale |
| --- | --- | --- |
| `sulfuric_acid_m_i` | `ACID_m` | 6.02214076e17/0.098079 |
| `organic_m_i` | `ORGANIC_m` | 6.02214076e17/0.15 |
| `black_carbon_m_i` | `BC_m` | 6.02214076e17/0.012011 |
| `num_m_i` | `N_m` | 1e-6 |

Bind gases `H2SO4` to `sulfuric_acid` and `ORGANIC_VAPOR` to `organic`, with one
vapor molecule per material unit. Read actual field names from `stateFields()`;
`m` in this table is a mode-name placeholder, not a literal field name.
`COAG_VOLUME_RATE` is a separate host coordinate, not a condensable binding.

For example, the YAML represents empirical nucleation directly as:

```yaml
- type: CUSTOM
  reactants: {H2SO4: 100}
  products: {ACID_aitken: 100, N_aitken: 1}
  parameters:
    rate_expression: 1.0e-8
    reaction_orders: {H2SO4: 1.5}
  stiff: true
```

The flux is J=1e-8*C_H2SO4^1.5, independent of the seed yield of 100 acid
molecules/particle. Concentrations use formula units/cm3 and particles/cm3.
An EX-aero empirical SI prefactor A_SI becomes A_cm=A_SI*10^(6*(p-1)); the
companion configuration uses p=1.5 and A_SI=1e-11. The seed radius is derived from
100 acid molecules at the configured density, including the receiving mode's
lognormal volume moment exp(4.5*ln(sigma)^2). No separate number source is needed.

The dynamic file also contains signed exchange fluxes and shared aging-frequency
definitions with empty-mode and positive-budget guards. Named expressions expand
symbolically, so MKPP differentiates their geometry, thermodynamics, and coating
dependence. The auxiliary fixed-role `COAG_VOLUME_RATE` has units m3 coating/m3
air/s; it is **not** a chemical concentration or an aerosol field. The host
populates it from selected gross gains/density/accepted interval and excludes it
from chemical material/number budgets.

Its illustrative exchange closure is the exact lognormal **continuum-limit**
uptake 4*pi*D*N*E[r], with a linear organic equilibrium concentration. It does
not reproduce the production transition-regime uptake integral returned by
`modalPhaseExchange`. The shell-aging moment follows the closure below, with
matching material masses and width. Use these laws as readable, compiler-checked
examples; thermodynamic/transport and MAM4 calibration remain scientific choices.
Do not evaluate the EX-aero exchange/aging callbacks on top of these YAML terms.
The callbacks remain the path to the existing EX-aero production closures until
the host or compiler explicitly connects them.

### Condensation and evaporation

```yaml
reactions:
  - type: ARRHENIUS
    reactants: {SOA_VAP_EQ: 1.0}
    products: {SOA_primary_carbon: 1.0}
    A: 1.0e-3
    stiff: true
  - type: ARRHENIUS
    reactants: {SOA_primary_carbon: 1.0}
    products: {SOA_VAP_EQ: 1.0}
    A: 2.0e-4
    stiff: true
```

This pair gives F=1e-3*Cgas-2e-4*Cparticle in material units/cm3/s, equivalent to
frozen uptake=1e-3/s and illustrative Ceq=0.2*Cparticle. Number is unchanged.
For the production closure, **replace both contributions** with
`modalPhaseExchange` RHS/Jacobian contributions at each stage. Supply current
geometry, equilibrium concentrations, and their derivatives. Do not add the
callback to these demonstration reactions or exchange will be counted twice.

### Empirical nucleation and number yield

```yaml
reactions:
  - type: ARRHENIUS
    reactants: {H2SO4_EQ: 1.0}
    products: {SO4_aitken: 1.0}
    A: 1.0e-3
    stiff: true
  - type: ARRHENIUS
    reactants: {H2SO4_EQ: 1.0}
    products: {H2SO4_EQ: 1.0, N_aitken: 0.01}
    A: 1.0e-3
    stiff: true
```

With gas/material in equivalent formula units/cm3 and number in particles/cm3,
these rows give J=1e-5*Cgas, material production=100*J and gas loss=100*J.
The second reaction is number bookkeeping: its gas reactant/product cancel.
The chosen seed has 100 material units/particle solely for illustration. For a
configured seed radius/composition, compute its material-unit yield from seed
mass and molecular weight; do not retain 100 independently of those parameters.
Putting 100 on the gas reactant would incorrectly create a 100th-order rate law.
`IModalCoupling::nucleationTerms` supplies configured source terms without
advancing state; use those contributions instead of these two rows for
the general empirical exponent/seed mixture. Host code must apply the matching
precursor sink and availability/positivity constraints once.

`H2SO4_EQ` counts SO4-equivalent units; it is not kg of physical H2SO4.
The gas mapping uses target SO4 molecular weight, while gas transport uses
physical H2SO4 molecular weight. Hydrogen/water and any neutralization partners
must be accounted for by the host mechanism. These rows demonstrate equivalent
sulfate and particle-number bookkeeping, not complete elemental chemistry.

### Aging transfers every material and number

```yaml
reactions:
  - type: ARRHENIUS
    reactants: {BC_primary_carbon: 1.0}
    products: {BC_accumulation: 1.0}
    A: 2.0e-5
    stiff: true
  - type: ARRHENIUS
    reactants: {N_primary_carbon: 1.0}
    products: {N_accumulation: 1.0}
    A: 2.0e-5
    stiff: true
```

The complete file repeats the same frequency for **all seven** configured
materials and number. Each row transfers a single quantity conservatively;
aging does not consume the coating material again. The constant frequency only
demonstrates stoichiometry. In the coupled model, replace all these rows with
`modalAgingTerms`, using the condensation RHS/Jacobian from `modalPhaseExchange`
and the gross incoming coagulation coating budget. Recompute the common
frequency and its derivatives at every chemical stage. Freezing a frequency in
YAML omits its geometry, coating, and thermodynamic Jacobian terms.

### Host callback schedule and processes outside the reaction list

At each stage, gather and convert the mechanism state to EX-aero SI units;
evaluate exchange, configured nucleation, and aging; then scatter their rates
and Jacobians with the unit factors described below. Add them to the ordinary
chemical reactions. The current MKPP generated solver needs host integration
work to call this combined RHS/Jacobian; the YAML alone does not install it.

Keep coagulation in EX-aero, passing its accepted gross gains to aging. Keep
`renameModalModes` after chemical acceptance with `transfer_owner: host`.
Neither operation belongs in the reaction list for this split schedule. Water
uptake/optical diagnostics also remain EX-aero diagnostics. Gas oxidation and
aqueous chemistry use ordinary balanced MKPP reactions as appropriate.
For nitrate, choose explicitly between ionic NO3 mass and complete NH4NO3 salt
mass and include NH3/HNO3 partners consistently; one bound vapor is not a
two-precursor ammonium-nitrate equilibrium model. Dust is carried by the same
aging transfer rules when configured, without an invented dust chemistry reaction.

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
using Kokkos 5.1.1. These checks validate the coupling and generic closure;
calibrated MAM4 aging equivalence, implicit generated-mechanism integration,
thermodynamic calibration, and GPU validation remain separate work. The
[dry-box comparison harness](modal-comparison.md) measures process differences
and can execute a generated MKPP mechanism with frozen host stage coefficients.
