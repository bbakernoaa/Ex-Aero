# Comparing configurable modal processes with MAM4xx

`exaero_modal_mam4_comparison` runs paired dry boxes with identical initial mass,
number, size distributions, gas, temperature and pressure. It compares the
configurable EX-aero interfaces with pinned MAM4xx kernels. This is a process
comparison with a common box integrator, not a full native `amicphys` integration.

## Experiments

| Experiment | Generic path | Native reference |
| --- | --- | --- |
| `coagulation_control` | Native coagulation through the generic adapter | Direct `mam_coag_1subarea` |
| `coagulation` | Brownian quadrature, all six mode pairs | Native selected mode pairs |
| `condensation` | Public dry transition-regime uptake | Native beta-adjusted two-point uptake quadrature |
| `nucleation_linear` | Configured first-order empirical helper | Wang PBL option 11 |
| `nucleation_quadratic` | Configured second-order empirical helper | Wang PBL option 12 |
| `nucleation_example` | Illustrative `J=1e-8*C**1.5` | Wang PBL option 11 |
| `aging` | Configurable sulfate shell coating | Native finite-step primary-carbon aging |
| `combined` | Brownian + uptake + illustrative nucleation + shell aging | Native kernels + Wang option 11 |
| `combined_wang_linear` | Same generic processes with matched first-order nucleation | Same native reference |

In these formulas C is molecules/cm3 and J is particles/cm3/s. The helper's SI
prefactor is converted explicitly. Matched Wang and native-coagulation controls
enforce parity; other cases measure differences without imposing an arbitrary
scientific acceptance limit.

The four modes and seven materials use preset widths/densities. Primary carbon
starts with POM and BC; other modes start with sulfate. Three cases use
(250 K, 60000 Pa, 1e5 acid molecules/cm3), (280 K, 90000 Pa, 1e6), and
(300 K, 101325 Pa, 1e7). Initial radii scale by 0.7, 1 and 1.5 respectively.
RH is zero. H2SO4 diffusivity and accommodation (0.65) are shared; quadrature
and molecular-speed prescriptions remain those of each implementation.

## Material and integration conventions

Both paths use a **controlled 100-H2SO4-molecule seed**. The generic nucleation
helper computes actual pure-acid mass; the adapter maps it to SO4-equivalent dry
mass. Public gas also uses SO4-equivalent kg/m3, with 0.09606 kg/mol of acid
equivalents. Vapor transport uses physical H2SO4, 0.098079 kg/mol. This compares
cluster formation J and material bookkeeping, **not native cluster survival and
growth into Aitken**. Wang's nominal 1 nm cluster is not silently equated to the
controlled seed.

Generic aging uses an exact lognormal shell of thickness `8*4.76e-10 m`.
Native aging uses eight sulfate monolayers and upstream equivalent-volume and
finite-step transfer rules. Its private mol/mol conversion explicitly uses its
effective molecular weights, including 115 g/mol sulfate; conversion back
preserves each species' mass. Native aging transfers coating species fully and
core mass/number by the aged fraction; generic aging transfers the same fraction
of every source quantity. Their results are expected to differ. The aging-only
case adds and budgets a 1e-15 kg/m3/s sulfate coating source.

Every split step applies coagulation, explicit condensation/nucleation, then
aging. Accepted gross coagulation gains feed aging. Runs use dt=1 and 0.5 s,
with output at 0, 1, 60, 300 and 600 s. Every step checks finite/nonnegative
state and per-species material conservation, including gas and coating sources.
First-step changes are effective tendencies; native finite updates are not
presented as instantaneous analytic derivatives.

## Generated MKPP execution

`EXAERO_COMPARISON_MKPP=ON` compiles
[`reactions.yaml`](../../tools/modal_comparison/reactions.yaml) using the MKPP
compiler from [MKPP PR #2](https://github.com/bbakernoaa/MKPP/pull/2). It explicitly
represents gas-to-mode transfer, balanced seed formation, and all mass/number
aging transfers. Generated kernels execute at every chemical stage and drive
the state update. Every physical rate must match an independent host RHS, and
fixed coefficients must have zero tendency. Material and particle coordinates
convert to/from units/cm3.

The host refreshes uptake, formation J and aging coefficients at each stage.
This benchmark freezes geometry within a stage; it does not exercise fully
symbolic geometry feedback, the dynamic three-mode educational example, or
MKPP's implicit solver. Those require a separate solver comparison. With MKPP
OFF, the same boxes use the host RHS and outputs explicitly say
`host_reference`; enabled builds say `mkpp_generated`.

## Build and run

Use Python with numpy, pandas and the EX-aero table-generator dependencies.
For generated chemistry, install the MKPP PR branch into that environment or
expose its `src` through `PYTHONPATH`. Kokkos 5.1.1, EKAT, yaml-cpp, Catch2 v2
and a Fortran compiler are required as in the main build.

```sh
cmake -S . -B build/comparison -DEXAERO_ENABLE_MAM4XX=ON \
  -DEXAERO_BUILD_CHEMISTRY=OFF -DEXAERO_COMPARISON_MKPP=ON \
  -DPython3_EXECUTABLE="$(command -v python)" \
  -DCMAKE_PREFIX_PATH="$CMAKE_PREFIX_PATH"
cmake --build build/comparison --target exaero_modal_mam4_comparison --parallel 2
ctest --test-dir build/comparison -R '^ExaeroModalMam4Comparison$' \
  --output-on-failure
export DATA="${DATA:?Set DATA to the comparison output directory}"
python tools/modal_comparison/report.py \
  --executable build/comparison/tests/exaero_modal_mam4_comparison \
  --mkpp-manifest build/comparison/tests/comparison_generated/modal_mam4_comparison_manifest.json
```

Use a clean comparison build when changing the MKPP compiler revision. The
normal MAM4xx-enabled CTest build includes the ten-second smoke comparison
without requiring MKPP. The report driver runs both full 600-second boxes
and writes these artifacts beneath DATA:

- `modal.comparison.csv` and `.half.csv`: paired state trajectories.
- `modal.comparison.first_step.csv`: effective first-step tendencies.
- `modal.comparison.summary.csv`: every final mass/number/gas/radius field,
  absolute and symmetric relative differences, and time-step sensitivity.
- `modal.comparison.md`: selected number, BC and gas results with scope notes.
- `modal.comparison.metadata.json`: source revisions, dirty-tree flag,
  executable/reaction checksums, cases, time steps and optional MKPP manifest.

Pass `--mkpp-source PATH_TO_MKPP_CHECKOUT` to additionally record the compiler's
Git revision and dirty-tree flag. The report also hashes the comparison sources
and generated rate header, so an uncommitted local run remains identifiable.

## Initial measured results

The initial GCC 14/Kokkos 5.1.1 OpenMP run executed generated MKPP rates and
passed all control, conservation and positivity checks across the 27 paired
experiments. At 280 K/90000 Pa, after 600 s with dt=0.5 s:

| Observable | Generic | Native kernel reference | Symmetric relative difference |
| --- | --- | --- | --- |
| Matched-law combined Aitken number [particles/m3] | 1.281264e9 | 1.281913e9 | 0.0507% |
| Matched-law combined acid vapor [SO4-equivalent kg/m3] | 1.017740e-13 | 1.018431e-13 | 0.0679% |
| Illustrative-law combined Aitken number [particles/m3] | 4.020309e9 | 1.281913e9 | 68.1% |
| Aging-only primary-carbon number [particles/m3] | 9.350116e7 | 9.271349e7 | 0.842% |

For the matched-law combined Aitken number, halving dt changed each path by
about 0.0057%; the acid-vapor changes were about 0.0080%. The isolated matched
linear/quadratic nucleation and coagulation-adapter controls agree to the
enforced numerical tolerances. These results describe the defined dry cases;
they do not establish broader atmospheric or full MAM4 equivalence.

Two time steps help distinguish integration sensitivity from formulation
differences but do not prove convergence. Wet condensation, organic partition
thermodynamics, native renaming, clouds, deposition, optics, and GPU execution
are not covered by these dry-box cases.
