# CAM-SIMA aerosol parity roadmap

The design under review is [CAM-SIMA aerosol parity](../specs/2026-10-28-cam-sima-aerosol-parity-design.md).
Implement it as two independently reviewable workstreams, using CAM-SIMA
`development` commit `52723a73a2a49e37fdd981f795d08084bb5bfda8` as the numerical
reference. Ex-Aero starts at `codex/native-aerosol-optics` commit `5b347cf`.
The comparison also covers the CAM-SIMA-pinned `ESCOMP/atmospheric_physics`
submodule at `6ee9b89afb294408529bea5b72c7bab47c0c55eb`, including its
wet-deposition core; its rate physics remains host-owned in Ex-Aero.

1. [Modal diagnostics, physics, and host-removal adapters](2026-10-29-modal-diagnostics-physics-parity.md)
   closes the reference crosswalk and tests number limits, activation,
   scavenging size, externally computed dry-deposition velocities and
   wet-scavenging coefficients specific to bulk species/bins and modal
   mass/number moments, with common constituent-mass removal within each mode,
   resuspension/deposition transforms, cloud/ice diagnostics, and bulk analogs
   where the representation supports them.
2. [Modal RRTMG band optics](2026-10-29-modal-rrtmg-band-optics.md) completes
   deterministic SW/LW band data for generic modal and MAM4xx where source
   indices cover the bands, without wavelength substitution or extrapolation.

Run the behavior-parity work first to define test fixtures and shared semantics.
Optics can then proceed independently against its pinned data recipe. Finish
with the complete CMake test suite, optional MAM4xx configuration, Fortran
interop tests, and an updated status matrix. Neither workstream changes MKPP's
ownership of chemistry, equilibrium, or integration.
