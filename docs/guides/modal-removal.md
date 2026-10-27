# Generic modal removal and number bounds

These are explicit host utilities for the configurable `modal` backend. They do not run during `executeMicrophysics` or a chemical solver stage. The fixed `mam4xx` package has a separate state layout and does not use these interfaces.

Each mode may set `dry_geometric_radius_min` and `dry_geometric_radius_max` in metres, in a `modes` entry or `mode_overrides`. Both values are required together; omitted bounds are disabled. After an accepted chemistry and removal update, call `enforceModalBounds` once. It preserves every material mass and adjusts only particle number to keep dry geometric radius in range. Record the signed returned `delta_number` as a correction budget. An empty mode remains empty; the host must provide both mass and number when injecting particles. The correction is non-differentiable at the bounds and is not a chemical tendency.

`inspectModalBounds` returns dry geometric radius, minimum permitted number, maximum permitted number, and a status for each mode. `modalScavengingSizes` returns dry and wet **number-median diameters** in metres. A host scavenging scheme chooses the size measure appropriate to its own collection law; these calls do not compute a removal rate.

Optional `resuspension_routes` name a dry material, source mode, and destination mode:

```yaml
resuspension_routes:
  - {species: DU, source: accumulation, destination: coarse}
```

`redistributeModalResuspension` accepts nonnegative returned-material **rates** in kg/m3/s. Number slots must be zero. Each route reads the original buffer once, so chains do not route a contribution twice. Unrouted material remains at its source. The host supplies a separate returned-number rate and applies both mass and number after accepting the step. Moving material into an empty destination without a positive number source is invalid ambient state.

`rebinModalDeposition` accepts time-integrated deposition **budgets**: kg/m2 for materials and particles/m2 for number. It reconstructs a dry fixed-width lognormal from the deposited total material volume and number, then returns selected-material kg/m2 in each requested dry geometric diameter bin. Two tail outputs preserve material below and above the requested edges. The deposited population can differ from ambient state after size-selective removal. Wet and aerodynamic size cuts require a separate model.

Typical host order: apply accepted chemistry and removal changes to ambient mass and number, apply optional host-owned growth transfer, enforce bounds once, record the number correction, and refresh properties. The host owns rollback and must not replay corrections or returned mass on a rejected step. The C and Fortran entry points use column-major `(cell,level,field)` buffers and report errors through `errmsg` and `errflg`.

For example, a host with an accepted removal step can route a returned dust mass rate with `redistributeModalResuspension`, then integrate that rate alongside its independently calculated returned-number rate. After updating both ambient fields, it calls `enforceModalBounds` and adds `delta_number` to its number-correction diagnostic. It calls `modalScavengingSizes` again before calculating the next removal step. A rejected step discards the proposed rates and does not call the bound correction. The host should size `errmsg` to at least 256 characters for each C or Fortran call; species selectors in C and C++ are zero-based. The status codes are exported as `EXAERO_MODAL_BOUNDS_*` in C and `exaero_MODAL_BOUNDS_*` in Fortran.
