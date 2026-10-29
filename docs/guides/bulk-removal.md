# Bulk removal coupling

See [host removal tendencies](host-removal.md) for host-supplied deposition
velocities and first-order wet-scavenging coefficients through C++, C and Fortran.

The GOCART bulk backend packs one mass field per configured species. These utilities provide size and removal budgets for a host model; they do not calculate scavenging, resuspension, or deposition rates. All operations are explicit and leave aerosol state unchanged.

`bulkFixedNumber` returns a diagnostic particle number concentration for each species, using the same fixed lognormal mass-to-number conversion as the existing GOCART diagnostics:

\[
N_s = \frac{M_s}{(\pi/6)\rho_s d_{g,s}^3\exp(4.5\ln^2\sigma_s)}.
\]

Here `M` is dry aerosol mass concentration, `rho` is configured dry density, `d_g` is configured lognormal geometric mean diameter, and `sigma` is configured geometric width. Number is particles/m3 and is not a prognostic field. There is no bulk analog of modal number-bound enforcement.

`bulkReferenceRadii` returns **physical dry radii** independent of GEOSmie or radiation tables. Optional per-species `dry_radius_lower` and `dry_radius_upper` must be supplied together in YAML or structured configuration; omitted bounds report `Unavailable` and zero. `dry_effective_radius` may be supplied explicitly; otherwise it is derived from the configured lognormal moments as `(lognormal_dg/2) * exp(2.5 * ln(lognormal_sigma)^2)` and reports `Derived`. The returned `(species,3)` values and statuses are lower, upper, effective, in metres. A supplied effective radius does not alter the existing mass-to-number formula.

`bulkScavengingSizes` reports configured dry particle diameter and kappa-equilibrium wet diameter by species and point. It uses configured hygroscopicity and host RH capped at 0.99, matching GOCART's analytical growth path; it does not query GEOSmie. An empty mass field has `BULK_SIZE_EMPTY` status and zero diameters. Populated fields have `BULK_SIZE_CONFIGURED` status. The configured `dry_particle_diameter` used here can differ from `lognormal_dg` used for number. The pre-existing `computeWetSizing` API remains a separate table-backed option for hosts that explicitly choose it.

`redistributeBulkResuspension` accepts nonnegative returned-material rates `(cell,level,species)` in kg/m3/s and explicit source/destination species-index pairs. Every source may route once, from the input snapshot. Unrouted mass stays in its original species field; cycles, self routes, and duplicate sources are rejected. A host decides when to apply the returned rate after accepting its step.

`rebinBulkDeposition` accepts nonnegative, time-integrated deposited budgets in kg/m2, selects one configured species, and reports dry-diameter bin mass and below/above tails in kg/m2. It uses the species' fixed configured lognormal distribution, or a single size when `lognormal_sigma` is exactly one. This assumes deposited particles follow the configured distribution, which may be inaccurate after size-selective removal. No mass is silently dropped outside the requested edges. Species indices on C/C++ boundaries are zero-based.

C and Fortran entry points mirror these operations with column-major arrays and `errmsg` buffers of at least 256 characters. In Fortran, `radii(nf,3)` and its status use lower, upper, effective order. The host owns rollback and only applies returned rates after accepting the step.
