# Test Fixtures — GEOSmie MIE Attribute Surfacing

Fixtures for `specs/001-geosmie-lut-attributes` (US2/US3 runtime-file paths, polish gates).
All fixtures are committed and deterministic; tests reference them by relative path from the
CTest working directory (`tests/`), never by absolute path.

## Conventions

- **Format**: the native portable text schema parsed by `MieFileLoader` (SPEC-CRTM-004
  stream pattern) — `#`/`!` comment lines, one `# @` metadata header, a dimension line, then
  whitespace-separated records. Fortran-style exponents (`1.0e-05`, `1.0d-05`) are tolerated.
- **Precision**: values are written with shortest round-trip `float64` repr so file-parsed
  values equal the source bits exactly (FR-005 grid-point tolerance 1e-7).
- **Units** are declared per attribute in the metadata header and MUST match the documented
  unit in `specs/001-geosmie-lut-attributes/data-model.md`, otherwise initialisation aborts
  with `FATAL ERROR:` (FR-009).
- **Golden reference values** are never duplicated here; tests compare against the pinned
  snapshot under `tools/geosmie_snapshot/` widened to `float64`.

## Files

| File | Species | Purpose |
|------|---------|---------|
| `mie_dust_monochromatic.txt` | DU | Valid file replacing the baked-in 30-band RRTMG λ axis with a 3-wavelength monochromatic set (0.355/0.55/1.33 µm) + polarized phase-function moments (rank-5, 6 elements × 3 moments; a file-only category). Proves file values *extend* baked bands (US2 T029, US3 T031/T034, invariant C7). |
| `mie_dust_override.txt` | DU | Valid file overriding `bext` (constant 1234.5) on the SAME 30-band axis; query must report `delivery_source=runtime-file` (invariant C7). Doubles as the file layer of the C13 precedence test. |
| `mie_nitrate_density.txt` | NI | Valid file supplying `wet_particle_density` + `growth_factor` (rank-2 radius,rh), which the pinned RRTMG band tables do not carry, proving config/file-supplied microphysics are surfaced (FR-001, US1). |
| `mie_dust7.txt` | DUST7 | Valid file for a species ABSENT from the baked set, declaring PRIME/non-source extents: 7 radius × 7 RH × 3 bands. Because the label is new, the file fully defines the curve — the no-hardcoding gate (invariant C12, T044). |
| `mie_corrupt_schema.txt` | DU | Unknown `# @bogus` directive → must abort with `FATAL ERROR:` and no silent fallback (invariant C6, FR-009). |
| `mie_corrupt_unit.txt` | DU | Declares `qext` in `m` instead of the canonical `1` → unit-mismatch abort (FR-009). |
| `mie_corrupt_version.txt` | DU | Missing `@source_version` → provenance abort (FR-016). |
| `mie_truncated.txt` | DU | Declares 2×2 axes but supplies 3 field values → truncation abort (FR-009). |

## Format notes

Files are self-describing: the loader derives each field's element count from the axes
declared **in the file** (`# @axis <name> <units> <n> [vals...]`), so valid fixtures emit
radius/rh even when identical to baked. The store merge keeps baked radius/rh and replaces
only what the file redefines; a file that changes the λ axis drops stale rank-3 fields it
does not re-supply (FR-012). The trailing count token on `# @field` lines is informational
readability only.

## Generating / refreshing

Regenerate the valid fixtures with `python3 tools/make_test_fixtures.py` (deterministic,
derived from the pinned snapshot). The corrupt fixtures are hand-authored and intentionally
invalid — do not "fix" them.
