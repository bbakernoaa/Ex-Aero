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
| `mie_dust_monochromatic.txt` | DU | Valid file extending baked-in RRTMG bands with monochromatic wavelengths + polarized phase-function moments (US2 T029, US3 T031/T034). Proves file values *extend* baked bands (invariant C7). |
| `mie_dust_override.txt` | DU | Valid file overriding a key that also exists baked-in; query must report `delivery_source=runtime-file` (invariant C7). |
| `mie_nitrate_rhop.txt` | NI | Valid file supplying `rhop` + `growth_factor`, which the pinned RRTMG tables do not carry, proving config/file-supplied microphysics are surfaced (FR-001, US1). |
| `mie_corrupt_schema.txt` | DU | Missing required variable/dimension ordering wrong → must abort with `FATAL ERROR:` and no silent fallback (invariant C6, FR-009). |
| `mie_corrupt_unit.txt` | DU | Declares `bext` in `m2 g-1` instead of `m2 (kg dry mass)-1` → unit-mismatch abort (FR-009). |
| `mie_corrupt_version.txt` | DU | Empty/absent `source_version` → provenance abort (FR-016). |
| `mie_truncated.txt` | DU | Declares more records than present → truncation abort (FR-009). |

## Generating / refreshing

Regenerate the valid fixtures with `python3 tools/make_test_fixtures.py` (deterministic,
derived from the pinned snapshot). The corrupt fixtures are hand-authored and intentionally
invalid — do not "fix" them.
