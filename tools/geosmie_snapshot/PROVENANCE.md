# GEOSmie Optics Table Snapshot — Provenance

**Pinned release**: UFS regression-test input data `input-data-20260617`, GOCART `p8c_5d`
experiment (NASA GOCART2G optics tables produced by the GEOSmie toolchain).

**Authoritative source**: `https://noaa-ufs-regtests-pds.s3.amazonaws.com/input-data-20260617/GOCART/p8c_5d/ExtData/`
(NOA/UFS public regression-test data bucket, S3, anonymous read).

**Generator lineage**: NASA GMAO GEOSmie (`github.com/bbakernoaa/geosmie` / `GEOS-ESM/GEOSmie`),
Mie/GSF table production per Kemppinen et al. 2022 (GMAO pubs docs `Kemppinen1447`).

**License**: NASA GOCART/GEOSmie tables — Apache-2.0 lineage; attribution retained per
Ex-Aero spec Assumptions (Licensing/provenance).

## Committed files (RRTMG band-averaged tables — the baked-in default source)

These are small (< 200 KB each) and committed directly so the build-time generator is
hermetic and offline-reproducible (FR-016, SC-007).

| File | Species (label) | SHA-256 |
|------|-----------------|---------|
| `opticsBands_BC.v1_3.RRTMG.nc` | Black carbon (v1.3) | `b90131256d0ad57b3175afbd684a16c851fde8b66bb4439fbd7ef30fb73defe0` |
| `opticsBands_DU.v15_3.RRTMG.nc` | Dust (v15.3) | `644928a0a0ec39bc98cdc5805616230cbefff973afed2f8c377229ef3c4693de` |
| `opticsBands_NI.v2_5.RRTMG.nc` | Nitrate (v2.5) | `6efdbcd64a31fdb3d52755969918ac1b1313c940bb587c7d325e74fa85ee5e6e` |
| `opticsBands_OC.v1_3.RRTMG.nc` | Organic carbon (v1.3) | `453f270d9c8cfce60771bb8c3c00329b2e13a6c0812bd909c0a1546efdc67699` |
| `opticsBands_SS.v3_3.RRTMG.nc` | Sea salt (v3.3) | `0aff4b0664d031c24d2bb0718a79d70aa4d223ff5309e65ed812bdbb507f7cc7` |
| `opticsBands_SU.v1_3.RRTMG.nc` | Sulfate (v1.3) | `00c2dead267b289b9e389fe9c05037e35b22ea85e1112ec817363f0f59984d7d` |

## NOT committed (monochromatic tables — runtime-file delivery path)

The monochromatic tables (full wavelength set + polarized phase-function moments) are
73 MB–317 MB each (~1 GB total) and are fetched on demand by `fetch_monochromatic.sh` into
this directory (the directory is git-ignored for `optics_*.nc`). They are consumed by the
optional runtime portable-file path (FR-012), not baked in (FR-017).

| File | Species | SHA-256 |
|------|---------|---------|
| `optics_BC.v1_3.nc` | Black carbon | `c395abf32c5c8aca8a6a1fb9d9e4fbad6edcfef50f54b27e6fa66b58fdecb73b` |
| `optics_BRC.v1_5.nc` | Brown carbon | `e357a56dfcbd673f29b6f7935ac0faa8de23bbc0ed953596cec20738d7a2b427` |
| `optics_DU.v15_3.nc` | Dust | `933cba60cb2921aaa382bb28f63bf4276512b10f0d51fa6f81a83166e7ffd558` |
| `optics_NI.v2_5.nc` | Nitrate | `45e4e2d3a5234b509d72b610bf87c880ed916b9e775a1ef4b9340e454334b318` |
| `optics_OC.v1_3.nc` | Organic carbon | `d3a2a2424e560fecd07ffc5898e1d53bf0c74c8018c7a49ea1711bec64849cdc` |
| `optics_SS.v3_3.nc` | Sea salt | `313e1e8e05fabc99de458c9b560ac8733c9b4a5a1b3457cfb9e07b3563edc682` |
| `optics_SU.v1_3.nc` | Sulfate | `d41c16bf8703d8c2b841cccc2da5d0f65141034a5b0cb8f831337c24061c96f8` |

## Verified content notes (from schema inspection, 2026-09-04)

- RRTMG band tables: dims `rh=36` (fraction 0–0.99), `lambda=30` (band **indices** 1–30:
  14 SW + 16 LW, edges from `bandaverage.py getBands('RRTMG')`), `radius` = 5 (DU/SS),
  3 (NI), 2 (BC/OC/SU). Variables: `rLow rUp rEff rMass qext qsca bext bsca bbck g
  refreal refimag` (+ `qname`, coords).
- Monochromatic tables: same base variables plus `pmom (nPol=6, nMom, radius, rh, lambda)`
  and `pback`, at real wavelengths (56 values for DU, 61 for others, 250 nm–40 µm);
  `nMom` = 301 (DU), 751 (BC/OC/SS/SU), 2001 (NI/BRC).
- `rhop` + `growth_factor` are present ONLY in the NI and BRC monochromatic tables;
  `volume`/`area`/`mass`/`qabs`/`ssa`/`lidar_ratio` are absent from both sets (surfaced as
  explicit "not available" or exact derivations — see feature data-model.md).
- **No RRTMG band file exists for BR/BRC** in this pinned release: BR's spectral data is
  monochromatic-only (delivered via the runtime file path).

## Upgrading the pinned version

Per FR-016: replace the files above, update the SHA-256 table, and re-run
`tools/generate_mie_tables.py`; CI diffs the regenerated header. The provenance version
string below is what every baked-in record carries.

`SOURCE_VERSION = ufs-regtests-input-data-20260617/GOCART/p8c_5d`
