/// @file GocartSolver.cpp
/// @brief Kokkos device kernels for the GOCART/GEOSmie aerosol package.
///
/// Contains the four hot-path kernels (diagnostics, optics, CCN
/// activation, emissions), the device solver state (species parameters +
/// the single flat curve pool), the public Kokkos environment lifecycle
/// wrappers, and the device-safe table-lookup helpers shared by the
/// kernels.
///
/// @section physics Physics summary
/// - **Diagnostics:** lognormal number-from-mass inversion, wet growth
///   (GEOSmie growth factor when curve-bound, else @f$ \kappa @f$-Köhler
///   @f$ g = (1 + \kappa \cdot \mathrm{RH}/(1-\mathrm{RH}))^{1/3} @f$),
///   SAD, aerosol liquid water, Stokes settling velocity, PM2.5/PM10 cuts.
/// - **Optics:** three resolution modes per species — Mode C
///   (GEOSmie spectral table), Mode B (band-integrated RH lookup), Mode A
///   (analytical Anomalous Diffraction Theory) — with
///   Henyey-Greenstein backscatter and Beer–Lambert column AOT
///   integration over @f$ \Delta z @f$.
/// - **CCN:** Köhler-theory critical supersaturation
///   @f$ S_c = \sqrt{4 A^3 / (27 \kappa D_{dry}^3)} @f$.
/// - **Emissions:** raw->target flux mapping with unit scaling and
///   modal mass-to-number conversion.
///
/// @section safeguards Numerical safeguards (applied in every kernel)
/// RH clamped to @f$ \le 0.99 @f$; interpolants never extrapolate beyond
/// declared axis points; guarded divisions with @c 1e-15/@c 1e-30 epsilon
/// denominators; NaN/Inf and negative fluxes zeroed at the boundary; the
/// small-@f$ \rho @f$ ADT branch uses a Taylor series to eliminate
/// floating-point cancellation.
///
/// @note Memory model: host arrays arrive as raw column-major pointers
/// and are wrapped in *unmanaged* Kokkos views (zero-copy). All table
/// data is device-resident in one pool uploaded at initialization — the
/// timestep loop performs no host-to-device transfers.
#include <Kokkos_Core.hpp>
#include <gocart/GocartSpeciesParams.hpp>

namespace exaero {

namespace {

// Device-safe bracket search on a strictly-increasing coordinate axis with
// clamp at both edges (never extrapolates beyond declared points).
// log_space != 0 blends in log coordinates (reference spectral interpolation,
//); otherwise linear. Mirrors the host MieTableStore locate rules.
/// @brief Locate the interpolation bracket @f$ [v_i, v_{i+1}] @f$ for @p v
/// on a strictly increasing axis.
/// @param axis Ascending coordinate axis of length @p n.
/// @param n Axis length (data, never a fixed literal).
/// @param v Query coordinate.
/// @param log_space Blend in log(@p axis) when non-zero (spectral axes);
///        0 blends linearly (RH axes).
/// @param[out] i0 Lower bracket index.
/// @param[out] i1 Upper bracket index (== i0 when clamped to an edge).
/// @param[out] w Blend weight in [0,1]; 0 at an edge clamp (no extrapolation).
KOKKOS_INLINE_FUNCTION void curve_locate(const double *axis, int n, double v,
                                         int log_space, int &i0, int &i1,
                                         double &w) {
  if (n <= 1) {
    i0 = 0;
    i1 = 0;
    w = 0.0;
    return;
  }
  if (v <= axis[0]) {
    i0 = 0;
    i1 = 0;
    w = 0.0;
    return;
  }
  if (v >= axis[n - 1]) {
    i0 = n - 1;
    i1 = n - 1;
    w = 0.0;
    return;
  }
  int lo = 0, hi = n - 1;
  while (hi - lo > 1) {
    const int mid = (lo + hi) / 2;
    if (axis[mid] <= v)
      lo = mid;
    else
      hi = mid;
  }
  i0 = lo;
  i1 = hi;
  const double v0 = axis[i0], v1 = axis[i1];
  if (log_space && v0 > 0.0 && v1 > 0.0 && v > 0.0) {
    const double l0 = Kokkos::log(v0);
    const double span = Kokkos::log(v1) - l0;
    w = (span > 0.0) ? (Kokkos::log(v) - l0) / span : 0.0;
  } else {
    const double span = v1 - v0;
    w = (span > 0.0) ? (v - v0) / span : 0.0;
  }
}

// GEOSmie table spectral read mass extinction [m^2/kg],
// single-scattering albedo and asymmetry at (rh, band coordinate) for the
// species' solver radius node. Block layout documented in GocartSpeciesParams;
// the extents are data, never literals. Returns false when the requested band
// coordinate lies OUTSIDE the curve's declared lambda domain (e.g. a physical
// wavelength in metres against a band-index axis): the caller must then fall
// back to its analytical path rather than silently clamping to an unrelated
// band (fail-loud intent, no silent wrong physics).
KOKKOS_INLINE_FUNCTION bool
species_spectral_table(const GocartSpeciesParams &params, const double *pool,
                       double rh, double band, double &ext_per_mass,
                       double &ssa_v, double &g_v) {
  const double *sb = pool + params.spec_offset;
  const int nH = params.n_spec_rh;
  const int nL = params.n_spec_lambda;
  const int nR = params.n_spec_radius;
  const double *srh = sb;
  const double *slam = sb + nH;
  if (!(band >= slam[0] && band <= slam[nL - 1]))
    return false;
  const long nfield = static_cast<long>(nR) * nH * nL;
  const double *f_ext = sb + nH + nL;
  const double *f_ssa = f_ext + nfield;
  const double *f_g = f_ssa + nfield;
  int rnode = params.solver_radius_node;
  if (rnode < 0)
    rnode = 0;
  if (rnode > nR - 1)
    rnode = nR - 1;
  int h0, h1, l0, l1;
  double wh, wl;
  curve_locate(srh, nH, rh, 0, h0, h1, wh);
  curve_locate(slam, nL, band, 1, l0, l1, wl);
  const long base = static_cast<long>(rnode) * nH;
  const double e00 = f_ext[(base + h0) * nL + l0];
  const double e10 = f_ext[(base + h1) * nL + l0];
  const double e01 = f_ext[(base + h0) * nL + l1];
  const double e11 = f_ext[(base + h1) * nL + l1];
  const double s00 = f_ssa[(base + h0) * nL + l0];
  const double s10 = f_ssa[(base + h1) * nL + l0];
  const double s01 = f_ssa[(base + h0) * nL + l1];
  const double s11 = f_ssa[(base + h1) * nL + l1];
  const double g00 = f_g[(base + h0) * nL + l0];
  const double g10 = f_g[(base + h1) * nL + l0];
  const double g01 = f_g[(base + h0) * nL + l1];
  const double g11 = f_g[(base + h1) * nL + l1];
  const double e_l = e00 + wh * (e10 - e00);
  const double e_h = e01 + wh * (e11 - e01);
  ext_per_mass = e_l + wl * (e_h - e_l);
  const double s_l = s00 + wh * (s10 - s00);
  const double s_h = s01 + wh * (s11 - s01);
  ssa_v = s_l + wl * (s_h - s_l);
  const double g_l = g00 + wh * (g10 - g00);
  const double g_h = g01 + wh * (g11 - g01);
  g_v = g_l + wl * (g_h - g_l);
  return true;
}

// Device read of one polarized phase-function moment. The
// (element, moment) pair is the documented encoding idx = moment*ELEMENT_STRIDE
// + element (data-model); the caller supplies element in [0,nP) and moment in
// [0,nM). Returns false when the species has no moment block or the band is
// outside the curve lambda domain (declines, never a silent wrong value). Value
// interpolates linearly in RH and takes the nearest grid band (matches the host
// store query).
KOKKOS_INLINE_FUNCTION bool
species_moment_table(const GocartSpeciesParams &params, const double *pool,
                     double rh, double band, int element, int moment,
                     double &value_out) {
  if (params.pmom_offset < 0)
    return false;
  const int nH = params.n_pmom_rh;
  const int nL = params.n_pmom_lambda;
  const int nR = params.n_pmom_radius;
  const int nP = params.n_pmom_pol;
  const int nM = params.n_pmom_moment;
  if (element < 0 || element >= nP || moment < 0 || moment >= nM)
    return false;
  const double *pb = pool + params.pmom_offset;
  const double *srh = pb;
  const double *slam = pb + nH;
  if (!(band >= slam[0] && band <= slam[nL - 1]))
    return false; // domain decline
  int rnode = params.solver_radius_node;
  if (rnode < 0)
    rnode = 0;
  if (rnode > nR - 1)
    rnode = nR - 1;
  const double *pmom = pb + nH + nL;
  // C-order ravel over (radius,rh,lambda,pol,moment): moment is the fastest
  // axis.
  auto at = [&](int hh, int ll) {
    return pmom[(((static_cast<long>(rnode) * nH + hh) * nL + ll) * nP +
                 element) *
                    nM +
                moment];
  };
  int h0, h1;
  double wh;
  curve_locate(srh, nH, rh, 0, h0, h1, wh);
  int l0, l1;
  double wl;
  curve_locate(slam, nL, band, 1, l0, l1, wl);
  const int l = (wl < 0.5) ? l0 : l1; // nearest grid band (grid read)
  const double v0 = at(h0, l);
  const double v1 = at(h1, l);
  value_out = v0 + wh * (v1 - v0);
  return true;
}

} // anonymous namespace

// Implement public environment lifecycle routines (declared in
// exaero/Environment.hpp)
// @copydoc exaero::initialize_environment
void initialize_environment() {
  if (!Kokkos::is_initialized()) {
    Kokkos::initialize();
  }
}

// @copydoc exaero::finalize_environment
void finalize_environment() {
  if (Kokkos::is_initialized()) {
    Kokkos::finalize();
  }
}

// @copydoc exaero::fence_environment
void fence_environment() { Kokkos::fence(); }

// Definition of our private solver state managing GPU-allocated memory views
/// @brief Device-resident solver state: species parameters + curve pool.
///
/// The only GPU allocations owned by the package; both are sized once at
/// initialization so the timestep kernels are pure readers (zero H2D
/// traffic, deterministic timing).
struct GocartSolverState {
  int num_species; ///< Species axis length (POD mirror of the host count).
  Kokkos::View<GocartSpeciesParams *, Kokkos::DefaultExecutionSpace>
      d_species_params; ///< Device copy of the per-species parameter blocks.
  // ONE flat device double pool holding every species' curve block.
  // Uploaded once here; never reallocated or H2D-copied inside the
  // timestep loop.
  Kokkos::View<double *, Kokkos::DefaultExecutionSpace>
      d_curve_pool; ///< Device-resident flat curve pool (CurvePoolDescriptor
                    ///< offsets index into it).
};

// Extern C/C++ helper functions declared in GocartPackage.cpp
/// @brief Build the device solver state from host arrays (one-time upload).
///
/// Deep-copies the species parameter blocks and the entire flat curve
/// pool to the default execution space. After this call the timestep
/// loop performs zero host-to-device transfers.
/// @param num_species Species count (length of @p params).
/// @param params Host species parameter blocks.
/// @param pool Flat device curve pool (all species' blocks concatenated).
/// @param pool_size Length of @p pool (0 => no pool allocation).
/// @return Newly allocated solver state; caller owns it (free via
///         free_solver_state()).
GocartSolverState *create_solver_state(int num_species,
                                       const GocartSpeciesParams *params,
                                       const double *pool, int pool_size) {
  auto *state = new GocartSolverState();
  state->num_species = num_species;

  // Allocate Device View
  state->d_species_params =
      Kokkos::View<GocartSpeciesParams *, Kokkos::DefaultExecutionSpace>(
          "d_species_params", num_species);

  // Allocate Host View Mirror
  auto h_view = Kokkos::create_mirror_view(state->d_species_params);
  for (int i = 0; i < num_species; ++i) {
    h_view(i) = params[i];
  }

  // Deep copy values to GPU Default Execution Space
  Kokkos::deep_copy(state->d_species_params, h_view);

  // Single H2D upload of the flat curve pool (zero transfers in the loop,
  //).
  if (pool_size > 0) {
    state->d_curve_pool = Kokkos::View<double *, Kokkos::DefaultExecutionSpace>(
        "d_curve_pool", pool_size);
    auto h_pool = Kokkos::create_mirror_view(state->d_curve_pool);
    for (int i = 0; i < pool_size; ++i)
      h_pool(i) = pool[i];
    Kokkos::deep_copy(state->d_curve_pool, h_pool);
  }

  return state;
}

/// @brief Release the device solver state (safe with nullptr).
void free_solver_state(GocartSolverState *state) {
  if (state) {
    delete state;
  }
}

/// @brief Launch the derived-diagnostics kernels (3D fields + column/surface).
///
/// Per (cell, level) the kernel sums, over species: mass-weighted PM2.5 /
/// PM10 cuts (dry diameter <= 2.5 / 10 um), number concentration from
/// the lognormal volume relation
/// @f$ N = M_{bulk} / [ (\pi/6) \rho_0 D_{gm}^3 e^{4.5 \ln^2 \sigma_g} ] @f$,
/// surface area density, aerosol liquid water from the wet-minus-dry
/// volume, and a mass-weighted Stokes settling velocity
/// @f$ v_g = \rho_w D_w^2 g / (18 \mu) @f$. Wet growth uses the GEOSmie
/// table when the species is curve-bound, else the @f$ \kappa @f$-Köhler
/// approximation. A second kernel integrates the column
/// @f$ \sum_l \rho_l \Delta z_l @f$ and extracts the surface layer
/// (level 0) into the 2D diagnostic slots.
/// @param state Device solver state (params + curve pool).
/// @param num_cells Horizontal cell count.
/// @param num_levels Vertical level count.
/// @param num_species Species axis length of @p state_ptr.
/// @param rh_ptr [fraction] relative humidity, column-major (cell, level).
/// @param thick_ptr [m] layer thickness, column-major (cell, level).
/// @param state_ptr [kg/m^3] aerosol mass, column-major (cell, level, species).
/// @param[out] diags_ptr Column-major (cell, level, NUM_DIAGNOSTICS);
///             slots per exaero::diagnostic_indices.
void run_gocart_diagnostics(GocartSolverState *state, int num_cells,
                            int num_levels, int num_species,
                            const double *rh_ptr, const double *thick_ptr,
                            const double *state_ptr, double *diags_ptr) {

  using MemSpace = typename Kokkos::DefaultExecutionSpace::memory_space;

  // Wrap incoming raw pointers into unmanaged Default Execution Space
  // column-major (LayoutLeft) memory Views (zero copy!)
  auto d_rh = Kokkos::View<const double **, Kokkos::LayoutLeft, MemSpace,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
      rh_ptr, num_cells, num_levels);
  auto d_thick = Kokkos::View<const double **, Kokkos::LayoutLeft, MemSpace,
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
      thick_ptr, num_cells, num_levels);
  auto d_state = Kokkos::View<const double ***, Kokkos::LayoutLeft, MemSpace,
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
      state_ptr, num_cells, num_levels, num_species);
  auto d_diags = Kokkos::View<double ***, Kokkos::LayoutLeft, MemSpace,
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
      diags_ptr, num_cells, num_levels, diagnostic_indices::NUM_DIAGNOSTICS);

  auto d_species_params = state->d_species_params;
  auto d_curve_pool = state->d_curve_pool;

  // Execute parallel diagnostics calculation on the GPU
  Kokkos::parallel_for(
      "GocartDiagnosticsKernel",
      Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {num_cells, num_levels}),
      KOKKOS_LAMBDA(int i_cell, int i_level) {
        double total_pm2_5_3d = 0.0;
        double total_pm10_3d = 0.0;
        double total_number_3d = 0.0;
        double total_sad_3d = 0.0;
        double total_alw_3d = 0.0;
        double weighted_vg_numerator = 0.0;
        double total_mass_vg_weight = 0.0;

        double current_rh = Kokkos::min(d_rh(i_cell, i_level), 0.99);

        for (int i_spec = 0; i_spec < num_species; ++i_spec) {
          const auto &params = d_species_params(i_spec);
          double mass_conc = d_state(i_cell, i_level, i_spec); // [kg/m³]

          if (mass_conc <= 0.0)
            continue;

          // 0. GEOSmie table-backed microphysics when this species is
          // curve-bound, growth factor and wet particle density are read
          // from the flat device pool (RH-linear, extent-driven scan — no
          // literals).
          double gf_table = -1.0; // wet/dry effective-radius growth factor
          double wd_table = -1.0; // wet particle density [kg/m^3]
          if (params.micro_offset >= 0) {
            const double *mb = d_curve_pool.data() + params.micro_offset;
            const int nH = params.n_micro_rh;
            const int nR = params.n_micro_radius;
            const double *mrh = mb;                // [nH]
            const double *mgf = mb + nH;           // [nR*nH]
            const double *mwd = mb + nH + nR * nH; // [nR*nH]
            int rnode =
                Kokkos::min(Kokkos::max(params.solver_radius_node, 0), nR - 1);
            int i_bin = 0;
            while (i_bin < nH - 2 && current_rh > mrh[i_bin + 1]) {
              i_bin++;
            }
            double w = (current_rh - mrh[i_bin]) /
                       (mrh[i_bin + 1] - mrh[i_bin] + 1e-15);
            w = Kokkos::min(Kokkos::max(w, 0.0),
                            1.0); // clamp: never extrapolate
            gf_table =
                mgf[rnode * nH + i_bin] +
                w * (mgf[rnode * nH + i_bin + 1] - mgf[rnode * nH + i_bin]);
            wd_table =
                mwd[rnode * nH + i_bin] +
                w * (mwd[rnode * nH + i_bin + 1] - mwd[rnode * nH + i_bin]);
          }

          // 1. Derived Number Concentration from Bulk Mass (lognormal
          // population)
          double ln_sig = Kokkos::log(params.lognormal_sigma);
          double vol_factor = (M_PI / 6.0) * params.dry_density *
                              Kokkos::pow(params.lognormal_dg, 3) *
                              Kokkos::exp(4.5 * ln_sig * ln_sig);
          double num_conc = mass_conc / vol_factor; // [particles/m³]
          total_number_3d += num_conc;

          // 2. Wet size calculation: GEOSmie growth factor (table-backed)
          // when available, else the kappa-Kohler approximation.
          double growth =
              (gf_table > 0.0)
                  ? gf_table
                  : Kokkos::pow(1.0 + params.hygroscopicity *
                                          (current_rh / (1.0 - current_rh)),
                                1.0 / 3.0);
          double wet_diameter = params.dry_particle_diameter * growth;

          // 3. Surface Area Density (SAD) [m²/m³]
          double sad_factor = M_PI * Kokkos::pow(params.lognormal_dg, 2) *
                              Kokkos::exp(2.0 * ln_sig * ln_sig);
          total_sad_3d += num_conc * sad_factor;

          // 4. Aerosol Liquid Water (ALW) Content [kg/m³] (density of liquid
          // water = 1000 kg/m³)
          double dry_vol =
              (M_PI / 6.0) * Kokkos::pow(params.dry_particle_diameter, 3);
          double wet_vol = (M_PI / 6.0) * Kokkos::pow(wet_diameter, 3);
          double alw_mass_spec =
              num_conc * (wet_vol - dry_vol) * 1000.0; // [kg/m³]
          total_alw_3d += alw_mass_spec;

          // 5. Gravitational Settling Fall Velocity (vg) [m/s] (Stokes
          // Settling) Wet density: GEOSmie table value when curve-bound,
          // else the dry-volume + condensate-mass mixture approximation.
          double dry_mass = dry_vol * params.dry_density;
          double water_mass = (wet_vol - dry_vol) * 1000.0;
          double wet_density =
              (wd_table > 0.0) ? wd_table : (dry_mass + water_mass) / wet_vol;

          double g_acc = 9.80665;     // [m/s²]
          double dyn_visc = 1.825e-5; // [kg/m-s] air dynamic viscosity
          double vg_spec = (wet_density * wet_diameter * wet_diameter * g_acc) /
                           (18.0 * dyn_visc);

          weighted_vg_numerator += mass_conc * vg_spec;
          total_mass_vg_weight += mass_conc;

          // 6. 3D PM2.5 and PM10 size cuts
          if (params.dry_particle_diameter <= 2.5e-6) {
            total_pm2_5_3d += mass_conc;
          }
          if (params.dry_particle_diameter <= 10.0e-6) {
            total_pm10_3d += mass_conc;
          }
        }

        // Write 3D grid cell diagnostics
        d_diags(i_cell, i_level, diagnostic_indices::MASS_CONCENTRATION) =
            d_state(i_cell, i_level, 0); // principal species mass
        d_diags(i_cell, i_level, diagnostic_indices::PM2_5_CONCENTRATION) =
            total_pm2_5_3d;
        d_diags(i_cell, i_level, diagnostic_indices::PM10_CONCENTRATION) =
            total_pm10_3d;
        d_diags(i_cell, i_level, diagnostic_indices::NUMBER_CONCENTRATION) =
            total_number_3d;
        d_diags(i_cell, i_level, diagnostic_indices::SURFACE_AREA_DENSITY) =
            total_sad_3d;
        d_diags(i_cell, i_level, diagnostic_indices::AEROSOL_LIQUID_WATER) =
            total_alw_3d;
        d_diags(i_cell, i_level,
                diagnostic_indices::GRAVITATIONAL_SETTLING_VELOCITY) =
            total_mass_vg_weight > 0.0
                ? (weighted_vg_numerator / total_mass_vg_weight)
                : 0.0;
      });
  Kokkos::fence();

  // 2D Column Integration and Surface extractions (k=0 surface layer,
  // integrated sum vertically)
  Kokkos::parallel_for(
      "GocartColumnIntegratorKernel",
      Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, num_cells),
      KOKKOS_LAMBDA(int i_cell) {
        double total_col_mass = 0.0;
        double total_col_pm25 = 0.0;

        for (int i_level = 0; i_level < num_levels; ++i_level) {
          double dz = d_thick(i_cell, i_level);

          // Sum vertical column mass
          for (int i_spec = 0; i_spec < num_species; ++i_spec) {
            double mass_conc = d_state(i_cell, i_level, i_spec);
            const auto &params = d_species_params(i_spec);

            if (mass_conc > 0.0) {
              total_col_mass += mass_conc * dz;
              if (params.dry_particle_diameter <= 2.5e-6) {
                total_col_pm25 += mass_conc * dz;
              }
            }
          }
        }

        // Surface layer is bottom layer (i_level = 0)
        double surface_mass = 0.0;
        double surface_pm25 = 0.0;
        for (int i_spec = 0; i_spec < num_species; ++i_spec) {
          double mass_conc = d_state(i_cell, 0, i_spec);
          const auto &params = d_species_params(i_spec);
          if (mass_conc > 0.0) {
            surface_mass += mass_conc;
            if (params.dry_particle_diameter <= 2.5e-6) {
              surface_pm25 += mass_conc;
            }
          }
        }

        // Store 2D column-integrated values in the bottom-most level slots of
        // output
        d_diags(i_cell, 0, diagnostic_indices::SURFACE_MASS) = surface_mass;
        d_diags(i_cell, 0, diagnostic_indices::COLUMN_MASS) = total_col_mass;
        d_diags(i_cell, 0, diagnostic_indices::SURFACE_PM2_5_MASS) =
            surface_pm25;
        d_diags(i_cell, 0, diagnostic_indices::COLUMN_PM2_5_MASS) =
            total_col_pm25;
      });
  Kokkos::fence();
}

/// @brief Launch the optical-properties kernels (3D coefficients + AOT).
///
/// For each (cell, level, band) the kernel resolves every species' mass
/// extinction and single-scattering albedo through one of three modes
/// (highest precedence first):
/// - **Mode C** — GEOSmie spectral table: bilinear in RH (linear) and
///   band coordinate (linear-in-log) at the solver radius node.
/// - **Mode B** — band-integrated RH lookup table, 1D linear in RH.
/// - **Mode A** — analytical Anomalous Diffraction Theory (ADT):
///   @f$ \rho = 2 x (n-1) @f$ with size parameter @f$ x = \pi D_w /\lambda @f$;
///   @f$ Q_{ext} = 2 - (4/\rho)\sin\rho + (4/\rho^2)(1-\cos\rho) @f$,
///   with a Taylor branch for @f$ \rho < 0.01 @f$ to avoid cancellation.
///   Absorption scales scattering as @f$ Q_{sca} = Q_{ext} e^{-2 x k} @f$.
/// Wet growth is @f$ \kappa @f$-Köhler in Mode A. Backscatter uses the
/// Henyey-Greenstein phase function at @f$ \theta = \pi @f$. A second
/// kernel integrates column AOTs (Beer–Lambert, @f$ \tau = \sum_l
/// \sigma_{ext,l} \Delta z_l @f$) including fine-mode and PM2.5 cuts.
/// @param state Device solver state (params + curve pool).
/// @param num_cells Horizontal cell count.
/// @param num_levels Vertical level count.
/// @param num_bands Wavelength band count.
/// @param num_species Species axis length of @p state_ptr.
/// @param wavelengths_ptr [m] queried free-space wavelengths, length num_bands.
/// @param rh_ptr [fraction] relative humidity, column-major (cell, level).
/// @param thick_ptr [m] layer thickness, column-major (cell, level).
/// @param state_ptr [kg/m^3] aerosol mass, column-major (cell, level, species).
/// @param[out] optics_ptr Column-major (cell, level, num_bands, NUM_OPTICS);
///             slots per exaero::optical_indices (AOTs at level 0).
void run_gocart_optics(GocartSolverState *state, int num_cells, int num_levels,
                       int num_bands, int num_species,
                       const double *wavelengths_ptr, const double *rh_ptr,
                       const double *thick_ptr, const double *state_ptr,
                       double *optics_ptr) {

  using MemSpace = typename Kokkos::DefaultExecutionSpace::memory_space;

  // Wrap raw pointers directly into unmanaged column-major (LayoutLeft) default
  // space Kokkos views (zero copy!)
  auto d_wavelengths =
      Kokkos::View<const double *, Kokkos::LayoutLeft, MemSpace,
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>(wavelengths_ptr,
                                                            num_bands);
  auto d_rh = Kokkos::View<const double **, Kokkos::LayoutLeft, MemSpace,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
      rh_ptr, num_cells, num_levels);
  auto d_thick = Kokkos::View<const double **, Kokkos::LayoutLeft, MemSpace,
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
      thick_ptr, num_cells, num_levels);
  auto d_state = Kokkos::View<const double ***, Kokkos::LayoutLeft, MemSpace,
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
      state_ptr, num_cells, num_levels, num_species);
  auto d_optics = Kokkos::View<double ****, Kokkos::LayoutLeft, MemSpace,
                               Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
      optics_ptr, num_cells, num_levels, num_bands,
      optical_indices::NUM_OPTICS);

  auto d_species_params = state->d_species_params;
  auto d_curve_pool = state->d_curve_pool;

  // 1. Calculate 3D Optical Coefficients (Extinction, Scattering, Lidar
  // Backscatter, Asymmetry)
  Kokkos::parallel_for(
      "GocartOptics3D_Kernel",
      Kokkos::MDRangePolicy<Kokkos::Rank<3>>(
          {0, 0, 0}, {num_cells, num_levels, num_bands}),
      KOKKOS_LAMBDA(int i_cell, int i_level, int i_band) {
        double total_ext_coeff = 0.0;
        double total_sca_coeff = 0.0;
        double weighted_asymmetry = 0.0;
        double total_lidar_backscatter = 0.0;

        double current_rh =
            Kokkos::min(d_rh(i_cell, i_level), 0.99); // clamp Relative Humidity
        double wavelength =
            d_wavelengths(i_band); // dynamic queried wavelength in meters

        for (int i_spec = 0; i_spec < num_species; ++i_spec) {
          const auto &params = d_species_params(i_spec);
          double mass_conc = d_state(i_cell, i_level, i_spec);

          if (mass_conc <= 0.0)
            continue;

          double ext_coeff_spec = 0.0;
          double sca_coeff_spec = 0.0;
          double asm_spec = 0.0;

          double tbl_ext = 0.0, tbl_ssa = 0.0, tbl_g = 0.0;
          const bool table_hit =
              params.spec_offset >= 0 &&
              species_spectral_table(params, d_curve_pool.data(), current_rh,
                                     wavelength, tbl_ext, tbl_ssa, tbl_g);
          if (table_hit) {
            // --- MODE C: GEOSmie table spectral read ---
            // Mass extinction per dry mass [m^2/kg] interpolated in RH (linear)
            // and band coordinate (linear-in-log), at the solver radius node.
            ext_coeff_spec = mass_conc * tbl_ext;
            sca_coeff_spec = ext_coeff_spec * tbl_ssa;
            asm_spec = tbl_g;
          } else if (params.has_optics_lookup) {
            // --- MODE B: RH Lookup Table 1D Linear Interpolation ---
            // Curve read from the flat device pool; loop bound from n_rh
            // extent, never a literal. Block layout:
            // {rh,ext,ssa,asm} x n_rh.
            const double *blk = d_curve_pool.data() + params.curve_offset;
            const int n_rh = params.n_rh;
            const double *rh_axis = blk + curve_slot::RH_AXIS * n_rh;
            const double *ext = blk + curve_slot::EXT * n_rh;
            const double *ssa = blk + curve_slot::SSA * n_rh;
            const double *asm_l = blk + curve_slot::ASM * n_rh;
            int i_bin = 0;
            while (i_bin < n_rh - 2 && current_rh > rh_axis[i_bin + 1]) {
              i_bin++;
            }

            double rh_lower = rh_axis[i_bin];
            double rh_upper = rh_axis[i_bin + 1];
            double weight =
                (current_rh - rh_lower) / (rh_upper - rh_lower + 1e-15);

            // Linear Interpolation
            double mee = ext[i_bin] + weight * (ext[i_bin + 1] - ext[i_bin]);
            double ssa_v = ssa[i_bin] + weight * (ssa[i_bin + 1] - ssa[i_bin]);
            asm_spec =
                asm_l[i_bin] + weight * (asm_l[i_bin + 1] - asm_l[i_bin]);

            ext_coeff_spec = mass_conc * mee * 1000.0;
            sca_coeff_spec = ext_coeff_spec * ssa_v;

          } else {
            // --- MODE A: Anomalous Diffraction Theory (ADT) Analytical Solver
            // ---
            // 1. Wet size growth
            double wet_diameter =
                params.dry_particle_diameter *
                Kokkos::pow(1.0 + params.hygroscopicity *
                                      (current_rh / (1.0 - current_rh)),
                            1.0 / 3.0);

            // 2. Size parameter x based on dynamic wavelength
            double x = (M_PI * wet_diameter) / wavelength;

            // 3. ADT Phase shift parameter rho
            double rho = 2.0 * x * (params.refractive_index_real - 1.0);
            rho = Kokkos::max(rho, 1e-12); // prevent division by zero

            // 4. Extinction efficiency Q_ext (Anomalous Diffraction Theory)
            double q_ext = 0.0;
            if (rho < 0.01) {
              // --- stable Taylor series expansion to eliminate small-particle
              // floating-point cancellation (Hole 3) ---
              q_ext = 0.5 * rho * rho - (4.0 / 45.0) * Kokkos::pow(rho, 4) +
                      (1.0 / 72.0) * Kokkos::pow(rho, 6);
            } else {
              // --- Standard ADT formula ---
              q_ext = 2.0 - (4.0 / rho) * Kokkos::sin(rho) +
                      (4.0 / (rho * rho)) * (1.0 - Kokkos::cos(rho));
            }

            // 5. Scattering efficiency Q_sca (scaled based on imaginary index
            // absorption)
            double q_sca =
                q_ext * Kokkos::exp(-2.0 * x * params.refractive_index_imag);

            // 6. Number Concentration
            double ln_sig = Kokkos::log(params.lognormal_sigma);
            double vol_factor = (M_PI / 6.0) * params.dry_density *
                                Kokkos::pow(params.lognormal_dg, 3) *
                                Kokkos::exp(4.5 * ln_sig * ln_sig);
            double num_conc = mass_conc / vol_factor;

            // 7. Coefficients calculations
            double cross_section = (M_PI / 4.0) * wet_diameter * wet_diameter;
            ext_coeff_spec = num_conc * cross_section * q_ext;
            sca_coeff_spec = num_conc * cross_section * q_sca;

            // Parameterize asymmetry as a smooth function of wet size parameter
            asm_spec =
                0.7 * (x / (x + 1.0)); // standard asymptotic growth curve
          }

          // Henyey-Greenstein (HG) backscatter phase function evaluated at 180
          // degrees (theta = pi, cos_theta = -1): P_HG(pi, g) = (1 - g^2) / (1
          // + g^2 + 2g)^1.5 = (1 - g) / (1 + g)^2
          double phase_hg_backscatter =
              (1.0 - asm_spec) / (4.0 * M_PI * Kokkos::pow(1.0 + asm_spec, 2));
          double backscatter_spec =
              sca_coeff_spec * phase_hg_backscatter; // [m⁻¹ sr⁻¹]

          total_ext_coeff += ext_coeff_spec;
          total_sca_coeff += sca_coeff_spec;
          weighted_asymmetry += sca_coeff_spec * asm_spec;
          total_lidar_backscatter += backscatter_spec;
        }

        d_optics(i_cell, i_level, i_band, optical_indices::EXTINCTION_COEFF) =
            total_ext_coeff;
        d_optics(i_cell, i_level, i_band, optical_indices::SCATTERING_COEFF) =
            total_sca_coeff;
        d_optics(i_cell, i_level, i_band, optical_indices::BACKSCATTER_COEFF) =
            total_ext_coeff * 0.05; // standard backscatter ratio
        d_optics(i_cell, i_level, i_band, optical_indices::ASYMMETRY_FACTOR) =
            total_sca_coeff > 0.0 ? (weighted_asymmetry / total_sca_coeff)
                                  : 0.0;
        d_optics(i_cell, i_level, i_band, optical_indices::LIDAR_BACKSCATTER) =
            total_lidar_backscatter;
      });
  Kokkos::fence();

  // 2. Calculate 2D Column-Integrated Optical AOTs
  Kokkos::parallel_for(
      "GocartOpticsAOT_Kernel",
      Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {num_cells, num_bands}),
      KOKKOS_LAMBDA(int i_cell, int i_band) {
        double total_ext_aot = 0.0;
        double total_sca_aot = 0.0;
        double finemode_ext_aot = 0.0;
        double finemode_sca_aot = 0.0;
        double pm25_ext_aot = 0.0;
        double pm25_sca_aot = 0.0;

        double wavelength = d_wavelengths(i_band);

        for (int i_level = 0; i_level < num_levels; ++i_level) {
          double dz = d_thick(i_cell, i_level);
          double current_rh = Kokkos::min(d_rh(i_cell, i_level), 0.99);

          for (int i_spec = 0; i_spec < num_species; ++i_spec) {
            double mass_conc = d_state(i_cell, i_level, i_spec);
            if (mass_conc <= 0.0)
              continue;

            const auto &params = d_species_params(i_spec);
            double ext_coeff_spec = 0.0;
            double sca_coeff_spec = 0.0;

            double tbl_ext = 0.0, tbl_ssa = 0.0, tbl_g = 0.0;
            const bool table_hit =
                params.spec_offset >= 0 &&
                species_spectral_table(params, d_curve_pool.data(), current_rh,
                                       wavelength, tbl_ext, tbl_ssa, tbl_g);
            if (table_hit) {
              // --- MODE C: GEOSmie table spectral read ---
              ext_coeff_spec = mass_conc * tbl_ext;
              sca_coeff_spec = ext_coeff_spec * tbl_ssa;
            } else if (params.has_optics_lookup) {
              const double *blk = d_curve_pool.data() + params.curve_offset;
              const int n_rh = params.n_rh;
              const double *rh_axis = blk + curve_slot::RH_AXIS * n_rh;
              const double *ext = blk + curve_slot::EXT * n_rh;
              const double *ssa = blk + curve_slot::SSA * n_rh;
              int i_bin = 0;
              while (i_bin < n_rh - 2 && current_rh > rh_axis[i_bin + 1]) {
                i_bin++;
              }
              double weight = (current_rh - rh_axis[i_bin]) /
                              (rh_axis[i_bin + 1] - rh_axis[i_bin] + 1e-15);
              double mee = ext[i_bin] + weight * (ext[i_bin + 1] - ext[i_bin]);
              double ssa_v =
                  ssa[i_bin] + weight * (ssa[i_bin + 1] - ssa[i_bin]);
              ext_coeff_spec = mass_conc * mee * 1000.0;
              sca_coeff_spec = ext_coeff_spec * ssa_v;
            } else {
              // ADT values
              double wet_diameter =
                  params.dry_particle_diameter *
                  Kokkos::pow(1.0 + params.hygroscopicity *
                                        (current_rh / (1.0 - current_rh)),
                              1.0 / 3.0);
              double x = (M_PI * wet_diameter) / wavelength;
              double rho = Kokkos::max(
                  2.0 * x * (params.refractive_index_real - 1.0), 1e-12);

              double q_ext = 0.0;
              if (rho < 0.01) {
                // --- stable Taylor series expansion to eliminate
                // small-particle floating-point cancellation (Hole 3) ---
                q_ext = 0.5 * rho * rho - (4.0 / 45.0) * Kokkos::pow(rho, 4) +
                        (1.0 / 72.0) * Kokkos::pow(rho, 6);
              } else {
                // --- Standard ADT formula ---
                q_ext = 2.0 - (4.0 / rho) * Kokkos::sin(rho) +
                        (4.0 / (rho * rho)) * (1.0 - Kokkos::cos(rho));
              }

              double q_sca =
                  q_ext * Kokkos::exp(-2.0 * x * params.refractive_index_imag);
              double ln_sig = Kokkos::log(params.lognormal_sigma);
              double vol_factor = (M_PI / 6.0) * params.dry_density *
                                  Kokkos::pow(params.lognormal_dg, 3) *
                                  Kokkos::exp(4.5 * ln_sig * ln_sig);
              double num_conc = mass_conc / vol_factor;
              double cross_section = (M_PI / 4.0) * wet_diameter * wet_diameter;
              ext_coeff_spec = num_conc * cross_section * q_ext;
              sca_coeff_spec = num_conc * cross_section * q_sca;
            }

            total_ext_aot += ext_coeff_spec * dz;
            total_sca_aot += sca_coeff_spec * dz;

            // Fine mode
            if (params.dry_particle_diameter <= 1.0e-6) {
              finemode_ext_aot += ext_coeff_spec * dz;
              finemode_sca_aot += sca_coeff_spec * dz;
            }

            // PM2.5
            if (params.dry_particle_diameter <= 2.5e-6) {
              pm25_ext_aot += ext_coeff_spec * dz;
              pm25_sca_aot += sca_coeff_spec * dz;
            }
          }
        }

        // Store 2D column AOT values in the bottom-most level slots of output
        d_optics(i_cell, 0, i_band, optical_indices::EXTINCTION_AOT) =
            total_ext_aot;
        d_optics(i_cell, 0, i_band, optical_indices::SCATTERING_AOT) =
            total_sca_aot;
        d_optics(i_cell, 0, i_band, optical_indices::FINE_MODE_EXTINCTION_AOT) =
            finemode_ext_aot;
        d_optics(i_cell, 0, i_band, optical_indices::FINE_MODE_SCATTERING_AOT) =
            finemode_sca_aot;
        d_optics(i_cell, 0, i_band, optical_indices::PM2_5_EXTINCTION_AOT) =
            pm25_ext_aot;
        d_optics(i_cell, 0, i_band, optical_indices::PM2_5_SCATTERING_AOT) =
            pm25_sca_aot;
      });
  Kokkos::fence();
}

/// @brief Launch the Köhler CCN-activation spectrum kernel.
///
/// For each (cell, level, supersaturation) sample, sums the number
/// concentration of every species whose critical supersaturation is
/// exceeded. From Köhler theory the analytical critical supersaturation
/// for a lognormal-bulk species is
/// @f$ S_c = \sqrt{ 4 A^3 / (27 \kappa D_{dry}^3) } @f$ with the Kelvin
/// parameter @f$ A = 2 \sigma M_w / (\rho_w R T) \approx 1.2\times10^{-9}/T @f$
/// [m]; the @f$ 10^{-30} @f$ guard prevents division by zero for
/// non-hygroscopic (@f$ \kappa = 0 @f$) species, which then never
/// activate. A species activates at the queried S when @f$ S \ge S_c @f$.
/// @param state Device solver state (species parameters).
/// @param num_cells Horizontal cell count.
/// @param num_levels Vertical level count.
/// @param num_ss Supersaturation sample count.
/// @param num_species Species axis length of @p state_ptr.
/// @param ss_ptr [fraction] queried supersaturations, length num_ss.
/// @param temp_ptr [K] temperature, column-major (cell, level).
/// @param rh_ptr [fraction] relative humidity, column-major (cell, level).
/// @param state_ptr [kg/m^3] aerosol mass, column-major (cell, level, species).
/// @param[out] ccn_ptr Activated CCN number concentration [particles/m^3],
///             column-major (cell, level, num_ss).
void run_gocart_ccn(GocartSolverState *state, int num_cells, int num_levels,
                    int num_ss, int num_species, const double *ss_ptr,
                    const double *temp_ptr, const double *rh_ptr,
                    const double *state_ptr, double *ccn_ptr) {

  using MemSpace = typename Kokkos::DefaultExecutionSpace::memory_space;

  auto d_ss =
      Kokkos::View<const double *, Kokkos::LayoutLeft, MemSpace,
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>(ss_ptr, num_ss);
  auto d_temp = Kokkos::View<const double **, Kokkos::LayoutLeft, MemSpace,
                             Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
      temp_ptr, num_cells, num_levels);
  auto d_rh = Kokkos::View<const double **, Kokkos::LayoutLeft, MemSpace,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
      rh_ptr, num_cells, num_levels);
  auto d_state = Kokkos::View<const double ***, Kokkos::LayoutLeft, MemSpace,
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
      state_ptr, num_cells, num_levels, num_species);
  auto d_ccn = Kokkos::View<double ***, Kokkos::LayoutLeft, MemSpace,
                            Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
      ccn_ptr, num_cells, num_levels, num_ss);

  auto d_species_params = state->d_species_params;

  // Parallel Cloud CCN Activation Spectrum Solver on the GPU
  Kokkos::parallel_for(
      "GocartCcnSpectrumKernel",
      Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0},
                                             {num_cells, num_levels, num_ss}),
      KOKKOS_LAMBDA(int i_cell, int i_level, int i_ss) {
        double total_activated_ccn = 0.0;

        double temp = d_temp(i_cell, i_level);
        double query_ss =
            d_ss(i_ss); // supersaturation (fraction, e.g. 0.005 for 0.5% SS)

        // Kelvin parameter: A_kelvin = 2 * sigma * M_w / (rho_water * R * T)
        // ≈ 1.2e-9 / T [m]
        double a_kelvin = 1.2e-9 / temp;

        for (int i_spec = 0; i_spec < num_species; ++i_spec) {
          const auto &params = d_species_params(i_spec);
          double mass_conc = d_state(i_cell, i_level, i_spec);

          if (mass_conc <= 0.0)
            continue;

          // Calculate analytical critical supersaturation S_c for liquid cloud
          // droplet activation S_c = sqrt((4 * A^3) / (27 * kappa * D_dry^3))
          // [fraction]
          double s_crit =
              Kokkos::sqrt((4.0 * Kokkos::pow(a_kelvin, 3)) /
                           (27.0 * params.hygroscopicity *
                                Kokkos::pow(params.dry_particle_diameter, 3) +
                            1e-30));

          // If queried supersaturation exceeds the critical activation
          // threshold, the species activates!
          if (query_ss >= s_crit) {
            double ln_sig = Kokkos::log(params.lognormal_sigma);
            double vol_factor = (M_PI / 6.0) * params.dry_density *
                                Kokkos::pow(params.lognormal_dg, 3) *
                                Kokkos::exp(4.5 * ln_sig * ln_sig);
            double num_conc = mass_conc / vol_factor; // [particles/m³]
            total_activated_ccn += num_conc;
          }
        }

        d_ccn(i_cell, i_level, i_ss) = total_activated_ccn;
      });
  Kokkos::fence();
}

/// @brief Launch the raw->target emissions mapping kernel.
///
/// Per (cell, level, target species): fetches the mapped raw CECE flux
/// (NaN/Inf/negative inputs are zeroed — defensive boundary check),
/// converts an area flux @f$ [kg\,m^{-2}\,s^{-1}] @f$ to a volumetric
/// rate by dividing by @f$ \Delta z @f$ (guarded against micro-layers),
/// applies the species mass split fraction, and — in modal mode —
/// converts emitted mass to emitted number via the lognormal single-mode
/// volume @f$ v = (\pi/6) \rho_0 D_e^3 e^{4.5\ln^2\sigma_g} @f$.
/// Inactive mappings write an explicit 0 (never uninitialized memory).
/// @param state Device solver state (species parameters + mappings).
/// @param num_cells Horizontal cell count.
/// @param num_levels Vertical level count.
/// @param num_raw_species Raw CECE species axis length of @p raw_emissions_ptr.
/// @param num_target_species Package species axis length of the output.
/// @param flux_type_code 0 = [kg/m^3/s] rate, 1 = [kg/m^2/s] area flux.
/// @param thick_ptr [m] layer thickness, column-major (cell, level).
/// @param raw_emissions_ptr Column-major (cell, level, num_raw_species).
/// @param[out] target_emissions_out_ptr Column-major
///             (cell, level, num_target_species).
void run_gocart_emissions(GocartSolverState *state, int num_cells,
                          int num_levels, int num_raw_species,
                          int num_target_species, int flux_type_code,
                          const double *thick_ptr,
                          const double *raw_emissions_ptr,
                          double *target_emissions_out_ptr) {

  using MemSpace = typename Kokkos::DefaultExecutionSpace::memory_space;

  auto d_thick = Kokkos::View<const double **, Kokkos::LayoutLeft, MemSpace,
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
      thick_ptr, num_cells, num_levels);
  auto d_raw = Kokkos::View<const double ***, Kokkos::LayoutLeft, MemSpace,
                            Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
      raw_emissions_ptr, num_cells, num_levels, num_raw_species);
  auto d_out = Kokkos::View<double ***, Kokkos::LayoutLeft, MemSpace,
                            Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
      target_emissions_out_ptr, num_cells, num_levels, num_target_species);

  auto d_species_params = state->d_species_params;

  Kokkos::parallel_for(
      "GocartEmissions_Kernel",
      Kokkos::MDRangePolicy<Kokkos::Rank<3>>(
          {0, 0, 0}, {num_cells, num_levels, num_target_species}),
      KOKKOS_LAMBDA(int i_cell, int i_level, int i_spec) {
        const auto &params = d_species_params(i_spec);
        const auto &map = params.emissions_mapping;

        if (!map.is_active || map.raw_cece_index < 0 ||
            map.raw_cece_index >= num_raw_species) {
          d_out(i_cell, i_level, i_spec) = 0.0;
          return;
        }

        // 1. Fetch raw flux and clamp negative or NaN values
        double raw_flux = d_raw(i_cell, i_level, map.raw_cece_index);
        if (Kokkos::isnan(raw_flux) || !Kokkos::isfinite(raw_flux) ||
            raw_flux < 0.0) {
          raw_flux = 0.0;
        }

        // 2. Perform Unit Scaling (Area Flux -> Volumetric)
        double vol_mass_flux = 0.0;
        if (flux_type_code == 1) { // AREA_FLUX
          double dz = d_thick(i_cell, i_level);
          if (dz > 1e-12) {
            vol_mass_flux = raw_flux / dz;
          } else {
            vol_mass_flux = 0.0; // Avoid division-by-zero on micro layers
          }
        } else { // MASS_CONCENTRATION_RATE
          vol_mass_flux = raw_flux;
        }

        // 3. Apply mass split fraction
        double final_mass_rate = vol_mass_flux * map.mass_split_fraction;

        // 4. Handle Modal Number Conversion if applicable
        if (map.is_modal_mode && map.emitted_particle_diameter > 1e-12) {
          double ln_sig = Kokkos::log(map.lognormal_sigma);
          double vol_factor = (M_PI / 6.0) * params.dry_density *
                              Kokkos::pow(map.emitted_particle_diameter, 3) *
                              Kokkos::exp(4.5 * ln_sig * ln_sig);

          if (vol_factor > 1e-30) {
            d_out(i_cell, i_level, i_spec) =
                final_mass_rate / vol_factor; // Emitted Number concentration
          } else {
            d_out(i_cell, i_level, i_spec) = 0.0;
          }
        } else {
          d_out(i_cell, i_level, i_spec) =
              final_mass_rate; // Standard mass rate
        }
      });
  Kokkos::fence();
}

} // namespace exaero
