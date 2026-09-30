!> @file exaero_interface.f90
!> @brief Fortran iso_c_binding interface to the EX-aero C surface.
!>
!> This module is the Fortran mirror of IAerosolPackage_C.h: it declares
!> the @c bind(c) interfaces (CCPP-standard errmsg/errflg signatures) and
!> the parameter constants that must stay numerically identical to the
!> C++ *_indices headers. A host Fortran/CCPP component @c USEs this
!> module and never sees a C++ type.
!>
!> @note ABI contract: every @c bind(c) scalar that the C callee writes
!> only on some paths (e.g. value_out on an available status) is declared
!> @c intent(inout), not @c intent(out) — see the rationale at
!> exaero_query_attribute. intent(out) would let an optimizing compiler
!> assume the callee always defines the dummy and discard the caller's
!> pre-call value.
module exaero_interface
  use iso_c_binding
  implicit none

  ! Public indices corresponding exactly to exaero/AerosolIndices.hpp
  integer(c_int), parameter :: exaero_PM2_5_CONCENTRATION = 1  !< [kg/m^3] PM2.5 mass
  integer(c_int), parameter :: exaero_COLUMN_MASS = 8          !< [kg/m^2] column mass
  integer(c_int), parameter :: exaero_EXTINCTION_COEFF = 0     !< [1/m] extinction coeff
  integer(c_int), parameter :: exaero_EXTINCTION_AOT = 5       !< column AOT (550 nm)

  ! --- GEOSmie MIE attribute surface indices (exaero/AttributeQuery.hpp) ---
  ! Categories (AttributeCategory)
  integer(c_int), parameter :: exaero_CAT_MICROPHYSICAL   = 0  !< (radius, rh)
  integer(c_int), parameter :: exaero_CAT_SPECTRAL_OPTICAL = 1 !< (radius, rh, lambda)
  integer(c_int), parameter :: exaero_CAT_POLARIZED_MOMENT = 2 !< (radius, rh, lambda, pol, moment)
  ! Availability status (AttributeStatus)
  integer(c_int), parameter :: exaero_STATUS_AVAILABLE        = 0  !< grid point, baked-in
  integer(c_int), parameter :: exaero_STATUS_AVAILABLE_FILE   = 1  !< grid point, runtime file
  integer(c_int), parameter :: exaero_STATUS_INTERPOLATED     = 2  !< between grid points
  integer(c_int), parameter :: exaero_STATUS_NOT_ACTIVATED    = 3  !< deselected by operator
  integer(c_int), parameter :: exaero_STATUS_NOT_IN_SOURCE    = 4  !< species lacks category
  integer(c_int), parameter :: exaero_STATUS_AVAILABLE_CONFIG = 5  !< config override curve
  ! Delivery source (DeliverySource)
  integer(c_int), parameter :: exaero_DELIVERY_BAKED_IN    = 0  !< build-time snapshot
  integer(c_int), parameter :: exaero_DELIVERY_RUNTIME_FILE = 1 !< optional data file
  integer(c_int), parameter :: exaero_DELIVERY_CONFIG      = 2  !< runtime config (wins)
  ! Microphysical attribute index codes (microphysical_indices)
  integer(c_int), parameter :: exaero_WET_PARTICLE_DENSITY = 0  !< [kg/m^3]
  integer(c_int), parameter :: exaero_GROWTH_FACTOR        = 1  !< wet/dry radius, >= 1
  integer(c_int), parameter :: exaero_EFFECTIVE_RADIUS     = 2  !< [m]
  integer(c_int), parameter :: exaero_MASS_MEAN_RADIUS     = 3  !< [m]
  integer(c_int), parameter :: exaero_BIN_LOWER_RADIUS     = 4  !< [m]
  integer(c_int), parameter :: exaero_BIN_UPPER_RADIUS     = 5  !< [m]
  integer(c_int), parameter :: exaero_VOLUME_PER_MASS      = 6  !< [m^3/kg]
  integer(c_int), parameter :: exaero_AREA_PER_MASS        = 7  !< [m^2/kg]
  integer(c_int), parameter :: exaero_PARTICLE_MASS        = 8  !< [kg]
  integer(c_int), parameter :: exaero_NUM_MICROPHYSICAL    = 9  !< slot count (extent)
  ! Spectral optical attribute index codes (spectral_optical_indices)
  integer(c_int), parameter :: exaero_QEXT       = 0   !< Q_ext dimensionless, >= 0
  integer(c_int), parameter :: exaero_QSCA       = 1   !< Q_sca in [0, Q_ext]
  integer(c_int), parameter :: exaero_QABS       = 2   !< Q_abs dimensionless, >= 0
  integer(c_int), parameter :: exaero_BEXT       = 3   !< mass extinction [m^2/kg]
  integer(c_int), parameter :: exaero_BSCA       = 4   !< mass scattering [m^2/kg]
  integer(c_int), parameter :: exaero_BBCK       = 5   !< mass backscatter [m^2/kg/sr]
  integer(c_int), parameter :: exaero_LIDAR_RATIO = 6  !< [sr], guarded division
  integer(c_int), parameter :: exaero_ASYM_FACTOR = 7  !< g in [-1, 1]
  integer(c_int), parameter :: exaero_SSA         = 8  !< single-scat. albedo in [0, 1]
  integer(c_int), parameter :: exaero_REFREAL     = 9  !< n (wet), > 0
  integer(c_int), parameter :: exaero_REFIMAG     = 10 !< k (wet), >= 0
  integer(c_int), parameter :: exaero_NUM_SPECTRAL = 11 !< slot count (extent)
  ! Polarized moment (polarized_moment_indices): single attribute, element ordering
  ! P11,P12,P33,P34,P22,P44; moment count is per-species DATA (never hardcoded).
  integer(c_int), parameter :: exaero_PMOM_ELEMENT_STRIDE = 6  !< P11,P12,P33,P34,P22,P44
  integer(c_int), parameter :: exaero_PHASE_FUNCTION_MOMENT = 0 !< addressed by (element, moment)

  !> @brief C-interoperable opaque handle to the EX-aero package instance.
  !>
  !> Wraps the C++ @c exaero_package_t (@c void*) so the Fortran side can
  !> hold it as a typed value and pass it @c value to the bind(c) routines
  !> without ever dereferencing the pointer.
  type, bind(c) :: exaero_package_t
    type(c_ptr) :: ptr = c_null_ptr  !< Opaque pointer to the C++ object.
  end type exaero_package_t

  !> @cond
  ! The bind(c) procedure contracts are documented once on the matching
  ! extern "C" declarations in IAerosolPackage_C.h (single source of truth);
  ! hiding them here avoids duplicate doc pages in the generated reference.
  !> @brief C-binding interfaces (CCPP-standard error signatures).
  !>
  !> Each @c bind(c) procedure maps one-to-one onto the extern "C" routine
  !> of the same name in IAerosolPackage_C.h. Multidimensional arrays are
  !> assumed-shape-free @c (*) arrays in Fortran column-major order, so a
  !> host (cell, level, species) array forwards zero-copy. Refer to the C
  !> header for the full per-argument unit/index contract.
  ! Public C-bindings interfaces (CCPP-standard error signatures)
  interface
    function exaero_create_gocart_package() result(pkg) bind(c, name="exaero_create_gocart_package")
      import :: c_ptr, exaero_package_t
      type(exaero_package_t) :: pkg
    end function exaero_create_gocart_package

    subroutine exaero_free_package(pkg) bind(c, name="exaero_free_package")
      import :: exaero_package_t
      type(exaero_package_t), value :: pkg
    end subroutine exaero_free_package

    subroutine exaero_initialize_package(pkg, config_yaml, errmsg, errflg) bind(c, name="exaero_initialize_package")
      import :: exaero_package_t, c_char, c_int
      type(exaero_package_t), value :: pkg
      character(kind=c_char), intent(in) :: config_yaml(*)
      character(kind=c_char), intent(out) :: errmsg(*)
      integer(c_int), intent(out) :: errflg
    end subroutine exaero_initialize_package

    subroutine exaero_compute_diagnostics(pkg, num_cells, num_levels, num_species, &
        temp_ptr, pres_ptr, dens_ptr, rh_ptr, thick_ptr, state_ptr, diags_ptr, &
        errmsg, errflg) bind(c, name="exaero_compute_diagnostics")
      import :: exaero_package_t, c_int, c_double, c_char
      type(exaero_package_t), value :: pkg
      integer(c_int), value :: num_cells, num_levels, num_species
      real(c_double), intent(in) :: temp_ptr(*), pres_ptr(*), dens_ptr(*), rh_ptr(*), thick_ptr(*)
      real(c_double), intent(in) :: state_ptr(*)
      real(c_double), intent(out) :: diags_ptr(*)
      character(kind=c_char), intent(out) :: errmsg(*)
      integer(c_int), intent(out) :: errflg
    end subroutine exaero_compute_diagnostics

    subroutine exaero_compute_emissions(pkg, num_cells, num_levels, num_raw, num_target, &
        flux_type_code, thick_ptr, raw_emissions_ptr, target_emissions_out_ptr, &
        errmsg, errflg) bind(c, name="exaero_compute_emissions")
      import :: exaero_package_t, c_int, c_double, c_char
      type(exaero_package_t), value :: pkg
      integer(c_int), value :: num_cells, num_levels, num_raw, num_target, flux_type_code
      real(c_double), intent(in) :: thick_ptr(*), raw_emissions_ptr(*)
      real(c_double), intent(out) :: target_emissions_out_ptr(*)
      character(kind=c_char), intent(out) :: errmsg(*)
      integer(c_int), intent(out) :: errflg
    end subroutine exaero_compute_emissions

    subroutine exaero_compute_optics(pkg, num_cells, num_levels, num_bands, num_species, &
        wavelengths_ptr, temp_ptr, pres_ptr, dens_ptr, rh_ptr, thick_ptr, state_ptr, optics_ptr, &
        errmsg, errflg) bind(c, name="exaero_compute_optics")
      import :: exaero_package_t, c_int, c_double, c_char
      type(exaero_package_t), value :: pkg
      integer(c_int), value :: num_cells, num_levels, num_bands, num_species
      real(c_double), intent(in) :: wavelengths_ptr(*)
      real(c_double), intent(in) :: temp_ptr(*), pres_ptr(*), dens_ptr(*), rh_ptr(*), thick_ptr(*)
      real(c_double), intent(in) :: state_ptr(*)
      real(c_double), intent(out) :: optics_ptr(*)
      character(kind=c_char), intent(out) :: errmsg(*)
      integer(c_int), intent(out) :: errflg
    end subroutine exaero_compute_optics

    subroutine exaero_compute_ccn(pkg, num_cells, num_levels, num_ss, num_species, &
        ss_ptr, temp_ptr, rh_ptr, state_ptr, ccn_ptr, &
        errmsg, errflg) bind(c, name="exaero_compute_ccn")
      import :: exaero_package_t, c_int, c_double, c_char
      type(exaero_package_t), value :: pkg
      integer(c_int), value :: num_cells, num_levels, num_ss, num_species
      real(c_double), intent(in) :: ss_ptr(*)
      real(c_double), intent(in) :: temp_ptr(*), rh_ptr(*)
      real(c_double), intent(in) :: state_ptr(*)
      real(c_double), intent(out) :: ccn_ptr(*)
      character(kind=c_char), intent(out) :: errmsg(*)
      integer(c_int), intent(out) :: errflg
    end subroutine exaero_compute_ccn

    function exaero_get_species_index(pkg, name) result(index) bind(c, name="exaero_get_species_index")
      import :: exaero_package_t, c_char, c_int
      type(exaero_package_t), value :: pkg
      character(kind=c_char), intent(in) :: name(*)
      integer(c_int) :: index
    end function exaero_get_species_index

    subroutine exaero_get_species_name(pkg, index, name_out, max_len) bind(c, name="exaero_get_species_name")
      import :: exaero_package_t, c_char, c_int
      type(exaero_package_t), value :: pkg
      integer(c_int), value :: index
      character(kind=c_char), intent(out) :: name_out(*)
      integer(c_int), value :: max_len
    end subroutine exaero_get_species_name

    ! --- GEOSmie MIE attribute surface (CCPP-standard error signatures) ---
    ! Array-of-pointer C arguments (const char* const*, const double* const*) bind to a
    ! single type(c_ptr), value; the caller passes c_loc of a target array of c_ptr.
    ! const int* arrays bind to intent(in) arrays (no VALUE, which conflicts with dimension).
    subroutine exaero_set_attribute_activation(pkg, species_names, num_species, &
        categories_mask, runtime_file_path, errmsg, errflg) &
        bind(c, name="exaero_set_attribute_activation")
      import :: exaero_package_t, c_char, c_int, c_ptr
      type(exaero_package_t), value :: pkg
      type(c_ptr), value :: species_names
      integer(c_int), value :: num_species
      integer(c_int), value :: categories_mask
      type(c_ptr), value :: runtime_file_path
      character(kind=c_char), intent(out) :: errmsg(*)
      integer(c_int), intent(out) :: errflg
    end subroutine exaero_set_attribute_activation

    subroutine exaero_set_species_curve_config(pkg, species_names, source_labels, &
        radius_nodes, num_radius_nodes, override_values, num_override_attributes, &
        override_attribute_ids, interpolate, num_curves, errmsg, errflg) &
        bind(c, name="exaero_set_species_curve_config")
      import :: exaero_package_t, c_char, c_int, c_ptr
      type(exaero_package_t), value :: pkg
      type(c_ptr), value :: species_names
      type(c_ptr), value :: source_labels
      type(c_ptr), value :: radius_nodes
      integer(c_int), intent(in) :: num_radius_nodes(*)
      type(c_ptr), value :: override_values
      integer(c_int), intent(in) :: num_override_attributes(*)
      type(c_ptr), value :: override_attribute_ids
      type(c_ptr), value :: interpolate
      integer(c_int), value :: num_curves
      character(kind=c_char), intent(out) :: errmsg(*)
      integer(c_int), intent(out) :: errflg
    end subroutine exaero_set_species_curve_config

    subroutine exaero_compute_attributes(pkg, num_cells, num_levels, &
        temp_ptr, pres_ptr, dens_ptr, rh_ptr, thick_ptr, state_ptr, &
        species_index, category, num_bands, wavelengths_ptr, &
        attributes_out, status_out, errmsg, errflg) &
        bind(c, name="exaero_compute_attributes")
      import :: exaero_package_t, c_int, c_double, c_char
      type(exaero_package_t), value :: pkg
      integer(c_int), value :: num_cells, num_levels
      real(c_double), intent(in) :: temp_ptr(*), pres_ptr(*), dens_ptr(*), rh_ptr(*), thick_ptr(*)
      real(c_double), intent(in) :: state_ptr(*)
      integer(c_int), value :: species_index, category
      integer(c_int), value :: num_bands
      real(c_double), intent(in) :: wavelengths_ptr(*)
      real(c_double), intent(out) :: attributes_out(*)
      integer(c_int), intent(out) :: status_out(*)
      character(kind=c_char), intent(out) :: errmsg(*)
      integer(c_int), intent(out) :: errflg
    end subroutine exaero_compute_attributes

    subroutine exaero_query_attribute(pkg, species_index, category, attribute_index, &
        rh, wavelength_m, value_out, unit_out, unit_max, version_out, version_max, &
        status_out, errmsg, errflg) &
        bind(c, name="exaero_query_attribute")
      import :: exaero_package_t, c_int, c_double, c_char
      type(exaero_package_t), value :: pkg
      integer(c_int), value :: species_index, category, attribute_index
      real(c_double), value :: rh, wavelength_m
      ! intent(inout), NOT intent(out): the C side writes value_out ONLY on an
      ! available/interpolated status (never a silent 0). intent(out) would let
      ! an optimizing compiler assume the callee always defines it and discard the
      ! caller's pre-call value on a not-in-source result (caught only in Release).
      real(c_double), intent(inout) :: value_out
      character(kind=c_char), intent(out) :: unit_out(*)
      integer(c_int), value :: unit_max
      character(kind=c_char), intent(out) :: version_out(*)
      integer(c_int), value :: version_max
      integer(c_int), intent(out) :: status_out
      character(kind=c_char), intent(out) :: errmsg(*)
      integer(c_int), intent(out) :: errflg
    end subroutine exaero_query_attribute

    ! Polarized-moment count query num_slots = num_pol * num_moment for the
    ! species' curve (data, never a literal). Both are 0 for a species without moments.
    function exaero_get_moment_counts(pkg, species_index, num_pol_out, num_moment_out, &
        errmsg, errflg) result(rc) bind(c, name="exaero_get_moment_counts")
      import :: exaero_package_t, c_int, c_char
      type(exaero_package_t), value :: pkg
      integer(c_int), value :: species_index
      integer(c_int), intent(out) :: num_pol_out, num_moment_out
      character(kind=c_char), intent(out) :: errmsg(*)
      integer(c_int), intent(out) :: errflg
      integer(c_int) :: rc
    end function exaero_get_moment_counts

    subroutine exaero_init_environment() bind(c, name="exaero_init_environment")
    end subroutine exaero_init_environment

    subroutine exaero_finalize_environment() bind(c, name="exaero_finalize_environment")
    end subroutine exaero_finalize_environment
  end interface
  !> @endcond

end module exaero_interface
