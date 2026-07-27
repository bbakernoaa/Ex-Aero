module exaero_interface
  use iso_c_binding
  implicit none

  ! Public indices corresponding exactly to exaero/AerosolIndices.hpp
  integer(c_int), parameter :: exaero_PM2_5_CONCENTRATION = 1
  integer(c_int), parameter :: exaero_COLUMN_MASS = 8
  integer(c_int), parameter :: exaero_EXTINCTION_COEFF = 0
  integer(c_int), parameter :: exaero_EXTINCTION_AOT = 5

  ! C-interoperable opaque derived type representing the EX-aero package instance
  type, bind(c) :: exaero_package_t
    type(c_ptr) :: ptr = c_null_ptr
  end type exaero_package_t

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

    subroutine exaero_init_environment() bind(c, name="exaero_init_environment")
    end subroutine exaero_init_environment

    subroutine exaero_finalize_environment() bind(c, name="exaero_finalize_environment")
    end subroutine exaero_finalize_environment
  end interface

end module exaero_interface
