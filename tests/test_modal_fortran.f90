!> @file test_modal_fortran.f90
!> @brief Generic modal host layout, units, and error-path interoperability.
!> @author EX-aero contributors
!> @details History: 2026-10. Usage: CTest. Input: inline YAML; no file I/O.
program test_modal_fortran
  use iso_c_binding
  use exaero_interface
  implicit none
  type(exaero_package_t) :: pkg
  integer(c_int) :: status
  character(kind=c_char) :: error(256)
  character(len=*), parameter :: yaml = &
    'backend: modal'//new_line('a')// &
    'coagulation: false'//new_line('a')// &
    'species: [{name: acid, density: 1700, molecular_weight: 0.098079, hygroscopicity: 0.507}]'//new_line('a')// &
    'modes: [{name: fine, geometric_radius: 3e-8, geometric_std_dev: 1.6}]'//new_line('a')// &
    'nucleation: {destination: fine, prefactor: 1e-12, exponent: 1, seed_radius: 1e-9, mass_fractions: [1]}'
  real(c_double), target :: jac(2,3,1,18)
  real(c_double) :: temp(2,3), rh(2,3), state(2,3,2), props(2,3,1,9)
  real(c_double) :: uptake(2,3,1,1), diff(1), mw(1), alpha(1)
  real(c_double) :: masses(1), number, loss, rate, mass_rate(1)
  integer(c_int) :: destination

  call exaero_init_environment()
  pkg = exaero_create_package('modal'//c_null_char, error, status)
  if (status /= 0) error stop 'FATAL ERROR: modal construction'
  call exaero_initialize_package(pkg, yaml//c_null_char, error, status)
  if (status /= 0) error stop 'FATAL ERROR: modal YAML'
  if (exaero_get_state_field_count(pkg) /= 2) error stop 'FATAL ERROR: modal field count'
  if (exaero_get_modal_mode_count(pkg) /= 1) error stop 'FATAL ERROR: modal mode count'
  masses = 1e-9_c_double
  call exaero_modal_emitted_number(pkg, 0_c_int, 1_c_int, masses, number, error, status)
  if (status /= 0 .or. number <= 0) error stop 'FATAL ERROR: emitted number'
  temp = 280.0_c_double; rh = 0.5_c_double
  state(:,:,1) = masses(1); state(:,:,2) = number
  call exaero_modal_properties(pkg, 2_c_int, 3_c_int, 2_c_int, 1_c_int, temp, rh, state, props, &
    c_loc(jac), error, status)
  if (status /= 0) error stop 'FATAL ERROR: modal properties'
  if (any(abs(props(:,:,1,1)-3e-8_c_double) > 1e-20_c_double)) error stop 'FATAL ERROR: dry radius units'
  if (any(props(:,:,1,2) <= props(:,:,1,1))) error stop 'FATAL ERROR: wet growth'
  if (any(jac(:,:,1,1) <= 0)) error stop 'FATAL ERROR: radius mass derivative'
  diff = 1e-5_c_double; mw = 0.098079_c_double; alpha = 1.0_c_double
  call exaero_modal_uptake(pkg, 2_c_int, 3_c_int, 2_c_int, 1_c_int, 1_c_int, temp, rh, state, &
    diff, mw, alpha, uptake, c_null_ptr, error, status)
  if (status /= 0 .or. any(uptake <= 0)) error stop 'FATAL ERROR: uptake units'
  call exaero_modal_nucleation(pkg, 1e-12_c_double, 1_c_int, destination, rate, loss, mass_rate, error, status)
  if (status /= 0 .or. destination /= 0 .or. loss /= mass_rate(1)) error stop 'FATAL ERROR: nucleation budget'
  call exaero_modal_properties(pkg, 2_c_int, 3_c_int, 2_c_int, 2_c_int, temp, rh, state, props, &
    c_null_ptr, error, status)
  if (status == 0) error stop 'FATAL ERROR: bad mode extent accepted'
  call exaero_free_package(pkg)
  call exaero_finalize_environment()
end program test_modal_fortran
