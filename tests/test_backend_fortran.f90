!> @file test_backend_fortran.f90
!> @brief Generic backend creation and modal metadata/compute through iso_c_binding.
program test_backend_fortran
  use iso_c_binding
  use exaero_interface
  implicit none
  type(exaero_package_t) :: pkg
  integer(c_int) :: status, count, idx, quantity, mode, cloud
  character(kind=c_char) :: errmsg(256)
  real(c_double), allocatable :: state(:,:,:)
  real(c_double) :: temp(2,3), pres(2,3), dens(2,3), rh(2,3), dz(2,3)
  if (exaero_package_available('mam4xx'//c_null_char) == 0) stop
  call exaero_init_environment()
  pkg = exaero_create_package('mam4xx'//c_null_char, errmsg, status)
  if (status /= 0 .or. .not. c_associated(pkg%ptr)) error stop 'FATAL ERROR: modal construction'
  if (exaero_get_representation(pkg) /= exaero_MODAL) error stop 'FATAL ERROR: modal representation'
  call exaero_initialize_package(pkg, '{}'//c_null_char, errmsg, status)
  if (status /= 0) error stop 'FATAL ERROR: modal initialization'
  count = exaero_get_state_field_count(pkg)
  if (count /= 64) error stop 'FATAL ERROR: modal state count'
  allocate(state(2,3,count))
  state = 0.0_c_double
  idx = exaero_get_species_index(pkg, 'SO4_accumulation_i'//c_null_char)
  if (idx < 0) error stop 'FATAL ERROR: mass field lookup'
  state(:,:,idx+1) = 1e-9_c_double
  idx = exaero_get_species_index(pkg, 'num_accumulation_i'//c_null_char)
  if (idx < 0) error stop 'FATAL ERROR: number field lookup'
  quantity = -1; mode = -1; cloud = -1
  status = exaero_get_state_field_info(pkg, idx, quantity, mode, cloud)
  if (status /= 0 .or. quantity /= exaero_STATE_NUMBER .or. mode /= 0 .or. cloud /= 0) &
      error stop 'FATAL ERROR: number metadata'
  state(:,:,idx+1) = 1e8_c_double
  temp = 280.0_c_double; pres = 90000.0_c_double; dens = 1.1_c_double
  rh = 0.0_c_double; dz = 100.0_c_double
  call exaero_execute_microphysics(pkg, 2_c_int, 3_c_int, count, temp, pres, dens, rh, dz, &
      state, 60.0_c_double, errmsg, status)
  if (status /= 0) error stop 'FATAL ERROR: modal Fortran microphysics'
  if (any(state(:,:,idx+1) <= 0.0_c_double)) error stop 'FATAL ERROR: invalid modal number'
  call exaero_free_package(pkg)
  call exaero_finalize_environment()
end program test_backend_fortran
