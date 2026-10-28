!> @file test_cloud_ice_fortran.f90
!> @brief Verify cloud/ice ISO_C_BINDING argument and error protocol.
program test_cloud_ice_fortran
  use iso_c_binding
  use exaero_interface
  implicit none
  type(exaero_package_t) :: package
  real(c_double) :: temperature(1), rh(1), state(1)
  real(c_double) :: ambient(12), cloud(12), supersaturation(1), fractions(2)
  integer(c_int) :: ambient_status(1), cloud_status(1), errflg
  character(kind=c_char) :: errmsg(256)

  package%ptr = c_null_ptr
  temperature = 280.0_c_double
  rh = 0.5_c_double
  state = 1.0e-12_c_double
  ambient = -9.0_c_double
  cloud = -9.0_c_double
  ambient_status = -9
  cloud_status = -9
  call exaero_bulk_cloud_ice(package,1_c_int,1_c_int,1_c_int,temperature,rh, &
      state,c_null_ptr,exaero_ICE_DUST,0_c_int,ambient,cloud, &
      ambient_status,cloud_status,errmsg,errflg)
  if (errflg /= 1_c_int) error stop "null bulk package accepted"
  if (any(ambient /= -9.0_c_double)) error stop "failed call changed ambient"

  supersaturation = 0.001_c_double
  fractions = -9.0_c_double
  call exaero_modal_activation_fractions(package,1_c_int,1_c_int,1_c_int, &
      1_c_int,temperature,rh,state,supersaturation,fractions,errmsg,errflg)
  if (errflg /= 1_c_int) error stop "null modal package accepted"
  if (any(fractions /= -9.0_c_double)) error stop "failed call changed fractions"
end program test_cloud_ice_fortran
