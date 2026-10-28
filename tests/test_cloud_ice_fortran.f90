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
  character(kind=c_char,len=512) :: yaml

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

  call exaero_init_environment()
  package = exaero_create_package('gocart'//c_null_char,errmsg,errflg)
  if (errflg /= 0_c_int) error stop "GOCART construction failed"
  yaml = "species:"//c_new_line// &
      "  - name: Dust"//c_new_line// &
      "    dry_density: 2600"//c_new_line// &
      "    molecular_weight: 100"//c_new_line// &
      "    dry_particle_diameter: 1.0e-6"//c_new_line// &
      "    hygroscopicity: 0.5"//c_new_line// &
      "    lognormal_sigma: 1.5"//c_new_line// &
      "    lognormal_dg: 8.0e-7"//c_new_line// &
      "    refractive_index_real: 1.55"//c_new_line// &
      "    refractive_index_imag: 0.002"//c_new_line// &
      "    ice_categories: [dust]"//c_new_line// &
      "    ice_rules: [{category: dust, kind: all}]"//c_new_line//c_null_char
  call exaero_initialize_package(package,yaml,errmsg,errflg)
  if (errflg /= 0_c_int) error stop "GOCART YAML initialization failed"
  state = 1.0e-12_c_double
  call exaero_bulk_cloud_ice(package,1_c_int,1_c_int,1_c_int,temperature,rh, &
      state,c_null_ptr,exaero_ICE_DUST,0_c_int,ambient,cloud, &
      ambient_status,cloud_status,errmsg,errflg)
  if (errflg /= 0_c_int) error stop "bulk cloud/ice binding failed"
  if (ambient_status(1) /= exaero_CLOUD_ICE_CONFIGURED) error stop "bad status"
  if (ambient(exaero_CLOUD_ICE_ELIGIBLE_NUMBER+1) <= 0) error stop "no ice candidates"
  supersaturation = 0.02_c_double
  call exaero_bulk_activation_fractions(package,1_c_int,1_c_int,1_c_int, &
      1_c_int,temperature,rh,state,supersaturation,fractions,errmsg,errflg)
  if (errflg /= 0_c_int) error stop "bulk activation binding failed"
  if (any(fractions /= 1.0_c_double)) error stop "bad activation fractions"
  call exaero_free_package(package)
  call exaero_finalize_environment()
end program test_cloud_ice_fortran
