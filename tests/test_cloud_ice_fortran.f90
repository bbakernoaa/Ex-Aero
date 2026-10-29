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
  call check_bulk_layout()
  call exaero_finalize_environment()
contains
  !> @brief Check multi-cell/level/species property slabs and absent reservoir.
  subroutine check_bulk_layout()
    implicit none
    integer(c_int),parameter :: nc=2,nl=3,nf=2
    type(exaero_package_t) :: pkg
    character(kind=c_char) :: error(256)
    integer(c_int) :: flag,astatus(nc,nl,nf),cstatus(nc,nl,nf)
    real(c_double) :: mass(nc,nl,nf),original(nc,nl,nf),temp(nc,nl),humidity(nc,nl)
    real(c_double) :: aout(nc,nl,nf,12),cout(nc,nl,nf,12),number(nc,nl,nf)
    integer :: cell,level,s
    character(len=*),parameter :: config = &
      'species:'//new_line('a')// &
      '  - {name: Dust, dry_density: 2600, molecular_weight: 100, dry_particle_diameter: 1e-6, hygroscopicity: 0, '// &
      'lognormal_sigma: 1.5, lognormal_dg: 8e-7, refractive_index_real: 1.55, refractive_index_imag: 0.002, '// &
      'ice_categories: [dust], ice_rules: [{category: dust, kind: all}]}'//new_line('a')// &
      '  - {name: SU, dry_density: 1800, molecular_weight: 100, dry_particle_diameter: 1e-7, hygroscopicity: 0.5, '// &
      'lognormal_sigma: 1.5, lognormal_dg: 1e-7, refractive_index_real: 1.55, refractive_index_imag: 0.002, '// &
      'ice_categories: [sulfate]}'
    pkg=exaero_create_package('gocart'//c_null_char,error,flag)
    if(flag/=0) error stop 'FATAL ERROR: bulk layout package'
    call exaero_initialize_package(pkg,config//c_null_char,error,flag)
    if(flag/=0) error stop 'FATAL ERROR: bulk layout initialization'
    do s=1,nf
      do level=1,nl
        do cell=1,nc
          mass(cell,level,s)=1e-12_c_double*(cell+2*level)*s
        end do
      end do
    end do
    original=mass; temp=280; humidity=0
    call exaero_bulk_fixed_number(pkg,nc,nl,nf,mass,number,error,flag)
    if(flag/=0) error stop 'FATAL ERROR: bulk layout number'
    call exaero_bulk_cloud_ice(pkg,nc,nl,nf,temp,humidity,mass,c_null_ptr,exaero_ICE_DUST, &
        0_c_int,aout,cout,astatus,cstatus,error,flag)
    if(flag/=0) error stop 'FATAL ERROR: bulk cloud layout'
    if(any(aout(:,:,:,exaero_CLOUD_ICE_MASS+1)/=mass)) error stop 'FATAL ERROR: bulk mass slab'
    if(any(abs(aout(:,:,:,exaero_CLOUD_ICE_NUMBER+1)-number)>1e-12_c_double*number)) &
        error stop 'FATAL ERROR: bulk number slab units'
    if(any(abs(aout(:,:,1,exaero_CLOUD_ICE_ELIGIBLE_NUMBER+1)-number(:,:,1))> &
        1e-12_c_double*number(:,:,1))) error stop 'FATAL ERROR: bulk candidate slab'
    if(any(aout(:,:,2,exaero_CLOUD_ICE_ELIGIBLE_NUMBER+1)/=0)) error stop 'FATAL ERROR: species role ordering'
    if(any(astatus(:,:,1)/=exaero_CLOUD_ICE_CONFIGURED)) error stop 'FATAL ERROR: configured slab'
    if(any(astatus(:,:,2)/=exaero_CLOUD_ICE_NOT_CONFIGURED)) error stop 'FATAL ERROR: missing rule slab'
    if(any(cstatus/=exaero_CLOUD_ICE_ABSENT).or.any(cout/=0)) error stop 'FATAL ERROR: absent cloud slab'
    if(any(mass/=original)) error stop 'FATAL ERROR: cloud diagnostics changed input'
    call exaero_free_package(pkg)
  end subroutine check_bulk_layout
end program test_cloud_ice_fortran
