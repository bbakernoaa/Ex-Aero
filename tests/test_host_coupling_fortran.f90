!> @file test_host_coupling_fortran.f90
!> @brief Verify population metadata, zero-based fields and atomic C errors.
program test_host_coupling_fortran
  use iso_c_binding
  use exaero_interface
  implicit none(type,external)
  type(exaero_package_t) :: pkg
  integer(c_int) :: flag,rep,nf,np,ns,idx,kind,field,pop,species,allowed,supported
  character(c_char) :: error(256),name(32)
  character(len=*),parameter :: bulk = &
    'species: [{name: A, dry_density: 1000, molecular_weight: 100, dry_particle_diameter: 2e-7, '// &
    'hygroscopicity: 0.6, lognormal_sigma: 1.5, lognormal_dg: 2e-7, '// &
    'refractive_index_real: 1.5, refractive_index_imag: 0.002}]'
  character(len=*),parameter :: modal = &
    'backend: modal'//new_line('a')// &
    'species: [{name: A, density: 1000, molecular_weight: 0.1, hygroscopicity: 0.6}]'//new_line('a')// &
    'modes: [{name: fine, geometric_radius: 1e-7, geometric_std_dev: 1.5}]'
  integer :: m
  call exaero_init_environment()
  do m=0,1
    if(m==0) then
      pkg=exaero_create_package('gocart'//c_null_char,error,flag)
      call exaero_initialize_package(pkg,bulk//c_null_char,error,flag)
    else
      pkg=exaero_create_package('modal'//c_null_char,error,flag)
      call exaero_initialize_package(pkg,modal//c_null_char,error,flag)
    end if
    if(flag/=0) error stop 'FATAL ERROR: metadata setup'
    call exaero_host_coupling_counts(pkg,rep,nf,np,ns,error,flag)
    if(flag/=0.or.rep/=m.or.nf/=m+1.or.np/=1.or.ns/=1) error stop 'FATAL ERROR: metadata counts'
    call exaero_host_population_info(pkg,0_c_int,idx,kind,field,name,32_c_int,error,flag)
    if(flag/=0.or.idx/=0.or.kind/=m.or.field/=2*m-1) error stop 'FATAL ERROR: metadata population'
    call exaero_host_constituent_info(pkg,0_c_int,pop,species,field,allowed,name,32_c_int,error,flag)
    if(flag/=0.or.pop/=0.or.species/=0.or.field/=0.or.allowed/=1.or.name(1)/='A') &
      error stop 'FATAL ERROR: metadata constituent'
    field=77; name='!'
    call exaero_host_population_info(pkg,0_c_int,idx,kind,field,name,1_c_int,error,flag)
    if(flag==0.or.field/=77.or.any(name/='!')) error stop 'FATAL ERROR: metadata atomic name error'
    supported=77
    call exaero_supports_host_operation(pkg,999_c_int,supported,error,flag)
    if(flag==0.or.supported/=77) error stop 'FATAL ERROR: metadata atomic operation error'
    call exaero_supports_host_operation(pkg,exaero_HOST_INSPECT_NUMBER_BOUNDS,supported,error,flag)
    if(flag/=0.or.supported/=m) error stop 'FATAL ERROR: metadata operation support'
    call exaero_free_package(pkg)
  end do
  call exaero_finalize_environment()
end program
