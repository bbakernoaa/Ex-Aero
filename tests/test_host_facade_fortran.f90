!> @file test_host_facade_fortran.f90
!> @brief Verify common host slabs and error atomicity for both representations.
program test_host_facade_fortran
  use iso_c_binding
  use exaero_interface
  implicit none(type,external)
  type(exaero_package_t) :: pkg
  integer(c_int) :: flag,nf,np,ncon,rep
  integer :: backend,c,k
  integer(c_int) :: sources(1),destinations(1),nroute
  real(c_double),allocatable :: returns(:,:,:),routed(:,:,:),bins(:,:,:),tails(:,:,:)
  real(c_double) :: edges(4)
  real(c_double),allocatable :: properties(:,:,:,:)
  integer(c_int),allocatable :: availability(:,:,:,:),pstatus(:,:,:)
  integer(c_int),parameter :: nc=2,nl=3
  character(c_char) :: error(256)
  real(c_double),allocatable,target :: state(:,:,:)
  real(c_double),allocatable :: a(:,:,:,:),b(:,:,:,:),fractions(:,:,:,:)
  integer(c_int),allocatable :: ast(:,:,:),bst(:,:,:)
  real(c_double) :: temp(nc,nl),rh(nc,nl),ss(3)
  character(len=*),parameter :: bulk= &
    'species:'// &
    new_line('a')// &
    '  - {name: A, dry_density: 1000, molecular_weight: 100, dry_particle_diameter: 8e-7, '// &
    'hygroscopicity: 0.6, lognormal_sigma: 1.5, lognormal_dg: 2e-7, refractive_index_real:'// &
    ' 1.55, refractive_index_imag: 0.002}'// &
    new_line('a')// &
    '  - {name: B, dry_density: 2000, molecular_weight: 100, dry_particle_diameter: 2e-6, '// &
    'hygroscopicity: 0, lognormal_sigma: 1.5, lognormal_dg: 4e-7, refractive_index_real: 1'// &
    '.55, refractive_index_imag: 0.002}'// &
    new_line('a')
  character(len=*),parameter :: modal= &
    'backend: modal'// &
    new_line('a')// &
    'coagulation: false'// &
    new_line('a')// &
    'species: [{name: A, density: 1000, molecular_weight: 0.1, hygroscopicity: 0.6}, {name'// &
    ': B, density: 2000, molecular_weight: 0.2, hygroscopicity: 0}]'// &
    new_line('a')// &
    'modes: [{name: fine, geometric_radius: 1e-7, geometric_std_dev: 1.5}, {name: coarse, '// &
    'geometric_radius: 2e-7, geometric_std_dev: 1.5, species: [A]}]'// &
    new_line('a')

  call exaero_init_environment()
  temp=280;rh=0.5;ss=[0.0_c_double,0.001_c_double,0.01_c_double]
  do backend=0,1
    if(backend==0) then
      pkg=exaero_create_package('gocart'//c_null_char,error,flag)
      call exaero_initialize_package(pkg,bulk//c_null_char,error,flag)
    else
      pkg=exaero_create_package('modal'//c_null_char,error,flag)
      call exaero_initialize_package(pkg,modal//c_null_char,error,flag)
    end if
    if(flag/=0) error stop 'FATAL ERROR: facade initialization'
    call exaero_host_coupling_counts(pkg,rep,nf,np,ncon,error,flag)
    if(flag/=0.or.np/=2) error stop 'FATAL ERROR: facade metadata'
    allocate(state(nc,nl,nf),a(nc,nl,np,12),b(nc,nl,np,12),fractions(nc,nl,np,6))
    allocate(ast(nc,nl,np),bst(nc,nl,np))
    state=0
    do k=1,nl
      do c=1,nc
        state(c,k,1)=2e-9_c_double*(c+k)
        state(c,k,2)=3e-9_c_double*(c+k)
        if(backend==1) then
          state(c,k,3)=5e7_c_double*(c+k)
          state(c,k,4)=1e-9_c_double*(c+k)
          state(c,k,6)=1e7_c_double*(c+k)
        end if
      end do
    end do
    call exaero_cloud_ice(pkg,nc,nl,nf,np,temp,rh,state,c_null_ptr, &
        exaero_ICE_DUST,0_c_int,a,b,ast,bst,error,flag)
    if(flag/=0.or.any(bst/=exaero_CLOUD_ICE_ABSENT).or.any(b/=0)) error stop 'FATAL ERROR: facade cloud status'
    if(backend==0) then
      if(any(a(:,:,:,exaero_CLOUD_ICE_MASS+1)/=state)) error stop 'FATAL ERROR: facade bulk mass slab'
    else
      if(any(a(:,:,1,exaero_CLOUD_ICE_MASS+1)/=state(:,:,1)+state(:,:,2))) &
          error stop 'FATAL ERROR: facade modal mass slab'
      if(any(a(:,:,2,exaero_CLOUD_ICE_NUMBER+1)/=state(:,:,6))) error stop 'FATAL ERROR: facade modal number slab'
    end if
    call exaero_activation_fractions(pkg,nc,nl,nf,np,3_c_int,temp,rh,state,ss,fractions,error,flag)
    if(flag/=0.or.any(fractions<0).or.any(fractions>1)) error stop 'FATAL ERROR: facade activation'
    if(any(fractions(:,:,:,1)/=0).or.any(fractions(:,:,:,2)/=0)) error stop 'FATAL ERROR: facade zero activation'
    a=77;b=88;ast=99;bst=100
    call exaero_cloud_ice(pkg,nc,nl,nf,np+1,temp,rh,state,c_null_ptr, &
        exaero_ICE_DUST,0_c_int,a,b,ast,bst,error,flag)
    if(flag==0.or.any(a/=77).or.any(b/=88).or.any(ast/=99).or.any(bst/=100)) &
        error stop 'FATAL ERROR: facade atomic cloud error'
    allocate(properties(nc,nl,np,exaero_POP_PROPERTY_COUNT),availability(nc,nl,np,exaero_POP_PROPERTY_COUNT))
    allocate(pstatus(nc,nl,np))
    call exaero_population_properties(pkg,nc,nl,nf,np,temp,rh,state,properties,availability,pstatus,error,flag)
    if(flag/=0.or.any(pstatus/=exaero_POPULATION_POPULATED)) error stop 'FATAL ERROR: facade properties status'
    if(any(properties(:,:,:,exaero_POP_NUMBER+1)<=0)) error stop 'FATAL ERROR: facade properties number'
    if(backend==0) then
      if(any(properties(:,:,:,exaero_POP_MASS+1)/=state)) error stop 'FATAL ERROR: facade bulk property mass'
    else
      if(any(properties(:,:,1,exaero_POP_MASS+1)/=state(:,:,1)+state(:,:,2))) &
        error stop 'FATAL ERROR: facade modal property mass'
      if(any(properties(:,:,2,exaero_POP_NUMBER+1)/=state(:,:,6))) error stop 'FATAL ERROR: facade property number'
    end if
    if(any(availability(:,:,:,exaero_POP_LOWER_DRY_RADIUS+1)/=exaero_PROPERTY_UNAVAILABLE)) &
      error stop 'FATAL ERROR: facade unavailable bound'
    properties=77;availability=88;pstatus=99
    call exaero_population_properties(pkg,nc,nl,nf,np+1,temp,rh,state,properties,availability,pstatus,error,flag)
    if(flag==0.or.any(properties/=77).or.any(availability/=88).or.any(pstatus/=99)) &
      error stop 'FATAL ERROR: facade atomic property error'
    deallocate(properties,availability,pstatus)
    allocate(returns(nc,nl,nf),routed(nc,nl,nf),bins(nc,nl,3),tails(nc,nl,2))
    returns=state
    sources=0;destinations=1
    if(backend==1) then
      returns(:,:,3)=0;returns(:,:,6)=0;destinations=3
    end if
    call exaero_resuspension_routes(pkg,0_c_int,nroute,c_null_ptr,c_null_ptr,error,flag)
    if(flag/=0.or.nroute/=0) error stop 'FATAL ERROR: facade route query'
    call exaero_redistribute_resuspension(pkg,nc,nl,nf,1_c_int,sources,destinations,returns,routed,error,flag)
    if(flag/=0.or.any(routed(:,:,1)/=0)) error stop 'FATAL ERROR: facade redistribution'
    if(any(abs(sum(routed,dim=3)-sum(returns,dim=3))>1e-20_c_double)) error stop 'FATAL ERROR: facade return budget'
    edges=[1e-12_c_double,2e-7_c_double,4e-7_c_double,1e-3_c_double]
    call exaero_rebin_deposition(pkg,nc,nl,nf,0_c_int,3_c_int,state,edges,bins,tails,error,flag)
    if(flag/=0) error stop 'FATAL ERROR: facade rebin'
    if(backend==0) then
      if(any(abs(sum(bins,dim=3)+sum(tails,dim=3)-state(:,:,1))>1e-20_c_double)) &
        error stop 'FATAL ERROR: facade bulk rebin budget'
    else
      if(any(abs(sum(bins,dim=3)+sum(tails,dim=3)-state(:,:,1)-state(:,:,4))>1e-20_c_double)) &
        error stop 'FATAL ERROR: facade modal rebin budget'
    end if
    bins=77;tails=88
    call exaero_rebin_deposition(pkg,nc,nl,nf,nf,3_c_int,state,edges,bins,tails,error,flag)
    if(flag==0.or.any(bins/=77).or.any(tails/=88)) error stop 'FATAL ERROR: facade atomic rebin error'
    deallocate(returns,routed,bins,tails)
    deallocate(state,a,b,fractions,ast,bst)
    call exaero_free_package(pkg)
  end do
  call exaero_finalize_environment()
end program
