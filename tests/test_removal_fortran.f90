!> @file test_removal_fortran.f90
!> @brief Physical sink signs, species/moment mapping and atomic ISO C errors.
!> @details Inline bulk/modal configurations; no I/O or accepted state update.
program test_removal_fortran
  use iso_c_binding
  use exaero_interface
  implicit none (type, external)
  call exaero_init_environment()
  call check_bulk()
  call check_modal()
  call exaero_finalize_environment()
contains
  !> @brief Compare rates relative to their own magnitude, including tiny mass sinks.
  pure elemental logical function near(a,b)
    implicit none
    real(c_double),intent(in) :: a,b
    near = abs(a-b) <= 1e-12_c_double*max(1e-300_c_double,abs(a),abs(b))
  end function near

  !> @brief Verify separate bulk species/bin coefficients and diagnostic number.
  subroutine check_bulk()
    implicit none (type, external)
    type(exaero_package_t) :: pkg
    integer(c_int),parameter :: nc=2,nl=3,nf=2
    integer(c_int) :: flag
    integer :: c,k,s
    character(kind=c_char) :: error(256)
    character(len=*),parameter :: yaml = &
      'species:'//new_line('a')// &
      '  - {name: Dust, dry_density: 1000, molecular_weight: 100, dry_particle_diameter: 8e-7, '// &
      'hygroscopicity: 0.5, lognormal_sigma: 1.5, lognormal_dg: 8e-7, '// &
      'refractive_index_real: 1.55, refractive_index_imag: 0.002}'//new_line('a')// &
      '  - {name: BC, dry_density: 2000, molecular_weight: 100, dry_particle_diameter: 2e-6, '// &
      'hygroscopicity: 0.0, lognormal_sigma: 1.0, lognormal_dg: 2e-6, '// &
      'refractive_index_real: 1.55, refractive_index_imag: 0.002}'
    real(c_double) :: state(nc,nl,nf),original(nc,nl,nf),rates(nc,nl,nf)
    real(c_double) :: number(nc,nl,nf),mass_sink(nc,nl,nf),number_sink(nc,nl,nf)
    real(c_double) :: velocity(nc,nl,nf),dz(nc,nl)
    pkg=exaero_create_package('gocart'//c_null_char,error,flag)
    if(flag/=0) error stop 'FATAL ERROR: bulk removal construction'
    call exaero_initialize_package(pkg,yaml//c_null_char,error,flag)
    if(flag/=0) error stop 'FATAL ERROR: bulk removal initialization'
    do s=1,nf
      do k=1,nl
        do c=1,nc
          state(c,k,s)=1e-9_c_double*(c+2*k)*s
          rates(c,k,s)=5e-4_c_double*(c+k)*s
          velocity(c,k,s)=0.01_c_double*(c+k)*s
          dz(c,k)=10.0_c_double*(c+k)
        end do
      end do
    end do
    original=state
    call exaero_bulk_fixed_number(pkg,nc,nl,nf,state,number,error,flag)
    if(flag/=0) error stop 'FATAL ERROR: bulk removal number diagnostic'
    call exaero_removal_rate_tendencies(pkg,nc,nl,nf,nf,state,rates,c_null_ptr,mass_sink,number_sink,error,flag)
    if(flag/=0) error stop 'FATAL ERROR: bulk removal rates'
    if(any(.not.near(mass_sink,-rates*state))) error stop 'FATAL ERROR: bulk mass rate layout'
    if(any(.not.near(number_sink,-rates*number))) error stop 'FATAL ERROR: bulk number rate layout'
    call exaero_removal_velocity_tendencies(pkg,nc,nl,nf,nf,state,velocity,dz, &
        mass_sink,number_sink,error,flag)
    if(flag/=0) error stop 'FATAL ERROR: bulk removal velocities'
    do s=1,nf
      if(any(.not.near(mass_sink(:,:,s),-velocity(:,:,s)/dz*state(:,:,s)))) &
          error stop 'FATAL ERROR: bulk velocity mass'
      if(any(.not.near(number_sink(:,:,s),-velocity(:,:,s)/dz*number(:,:,s)))) &
          error stop 'FATAL ERROR: bulk velocity number'
    end do
    mass_sink=77.0_c_double; number_sink=88.0_c_double
    rates(nc,nl,nf)=-1.0_c_double
    call exaero_bulk_removal_rate_tendencies(pkg,nc,nl,nf,state,rates,mass_sink,number_sink,error,flag)
    if(flag==0.or.any(mass_sink/=77).or.any(number_sink/=88)) error stop 'FATAL ERROR: bulk atomic rate error'
    if(any(state/=original)) error stop 'FATAL ERROR: bulk removal mutated state'
    call exaero_free_package(pkg)
  end subroutine check_bulk

  !> @brief Verify common constituent-mass loss and independent mode-number loss.
  subroutine check_modal()
    implicit none (type, external)
    type(exaero_package_t) :: pkg
    integer(c_int),parameter :: nc=2,nl=3,nf=6,nm=2
    integer(c_int) :: flag
    integer :: c,k,m,s,f
    character(kind=c_char) :: error(256)
    character(len=*),parameter :: yaml = &
      'backend: modal'//new_line('a')// &
      'coagulation: false'//new_line('a')// &
      'species: [{name: A, density: 1000, molecular_weight: 0.1, hygroscopicity: 0.6}, '// &
      '{name: B, density: 2000, molecular_weight: 0.2, hygroscopicity: 0.0}]'//new_line('a')// &
      'modes: [{name: fine, geometric_radius: 1e-7, geometric_std_dev: 1.5}, '// &
      '{name: coarse, geometric_radius: 2e-7, geometric_std_dev: 1.5}]'
    real(c_double) :: state(nc,nl,nf),original(nc,nl,nf),lm(nc,nl,nm),field_rates(nc,nl,nf)
    real(c_double),target :: ln(nc,nl,nm)
    real(c_double) :: mass_sink(nc,nl,nf),number_sink(nc,nl,nm),velocity(nc,nl,nm),dz(nc,nl)
    pkg=exaero_create_package('modal'//c_null_char,error,flag)
    if(flag/=0) error stop 'FATAL ERROR: modal removal construction'
    call exaero_initialize_package(pkg,yaml//c_null_char,error,flag)
    if(flag/=0) error stop 'FATAL ERROR: modal removal initialization'
    do m=1,nm
      do k=1,nl
        do c=1,nc
          state(c,k,(m-1)*3+1)=2e-9_c_double*(c+k)*m
          state(c,k,(m-1)*3+2)=3e-9_c_double*(c+k)*m
          state(c,k,m*3)=5e7_c_double*(c+k)*m
          lm(c,k,m)=5e-4_c_double*m
          ln(c,k,m)=1e-4_c_double*(3-m)*(c+k)
          velocity(c,k,m)=0.01_c_double*m
          dz(c,k)=10.0_c_double*(c+k)
        end do
      end do
    end do
    original=state
    call exaero_modal_removal_rate_tendencies(pkg,nc,nl,nf,nm,state,lm,ln,mass_sink,number_sink,error,flag)
    if(flag/=0) error stop 'FATAL ERROR: modal removal rates'
    do m=1,nm
      do s=1,2
        f=(m-1)*3+s
        if(any(.not.near(mass_sink(:,:,f),-lm(:,:,m)*state(:,:,f)))) &
            error stop 'FATAL ERROR: modal mixed mass rate'
      end do
      if(any(mass_sink(:,:,m*3)/=0)) error stop 'FATAL ERROR: modal reserved number slot'
      if(any(.not.near(number_sink(:,:,m),-ln(:,:,m)*state(:,:,m*3)))) &
          error stop 'FATAL ERROR: modal independent number rate'
    end do
    field_rates=0
    do m=1,nm
      do s=1,2
        f=(m-1)*3+s
        field_rates(:,:,f)=lm(:,:,m)*s
      end do
    end do
    call exaero_removal_rate_tendencies(pkg,nc,nl,nf,nm,state,field_rates,c_loc(ln), &
        mass_sink,number_sink,error,flag)
    if(flag/=0) error stop 'FATAL ERROR: common modal species rates'
    do m=1,nm
      do s=1,2
        f=(m-1)*3+s
        if(any(.not.near(mass_sink(:,:,f),-field_rates(:,:,f)*state(:,:,f)))) &
          error stop 'FATAL ERROR: common modal constituent layout'
      end do
      if(any(mass_sink(:,:,m*3)/=0)) error stop 'FATAL ERROR: common modal number mass slot'
      if(any(.not.near(number_sink(:,:,m),-ln(:,:,m)*state(:,:,m*3)))) &
          error stop 'FATAL ERROR: common modal number loss'
    end do
    call exaero_removal_velocity_tendencies(pkg,nc,nl,nf,nm,state,velocity,dz, &
        mass_sink,number_sink,error,flag)
    if(flag/=0) error stop 'FATAL ERROR: modal removal velocities'
    do m=1,nm
      do s=1,2
        f=(m-1)*3+s
        if(any(.not.near(mass_sink(:,:,f),-velocity(:,:,m)/dz*state(:,:,f)))) &
            error stop 'FATAL ERROR: modal velocity mass layout'
      end do
      if(any(.not.near(number_sink(:,:,m),-velocity(:,:,m)/dz*state(:,:,m*3)))) &
          error stop 'FATAL ERROR: modal velocity number layout'
    end do
    mass_sink=77.0_c_double; number_sink=88.0_c_double; dz(nc,nl)=0
    call exaero_modal_removal_velocity_tendencies(pkg,nc,nl,nf,nm,state,velocity,dz, &
        mass_sink,number_sink,error,flag)
    if(flag==0.or.any(mass_sink/=77).or.any(number_sink/=88)) error stop 'FATAL ERROR: modal atomic depth error'
    call exaero_modal_removal_rate_tendencies(pkg,nc,nl,nf,1_c_int,state,lm,ln,mass_sink,number_sink,error,flag)
    if(flag==0.or.any(mass_sink/=77).or.any(number_sink/=88)) error stop 'FATAL ERROR: modal atomic extent error'
    mass_sink=77;number_sink=88
    call exaero_removal_rate_tendencies(pkg,nc,nl,nf,nm,state,field_rates,c_null_ptr, &
        mass_sink,number_sink,error,flag)
    if(flag==0.or.any(mass_sink/=77).or.any(number_sink/=88)) &
      error stop 'FATAL ERROR: common modal missing number rate'
    call exaero_removal_rate_tendencies(pkg,nc,nl,nf,nm+1,state,field_rates,c_loc(ln), &
        mass_sink,number_sink,error,flag)
    if(flag==0.or.any(mass_sink/=77).or.any(number_sink/=88)) &
      error stop 'FATAL ERROR: common modal population count'
    if(any(state/=original)) error stop 'FATAL ERROR: modal removal mutated state'
    call exaero_free_package(pkg)
  end subroutine check_modal
end program test_removal_fortran
