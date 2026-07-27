program test_fortran_interface
  use exaero_interface
  use iso_c_binding
  implicit none

  write(*,*) "Fortran Interface Integration Test Started..."

  ! 1. Initialize the C++/Kokkos environment
  call exaero_init_environment()
  write(*,*) "C++/Kokkos Environment Initialized."

  ! 2. Execute the test logic inside an isolated subroutine scope to trigger RAII finalization
  call run_test()
  write(*,*) "run_test Subroutine completed successfully."

  ! 3. Finalize the C++/Kokkos environment safely (all package views are already deallocated!)
  call exaero_finalize_environment()
  write(*,*) "C++/Kokkos Environment Finalized."

  write(*,*) "Fortran Interface Integration Test: PASS"

contains

  subroutine run_test()
    type(exaero_package_t) :: pkg
    character(len=512, kind=c_char) :: yaml_string
    character(len=32, kind=c_char) :: queried_name
    integer(c_int) :: resolved_idx

    ! CCPP-standard error variables
    character(len=256, kind=c_char) :: errmsg
    integer(c_int) :: errflg

    ! Mock dimensions
    integer, parameter :: num_cells = 1
    integer, parameter :: num_levels = 1
    integer, parameter :: num_species = 1
    integer, parameter :: num_bands = 2
    integer, parameter :: num_ss = 3

    ! Environmental states and aerosol state variables
    real(c_double) :: temp(num_cells, num_levels)
    real(c_double) :: pres(num_cells, num_levels)
    real(c_double) :: dens(num_cells, num_levels)
    real(c_double) :: rh(num_cells, num_levels)
    real(c_double) :: thick(num_cells, num_levels)
    real(c_double) :: mass_state(num_cells, num_levels)

    ! Outputs
    real(c_double) :: diags(num_cells, num_levels, 11) ! NUM_DIAGNOSTICS is 11
    real(c_double) :: optics(num_cells, num_levels, num_bands, 11) ! NUM_OPTICS is 11
    real(c_double) :: ccn(num_cells, num_levels, num_ss)
    real(c_double) :: raw_emit(num_cells, num_levels, 2)
    real(c_double) :: target_out(num_cells, num_levels, num_species)

    ! Dynamic queries
    real(c_double) :: wavelengths(num_bands)
    real(c_double) :: supersaturations(num_ss)

    ! Create the GOCART Package instance
    write(*,*) "Creating package instance..."
    pkg = exaero_create_gocart_package()
    if (.not. c_associated(pkg%ptr)) then
      write(*,*) "Error: Failed to create GocartPackage from Fortran!"
      call exit(1)
    end if
    write(*,*) "Package instance created successfully."

    ! Prepare YAML string - written as a single, flat, continuous line of kind=c_char with no concatenations or trims!
    yaml_string = c_char_"species: [{name: 'Dust', dry_density: 2600.0, molecular_weight: 100.0, dry_particle_diameter: 0.15e-6, hygroscopicity: 0.1, lognormal_sigma: 1.5, lognormal_dg: 0.1e-6, refractive_index_real: 1.55, refractive_index_imag: 0.002}]" // c_null_char

    ! Initialize package from Fortran
    write(*,*) "Initializing package from YAML..."
    errmsg = ""
    errflg = 0
    call exaero_initialize_package(pkg, yaml_string, errmsg, errflg)
    if (errflg /= 0) then
      write(*,*) "Error: Failed to initialize package: ", errmsg
      call exit(1)
    end if
    write(*,*) "Package initialized successfully."

    ! --- Dynamic Metadata Queries Validation (Hole 4) ---
    write(*,*) "Running dynamic metadata queries..."
    resolved_idx = exaero_get_species_index(pkg, "Dust" // char(0))
    if (resolved_idx /= 0) then
      write(*,*) "Error: Failed to dynamically map Dust species to index offset!"
      call exit(1)
    end if

    call exaero_get_species_name(pkg, 0, queried_name, 32)
    if (queried_name(1:4) /= "Dust") then
      write(*,*) "Error: Dynamically queried species name mismatch in Fortran!"
      call exit(1)
    end if
    write(*,*) "Dynamic metadata queries completed successfully."

    ! Populate mock inputs
    temp = 298.0d0
    pres = 101325.0d0
    dens = 1.2d0
    rh = 0.50d0
    thick = 100.0d0
    mass_state = 1.0d-6

    wavelengths(1) = 550.0d-9
    wavelengths(2) = 870.0d-9
    supersaturations(1) = 0.0005d0
    supersaturations(2) = 0.001d0
    supersaturations(3) = 0.005d0

    ! Clear output arrays
    diags = 0.0d0
    optics = 0.0d0
    ccn = 0.0d0

    ! Run diagnostics
    write(*,*) "Computing diagnostics..."
    errmsg = ""
    errflg = 0
    call exaero_compute_diagnostics(pkg, num_cells, num_levels, num_species, &
        temp, pres, dens, rh, thick, mass_state, diags, errmsg, errflg)
    if (errflg /= 0) then
      write(*,*) "Error: Failed to compute diagnostics: ", errmsg
      call exit(1)
    end if
    write(*,*) "Diagnostics completed successfully."

    ! Assert 3D PM2.5 and Column Mass calculations are correct in Fortran space!
    if (abs(diags(1, 1, exaero_COLUMN_MASS + 1) - 1.0d-4) > 1e-12) then
      write(*,*) "Error: Column Mass calculation mismatch in Fortran!"
      call exit(1)
    end if
    if (diags(1, 1, exaero_PM2_5_CONCENTRATION + 1) /= 1.0d-6) then
      write(*,*) "Error: PM2.5 Concentration mismatch in Fortran!"
      call exit(1)
    end if

    ! Run optics
    write(*,*) "Computing optics..."
    errmsg = ""
    errflg = 0
    call exaero_compute_optics(pkg, num_cells, num_levels, num_bands, num_species, &
        wavelengths, temp, pres, dens, rh, thick, mass_state, optics, errmsg, errflg)
    if (errflg /= 0) then
      write(*,*) "Error: Failed to compute optics: ", errmsg
      call exit(1)
    end if
    write(*,*) "Optics completed successfully."

    if (optics(1, 1, 1, exaero_EXTINCTION_COEFF + 1) <= 0.0) then
      write(*,*) "Error: Extinction coefficient is zero or negative!"
      call exit(1)
    end if

    ! Run CCN activation spectrum
    write(*,*) "Computing cloud CCN activation..."
    errmsg = ""
    errflg = 0
    call exaero_compute_ccn(pkg, num_cells, num_levels, num_ss, num_species, &
        supersaturations, temp, rh, mass_state, ccn, errmsg, errflg)
    if (errflg /= 0) then
      write(*,*) "Error: Failed to compute CCN: ", errmsg
      call exit(1)
    end if
    write(*,*) "Cloud CCN activation completed successfully."

    if (ccn(1, 1, 3) < ccn(1, 1, 2) .or. ccn(1, 1, 2) < ccn(1, 1, 1)) then
      write(*,*) "Error: Monotonicity activation spectrum failed in Fortran!"
      call exit(1)
    end if

    ! Run Emissions mapping stub
    write(*,*) "Computing emissions mapping stub..."
    raw_emit = 1.0d-6
    target_out = 0.0d0
    errmsg = ""
    errflg = 0
    call exaero_compute_emissions(pkg, num_cells, num_levels, 2, num_species, &
        0, thick, raw_emit, target_out, errmsg, errflg)
    if (errflg /= 0) then
      write(*,*) "Error: Failed to compute emissions: ", errmsg
      call exit(1)
    end if
    write(*,*) "Emissions mapping stub completed successfully."

    ! Manually free Gocart Package instance before Kokkos finalizes
    write(*,*) "Freeing package instance..."
    call exaero_free_package(pkg)
    write(*,*) "Package instance freed successfully."

  end subroutine run_test

end program test_fortran_interface
