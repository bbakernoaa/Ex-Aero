program test_fortran_interface
  use exaero_interface
  use iso_c_binding
  implicit none

  type(exaero_package_t) :: pkg
  character(len=512) :: yaml_string
  character(len=32, kind=c_char) :: queried_name
  integer(c_int) :: resolved_idx

  ! Mock dimensions (declared as standard default-kind integer parameters for constant bounds)
  integer, parameter :: num_cells = 1
  integer, parameter :: num_levels = 1
  integer, parameter :: num_species = 1
  integer, parameter :: num_bands = 2
  integer, parameter :: num_ss = 3

  ! Environmental states and aerosol state variables grouped together
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

  ! Dynamic queries
  real(c_double) :: wavelengths(num_bands)
  real(c_double) :: supersaturations(num_ss)

  write(*,*) "Fortran Interface Integration Test Started..."

  ! 1. Initialize the C++/Kokkos environment
  call exaero_init_environment()

  ! 2. Create the GOCART Package instance
  pkg = exaero_create_gocart_package()
  if (.not. c_associated(pkg%ptr)) then
    write(*,*) "Error: Failed to create GocartPackage from Fortran!"
    call exit(1)
  end if

  ! 3. Prepare YAML string (using simple sequential appends to avoid Fortran string continuation bugs)
  yaml_string = "species: [{name: 'Dust', dry_density: 2600.0, "
  yaml_string = trim(yaml_string) // "molecular_weight: 100.0, dry_particle_diameter: 0.15e-6, "
  yaml_string = trim(yaml_string) // "hygroscopicity: 0.1, lognormal_sigma: 1.5, "
  yaml_string = trim(yaml_string) // "lognormal_dg: 0.1e-6, refractive_index_real: 1.55, "
  yaml_string = trim(yaml_string) // "refractive_index_imag: 0.002}]" // char(0)

  ! 4. Initialize package from Fortran
  call exaero_initialize_package(pkg, yaml_string)

  ! --- Dynamic Metadata Queries Validation (Hole 4) ---
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

  ! 5. Populate mock inputs as standard-compliant whole-array slice assignments
  temp(:, :) = 298.0d0
  pres(:, :) = 101325.0d0
  dens(:, :) = 1.2d0
  rh(:, :) = 0.50d0
  thick(:, :) = 100.0d0
  mass_state(:, :) = 1.0d-6

  wavelengths(1) = 550.0d-9
  wavelengths(2) = 870.0d-9
  supersaturations(1) = 0.0005d0
  supersaturations(2) = 0.001d0
  supersaturations(3) = 0.005d0

  ! Clear output arrays using standard slice assignments
  diags(:, :, :) = 0.0d0
  optics(:, :, :, :) = 0.0d0
  ccn(:, :, :) = 0.0d0

  ! 6. Run diagnostics
  call exaero_compute_diagnostics(pkg, num_cells, num_levels, num_species, &
      temp, pres, dens, rh, thick, mass_state, diags)

  ! Assert 3D PM2.5 and Column Mass calculations are correct in Fortran space!
  if (abs(diags(1, 1, exaero_COLUMN_MASS + 1) - 1.0d-4) > 1e-12) then
    write(*,*) "Error: Column Mass calculation mismatch in Fortran!"
    call exit(1)
  end if
  if (diags(1, 1, exaero_PM2_5_CONCENTRATION + 1) /= 1.0d-6) then
    write(*,*) "Error: PM2.5 Concentration mismatch in Fortran!"
    call exit(1)
  end if

  ! 7. Run optics
  call exaero_compute_optics(pkg, num_cells, num_levels, num_bands, num_species, &
      wavelengths, temp, pres, dens, rh, thick, mass_state, optics)

  if (optics(1, 1, 1, exaero_EXTINCTION_COEFF + 1) <= 0.0) then
    write(*,*) "Error: Extinction coefficient is zero or negative!"
    call exit(1)
  end if

  ! 8. Run CCN activation spectrum
  call exaero_compute_ccn(pkg, num_cells, num_levels, num_ss, num_species, &
      supersaturations, temp, rh, mass_state, ccn)

  if (ccn(1, 1, 3) < ccn(1, 1, 2) .or. ccn(1, 1, 2) < ccn(1, 1, 1)) then
    write(*,*) "Error: Monotonicity activation spectrum failed in Fortran!"
    call exit(1)
  end if

  ! 9. Free the Gocart Package instance
  call exaero_free_package(pkg)

  ! 10. Finalize the C++/Kokkos environment
  call exaero_finalize_environment()

  write(*,*) "Fortran Interface Integration Test: PASS"

end program test_fortran_interface
