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

    ! GEOSmie MIE attribute round-trip (T022)
    real(c_double) :: mie_value
    character(len=32, kind=c_char) :: mie_unit
    character(len=128, kind=c_char) :: mie_version
    integer(c_int) :: mie_status

    ! GEOSmie polarized-moment round-trip (T034, FR-003)
    type(exaero_package_t) :: mpkg
    character(len=4096, kind=c_char) :: mom_yaml
    character(len=4096) :: data_dir
    integer :: env_status
    real(c_double) :: mom_value
    integer(c_int) :: mom_status, n_pol, n_mom, rc
    real(c_double) :: mom_attrs(1, 1, 18)
    integer(c_int) :: mom_stat(1, 1, 18)
    real(c_double) :: mom_wl(1)

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
    yaml_string = c_char_"species: [{name: 'Dust', dry_density: 2600.0, molecular_weight: 100.0, dry_particle_diameter: 0.15e-6, hygroscopicity: 0.1, lognormal_sigma: 1.5, lognormal_dg: 0.1e-6, refractive_index_real: 1.55, refractive_index_imag: 0.002, mie_table: {source: DU}}]" // c_null_char

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

    ! --- GEOSmie MIE microphysical round-trip (T022): query effective radius of the
    !     DU-bound Dust species at the dry grid point through the C boundary. ---
    write(*,*) "Querying MIE microphysical attribute..."
    errmsg = ""
    errflg = 0
    mie_value = -1.0d0
    call exaero_query_attribute(pkg, 0, exaero_CAT_MICROPHYSICAL, exaero_EFFECTIVE_RADIUS, &
        0.0d0, 0.0d0, mie_value, mie_unit, 32, mie_version, 128, mie_status, errmsg, errflg)
    if (errflg /= 0) then
      write(*,*) "Error: MIE query failed: ", errmsg
      call exit(1)
    end if
    if (mie_status /= exaero_STATUS_AVAILABLE) then
      write(*,*) "Error: MIE status not Available:", mie_status
      call exit(1)
    end if
    ! DU bin0 dry effective radius = 6.358845325848961e-07 m (pinned snapshot).
    if (abs(mie_value - 6.358845325848961d-07) > 1d-13) then
      write(*,*) "Error: MIE effective radius mismatch:", mie_value
      call exit(1)
    end if
    if (mie_unit(1:1) /= 'm') then
      write(*,*) "Error: MIE unit not meters:", mie_unit
      call exit(1)
    end if
    write(*,*) "MIE microphysical round-trip completed successfully."

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
    ! --- GEOSmie polarized-moment round-trip (T034, FR-003/C1): a second package loads
    !     the DU monochromatic runtime file, which carries the rank-5 pmom array over
    !     (radius, rh, lambda, pol, moment) with the documented element ordering
    !     P11,P12,P33,P34,P22,P44. Counts are curve DATA (6 x 3 here), never literals. ---
    write(*,*) "Running polarized-moment round-trip..."
    call get_environment_variable("EXAERO_TEST_DATA_DIR", data_dir, status=env_status)
    if (env_status /= 0) then
      write(*,*) "Error: EXAERO_TEST_DATA_DIR not set for moments fixture"
      call exit(1)
    end if

    mom_yaml = c_char_"{species: [{name: 'Dust', dry_density: 2600.0, molecular_weight: 100.0, dry_particle_diameter: 0.15e-6, hygroscopicity: 0.1, lognormal_sigma: 1.5, lognormal_dg: 0.1e-6, refractive_index_real: 1.55, refractive_index_imag: 0.002, mie_table: {source: DU}}, {name: 'SeaSalt', dry_density: 2600.0, molecular_weight: 100.0, dry_particle_diameter: 0.15e-6, hygroscopicity: 0.1, lognormal_sigma: 1.5, lognormal_dg: 0.1e-6, refractive_index_real: 1.55, refractive_index_imag: 0.002}], activation: {categories: [microphysical, spectral, polarized], data_file: '" &
        // trim(data_dir) // "/mie_dust_monochromatic.txt'}}" // c_null_char

    mpkg = exaero_create_gocart_package()
    if (.not. c_associated(mpkg%ptr)) then
      write(*,*) "Error: Failed to create moments package!"
      call exit(1)
    end if
    errmsg = ""
    errflg = 0
    call exaero_initialize_package(mpkg, mom_yaml, errmsg, errflg)
    if (errflg /= 0) then
      write(*,*) "Error: moments package init failed: ", errmsg
      call exit(1)
    end if

    ! (1) Counts are data: the fixture declares 6 elements x 3 moments (FR-003).
    errmsg = ""
    errflg = 0
    rc = exaero_get_moment_counts(mpkg, 0, n_pol, n_mom, errmsg, errflg)
    if (errflg /= 0 .or. rc /= 0) then
      write(*,*) "Error: moment counts failed: ", errmsg
      call exit(1)
    end if
    if (n_pol /= 6 .or. n_mom /= 3) then
      write(*,*) "Error: moment counts wrong:", n_pol, n_mom
      call exit(1)
    end if

    ! (2) Scalar query: element P11 (0), moment 0 at the file's central band and the
    !     rh=0.5 grid point -> DU qext[bin0,rh10,band2] = 1.9339340925216675, file-delivered.
    mom_value = -1.0d0
    call exaero_query_attribute(mpkg, 0, exaero_CAT_POLARIZED_MOMENT, &
        exaero_PHASE_FUNCTION_MOMENT, 0.5d0, 5.5d-7, mom_value, &
        mie_unit, 32, mie_version, 128, mom_status, errmsg, errflg)
    if (errflg /= 0) then
      write(*,*) "Error: moment query failed: ", errmsg
      call exit(1)
    end if
    if (mom_status /= exaero_STATUS_AVAILABLE_FILE) then
      write(*,*) "Error: moment status not AvailableFile:", mom_status
      call exit(1)
    end if
    if (abs(mom_value - 1.9339340925216675d0) > 1d-13) then
      write(*,*) "Error: P11 moment-0 value mismatch:", mom_value
      call exit(1)
    end if

    ! (3) Documented ordering: element P22 (e=4, m=0) -> qe00 * 1.4.
    mom_value = -1.0d0
    call exaero_query_attribute(mpkg, 0, exaero_CAT_POLARIZED_MOMENT, &
        0 * exaero_PMOM_ELEMENT_STRIDE + 4, 0.5d0, 5.5d-7, mom_value, &
        mie_unit, 32, mie_version, 128, mom_status, errmsg, errflg)
    if (errflg /= 0 .or. mom_status /= exaero_STATUS_AVAILABLE_FILE) then
      write(*,*) "Error: P22 query failed:", mom_status
      call exit(1)
    end if
    if (abs(mom_value - 2.7075077295303345d0) > 1d-13) then
      write(*,*) "Error: P22 moment-0 value mismatch:", mom_value
      call exit(1)
    end if

    ! (4) Moment 1 (idx = 1*STRIDE + 0) halves the value: 0.9669670462608337.
    mom_value = -1.0d0
    call exaero_query_attribute(mpkg, 0, exaero_CAT_POLARIZED_MOMENT, &
        1 * exaero_PMOM_ELEMENT_STRIDE, 0.5d0, 5.5d-7, mom_value, &
        mie_unit, 32, mie_version, 128, mom_status, errmsg, errflg)
    if (errflg /= 0 .or. mom_status /= exaero_STATUS_AVAILABLE_FILE) then
      write(*,*) "Error: moment-1 query failed:", mom_status
      call exit(1)
    end if
    if (abs(mom_value - 0.9669670462608337d0) > 1d-13) then
      write(*,*) "Error: P11 moment-1 value mismatch:", mom_value
      call exit(1)
    end if

    ! (5) Out-of-range moment (M=3, idx 18): explicit NotInSource, value untouched.
    mom_value = -777.0d0
    call exaero_query_attribute(mpkg, 0, exaero_CAT_POLARIZED_MOMENT, &
        3 * exaero_PMOM_ELEMENT_STRIDE, 0.5d0, 5.5d-7, mom_value, &
        mie_unit, 32, mie_version, 128, mom_status, errmsg, errflg)
    if (errflg /= 0 .or. mom_status /= exaero_STATUS_NOT_IN_SOURCE) then
      write(*,*) "Error: out-of-range moment status:", mom_status
      call exit(1)
    end if
    if (mom_value /= -777.0d0) then
      write(*,*) "Error: out-of-range moment clobbered value!"
      call exit(1)
    end if

    ! (6) A spherical species (SeaSalt, no moments): explicit NotInSource (FR-008).
    n_pol = -1
    n_mom = -1
    rc = exaero_get_moment_counts(mpkg, 1, n_pol, n_mom, errmsg, errflg)
    if (errflg /= 0 .or. rc /= 0 .or. n_pol /= 0 .or. n_mom /= 0) then
      write(*,*) "Error: SeaSalt counts should be 0:", n_pol, n_mom
      call exit(1)
    end if

    ! (7) Bulk path: compute_attributes fills n_pol*n_moment = 18 slots, dense
    !     row-major (element fastest); slot 5 (e=5,m=0) = qe00*1.5 = 2.900901138782501.
    mom_attrs = -777.0d0
    mom_stat = -1
    mom_wl(1) = 5.5d-7
    call exaero_compute_attributes(mpkg, 1, 1, temp, pres, dens, rh, thick, &
        mass_state, 0, exaero_CAT_POLARIZED_MOMENT, 1, mom_wl, &
        mom_attrs, mom_stat, errmsg, errflg)
    if (errflg /= 0) then
      write(*,*) "Error: moment bulk compute failed: ", errmsg
      call exit(1)
    end if
    if (mom_stat(1, 1, 1) /= exaero_STATUS_AVAILABLE_FILE) then
      write(*,*) "Error: bulk slot 0 status:", mom_stat(1, 1, 1)
      call exit(1)
    end if
    if (abs(mom_attrs(1, 1, 1) - 1.9339340925216675d0) > 1d-13) then
      write(*,*) "Error: bulk slot 0 value:", mom_attrs(1, 1, 1)
      call exit(1)
    end if
    if (abs(mom_attrs(1, 1, 6) - 2.900901138782501d0) > 1d-13) then
      write(*,*) "Error: bulk slot 5 (P44,m0) value:", mom_attrs(1, 1, 6)
      call exit(1)
    end if
    if (abs(mom_attrs(1, 1, 18) - 0.7252252846956253d0) > 1d-13) then
      write(*,*) "Error: bulk slot 17 (P44,m2) value:", mom_attrs(1, 1, 18)
      call exit(1)
    end if

    call exaero_free_package(mpkg)
    write(*,*) "Polarized-moment round-trip completed successfully."

    ! Manually free Gocart Package instance before Kokkos finalizes
    write(*,*) "Freeing package instance..."
    call exaero_free_package(pkg)
    write(*,*) "Package instance freed successfully."

  end subroutine run_test

end program test_fortran_interface
