# EX-aero: Exascale Aerosol Interface

**EX-aero** (**EX**ascale **AERO**sol Interface) is an independent, zero-dependency, layout-agnostic C++20 and Kokkos interface library. It is designed to decouple host atmospheric dynamical and physical cores (such as CATChem, CECE, or CCPP) from advanced aerosol microphysics and optics engines.

EX-aero serves as a high-performance, GPU-portable aerosol state and microphysics boundary. This repository contains the complete software scaffolding, a privately vendored single-header `mdspan` backport, and a fully functional, dynamic-YAML-driven **GOCART (bulk aerosol) proof of concept (PoC)** solver.

---

## 1. Core Architecture & Separation of Concerns

EX-aero implements a target-level compiler isolation model to prevent template, header guard, and namespace collisions between C++20 `std::experimental::mdspan` and `<Kokkos_Core.hpp>` (such as known adjacent template parameter pack parsing bugs on Apple Clang):

*   **Public API Layer (`src/exaero/`):** Defines clean abstract interfaces (like `IAerosolPackage`) using layout-agnostic `mdspan` views. It has absolutely **zero** dependencies on Kokkos headers, keeping host builds lightweight and clean.
*   **GPU Solver Layer (`exaero_solver` target):** Contains raw C++ and Kokkos GPU parallel kernels. It includes `<Kokkos_Core.hpp>` but has **zero** access to public `mdspan` headers, keeping compilation isolated and compile-safe.
*   **Interface Adapter (`exaero_impl` target):** Bridges the public API and Kokkos solvers. It extracts raw data pointers (`.data_handle()`) on the `mdspan` boundary and routes them zero-copy to unmanaged `Kokkos::View` structures on the Default Execution Space (GPU/OpenMP) inside the solvers.

---

## 2. Directory Layout

The repository is organized to strictly enforce compilation isolation:

```text
Ex-Aero/
├── CMakeLists.txt                  # Top-level C++20 build system
├── README.md                       # Documentation
├── docs/                           # Specifications and Charters
└── src/
    ├── CMakeLists.txt              # Shared public exaero library target
    ├── exaero/                     # PUBLIC API HEADERS (Namespace: exaero)
    │   ├── IAerosolPackage.hpp     # Abstract aerosol package interface
    │   ├── AerosolIndices.hpp      # Public diagnostic & optical enums
    │   └── Environment.hpp         # Public environment lifecycle helpers
    └── src_impl/                   # PRIVATE IMPLEMENTATIONS (Internal only)
        ├── CMakeLists.txt          # Defines exaero_impl and exaero_solver targets
        ├── gocart/                 # GOCART diagnostics & optical solvers
        │   ├── GocartPackage.hpp
        │   ├── GocartPackage.cpp   # YAML parsing & raw pointer routing
        │   ├── GocartSolver.cpp    # Kokkos GPU parallel kernels
        │   └── GocartSpeciesParams.hpp # Isolated species structures
        └── third_party/            # Privately vendored dependencies
            └── mdspan/             # Isolated single-header C++20 mdspan
```

---

## 3. Supported Physical & Optical Solvers

The GOCART proof of concept implements advanced exascale aerosol solvers executing performance-portably on GPUs:

### **A. Mass Sizing & ESM Diagnostics (Task 5)**
*   **$\kappa$-Kohler Wet Sizing:** Dynamically calculates wet aerosol particle sizes under varying relative humidities on the default GPU memory space:
    $$D_{\text{wet}} = D_{\text{dry}} \left(1.0 + \kappa \frac{RH}{1.0 - RH}\right)^{1/3}$$
*   **Derived Lognormal Number Concentration:** Automatically derives 3D lognormal particle concentrations from GOCART's prognostic bulk mass concentrations:
    $$N_i = \frac{C_i}{\frac{\pi}{6} \rho_i D_{g,i}^3 \exp\left(\frac{9}{2} \ln^2 \sigma_{g,i}\right)}$$
*   **Aerosol Liquid Water (ALW):** Evaluates wet water mass absorbed via hygroscopic wet-to-dry volume differences [kg/m³].
*   **Gravitational Settling Terminal Velocity ($v_g$):** Computes mass-weighted Stokes settling terminal velocities [m/s] for vertical downward transport.
*   **3D Size Cuts & SAD:** Calculates 3D PM2.5 and PM10 concentrations, and total Surface Area Density (SAD).
*   **Vertical Column Integrator:** Integrates 3D mass concentrations vertically over layer thicknesses $\Delta z$ [m] to return 2D Column Mass [kg/m²] and Surface Mass.

### **B. Dual-Mode Optical Solver (Task 6)**
Supports multi-band radiative transfer queries (such as RRTMG) across any arbitrary number of wavelength bands simultaneously in a single parallel GPU kernel launch. It supports:
1.  **ADT (Anomalous Diffraction Theory) Mode:** Evaluates analytical Mie scattering efficiencies ($Q_{\text{ext}}, Q_{\text{sca}}, g$) on-the-fly on the GPU default execution space based on size parameter $x = \frac{\pi D_{\text{wet}}}{\lambda}$ and species real/imaginary refractive indices ($n, k$):
    $$Q_{\text{ext}}(\rho) = 2 - \frac{4}{\rho}\sin\rho + \frac{4}{\rho^2}(1 - \cos\rho)$$
2.  **RH Lookup Table Mode:** Performs high-performance 1D linear relative humidity interpolation directly on the GPU default space based on 8 standard species-specific RH bins (0% to 99%) loaded from YAML.
3.  **Column AOT Integrator:** Vertically integrates 3D extinction coefficients over layer heights to return Column Aerosol Optical Thickness (AOT):
    $$\tau_{\text{band}} = \sum_k \beta_{\text{ext}, k}(\text{band}) \Delta z_k$$

### **C. Cloud droplet CCN Activation Spectra**
*   **Physics:** Evaluates Kohler analytical critical activation thresholds $S_c$ on-the-fly inside GPU parallel kernels across a spectrum of standard supersaturations (S_index) simultaneously:
    $$S_c = \left( \frac{4 A^3}{27 \kappa D_{\text{dry}}^3} \right)^{1/2}$$
    If $S_{\text{query}} \ge S_c$, the species activates, adding $N_i$ to the total activated cloud droplet count.

---

## 4. How to Build & Run Tests

### **A. Prerequisites**
Ensure that the following dependencies are installed and discoverable on your platform (e.g. via Homebrew or standard environment modules):
*   CMake (VERSION >= 3.18)
*   Kokkos (built with your active execution space, e.g. CUDA, HIP, or OpenMP)
*   `yaml-cpp`
*   OpenMP (if compiling for CPU targets)

### **B. Configure & Compile**

#### **macOS (Apple Clang + Homebrew libomp):**
Configure passing the Homebrew `libomp` compiler flags explicitly to enable OpenMP:
```bash
cmake -S . -B build \
  -DCMAKE_CXX_STANDARD=20 \
  -DOpenMP_CXX_FLAGS="-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include" \
  -DOpenMP_CXX_LIB_NAMES="omp" \
  -DOpenMP_omp_LIBRARY=/opt/homebrew/opt/libomp/lib/libomp.dylib

cmake --build build
```

#### **Linux (GCC + Native OpenMP):**
```bash
cmake -S . -B build -DCMAKE_CXX_STANDARD=20
cmake --build build
```

---

## 5. Running the Test Suites

EX-aero contains two distinct test runners verifying library correctness and mathematical/crash robustness:

### **A. Standard Unit & Integration Tests**
Executes deterministic physical and optical validations:
```bash
./build/tests/exaero_test_runner
```
*   **YAML Parsing Unit Test:** Confirms YAML dynamic species loading and deep-copying configurations to GPU default execution space views.
*   **Memory Mapping Zero-Copy Test:** Confirms raw pointers wrapped in C++20 `mdspan` are mapped to unmanaged Kokkos Views, integrated over vertical layer thicknesses, and written back to raw host memory in-place zero-copy.
*   **GOCART Optics Dual-Mode Test:** Asserts that 1D linear relative humidity lookups and ADT analytical Mie solvers calculate and scale correctly across multiple bands (550nm and 870nm) on the GPU.
*   **GOCART Cloud CCN Spectra Test:** Verifies that cloud droplet activation matches Kohler thresholds and is physically monotonic across standard supersaturations.

### **B. Property-Based and Fuzz/Crash Resilience Tests**
Executes extensive randomized property-based validations and extreme crash fuzzer runs:
```bash
./build/tests/exaero_fuzz_test
```
*   **Property-Based Invariants (1000 Trials):** Runs 1,000 pseudo-random trials (with a fixed, reproducible seed) asserting key physical invariants (mass conservation, wet-size monotonicity, and $0 \le \text{SSA} \le 1$).
*   **Fuzz/Crash Resilience Boundary Testing:** Throws highly malformed, boundary-pushing, and invalid inputs (such as $RH \ge 1.0$ or $RH < 0.0$, negative layer heights, subnormal densities, quiet **NaNs**, and positive/negative **Infinities**) to verify the GPU parallel kernels defensively clamp them, gracefully handle errors, and **never trigger segmentation faults, abort traps, or GPU memory page faults**.

---

## 6. Disclaimer

This project is part of NOAA-EMC Ecosystem. Code is provided on an "as is" basis, and the user assumes responsibility for its use. See LICENSE and DISCLAIMER for details.
