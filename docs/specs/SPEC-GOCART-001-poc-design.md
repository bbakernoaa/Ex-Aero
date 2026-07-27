# Design Spec: EX-aero (Exascale Aerosol Interface) GOCART Proof of Concept

**Date:** July 24, 2026  
**Status:** Approved Draft  
**Target Architecture:** C++20, Kokkos, Private Header-Only `mdspan` Backport

---

## 1. Executive Summary & Architecture Goals

This document specifies the initial implementation design for **EX-aero**, an independent, zero-dependency, layout-agnostic C++20 interface library designed to decouple host atmospheric models from aerosol microphysics solvers. 

Our immediate goal is to establish a **GOCART (bulk aerosol) proof of concept (PoC)**. In keeping with strict Separation of Concerns:
1. **Dynamic GOCART Characterization:** GOCART is a passive bulk scheme. Unlike complex modal models (such as MAM4xx), GOCART does not execute dynamic internal microphysics solvers (like coagulation, condensation, or nucleation).
2. **Separation of Concerns:** The GOCART package's primary role in EX-aero is to calculate and prepare aerosol state diagnostics (e.g., wet/dry particle diameters, PM2.5/PM10 concentrations, surface area density (SAD), and hygroscopic properties) and hand them back to the host's actual physics and cloud microphysics schemes (e.g., CCPP Thompson cloud microphysics, RRTMG radiation).
3. **Execution Portability:** All calculations are executed performance-portably on GPUs using **Kokkos**.
4. **Zero-Copy Memory Boundary:** Host raw pointers are wrapped into layout-agnostic C++20 `std::experimental::mdspan` view abstractions, which map zero-copy to unmanaged `Kokkos::View`s inside the GPU kernels.
5. **Zero-Dependency Core:** Vendors the header-only `mdspan` backport privately to maintain a completely clean public API.

---

## 2. Directory Layout & Vendoring Strategy

The repository structure isolates public API definitions from internal GOCART implementations and third-party dependencies.

```text
/Users/barry/Documents/Ex-Aero/
├── CMakeLists.txt                  # Top-level build configuration
├── README.md
├── docs/
│   ├── Project-Charter.md
│   └── superpowers/specs/2026-07-24-exaero-gocart-poc-design.md
└── src/
    ├── CMakeLists.txt              # Primary targets: exaero (public API) and exaero_impl
    ├── exaero/                     # PUBLIC API HEADERS (Namespace: exaero)
    │   ├── IAerosolPackage.hpp     # Abstract package interface
    │   ├── EnvironmentalState.hpp  # Host environment descriptors
    │   └── CcppConstituentAdapter.hpp # CCPP 3D array pointer adapter
    └── src_impl/                   # PRIVATE IMPLEMENTATION (Internal only)
        ├── CMakeLists.txt
        ├── gocart/                 # GOCART proof-of-concept diagnostics package
        │   ├── GocartPackage.hpp
        │   └── GocartPackage.cpp
        ├── mam4xx_wrapper/         # Shim wrapping externally pre-built mam4xx
        └── third_party/            # Privately vendored dependencies
            └── mdspan/             # Header-only Kokkos/std-backport mdspan
```

---

## 3. Zero-Copy Memory Boundary & Pointer Adapters

EX-aero bridges C++ and Fortran host memory models zero-copy. On the API boundary, we wrap raw pointers into `std::experimental::mdspan` instances. Inside GPU parallel regions, we wrap these same pointers into unmanaged `Kokkos::View` objects.

### 3.1. Raw Pointer wrapping into `mdspan`
Multi-dimensional host arrays (e.g., density, concentration state) are wrapped directly into layout-agnostic `mdspan` views matching the host's native memory layout:

```cpp
namespace exaero {
    // 3D View alias: (grid_cells, vertical_levels, species)
    template <typename T, typename Layout = std::experimental::layout_right>
    using mdspan_3d = std::experimental::mdspan<
        T,
        std::experimental::extents<size_t, std::experimental::dynamic_extent, std::experimental::dynamic_extent, std::experimental::dynamic_extent>,
        Layout
    >;

    // 2D View alias: (grid_cells, vertical_levels)
    template <typename T, typename Layout = std::experimental::layout_right>
    using mdspan_2d = std::experimental::mdspan<
        T,
        std::experimental::extents<size_t, std::experimental::dynamic_extent, std::experimental::dynamic_extent>,
        Layout
    >;
}
```

---

## 4. Dynamic GOCART YAML & GPU Property Mapping

To fulfill the dynamic specification, the active GOCART species list and physical properties are configured via YAML at runtime. Species parameters are parsed at startup and mapped to a flat, GPU-cached Kokkos View.

### 4.1. Trivially Copyable Parameter Struct
To support PM2.5/PM10 calculations, **derived number concentration**, and our **Dual-Mode Optical Solver (Anomalous Diffraction Theory & RH Lookup Tables)**, our species parameter structure includes standard lognormal distribution variables, refractive index parameters, and fixed-size RH lookup arrays:

```cpp
namespace exaero {
    struct GocartSpeciesParams {
        double dry_density;             // [kg/m³]
        double molecular_weight;        // [g/mol]
        double dry_particle_diameter;   // [m]
        double hygroscopicity;          // kappa value for wet sizing calculations
        
        // Assumed lognormal parameters for PM and Number Concentration derivations
        double lognormal_sigma;         // Geometric standard deviation (GSD)
        double lognormal_dg;            // Geometric mean diameter (GMD) [m]
        
        // Refractive indices for analytical ADT (Anomalous Diffraction Theory) optical calculations
        double refractive_index_real;   // Real part of the refractive index (n)
        double refractive_index_imag;   // Imaginary part of the refractive index (k)
        
        // Tabulated RH Lookup Table Parameters (Optional, enabled per species)
        bool has_optics_lookup;         // Flag enabling lookups instead of ADT
        double rh_bins[8];              // Standard relative humidity bins (e.g. 0.0, 0.5, 0.7, 0.8, 0.9, 0.95, 0.98, 0.99)
        double ext_lookup[8];           // Pre-tabulated Mass Extinction Efficiency [m²/g] at each RH bin
        double ssa_lookup[8];           // Pre-tabulated Single Scattering Albedo [fraction] at each RH bin
        double asm_lookup[8];           // Pre-tabulated Asymmetry Parameter [g] at each RH bin
    };
}
```

### 4.2. Runtime Initialization Procedure
During `GocartPackage::initialize(const std::string& config_yaml)`:
1.  **YAML Parsing:** Parse active species from the YAML file.
2.  **State Sizing:** Cache `num_species = species_list.size()`.
3.  **View Allocation & Copy:** Allocate a 1D default execution space `Kokkos::View<GocartSpeciesParams*, MemorySpace>` and copy the parsed parameters from a host mirrors to the device.

---

## 5. Core Public API Definitions

The public-facing `IAerosolPackage` keeps a strict separation of concerns, offering separate methods for microphysics evolution, diagnostics computations, and optical computations.

### 5.1. `src/exaero/IAerosolPackage.hpp`
```cpp
#pragma once
#include <string>

// Force backport to compile under exaero_mdspan namespace to prevent redefinition conflicts with system Kokkos
#define MDSPAN_IMPL_STANDARD_NAMESPACE exaero_mdspan
#include <exaero_experimental/mdspan>

namespace exaero {

    template <typename T>
    using View3D = exaero_mdspan::mdspan<T, exaero_mdspan::extents<size_t, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent>>;

    template <typename T>
    using View2D = exaero_mdspan::mdspan<T, exaero_mdspan::extents<size_t, std::dynamic_extent, std::dynamic_extent>>;

    struct EnvironmentalStateView {
        View2D<const double> temperature;       // [K] (cell, level)
        View2D<const double> pressure;          // [Pa] (cell, level)
        View2D<const double> air_density;       // [kg/m³] (cell, level)
        View2D<const double> relative_humidity;  // [fraction 0-1] (cell, level)
        View2D<const double> layer_thickness;   // [m] grid layer vertical thickness Delta_z (cell, level)
    };

    namespace diagnostic_indices {
        // 3D Mass and Size Diagnostics (cell, level, index)
        constexpr int MASS_CONCENTRATION = 0;       // 3D mass concentration per species [kg/m³]
        constexpr int PM2_5_CONCENTRATION = 1;      // 3D PM2.5 mass concentration [kg/m³]
        constexpr int PM10_CONCENTRATION = 2;       // 3D PM10 mass concentration [kg/m³]
        constexpr int NUMBER_CONCENTRATION = 3;     // Derived number concentration [particles/m³]
        constexpr int SURFACE_AREA_DENSITY = 4;     // SAD [m²/m³]
        
        // 2D Column-Integrated and Surface Mass Diagnostics (cell, level=0, index)
        constexpr int SURFACE_MASS = 5;             // Surface layer mass concentration [kg/m³]
        constexpr int COLUMN_MASS = 6;              // Column-integrated mass density [kg/m²]
        constexpr int SURFACE_PM2_5_MASS = 7;       // Surface layer PM2.5 mass concentration [kg/m³]
        constexpr int COLUMN_PM2_5_MASS = 8;        // Column-integrated PM2.5 mass density [kg/m²]
        
        constexpr int NUM_DIAGNOSTICS = 9;
    }

    namespace optical_indices {
        // 3D Optical Properties (cell, level, index)
        constexpr int EXTINCTION_COEFF = 0;         // 3D extinction coefficient [m⁻¹]
        constexpr int SCATTERING_COEFF = 1;         // 3D scattering coefficient [m⁻¹]
        constexpr int BACKSCATTER_COEFF = 2;        // 3D backscatter coefficient [m⁻¹ sr⁻¹]
        constexpr int ASYMMETRY_FACTOR = 3;         // g asymmetry parameter [fraction -1 to 1]
        
        // 2D Column-Integrated Optical Diagnostics (Column AOT / AOD at 550nm) (cell, level=0, index)
        constexpr int EXTINCTION_AOT = 4;           // Total Extinction Aerosol Optical Thickness (AOT)
        constexpr int SCATTERING_AOT = 5;           // Total Scattering Aerosol Optical Thickness (AOT)
        constexpr int FINE_MODE_EXTINCTION_AOT = 6;  // Fine mode (sub-micron) extinction AOT
        constexpr int FINE_MODE_SCATTERING_AOT = 7;  // Fine mode (sub-micron) scattering AOT
        constexpr int PM2_5_EXTINCTION_AOT = 8;     // PM2.5 extinction AOT
        constexpr int PM2_5_SCATTERING_AOT = 9;     // PM2.5 scattering AOT
        
        constexpr int NUM_OPTICS = 10;
    }

    class IAerosolPackage {
    public:
        virtual ~IAerosolPackage() = default;

        // Dynamic Initialization (Loads YAML config)
        virtual void initialize(const std::string& config_yaml) = 0;

        // Passive microphysics step
        virtual void executeMicrophysics(
            const EnvironmentalStateView& env,
            View3D<double>& state,          // (cell, level, species)
            double delta_time_sec) = 0;

        // Active Diagnostics Step (PM2.5, PM10, Derived Number, Column Mass, etc.)
        virtual void computeDerivedDiagnostics(
            const EnvironmentalStateView& env,
            const View3D<const double>& state,    // (cell, level, species)
            View3D<double>& diagnostics_out) = 0; // (cell, level, diagnostic_index)

        // Active Optical Properties Step (Extinction, Scattering, AOT, etc.)
        virtual void computeOptics(
            const EnvironmentalStateView& env,
            const View3D<const double>& state,    // (cell, level, species)
            View3D<double>& optics_out) = 0;      // (cell, level, optics_index)
    };

} // namespace exaero
```

### 5.2. Private GOCART Diagnostic & Optical Kernels (`src_impl/gocart/GocartSolver.cpp`)

For GOCART, `executeMicrophysics` remains a passive pass-through. The core GOCART implementations are the GPU-accelerated diagnostics and optical solvers:

#### A. Derived Number Concentration & Column-Mass Physics
GOCART tracks dynamic mass concentration $C_i$ [kg/m³]. Number Concentration $N_i$ [particles/m³] is derived by representing the bulk mass as a lognormal distribution:
$$N_i = \frac{C_i}{\frac{\pi}{6} \rho_i D_{g,i}^3 \exp\left(\frac{9}{2} \ln^2 \sigma_{g,i}\right)}$$
Surface mass represents the bottom vertical layer ($k = 0$). Column-integrated mass density is calculated by integrating concentrations vertically over grid layer thicknesses $\Delta z$:
$$M_{\text{col}} = \sum_{k} C_{k} \Delta z_{k}$$

#### B. Dual-Mode Aerosol Optics (`computeOptics`)
Our GPU solver evaluates optical characteristics per grid cell and level, dynamically selecting the mode based on `params.has_optics_lookup`:

##### 1. ADT (Anomalous Diffraction Theory) Mode (Analytical)
We calculate the extinction efficiency $Q_{\text{ext}}$ dynamically using Van de Hulst's formulation based on phase shift $\rho = 2x(n - 1)$ where $x = \frac{\pi D_{\text{wet}}}{\lambda}$ is the wet size parameter, $n$ is the real refractive index:
$$Q_{\text{ext}} = 2 - \frac{4}{\rho}\sin\rho + \frac{4}{\rho^2}(1 - \cos\rho)$$
Scattering efficiency is scaled based on imaginary index absorption $k$, and asymmetry $g$ is approximated as a smooth function of the size parameter (Rayleigh $g \approx 0$ to geometric $g \approx 0.7$). Extinction and scattering coefficients [m⁻¹] are derived by multiplying cross-section by efficiencies:
$$\beta_{\text{ext}, i} = N_i \cdot \frac{\pi D_{\text{wet}, i}^2}{4} \cdot Q_{\text{ext}, i}$$

##### 2. RH Lookup Table Mode (Tabulated Interpolation)
If lookups are enabled, we read pre-tabulated Mass Extinction Efficiency ($MEE_i$) [m²/g], Single Scattering Albedo ($SSA_i$), and asymmetry ($g_i$) at 8 standard RH bins (from $0\%$ to $99\%$). The GPU kernel executes a fast 1D linear interpolation to obtain the local optical variables at the cell's relative humidity:
$$\beta_{\text{ext}, i} = C_i \cdot MEE_i(RH) \cdot 10^3$$

The 3D coefficients are vertically integrated over layer thicknesses $\Delta z$ to return the standard Column Aerosol Optical Thickness (AOT):
$$\tau = \sum_{k} \beta_{\text{ext}, k} \Delta z_{k}$$

---

## 6. Phased Implementation Plan

### Phase 1: Foundational Scaffolding & Vendoring
*   Configure CMake build targets.
*   Vendor the `mdspan` stable backport into a single-header `src/src_impl/third_party/mdspan/include/experimental/mdspan`.
*   Establish public `IAerosolPackage` and `EnvironmentalStateView` structures.
*   Verify basic compilation with a simple CPU target.

### Phase 2: YAML Parsing & Dynamic Property Upload
*   Integrate `yaml-cpp` as a build dependency.
*   Implement GOCART `GocartSpeciesParams` mapping (incorporating ADT refractive indices and the 8-bin RH lookup arrays).
*   Implement `GocartPackage::initialize(...)` to read dynamic species arrays and upload them to the default GPU execution space.
*   Verify `make_unmanaged_kokkos_view` compilation on GPU space.

### Phase 3: GOCART Passive Diagnostics Proof-of-Concept
*   Implement `GocartPackage::executeMicrophysics` as an intentional passive NO-OP.
*   Implement `GocartPackage::computeDerivedDiagnostics` on the GPU, calculating wet size growth factors, Surface Area Density (SAD), derived number concentration, and column-integrated masses.
*   Add unit tests verifying precise diagnostic and column integration calculations.

### Phase 4: GOCART Dual-Mode Optical Solver
*   Implement `GocartPackage::computeOptics` on the GPU.
*   Program the **Anomalous Diffraction Theory (ADT)** analytical solver kernels.
*   Program the **RH 1D Linear Interpolation Lookup Table** solver kernels.
*   Implement 2D Column AOT integrations (vertical sum over $\beta_{\text{ext}} \cdot \Delta z$).
*   Add comprehensive unit tests verifying both ADT and RH table lookup calculations.

### Phase 5: MAM4xx Linking & Validation
*   Configure EX-aero to link against an external pre-built `mam4xx` library.
*   Incorporate MAM4xx under the unified `IAerosolPackage` interface.
*   Run validation suite and performance benchmarks verifying zero-copy overhead is `< 0.1%`.
