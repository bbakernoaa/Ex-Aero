# Project Charter: EX-aero (Exascale Aerosol Interface)

**Project Name:** EX-aero (**EX**ascale **AERO**sol Interface)

**Target Artifacts:** `libexaero.so` / `libexaero.a`

**Namespace:** `exaero::`

**Header Directory:** `<exaero/...>`

**Language Standard:** C++20 / Kokkos

**License:** Open Source (Apache 2.0 / BSD-3-Clause compatible)

---

## 1. Executive Summary & Vision

The transition of Earth System Models (ESMs) toward exascale, GPU-accelerated computing requires breaking down monolithic atmospheric composition codes. Within NOAA's Unified Forecast System (UFS), the shift of the **CATChem** (Configurable ATmospheric Chemistry) core to C++ and the deployment of **CECE** (Community Emissions Computing Engine) in C++20/Kokkos present an opportunity to unify atmospheric composition modeling.

**EX-aero** serves as an **independent, zero-dependency, layout-agnostic C++20 interface library**. Adhering to the principles established by the community-wide Generalized Aerosol/Chemistry Interface (GIANT) initiative, **EX-aero** decouples host dynamics and physics schemes from specific aerosol microphysics implementations.

Using standard `mdspan` view abstractions and unmanaged Kokkos wrappers, **EX-aero** enables zero-copy GPU memory exchange among chemistry cores (CATChem), emissions processors (CECE), and physics drivers (CCPP schemes like Thompson cloud microphysics and RRTMG radiation).

---

## 2. Scope & Separation of Concerns

To prevent architectural bloat and maintain high performance, **EX-aero** strictly enforces Separation of Concerns. It is an *aerosol state and microphysics engine*, not a complete atmospheric model.

**What EX-aero DOES:**

* **Aerosol Microphysics:** Condensation, coagulation, nucleation, and gas-particle partitioning.


* **Aerosol Optics:** Computes Extinction, Single Scattering Albedo (SSA), and Asymmetry parameters based on size and refractive indices.


* **Aerosol Diagnostics:** Derives PM2.5, PM10, Number Concentration, and Cloud Condensation Nuclei (CCN) / Ice Nucleating Particles (INP) from the aerosol state.


* **Deposition Physics:** Calculates dry deposition and wet scavenging removal tendencies.



**What EX-aero DOES NOT do (Handled by Host/CATChem):**

* **Advection/Transport:** The dynamical core (e.g., FV3) moves tracers.
* **Gas-Phase & Aqueous Chemistry:** CATChem's solvers handle chemical kinetics (using EX-aero's Surface Area Density and aerosol liquid water).


* **Radiative Transfer & Photolysis:** RRTMG and Fast-J use EX-aero's optics to compute heating rates and actinic fluxes.



---

## 3. High-Level Objectives

1. **Zero-Dependency Architecture:** Implement `libexaero` as an independent software library free of CATChem, CECE, CCPP, or host model headers.
2. **Zero-Copy Memory Boundary:** Utilize C++20 `std::experimental::mdspan` view abstractions to wrap memory pointers zero-copy across C++ and Fortran boundaries.


3. **GPU Performance Portability:** Leverage unmanaged `Kokkos::View` structures inside aerosol package solvers to execute GPU parallel kernels directly on host-allocated device memory.


4. **Dual Launch Packages:**
* **MAM4xx:** Modern 4-mode lognormal modal aerosol package in C++/Kokkos tracking mass and number concentrations.


* **GOCART:** Fully configurable bulk aerosol package driven dynamically via YAML runtime files.


5. **CCPP Constituent Array Interoperability:** Include a native `CcppConstituentAdapter` that allows CCPP physics schemes to slice into CCPP's 3D tracer array zero-copy.



---

## 4. Dimensional Units & Standards

To ensure clarity across chemistry, aerosol, and cloud microphysics modules:

* **Aerosol Phase Species:** Expressed in volumetric concentration: **μm³/m³** (cubic micrometers of aerosol volume per cubic meter of air).
* **Aerosol Number Concentration:** Expressed in numerical concentration: **particles/m³**.
* **Surface Emission Fluxes:** Expressed in volumetric flux: **μm³/m²/s**.
* **Gas-to-Particle Conversion:** Driven by the host's chemistry solver using species molecular weight and dry particle density to convert from **mol/m³** (gas) to **μm³/m³** (aerosol).

---

## 5. Generic Physics Engines

Instead of duplicating logic for bulk and modal schemes, EX-aero relies on three generic, Kokkos-accelerated engines to standardize physical processes.

### A. Emissions Mapping Engine

Surface emission schemes (CECE, dust/sea salt models) output fluxes across fixed physical size bands.

* **Function:** Applies a GPU-accelerated fractional overlap matrix to re-bin raw emission size ranges into package-specific modes or bins.
* **Number Generation:** For modal packages (MAM4xx), automatically converts incoming volume fluxes into Emitted Number Fluxes using characteristic emitted particle volumes.

### B. Generic Deposition Engine

Calculates dry deposition and wet scavenging by treating bulk schemes as a special case of a modal distribution.

* **Function:** Computes size-resolved physics (Stokes drag, Brownian diffusion, precipitation collision efficiency) based purely on particle diameter and density.


* **Mechanism:** Integrates across a size distribution curve. For modal (MAM4xx), it integrates across the lognormal curve (geometric standard deviation > 1.0) yielding distinct removal rates for mass and number. For bulk (GOCART), it evaluates a monodisperse size (geometric standard deviation = 1.0).

### C. Derived Parameter & Diagnostics Engine

Bridges bulk models with advanced cloud microphysics and data assimilation.

* **Function:** Calculates Number Concentration, Ice-Friendly Aerosols (INP), Water-Friendly Aerosols (CCN), Surface Area Density (SAD), PM2.5, and PM10.


* **Mechanism:** Derives these metrics for bulk schemes (GOCART) using assumed lognormal parameters loaded via YAML. For modal schemes (MAM4xx), computes dynamic metrics directly from prognostic mass and number.



---

## 6. Core Software API

The API enforces a minimal, strictly scoped interaction model during the host's timestep.

```cpp
namespace exaero {
    class IAerosolPackage {
    public:
        virtual ~IAerosolPackage() = default;

        // 1. Initialization
        virtual void initialize(const std::string& config_yaml) = 0;

        // 2. Microphysics Evolution
        virtual void executeMicrophysics(
            const EnvironmentalStateView<>& env,
            const SurfaceFluxView<>& emissions,
            AerosolStateView<>& state, // Configured for CCPP adapter or Native C++
            double delta_time_sec) = 0;

        // 3. Optics (For RRTMG / Photolysis)
        virtual void computeOptics(
            const EnvironmentalStateView<>& env,
            const AerosolStateView<>& state,
            OpticalPropertyView<>& optics_out) = 0;

        // 4. Deposition (For Land/Ocean coupling and Atmospheric sink)
        virtual void computeDeposition(
            const EnvironmentalStateView<>& env,
            const SurfaceEnvironmentView<>& sfc_env,
            const PrecipitationStateView<>& precip,
            AerosolStateView<>& state_tendencies,
            SurfaceDepositionFluxView<>& sfc_fluxes_out) = 0;

        // 5. Diagnostics (For Cloud Microphysics & Data Assimilation)
        virtual void computeDerivedDiagnostics(
            const EnvironmentalStateView<>& env,
            const AerosolStateView<>& state,
            DerivedDiagnosticViews<>& out_diagnostics) = 0;
    };
}

```

---

## 7. Success Criteria & Verification Framework

1. **Strict Decoupling:** Compilation of `libexaero` requires zero external host headers (only C++20 standard library, `mdspan`, and Kokkos).
2. **Unit Test Coverage:** Achieve ≥ 80% code coverage on all non-GPU host routines, view adapters, and generic physics engines.


3. **Zero Overhead Verification:** Validate that wrapping host memory into `mdspan` and unmanaged `Kokkos::View` incurs < 0.1% performance overhead compared to raw pointer operations.
4. **Data Locality:** Confirm via profiling tools (NVIDIA Nsight Systems / AMD ROCm) that zero PCI-e host-device transfers occur during CECE → CATChem → EX-aero execution.
5. **Multi-Platform CI/CD:** Maintain continuous integration testing on CPU (OpenMP), NVIDIA GPUs (CUDA), AMD GPUs (HIP), and Intel GPUs (SYCL).