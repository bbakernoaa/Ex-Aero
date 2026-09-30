#pragma once
/// @file Environment.hpp
/// @brief Kokkos runtime lifecycle wrapped behind a Kokkos-free header.
///
/// Host applications (and the C/Fortran CCPP layer) start and stop the
/// performance-portability runtime through these three functions without
/// ever including <Kokkos_Core.hpp>, keeping the public interface
/// dependency-free (ADR-001).
namespace exaero {

/// @brief Initialize the Kokkos runtime (idempotent guard is the
/// implementation's responsibility).
///
/// Must be called once before any compute step that launches device work.
/// Wraps @c Kokkos::initialize without exposing Kokkos headers to callers.
void initialize_environment();

/// @brief Finalize the Kokkos runtime and release device resources.
///
/// Call once at teardown after all compute steps have completed. Wraps
/// @c Kokkos::finalize without exposing Kokkos headers to callers.
void finalize_environment();

/// @brief Block until all previously launched device work has completed.
///
/// Wraps @c Kokkos::fence without exposing Kokkos headers. Needed by
/// timing/hot-path consumers that must not read output buffers until the
/// asynchronous device kernels have drained.
void fence_environment();

} // namespace exaero
