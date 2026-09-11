#pragma once

namespace exaero {

// Public functions to initialize and finalize the Kokkos runtime
// without exposing any Kokkos headers to host applications.
void initialize_environment();
void finalize_environment();

// Block until all previously launched device work has completed
// (Kokkos::fence), without exposing Kokkos headers. Needed by timing/hot-path
// consumers.
void fence_environment();

} // namespace exaero
