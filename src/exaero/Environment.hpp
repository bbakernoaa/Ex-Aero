#pragma once

namespace exaero {

    // Public functions to initialize and finalize the Kokkos runtime
    // without exposing any Kokkos headers to host applications.
    void initialize_environment();
    void finalize_environment();

} // namespace exaero
