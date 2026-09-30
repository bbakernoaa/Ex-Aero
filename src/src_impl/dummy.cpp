/// @file dummy.cpp
/// @brief Placeholder translation unit for the @c exaero_impl target.
///
/// Historically needed to satisfy CMake until GocartPackage was
/// implemented; kept as a stable link anchor for the implementation library.

/// @brief No-op anchor symbol; guarantees the library has at least one
/// exported definition even when built without optional components.
void exaero_impl_dummy_function() {}
