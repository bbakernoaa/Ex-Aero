#include <exaero/IAerosolPackage.hpp>
#include <exaero/MemoryMapper.hpp>
#include <exaero/Environment.hpp>
#include <cassert>
#include <iostream>

void test_zero_copy_mapping() {
    double data_raw[2 * 2 * 2] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0 };
    
    // Wrap raw array in public C++ View3D mdspan (zero-copy)
    exaero::View3D<double> span(data_raw, 2, 2, 2);
    
    // Map standard mdspan to unmanaged Kokkos View
    auto view = exaero::make_unmanaged_kokkos_view(span);
    
    // 1. Verify indices and extents match
    assert(view.extent(0) == 2);
    assert(view.extent(1) == 2);
    assert(view.extent(2) == 2);
    
    // 2. Verify values are exactly mapped
    assert(view(0, 0, 0) == 1.0);
    assert(view(1, 1, 1) == 8.0);
    
    // 3. Verify zero-copy: modifications inside the View write back to raw host memory in-place
    view(0, 1, 0) = 999.0;
    assert(data_raw[2] == 999.0); // Modified in-place!
    
    std::cout << "Memory Mapping Zero-Copy Unit Test: PASS" << std::endl;
}

int main() {
    exaero::initialize_environment();
    
    test_zero_copy_mapping();
    
    exaero::finalize_environment();
    return 0;
}
