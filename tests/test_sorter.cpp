#include <gtest/gtest.h>
#include "exaero/sorter.hpp"
#include <vector>

static size_t dynamic_allocations = 0;

extern "C" void kokkos_profiling_allocate_data(const Kokkos::Profiling::SpaceHandle, const char*, const void*, const uint64_t) {
    dynamic_allocations++;
}

extern "C" void kokkos_profiling_free_data(const Kokkos::Profiling::SpaceHandle, const char*, const void*, const uint64_t) {}

TEST(SorterTest, ZeroAllocationDuringSort) {
    std::vector<double> dummy_conc(5 * 100);
    std::vector<double> dummy_met(10 * 100);

    // Put reverse sorted SZA to force it to work
    for(int i=0; i<100; ++i) {
        dummy_met[i] = 100.0 - i;
    }

    exaero::UnmanagedDeviceState state(dummy_conc.data(), dummy_met.data(), 5, 100, 1, 1);
    exaero::ExaeroContext ctx("ChapmanCycle", 100);
    
    // Reset counter before tracking the integration loop
    dynamic_allocations = 0;

    exaero::SZA_Sorter::sort_workload(ctx, state);

    // Track memory (T016)
    EXPECT_EQ(dynamic_allocations, 0) << "Sorting triggered dynamic memory allocations!";

    auto h_indices = Kokkos::create_mirror_view(ctx.sza_sorted_indices);
    Kokkos::deep_copy(h_indices, ctx.sza_sorted_indices);

    EXPECT_EQ(h_indices(0), 99);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    Kokkos::initialize(argc, argv);
    int result = RUN_ALL_TESTS();
    Kokkos::finalize();
    return result;
}
