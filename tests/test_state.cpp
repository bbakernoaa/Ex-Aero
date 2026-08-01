#include <gtest/gtest.h>
#include "exaero/state.hpp"
#include <vector>

TEST(ExaeroStateTest, Instantiation) {
    Kokkos::initialize();
    {
        exaero::ExaeroContext ctx("ChapmanCycle", 100);
        EXPECT_EQ(ctx.active_mechanism, "ChapmanCycle");
        EXPECT_EQ(ctx.total_grid_cells, 100);
        EXPECT_EQ(ctx.sza_sorted_indices.extent(0), 100);

        std::vector<double> dummy_conc(5 * 10 * 10 * 1);
        std::vector<double> dummy_met(10 * 10 * 10 * 1);

        exaero::UnmanagedDeviceState state(dummy_conc.data(), dummy_met.data(), 5, 10, 10, 1);
        
        EXPECT_EQ(state.concentrations.extent(0), 5);
        EXPECT_EQ(state.concentrations.extent(1), 10);
        EXPECT_EQ(state.meteorology.extent(0), 10);
    }
    Kokkos::finalize();
}
