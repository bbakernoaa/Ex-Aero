#include <gtest/gtest.h>
#include "auxiliary_engines.hpp"
#include <vector>

TEST(AuxiliaryEnginesTest, ZeroCopyExtractionDispatch) {
    Kokkos::initialize();
    {
        // Allocate raw host memory
        std::vector<double> dummy_meteo(10 * 10 * 10 * 1);
        std::vector<double> dummy_conc(5 * 10 * 10 * 1);
        
        // Populate dummy array to test pointer extraction mapping
        dummy_meteo[0] = 288.15; // e.g. T
        dummy_conc[0]  = 0.01;   // e.g. O3

        exaero::UnmanagedDeviceState state(dummy_conc.data(), dummy_meteo.data(), 5, 10, 10, 1);
        exaero::ExaeroContext ctx("ChapmanCycle", 100);

        EXPECT_NO_THROW({
            exaero::AuxiliaryEngines::compute_photolysis(ctx, state, 0.1);
            exaero::AuxiliaryEngines::compute_thermodynamics(ctx, state, 0.1);
        });

        // The underlying Kokkos LayoutLeft unmanaged view's .data() pointer MUST point to exactly
        // the same memory address as the host std::vector data(), guaranteeing zero-copy mapping (SC-003).
        EXPECT_EQ(state.meteorology.data(), dummy_meteo.data());
        EXPECT_EQ(state.concentrations.data(), dummy_conc.data());
    }
    Kokkos::finalize();
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
    return result;
}
