#include <gtest/gtest.h>
#include "exaero/registry.hpp"
#include <stdexcept>
#include <vector>

// Note: Kokkos initialization/finalization is handled by the test runner's global environment setup

TEST(RegistryTest, UnknownMechanismThrows) {
    std::vector<double> dummy_buffer(10);
    exaero::UnmanagedDeviceState dummy_state(dummy_buffer.data(), dummy_buffer.data(), 1, 1, 1, 1);
    exaero::ExaeroContext ctx("UnknownMechanism", 100);
    
    EXPECT_FALSE(exaero::Registry::has_mechanism("UnknownMechanism"));
    EXPECT_THROW(exaero::Registry::dispatch(ctx, dummy_state, 0.1), std::runtime_error);
}

TEST(RegistryTest, KnownMechanismResolves) {
    std::vector<double> dummy_buffer(10);
    exaero::UnmanagedDeviceState dummy_state(dummy_buffer.data(), dummy_buffer.data(), 1, 1, 1, 1);
    exaero::ExaeroContext ctx("ChapmanCycle", 100);
    
    EXPECT_TRUE(exaero::Registry::has_mechanism("ChapmanCycle"));
    EXPECT_NO_THROW(exaero::Registry::dispatch(ctx, dummy_state, 0.1));
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    Kokkos::initialize(argc, argv);
    int result = RUN_ALL_TESTS();
    Kokkos::finalize();
    return result;
}
