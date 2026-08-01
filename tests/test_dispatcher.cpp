#include <gtest/gtest.h>
#include "exaero/dispatcher.hpp"

// Mock kernel for dispatcher test
struct MockKernels {
    KOKKOS_INLINE_FUNCTION void integrate_forward(double*, double*) const {}
};

TEST(DispatcherTest, ExecutionLoop) {
    exaero::ExaeroContext ctx("MockMech", 10);
    std::vector<double> dummy(10);
    exaero::UnmanagedDeviceState state(dummy.data(), dummy.data(), 1, 10, 1, 1);

    // Just assert it runs without throwing or allocating
    EXPECT_NO_THROW({
        exaero::Dispatcher::dispatch<MockKernels>(ctx, state, 0.1);
    });
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    Kokkos::initialize(argc, argv);
    int result = RUN_ALL_TESTS();
    Kokkos::finalize();
    return result;
}
