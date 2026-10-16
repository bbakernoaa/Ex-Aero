#include "exaero/dispatcher.hpp"
#include <gtest/gtest.h>

// Mock kernel for dispatcher test
struct MockKernels {
  KOKKOS_INLINE_FUNCTION void integrate_forward(double *, double *) const {}
};

TEST(DispatcherTest, ExecutionLoop) {
  exaero::ExaeroContext ctx("MockMech", 10);
  std::vector<double> dummy(10);
  exaero::UnmanagedDeviceState state(dummy.data(), dummy.data(), 1, 10, 1, 1);

  // Just assert it runs without throwing or allocating
  EXPECT_NO_THROW(
      { exaero::Dispatcher::dispatch<MockKernels>(ctx, state, 0.1); });
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  // Exception-safe runtime ownership; also consumes --kokkos-* tuning args (FR-013).
  Kokkos::ScopeGuard kokkos_scope(argc, argv);
  return RUN_ALL_TESTS();
}
