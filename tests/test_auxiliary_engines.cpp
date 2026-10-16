#include "auxiliary_engines.hpp"
#include <gtest/gtest.h>
#include <stdexcept>
#include <vector>

// The auxiliary engines are host-confined and optional. In the default build
// (EXAERO_WITH_CLOUDJ / EXAERO_WITH_ISORROPIALITE OFF) the wrapped Fortran
// kernels are not linked, so the entry points MUST fail loud rather than
// launch an empty team and silently pretend to have computed physics. When an
// engine IS compiled in, the same call dispatches on the host execution space
// and must not throw for a well-formed input.
TEST(AuxiliaryEnginesTest, PhotolysisIsNotASilentNoOp) {
  std::vector<double> dummy_meteo(10 * 10 * 10 * 1, 288.15);
  std::vector<double> dummy_conc(5 * 10 * 10 * 1, 0.01);
  exaero::UnmanagedDeviceState state(dummy_conc.data(), dummy_meteo.data(), 5,
                                     10, 10, 1);
  exaero::ExaeroContext ctx("ChapmanCycle", 100);

#ifndef EXAERO_WITH_CLOUDJ
  EXPECT_THROW(exaero::AuxiliaryEngines::compute_photolysis(ctx, state, 0.1),
               std::runtime_error);
#endif
}

TEST(AuxiliaryEnginesTest, ThermodynamicsIsNotASilentNoOp) {
  std::vector<double> dummy_meteo(10 * 10 * 10 * 1, 288.15);
  std::vector<double> dummy_conc(5 * 10 * 10 * 1, 0.01);
  exaero::UnmanagedDeviceState state(dummy_conc.data(), dummy_meteo.data(), 5,
                                     10, 10, 1);
  exaero::ExaeroContext ctx("ChapmanCycle", 100);

#ifndef EXAERO_WITH_ISORROPIALITE
  EXPECT_THROW(
      exaero::AuxiliaryEngines::compute_thermodynamics(ctx, state, 0.1),
      std::runtime_error);
#endif
}

TEST(AuxiliaryEnginesTest, ZeroCopyPointerMapping) {
  // Allocate raw host memory and populate it so the pointer-extraction
  // mapping can be inspected.
  std::vector<double> dummy_meteo(10 * 10 * 10 * 1);
  std::vector<double> dummy_conc(5 * 10 * 10 * 1);
  dummy_meteo[0] = 288.15; // e.g. T
  dummy_conc[0] = 0.01;    // e.g. O3

  exaero::UnmanagedDeviceState state(dummy_conc.data(), dummy_meteo.data(), 5,
                                     10, 10, 1);

  // The underlying Kokkos LayoutLeft unmanaged view's.data() pointer MUST
  // point to exactly the same memory address as the host std::vector data(),
  // guaranteeing zero-copy mapping.
  EXPECT_EQ(state.meteorology.data(), dummy_meteo.data());
  EXPECT_EQ(state.concentrations.data(), dummy_conc.data());
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  // Exception-safe runtime ownership; also consumes --kokkos-* tuning args.
  Kokkos::ScopeGuard kokkos_scope(argc, argv);
  return RUN_ALL_TESTS();
}
