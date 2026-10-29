#pragma once
// Test-only initialized packages and error assertions for host-coupling tests.
#include <algorithm>
#include <cassert>
#include <cmath>
#include <exaero/PackageFactory.hpp>
#include <exaero/IAerosolPackage_C.h>
#include <stdexcept>
namespace host_test {
inline const char *bulk_yaml = R"(species:
  - {name: A, dry_density: 1000, molecular_weight: 100, dry_particle_diameter: 8e-7, hygroscopicity: 0.6, lognormal_sigma: 1.5, lognormal_dg: 2e-7, refractive_index_real: 1.55, refractive_index_imag: 0.002}
  - {name: B, dry_density: 2000, molecular_weight: 100, dry_particle_diameter: 2e-6, hygroscopicity: 0, lognormal_sigma: 1.5, lognormal_dg: 4e-7, refractive_index_real: 1.55, refractive_index_imag: 0.002}
)";
inline const char *modal_yaml = R"(backend: modal
coagulation: false
species: [{name: A, density: 1000, molecular_weight: 0.1, hygroscopicity: 0.6}, {name: B, density: 2000, molecular_weight: 0.2, hygroscopicity: 0}]
modes: [{name: fine, geometric_radius: 1e-7, geometric_std_dev: 1.5}, {name: coarse, geometric_radius: 2e-7, geometric_std_dev: 1.5, species: [A]}]
)";
inline auto make(bool modal) {
  auto p=exaero::createPackage(modal ? "modal":"gocart");
  p->initialize(modal ? modal_yaml : bulk_yaml); return p;
}
inline bool near(double a,double b) {
  return std::abs(a-b)<=1e-11*std::max({1e-300,std::abs(a),std::abs(b)});
}
template<class F> void rejects(F f) { bool failed=false; try { f(); } catch(const std::exception&) { failed=true; } assert(failed); }
struct Runtime { Runtime(){exaero_init_environment();} ~Runtime(){exaero_finalize_environment();} };
}
