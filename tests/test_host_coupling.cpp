#ifdef NDEBUG
#undef NDEBUG
#endif
#include "../src/src_impl/HostCouplingInternal.hpp"
#include "host_coupling_fixture.hpp"
#include <cstring>
#include <exaero/GocartConfig.hpp>
#include <exaero/HostCoupling.hpp>
#include <exaero/HostCoupling_C.h>
#include <exaero/ModalCoupling.hpp>
int main() {
  host_test::Runtime runtime;
  for (bool modal : {false, true}) {
    auto p = host_test::make(modal);
    auto info = exaero::hostCouplingInfo(*p);
    assert(info.populations.size() == 2);
    // Extent validation must account for each output axis once, without
    // allocation.
    (void)exaero::host_detail::dimensions(*p, 10000000, 1, info.field_count, 2,
                                          15);
    host_test::rejects([&] {
      exaero::host_detail::dimensions(*p, std::numeric_limits<int>::max(), 1,
                                      info.field_count, 2, 15);
    });
    assert(info.field_count == (modal ? 6 : 2));
    assert(info.constituents.size() == (modal ? 4 : 2));
    for (int i = 0; i < 2; ++i) {
      assert(info.populations[i].index == i);
      assert(info.populations[i].number_field == (modal ? i * 3 + 2 : -1));
      assert(info.populations[i].number_kind ==
             (modal ? exaero::NumberKind::Prognostic
                    : exaero::NumberKind::Diagnosed));
    }
    assert(info.constituents.back().mass_field == (modal ? 4 : 1));
    assert(info.constituents.back().permitted == !modal);
    assert(exaero::supportsHostOperation(
               *p, exaero::HostOperation::InspectNumberBounds) == modal);
    int rep = -9, nf = -9, np = -9, ns = -9, flag = -1;
    char error[256];
    exaero_host_coupling_counts(p.get(), &rep, &nf, &np, &ns, error, &flag);
    assert(flag == 0 && rep == int(modal) && nf == info.field_count &&
           np == 2 && ns == int(info.constituents.size()));
    int idx = -9, kind = -9, field = -9;
    char name[32];
    exaero_host_population_info(p.get(), 1, &idx, &kind, &field, name, 32,
                                error, &flag);
    assert(flag == 0 && idx == 1 && field == (modal ? 5 : -1));
    const auto old = field;
    name[0] = '!';
    exaero_host_population_info(p.get(), 1, &idx, &kind, &field, name, 1, error,
                                &flag);
    assert(flag != 0 && field == old && name[0] == '!');
    int pop = -9, species = -9, mass = -9, allowed = -9;
    exaero_host_constituent_info(p.get(), ns - 1, &pop, &species, &mass,
                                 &allowed, name, 32, error, &flag);
    assert(flag == 0 && pop == 1 && species == 1 && allowed == int(!modal));
    assert(mass == (modal ? 4 : 1));
    exaero_host_constituent_info(p.get(), ns, &pop, &species, &mass, &allowed,
                                 name, 32, error, &flag);
    assert(flag != 0 && mass == (modal ? 4 : 1));
    int supported = 77;
    exaero_supports_host_operation(p.get(), 999, &supported, error, &flag);
    assert(flag != 0 && supported == 77);
    p->initialize(modal ? host_test::modal_yaml : host_test::bulk_yaml);
    assert(exaero::hostCouplingInfo(*p).populations.size() == 2);
    if (modal) {
      auto cfg = exaero::modalCoupling(*p).modalConfig();
      std::reverse(cfg.species.begin(), cfg.species.end());
      std::reverse(cfg.modes.begin(), cfg.modes.end());
      exaero::initializeModal(*p, cfg);
    } else {
      const std::string yaml = host_test::bulk_yaml;
      const auto b = yaml.find("  - {name: B");
      p->initialize("species:\n" + yaml.substr(b) +
                    yaml.substr(yaml.find("  - {name: A"),
                                b - yaml.find("  - {name: A")));
    }
    const auto reordered = exaero::hostCouplingInfo(*p);
    assert(reordered.constituents.front().name == "B" &&
           reordered.constituents.front().mass_field == 0);
    assert(reordered.populations.front().name == (modal ? "coarse" : "B"));
    assert(reordered.constituents.front().permitted == !modal);
  }
  {
    auto p = exaero::createPackage("gocart");
    std::string invalid = host_test::bulk_yaml;
    invalid.replace(invalid.find("name: B,"), 8,
                    "name: B, dry_radius_lower: 1e-7,");
    for (bool prior_success : {false, true}) {
      if (prior_success)
        p->initialize(host_test::bulk_yaml);
      host_test::rejects([&] { p->initialize(invalid); });
      host_test::rejects([&] { exaero::hostCouplingInfo(*p); });
      for (int op = 0; op <= 8; ++op)
        assert(!exaero::supportsHostOperation(
            *p, static_cast<exaero::HostOperation>(op)));
      int rep = 77, nf = 77, np = 77, ns = 77, flag;
      char error[256];
      exaero_host_coupling_counts(p.get(), &rep, &nf, &np, &ns, error, &flag);
      assert(flag != 0 && rep == 77 && nf == 77 && np == 77 && ns == 77);
      p->initialize(host_test::bulk_yaml);
      assert(exaero::hostCouplingInfo(*p).populations.size() == 2);
    }
    host_test::rejects([&] { p->initialize("species: ["); });
    assert(!exaero::supportsHostOperation(*p,
                                          exaero::HostOperation::RemovalRates));
    p->initialize(host_test::bulk_yaml);
    exaero::GocartConfig cfg;
    exaero::GocartSpeciesConfig species;
    species.name = "A";
    species.dry_density = 1000;
    species.molecular_weight = 100;
    species.dry_particle_diameter = 2e-7;
    species.lognormal_dg = 2e-7;
    species.lognormal_sigma = 1.5;
    species.refractive_index_real = 1.5;
    species.refractive_index_imag = .002;
    species.dry_radius_lower = 1e-7;
    cfg.species.push_back(species);
    host_test::rejects([&] { p->initialize(cfg); });
    assert(!exaero::supportsHostOperation(*p,
                                          exaero::HostOperation::RemovalRates));
    p->initialize(host_test::bulk_yaml);
    assert(exaero::hostCouplingInfo(*p).populations.size() == 2);
  }
  {
    std::string yaml = "species:\n", line = host_test::bulk_yaml;
    line = line.substr(line.find("  - {name: B"));
    for (int i = 0; i < 100; ++i) {
      auto species = line;
      species.replace(species.find("name: B"), 7,
                      "name: bin" + std::to_string(i));
      yaml += species;
    }
    auto p = exaero::createPackage("gocart");
    p->initialize(yaml);
    // 30 million actual elements, not three billion: rebin has no population
    // output axis.
    exaero::host_detail::rebinDimensions(*p, 300000, 1, 100, 100);
    host_test::rejects([&] {
      exaero::host_detail::rebinDimensions(*p, std::numeric_limits<int>::max(),
                                           1, 100, 100);
    });
  }
  auto uninitialized = exaero::createPackage("modal");
  assert(!exaero::supportsHostOperation(*uninitialized,
                                        exaero::HostOperation::RemovalRates));
  host_test::rejects([&] { exaero::hostCouplingInfo(*uninitialized); });
  if (exaero::packageAvailable("mam4xx")) {
    auto native = exaero::createPackage("mam4xx");
    native->initialize("backend: mam4xx");
    assert(!exaero::supportsHostOperation(*native,
                                          exaero::HostOperation::RemovalRates));
    host_test::rejects([&] { exaero::hostCouplingInfo(*native); });
  }
}
