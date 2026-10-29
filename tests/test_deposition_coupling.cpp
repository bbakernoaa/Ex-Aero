#ifdef NDEBUG
#undef NDEBUG
#endif
#include "host_coupling_fixture.hpp"
#include <exaero/BulkRemoval.hpp>
#include <exaero/DepositionCoupling.hpp>
#include <exaero/DepositionCoupling_C.h>
#include <exaero/ModalCoupling.hpp>
#include <exaero/ModalRemoval.hpp>
#include <limits>
int main() {
  host_test::Runtime runtime;
  using namespace exaero;
  for (bool modal : {false, true}) {
    auto p = host_test::make(modal);
    if (modal) {
      auto cfg = modalCoupling(*p).modalConfig();
      auto third = cfg.modes.back();
      third.name = "largest";
      cfg.modes.push_back(third);
      cfg.resuspension_routes = {{"A", "fine", "coarse"},
                                 {"A", "coarse", "largest"}};
      initializeModal(*p, cfg);
    } else {
      std::string yaml = host_test::bulk_yaml;
      auto extra = std::string(host_test::bulk_yaml);
      const auto pos = extra.find("  - {name: B");
      extra = extra.substr(pos);
      auto c = extra;
      c.replace(c.find("name: B"), 7, "name: C");
      auto d = extra;
      d.replace(d.find("name: B"), 7, "name: D");
      p->initialize(yaml + c + d);
    }
    auto info = hostCouplingInfo(*p);
    const int nf = info.field_count, nc = 2, nl = 3;
    const int fa = 0, fb = modal ? 3 : 1, fc = modal ? 6 : 2,
              fd = modal ? 1 : 3;
    std::vector<double> returns(nc * nl * nf), out(returns.size()),
        reference(out.size());
    View3D<double> x(returns.data(), nc, nl, nf), y(out.data(), nc, nl, nf),
        ref(reference.data(), nc, nl, nf);
    for (int k = 0; k < nl; ++k)
      for (int c = 0; c < nc; ++c) {
        x(c, k, fa) = 10 * (c + k + 1);
        x(c, k, fb) = 3 * (c + k + 1);
        x(c, k, fd) = 4 * (c + k + 1);
      }
    View3D<const double> input(returns.data(), nc, nl, nf);
    std::vector<ResuspensionFieldRoute> routes{
        {std::size_t(fa), std::size_t(fb)}, {std::size_t(fb), std::size_t(fc)}};
    redistributeResuspension(*p, input, routes, y);
    assert(y(0, 0, fa) == 0 && y(0, 0, fb) == 10 && y(0, 0, fc) == 3 &&
           y(0, 0, fd) == 4);
    if (modal) {
      assert(configuredResuspensionRoutes(*p).size() == 2);
      redistributeModalResuspension(*p, input, ref);
    } else {
      assert(configuredResuspensionRoutes(*p).empty());
      std::vector<BulkResuspensionRoute> br{{0, 1}, {1, 2}};
      redistributeBulkResuspension(*p, input, br, ref);
    }
    assert(out == reference);
    int sources[2]{fa, fb}, destinations[2]{fb, fc}, flag;
    char error[256];
    exaero_redistribute_resuspension(p.get(), nc, nl, nf, 2, sources,
                                     destinations, returns.data(), out.data(),
                                     error, &flag);
    assert(flag == 0 && out == reference);
    int required = 77;
    exaero_resuspension_routes(p.get(), 0, &required, nullptr, nullptr, error,
                               &flag);
    assert(flag == 0 && required == (modal ? 2 : 0));
    if (modal) {
      int srcs[2], dsts[2];
      exaero_resuspension_routes(p.get(), 2, &required, srcs, dsts, error,
                                 &flag);
      assert(flag == 0 && srcs[0] == 0 && dsts[0] == 3 && srcs[1] == 3 &&
             dsts[1] == 6);
      required = 77;
      int src = 88, dst = 99;
      exaero_resuspension_routes(p.get(), 1, &required, &src, &dst, error,
                                 &flag);
      assert(flag != 0 && required == 77 && src == 88 && dst == 99);
    }
    std::vector<std::vector<ResuspensionFieldRoute>> invalid{
        {{0, 0}},
        {{0, std::size_t(fb)}, {0, std::size_t(fc)}},
        {{0, std::size_t(fb)}, {std::size_t(fb), 0}},
        {{0, std::size_t(nf)}}};
    if (modal) {
      invalid.push_back({{0, 2}});
      invalid.push_back({{0, 1}});
      invalid.push_back({{1, 4}});
    }
    for (const auto &r : invalid) {
      std::fill(out.begin(), out.end(), 77);
      host_test::rejects([&] { redistributeResuspension(*p, input, r, y); });
      assert(out.front() == 77 && out.back() == 77);
    }
    // Routed intermediate originals leave before incoming material is
    // accumulated.
    x(1, 2, fa) = std::numeric_limits<double>::max();
    x(1, 2, fb) = std::numeric_limits<double>::max();
    redistributeResuspension(*p, input, routes, y);
    assert(y(1, 2, fb) == std::numeric_limits<double>::max() &&
           y(1, 2, fc) == std::numeric_limits<double>::max());
    x(1, 2, fc) = std::numeric_limits<double>::max();
    std::fill(out.begin(), out.end(), 77);
    host_test::rejects([&] { redistributeResuspension(*p, input, routes, y); });
    assert(out.front() == 77 && out.back() == 77);
  }
  for (bool modal : {false, true}) {
    auto p = host_test::make(modal);
    host_test::Grid g(modal);
    std::fill(g.state.begin(), g.state.end(), 0);
    View3D<double> x(g.state.data(), g.nc, g.nl, g.nf);
    for (int k = 0; k < g.nl; ++k)
      for (int c = 0; c < g.nc; ++c) {
        x(c, k, 0) = 2e-9 * (c + k + 1);
        if (modal)
          x(c, k, 2) = 5e7 * (c + k + 1);
      }
    double edge[4]{1e-12, 2e-7, 4e-7, 1e-3};
    View1D<const double> edges(edge, 4);
    std::vector<double> bins(g.nc * g.nl * 3), tails(g.nc * g.nl * 2),
        rb(bins.size()), rt(tails.size());
    View3D<double> b(bins.data(), g.nc, g.nl, 3),
        t(tails.data(), g.nc, g.nl, 2), br(rb.data(), g.nc, g.nl, 3),
        tr(rt.data(), g.nc, g.nl, 2);
    rebinDeposition(*p, g.input(), 0, edges, b, t);
    if (modal)
      rebinModalDeposition(*p, g.input(), 0, edges, br, tr);
    else
      rebinBulkDeposition(*p, g.input(), 0, edges, br, tr);
    assert(bins == rb && tails == rt);
    const double ln = std::log(1.5),
                 dg =
                     modal
                         ? std::cbrt(6 * 2e-9 / (std::acos(-1.) * 1000 * 5e7)) *
                               std::exp(-1.5 * ln * ln)
                         : 2e-7;
    auto cdf = [&](double d) {
      return .5 *
             std::erfc((std::log(dg / d) + 3 * ln * ln) / (std::sqrt(2.) * ln));
    };
    for (int j = 0; j < 3; ++j)
      assert(host_test::near(b(0, 0, j),
                             2e-9 * (cdf(edge[j + 1]) - cdf(edge[j]))));
    assert(host_test::near(
        t(0, 0, 0) + t(0, 0, 1) + b(0, 0, 0) + b(0, 0, 1) + b(0, 0, 2), 2e-9));
    int flag;
    char error[256];
    exaero_rebin_deposition(p.get(), g.nc, g.nl, g.nf, 0, 3, g.state.data(),
                            edge, bins.data(), tails.data(), error, &flag);
    assert(flag == 0 && bins == rb && tails == rt);
    if (modal) {
      const auto original = bins;
      for (auto &v : g.state)
        v *= 2;
      rebinDeposition(*p, g.input(), 0, edges, b, t);
      for (std::size_t i = 0; i < bins.size(); ++i)
        assert(host_test::near(bins[i], 2 * original[i]));
      x(0, 0, 2) *= 4;
      rebinDeposition(*p, g.input(), 0, edges, b, t);
      assert(!host_test::near(b(0, 0, 0), 2 * original[0]));
    }
    if (modal) {
      x(0, 0, 1) = 8e-9;
      rebinDeposition(*p, g.input(), 0, edges, b, t);
      assert(host_test::near(t(0, 0, 0) + t(0, 0, 1) + b(0, 0, 0) + b(0, 0, 1) +
                                 b(0, 0, 2),
                             x(0, 0, 0)));
    }
    std::fill(g.state.begin(), g.state.end(), 0);
    rebinDeposition(*p, g.input(), 0, edges, b, t);
    assert(std::all_of(bins.begin(), bins.end(),
                       [](double v) { return v == 0; }) &&
           std::all_of(tails.begin(), tails.end(),
                       [](double v) { return v == 0; }));
    std::fill(bins.begin(), bins.end(), 77);
    std::fill(tails.begin(), tails.end(), 88);
    edge[2] = edge[1];
    host_test::rejects([&] { rebinDeposition(*p, g.input(), 0, edges, b, t); });
    assert(bins.front() == 77 && tails.back() == 88);
  }
}
