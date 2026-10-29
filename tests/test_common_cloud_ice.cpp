#ifdef NDEBUG
#undef NDEBUG
#endif
#include "host_coupling_fixture.hpp"
#include <exaero/CloudIceCoupling.hpp>
#include <exaero/CloudIceCoupling_C.h>
#include <exaero/HostCoupling.hpp>
#include <limits>
int main() {
  host_test::Runtime runtime;
  using namespace exaero;
  for (bool modal : {false, true}) {
    auto p = host_test::make(modal);
    host_test::Grid g(modal);
    const int n = g.nc * g.nl * g.np, props = cloud_ice_properties::COUNT;
    std::vector<double> a(n * props), b(a.size()), refa(a.size()),
        refb(a.size()), cloud = g.state;
    for (auto &v : cloud)
      v *= .3;
    View3D<const double> cv(cloud.data(), g.nc, g.nl, g.nf);
    std::vector<int> as(n), bs(n), ras(n), rbs(n);
    View4D<double> av(a.data(), g.nc, g.nl, g.np, props),
        bv(b.data(), g.nc, g.nl, g.np, props),
        rav(refa.data(), g.nc, g.nl, g.np, props),
        rbv(refb.data(), g.nc, g.nl, g.np, props);
    View3D<int> asi(as.data(), g.nc, g.nl, g.np),
        bsi(bs.data(), g.nc, g.nl, g.np), rasi(ras.data(), g.nc, g.nl, g.np),
        rbsi(rbs.data(), g.nc, g.nl, g.np);
    auto input = g.input();
    const auto env = g.env();
    cloudIce(*p, env, input, &cv, IceCategory::Dust, false, av, bv, asi, bsi);
    if (modal)
      modalCloudIce(*p, env, input, &cv, IceCategory::Dust, false, rav, rbv,
                    rasi, rbsi);
    else
      bulkCloudIce(*p, env, input, &cv, IceCategory::Dust, false, rav, rbv,
                   rasi, rbsi);
    assert(a == refa && b == refb && as == ras && bs == rbs &&
           host_test::near(bv(0, 0, 0, 0), .3 * av(0, 0, 0, 0)));
    cloudIce(*p, env, input, nullptr, IceCategory::Dust, false, av, bv, asi,
             bsi);
    assert(std::all_of(bs.begin(), bs.end(),
                       [](int v) { return v == int(CloudIceStatus::Absent); }));
    std::fill(cloud.begin(), cloud.end(), 0);
    cloudIce(*p, env, input, &cv, IceCategory::Dust, false, av, bv, asi, bsi);
    assert(std::all_of(bs.begin(), bs.end(),
                       [](int v) { return v == int(CloudIceStatus::Empty); }));
    double ss[3]{0, .001, .01};
    std::vector<double> f(n * 6), rf(f.size());
    View1D<const double> ssv(ss, 3);
    View4D<double> fv(f.data(), g.nc, g.nl, g.np, 6),
        rfv(rf.data(), g.nc, g.nl, g.np, 6);
    activationFractions(*p, env, input, ssv, fv);
    if (modal)
      modalActivationFractions(*p, env, input, ssv, rfv);
    else
      bulkActivationFractions(*p, env, input, ssv, rfv);
    assert(f == rf);
    int flag;
    char error[256];
    exaero_activation_fractions(p.get(), g.nc, g.nl, g.nf, g.np, 3,
                                g.temp.data(), g.rh.data(), g.state.data(), ss,
                                f.data(), error, &flag);
    assert(flag == 0 && f == rf);
    std::fill(a.begin(), a.end(), 77);
    std::fill(b.begin(), b.end(), 88);
    exaero_cloud_ice(p.get(), g.nc, g.nl, g.nf, g.np, g.temp.data(),
                     g.rh.data(), g.state.data(), nullptr, 999, 0, a.data(),
                     b.data(), as.data(), bs.data(), error, &flag);
    assert(flag != 0 && a.front() == 77 && b.back() == 88);
    host_test::rejects([&] {
      cloudIce(*p, env, input, nullptr, IceCategory::Dust, false, av, av, asi,
               bsi);
    });
    assert(a.front() == 77 && b.back() == 88);
    assert(supportsHostOperation(*p, HostOperation::CloudIce) &&
           supportsHostOperation(*p, HostOperation::Activation));
  }
}
