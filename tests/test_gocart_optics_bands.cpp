#ifdef NDEBUG
#undef NDEBUG
#endif
#include <exaero/AerosolIndices.hpp>
#include <exaero/Environment.hpp>
#include <exaero/OpticsBandCoupling.hpp>
#include <gocart/GocartPackage.hpp>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>

namespace {
exaero::GocartSpeciesConfig species(const char *name) {
  exaero::GocartSpeciesConfig s;
  s.name = name;
  s.dry_density = 1800.; s.molecular_weight = 98.;
  s.dry_particle_diameter = 2.e-7; s.hygroscopicity = .2;
  s.lognormal_sigma = 1.5; s.lognormal_dg = 1.5e-7;
  s.refractive_index_real = 1.5; s.refractive_index_imag = .01;
  return s;
}
exaero::EnvironmentalStateView environment(double &t, double &p, double &rho,
                                           double &rh, double &dz) {
  return {exaero::View2D<const double>(&t, 1, 1),
          exaero::View2D<const double>(&p, 1, 1),
          exaero::View2D<const double>(&rho, 1, 1),
          exaero::View2D<const double>(&rh, 1, 1),
          exaero::View2D<const double>(&dz, 1, 1)};
}
void close(double actual, double expected) {
  assert(std::abs(actual - expected) <= 1.e-10 * std::max(1., std::abs(expected)));
}
}

int main() {
  using namespace exaero;
  initialize_environment();
  // Native band values for one source species must agree with the table query.
  {
    GocartConfig cfg; cfg.species.push_back(species("DU"));
    GocartPackage pkg; pkg.initialize(cfg);
    auto &bands = opticsBandCoupling(pkg);
    const int selected[] = {2};
    double t=280., p=90000., rho=1., rh=.5, dz=120., mass=1.e-6;
    auto env=environment(t,p,rho,rh,dz);
    View3D<const double> state(&mass,1,1,1);
    double result[optical_indices::NUM_OPTICS]{};
    View4D<double> out(result,1,1,1,optical_indices::NUM_OPTICS);
    bands.computeBandOptics(env,state,RadiationBandSet::RRTMG_30,selected,out);
    double bext=0., bsca=0., bbck=0., g=0.; ProvenanceInfo provenance{};
    pkg.queryAttribute(0,AttributeCategory::SpectralOptical,
      spectral_optical_indices::MASS_EXTINCTION,rh,3.,&bext,&provenance);
    pkg.queryAttribute(0,AttributeCategory::SpectralOptical,
      spectral_optical_indices::MASS_SCATTERING,rh,3.,&bsca,&provenance);
    pkg.queryAttribute(0,AttributeCategory::SpectralOptical,
      spectral_optical_indices::MASS_BACKSCATTER,rh,3.,&bbck,&provenance);
    pkg.queryAttribute(0,AttributeCategory::SpectralOptical,
      spectral_optical_indices::ASYMMETRY_FACTOR,rh,3.,&g,&provenance);
    close(out(0,0,0,optical_indices::EXTINCTION_COEFF),mass*bext);
    close(out(0,0,0,optical_indices::SCATTERING_COEFF),mass*bsca);
    close(out(0,0,0,optical_indices::BACKSCATTER_COEFF),mass*bbck);
    close(out(0,0,0,optical_indices::LIDAR_BACKSCATTER),mass*bbck);
    close(out(0,0,0,optical_indices::ASYMMETRY_FACTOR),g);
    close(out(0,0,0,optical_indices::EXTINCTION_AOT),mass*bext*dz);
    double wavelength=550.e-9;
    View1D<const double> wavelengths(&wavelength,1);
    double spectral_result[optical_indices::NUM_OPTICS]{};
    View4D<double> spectral_out(spectral_result,1,1,1,optical_indices::NUM_OPTICS);
    pkg.computeOptics(env,state,wavelengths,spectral_out);
    assert(std::isfinite(spectral_out(0,0,0,optical_indices::EXTINCTION_COEFF)));
  }
  // All baked band species are supported; brown carbon is explicit missing data.
  {
    GocartConfig cfg;
    for (const char *name : {"DU","SS","SU","BC","OC","NI","BR"})
      cfg.species.push_back(species(name));
    GocartPackage pkg; pkg.initialize(cfg);
    const auto availability=opticsBandCoupling(pkg).bandAvailability(RadiationBandSet::RRTMG_30);
    assert(availability.size()==7);
    for (std::size_t i=0;i<6;++i)
      for (auto status:availability[i].status) assert(status==RadiationBandStatus::Available);
    for (auto status:availability[6].status) assert(status==RadiationBandStatus::NotInSource);
  }
  // Unsupported positive mass throws before touching any output element.
  {
    GocartConfig cfg; cfg.species={species("DU"),species("BR")};
    GocartPackage pkg; pkg.initialize(cfg);
    double t[]={280.,280.}, p[]={90000.,90000.}, rho[]={1.,1.}, rh[]={.5,.5}, dz[]={120.,120.};
    double mass[]={1.e-6,0.,0.,1.e-6};
    EnvironmentalStateView env{View2D<const double>(t,2,1),View2D<const double>(p,2,1),
      View2D<const double>(rho,2,1),View2D<const double>(rh,2,1),View2D<const double>(dz,2,1)};
    View3D<const double> state(mass,2,1,2);
    int selected[]={2};
    std::vector<double> result(2*optical_indices::NUM_OPTICS,123.); View4D<double> out(result.data(),2,1,1,optical_indices::NUM_OPTICS);
    bool threw=false;
    try { opticsBandCoupling(pkg).computeBandOptics(env,state,RadiationBandSet::RRTMG_30,selected,out); }
    catch (const std::runtime_error &) { threw=true; }
    assert(threw); for (double x:result) assert(x==123.);
  }
  // Shape and selector validation happen before any caller output write.
  {
    GocartConfig cfg; cfg.species.push_back(species("DU"));
    GocartPackage pkg; pkg.initialize(cfg);
    double t=280., p=90000., rho=1., rh=.5, dz=120., mass=1.e-6;
    auto env=environment(t,p,rho,rh,dz); View3D<const double> state(&mass,1,1,1);
    std::vector<double> result(optical_indices::NUM_OPTICS,321.);
    View4D<double> wrong_shape(result.data(),1,1,1,optical_indices::NUM_OPTICS-1);
    bool threw=false;
    try { opticsBandCoupling(pkg).computeBandOptics(env,state,RadiationBandSet::RRTMG_30,{},wrong_shape); }
    catch (const std::invalid_argument &) { threw=true; }
    assert(threw); for (double x:result) assert(x==321.);
  }
  // Unsupported zero-mass BR does not block DU and its source result.
  {
    GocartConfig cfg; cfg.species={species("DU"),species("BR")};
    GocartPackage pkg; pkg.initialize(cfg);
    double t=280., p=90000., rho=1., rh=.5, dz=120.; double mass[]={1.e-6,0.};
    auto env=environment(t,p,rho,rh,dz); View3D<const double> state(mass,1,1,2);
    int selected[]={2}; double result[optical_indices::NUM_OPTICS]{}; View4D<double> out(result,1,1,1,optical_indices::NUM_OPTICS);
    opticsBandCoupling(pkg).computeBandOptics(env,state,RadiationBandSet::RRTMG_30,selected,out);
    assert(std::isfinite(out(0,0,0,optical_indices::EXTINCTION_COEFF)) && result[0]>0.);
  }
  // A second package initialization must not mutate the first package's curves.
  {
    GocartConfig first; first.species.push_back(species("DU"));
    GocartPackage pkg_a; pkg_a.initialize(first);
    GocartConfig second; second.species.push_back(species("SS"));
    GocartPackage pkg_b; pkg_b.initialize(second);
    const auto a = opticsBandCoupling(pkg_a).bandAvailability(RadiationBandSet::RRTMG_30);
    assert(a.size()==1 && a[0].status[1]==RadiationBandStatus::Available);
    double t=280., p=90000., rho=1., rh=.5, dz=120., mass=1.e-6;
    auto env=environment(t,p,rho,rh,dz); View3D<const double> state(&mass,1,1,1);
    int selected[]={2}; double result[optical_indices::NUM_OPTICS]{};
    View4D<double> out(result,1,1,1,optical_indices::NUM_OPTICS);
    opticsBandCoupling(pkg_a).computeBandOptics(env,state,RadiationBandSet::RRTMG_30,selected,out);
    assert(result[optical_indices::EXTINCTION_COEFF]>0.);
  }
  finalize_environment();
}
