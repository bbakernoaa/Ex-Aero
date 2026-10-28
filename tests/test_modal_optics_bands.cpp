#include <algorithm>
#include <exaero/Environment.hpp>
#include <exaero/ModalCoupling.hpp>
#include <exaero/OpticsBandCoupling.hpp>
#include <exaero/PackageFactory.hpp>
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace {
using namespace exaero;

void check(bool ok, const char *message) {
  if (!ok) throw std::runtime_error(message);
}

class BandProvider final : public IModalBandOpticsProvider {
public:
  bool supports(const ModalConfig &, std::size_t,
                const RadiationBandInfo &band) const override {
    return band.source_id == 3;
  }
  ModalOpticalResult evaluateBand(
      const ModalConfig &, std::size_t, const RadiationBandInfo &band,
      double, double, std::span<const double> masses,
      std::span<const double>) const override {
    if (band.source_id != 3)
      throw std::invalid_argument("unsupported band");
    double total = 0;
    for (double mass : masses) total += mass;
    return {10 * total, 5 * total, .4, .1 * total, .1 * total,
            8 * total, 4 * total, 7 * total, 3.5 * total};
  }
};

ModalConfig config(bool internal) {
  ModalConfig cfg;
  cfg.coagulation = false;
  cfg.species = {{"acid", 1800., .098, .6, {}}};
  cfg.modes = {{"mode", 1.e-7, 1.5, {},
                internal ? ModalOpticalTreatment::InternalMixture
                         : ModalOpticalTreatment::ComponentTables}};
  return cfg;
}
}

void modal_native_band_provider_reports_support_and_computes_atomically() {
  auto package = createPackage("modal");
  initializeModal(*package, config(true));
  auto &coupling = modalCoupling(*package);
  coupling.setModalBandOpticsProvider(std::make_shared<BandProvider>());
  auto &bands = opticsBandCoupling(*package);
  const auto availability = bands.bandAvailability(RadiationBandSet::RRTMG_30);
  check(availability.size() == 1, "component availability size");
  check(availability[0].status[2] == RadiationBandStatus::Available,
        "supported band status");
  check(availability[0].status[1] == RadiationBandStatus::NotInSource,
        "unsupported band status");
  check(package->capabilities() & OpticsBands, "optics band capability");

  double t=280., p=90000., rho=1., rh=.5, dz=100.;
  EnvironmentalStateView env{View2D<const double>(&t,1,1),
    View2D<const double>(&p,1,1),View2D<const double>(&rho,1,1),
    View2D<const double>(&rh,1,1),View2D<const double>(&dz,1,1)};
  const std::array<double,1> masses{1.e-6};
  const double number=modalCoupling(*package).emittedNumber(0,masses);
  const double state_data[]={masses[0],number};
  View3D<const double> state(state_data,1,1,2);
  const int selected[]={2};
  std::array<double,optical_indices::NUM_OPTICS> result{};
  View4D<double> out(result.data(),1,1,1,optical_indices::NUM_OPTICS);
  bands.computeBandOptics(env,state,RadiationBandSet::RRTMG_30,selected,out);
  check(std::abs(out(0,0,0,optical_indices::EXTINCTION_COEFF)-1.e-5)<1.e-15,
        "integrated extinction");
  check(std::abs(out(0,0,0,optical_indices::EXTINCTION_AOT)-1.e-3)<1.e-13,
        "integrated AOT");

  std::fill(result.begin(),result.end(),123.);
  const int unsupported[]={1};
  bool threw=false;
  try { bands.computeBandOptics(env,state,RadiationBandSet::RRTMG_30,
                                unsupported,out); }
  catch (const std::exception &) { threw=true; }
  check(threw,"unsupported band must fail");
  for (double value : result) check(value == 123.,"failure must be atomic");
}

void modal_component_table_without_native_bands_is_unavailable() {
  auto package = createPackage("modal");
  auto cfg=config(false);
  ModalTableBinding binding;
  binding.kind=ModalTableKind::LognormalAverage;
  binding.geometric_std_dev=1.5;
  binding.curve.source_label="MODAL_TEST";
  binding.data_file=std::string(EXAERO_TEST_DATA_DIR)+"/mie_modal_reference.txt";
  binding.wavelength_nodes={5.e-7,1.e-6};
  cfg.species[0].optics=binding;
  initializeModal(*package, cfg);
  auto availability=opticsBandCoupling(*package).bandAvailability(
      RadiationBandSet::RRTMG_30);
  check(availability.size()==1,"component availability size");
  for (auto status:availability[0].status)
    check(status==RadiationBandStatus::NotInSource,"no generated band data");
  check(!(package->capabilities() & OpticsBands),"capability requires data");
}

void modal_native_band_table_uses_source_ids_and_native_coefficients() {
  auto cfg=config(false);
  ModalTableBinding binding;
  binding.kind=ModalTableKind::LognormalAverage;
  binding.geometric_std_dev=1.5;
  binding.curve.source_label="MODAL_TEST";
  binding.data_file=std::string(EXAERO_TEST_DATA_DIR)+"/mie_modal_band_reference.txt";
  cfg.species[0].optics=binding;
  auto package=createPackage("modal");
  initializeModal(*package,cfg);
  auto &bands=opticsBandCoupling(*package);
  const auto availability=bands.bandAvailability(RadiationBandSet::RRTMG_30);
  check(availability.size()==1,"native table component count");
  check(availability[0].status[2]==RadiationBandStatus::Available,
        "source band 3 availability");
  check(availability[0].status[4]==RadiationBandStatus::Available,
        "source band 5 availability");
  check(availability[0].status[3]==RadiationBandStatus::NotInSource,
        "absent source band unavailable");
  const std::array<double,1> masses{1.e-6};
  const double number=modalCoupling(*package).emittedNumber(0,masses);
  const double state_data[]={masses[0],number};
  double t=280.,p=90000.,rho=1.,rh=.5,dz=100.;
  EnvironmentalStateView env{View2D<const double>(&t,1,1),View2D<const double>(&p,1,1),
    View2D<const double>(&rho,1,1),View2D<const double>(&rh,1,1),View2D<const double>(&dz,1,1)};
  View3D<const double> state(state_data,1,1,2);
  const int selected[]={2};
  std::array<double,optical_indices::NUM_OPTICS> result{};
  View4D<double> out(result.data(),1,1,1,optical_indices::NUM_OPTICS);
  bands.computeBandOptics(env,state,RadiationBandSet::RRTMG_30,selected,out);
  check(std::abs(out(0,0,0,optical_indices::EXTINCTION_COEFF)-1.e-5)<1.e-15,
        "native component extinction");
  check(std::abs(out(0,0,0,optical_indices::SCATTERING_COEFF)-5.e-6)<1.e-15,
        "native component scattering");
  check(std::abs(out(0,0,0,optical_indices::ASYMMETRY_FACTOR)-.4)<1.e-14,
        "native component asymmetry");

  binding.kind=ModalTableKind::ParticleResolved;
  cfg.species[0].optics=binding;
  auto particle_package=createPackage("modal");
  initializeModal(*particle_package,cfg);
  const double particle_number=modalCoupling(*particle_package).emittedNumber(0,masses);
  const double particle_state_data[]={masses[0],particle_number};
  View3D<const double> particle_state(particle_state_data,1,1,2);
  std::array<double,optical_indices::NUM_OPTICS> particle_result{};
  View4D<double> particle_out(particle_result.data(),1,1,1,
                              optical_indices::NUM_OPTICS);
  opticsBandCoupling(*particle_package).computeBandOptics(
      env,particle_state,RadiationBandSet::RRTMG_30,selected,particle_out);
  check(std::abs(particle_out(0,0,0,optical_indices::EXTINCTION_COEFF)-1.e-5)<1.e-12,
        "particle-resolved native band integral");
}

int main() {
  initialize_environment();
  modal_native_band_provider_reports_support_and_computes_atomically();
  modal_component_table_without_native_bands_is_unavailable();
  modal_native_band_table_uses_source_ids_and_native_coefficients();
  finalize_environment();
}
