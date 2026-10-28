#ifdef NDEBUG
#undef NDEBUG
#endif
#include <exaero/OpticsBandCoupling.hpp>
#include <exaero/RadiationBands.hpp>
#include <cassert>
#include <cmath>
#include <stdexcept>

int main() {
  const auto bands = exaero::radiationBandInfo(exaero::RadiationBandSet::RRTMG_30);
  assert(bands.size() == 30);
  for (int i = 0; i < 30; ++i) {
    assert(bands[static_cast<std::size_t>(i)].source_id == i + 1);
    assert(bands[static_cast<std::size_t>(i)].region ==
           (i < 14 ? exaero::RadiationBandRegion::Shortwave
                   : exaero::RadiationBandRegion::Longwave));
    assert(std::isfinite(bands[static_cast<std::size_t>(i)].lower_wavenumber_cm1));
    assert(std::isfinite(bands[static_cast<std::size_t>(i)].upper_wavenumber_cm1));
    assert(0 < bands[static_cast<std::size_t>(i)].lower_wavenumber_cm1);
    assert(bands[static_cast<std::size_t>(i)].lower_wavenumber_cm1 <
           bands[static_cast<std::size_t>(i)].upper_wavenumber_cm1);
  }
  const auto all = exaero::selectRadiationBands(exaero::RadiationBandSet::RRTMG_30, {});
  assert(all.size() == 30);
  const int requested[] = {14, 2};
  const auto selected = exaero::selectRadiationBands(
      exaero::RadiationBandSet::RRTMG_30, requested);
  assert(selected.size() == 2 && selected[0].source_id == 15 && selected[1].source_id == 3);
  for (const auto invalid : {std::vector<int>{-1}, std::vector<int>{30},
                             std::vector<int>{2, 2}}) {
    bool threw = false;
    try { (void)exaero::selectRadiationBands(exaero::RadiationBandSet::RRTMG_30, invalid); }
    catch (const std::invalid_argument &) { threw = true; }
    assert(threw);
  }
  bool threw = false;
  try { (void)exaero::radiationBandInfo(static_cast<exaero::RadiationBandSet>(99)); }
  catch (const std::invalid_argument &) { threw = true; }
  assert(threw);
}
