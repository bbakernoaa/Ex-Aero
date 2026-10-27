#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <cmath>
#include <exaero/ModalRemovalMath.hpp>

int main() {
  const double sigma = 1.5;
  const double number = 1e8;
  const double radius = 1e-7;
  const double volume = exaero::modalNumberVolume(radius, sigma) * number;
  assert(std::abs(exaero::modalRadius(volume, number, sigma) / radius - 1) <
         1e-14);
  assert(
      std::abs(exaero::modalNumberAtRadius(volume, 2e-7, sigma) / (number / 8) -
               1) < 1e-14);
  assert(std::abs(exaero::modalMassFractionBelow(
                      2 * radius * std::exp(3 * std::pow(std::log(sigma), 2)),
                      radius, sigma) -
                  .5) < 1e-14);
  assert(exaero::modalMassFractionBelow(0, radius, sigma) == 0);
  assert(exaero::modalMassFractionBelow(INFINITY, radius, sigma) == 1);
  assert(exaero::modalMassFractionAbove(INFINITY, radius, sigma) == 0);
  assert(exaero::modalMassFractionAbove(0, radius, sigma) == 1);
  assert(exaero::modalMassFractionAbove(100 * radius, radius, sigma) > 0);
}
