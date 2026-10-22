/// @file ModalPresets.cpp
/// @brief Editable presets based on pinned MAM4xx mode/material metadata.
/// @author EX-aero contributors
/// @details History: 2026-10. In-memory configuration only; no file/device I/O.
#include <exaero/ModalPresets.hpp>

namespace exaero {
ModalConfig mam4ModalConfig() {
  ModalConfig c;
  // Preserve the pinned native constants when comparing diagnostics/CCN.
  c.constants.gravity = 9.80616;
  c.constants.water_molecular_weight = .018016;
  c.constants.gas_constant = 6.022214076e23 * 1.380649e-23;
  const char *names[] = {"SOA", "SO4", "POM", "BC", "NaCl", "DST", "MOM"};
  const char *labels[] = {"MAM4_SOA",  "MAM4_SO4", "MAM4_POM", "MAM4_BC",
                          "MAM4_NACL", "MAM4_DST", "MAM4_MOM"};
  const double density[] = {1000, 1770, 1000, 1700, 1900, 2600, 1601};
  const double kappa[] = {.1, .507, 1e-10, 1e-10, 1.16, .14, .1};
  // Explicit chemistry surrogates in kg/mol. MOM is 250.093 g/mol; native
  // coagulation/water do not consume molecular weights. MKPP owns chemistry.
  const double mw[] = {.0120107, .09606,  .0120107, .0120107,
                       .05844,   .135065, .250093};
  for (int s = 0; s < 7; ++s) {
    ModalTableBinding binding;
    binding.curve.source_label = labels[s];
    c.species.push_back({names[s], density[s], mw[s], kappa[s], binding});
  }
  c.modes = {{"accumulation", 5.5e-8, 1.8, {}},
             {"aitken", 1.3e-8, 1.6, {}},
             {"coarse", 1e-6, 1.8, {}},
             {"primary_carbon", 2.5e-8, 1.60000002384186, {}}};
  c.coagulation_routes = {{"accumulation", "aitken", "accumulation"},
                          {"accumulation", "coarse", "coarse"},
                          {"accumulation", "primary_carbon", "accumulation"},
                          {"aitken", "coarse", "coarse"},
                          {"aitken", "primary_carbon", "primary_carbon"},
                          {"coarse", "primary_carbon", "coarse"}};
  return c;
}
ModalConfig mam4xxReferenceModalConfig() {
  auto c = mam4ModalConfig();
  c.coagulation_model = ModalCoagulationModel::Mam4xx;
  c.water_uptake = ModalWaterUptakeModel::Mam4xx;
  c.coagulation_routes.clear();
  return c;
}
} // namespace exaero
