#pragma once
/// @file ModalPresets.hpp
/// @brief Editable MAM4 modal profiles with explicit process choices.
/// @author EX-aero contributors
/// @details History: 2026-10. Return host configuration only; no file/device
/// I/O.
#include <exaero/ModalConfig.hpp>

namespace exaero {
/// @brief Four native MAM4 modes/materials with generic Brownian/kappa physics.
/// All seven mass slots permit coatings; this is not a complete MAM4 scheme.
/// Available without MAM4xx. The returned configuration is freely editable.
ModalConfig mam4ModalConfig();
/// @brief Native coagulation/water reference profile, 32 interstitial fields.
/// Initialization requires EXAERO_ENABLE_MAM4XX. Native coagulation fixes mode
/// widths/material densities/kappa; generic processes permit arbitrary changes.
ModalConfig mam4xxReferenceModalConfig();
} // namespace exaero
