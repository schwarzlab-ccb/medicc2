// Saturating mandatory WGD. X (symbol 9) is handled separately by callers.
#pragma once
namespace wgd_model {
constexpr int double_cn(int cn) { return cn <= 4 ? 2 * cn : 8; }
} // namespace wgd_model
