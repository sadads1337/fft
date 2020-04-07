#pragma once

#include "Types.h"

namespace utils
{

auto get_idx(const Precision l_edge, const Precision r_edge, const Precision value, const size_t grid_size);

} // namespace utils

#include "IdxUtils.inl"
