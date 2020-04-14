#pragma once

#include <Core/Types.h>

namespace utils {

template <bool vectorized>
auto apply_operation(
    const Grid1D& lhs, const Grid1D& rhs,
    const std::function<Precision(Precision, Precision)>& function);

template <bool vectorized>
auto apply_conv_factor(const Grid1D& input, Precision mult);

template <bool vectorized>
auto apply_corr_factor(const Grid1D& input, Precision mult);

template <bool vectorized>
auto summ_real(const Grid1D& x, const Grid1D& y);

}  // namespace utils

#include "IdxUtils.inl"
