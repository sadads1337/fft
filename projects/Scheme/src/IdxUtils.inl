#pragma once

#include <Core/MKL/Utils.h>
#include <Core/MakeWithCapacity.h>

#include "OpenMP.h"

namespace utils {

template <bool vectorized>
inline auto apply_operation(
    const Grid1D& lhs, const Grid1D& rhs,
    const std::function<Precision(Precision, Precision)>& function) {
  assert(lhs.size() && rhs.size());
  assert(lhs.size() == rhs.size());
  Grid1D result(rhs.size());
  if constexpr (vectorized) {
    FFT_OMP_PRAGMA("omp simd")
    for (auto idx = 0u; idx < lhs.size(); ++idx) {
      result[idx] = function(lhs[idx], rhs[idx]);
    }
  } else {
    for (auto idx = 0u; idx < lhs.size(); ++idx) {
      result[idx] = function(lhs[idx], rhs[idx]);
    }
  }
  return result;
}

template <bool vectorized>
inline auto apply_conv_factor(const Grid1D& input, const Precision mult) {
  //! \todo: remove redudant allocation
  Grid1D result(input.size());
  if constexpr (vectorized) {
    FFT_OMP_PRAGMA("omp simd")
    for (auto k_idx = 0u; k_idx < input.size(); ++k_idx) {
      result[k_idx] = mult * static_cast<Precision>(k_idx) * input[k_idx];
    }
  } else {
    for (auto k_idx = 0u; k_idx < input.size(); ++k_idx) {
      result[k_idx] = mult * static_cast<Precision>(k_idx) * input[k_idx];
    }
  }
  return result;
}

template <bool vectorized>
inline auto apply_corr_factor(const Grid1D& input, const Precision mult) {
  //! \todo: remove redudant allocation
  Grid1D result(input.size());
  if constexpr (vectorized) {
    FFT_OMP_PRAGMA("omp simd")
    for (auto k_idx = 0u; k_idx < input.size(); ++k_idx) {
      result[k_idx] =
          mult * static_cast<Precision>(input.size() - k_idx) * input[k_idx];
    }
  } else {
    for (auto k_idx = 0u; k_idx < input.size(); ++k_idx) {
      result[k_idx] =
          mult * static_cast<Precision>(input.size() - k_idx) * input[k_idx];
    }
  }
  return result;
}

template <bool vectorized>
auto inline summ_real(const Grid1D& x, const Grid1D& y) {
  assert(x.size() == y.size());
  //! \todo: remove redudant allocation
  Grid1D result(x.size());

  if constexpr (vectorized) {
    FFT_OMP_PRAGMA("omp simd")
    for (auto idx = 0u; idx < x.size(); ++idx) {
      result[idx] = x[idx] + y[idx];
    }
  } else {
    for (auto idx = 0u; idx < x.size(); ++idx) {
      result[idx] = x[idx] + y[idx];
    }
  }
  return result;
}

}  // namespace utils
