#pragma once

#include <Core/MKL/Utils.h>
#include <Core/MakeWithCapacity.h>

#include "OpenMP.h"

#if FFT_ENABLE_MANUAL_VECT
#include "IdxUtilsAVX.inl"
#endif

namespace utils {

inline namespace omp {

inline auto summ_real_impl(const Grid1D& x, const Grid1D& y, Grid1D& result) {
  FFT_OMP_PRAGMA("omp simd")
  for (auto idx = 0u; idx < x.size(); ++idx) {
    result[idx] = x[idx] + y[idx];
  }
}

inline auto apply_corr_factor_impl(const Grid1D& input, const Precision mult,
                                   const size_t offset, const size_t fg_count,
                                   Grid1D& result) {
  FFT_OMP_PRAGMA("omp simd")
  for (auto k_idx = 0u; k_idx < input.size(); ++k_idx) {
    const auto k_idx_real = offset + k_idx;
    const auto size_real = fg_count * input.size();
    result[k_idx] =
        mult * static_cast<Precision>(size_real - k_idx_real) * input[k_idx];
  }
}

inline auto apply_conv_factor_impl(const Grid1D& input, const Precision mult,
                                   const size_t offset, Grid1D& result) {
  FFT_OMP_PRAGMA("omp simd")
  for (auto k_idx = 0u; k_idx < input.size(); ++k_idx) {
    const auto k_idx_real = offset + k_idx;
    result[k_idx] = mult * static_cast<Precision>(k_idx_real) * input[k_idx];
  }
}

inline auto apply_operation_impl(
    const Grid1D& lhs, const Grid1D& rhs,
    const std::function<Precision(Precision, Precision)>& function,
    Grid1D& result) {
  FFT_OMP_PRAGMA("omp simd")
  for (auto idx = 0u; idx < lhs.size(); ++idx) {
    result[idx] = function(lhs[idx], rhs[idx]);
  }
  return result;
}

}  // namespace omp

inline namespace normal {

inline auto conv_real_impl(const Grid1D &x, const Grid1D &y, Grid1D &result) {
  assert(x.size() + y.size() - 1u == result.size());
  for (auto i = 0u; i < x.size() + y.size() - 1u; i++) {
    result[i] = static_cast<Precision>(0.);
    const auto start = i >= x.size() ? i - x.size() + 1u : 0u;
    const auto end = i < y.size() ? i : y.size() - 1u;
    for (auto k = start; k <= end; k++) {
      result[i] += x[i - k] * y[k];
    }
  }
}

inline auto corr_real_impl(const Grid1D &x, const Grid1D &y, Grid1D &result) {
  assert(x.size() + y.size() - 1u == result.size());
  for (auto i = 0u; i < result.size(); i++) {
    result[i] = static_cast<Precision>(0.);
    const auto start = i >= x.size() ? i - x.size() + 1u : 0u;
    const auto end = i < y.size() ? i : y.size() - 1u;
    for (auto k = start; k <= end; k++) {
      result[i] += x[x.size() - 1u - (i - k)] * y[k];
    }
  }
}

inline auto summ_real_impl(const Grid1D& x, const Grid1D& y, Grid1D& result) {
  for (auto idx = 0u; idx < x.size(); ++idx) {
    result[idx] = x[idx] + y[idx];
  }
}

inline auto apply_corr_factor_impl(const Grid1D& input, const Precision mult,
                                   const size_t offset, const size_t fg_count,
                                   Grid1D& result) {
  for (auto k_idx = 0u; k_idx < input.size(); ++k_idx) {
    const auto k_idx_real = offset + k_idx;
    const auto size_real = fg_count * input.size();
    result[k_idx] =
        mult * static_cast<Precision>(size_real - k_idx_real) * input[k_idx];
  }
}

inline auto apply_conv_factor_impl(const Grid1D& input, const Precision mult,
                                   const size_t offset, Grid1D& result) {
  for (auto k_idx = 0u; k_idx < input.size(); ++k_idx) {
    const auto k_idx_real = offset + k_idx;
    result[k_idx] = mult * static_cast<Precision>(k_idx_real) * input[k_idx];
  }
}

inline auto apply_operation_impl(
    const Grid1D& lhs, const Grid1D& rhs,
    const std::function<Precision(Precision, Precision)>& function,
    Grid1D& result) {
  for (auto idx = 0u; idx < lhs.size(); ++idx) {
    result[idx] = function(lhs[idx], rhs[idx]);
  }
}

}  // namespace normal

template <bool vectorized>
inline auto apply_operation(
    const Grid1D& lhs, const Grid1D& rhs,
    const std::function<Precision(Precision, Precision)>& function) {
  assert(lhs.size() && rhs.size());
  assert(lhs.size() == rhs.size());
  Grid1D result(rhs.size());
  if constexpr (vectorized) {
#if FFT_ENABLE_MANUAL_VECT
    avx::apply_operation_impl(lhs, rhs, function, result);
#else
    omp::apply_operation_impl(lhs, rhs, function, result);
#endif
  } else {
    normal::apply_operation_impl(lhs, rhs, function, result);
  }
  return result;
}

template <bool vectorized>
inline auto apply_conv_factor(const Grid1D& input, const Precision mult,
                              const size_t offset) {
  //! \todo: remove redudant allocation
  Grid1D result(input.size());
  if constexpr (vectorized) {
#if FFT_ENABLE_MANUAL_VECT
    avx::apply_conv_factor_impl(input, mult, offset, result);
#else
    omp::apply_conv_factor_impl(input, mult, offset, result);
#endif
  } else {
    normal::apply_conv_factor_impl(input, mult, offset, result);
  }
  return result;
}

template <bool vectorized>
inline auto apply_corr_factor(const Grid1D& input, const Precision mult,
                              const size_t offset, const size_t fg_count) {
  //! \todo: remove redudant allocation
  Grid1D result(input.size());
  if constexpr (vectorized) {
#if FFT_ENABLE_MANUAL_VECT
    avx::apply_corr_factor_impl(input, mult, offset, fg_count, result);
#else
    omp::apply_corr_factor_impl(input, mult, offset, fg_count, result);
#endif
  } else {
    normal::apply_corr_factor_impl(input, mult, offset, fg_count, result);
  }
  return result;
}

template <bool vectorized>
auto inline conv_real(const Grid1D& x, const Grid1D& y) {
  Grid1D result(x.size() + y.size() - 1u);
  if constexpr (vectorized) {
#if FFT_ENABLE_MANUAL_VECT
    avx::conv_real_impl(x, y, result);
#else
    assert(false);
#endif
  }
  else {
    normal::conv_real_impl(x, y, result);
  }
  return result;
}

template <bool vectorized>
auto inline corr_real(const Grid1D& x, const Grid1D& y) {
  Grid1D result(x.size() + y.size() - 1u);
  if constexpr (vectorized) {
#if FFT_ENABLE_MANUAL_VECT
    avx::corr_real_impl(x, y, result);
#else
    assert(false);
#endif
  }
  else {
    normal::corr_real_impl(x, y, result);
  }
  return result;
}

template <bool vectorized>
auto inline summ_real(const Grid1D& x, const Grid1D& y) {
  assert(x.size() == y.size());
  //! \todo: remove redudant allocation
  Grid1D result(x.size());

  if constexpr (vectorized) {
#if FFT_ENABLE_MANUAL_VECT
    avx::summ_real_impl(x, y, result);
#else
    omp::summ_real_impl(x, y, result);
#endif
  } else {
    normal::summ_real_impl(x, y, result);
  }
  return result;
}

}  // namespace utils
