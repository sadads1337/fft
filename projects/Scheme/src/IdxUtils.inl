#pragma once

#include <Core/MKL/Utils.h>
#include <Core/MakeWithCapacity.h>

#include "OpenMP.h"

#if FFT_ENABLE_MANUAL_VECT
#include <immintrin.h>
#endif

namespace utils {

inline namespace avx
{

inline auto summ_real_impl(const Grid1D& x, const Grid1D& y, Grid1D & result)
{
  static_assert(sizeof(Grid1D::value_type) == 8u);
  for (auto idx = 0u; idx < x.size() / 4u; ++idx)
  {
    inline auto x_data = _mm256_load_pd(x.data() + 4u * idx);
    inline auto y_data = _mm256_load_pd(y.data() + 4u * idx);

    inline const auto sum = _mm256_add_pd(x_data, y_data);
    _mm256_store_pd(result.data() + 4u * idx, sum);
  }
  for (auto idx = x.size() - x.size() % 4u; idx < x.size(); ++idx)
  {
    result[idx] = x[idx] + y[idx];
  }
}

inline auto apply_corr_factor_impl(const Grid1D& input, const Precision mult, Grid1D & result)
{
  static_assert(sizeof(Grid1D::value_type) == 8u);
  for (auto idx = 0u; idx < input.size() / 4u; ++idx)
  {
    const auto k_idx = idx * 4u;

    inline auto input_data = _mm256_load_pd(input.data() + k_idx);
    inline const auto factors_data = _mm256_set_pd(
        mult * static_cast<Precision>(input.size() - k_idx - 3u),
        mult * static_cast<Precision>(input.size() - k_idx - 2u),
        mult * static_cast<Precision>(input.size() - k_idx - 1u),
        mult * static_cast<Precision>(input.size() - k_idx));

    inline const auto factorized = _mm256_mul_pd(input_data, factors_data);
    _mm256_store_pd(result.data() + k_idx, factorized);
  }

  for (auto k_idx = input.size() - input.size() % 4u; k_idx < input.size(); ++k_idx) {
    result[k_idx] =
        mult * static_cast<Precision>(input.size() - k_idx) * input[k_idx];
  }
}

inline auto apply_conv_factor_impl(const Grid1D& input, const Precision mult, Grid1D & result)
{
  static_assert(sizeof(Grid1D::value_type) == 8u);
  for (auto idx = 0u; idx < input.size() / 4u; ++idx)
  {
    inline const auto k_idx = idx * 4u;

    inline auto input_data = _mm256_load_pd(input.data() + k_idx);
    inline const auto factors_data = _mm256_set_pd(
        mult * static_cast<Precision>(k_idx + 3u),
        mult * static_cast<Precision>(k_idx + 2u),
        mult * static_cast<Precision>(k_idx + 1u),
        mult * static_cast<Precision>(k_idx));

    inline const auto factorized = _mm256_mul_pd(input_data, factors_data);
    _mm256_store_pd(result.data() + k_idx, factorized);
  }

  for (auto k_idx = input.size() - input.size() % 4u; k_idx < input.size(); ++k_idx) {
    result[k_idx] =
        mult * static_cast<Precision>(k_idx) * input[k_idx];
  }
}

inline auto apply_operation_impl(
    const Grid1D& lhs, const Grid1D& rhs,
    const std::function<Precision(Precision, Precision)>& function,
    Grid1D & result) {
  static_assert(sizeof(Grid1D::value_type) == 8u);
  for (auto idx = 0u; idx < lhs.size() / 4u; ++idx)
  {
    inline const auto function_apply_data = _mm256_set_pd(
        function(lhs[idx * 4u + 3u], rhs[idx * 4u + 3u]),
        function(lhs[idx * 4u + 2u], rhs[idx * 4u + 2u]),
        function(lhs[idx * 4u + 1u], rhs[idx * 4u + 1u]),
        function(lhs[idx * 4u], rhs[idx * 4u]));
    _mm256_store_pd(result.data() + 4u * idx, function_apply_data);
  }
  for (auto idx = lhs.size() - lhs.size() % 4u; idx < lhs.size(); ++idx)
  {
    result[idx] = function(lhs[idx], rhs[idx]);
  }
}

} // inline namespace avx

inline namespace omp
{

inline auto summ_real_impl(const Grid1D& x, const Grid1D& y, Grid1D & result)
{
  FFT_OMP_PRAGMA("omp simd")
  for (auto idx = 0u; idx < x.size(); ++idx)
  {
    result[idx] = x[idx] + y[idx];
  }
}

inline auto apply_corr_factor_impl(const Grid1D& input, const Precision mult, Grid1D & result) {
  FFT_OMP_PRAGMA("omp simd")
  for (auto k_idx = 0u; k_idx < input.size(); ++k_idx) {
    result[k_idx] =
        mult * static_cast<Precision>(input.size() - k_idx) * input[k_idx];
  }
}

inline auto apply_conv_factor_impl(const Grid1D& input, const Precision mult, Grid1D & result) {
  FFT_OMP_PRAGMA("omp simd")
  for (auto k_idx = 0u; k_idx < input.size(); ++k_idx) {
    result[k_idx] =
        mult * static_cast<Precision>(k_idx) * input[k_idx];
  }
}

inline auto apply_operation_impl(
    const Grid1D& lhs, const Grid1D& rhs,
    const std::function<Precision(Precision, Precision)>& function,
    Grid1D & result) {
  FFT_OMP_PRAGMA("omp simd")
  for (auto idx = 0u; idx < lhs.size(); ++idx) {
    result[idx] = function(lhs[idx], rhs[idx]);
  }
  return result;
}

} // inline namespace omp

inline namespace normal
{

inline auto summ_real_impl(const Grid1D& x, const Grid1D& y, Grid1D & result)
{
  for (auto idx = 0u; idx < x.size(); ++idx)
  {
    result[idx] = x[idx] + y[idx];
  }
}

inline auto apply_corr_factor_impl(const Grid1D& input, const Precision mult, Grid1D & result) {
  for (auto k_idx = 0u; k_idx < input.size(); ++k_idx) {
    result[k_idx] =
        mult * static_cast<Precision>(input.size() - k_idx) * input[k_idx];
  }
}

inline auto apply_conv_factor_impl(const Grid1D& input, const Precision mult, Grid1D & result) {
  for (auto k_idx = 0u; k_idx < input.size(); ++k_idx) {
    result[k_idx] =
        mult * static_cast<Precision>(k_idx) * input[k_idx];
  }
}

inline auto apply_operation_impl(
    const Grid1D& lhs, const Grid1D& rhs,
    const std::function<Precision(Precision, Precision)>& function,
    Grid1D & result) {
  for (auto idx = 0u; idx < lhs.size(); ++idx) {
    result[idx] = function(lhs[idx], rhs[idx]);
  }
}

} // inline namespace normal

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
inline auto apply_conv_factor(const Grid1D& input, const Precision mult) {
  //! \todo: remove redudant allocation
  Grid1D result(input.size());
  if constexpr (vectorized) {
#if FFT_ENABLE_MANUAL_VECT
    avx::apply_conv_factor_impl(input, mult, result);
#else
    omp::apply_conv_factor_impl(input, mult, result);
#endif
  } else {
    normal::apply_conv_factor_impl(input, mult, result);
  }
  return result;
}

template <bool vectorized>
inline auto apply_corr_factor(const Grid1D& input, const Precision mult) {
  //! \todo: remove redudant allocation
  Grid1D result(input.size());
  if constexpr (vectorized) {
#if FFT_ENABLE_MANUAL_VECT
    avx::apply_corr_factor_impl(input, mult, result);
#else
    omp::apply_corr_factor_impl(input, mult, result);
#endif
  } else {
    normal::apply_corr_factor_impl(input, mult, result);
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
