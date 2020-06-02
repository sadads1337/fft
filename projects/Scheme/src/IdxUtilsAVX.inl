#pragma once

#include <Core/Types.h>

#if FFT_ENABLE_MANUAL_VECT
#include <immintrin.h>
#else
#error "Must not be included if FFT_ENABLE_MANUAL_VECT is 0"
#endif

namespace utils {

inline namespace avx {

inline auto summ_real_impl(const Grid1D &x, const Grid1D &y, Grid1D &result) {
  static_assert(sizeof(Grid1D::value_type) == 8u);
  for (auto idx = 0u; idx < x.size() / 4u; ++idx) {
    auto x_data = _mm256_load_pd(x.data() + 4u * idx);
    auto y_data = _mm256_load_pd(y.data() + 4u * idx);

    const auto sum = _mm256_add_pd(x_data, y_data);
    _mm256_store_pd(result.data() + 4u * idx, sum);
  }
  for (auto idx = x.size() - x.size() % 4u; idx < x.size(); ++idx) {
    result[idx] = x[idx] + y[idx];
  }
}

inline auto apply_corr_factor_impl(const Grid1D &input, const Precision mult,
                                   const size_t offset, const size_t fg_count,
                                   Grid1D &result) {
  static_assert(sizeof(Grid1D::value_type) == 8u);
  for (auto idx = 0u; idx < input.size() / 4u; ++idx) {
    const auto k_idx = idx * 4u;

    auto input_data = _mm256_load_pd(input.data() + k_idx);
    const auto k_idx_real = offset + k_idx;
    const auto size_real = fg_count * input.size();
    const auto factors_data = _mm256_set_pd(
        mult * static_cast<Precision>(size_real - k_idx_real - 3u),
        mult * static_cast<Precision>(size_real - k_idx_real - 2u),
        mult * static_cast<Precision>(size_real - k_idx_real - 1u),
        mult * static_cast<Precision>(size_real - k_idx_real));

    const auto factorized = _mm256_mul_pd(input_data, factors_data);
    _mm256_store_pd(result.data() + k_idx, factorized);
  }

  for (auto k_idx = input.size() - input.size() % 4u; k_idx < input.size();
       ++k_idx) {
    const auto k_idx_real = offset + k_idx;
    const auto size_real = fg_count * input.size();
    result[k_idx] =
        mult * static_cast<Precision>(size_real - k_idx_real) * input[k_idx];
  }
}

inline auto apply_conv_factor_impl(const Grid1D &input, const Precision mult,
                                   const size_t offset, Grid1D &result) {
  static_assert(sizeof(Grid1D::value_type) == 8u);
  for (auto idx = 0u; idx < input.size() / 4u; ++idx) {
    const auto k_idx = idx * 4u;
    const auto k_idx_real = offset + k_idx;
    auto input_data = _mm256_load_pd(input.data() + k_idx);
    const auto factors_data =
        _mm256_set_pd(mult * static_cast<Precision>(k_idx_real + 3u),
                      mult * static_cast<Precision>(k_idx_real + 2u),
                      mult * static_cast<Precision>(k_idx_real + 1u),
                      mult * static_cast<Precision>(k_idx_real));

    const auto factorized = _mm256_mul_pd(input_data, factors_data);
    _mm256_store_pd(result.data() + k_idx, factorized);
  }

  for (auto k_idx = input.size() - input.size() % 4u; k_idx < input.size();
       ++k_idx) {
    const auto k_idx_real = offset + k_idx;
    result[k_idx] = mult * static_cast<Precision>(k_idx_real) * input[k_idx];
  }
}

inline auto apply_operation_impl(
    const Grid1D &lhs, const Grid1D &rhs,
    const std::function<Precision(Precision, Precision)> &function,
    Grid1D &result) {
  static_assert(sizeof(Grid1D::value_type) == 8u);
  for (auto idx = 0u; idx < lhs.size() / 4u; ++idx) {
    const auto function_apply_data =
        _mm256_set_pd(function(lhs[idx * 4u + 3u], rhs[idx * 4u + 3u]),
                      function(lhs[idx * 4u + 2u], rhs[idx * 4u + 2u]),
                      function(lhs[idx * 4u + 1u], rhs[idx * 4u + 1u]),
                      function(lhs[idx * 4u], rhs[idx * 4u]));
    _mm256_store_pd(result.data() + 4u * idx, function_apply_data);
  }
  for (auto idx = lhs.size() - lhs.size() % 4u; idx < lhs.size(); ++idx) {
    result[idx] = function(lhs[idx], rhs[idx]);
  }
}

}  // namespace avx

}  // namespace utils