#include <Core/MKL/Utils.h>
#include <Core/TimeGuard.h>
#include <Scheme/Scheme.h>

#include "IdxUtils.h"
#include "OpenMP.h"

namespace scheme {

Params create_params(const size_t num_threads, const Precision z_limit_value, const size_t z_grid_size,
                     const size_t t_limit, const Precision t_grid_step, const size_t k_limit) {
  const auto z_grid_step = z_limit_value / static_cast<Precision>(z_grid_size);

  /*assert(z_grid_step <
                static_cast<Precision>(1) / (static_cast<Precision>(2.) *
                                              static_cast<Precision>(k_limit)));
  assert(t_grid_step <= z_grid_step / sqrt(2.));*/
  return {
      num_threads,
      z_limit_value,
      z_grid_size,
      z_grid_step,
      t_limit,
      t_grid_step,
      k_limit,
  };
}

//! \todo: fixme скопипащенный legacy c
Grid1D source(int IG, Precision WN7, Precision DT, Precision DZ, int K8) {
  auto F = Grid1D(static_cast<size_t>(K8));
  Precision C, B, T, T1, T2, T3, T4, PII, X, A, S;
  int N, GI;
  PII = std::acos(static_cast<Precision>(-1.));
  int NPT1 = K8;
  GI = IG;
  C = DT / DZ;
  B = static_cast<Precision>(0.);
  T1 = static_cast<Precision>(1.5);
  T2 = static_cast<Precision>(3.5);
  T3 = static_cast<Precision>(5.5);
  T = static_cast<Precision>(2.) * PII * WN7;

  if (GI <= static_cast<Precision>(2.01))
    T4 = T1;
  else {
    if ((GI <= static_cast<Precision>(4.01)) &&
        (GI >= static_cast<Precision>(2.01)))
      T4 = T2;
    else
      T4 = T3;
  }
  N = int(T4 / (DT * WN7));
  X = T4 / (static_cast<Precision>(2.) * WN7);
  A = (T / GI) * (T / GI);
  S = -X;
  for (auto J = 0; J < NPT1; ++J) {
    if (J + 1 <= N + 1)
      F[J] = C * exp((-A) * (S * S)) * cos(T * S + B);
    else
      F[J] = static_cast<Precision>(0.);
    S = S + DT;
  }
  for (auto i = 0; i < NPT1; ++i) F[i] = F[i] * DT / (DZ);

  return F;
}

Precision u_func(const Grid2D& u, const size_t x_idx, const size_t z_idx, const Params& params) {
  const auto b = params.z_limit_value;
  auto result = static_cast<Precision>(0.);
  for (auto k_idx = 0u; k_idx < params.k_limit; ++k_idx) {
    const auto x = static_cast<Precision>(x_idx) * params.z_grid_step;
    result += u[z_idx][k_idx] *
              std::sin(static_cast<Precision>(k_idx) * g_pi / b * x);
  }
  return result * (static_cast<Precision>(2.) / b);
}

template <typename... T>
auto inline summ_real(T&&... args) {
  constexpr bool vectorized = false;
  return utils::summ_real<vectorized>(std::forward<T>(args)...);
}

template <typename... T>
auto inline apply_conv_factor(T&&... args) {
  constexpr bool vectorized = false;
  return utils::apply_conv_factor<vectorized>(std::forward<T>(args)...);
}

template <typename... T>
auto inline apply_corr_factor(T&&... args) {
  constexpr bool vectorized = false;
  return utils::apply_corr_factor<vectorized>(std::forward<T>(args)...);
}

template <typename... T>
auto inline apply_operation(T&&... args) {
  constexpr bool vectorized = false;
  return utils::apply_operation<vectorized>(std::forward<T>(args)...);
}

template <typename... T>
auto inline conv_real(T&&... args) {
  constexpr bool vectorized = false;
  return utils::conv_real<vectorized>(std::forward<T>(args)...);
}

template <typename... T>
auto inline corr_real(T&&... args) {
  constexpr bool vectorized = false;
  return utils::corr_real<vectorized>(std::forward<T>(args)...);
}

void calculate_one_step(const Values& prev_values, Values& values,
                        const Env& env, const size_t t_idx, const size_t fg_num,
                        const size_t fg_size, const size_t fg_count,
                        const Params& params) {
  const auto& prev_u = prev_values.u;
  const auto& prev_w = prev_values.w;
  const auto& prev_p = prev_values.p;
  const auto& prev_q = prev_values.q;
  const auto& prev_s = prev_values.s;

  auto& u = values.u;
  auto& w = values.w;
  auto& p = values.p;
  auto& q = values.q;
  auto& s = values.s;

  const auto& rho = env.rho;
  const auto& lambda = env.lambda;
  const auto& mu = env.mu;
  const auto& f = env.f;

  const auto z_idx_start = fg_num == 0u ? 1u : 0u;
  const auto z_idx_end = fg_size - 1u;

  const auto offset = fg_num * fg_size;

  assert(z_idx_start >= 0);
  assert(z_idx_end < params.z_grid_size);
  assert(fg_size * fg_count == params.z_grid_size);

#if FFT_ENABLE_OPENMP
  omp_set_dynamic(0);
  omp_set_num_threads(params.num_threads);
#endif

  FFT_OMP_PRAGMA("omp parallel for")
  for (auto z_idx = z_idx_start; z_idx < z_idx_end; ++z_idx) {
    //! For u
    const auto der_q = apply_operation(prev_q[z_idx + 1u], prev_q[z_idx],
                                       [&](const auto& lhs, const auto& rhs) {
                                         return (lhs - rhs) / params.z_grid_step;
                                       });
    const auto q_rho_for_u = summ_real(
        conv_real(apply_conv_factor(der_q, g_pi / params.z_limit_value, offset),
                       rho),
        corr_real(
            apply_corr_factor(der_q, g_pi / params.z_limit_value, offset, fg_count),
            rho));
    const auto p_rho_for_u = summ_real(
        conv_real(
            apply_conv_factor(prev_p[z_idx], g_pi / params.z_limit_value, offset),
            rho),
        corr_real(apply_corr_factor(prev_p[z_idx], g_pi / params.z_limit_value,
                                         offset, fg_count),
                       rho));

    //! For w
    const auto q_rho_for_w = summ_real(
        conv_real(apply_conv_factor(prev_q[z_idx + 1u],
                                         g_pi / params.z_limit_value, offset),
                       rho),
        corr_real(
            apply_corr_factor(prev_q[z_idx + 1u], g_pi / params.z_limit_value,
                              offset, fg_count),
            rho));
    const auto der_s = apply_operation(prev_s[z_idx + 1u], prev_s[z_idx],
                                       [&](const auto& lhs, const auto& rhs) {
                                         return (lhs - rhs) / params.z_grid_step;
                                       });
    const auto s_rho_for_w = summ_real(
        conv_real(apply_conv_factor(der_s, g_pi / params.z_limit_value, offset),
                       rho),
        corr_real(
            apply_corr_factor(der_s, g_pi / params.z_limit_value, offset, fg_count),
            rho));

    //! For p
    const auto der_w = apply_operation(prev_w[z_idx + 1u], prev_w[z_idx],
                                       [&](const auto& lhs, const auto& rhs) {
                                         return (lhs - rhs) / params.z_grid_step;
                                       });
    const auto w_lambda_for_p = summ_real(
        conv_real(apply_conv_factor(der_w, g_pi / params.z_limit_value, offset),
                       lambda),
        corr_real(
            apply_corr_factor(der_w, g_pi / params.z_limit_value, offset, fg_count),
            lambda));
    const auto summ_lambda_mu =
        apply_operation(lambda, mu, [](const auto& lhs, const auto& rhs) {
          return lhs + static_cast<Precision>(2.) * rhs;
        });
    const auto u_summ_lambda_mu_for_p = summ_real(
        conv_real(
            apply_conv_factor(prev_u[z_idx], g_pi / params.z_limit_value, offset),
            summ_lambda_mu),
        corr_real(apply_corr_factor(prev_u[z_idx], g_pi / params.z_limit_value,
                                         offset, fg_count),
                       summ_lambda_mu));

    //! For q
    const auto der_u = apply_operation(prev_u[z_idx + 1u], prev_u[z_idx],
                                       [&](const auto& lhs, const auto& rhs) {
                                         return (lhs - rhs) / params.z_grid_step;
                                       });
    const auto u_mu_for_q = summ_real(
        conv_real(apply_conv_factor(der_u, g_pi / params.z_limit_value, offset),
                       mu),
        corr_real(
            apply_corr_factor(der_u, g_pi / params.z_limit_value, offset, fg_count),
            mu));
    const auto w_mu_for_q = summ_real(
        conv_real(apply_conv_factor(prev_w[z_idx + 1u],
                                         g_pi / params.z_limit_value, offset),
                       mu),
        corr_real(
            apply_corr_factor(prev_w[z_idx + 1u], g_pi / params.z_limit_value,
                              offset, fg_count),
            mu));

    //! For s
    const auto w_summ_lambda_mu_for_s = summ_real(
        conv_real(apply_conv_factor(der_w, g_pi / params.z_limit_value, offset),
                       summ_lambda_mu),
        corr_real(
            apply_corr_factor(der_w, g_pi / params.z_limit_value, offset, fg_count),
            summ_lambda_mu));
    const auto u_lambda_for_s = summ_real(
        conv_real(
            apply_conv_factor(prev_u[z_idx], g_pi / params.z_limit_value, offset),
            lambda),
        corr_real(apply_corr_factor(prev_u[z_idx], g_pi / params.z_limit_value,
                                         offset, fg_count),
                       lambda));

    //! For f
    const auto z_0_idx = params.z_grid_size / 2u;
    const auto z_0 = static_cast<Precision>(z_0_idx) * params.z_grid_step;
    const auto f_x_h = (fg_num * fg_size + z_idx) == z_0_idx ? f[t_idx] : 0.;

    for (auto k_idx = 0u; k_idx < params.k_limit; ++k_idx) {
      constexpr auto half = static_cast<Precision>(0.5);

      u[z_idx][k_idx] =
          prev_u[z_idx][k_idx] +
          half * params.t_grid_step * (q_rho_for_u[k_idx] - p_rho_for_u[k_idx]);
      w[z_idx + 1u][k_idx] =
          prev_w[z_idx + 1u][k_idx] +
          half * params.t_grid_step * (q_rho_for_w[k_idx] + s_rho_for_w[k_idx]);
      p[z_idx][k_idx] =
          prev_p[z_idx][k_idx] +
          half * params.t_grid_step *
              (w_lambda_for_p[k_idx] - u_summ_lambda_mu_for_p[k_idx]) +
          f_x_h * std::cos(static_cast<Precision>(offset + k_idx) * g_pi /
                           params.z_limit_value * static_cast<Precision>(z_0));
      q[z_idx + 1u][k_idx] =
          prev_q[z_idx + 1u][k_idx] +
          half * params.t_grid_step * (u_mu_for_q[k_idx] - w_mu_for_q[k_idx]);
      s[z_idx][k_idx] =
          prev_s[z_idx][k_idx] +
          half * params.t_grid_step *
              (w_summ_lambda_mu_for_s[k_idx] - u_lambda_for_s[k_idx]) +
          f_x_h * std::cos(static_cast<Precision>(offset + k_idx) * g_pi /
                           params.z_limit_value * static_cast<Precision>(z_0));
    }
  }
}

void main_loop_for_t(Values& prev_values, Values& values, const Env& env,
                     const std::function<void(const Values&)>& callback,
                     const Params& params) {
  utils::TimeGuard time_guard{};
  for (auto t_idx = 0u; t_idx < params.t_limit; ++t_idx) {
    constexpr auto fg_num = 0u;
    const auto fg_size = params.z_grid_size;
    constexpr auto fg_count = 1u;
    calculate_one_step(prev_values, values, env, t_idx, fg_num, fg_size,
                       fg_count, params);

    //! new in values_2 now; we don't need values_2, swap here, do not copy
    std::swap(prev_values, values);

    if (callback) {
      callback(prev_values);
    }
  }
}

}  // namespace scheme
