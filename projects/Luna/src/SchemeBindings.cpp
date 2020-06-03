#include <Core/MoveOnly.h>
#include <Core/Types.h>
#include <Luna/SchemeBindings.h>

#include <cassert>
#include <iostream>

namespace {

inline auto calculate_source() {
  return scheme::source(4, static_cast<Precision>(1.), scheme::g_t_grid_step,
                        scheme::g_z_grid_step, scheme::g_t_grid_size);
}

inline auto operator==(const scheme::Values& lhs,
                       const scheme::Values& rhs) noexcept {
  return lhs.u == rhs.u && lhs.w == rhs.w && lhs.p == rhs.p && lhs.q == rhs.q &&
         lhs.s == rhs.s;
}

inline auto operator!=(const scheme::Values& lhs,
                       const scheme::Values& rhs) noexcept {
  return !(lhs == rhs);
}

template <typename F, typename... T>
decltype(auto) no_omp_wrapper(F&& function, T&&... args) {
#if FFT_ENABLE_OPENMP
#error "OPENMP must be disabled while running with LuNA RTS"
#endif
  return std::forward<F>(function)(std::forward<T>(args)...);
}

}  // namespace

void init_env([[maybe_unused]] const std::int32_t fg_num,
              [[maybe_unused]] const std::int32_t fg_size,
              luna::ucenv::OutputDF& out_env) {
  //! Since our environment is constant \b fg_num and \b fg_size unused.
  //! \todo: Fix it later, when environment become non constant.

  const scheme::Env env{
      Grid1D(scheme::g_k_limit, static_cast<Precision>(1.)),
      Grid1D(scheme::g_k_limit, static_cast<Precision>(1.)),
      Grid1D(scheme::g_k_limit, static_cast<Precision>(1.)),
      utils::move_only(calculate_source()),
  };

  out_env.set(env);
}

void check_env(const luna::ucenv::InputDF& env,
               [[maybe_unused]] const std::int32_t fg_num,
               [[maybe_unused]] const std::int32_t fg_size) {
  const auto env_data = env.getData<scheme::Env>();
  assert(env_data);
  assert(env_data->rho.size() == scheme::g_k_limit);
  assert(env_data->lambda.size() == scheme::g_k_limit);
  assert(env_data->mu.size() == scheme::g_k_limit);
  assert(env_data->f.size() == scheme::g_t_grid_size);

  const auto check_values = [](const auto& values, const auto& expected_value) {
    for (const auto& value : values) {
      if (expected_value != value) {
        std::cout << "Failed check in: " << typeid(decltype(values)).name();
      }
    }
  };

  check_values(env_data->rho, static_cast<Precision>(1.));
  check_values(env_data->lambda, static_cast<Precision>(1.));
  check_values(env_data->mu, static_cast<Precision>(1.));

  if (env_data->f != utils::move_only(calculate_source())) {
    std::cout << std::string{__func__} + " failed check in: "
              << typeid(env_data->f).name() << std::endl;
  }
}

void init_values([[maybe_unused]] const std::int32_t fg_num,
                 const std::int32_t fg_size, luna::ucenv::OutputDF& values) {
  scheme::Values values_data{
      make_grid_2d(fg_size, scheme::g_k_limit),
      make_grid_2d(fg_size, scheme::g_k_limit),
      make_grid_2d(fg_size, scheme::g_k_limit),
      make_grid_2d(fg_size, scheme::g_k_limit),
      make_grid_2d(fg_size, scheme::g_k_limit),
  };

  values.set(utils::move_only(values_data));
}

//! import check_prev_values(value #prev_values, int #fgnum, int #size) as
//! check_prev_values; import check_values(value #values, int #fgnum, int #size)
//! as check_values; \note scheme::g_z_grid_size = fg_num * fg_size; \warning
//! only for debug purposes;
void check_values(const luna::ucenv::InputDF& values,
                  [[maybe_unused]] const std::int32_t fg_num,
                  [[maybe_unused]] const std::int32_t fg_size) {
  const auto values_data = values.getData<scheme::Values>();
  assert(values_data);
  assert(values_data->u.size() == static_cast<size_t>(fg_size));
  assert(values_data->w.size() == static_cast<size_t>(fg_size));
  assert(values_data->p.size() == static_cast<size_t>(fg_size));
  assert(values_data->q.size() == static_cast<size_t>(fg_size));
  assert(values_data->s.size() == static_cast<size_t>(fg_size));

  const auto check_values = [](const auto& values, const auto& expected_value) {
    for (const auto& value : values) {
      for (const auto& el : value) {
        if (expected_value != el) {
          std::cout << std::string{__func__} + " failed check in: "
                    << typeid(decltype(values)).name() << std::endl;
        }
      }
    }
  };

  check_values(values_data->u, static_cast<Precision>(0.));
  check_values(values_data->w, static_cast<Precision>(0.));
  check_values(values_data->p, static_cast<Precision>(0.));
  check_values(values_data->q, static_cast<Precision>(0.));
  check_values(values_data->s, static_cast<Precision>(0.));
}

void calculate_one_step(const luna::ucenv::InputDF& prev_values,
                        const luna::ucenv::InputDF& env,
                        const std::int32_t t_idx, const std::int32_t fg_num,
                        const std::int32_t fg_size, const std::int32_t fg_count,
                        luna::ucenv::OutputDF& values) {
  const auto prev_values_data = prev_values.getData<scheme::Values>();
  assert(prev_values_data);
  assert(prev_values_data->u.size() == static_cast<size_t>(fg_size));
  assert(prev_values_data->w.size() == static_cast<size_t>(fg_size));
  assert(prev_values_data->p.size() == static_cast<size_t>(fg_size));
  assert(prev_values_data->q.size() == static_cast<size_t>(fg_size));
  assert(prev_values_data->s.size() == static_cast<size_t>(fg_size));
  const auto env_data = env.getData<scheme::Env>();
  assert(env_data);
  assert(env_data->rho.size() == scheme::g_k_limit);
  assert(env_data->lambda.size() == scheme::g_k_limit);
  assert(env_data->mu.size() == scheme::g_k_limit);
  assert(env_data->f.size() == scheme::g_t_grid_size);

  //! \todo: fix redudant copy
  scheme::Values values_data{
      make_grid_2d(fg_size, scheme::g_k_limit),
      make_grid_2d(fg_size, scheme::g_k_limit),
      make_grid_2d(fg_size, scheme::g_k_limit),
      make_grid_2d(fg_size, scheme::g_k_limit),
      make_grid_2d(fg_size, scheme::g_k_limit),
  };

  no_omp_wrapper(scheme::calculate_one_step, *prev_values_data, values_data,
                 *env_data, t_idx, fg_num, fg_size, fg_count);

  assert(values_data != *prev_values_data);

  values.set(values_data);
}

void check_one_step(const luna::ucenv::InputDF& prev_values,
                    const luna::ucenv::InputDF& values,
                    const luna::ucenv::InputDF& env, const std::int32_t t_idx,
                    const std::int32_t fg_num, const std::int32_t fg_size,
                    const std::int32_t fg_count) {
  scheme::Values expected_values_data{
      make_grid_2d(fg_size, scheme::g_k_limit),
      make_grid_2d(fg_size, scheme::g_k_limit),
      make_grid_2d(fg_size, scheme::g_k_limit),
      make_grid_2d(fg_size, scheme::g_k_limit),
      make_grid_2d(fg_size, scheme::g_k_limit),
  };

  const auto prev_values_data = prev_values.getData<scheme::Values>();
  assert(prev_values_data);
  assert(prev_values_data->u.size() == static_cast<size_t>(fg_size));
  assert(prev_values_data->w.size() == static_cast<size_t>(fg_size));
  assert(prev_values_data->p.size() == static_cast<size_t>(fg_size));
  assert(prev_values_data->q.size() == static_cast<size_t>(fg_size));
  assert(prev_values_data->s.size() == static_cast<size_t>(fg_size));
  const auto values_data = values.getData<scheme::Values>();
  assert(values_data);
  assert(values_data->u.size() == static_cast<size_t>(fg_size));
  assert(values_data->w.size() == static_cast<size_t>(fg_size));
  assert(values_data->p.size() == static_cast<size_t>(fg_size));
  assert(values_data->q.size() == static_cast<size_t>(fg_size));
  assert(values_data->s.size() == static_cast<size_t>(fg_size));
  const auto env_data = env.getData<scheme::Env>();
  assert(env_data);
  assert(env_data->rho.size() == scheme::g_k_limit);
  assert(env_data->lambda.size() == scheme::g_k_limit);
  assert(env_data->mu.size() == scheme::g_k_limit);
  assert(env_data->f.size() == scheme::g_t_grid_size);

  //! Since \b prev_values are correct, we check only next step. Induction works
  //! here.
  no_omp_wrapper(scheme::calculate_one_step, *prev_values_data,
                 expected_values_data, *env_data, t_idx, fg_num, fg_size,
                 fg_count);

  if (*values_data != expected_values_data) {
    std::cout << std::string{__func__} + " failed check in: "
              << typeid(decltype(values)).name() << " step: " << t_idx
              << std::endl;
  }
}
