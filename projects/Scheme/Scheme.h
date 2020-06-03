#pragma once

#include <Core/GridUtils.h>
#include <Core/Types.h>

#include <cmath>
#include <functional>
#include <type_traits>

namespace scheme {

//! \todo: Сейчас все пресеты захардкожены, нужны вытащить их в параметры

//! \todo: Сейчас сетка от 0. до g_*_limit_value, поправить для произвольной. И
//! ниже.
inline static constexpr auto g_z_limit_value = static_cast<Precision>(1.);
inline static constexpr auto g_t_limit_value = static_cast<Precision>(1.);

//! Схема не является безусловно устойчивой -> ограничения должны удовлетворять
//! условиям Куранта. \todo: fixme!!! constexpr check with static assertion
inline static constexpr auto g_t_grid_size = 2048u;
inline static constexpr auto g_z_grid_size = 1024u;

inline static constexpr auto g_k_limit = 32u;

static_assert(
    std::is_same_v<std::remove_cv_t<decltype(g_z_limit_value)>, double>,
    "g_z_grid_size must be double");
static_assert(
    std::is_same_v<std::remove_cv_t<decltype(g_t_limit_value)>, double>,
    "g_t_limit_value must be double");

//! Care double value might be out of range.
inline static constexpr auto g_z_grid_step =
    g_z_limit_value / static_cast<Precision>(g_z_grid_size);
inline static constexpr auto g_t_grid_step =
    g_t_limit_value / static_cast<Precision>(g_t_grid_size);

static_assert(g_z_grid_step <
              static_cast<Precision>(1.) / (static_cast<Precision>(2.) *
                                            static_cast<Precision>(g_k_limit)));

inline static const auto g_pi =
    std::atan(static_cast<Precision>(1.)) * static_cast<Precision>(4.);

struct Values final {
  Grid2D u = make_grid_2d(g_z_grid_size, g_k_limit);
  Grid2D w = make_grid_2d(g_z_grid_size, g_k_limit);
  Grid2D p = make_grid_2d(g_z_grid_size, g_k_limit);
  Grid2D q = make_grid_2d(g_z_grid_size, g_k_limit);
  Grid2D s = make_grid_2d(g_z_grid_size, g_k_limit);
};

struct Env final {
  Grid1D rho;
  Grid1D lambda;
  Grid1D mu;
  Grid1D f;

  //! \todo: maybe constant initilized, but there is a bug in LuNA df
  //! initialization code.
};

//! legacy c
Grid1D source(int IG, Precision WN7, Precision DT, Precision DZ, int K8);

Precision u_func(const Grid2D& u, size_t x_idx, size_t z_idx);

void calculate_one_step(const Values& prev_values, Values& values,
                        const Env& env, const size_t t_idx, const size_t fg_num,
                        const size_t fg_size, const size_t fg_count,
                        size_t num_threads);

void main_loop_for_t(Values& prev_values, Values& values, const Env& env,
                     size_t t_idx_limit,
                     const std::function<void(const Values&)>& callback,
                     size_t num_threads);

}  // namespace scheme
