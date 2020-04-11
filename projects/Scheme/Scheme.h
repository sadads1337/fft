#pragma once

#include <Core/GridUtils.h>
#include <Core/Types.h>

#include <type_traits>
#include <cmath>

namespace scheme
{

//! \todo: Сейчас все пресеты захардкожены, нужны вытащить их в параметры

//! \todo: Сейчас сетка от 0. до g_*_limit_value, поправить для произвольной. И ниже.
inline static constexpr auto g_z_limit_value = static_cast<Precision>(1.);
inline static constexpr auto g_t_limit_value = static_cast<Precision>(1.);

//! Схема не является безусловно устойчивой -> ограничения должны удовлетворять условиям Куранта.
//! \todo: fixme!!! constexpr check with static assertion
inline static constexpr auto g_t_grid_size = 1000u;
inline static constexpr auto g_z_grid_size = 500u;

inline static constexpr auto g_k_limit = 30u;

static_assert(std::is_same_v<std::remove_cv_t<decltype(g_z_limit_value)>, double>, "g_z_grid_size must be double");
static_assert(std::is_same_v<std::remove_cv_t<decltype(g_t_limit_value)>, double>, "g_t_limit_value must be double");

//! Care double value might be out of range.
inline static constexpr auto g_z_grid_step = g_z_limit_value / static_cast<Precision>(g_z_grid_size);
inline static constexpr auto g_t_grid_step = g_t_limit_value / static_cast<Precision>(g_t_grid_size);

static_assert(g_z_grid_step < static_cast<Precision>(1.) / (static_cast<Precision>(2.) * static_cast<Precision>(g_k_limit)));

inline static const auto g_pi = std::atan(static_cast<Precision>(1.)) * static_cast<Precision>(4.);

struct Values final
{
	Grid2D u = make_grid_2d(g_z_grid_size, g_k_limit);
	Grid2D w = make_grid_2d(g_z_grid_size, g_k_limit);
	Grid2D p = make_grid_2d(g_z_grid_size, g_k_limit);
	Grid2D q = make_grid_2d(g_z_grid_size, g_k_limit);
	Grid2D s = make_grid_2d(g_z_grid_size, g_k_limit);
};

struct Env final
{
	const Grid1D rho;
	const Grid1D lambda;
	const Grid1D mu;
	const Grid1D f;

	Env() = delete;
};

//! legacy c
Grid1D source(int IG, float WN7, float DT, float DZ, int K8);

Precision u_func(const Grid2D & u, const size_t x_idx, const size_t z_idx);

void calculate_one_step(const Values & prev_values, Values & values, const Env & env, const size_t t_idx);

} // namespace scheme
