#include <Core/MKL/Utils.h>
#include <Core/MakeWithCapacity.h>
#include <Core/MoveAndClear.h>
#include <Core/MoveOnly.h>

#include <iostream>
#include <array>
#include <type_traits>
#include <functional>

//! \todo: Студия не дружит с cassert, но очень нужно обмазать код assert'ми

using Precision = double;
static_assert(std::is_floating_point_v<Precision>, "Precision must be float or double");
using Grid1D = std::vector<double>;
using Grid2D = std::vector<Grid1D>;
using Grid3D = std::vector<Grid2D>;

auto make_grid(const size_t t, const size_t z, const size_t k, const double default_value = 0.)
{
	return Grid3D(t, Grid2D(z, Grid1D(k, default_value)));
}

//! \todo: Сейчас сетка от 0. до g_*_limit_value, поправить для произвольной. И ниже.
constexpr const auto g_z_limit_value = 1.;
constexpr const auto g_t_limit_value = 1.;

//! Схема не является безусловно устойчивой -> ограничения должны удовлетворять условиям Куранта.
//! \todo: fixme!!!
constexpr const auto g_t_grid_size = 10u;
constexpr const auto g_z_grid_size = 10u;

constexpr const auto g_k_limit = 10u;

static_assert(std::is_same<std::remove_cv_t<decltype(g_z_limit_value)>, double>::value, "g_z_grid_size must be double");
static_assert(std::is_same<std::remove_cv_t<decltype(g_t_limit_value)>, double>::value, "g_t_limit_value must be double");
//! Care double value might be out of range.
constexpr const auto g_z_grid_step = g_z_limit_value / static_cast<double>(g_z_grid_size);
constexpr const auto g_t_grid_step = g_t_limit_value / static_cast<double>(g_t_grid_size);

auto apply_conv_factor(const Grid1D & input, const size_t k)
{
	auto result = input;
	for (auto idx = 0u; idx <= input.size(); ++idx)
	{
		if (0u <= k - idx < input.size())
		{
			result[k - idx] = (k - idx) * result[k - idx];
		}
	}
	return utils::move_and_clear(result);
}

auto apply_corr_factor(const Grid1D & input, const size_t k)
{
	auto result = input;
	for (auto idx = 0u; idx <= input.size(); ++idx)
	{
		if (0u <= k - idx < input.size())
		{
			result[k - idx] = (k + idx) * result[k - idx];
		}
	}
	return utils::move_and_clear(result);
}

auto apply_operation(
	const Grid1D & lhs,
	const Grid1D & rhs,
	const std::function<Precision(Precision, Precision)> & function)
{
	auto result = utils::make_with_capacity<fft::RealContainer>(lhs.size());
	auto it_x = lhs.begin();
	auto it_y = rhs.begin();
	for (; it_x != lhs.end(); ++it_x, ++it_y)
	{
		result.emplace_back(function(*it_x, *it_y));
	}
	return result;
}

//! legacy c
auto f(int IG, float WN7, float DT, float DZ, int K8)
{
	auto F = Grid1D(static_cast<size_t>(K8));
	float C, B, T, T1, T2, T3, T4, PII, X, A, S, DX;
	int J, N, GI;
	float t0;
	PII = acos(-1.0);
	int NPT1 = K8;
	GI = IG;
	C = DT / DZ;
	B = 0.00;
	T1 = 1.500;
	T2 = 3.500;
	T3 = 5.500;
	T = 2.00*PII*WN7;

	if (GI <= 2.0100) T4 = T1;
	else {
		if ((GI <= 4.0100) && (GI >= 2.0100)) T4 = T2;
		else T4 = T3;
	}
	N = int(T4 / (DT*WN7));
	X = T4 / (2.00*WN7);
	A = (T / GI)*(T / GI);
	S = -X;
	for (J = 0; J<NPT1; J++) {
		if (J + 1 <= N + 1) F[J] = C * exp((-A)*(S*S))*cos(T*S + B);
		else F[J] = 0.00;
		S = S + DT;
	}
	for (int i = 0; i<NPT1; i++)
		F[i] = F[i] * DT / (DZ);

	return F;
}

int main() try
{
	auto u = make_grid(g_t_grid_size, g_z_grid_size, g_k_limit);
	auto w = make_grid(g_t_grid_size, g_z_grid_size, g_k_limit);
	auto p = make_grid(g_t_grid_size, g_z_grid_size, g_k_limit);
	auto q = make_grid(g_t_grid_size, g_z_grid_size, g_k_limit);
	auto s = make_grid(g_t_grid_size, g_z_grid_size, g_k_limit);

	const auto rho = Grid1D(g_z_grid_size, 1.);
	const auto lambda = Grid1D(g_z_grid_size, 1.);
	const auto mu = Grid1D(g_z_grid_size, 1.);

	for(auto t_idx = 0u; t_idx < g_t_grid_size; t_idx += 2)
	{
		for (auto z_idx = 2u; z_idx < g_z_grid_size; z_idx += 2)
		{
			for (auto k_idx = 0u; k_idx < g_k_limit; k_idx += 2)
			{
				const auto diffence_q = apply_operation(
					q[t_idx][z_idx + 1],
					q[t_idx][z_idx - 1],
					[](const auto & lhs, const auto & rhs) { return (lhs - rhs) / g_z_grid_step; });
				const auto q_rho_for_u = fft::summ_real(
					fft::conv_real(apply_conv_factor(diffence_q, k_idx), rho),
					fft::corr_real(apply_corr_factor(diffence_q, k_idx), rho));
				const auto p_rho_for_u = fft::summ_real(
					fft::conv_real(apply_conv_factor(p[t_idx][z_idx], k_idx), rho),
					fft::corr_real(apply_corr_factor(p[t_idx][z_idx], k_idx), rho));
				u[t_idx + 2u][z_idx][k_idx] = g_t_grid_step * u[t_idx][z_idx][k_idx] 
					+ g_t_grid_step * 0.5 * (q_rho_for_u[k_idx] - p_rho_for_u[k_idx]);

				const auto q_rho_for_w = fft::summ_real(
					fft::conv_real(apply_conv_factor(q[t_idx][z_idx + 1], k_idx), rho),
					fft::corr_real(apply_corr_factor(q[t_idx][z_idx + 1], k_idx), rho));
				const auto diffence_s = apply_operation(
					s[t_idx][z_idx + 2],
					s[t_idx][z_idx],
					[](const auto & lhs, const auto & rhs) { return (lhs - rhs) / g_z_grid_step; });
				const auto s_rho_for_u = fft::summ_real(
					fft::conv_real(apply_conv_factor(diffence_s, k_idx), rho),
					fft::corr_real(apply_corr_factor(diffence_s, k_idx), rho));
				w[t_idx + 2u][z_idx + 1u][k_idx] = g_t_grid_step * w[t_idx][z_idx + 1u][k_idx]
					+ g_t_grid_step * 0.5 * (s_rho_for_u[k_idx] + s_rho_for_u[k_idx]);

				const auto diffence_w = apply_operation(
					s[t_idx][z_idx + 1],
					s[t_idx][z_idx - 1],
					[](const auto & lhs, const auto & rhs) { return (lhs - rhs) / g_z_grid_step; });
				const auto w_lambda_for_p = fft::summ_real(
					fft::conv_real(apply_conv_factor(diffence_w, k_idx), rho),
					fft::corr_real(apply_corr_factor(diffence_w, k_idx), rho));
				const auto summ_lambda_mu = apply_operation(
					lambda,
					mu,
					[](const auto & lhs, const auto & rhs) { return lhs + 2. * rhs; });
				const auto u_summ_lambda_mu_for_p = fft::summ_real(
					fft::conv_real(apply_conv_factor(u[t_idx][z_idx], k_idx), summ_lambda_mu),
					fft::corr_real(apply_corr_factor(u[t_idx][z_idx], k_idx), summ_lambda_mu));
				const auto f_x_h = z_idx == 2u
					? f(4, 0.1, g_t_grid_step, g_z_grid_step, g_k_limit)
					: Grid1D(g_k_limit, 0.);
				p[t_idx + 2u][z_idx][k_idx] = g_t_grid_step * p[t_idx][z_idx][k_idx]
					+ g_t_grid_step * 0.5 * (w_lambda_for_p[k_idx] - u_summ_lambda_mu_for_p[k_idx])
					+ f_x_h[k_idx];

				const auto diffence_u = apply_operation(
					u[t_idx][z_idx + 2],
					u[t_idx][z_idx],
					[](const auto & lhs, const auto & rhs) { return (lhs - rhs) / g_z_grid_step; });
				const auto u_mu_for_q = fft::summ_real(
					fft::conv_real(apply_conv_factor(diffence_u, k_idx), mu),
					fft::corr_real(apply_corr_factor(diffence_u, k_idx), mu));
				const auto w_mu_for_q = fft::summ_real(
					fft::conv_real(apply_conv_factor(w[t_idx][z_idx + 1], k_idx), mu),
					fft::corr_real(apply_corr_factor(w[t_idx][z_idx + 1], k_idx), mu));
				q[t_idx + 2u][z_idx + 1u][k_idx] = g_t_grid_step * q[t_idx][z_idx + 1u][k_idx]
					+ g_t_grid_step * 0.5 * (u_mu_for_q[k_idx] - w_mu_for_q[k_idx]);

				const auto w_lambda_for_s = fft::summ_real(
					fft::conv_real(apply_conv_factor(diffence_w, k_idx), u_summ_lambda_mu_for_p),
					fft::corr_real(apply_corr_factor(diffence_w, k_idx), u_summ_lambda_mu_for_p));
				const auto u_lambda_for_s = fft::summ_real(
					fft::conv_real(apply_conv_factor(u[t_idx][z_idx], k_idx), lambda),
					fft::corr_real(apply_corr_factor(u[t_idx][z_idx], k_idx), lambda));
				s[t_idx + 2u][z_idx][k_idx] = g_t_grid_step * s[t_idx][z_idx][k_idx]
					+ g_t_grid_step * 0.5 * (w_lambda_for_s[k_idx] - u_lambda_for_s[k_idx])
					+ f_x_h[k_idx];
			}
		}
	}
}
catch(const utils::MKLException & exception)
{
	std::cout << "MKL exception happend: " << exception.what();
}
catch(const std::exception & exception)
{
	std::cout << "Another exception happend: " << exception.what();
}