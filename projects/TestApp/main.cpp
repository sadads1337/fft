#include <Core/MKL/Utils.h>
#include <Core/MakeWithCapacity.h>
#include <Core/MoveAndClear.h>
#include <Core/MoveOnly.h>
#include <Core/matplotlib-cpp/matplotlibcpp.h>
#include <Core/Types.h>

#include <iostream>
#include <array>
#include <type_traits>
#include <functional>

//! \todo: Студия не дружит с cassert, но очень нужно обмазать код assert'ми
#include <cassert>

auto make_grid(const size_t t, const size_t z, const size_t k, const double default_value = 0.)
{
	assert(t && z && k);
	return Grid3D(t, Grid2D(z, Grid1D(k, default_value)));
}

//! \todo: Сейчас сетка от 0. до g_*_limit_value, поправить для произвольной. И ниже.
constexpr const auto g_z_limit_value = 1.;
constexpr const auto g_t_limit_value = 1.;

//! Схема не является безусловно устойчивой -> ограничения должны удовлетворять условиям Куранта.
//! \todo: fixme!!! constexpr check with static assertion
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
	auto result = Grid1D(input.size(), 0.);
	for (auto idx = 0u; idx <= input.size(); ++idx)
	{
		if (0u <= k - idx && k - idx < input.size())
		{
			result[k - idx] = (k - idx) * input[k - idx];
		}
	}
	return result;
}

auto apply_corr_factor(const Grid1D & input, const size_t k)
{
	auto result = Grid1D(input.size(), 0.);
	for (auto idx = 0u; idx <= input.size(); ++idx)
	{
		if (0u <= k - idx && k - idx < input.size())
		{
			result[k - idx] = (k + idx) * input[k - idx];
		}
	}
	return result;
}

auto apply_operation(
	const Grid1D & lhs,
	const Grid1D & rhs,
	const std::function<Precision(Precision, Precision)> & function)
{
	assert(lhs.size() && rhs.size());
	auto result = utils::make_with_capacity<fft::RealContainer>(lhs.size());
	auto it_x = lhs.begin();
	auto it_y = rhs.begin();
	assert(lhs.size() == rhs.size());
	for (; it_x != lhs.end(); ++it_x, ++it_y)
	{
		result.emplace_back(function(*it_x, *it_y));
	}
	assert(result.size());
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

auto u_func(const Grid3D & u, const size_t x_idx, const size_t z_idx, const size_t t_idx)
{
	const auto b = 1.;
	static_assert(std::is_same_v<std::remove_cv_t<decltype(b)>, Precision>, "b must be the same type with main precision type");
	auto result = 2. / b;
	for (auto k_idx = 0; k_idx < g_k_limit; ++k_idx)
	{
		result += u[t_idx][z_idx][k_idx]
			* std::sin(static_cast<float>(k_idx) * x_idx * g_z_grid_step);
	}
	return result;
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

	for(auto t_idx = 0u; t_idx < g_t_grid_size - 2; t_idx += 2)
	{
		for (auto z_idx = 2u; z_idx < g_z_grid_size - 2u; z_idx += 2)
		{
			for (auto k_idx = 0u; k_idx < g_k_limit; ++k_idx)
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
				u[t_idx + 2u][z_idx][k_idx] = u[t_idx][z_idx][k_idx]
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
				w[t_idx + 2u][z_idx + 1u][k_idx] = w[t_idx][z_idx + 1u][k_idx]
					+ g_t_grid_step * 0.5 * (s_rho_for_u[k_idx] + s_rho_for_u[k_idx]);

				const auto diffence_w = apply_operation(
					s[t_idx][z_idx + 1],
					s[t_idx][z_idx - 1],
					[](const auto & lhs, const auto & rhs) { return (lhs - rhs) / g_z_grid_step; });
				const auto w_lambda_for_p = fft::summ_real(
					fft::conv_real(apply_conv_factor(diffence_w, k_idx), lambda),
					fft::corr_real(apply_corr_factor(diffence_w, k_idx), lambda));
				const auto summ_lambda_mu = apply_operation(
					lambda,
					mu,
					[](const auto & lhs, const auto & rhs) { return lhs + 2. * rhs; });
				const auto u_summ_lambda_mu_for_p = fft::summ_real(
					fft::conv_real(apply_conv_factor(u[t_idx][z_idx], k_idx), summ_lambda_mu),
					fft::corr_real(apply_corr_factor(u[t_idx][z_idx], k_idx), summ_lambda_mu));
				const auto f_x_h = z_idx == (g_z_grid_size / 2)
					? f(4, 10., g_t_grid_step, g_z_grid_step, g_t_grid_size)
					: Grid1D(g_t_grid_size, 0.);
				p[t_idx + 2u][z_idx][k_idx] = p[t_idx][z_idx][k_idx]
					+ g_t_grid_step * 0.5 * (w_lambda_for_p[k_idx] - u_summ_lambda_mu_for_p[k_idx])
					+ f_x_h[t_idx] * std::cos(k_idx * g_z_limit_value / 2.);

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
				q[t_idx + 2u][z_idx + 1u][k_idx] = q[t_idx][z_idx + 1u][k_idx]
					+ g_t_grid_step * 0.5 * (u_mu_for_q[k_idx] - w_mu_for_q[k_idx]);

				const auto w_lambda_for_s = fft::summ_real(
					fft::conv_real(apply_conv_factor(diffence_w, k_idx), u_summ_lambda_mu_for_p),
					fft::corr_real(apply_corr_factor(diffence_w, k_idx), u_summ_lambda_mu_for_p));
				const auto u_lambda_for_s = fft::summ_real(
					fft::conv_real(apply_conv_factor(u[t_idx][z_idx], k_idx), lambda),
					fft::corr_real(apply_corr_factor(u[t_idx][z_idx], k_idx), lambda));
				s[t_idx + 2u][z_idx][k_idx] = s[t_idx][z_idx][k_idx]
					+ g_t_grid_step * 0.5 * (w_lambda_for_s[k_idx] - u_lambda_for_s[k_idx])
					+ f_x_h[t_idx] * std::cos(k_idx * g_z_limit_value / 2.);
			}
		}
	}

	std::vector<std::vector<double>> x, y, z;
	for(auto x_idx = 2u; x_idx < g_z_grid_size - 2; x_idx += 2)
	{
		std::vector<double> x_row, y_row, z_row;
		for(auto z_idx = 2u; z_idx < g_z_grid_size - 2; z_idx += 2)
		{
			const auto t_idx = 4;
			x_row.push_back(x_idx * g_z_grid_step);
			y_row.push_back(z_idx * g_z_grid_step);
			const auto f_value = u_func(u, x_idx, z_idx, t_idx);
			z_row.push_back(f_value);
		}
		x.push_back(x_row);
		y.push_back(y_row);
		z.push_back(z_row);
	}

	namespace plt = matplotlibcpp;

	plt::plot_surface(x, y, z);
	plt::show();
}
catch(const utils::MKLException & exception)
{
	std::cout << "MKL exception happend: " << exception.what();
}
catch(const std::exception & exception)
{
	std::cout << "Another exception happend: " << exception.what();
}
