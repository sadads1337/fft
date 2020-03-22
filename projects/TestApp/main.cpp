#include <Core/MKL/Utils.h>
#include <Core/MakeWithCapacity.h>
#include <Core/MoveAndClear.h>
#include <Core/MoveOnly.h>
#include <Core/matplotlib-cpp/matplotlibcpp.h>
#include <Core/Types.h>
#include <Core/IdxUtils.h>

#include <iostream>
#include <array>
#include <type_traits>
#include <functional>

//! \todo: Студия не дружит с cassert, но очень нужно обмазать код assert'ми
#include <cassert>

auto make_grid_2d(const size_t z, const size_t k, const Precision default_value = 0.)
{
	assert(z && k);
	return Grid2D(z, Grid1D(k, default_value));
}

auto make_grid_3d(const size_t t, const size_t z, const size_t k, const Precision default_value = 0.)
{
	assert(t && z && k);
	return Grid3D(t, Grid2D(z, Grid1D(k, default_value)));
}

//! \todo: Сейчас сетка от 0. до g_*_limit_value, поправить для произвольной. И ниже.
constexpr auto g_z_limit_value = 1.;
constexpr auto g_t_limit_value = 1.;

//! Схема не является безусловно устойчивой -> ограничения должны удовлетворять условиям Куранта.
//! \todo: fixme!!! constexpr check with static assertion
constexpr auto g_t_grid_size = 101u;
constexpr auto g_z_grid_size = 51u;

constexpr auto g_k_limit = 10u;

static_assert(std::is_same_v<std::remove_cv_t<decltype(g_z_limit_value)>, double>, "g_z_grid_size must be double");
static_assert(std::is_same_v<std::remove_cv_t<decltype(g_t_limit_value)>, double>, "g_t_limit_value must be double");
//! Care double value might be out of range.
constexpr auto g_z_grid_step = g_z_limit_value / static_cast<Precision>(g_z_grid_size);
constexpr auto g_t_grid_step = g_t_limit_value / static_cast<Precision>(g_t_grid_size);

inline const auto g_pi = std::atan(1.) * 4.;

//! legacy c
auto source(int IG, float WN7, float DT, float DZ, int K8)
{
	auto F = Grid1D(static_cast<size_t>(K8));
    float C,B,T,T1,T2,T3,T4,PII,X,A,S;
    int J,N,GI;
	 PII=acos(-1.0);
	 int NPT1=K8;
	 GI=IG;
      C=DT/DZ;
       B=0.00;
       T1=1.500;
       T2=3.500;
       T3=5.500;
       T=2.00*PII*WN7;
    
       if(GI<=2.0100) T4=T1;
    else {
     if((GI<=4.0100)&&(GI>=2.0100)) T4=T2;
        else T4=T3;
    }
       N=int(T4/(DT*WN7));
       X=T4/(2.00*WN7);
       A=(T/GI)*(T/GI);
       S=-X;
    for(J=0; J<NPT1; J++){
     if(J+1<=N+1) F[J]=C*exp((-A)*(S*S))*cos(T*S+B);
     else F[J]=0.00;
     S=S+DT;
    }
   for(int i=0; i<NPT1; i++) F[i]=F[i]*DT/(DZ);
	
	return F;
}

auto u_func(const Grid2D & u, const size_t x_idx, const size_t z_idx)
{
	constexpr auto b = g_z_limit_value;
	static_assert(std::is_same_v<std::remove_cv_t<decltype(b)>, Precision>, "b must be the same type with main precision type");
	auto result = 0.;
	for (auto k_idx = 1u; k_idx < g_k_limit; ++k_idx)
	{
		const auto x = static_cast<Precision>(x_idx) * g_z_grid_step;

		result += u[z_idx][k_idx]
			* std::sin(static_cast<Precision>(k_idx) * g_pi / b * x);
	}
	return result * (2. / b);
}

/*auto w_func(const Grid3D & w, const Precision x, const Precision z, const Precision t)
{
	const auto b = 1.;
	static_assert(std::is_same_v<std::remove_cv_t<decltype(b)>, Precision>, "b must be the same type with main precision type");
	auto result = 1. / b * w[get_idx(0., g_t_limit_value, t, g_t_grid_step)][get_idx(0., g_z_limit_value, z, g_z_grid_step)][0];
	auto sum = 0.;
	for (auto k_idx = 2; k_idx < g_k_limit; k_idx+=2)
	{
		sum += w[get_idx(0., g_t_limit_value, t, g_t_grid_step)][get_idx(0., g_z_limit_value, z, g_z_grid_step)][k_idx]
			* std::cos(static_cast<float>(k_idx) * x);
	}
	return result + 2. / b * sum;
}*/

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
	const Grid1D rho = Grid1D(g_z_grid_size, 1.);
	const Grid1D lambda = Grid1D(g_z_grid_size, 1.);
	const Grid1D mu = Grid1D(g_z_grid_size, 1.);
};

void calculate_one_step(Values & prev_values, Values & values, const Env & env, const size_t t_idx)
{
	auto & prev_u = prev_values.u;
	auto & prev_w = prev_values.w;
	auto & prev_p = prev_values.p;
	auto & prev_q = prev_values.q;
	auto & prev_s = prev_values.s;

	auto & u = values.u;
	auto & w = values.w;
	auto & p = values.p;
	auto & q = values.q;
	auto & s = values.s;
	
	auto & rho = env.rho;
	auto & lambda = env.lambda;
	auto & mu = env.mu;

	for (auto z_idx = 1u; z_idx < g_z_grid_size - 1u; ++z_idx)
	{
		//! For u
		const auto der_q = utils::apply_operation(
			prev_q[z_idx],
			prev_q[z_idx - 1u],
			[](const auto & lhs, const auto & rhs) { return (lhs - rhs) / g_z_grid_step; });
		const auto q_rho_for_u = fft::summ_real(
			fft::conv_real(utils::apply_conv_factor(der_q, g_pi / g_z_limit_value), rho),
			fft::corr_real(utils::apply_corr_factor(der_q, g_pi / g_z_limit_value), rho));
		const auto p_rho_for_u = fft::summ_real(
			fft::conv_real(utils::apply_conv_factor(prev_p[z_idx], g_pi / g_z_limit_value), rho),
			fft::corr_real(utils::apply_corr_factor(prev_p[z_idx], g_pi / g_z_limit_value), rho));

		//! For w
		const auto q_rho_for_w = fft::summ_real(
			fft::conv_real(utils::apply_conv_factor(prev_q[z_idx], g_pi / g_z_limit_value), rho),
			fft::corr_real(utils::apply_corr_factor(prev_q[z_idx], g_pi / g_z_limit_value), rho));
		const auto der_s = utils::apply_operation(
			prev_s[z_idx + 1u],
			prev_s[z_idx],
			[](const auto & lhs, const auto & rhs) { return (lhs - rhs) / g_z_grid_step; });
		const auto s_rho_for_w = fft::summ_real(
			fft::conv_real(utils::apply_conv_factor(der_s, g_pi / g_z_limit_value), rho),
			fft::corr_real(utils::apply_corr_factor(der_s, g_pi / g_z_limit_value), rho));

		//! For p
		const auto der_w = utils::apply_operation(
			prev_w[z_idx],
			prev_w[z_idx - 1u],
			[](const auto & lhs, const auto & rhs) { return (lhs - rhs) / g_z_grid_step; });
		const auto w_lambda_for_p = fft::summ_real(
			fft::conv_real(utils::apply_conv_factor(der_w, g_pi / g_z_limit_value), lambda),
			fft::corr_real(utils::apply_corr_factor(der_w, g_pi / g_z_limit_value), lambda));
		const auto summ_lambda_mu = utils::apply_operation(
			lambda,
			mu,
			[](const auto & lhs, const auto & rhs) { return lhs + 2. * rhs; });
		const auto u_summ_lambda_mu_for_p = fft::summ_real(
			fft::conv_real(utils::apply_conv_factor(prev_u[z_idx], g_pi / g_z_limit_value), summ_lambda_mu),
			fft::corr_real(utils::apply_corr_factor(prev_u[z_idx], g_pi / g_z_limit_value), summ_lambda_mu));

		//! For q
		const auto der_u = utils::apply_operation(
			prev_u[z_idx + 1u],
			prev_u[z_idx],
			[](const auto & lhs, const auto & rhs) { return (lhs - rhs) / g_z_grid_step; });
		const auto u_mu_for_q = fft::summ_real(
			fft::conv_real(utils::apply_conv_factor(der_u, g_pi / g_z_limit_value), mu),
			fft::corr_real(utils::apply_corr_factor(der_u, g_pi / g_z_limit_value), mu));
		const auto w_mu_for_q = fft::summ_real(
			fft::conv_real(utils::apply_conv_factor(prev_w[z_idx], g_pi / g_z_limit_value), mu),
			fft::corr_real(utils::apply_corr_factor(prev_w[z_idx], g_pi / g_z_limit_value), mu));

		//! For s
		const auto w_summ_lambda_mu_for_s = fft::summ_real(
			fft::conv_real(utils::apply_conv_factor(der_w, g_pi / g_z_limit_value), summ_lambda_mu),
			fft::corr_real(utils::apply_corr_factor(der_w, g_pi / g_z_limit_value), summ_lambda_mu));
		const auto u_lambda_for_s = fft::summ_real(
			fft::conv_real(utils::apply_conv_factor(prev_u[z_idx], g_pi / g_z_limit_value), lambda),
			fft::corr_real(utils::apply_corr_factor(prev_u[z_idx], g_pi / g_z_limit_value), lambda));

		//! For f
		constexpr auto z_0_idx = (g_z_grid_size / 2);
		const auto f_x_h = z_idx == z_0_idx
			? source(4, 1., g_t_grid_step, g_z_grid_step, g_t_grid_size)
			: Grid1D(g_t_grid_size, 0.);

		for (auto k_idx = 0u; k_idx < g_k_limit; ++k_idx)
		{
			constexpr auto half = 0.5;
			static_assert(std::is_same_v<std::remove_cv_t<decltype(half)>, double>, "g_z_grid_size must be double");

			const auto z_0 = static_cast<Precision>(z_0_idx) * g_z_grid_step;
			u[z_idx][k_idx] = static_cast<Precision>(z_idx == z_0_idx) * std::cos(static_cast<Precision>(k_idx * g_pi / g_z_limit_value  * z_0));

			/*u[z_idx][k_idx] = prev_u[z_idx][k_idx]
				+ half * g_t_grid_step * (q_rho_for_u[k_idx] - p_rho_for_u[k_idx]);
			w[z_idx][k_idx] = prev_w[z_idx - 1u][k_idx]
				+ half * g_t_grid_step * (q_rho_for_w[k_idx] + s_rho_for_w[k_idx]);
			p[z_idx][k_idx] = prev_p[z_idx][k_idx]
				+ half * g_t_grid_step * (w_lambda_for_p[k_idx] - u_summ_lambda_mu_for_p[k_idx])
				+ f_x_h[t_idx] * std::cos(static_cast<Precision>(k_idx) * g_pi / g_z_limit_value * static_cast<Precision>(z_0) * g_z_grid_step);
			q[z_idx][k_idx] = prev_q[z_idx - 1u][k_idx]
				+ half * g_t_grid_step * (u_mu_for_q[k_idx] - w_mu_for_q[k_idx]);
			s[z_idx][k_idx] = prev_s[z_idx][k_idx]
				+ half * g_t_grid_step * (w_summ_lambda_mu_for_s[k_idx] - u_lambda_for_s[k_idx])
				+ f_x_h[t_idx] * std::cos(static_cast<Precision>(k_idx * g_pi / g_z_limit_value  * static_cast<Precision>(z_0) * g_z_grid_step));*/
		}
	}
}

int main() try
{
	Values values_1{};
	Values values_2{};
	const Env env{};
	
	for(auto t_idx = 0u; t_idx < 2u/*g_t_grid_size*/; ++t_idx)
	{
		calculate_one_step(values_1, values_2, env, t_idx);

		//! new in values_2 now; we don't need values_2, swap here, do not copy
		std::swap(values_1, values_2);

		//! draw values_1;
		std::vector<std::vector<double>> x, y, z;
		for(auto x_idx = 0u; x_idx < g_z_grid_size; ++x_idx)
		{
			std::vector<double> x_row, y_row, z_row;
			for(auto z_idx = 0u; z_idx < g_z_grid_size; ++z_idx)
			{
				x_row.push_back(x_idx * g_z_grid_step);
				y_row.push_back(z_idx * g_z_grid_step);
				const auto f_value = u_func(values_1.u, x_idx, z_idx);
				z_row.push_back(f_value);
			}
			x.push_back(x_row);
			y.push_back(y_row);
			z.push_back(z_row);
		}

		namespace plt = matplotlibcpp;
		
		plt::plot_surface(x, y, z);
		plt::xlabel("x");
		plt::ylabel("z");
		plt::show();
	}
	
	/*auto u = make_grid(g_t_grid_size, g_z_grid_size, g_k_limit);
	auto w = make_grid(g_t_grid_size, g_z_grid_size, g_k_limit);
	auto p = make_grid(g_t_grid_size, g_z_grid_size, g_k_limit);
	auto q = make_grid(g_t_grid_size, g_z_grid_size, g_k_limit);
	auto s = make_grid(g_t_grid_size, g_z_grid_size, g_k_limit);

	const auto rho = Grid1D(g_z_grid_size, 1.);
	const auto lambda = Grid1D(g_z_grid_size, 1.);
	const auto mu = Grid1D(g_z_grid_size, 1.);

	for(auto t_idx = 0u; t_idx < g_t_grid_size - 1; ++t_idx)
	{
		//! Edge equations
		const auto s_edge = utils::apply_operation(
			s[t_idx][2],
			s[t_idx][2],
			[](const auto & lhs, const auto & rhs) { return (lhs + rhs) / g_z_grid_step; });
		const auto s_edge_rho_for_w = fft::summ_real(
			fft::conv_real(utils::apply_conv_factor(s_edge), rho),
			fft::corr_real(utils::apply_corr_factor(s_edge), rho));

		for (auto k_idx = 0u; k_idx < g_k_limit; ++k_idx)
		{
			w[t_idx + 1u][1u][k_idx] = g_t_grid_step * 0.5 * s_edge_rho_for_w[k_idx]
				+ w[t_idx][1u][k_idx];
		}

		for (auto z_idx = 2u; z_idx < g_z_grid_size - 2u; z_idx += 2)
		{
			//! For u
			const auto diffence_q = utils::apply_operation(
				q[t_idx][z_idx + 1],
				q[t_idx][z_idx - 1],
				[](const auto & lhs, const auto & rhs) { return (lhs - rhs) / g_z_grid_step; });
			const auto q_rho_for_u = fft::summ_real(
				fft::conv_real(utils::apply_conv_factor(diffence_q), rho),
				fft::corr_real(utils::apply_corr_factor(diffence_q), rho));
			const auto p_rho_for_u = fft::summ_real(
				fft::conv_real(utils::apply_conv_factor(p[t_idx][z_idx]), rho),
				fft::corr_real(utils::apply_corr_factor(p[t_idx][z_idx]), rho));

			//! For w
			const auto q_rho_for_w = fft::summ_real(
				fft::conv_real(utils::apply_conv_factor(q[t_idx][z_idx + 1]), rho),
				fft::corr_real(utils::apply_corr_factor(q[t_idx][z_idx + 1]), rho));
			const auto diffence_s = utils::apply_operation(
				s[t_idx][z_idx + 2],
				s[t_idx][z_idx],
				[](const auto & lhs, const auto & rhs) { return (lhs - rhs) / g_z_grid_step; });
			const auto s_rho_for_w = fft::summ_real(
				fft::conv_real(utils::apply_conv_factor(diffence_s), rho),
				fft::corr_real(utils::apply_corr_factor(diffence_s), rho));

			//! For p
			const auto diffence_w = utils::apply_operation(
				s[t_idx][z_idx + 1],
				s[t_idx][z_idx - 1],
				[](const auto & lhs, const auto & rhs) { return (lhs - rhs) / g_z_grid_step; });
			const auto w_lambda_for_p = fft::summ_real(
				fft::conv_real(utils::apply_conv_factor(diffence_w), lambda),
				fft::corr_real(utils::apply_corr_factor(diffence_w), lambda));
			const auto summ_lambda_mu = utils::apply_operation(
				lambda,
				mu,
				[](const auto & lhs, const auto & rhs) { return lhs + 2. * rhs; });
			const auto u_summ_lambda_mu_for_p = fft::summ_real(
				fft::conv_real(utils::apply_conv_factor(u[t_idx][z_idx]), summ_lambda_mu),
				fft::corr_real(utils::apply_corr_factor(u[t_idx][z_idx]), summ_lambda_mu));

			//! For q
			const auto diffence_u = utils::apply_operation(
				u[t_idx][z_idx + 2],
				u[t_idx][z_idx],
				[](const auto & lhs, const auto & rhs) { return (lhs - rhs) / g_z_grid_step; });
			const auto u_mu_for_q = fft::summ_real(
				fft::conv_real(utils::apply_conv_factor(diffence_u), mu),
				fft::corr_real(utils::apply_corr_factor(diffence_u), mu));
			const auto w_mu_for_q = fft::summ_real(
				fft::conv_real(utils::apply_conv_factor(w[t_idx][z_idx + 1]), mu),
				fft::corr_real(utils::apply_corr_factor(w[t_idx][z_idx + 1]), mu));

			//! For s
			const auto w_lambda_for_s = fft::summ_real(
				fft::conv_real(utils::apply_conv_factor(diffence_w), u_summ_lambda_mu_for_p),
				fft::corr_real(utils::apply_corr_factor(diffence_w), u_summ_lambda_mu_for_p));
			const auto u_lambda_for_s = fft::summ_real(
				fft::conv_real(utils::apply_conv_factor(u[t_idx][z_idx]), lambda),
				fft::corr_real(utils::apply_corr_factor(u[t_idx][z_idx]), lambda));

			//! f
			//const auto f_x_h = f(4, 10., g_t_grid_step, g_z_grid_step, g_k_limit);
			constexpr auto z_0 = (g_z_grid_size / 2);
			const auto f_x_h = (z_idx == z_0)
				? source(4, 10., g_t_grid_step, g_z_grid_step, g_t_grid_size)
				: Grid1D(g_t_grid_size, 0.);

			for (auto k_idx = 0u; k_idx < g_k_limit; ++k_idx)
			{
				u[t_idx + 1u][z_idx][k_idx] = g_t_grid_step * 0.5 * (q_rho_for_u[k_idx] - p_rho_for_u[k_idx])
					+ u[t_idx][z_idx][k_idx];

				w[t_idx + 1u][z_idx + 1u][k_idx] = g_t_grid_step * 0.5 * (q_rho_for_w[k_idx] + s_rho_for_w[k_idx])
					+ w[t_idx][z_idx + 1u][k_idx];

				constexpr auto x_0 = g_z_limit_value / 2.;
				p[t_idx + 1u][z_idx][k_idx] = g_t_grid_step * 0.5 * (w_lambda_for_p[k_idx] - u_summ_lambda_mu_for_p[k_idx])
					+ p[t_idx][z_idx][k_idx]
					+ f_x_h[t_idx] * std::cos(k_idx * x_0);

				q[t_idx + 1u][z_idx + 1u][k_idx] = g_t_grid_step * 0.5 * (u_mu_for_q[k_idx] - w_mu_for_q[k_idx])
					+ q[t_idx][z_idx + 1u][k_idx];

				constexpr auto z_0 = x_0;
				s[t_idx + 1u][z_idx][k_idx] = g_t_grid_step * 0.5 * (w_lambda_for_s[k_idx] - u_lambda_for_s[k_idx])
					+ s[t_idx][z_idx][k_idx]
					+ f_x_h[t_idx] * std::cos(k_idx * z_0);
			}
		}
	}

	std::vector<std::vector<double>> x, y, z;
	for(auto x_idx = 2u; x_idx < g_z_grid_size - 2; x_idx += 2)
	{
		std::vector<double> x_row, y_row, z_row;
		for(auto z_idx = 2u; z_idx < g_z_grid_size - 2; z_idx += 2)
		{
			const auto t_idx = 2;
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
	plt::show();*/
}
catch(const utils::MKLException & exception)
{
	std::cout << "MKL exception happend: " << exception.what();
}
catch(const std::exception & exception)
{
	std::cout << "Another exception happend: " << exception.what();
}
