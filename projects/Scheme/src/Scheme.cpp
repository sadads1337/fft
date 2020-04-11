#include <Scheme/Scheme.h>


#include <Core/IdxUtils.h>
#include <Core/MKL/Utils.h>

#include "OpenMP.h"

namespace scheme
{

//! \todo: fixme скопипащенный legacy c
Grid1D source(int IG, float WN7, float DT, float DZ, int K8)
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

Precision u_func(const Grid2D & u, const size_t x_idx, const size_t z_idx)
{
	constexpr auto b = g_z_limit_value;
	static_assert(std::is_same_v<std::remove_cv_t<decltype(b)>, Precision>, "b must be the same type with main precision type");
	auto result = static_cast<Precision>(0.);
	for (auto k_idx = 0u; k_idx < g_k_limit; ++k_idx)
	{
		const auto x = static_cast<Precision>(x_idx) * g_z_grid_step;
		result += u[z_idx][k_idx] * std::sin(static_cast<Precision>(k_idx) * g_pi / b * x);
	}
	return result * (static_cast<Precision>(2.) / b);
}

void calculate_one_step(const Values & prev_values, Values & values, const Env & env, const size_t t_idx)
{
	const auto & prev_u = prev_values.u;
	const auto & prev_w = prev_values.w;
	const auto & prev_p = prev_values.p;
	const auto & prev_q = prev_values.q;
	const auto & prev_s = prev_values.s;

	auto & u = values.u;
	auto & w = values.w;
	auto & p = values.p;
	auto & q = values.q;
	auto & s = values.s;

	const auto & rho = env.rho;
	const auto & lambda = env.lambda;
	const auto & mu = env.mu;

	for (auto z_idx = 1u; z_idx < g_z_grid_size - 1u; ++z_idx)
	{
		//! For u
		const auto der_q = utils::apply_operation(
			prev_q[z_idx + 1u],
			prev_q[z_idx],
			[](const auto & lhs, const auto & rhs) { return (lhs - rhs) / g_z_grid_step; });
		const auto q_rho_for_u = fft::summ_real(
			fft::conv_real(utils::apply_conv_factor(der_q, g_pi / g_z_limit_value), rho),
			fft::corr_real(utils::apply_corr_factor(der_q, g_pi / g_z_limit_value), rho));
		const auto p_rho_for_u = fft::summ_real(
			fft::conv_real(utils::apply_conv_factor(prev_p[z_idx], g_pi / g_z_limit_value), rho),
			fft::corr_real(utils::apply_corr_factor(prev_p[z_idx], g_pi / g_z_limit_value), rho));

		//! For w
		const auto q_rho_for_w = fft::summ_real(
			fft::conv_real(utils::apply_conv_factor(prev_q[z_idx + 1u], g_pi / g_z_limit_value), rho),
			fft::corr_real(utils::apply_corr_factor(prev_q[z_idx + 1u], g_pi / g_z_limit_value), rho));
		const auto der_s = utils::apply_operation(
			prev_s[z_idx + 1u],
			prev_s[z_idx],
			[](const auto & lhs, const auto & rhs) { return (lhs - rhs) / g_z_grid_step; });
		const auto s_rho_for_w = fft::summ_real(
			fft::conv_real(utils::apply_conv_factor(der_s, g_pi / g_z_limit_value), rho),
			fft::corr_real(utils::apply_corr_factor(der_s, g_pi / g_z_limit_value), rho));

		//! For p
		const auto der_w = utils::apply_operation(
			prev_w[z_idx + 1u],
			prev_w[z_idx],
			[](const auto & lhs, const auto & rhs) { return (lhs - rhs) / g_z_grid_step; });
		const auto w_lambda_for_p = fft::summ_real(
			fft::conv_real(utils::apply_conv_factor(der_w, g_pi / g_z_limit_value), lambda),
			fft::corr_real(utils::apply_corr_factor(der_w, g_pi / g_z_limit_value), lambda));
		const auto summ_lambda_mu = utils::apply_operation(
			lambda,
			mu,
			[](const auto & lhs, const auto & rhs) { return lhs + static_cast<Precision>(2.) * rhs; });
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
			fft::conv_real(utils::apply_conv_factor(prev_w[z_idx + 1u], g_pi / g_z_limit_value), mu),
			fft::corr_real(utils::apply_corr_factor(prev_w[z_idx + 1u], g_pi / g_z_limit_value), mu));

		//! For s
		const auto w_summ_lambda_mu_for_s = fft::summ_real(
			fft::conv_real(utils::apply_conv_factor(der_w, g_pi / g_z_limit_value), summ_lambda_mu),
			fft::corr_real(utils::apply_corr_factor(der_w, g_pi / g_z_limit_value), summ_lambda_mu));
		const auto u_lambda_for_s = fft::summ_real(
			fft::conv_real(utils::apply_conv_factor(prev_u[z_idx], g_pi / g_z_limit_value), lambda),
			fft::corr_real(utils::apply_corr_factor(prev_u[z_idx], g_pi / g_z_limit_value), lambda));

		//! For f
		constexpr auto z_0_idx = g_z_grid_size / 2;
		const auto f_x_h = z_idx == z_0_idx
			? source(4, static_cast<Precision>(1.), g_t_grid_step, g_z_grid_step, g_t_grid_size)
			: Grid1D(g_t_grid_size, static_cast<Precision>(0.));

		for (auto k_idx = 0u; k_idx < g_k_limit; ++k_idx)
		{
			constexpr auto half = static_cast<Precision>(0.5);
			static_assert(std::is_same_v<std::remove_cv_t<decltype(half)>, double>, "g_z_grid_size must be double");

			const auto z_0 = static_cast<Precision>(z_0_idx) * g_z_grid_step;

			u[z_idx][k_idx] = prev_u[z_idx][k_idx]
				+ half * g_t_grid_step * (q_rho_for_u[k_idx] - p_rho_for_u[k_idx]);
			w[z_idx + 1u][k_idx] = prev_w[z_idx + 1u][k_idx]
				+ half * g_t_grid_step * (q_rho_for_w[k_idx] + s_rho_for_w[k_idx]);
			p[z_idx][k_idx] = prev_p[z_idx][k_idx]
				+ half * g_t_grid_step * (w_lambda_for_p[k_idx] - u_summ_lambda_mu_for_p[k_idx])
				+ f_x_h[t_idx] * std::cos(static_cast<Precision>(k_idx) * g_pi / g_z_limit_value * static_cast<Precision>(z_0));
			q[z_idx + 1u][k_idx] = prev_q[z_idx + 1u][k_idx]
				+ half * g_t_grid_step * (u_mu_for_q[k_idx] - w_mu_for_q[k_idx]);
			s[z_idx][k_idx] = prev_s[z_idx][k_idx]
				+ half * g_t_grid_step * (w_summ_lambda_mu_for_s[k_idx] - u_lambda_for_s[k_idx])
				+ f_x_h[t_idx] * std::cos(static_cast<Precision>(k_idx) * g_pi / g_z_limit_value  * static_cast<Precision>(z_0));
		}
	}
}

} // namespace scheme
