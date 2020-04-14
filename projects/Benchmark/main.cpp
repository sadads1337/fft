#include <Scheme/Scheme.h>

#include <Core/TimeGuard.h>
#include <Core/Types.h>
#include <Core/MKL/Utils.h>

#include <iostream>
#include <chrono>

int main() try
{
	scheme::Values values_1{};
	scheme::Values values_2{};
	const scheme::Env env{
		Grid1D(scheme::g_k_limit, static_cast<Precision>(1.)),
		Grid1D(scheme::g_k_limit, static_cast<Precision>(1.)),
		Grid1D(scheme::g_k_limit, static_cast<Precision>(1.)),
		scheme::source(
			4,
			static_cast<Precision>(1.),
			scheme::g_t_grid_step,
			scheme::g_z_grid_step,
			scheme::g_t_grid_size),
	};

	{
		utils::TimeGuard time_guard{};
		for (auto t_idx = 0u; t_idx < scheme::g_t_grid_size; ++t_idx) {
			calculate_one_step(values_1, values_2, env, t_idx);

			//! new in values_2 now; we don't need values_2, swap here, do not copy
			std::swap(values_1, values_2);
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
