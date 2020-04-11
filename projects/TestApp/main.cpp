#include <matplotlibcpp/matplotlibcpp.h>

#include <Scheme/Scheme.h>

#include <Core/Types.h>
#include <Core/MKL/Utils.h>

#include <iostream>

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
	
	for(auto t_idx = 0u; t_idx < 2u/*g_t_grid_size*/; ++t_idx)
	{
		calculate_one_step(values_1, values_2, env, t_idx);

		//! new in values_2 now; we don't need values_2, swap here, do not copy
		std::swap(values_1, values_2);

		//! draw values_1;
		std::vector<double> x, y;
		for(auto x_idx = 0u; x_idx < scheme::g_z_grid_size; ++x_idx)
		{
			x.push_back(static_cast<Precision>(x_idx) * scheme::g_z_grid_step);
			constexpr auto z_0 = scheme::g_z_grid_size / 2;
			const auto f_value = scheme::u_func(values_1.u, x_idx, z_0);
			y.push_back(f_value);
		}
		namespace plt = matplotlibcpp;
		plt::plot(x, y);
		plt::show();

		/*std::vector<double> x, y;
		for(auto x_idx = 0u; x_idx < g_z_grid_size; ++x_idx)
		{
			x.push_back(static_cast<Precision>(x_idx) * g_z_grid_step);
			constexpr auto z_0 = g_z_grid_size / 2;
			const auto f_value = u_func(values_1.u, z_0, x_idx);
			y.push_back(f_value);
		}
		namespace plt = matplotlibcpp;
		plt::plot(x, y);
		plt::show();*/

		/*{
			std::vector<std::vector<double>> x, y, z;
			for (auto x_idx = 0u; x_idx < g_z_grid_size; ++x_idx) {
				std::vector<double> x_row, y_row, z_row;
				for (auto z_idx = 0u; z_idx < g_z_grid_size; ++z_idx) {
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
		}*/
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
