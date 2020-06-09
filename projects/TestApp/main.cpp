//#include <Core/MKL/Utils.h>
#include <Core/Types.h>
#include <Scheme/Scheme.h>
#include <matplotlibcpp/matplotlibcpp.h>

#include <iostream>
#include <thread>

int main() try {
  constexpr auto t_idx_limit = 2u;
  const auto params = scheme::create_params(std::thread::hardware_concurrency(), 5.115, 1u << 10u,
                                            t_idx_limit, 1. / static_cast<Precision>(1u << 11u), 32u);
  scheme::Values values_1{
      make_grid_2d(params.z_grid_size, params.k_limit),
      make_grid_2d(params.z_grid_size, params.k_limit),
      make_grid_2d(params.z_grid_size, params.k_limit),
      make_grid_2d(params.z_grid_size, params.k_limit),
      make_grid_2d(params.z_grid_size, params.k_limit),
  };
  scheme::Values values_2{
      make_grid_2d(params.z_grid_size, params.k_limit),
      make_grid_2d(params.z_grid_size, params.k_limit),
      make_grid_2d(params.z_grid_size, params.k_limit),
      make_grid_2d(params.z_grid_size, params.k_limit),
      make_grid_2d(params.z_grid_size, params.k_limit),
  };
  const scheme::Env env{
      Grid1D(params.k_limit, static_cast<Precision>(1.)),
      Grid1D(params.k_limit, static_cast<Precision>(1.)),
      Grid1D(params.k_limit, static_cast<Precision>(1.)),
      scheme::source(4, static_cast<Precision>(1.), params.t_grid_step,
                     params.z_grid_step, params.t_limit),
  };

  const auto draw_plot = [&](const auto& values) {
    //! draw values_1;
    std::vector<double> x, y;
    for (auto x_idx = 0u; x_idx < params.z_grid_size; ++x_idx) {
      x.push_back(static_cast<Precision>(x_idx) * params.z_grid_step);
      const auto z_0 = params.z_grid_size / 2u;
      const auto f_value = scheme::u_func(values.u, x_idx, z_0, params);
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
                            const auto f_value = u_func(values_1.u, x_idx,
    z_idx); z_row.push_back(f_value);
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
  };

  main_loop_for_t(values_1, values_2, env, draw_plot,
                  params);

} catch (const utils::MKLException& exception) {
  std::cout << "MKL exception happend: " << exception.what();
} catch (const std::exception& exception) {
  std::cout << "Another exception happend: " << exception.what();
}
