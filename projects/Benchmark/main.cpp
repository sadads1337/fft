#include <Core/MKL/Utils.h>
#include <Core/TimeGuard.h>
#include <Core/Types.h>
#include <Scheme/Scheme.h>

#include <chrono>
#include <cstdlib>
#include <iostream>
#include <thread>

int main([[maybe_unused]] int argc, [[maybe_unused]] char** argv) try {
  const auto num_threads = std::atoi(argv[1]);
  const auto z_limit_value = std::atof(argv[2]);
  const auto z_grid_size = static_cast<size_t>(std::atoi(argv[3]));
  const auto t_limit = static_cast<size_t>(std::atoi(argv[4]));
  const auto t_grid_step = std::atof(argv[5]);
  const auto k_limit = static_cast<size_t>(std::atoi(argv[6]));

  const auto params = scheme::create_params(num_threads, z_limit_value, z_grid_size,
                                            t_limit, t_grid_step, k_limit);

  std::cout << "num_threads: " << num_threads
            << " z_limit: " << params.z_limit_value
            << " z_grid_size: " << params.z_grid_size
            << " t_limit: " << params.t_limit
            << " t_grid_step: " << params.t_grid_step
            << " k_limit: " << params.k_limit << std::endl;

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
      scheme::source(4, static_cast<Precision>(5.), params.t_grid_step,
                     params.z_grid_step, params.t_limit),
  };

  main_loop_for_t(values_1, values_2, env, {}, params);
} catch (const utils::MKLException& exception) {
  std::cout << "MKL exception happend: " << exception.what();
} catch (const std::exception& exception) {
  std::cout << "Another exception happend: " << exception.what();
}
