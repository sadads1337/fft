#include <Core/MKL/Utils.h>
#include <Core/TimeGuard.h>
#include <Core/Types.h>
#include <Scheme/Scheme.h>

#include <chrono>
#include <cstdlib>
#include <iostream>
#include <thread>

int main([[maybe_unused]] int argc, [[maybe_unused]] char** argv) try {
  scheme::Values values_1{};
  scheme::Values values_2{};
  const scheme::Env env{
      Grid1D(scheme::g_k_limit, static_cast<Precision>(1.)),
      Grid1D(scheme::g_k_limit, static_cast<Precision>(1.)),
      Grid1D(scheme::g_k_limit, static_cast<Precision>(1.)),
      scheme::source(4, static_cast<Precision>(1.), scheme::g_t_grid_step,
                     scheme::g_z_grid_step, scheme::g_t_grid_size),
  };

  const auto num_threads = argc <= 1 ? std::thread::hardware_concurrency()
                                     : static_cast<size_t>(std::atoi(argv[1]));

  main_loop_for_t(values_1, values_2, env, scheme::g_t_grid_size, {},
                  num_threads);
} catch (const utils::MKLException& exception) {
  std::cout << "MKL exception happend: " << exception.what();
} catch (const std::exception& exception) {
  std::cout << "Another exception happend: " << exception.what();
}
