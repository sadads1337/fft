#pragma once

#include <Core/GridUtils.h>
#include <Core/Types.h>

#include <cmath>
#include <functional>
#include <type_traits>

namespace scheme {

struct Params final {
  size_t num_threads;

  Precision z_limit_value;
  size_t z_grid_size;
  Precision z_grid_step;

  size_t t_limit;
  Precision t_grid_step;

  size_t k_limit;
};

Params create_params(size_t num_threads, Precision z_limit_value, size_t z_grid_size,
                     size_t t_limit, Precision t_grid_step, size_t k_limit);

inline static const auto g_pi =
    std::atan(static_cast<Precision>(1.)) * static_cast<Precision>(4.);

struct Values final {
  Grid2D u;
  Grid2D w;
  Grid2D p;
  Grid2D q;
  Grid2D s;
};

struct Env final {
  Grid1D rho;
  Grid1D lambda;
  Grid1D mu;
  Grid1D f;

  //! \todo: maybe constant initilized, but there is a bug in LuNA df
  //! initialization code.
};

//! legacy c
Grid1D source(int IG, Precision WN7, Precision DT, Precision DZ, int K8);

Precision u_func(const Grid2D& u, size_t x_idx, size_t z_idx, const Params& params);

void calculate_one_step(const Values& prev_values, Values& values,
                        const Env& env, size_t t_idx, size_t fg_num,
                        size_t fg_size, size_t fg_count,
                        const Params& params);

void main_loop_for_t(Values& prev_values, Values& values, const Env& env,
                     const std::function<void(const Values&)>& callback,
                     const Params& params);

}  // namespace scheme
