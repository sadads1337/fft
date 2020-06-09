#pragma once

#include <Scheme/Scheme.h>
#include <ucenv/ucenv.h>

#include <cstdint>

extern "C" {

//! import init_params(
//!     real #z_limit_value,
//!     int #z_grid_size,
//!     int #t_limit,
//!     real #t_grid_step,
//!     int #k_limit,
//!     name #params) as init_params;
//! \note scheme::g_z_grid_size = fg_num * fg_size;
void init_params(Precision z_limit_value, std::int32_t z_grid_size,
                 std::int32_t t_limit, Precision t_grid_step, std::int32_t k_limit,
                 luna::ucenv::OutputDF& params);

//! import init_env(value #params, int #fgnum, int #fgsize, name) as init_env;
//! \note scheme::g_z_grid_size = fg_num * fg_size;
void init_env(const luna::ucenv::InputDF& params, std::int32_t fg_num,
              std::int32_t fg_size, luna::ucenv::OutputDF& env);

//! import check_env(value #params, value #env, int #fgnum, int #size) as check_env;
//! \note scheme::g_z_grid_size = fg_num * fg_size;
//! \warning only for debug purposes;
void check_env(const luna::ucenv::InputDF& params, const luna::ucenv::InputDF& env,
               std::int32_t fg_num, std::int32_t fg_size);

//! import init_values(value #params, int #fgnum, int #fgsize, name) as init_values;
//! \note scheme::g_z_grid_size = fg_num * fg_size;
void init_values(const luna::ucenv::InputDF& params, std::int32_t fg_num,
                 std::int32_t fg_size, luna::ucenv::OutputDF& values);

//! import check_values(value #values, int #fgnum, int #size) as check_values;
//! \note scheme::g_z_grid_size = fg_num * fg_size;
//! \warning only for debug purposes;
void check_values(const luna::ucenv::InputDF& values, std::int32_t fg_num,
                  std::int32_t fg_size);

//! import calculate_one_step(
//!     value #params,
//!     value #prev_values,
//!     value #env,
//!     int #t_idx,
//!     int #fgnum,
//!     int #fgsize,
//!     int #fgcount,
//!     name #values) as calculate_one_step;
void calculate_one_step(const luna::ucenv::InputDF& params,
                        const luna::ucenv::InputDF& prev_values,
                        const luna::ucenv::InputDF& env, std::int32_t t_idx,
                        std::int32_t fg_num, std::int32_t fg_size,
                        std::int32_t fg_count, luna::ucenv::OutputDF& values);

//! import check_one_step(
//!     value #prev_values,
//!     value #values,
//!     value #env,
//!     int #t_idx,
//!     int #fgnum,
//!     int #fgsize
//!     int #fgcount) as check_one_step;
//! \warning only for debug purposes;
void check_one_step(const luna::ucenv::InputDF& prev_values,
                    const luna::ucenv::InputDF& values,
                    const luna::ucenv::InputDF& env, std::int32_t t_idx,
                    std::int32_t fg_num, std::int32_t fg_size,
                    std::int32_t fg_count);

}  // extern "C"
