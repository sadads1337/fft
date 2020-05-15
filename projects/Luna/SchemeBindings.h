#pragma once

#include <cstdint>

#include <ucenv/ucenv.h>

#include <Scheme/Scheme.h>

extern "C" {

//! import init_env(int #fgnum, int #fgsize, name) as init_env;
//! \note scheme::g_z_grid_size = fg_num * fg_size;
void init_env(
    std::int32_t fg_num,
    std::int32_t fg_size,
    luna::ucenv::OutputDF & env);

//! import check_env(value #env, int #fgnum, int #size) as check_env;
//! \note scheme::g_z_grid_size = fg_num * fg_size;
//! \warning only for debug purposes;
void check_env(
    const luna::ucenv::InputDF & env,
    std::int32_t fg_num,
    std::int32_t fg_size);

//! import init_values(int #fgnum, int #fgsize, name) as init_values;
//! \note scheme::g_z_grid_size = fg_num * fg_size;
void init_values(
    std::int32_t fg_num,
    std::int32_t fg_size,
    luna::ucenv::OutputDF & values);

//! import check_values(value #values, int #fgnum, int #size) as check_values;
//! \note scheme::g_z_grid_size = fg_num * fg_size;
//! \warning only for debug purposes;
void check_values(
    const luna::ucenv::InputDF & values,
    std::int32_t fg_num,
    std::int32_t fg_size);

//! import calculate_one_step(value #prev_values, value #env, int #t_idx, name #values) as calculate_one_step;
void calculate_one_step(
    const luna::ucenv::InputDF & prev_values,
    const luna::ucenv::InputDF & env,
    std::int32_t t_idx,
    luna::ucenv::OutputDF & values);

//! import check_one_step(
//!     value #prev_values,
//!     value #values,
//!     value #env,
//!     int #t_idx,
//!     int #fgnum,
//!     int #fgsize) as check_one_step;
//! \warning only for debug purposes;
void check_one_step(
    const luna::ucenv::InputDF & prev_values,
    const luna::ucenv::InputDF & values,
    const luna::ucenv::InputDF & env,
    std::int32_t t_idx,
    std::int32_t fg_num,
    std::int32_t fg_size);

} // extern "C"

