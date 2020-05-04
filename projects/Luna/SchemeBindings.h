#pragma once

#include <vector>

#include <ucenv/ucenv.h>

extern "C" {

//! import init_vector(int #fgnum, int #fgsize, name) as init_vector_x;
//! import init_vector(int #fgnum, int #fgsize, name) as init_vector_y;
void init_vector(
    std::int32_t fg_num,
    std::int32_t fg_size,
    luna::ucenv::OutputDF & vec);

//! import sum_vectors(value #vec1, value #vec2, int #fgsize, name #result) as sum;
void sum_vectors(
    const luna::ucenv::InputDF & x,
    const luna::ucenv::InputDF & y,
    std::int32_t fg_size,
    luna::ucenv::OutputDF & result);

//! import check_vector(value #vec, int #fgnum, int #size) as check_vector;
void check_vector(
    const luna::ucenv::InputDF & x,
    std::int32_t fg_num,
    std::int32_t fg_size);

} // extern "C"

