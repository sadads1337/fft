#include <Luna/SchemeBindings.h>

#include <cstdint>
#include <cassert>
#include <iostream>

void init_vector(
    const std::int32_t fg_num,
    const std::int32_t fg_size,
    luna::ucenv::OutputDF & vec)
{
  vec.create(fg_size * sizeof(std::int32_t));
  auto * const vec_array = vec.getData<std::int32_t>();
  for (auto idx = 0; idx < fg_size; ++idx)
  {
    vec_array[idx] = fg_num;
  }
}

void sum_vectors(
    const luna::ucenv::InputDF & x,
    const luna::ucenv::InputDF & y,
    const std::int32_t fg_size,
    luna::ucenv::OutputDF & result)
{
  assert(x.getArraySize<std::int32_t>() == y.getArraySize<std::int32_t>());
  assert(static_cast<std::size_t>(fg_size) == x.getArraySize<std::int32_t>());

  result.setArray(fg_size, std::int32_t{ 0 });
  result.create(fg_size * sizeof(std::int32_t));
  auto * const result_array = result.getData<std::int32_t>();
  const auto * const x_array = x.getData<std::int32_t>();
  const auto * const y_array = y.getData<std::int32_t>();

  for (auto idx = 0; idx < fg_size; ++idx)
  {
    result_array[idx] = x_array[idx] + y_array[idx];
  }

}

void check_vector(
    const luna::ucenv::InputDF & x,
    const std::int32_t fg_num,
    const std::int32_t fg_size)
{
  assert(static_cast<std::size_t>(fg_size) == x.getArraySize<std::int32_t>());

  const auto * const x_array = x.getData<std::int32_t>();
  const auto correct_value = static_cast<std::int32_t>(fg_num * 2);
  for (auto idx = 0; idx < fg_size; ++idx)
  {
    if (correct_value != x_array[idx])
    {
      std::cout << "failed in idx: " << idx;
    }
  }
}
