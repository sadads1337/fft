#pragma once
#include "MoveAndClear.h"
#include "MoveOnly.h"

namespace utils {

template <typename ContainerType>
std::remove_reference_t<ContainerType> move_and_clear(ContainerType&& val) {
  static_assert(
      std::is_move_constructible<std::remove_reference_t<ContainerType>>::value,
      "Need move constructor");
  std::remove_reference_t<ContainerType> temp(move_only(val));
  val.clear();
  return temp;
}

}  // namespace utils
