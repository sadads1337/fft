#pragma once

#include <type_traits>

namespace utils
{

template <typename ContainerType>
std::remove_reference_t<ContainerType> move_and_clear(ContainerType && val);

} // namespace utils

#include "MoveAndClear.inl"
