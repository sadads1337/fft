#pragma once

#include <type_traits>

namespace fft
{

template <typename ContainerType>
std::remove_reference_t<ContainerType> move_and_clear(ContainerType && val);

} // namespace fft

#include "MoveAndClear.inl"
