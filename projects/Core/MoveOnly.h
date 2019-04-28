#pragma once

#include <type_traits>

namespace fft
{

template <typename Type>
typename std::remove_reference<Type>::type && move_only(Type && val);

template <class InputIt, class OutputIt>
OutputIt move_only(InputIt first, InputIt last, OutputIt d_first);

} // namespace fft

#include "MoveOnly.inl"
