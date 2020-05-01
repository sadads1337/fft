#pragma once
#include <ciso646>
#include <memory>

#include "MoveOnly.h"

namespace utils {

template <typename Type>
typename std::remove_reference<Type>::type&& move_only(Type&& val) {
  static_assert(
      not std::is_const<typename std::remove_reference<Type>::type>::value,
      "value is const");
  return std::move(val);
}

template <class InputIt, class OutputIt>
OutputIt move_only(const InputIt first, const InputIt last,
                   const OutputIt d_first) {
  static_assert(
      not std::is_const<
          typename std::remove_reference<decltype(*first)>::type>::value,
      "value is const");
  return std::move(first, last, d_first);
}

}  // namespace utils
