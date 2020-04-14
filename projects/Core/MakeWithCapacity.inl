#pragma once

#include "MakeWithCapacity.h"

namespace utils {

template <typename VectorClass>
VectorClass make_with_capacity(typename VectorClass::size_type capacity) {
  VectorClass result;
  result.reserve(capacity);
  return result;
}

}  // namespace utils
