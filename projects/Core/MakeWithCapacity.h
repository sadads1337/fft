#pragma once
#include "MakeWithCapacity.h"

namespace utils {

template <typename VectorClass>
VectorClass make_with_capacity(typename VectorClass::size_type capacity);

}  // namespace utils

#include "MakeWithCapacity.inl"
