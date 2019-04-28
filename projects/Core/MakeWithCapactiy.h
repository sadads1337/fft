#pragma once
#include "MakeWithCapactiy.h"

namespace fft
{

template<typename VectorClass>
VectorClass make_with_capacity(typename VectorClass::size_type capacity);

} // namespace fft

#include "MakeWithCapactiy.inl"
