#pragma once

#include "Types.h"

//! \todo: Студия не дружит с cassert, но очень нужно обмазать код assert'ми
#include <cassert>

inline auto make_grid_2d(size_t z, size_t k, Precision default_value = 0.) {
  assert(z && k);
  return Grid2D(z, Grid1D(k, default_value));
}

inline auto make_grid_3d(size_t t, size_t z, size_t k,
                         Precision default_value = 0.) {
  assert(t && z && k);
  return Grid3D(t, Grid2D(z, Grid1D(k, default_value)));
}
