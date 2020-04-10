#pragma once

#include "Types.h"

//! \todo: Студия не дружит с cassert, но очень нужно обмазать код assert'ми
#include <cassert>

inline auto make_grid_2d(const size_t z, const size_t k, const Precision default_value = 0.)
{
	assert(z && k);
	return Grid2D(z, Grid1D(k, default_value));
}

inline auto make_grid_3d(const size_t t, const size_t z, const size_t k, const Precision default_value = 0.)
{
	assert(t && z && k);
	return Grid3D(t, Grid2D(z, Grid1D(k, default_value)));
}
