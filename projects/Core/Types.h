#pragma once

#include <aligned/aligned_vector>
#include <type_traits>
#include <vector>

using Precision = double;
static_assert(std::is_floating_point_v<Precision>,
              "Precision must be float or double");

using Grid1D = aligned::aligned_vector<Precision, aligned::alignment::avx>;
using Grid2D = std::vector<Grid1D>;
using Grid3D = std::vector<Grid2D>;
