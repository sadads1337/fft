#pragma once

namespace utils
{

inline auto get_idx(const Precision l_edge, const Precision r_edge, const Precision value, const size_t grid_size)
{
	assert(l_edge < r_edge);
	assert(l_edge <= value && value <= r_edge);
	const auto idx = static_cast<size_t>(std::floor(value / ((r_edge - l_edge) / grid_size)));
	assert(0 <= idx && idx <= grid_size);
	return idx;
}

} // namespace utils
