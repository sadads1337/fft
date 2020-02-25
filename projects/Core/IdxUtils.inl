#pragma once

#include <cstdint>

#include <Core/MakeWithCapacity.h>
#include <Core/MKL/Utils.h>

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

inline auto apply_operation(
	const Grid1D & lhs,
	const Grid1D & rhs,
	const std::function<Precision(Precision, Precision)> & function)
{
	assert(lhs.size() && rhs.size());
	auto result = utils::make_with_capacity<fft::RealContainer>(lhs.size());
	auto it_x = lhs.begin();
	auto it_y = rhs.begin();
	assert(lhs.size() == rhs.size());
	for (; it_x != lhs.end(); ++it_x, ++it_y)
	{
		result.push_back(function(*it_x, *it_y));
	}
	assert(result.size() == lhs.size());
	return result;
}

inline auto apply_conv_factor(const Grid1D & input)
{
	auto result = Grid1D(input.size(), 0.);
	for (auto k_idx = 0u; k_idx < input.size(); ++k_idx)
	{
		//! \todo pi / b
		result[k_idx] = k_idx * input[k_idx];
	}
	return result;
}

inline auto apply_corr_factor(const Grid1D & input)
{
	auto result = Grid1D(input.size(), 0.);
	for (auto k_idx = 0u; k_idx < input.size(); ++k_idx)
	{
		result[k_idx] = k_idx * input[k_idx];
	}
	return result;
}

} // namespace utils
