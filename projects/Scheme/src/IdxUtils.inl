#pragma once

#include <Core/MakeWithCapacity.h>
#include <Core/MKL/Utils.h>

#include "OpenMP.h"

namespace utils
{

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

inline auto apply_conv_factor(const Grid1D & input, const Precision mult)
{
	//! \todo: remove redudant allocation
	auto result = Grid1D(input.size(), static_cast<Precision>(0.));
	for (auto k_idx = 0u; k_idx < input.size(); ++k_idx)
	{
		result[k_idx] = mult * static_cast<Precision>(k_idx) * input[k_idx];
	}
	return result;
}

inline auto apply_corr_factor(const Grid1D & input, const Precision mult)
{
	//! \todo: remove redudant allocation
	auto result = Grid1D(input.size(), static_cast<Precision>(0.));
	for (auto k_idx = 0u; k_idx < input.size(); ++k_idx)
	{
		result[k_idx] = mult * static_cast<Precision>(input.size() - k_idx) * input[k_idx];
	}
	return result;
}

template<bool vectorized>
auto inline summ_real(const Grid1D & x, const Grid1D & y)
{
	assert(x.size() == y.size());
	//! \todo: remove redudant allocation
	auto result = utils::make_with_capacity<Grid1D>(x.size());

	if constexpr (vectorized)
	{
		FFT_OMP_PRAGMA("omp parallel for")
		for (auto idx = 0u; idx < x.size(); ++idx) {
			result.emplace_back(x[idx] + y[idx]);
		}
	}
	else
	{
		for (auto idx = 0u; idx < x.size(); ++idx) {
			result.emplace_back(x[idx] + y[idx]);
		}
	}
	return result;
}

} // namespace utils
