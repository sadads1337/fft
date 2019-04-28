#pragma once

#include <vector>
#include <complex>

#include <mkl_types.h>


namespace utils
{

struct MKLException final : std::runtime_error
{
	explicit MKLException(const std::string & message)
		: std::runtime_error(message)
	{
	}
	virtual ~MKLException() override = default;
};

template<typename Function, typename... Args>
MKL_LONG save_mkl_call(Function & function, Args && ...args);
	
} // namespace utils

namespace fft
{

using ComplexContainer = std::vector<std::complex<float>>;
using RealContainer = std::vector<float>;

ComplexContainer fft_complex(ComplexContainer & in);

ComplexContainer fft_real(const RealContainer & in_real);

RealContainer ifft_real(ComplexContainer & in);
	
} // namespace fft

#include "Utils.inl"
