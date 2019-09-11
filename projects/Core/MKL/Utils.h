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

using ComplexContainer = std::vector<std::complex<double>>;
using RealContainer = std::vector<double>;

ComplexContainer fft_complex(ComplexContainer & in);

ComplexContainer fft_real(const RealContainer & in_real);

RealContainer ifft_real(ComplexContainer & in);

RealContainer conv_real(const RealContainer & x, const RealContainer & y);

RealContainer corr_real(const RealContainer & x, const RealContainer & y);

RealContainer summ_real(const RealContainer & x, const RealContainer & y);

} // namespace fft

#include "Utils.inl"
