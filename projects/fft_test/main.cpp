#include <iostream>
#include <vector>
#include <complex>
#include <algorithm>
#include <numeric>

#include <mkl.h>

namespace utils
{

template<typename VectorClass>
VectorClass make_with_capacity(const typename VectorClass::size_type capacity)
{
	VectorClass result;
	result.reserve(capacity);
	return result;
}

struct MKLException final : std::runtime_error
{
	explicit MKLException(const std::string & message)
		: std::runtime_error(message)
	{
	}
	virtual ~MKLException() override = default;
};

template<typename Function, typename... Args>
MKL_LONG save_mkl_call(Function & function, Args && ...args)
{
	//! \todo: MKL принимает входные данные не по константному указателю, и в теории может поменять входные данные
	//! \todo: WTF???, нужно это как-то обойти.
	const auto status = function(std::forward<Args>(args)...);
	if (status)
	{
		throw MKLException(DftiErrorMessage(status));
	}
	return status;
}

} // namespace utils

using ComplexContainer = std::vector<std::complex<float>>;
using RealContainer = std::vector<float>;

RealContainer generate_grid(const float start, const float end, const std::uint32_t n)
{
	auto result = utils::make_with_capacity<RealContainer>(n);
	for (auto i = 0u; i <= n; ++i)
	{
		const auto value = start + (end - start) / n * i;
		result.emplace_back(value);
	}
	return result;
}

float gauss_function(const float argument) noexcept
{
	return std::exp(argument * argument * -2.f);
}

float step_function(const float argument) noexcept
{
	return argument < 0.f ? 0.f : 1.f;
}

ComplexContainer fft_complex(ComplexContainer & in)
{
	ComplexContainer out(in.size());

	DFTI_DESCRIPTOR_HANDLE descriptor;
	utils::save_mkl_call(DftiCreateDescriptor, &descriptor, DFTI_SINGLE, DFTI_COMPLEX, 1, in.size());
	utils::save_mkl_call(DftiSetValue, descriptor, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
	utils::save_mkl_call(DftiCommitDescriptor, descriptor);
	utils::save_mkl_call(DftiComputeForward, descriptor, in.data(), out.data());
	utils::save_mkl_call(DftiFreeDescriptor, &descriptor);

	return out;
}

ComplexContainer fft_real(const RealContainer & in_real)
{
	auto in = utils::make_with_capacity<ComplexContainer>(in_real.size());
	//! \todo: потенциально тяжелая операцция копирования, fixme!
	std::copy(in_real.begin(), in_real.end(), std::back_inserter(in));
	return fft_complex(in);
}

RealContainer ifft_real(ComplexContainer & in)
{
	std::vector<std::complex<float>> out(in.size());

	DFTI_DESCRIPTOR_HANDLE descriptor;
	utils::save_mkl_call(DftiCreateDescriptor, &descriptor, DFTI_SINGLE, DFTI_COMPLEX, 1, in.size());
	utils::save_mkl_call(DftiSetValue, descriptor, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
	utils::save_mkl_call(DftiSetValue, descriptor, DFTI_BACKWARD_SCALE, 1.0f / in.size());
	utils::save_mkl_call(DftiCommitDescriptor, descriptor);
	utils::save_mkl_call(DftiComputeBackward, descriptor, in.data(), out.data());
	utils::save_mkl_call(DftiFreeDescriptor, &descriptor);

	auto result = utils::make_with_capacity<std::vector<float>>(out.size());
	for(const auto & value : out)
	{
		result.emplace_back(value.real());
	}
	return result;
}

int main() try
{
	auto grid = generate_grid(-3.14f, 3.14f, 512);
	std::transform(
		grid.begin(),
		grid.end(),
		grid.begin(),
		[](const auto & value)
		{
			//return gauss_function(value);
			return step_function(value);
		});

	auto fft_result = fft_real(grid);
	const auto ifft_result = ifft_real(fft_result);

	auto it1 = grid.begin();
	auto it2 = ifft_result.begin();
	auto max_error = 0.f;

	for (; it1 != grid.end(); ++it1, ++it2)
	{
		const auto error = std::abs(*it1 - *it2);
		if (error > max_error)
		{
			max_error = error;
		}
	}
	std::cout << "Max error: " << max_error << std::endl;
}
catch(const utils::MKLException & mkl_exception)
{
	std::cout << "MKL Exception happend: " << mkl_exception.what();
}
catch(const std::exception & exception)
{
	std::cout << "Exception happend: " << exception.what();
}
