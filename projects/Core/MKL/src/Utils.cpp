#include <Core/MKL/Utils.h>
#include <Core/MakeWithCapactiy.h>

#include <mkl.h>

namespace fft
{

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
	for (const auto & value : out)
	{
		result.emplace_back(value.real());
	}
	return result;
}

} // namespace fft
