#include <Core/MKL/Utils.h>
#include <Core/MakeWithCapacity.h>
#include <mkl.h>

namespace fft {

ComplexContainer fft_complex(ComplexContainer& in) {
  ComplexContainer out(in.size());

  DFTI_DESCRIPTOR_HANDLE descriptor;
  utils::save_mkl_call(DftiCreateDescriptor, &descriptor, DFTI_SINGLE,
                       DFTI_COMPLEX, 1, in.size());
  utils::save_mkl_call(DftiSetValue, descriptor, DFTI_PLACEMENT,
                       DFTI_NOT_INPLACE);
  utils::save_mkl_call(DftiCommitDescriptor, descriptor);
  utils::save_mkl_call(DftiComputeForward, descriptor, in.data(), out.data());
  utils::save_mkl_call(DftiFreeDescriptor, &descriptor);

  return out;
}

ComplexContainer fft_real(const RealContainer& in_real) {
  auto in = utils::make_with_capacity<ComplexContainer>(in_real.size());
  //! \todo: потенциально тяжелая операцция копирования, fixme!
  std::copy(in_real.begin(), in_real.end(), std::back_inserter(in));
  return fft_complex(in);
}

RealContainer ifft_real(ComplexContainer& in) {
  ComplexContainer out(in.size());

  DFTI_DESCRIPTOR_HANDLE descriptor;
  utils::save_mkl_call(DftiCreateDescriptor, &descriptor, DFTI_SINGLE,
                       DFTI_COMPLEX, 1, in.size());
  utils::save_mkl_call(DftiSetValue, descriptor, DFTI_PLACEMENT,
                       DFTI_NOT_INPLACE);
  utils::save_mkl_call(DftiSetValue, descriptor, DFTI_BACKWARD_SCALE,
                       1.0f / in.size());
  utils::save_mkl_call(DftiCommitDescriptor, descriptor);
  utils::save_mkl_call(DftiComputeBackward, descriptor, in.data(), out.data());
  utils::save_mkl_call(DftiFreeDescriptor, &descriptor);

  auto result = utils::make_with_capacity<RealContainer>(out.size());
  for (const auto& value : out) {
    result.emplace_back(value.real());
  }
  return result;
}

RealContainer conv_real(const RealContainer& x, const RealContainer& y) {
  RealContainer z(x.size() + y.size() - 1);
  VSLConvTaskPtr task;
  utils::save_mkl_call(vsldConvNewTask1D, &task, VSL_CONV_MODE_AUTO,
                       static_cast<int>(x.size()), static_cast<int>(y.size()),
                       static_cast<int>(z.size()));
  utils::save_mkl_call(vsldConvExec1D, task, x.data(), 1, y.data(), 1, z.data(),
                       1);
  utils::save_mkl_call(vslConvDeleteTask, &task);
  return z;
}

RealContainer corr_real(const RealContainer& x, const RealContainer& y) {
  RealContainer z(x.size() + y.size() - 1);
  VSLCorrTaskPtr task;
  utils::save_mkl_call(vsldCorrNewTask1D, &task, VSL_CORR_MODE_AUTO,
                       static_cast<int>(x.size()), static_cast<int>(y.size()),
                       static_cast<int>(z.size()));
  utils::save_mkl_call(vsldCorrExec1D, task, x.data(), 1, y.data(), 1, z.data(),
                       1);
  utils::save_mkl_call(vslCorrDeleteTask, &task);
  return z;
}

}  // namespace fft
