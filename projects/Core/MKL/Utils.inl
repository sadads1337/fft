#pragma once

#include "Utils.h"

namespace utils
{

template<typename Function, typename... Args>
inline MKL_LONG save_mkl_call(Function & function, Args && ...args)
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
