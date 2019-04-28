#include <iostream>
#include <vector>
#include <complex>
#include <algorithm>

#include <Core/MakeWithCapactiy.h>
#include <Core/MKL/Utils.h>

fft::RealContainer generate_grid(const float start, const float end, const std::uint32_t n)
{
	auto result = utils::make_with_capacity<fft::RealContainer>(n);
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

	auto fft_result = fft::fft_real(grid);
	const auto ifft_result = fft::ifft_real(fft_result);

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
