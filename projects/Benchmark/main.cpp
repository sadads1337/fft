#include <Scheme/Scheme.h>

#include <Core/Types.h>
#include <Core/MKL/Utils.h>

#include <iostream>
#include <chrono>

namespace
{

class TimeGuard
{
public:
	TimeGuard()
		: start_{std::chrono::high_resolution_clock::now()}
	{
	}

	~TimeGuard()
	{
		const auto stop = std::chrono::high_resolution_clock::now();
		const auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(stop - start_);
		std::cout << "Elapsed in " << elapsed.count() << "seconds";
	}

private:
	std::chrono::time_point<std::chrono::high_resolution_clock> start_;
};

} // namespace

int main() try
{
	scheme::Values values_1{};
	scheme::Values values_2{};
	const scheme::Env env{};

	{
		TimeGuard time_guard{};
		for (auto t_idx = 0u; t_idx < scheme::g_t_grid_size; ++t_idx) {
			calculate_one_step(values_1, values_2, env, t_idx);

			//! new in values_2 now; we don't need values_2, swap here, do not copy
			std::swap(values_1, values_2);
		}
	}
}
catch(const utils::MKLException & exception)
{
	std::cout << "MKL exception happend: " << exception.what();
}
catch(const std::exception & exception)
{
	std::cout << "Another exception happend: " << exception.what();
}
