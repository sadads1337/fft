#include <benchmark/benchmark.h>

#include <Core/Types.h>
#include <Scheme/src/IdxUtils.h>

#include <random>
#include <algorithm>
#include <limits>

namespace
{

class Random final
{
public:
	explicit Random()
		: generator_{rd_()}
	{
	}

public:
	[[nodiscard]] Precision Next(const Precision from_inclusive, const Precision to_exclusive)
	{
		return std::uniform_real_distribution<Precision>(from_inclusive, to_exclusive)(generator_);
	}

private:
	std::random_device rd_;
	std::mt19937 generator_;
};

inline auto InitializeWithRandomNumbers(const size_t size)
{
	auto array = Grid1D(size);
	Random random{};
	std::generate(
		array.begin(),
		array.end(),
		[&random] {
			return random.Next(static_cast<Precision>(-1.), static_cast<Precision>(1.));});
	return array;
}

} // namespace

static void summ_real(benchmark::State& state)
{
	const auto size = static_cast<size_t>(state.range(0u));
	const auto lhs = InitializeWithRandomNumbers(size);
	const auto rhs = InitializeWithRandomNumbers(size);
	for ([[maybe_unused]] auto _ : state)
	{
		benchmark::DoNotOptimize(utils::summ_real<false>(lhs, rhs));
	}
}

static void summ_real_vectorized(benchmark::State& state)
{
	const auto size = static_cast<size_t>(state.range(0u));
	const auto lhs = InitializeWithRandomNumbers(size);
	const auto rhs = InitializeWithRandomNumbers(size);
	for ([[maybe_unused]] auto _ : state)
	{
		benchmark::DoNotOptimize(utils::summ_real<true>(lhs, rhs));
	}
}

static auto customize_benchmark(benchmark::internal::Benchmark * const benchmark)
{
	for (auto size = (1ull << 10ull); size <= (1ull << 24ull); size <<= 1ull)
	{
		benchmark->Arg(size);
	}
}

BENCHMARK(summ_real)->Apply(customize_benchmark);
BENCHMARK(summ_real_vectorized)->Apply(customize_benchmark);

BENCHMARK_MAIN();