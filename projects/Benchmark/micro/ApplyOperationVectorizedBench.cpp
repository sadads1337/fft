#include <benchmark/benchmark.h>

#include <Core/Types.h>
#include <Scheme/src/IdxUtils.h>

#include "Utils.h"

namespace
{

auto f(const Precision lhs, const Precision rhs)
{
	constexpr auto factor = static_cast<Precision>(2.);
	return (lhs + static_cast<Precision>(rhs)) / factor;
}

} // namespace

static void apply_operation(benchmark::State& state)
{
	const auto size = static_cast<size_t>(state.range(0u));
	const auto lhs = InitializeWithRandomNumbers(size);
	const auto rhs = InitializeWithRandomNumbers(size);
	for ([[maybe_unused]] auto _ : state)
	{
		benchmark::DoNotOptimize(utils::apply_operation<false>(lhs, rhs, f));
	}
}

static void apply_operation_vectorized(benchmark::State& state)
{
	const auto size = static_cast<size_t>(state.range(0u));
	const auto lhs = InitializeWithRandomNumbers(size);
	const auto rhs = InitializeWithRandomNumbers(size);
	for ([[maybe_unused]] auto _ : state)
	{
		benchmark::DoNotOptimize(utils::apply_operation<true>(lhs, rhs, f));
	}
}

static auto customize_benchmark(benchmark::internal::Benchmark * const benchmark)
{
	for (auto size = (1ull << 10ull); size <= (1ull << 24ull); size <<= 1ull)
	{
		benchmark->Arg(size);
	}
}

BENCHMARK(apply_operation)->Apply(customize_benchmark);
BENCHMARK(apply_operation_vectorized)->Apply(customize_benchmark);

BENCHMARK_MAIN();