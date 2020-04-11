#include <benchmark/benchmark.h>

#include <Core/Types.h>
#include <Scheme/src/IdxUtils.h>

#include "Utils.h"

static void apply_conv_factor(benchmark::State& state)
{
	const auto size = static_cast<size_t>(state.range(0u));
	const auto lhs = InitializeWithRandomNumbers(size);
	constexpr auto factor = static_cast<Precision>(1.5);
	for ([[maybe_unused]] auto _ : state)
	{
		benchmark::DoNotOptimize(utils::apply_conv_factor<false>(lhs, factor));
	}
}

static void apply_conv_factor_vectorized(benchmark::State& state)
{
	const auto size = static_cast<size_t>(state.range(0u));
	const auto lhs = InitializeWithRandomNumbers(size);
	constexpr auto factor = static_cast<Precision>(1.5);
	for ([[maybe_unused]] auto _ : state)
	{
		benchmark::DoNotOptimize(utils::apply_conv_factor<true>(lhs, factor));
	}
}

static auto customize_benchmark(benchmark::internal::Benchmark * const benchmark)
{
	for (auto size = (1ull << 10ull); size <= (1ull << 24ull); size <<= 1ull)
	{
		benchmark->Arg(size);
	}
}

BENCHMARK(apply_conv_factor)->Apply(customize_benchmark);
BENCHMARK(apply_conv_factor_vectorized)->Apply(customize_benchmark);

BENCHMARK_MAIN();