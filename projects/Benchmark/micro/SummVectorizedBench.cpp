#include <Core/Types.h>
#include <Scheme/src/IdxUtils.h>
#include <benchmark/benchmark.h>

#include "Utils.h"

static void summ_real(benchmark::State& state) {
  const auto size = static_cast<size_t>(state.range(0u));
  const auto lhs = InitializeWithRandomNumbers(size);
  const auto rhs = InitializeWithRandomNumbers(size);
  for ([[maybe_unused]] auto _ : state) {
    benchmark::DoNotOptimize(utils::summ_real<false>(lhs, rhs));
  }
}

static void summ_real_vectorized(benchmark::State& state) {
  const auto size = static_cast<size_t>(state.range(0u));
  const auto lhs = InitializeWithRandomNumbers(size);
  const auto rhs = InitializeWithRandomNumbers(size);
  for ([[maybe_unused]] auto _ : state) {
    benchmark::DoNotOptimize(utils::summ_real<true>(lhs, rhs));
  }
}

static auto customize_benchmark(
    benchmark::internal::Benchmark* const benchmark) {
  for (auto size = (1ull << 10ull); size <= (1ull << 24ull); size <<= 1ull) {
    benchmark->Arg(size);
  }
}

BENCHMARK(summ_real)->Apply(customize_benchmark);
BENCHMARK(summ_real_vectorized)->Apply(customize_benchmark);

BENCHMARK_MAIN();