#pragma once

#include <Core/Types.h>

#include <algorithm>
#include <limits>
#include <random>

class Random final {
 public:
  explicit Random() : generator_{rd_()} {}

 public:
  [[nodiscard]] Precision Next(const Precision from_inclusive,
                               const Precision to_exclusive) {
    return std::uniform_real_distribution<Precision>(from_inclusive,
                                                     to_exclusive)(generator_);
  }

 private:
  std::random_device rd_;
  std::mt19937 generator_;
};

inline auto InitializeWithRandomNumbers(const size_t size) {
  auto array = Grid1D(size);
  Random random{};
  std::generate(array.begin(), array.end(), [&random] {
    return random.Next(static_cast<Precision>(-1.), static_cast<Precision>(1.));
  });
  return array;
}
