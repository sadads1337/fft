#pragma once

#include <chrono>

namespace utils {

class TimeGuard {
 public:
  TimeGuard();
  ~TimeGuard();

 private:
  std::chrono::time_point<std::chrono::high_resolution_clock> start_;
};

}  // namespace utils
