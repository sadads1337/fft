#include <Core/TimeGuard.h>

#include <iostream>

namespace utils {

TimeGuard::TimeGuard() : start_{std::chrono::high_resolution_clock::now()} {}

TimeGuard::~TimeGuard() {
  const auto stop = std::chrono::high_resolution_clock::now();
  const auto elapsed =
      std::chrono::duration_cast<std::chrono::seconds>(stop - start_);
  std::cout << "Elapsed in " << elapsed.count() << "seconds." << std::endl;
}

}  // namespace utils
