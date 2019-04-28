cmake_minimum_required (VERSION 3.10)

project (FFT C CXX)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

add_subdirectory(projects/Core)
add_subdirectory(projects/fft_test)
