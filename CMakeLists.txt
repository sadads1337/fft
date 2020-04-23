cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

project(FFT C CXX)

option(FFT_SANITIZE "Enable builds with sanitizers support" OFF)
option(FFT_ENABLE_TESTS "Enable tests build" OFF)
option(FFT_ENABLE_OPENMP "Enable OpenMP support" ON)
option(FFT_ENABLE_MANUAL_VECT "Enable manual vectorization" ON)
option(FFT_ENABLE_VECT_REPORT "Enable additional vectorization report while compiling" ON)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(common)


add_subdirectory(projects/Core)
add_subdirectory(projects/Scheme)
add_subdirectory(projects/TestApp)
add_subdirectory(projects/Benchmark)
