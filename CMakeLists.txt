cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

project(FFT C CXX)

option(FFT_SANITIZE "Enable builds with sanitizers support" OFF)
option(FFT_ENABLE_TESTS "Enable tests build" ON)
option(FFT_ENABLE_OPENMP "Enable OpenMP support" ON)
option(FFT_ENABLE_LUNA "Enable LuNA support" OFF)
option(FFT_ENABLE_MANUAL_VECT "Enable manual vectorization" OFF)
option(FFT_ENABLE_VECT_REPORT "Enable additional vectorization report while compiling" OFF)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(common)

add_subdirectory(projects/Core)
add_subdirectory(projects/Scheme)
add_subdirectory(projects/TestApp)
add_subdirectory(projects/Benchmark)
if (FFT_ENABLE_LUNA)
    add_subdirectory(projects/Luna)
endif ()

install(TARGETS
        Benchmark
        TestApp
        Scheme
        ConvVectorizedBenchmark
        CorrVectorizedBenchmark
        ApplyOperationVectorizedBenchmark
        )

# Package generation
set(CPACK_SOURCE_GENERATOR "TGZ;ZIP")
set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})
include(CPack)
