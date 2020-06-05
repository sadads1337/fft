add_executable(AFSap main.cpp)

if (FFT_ENABLE_OPENMP)
    message(STATUS "Trying to find OpenMP")
    if (APPLE)
        if (CMAKE_C_COMPILER_ID MATCHES "Clang")
            set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp")
            set(OpenMP_C_LIB_NAMES "omp")
            set(OpenMP_omp_LIBRARY omp)
        endif ()
        if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
            set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp")
            set(OpenMP_CXX_LIB_NAMES "omp")
            set(OpenMP_omp_LIBRARY omp)
        endif ()
        target_include_directories(AFSap PUBLIC /usr/local/opt/libomp/include)
        target_link_directories(AFSap PUBLIC /usr/local/opt/libomp/lib)
    endif ()
    find_package(OpenMP REQUIRED)
    target_link_libraries(AFSap
            PUBLIC
            OpenMP::OpenMP_CXX
            )
endif ()
