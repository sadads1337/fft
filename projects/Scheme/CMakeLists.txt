set(MODULE_NAME Scheme)

set(SRC_LIST
        Scheme.h
        src/IdxUtils.h
        src/IdxUtils.inl
        src/IdxUtilsAVX.inl
        src/OpenMP.h
        src/Scheme.cpp
        )

add_library(${MODULE_NAME} ${SRC_LIST})

add_vectorization_report(${MODULE_NAME})

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
        target_include_directories(${MODULE_NAME} PUBLIC /usr/local/opt/libomp/include)
        target_link_directories(${MODULE_NAME} PUBLIC /usr/local/opt/libomp/lib)
    endif ()
    find_package(OpenMP REQUIRED)
    target_link_libraries(${MODULE_NAME}
            PUBLIC
            OpenMP::OpenMP_CXX
            )
    target_compile_definitions(${MODULE_NAME} PUBLIC FFT_ENABLE_OPENMP)
endif ()

if (FFT_ENABLE_MANUAL_VECT)
    find_package(SSE REQUIRED)
    target_link_libraries(${MODULE_NAME}
            PUBLIC
            FFT::AVX2
            )
    target_compile_definitions(${MODULE_NAME} PUBLIC FFT_ENABLE_MANUAL_VECT)
endif ()

target_link_libraries(${MODULE_NAME}
        PUBLIC
        FFT::Core
        )

add_library(${PROJECT_NAME}::${MODULE_NAME} ALIAS ${MODULE_NAME})

if (FFT_ENABLE_TESTS)
    set(TEST_SRC_LIST
            test/IdxUtilsTest.cpp
            )

    add_executable(${MODULE_NAME}Test ${TEST_SRC_LIST})

    find_package(gtest REQUIRED)

    target_link_libraries(${MODULE_NAME}Test
            PRIVATE
            gtest::gtest
            ${MODULE_NAME}
            )
    add_test(NAME ${MODULE_NAME}TestRunner COMMAND ${MODULE_NAME}Test)
endif ()
