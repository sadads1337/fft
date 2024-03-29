set(MODULE_NAME Core)

add_subdirectory(aligned)

set(SRC_LIST
        MakeWithCapacity.h
        MakeWithCapacity.inl
        MoveOnly.h
        MoveOnly.inl
        MoveAndClear.h
        MoveAndClear.inl
        GridUtils.h
        TimeGuard.h
        Types.h

        src/TimeGuard.cpp

        MKL/Utils.h
        MKL/Utils.inl
        MKL/src/Utils.cpp
        )

add_library(${MODULE_NAME} ${SRC_LIST})

set_target_properties(${MODULE_NAME} PROPERTIES LINKER_LANGUAGE CXX)

find_package(MKL REQUIRED)

target_link_libraries(${MODULE_NAME}
        PUBLIC
        MKL::MKL
        aligned
        )

if (FFT_ENABLE_LUNA)
    target_compile_options(${MODULE_NAME} PRIVATE -fPIC)
endif ()

add_library(${PROJECT_NAME}::${MODULE_NAME} ALIAS ${MODULE_NAME})

if (FFT_ENABLE_TESTS)
    set(TEST_SRC_LIST
            test/ConvTest.cpp
            )

    add_executable(CoreTest ${TEST_SRC_LIST})

    find_package(gtest REQUIRED)

    target_link_libraries(${MODULE_NAME}Test
            PRIVATE
            gtest::gtest
            ${MODULE_NAME}
            )
    add_test(NAME ${MODULE_NAME}TestRunner COMMAND ${MODULE_NAME}Test)
endif ()

add_subdirectory(matplotlib-cpp)

add_library(${PROJECT_NAME}::matplotlibcpp ALIAS matplotlibcpp)

set(BENCHMARK_ENABLE_TESTING OFF)
set(BENCHMARK_ENABLE_GTEST_TESTS OFF)
set(BENCHMARK_ENABLE_INSTALL OFF)

add_subdirectory(benchmark)

if (FFT_ENABLE_LUNA)
    add_subdirectory(luna)
endif ()