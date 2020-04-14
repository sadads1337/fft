set(MODULE_NAME Scheme)

set(SRC_LIST
        Scheme.h
        src/IdxUtils.h
        src/IdxUtils.inl
        src/OpenMP.h
        src/Scheme.cpp
        )

add_library(${MODULE_NAME} ${SRC_LIST})

if (FFT_ENABLE_OPENMP)
    message(STATUS "Trying to find OpenMP")
    find_package(OpenMP REQUIRED)
    target_link_libraries(${MODULE_NAME}
            PUBLIC
            OpenMP::OpenMP_CXX
            )
    target_compile_definitions(${MODULE_NAME} PUBLIC FFT_ENABLE_OPENMP)
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
