set(CMAKE_CXX_STANDARD 17)

if (NOT MSVC)
    add_compile_options(-Werror=all)

    if (CMAKE_BUILD_TYPE EQUAL "RELEASE")
        # Since our compiler is apple clang or icc
        # Enable O3 in release, hope there is no UB in thirdparty libs
        set (FFT_COMPILATION_FLAGS -O3 -march=native)

        if (FFT_ENABLE_OPENMP)
            list(APPEND FFT_COMPILATION_FLAGS -qopenmp-simd -ffast-math)
        endif ()

        add_compile_options(${FFT_COMPILATION_FLAGS})
    endif ()
else ()
    add_compile_options(/W4)

    if (CMAKE_BUILD_TYPE EQUAL "RELEASE")
        # Enable O3 in release, hope there is no UB in thirdparty libs
        add_compile_options(/O3)
    endif ()
endif ()

include_directories(${PROJECT_SOURCE_DIR}/Projects)

if (FFT_ENABLE_TESTS)
    enable_testing()
endif ()

if (FFT_SANITIZE)
    message(STATUS "Trying to check available sanitizers")

    include(CheckCXXCompilerFlag)

    set(CMAKE_REQUIRED_LIBRARIES "-fsanitize=undefined")
    check_cxx_compiler_flag("-fsanitize=undefined -fno-sanitize-recover=all" UBSAN_FOUND)

    set(CMAKE_REQUIRED_LIBRARIES "-fsanitize=thread")
    check_cxx_compiler_flag("-fsanitize=thread -fno-sanitize-recover=all" TSAN_FOUND)

    set(CMAKE_REQUIRED_LIBRARIES "-fsanitize=address")
    check_cxx_compiler_flag("-fsanitize=address -fno-sanitize-recover=all" ASAN_FOUND)

    unset(CMAKE_REQUIRED_LIBRARIES)

    if (UBSAN_FOUND)
        message(STATUS "UB sanitizer available. Adding UBsan build type")
        list(APPEND CMAKE_CONFIGURATION_TYPES UBsan)
        set(CMAKE_C_FLAGS_UBSAN "${CMAKE_C_FLAGS_DEBUG} -fsanitize=undefined -fno-sanitize-recover=all")
        set(CMAKE_CXX_FLAGS_UBSAN "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=undefined -fno-sanitize-recover=all")
        set(CMAKE_EXE_LINKER_FLAGS_UBSAN "${CMAKE_EXE_LINKER_FLAGS_DEBUG} -fsanitize=undefined -fno-sanitize-recover=all")
    endif ()

    if (TSAN_FOUND)
        message(STATUS "Thread sanitizer available. Adding Tsan build type")
        list(APPEND CMAKE_CONFIGURATION_TYPES Tsan)
        set(CMAKE_C_FLAGS_TSAN "${CMAKE_C_FLAGS_DEBUG} -fsanitize=thread -fno-sanitize-recover=all")
        set(CMAKE_CXX_FLAGS_TSAN "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=thread -fno-sanitize-recover=all")
        set(CMAKE_EXE_LINKER_FLAGS_TSAN "${CMAKE_EXE_LINKER_FLAGS_DEBUG} -fsanitize=thread -fno-sanitize-recover=all")
    endif ()

    if (ASAN_FOUND)
        message(STATUS "Address sanitizer available. Adding Asan build type")
        list(APPEND CMAKE_CONFIGURATION_TYPES Asan)
        set(CMAKE_C_FLAGS_ASAN "${CMAKE_C_FLAGS_DEBUG} -fsanitize=address -fno-sanitize-recover=all")
        set(CMAKE_CXX_FLAGS_ASAN "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address -fno-sanitize-recover=all")
        set(CMAKE_EXE_LINKER_FLAGS_ASAN "${CMAKE_EXE_LINKER_FLAGS_DEBUG} -fsanitize=address -fno-sanitize-recover=all")
    endif ()
endif ()

function (add_vectorization_report MODULE_NAME)
    if (FFT_ENABLE_VECT_REPORT)
        if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
            SET(FFT_COMPILATION_FLAGS_FOR_REPORT
                    -Rpass="loop|vect"
                    -Rpass-missed="loop|vect"
                    -Rpass-analysis="loop|vect"
                    )
        elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
            message(STATUS "GNU Compiler report for vectorization is not supported")
        elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
            SET(FFT_COMPILATION_FLAGS_FOR_REPORT
                    -qopt-report-file=stdout
                    -qopt-report-format=vs
                    -qopt-report=5
                    -qopt-report-phase=loop,vec
                    )
        elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
            message(STATUS "MSVC Compiler report for vectorization is not supported")
        endif()

        if (FFT_COMPILATION_FLAGS_FOR_REPORT)
            target_compile_options(${MODULE_NAME}
                    PRIVATE
                    ${FFT_COMPILATION_FLAGS_FOR_REPORT}
                    )
        endif ()
    endif ()
endfunction ()
