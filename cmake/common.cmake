set(CMAKE_CXX_STANDARD 17)

if (NOT MSVC)
    # Since our compiler is apple-clang
    add_compile_options(-Werror=all)
endif()

include_directories(${PROJECT_SOURCE_DIR}/Projects)

enable_testing()

option(FFT_SANITIZE On "Enable builds with sanitizers support")

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
    endif()

    if (TSAN_FOUND)
        message(STATUS "Thread sanitizer available. Adding Tsan build type")
        list(APPEND CMAKE_CONFIGURATION_TYPES Tsan)
        set(CMAKE_C_FLAGS_TSAN "${CMAKE_C_FLAGS_DEBUG} -fsanitize=thread -fno-sanitize-recover=all")
        set(CMAKE_CXX_FLAGS_TSAN "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=thread -fno-sanitize-recover=all")
        set(CMAKE_EXE_LINKER_FLAGS_TSAN "${CMAKE_EXE_LINKER_FLAGS_DEBUG} -fsanitize=thread -fno-sanitize-recover=all")
    endif()

    if (ASAN_FOUND)
        message(STATUS "Address sanitizer available. Adding Asan build type")
        list(APPEND CMAKE_CONFIGURATION_TYPES Asan)
        set(CMAKE_C_FLAGS_ASAN "${CMAKE_C_FLAGS_DEBUG} -fsanitize=address -fno-sanitize-recover=all")
        set(CMAKE_CXX_FLAGS_ASAN "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address -fno-sanitize-recover=all")
        set(CMAKE_EXE_LINKER_FLAGS_ASAN "${CMAKE_EXE_LINKER_FLAGS_DEBUG} -fsanitize=address -fno-sanitize-recover=all")
    endif()
endif()
