message(STATUS "Trying to check available SSE/AVX instructions")

include(CheckCXXCompilerFlag)

set(CMAKE_REQUIRED_LIBRARIES "-mavx2")
check_cxx_compiler_flag("-mavx2" AVX2_FOUND)

unset(CMAKE_REQUIRED_LIBRARIES)

if (AVX2_FOUND)
    message(STATUS "Found AVX2")

    add_library(AVX2 INTERFACE IMPORTED)
    set_target_properties(AVX2 PROPERTIES IMPORTED_GLOBAL TRUE)
    if(MSVC)
        target_compile_options(AVX2 INTERFACE /arch:AVX2)
    else ()
        target_compile_options(AVX2 INTERFACE -mavx)
    endif ()
    add_library(${PROJECT_NAME}::AVX2 ALIAS AVX2)
    mark_as_advanced(AVX2_FOUND)
endif ()

# todo: find SSE4.1 and AVX-512

