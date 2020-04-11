set(CMAKE_CXX_STANDARD 17)

if (NOT MSVC)
    # Since our compiler is apple-clang
    add_compile_options(-Werror=all)
endif()

include_directories(${PROJECT_SOURCE_DIR}/Projects)

enable_testing()
