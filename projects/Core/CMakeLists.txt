cmake_minimum_required (VERSION 3.10)

set(MODULE_NAME Core)

set(SRC_LIST 
	"MakeWithCapacity.h"
	"MakeWithCapacity.inl"
	"MoveOnly.h"
	"MoveOnly.inl"
	"MoveAndClear.h"
	"MoveAndClear.inl"
	
	"MKL/Utils.h"
	"MKL/Utils.inl"
	"MKL/src/Utils.cpp"
)


add_library(${MODULE_NAME} ${SRC_LIST})
set_target_properties(${MODULE_NAME} PROPERTIES LINKER_LANGUAGE CXX)

find_package(MKL REQUIRED)

target_link_libraries(${MODULE_NAME}
	PUBLIC
		MKL::MKL
)

target_include_directories(${MODULE_NAME}
	PUBLIC
		${PROJECT_SOURCE_DIR}/Projects/
)

add_library(${PROJECT_NAME}::${MODULE_NAME} ALIAS ${MODULE_NAME})


set(TEST_SRC_LIST 
	"test/conv_test.cpp"
)

add_executable(conv_test ${TEST_SRC_LIST})
target_link_libraries(conv_test 
	gtest_main
	${MODULE_NAME}
)
add_test(NAME conv_test_runner COMMAND conv_test)
