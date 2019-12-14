cmake_minimum_required (VERSION 3.10)

set(MODULE_NAME Core)

set(SRC_LIST 
	"MakeWithCapacity.h"
	"MakeWithCapacity.inl"
	"MoveOnly.h"
	"MoveOnly.inl"
	"MoveAndClear.h"
	"MoveAndClear.inl"
	"IdxUtils.h"
	"IdxUtils.inl"
	"Types.h"
	
	"MKL/Utils.h"
	"MKL/Utils.inl"
	"MKL/src/Utils.cpp"
	
	"matplotlib-cpp/matplotlibcpp.h"
)


add_library(${MODULE_NAME} ${SRC_LIST})
set_target_properties(${MODULE_NAME} PROPERTIES LINKER_LANGUAGE CXX)

find_package(MKL REQUIRED)
find_package(Python2 COMPONENTS Development NumPy REQUIRED)

target_link_libraries(${MODULE_NAME}
	PUBLIC
		MKL::MKL
		Python2::Python
		Python2::NumPy
)

target_include_directories(${MODULE_NAME}
	PUBLIC
		${PROJECT_SOURCE_DIR}/Projects/
		${Python2_INCLUDE_DIRS}
		${Python2_NumPy_INCLUDE_DIRS}
)

add_library(${PROJECT_NAME}::${MODULE_NAME} ALIAS ${MODULE_NAME})


set(TEST_SRC_LIST 
	"test/ConvTest.cpp"
	"test/IdxUtilsTest.cpp"
)

add_executable(CoreTest ${TEST_SRC_LIST})
target_link_libraries(CoreTest
	gtest_main
	${MODULE_NAME}
)
add_test(NAME CoreTestRunner COMMAND CoreTest)
