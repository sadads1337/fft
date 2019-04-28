cmake_minimum_required (VERSION 3.10)

set(MODULE_NAME Core)

set(SRC_LIST 
	"MakeWithCapactiy.h"
	"MakeWithCapactiy.inl"
	"MoveOnly.h"
	"MoveOnly.inl"
	"MoveAndClear.h"
	"MoveAndClear.inl"
)

enable_language(CXX)

add_library(${MODULE_NAME} ${SRC_LIST})
set_target_properties(${MODULE_NAME} PROPERTIES LINKER_LANGUAGE CXX)

add_library(${PROJECT_NAME}::${MODULE_NAME} ALIAS ${MODULE_NAME})
