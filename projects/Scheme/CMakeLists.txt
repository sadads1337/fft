set(MODULE_NAME Scheme)

set(SRC_LIST
	Scheme.h
	src/Scheme.cpp
)

add_library(${MODULE_NAME} ${SRC_LIST})

target_link_libraries(${MODULE_NAME}
	PUBLIC
		FFT::Core
)

add_library(${PROJECT_NAME}::${MODULE_NAME} ALIAS ${MODULE_NAME})
