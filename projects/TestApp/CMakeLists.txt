cmake_minimum_required (VERSION 3.10)

add_executable (TestApp "main.cpp")

target_link_libraries(TestApp 
	PUBLIC
		FFT::Core
)

# TODO: Add tests and install targets if needed.
