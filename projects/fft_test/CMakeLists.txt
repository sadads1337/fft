cmake_minimum_required (VERSION 3.10)

add_executable (ffttest "main.cpp")

target_link_libraries(ffttest 
	PUBLIC
		FFT::Core
)

# TODO: Add tests and install targets if needed.
