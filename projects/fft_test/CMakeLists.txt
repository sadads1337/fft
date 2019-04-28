cmake_minimum_required (VERSION 3.10)

add_executable (ffttest "main.cpp")

find_package(MKL REQUIRED)

target_link_libraries(ffttest 
	PRIVATE
		#FFT::Core
		MKL::MKL
)

# TODO: Add tests and install targets if needed.
