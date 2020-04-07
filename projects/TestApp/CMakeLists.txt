add_executable(TestApp main.cpp)

target_link_libraries(TestApp 
	PUBLIC
		FFT::Core
		matplotlibcpp
)

# TODO: Add tests and install targets if needed.
