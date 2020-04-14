add_executable(TestApp main.cpp)

target_link_libraries(TestApp
        PRIVATE
        FFT::Scheme
        FFT::matplotlibcpp
        )

# TODO: Add tests and install targets if needed.
