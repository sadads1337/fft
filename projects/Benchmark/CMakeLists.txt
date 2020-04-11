add_executable(Benchmark main.cpp)

target_link_libraries(Benchmark
	PRIVATE
		FFT::Scheme
)

add_executable(SummVectorizedBenchmark micro/SummVectorizedBench.cpp)

target_link_libraries(SummVectorizedBenchmark
	PUBLIC
		FFT::Scheme
		benchmark::benchmark
)
