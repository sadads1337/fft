add_executable(Benchmark main.cpp)

target_link_libraries(Benchmark
        PRIVATE
        FFT::Scheme
        )

#todo: Add function to common.cmake

add_executable(SummVectorizedBenchmark micro/SummVectorizedBench.cpp micro/Utils.h)

target_link_libraries(SummVectorizedBenchmark
        PUBLIC
        FFT::Scheme
        benchmark::benchmark
        )

add_executable(ConvVectorizedBenchmark micro/ConvVectorizedBench.cpp micro/Utils.h)

target_link_libraries(ConvVectorizedBenchmark
        PUBLIC
        FFT::Scheme
        benchmark::benchmark
        )

add_executable(CorrVectorizedBenchmark micro/CorrVectorizedBench.cpp micro/Utils.h)

target_link_libraries(CorrVectorizedBenchmark
        PUBLIC
        FFT::Scheme
        benchmark::benchmark
        )

add_executable(ApplyOperationVectorizedBenchmark micro/ApplyOperationVectorizedBench.cpp micro/Utils.h)

target_link_libraries(ApplyOperationVectorizedBenchmark
        PUBLIC
        FFT::Scheme
        benchmark::benchmark
        )
