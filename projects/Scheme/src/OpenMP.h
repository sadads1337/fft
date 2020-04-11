#pragma once

#ifndef FFT_ENABLE_OPENMP
	#define FFT_ENABLE_OPENMP 0
#endif // #ifndef FFT_ENABLE_OPENMP

#ifdef FFT_OMP
	#undef FFT_OMP
#endif // #ifdef FFT_OMP

#if FFT_ENABLE_OPENMP
	#include <omp.h>
	#define FFT_OMP(omp_directive) omp_directive
#else
	#define FFT_OMP(omp_directive)
#endif // #ifdef FFT_ENABLE_OPENMP
