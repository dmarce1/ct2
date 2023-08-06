/*
 * cosmictiger.hpp
 *
 *  Created on: Jul 26, 2023
 *      Author: dmarce1
 */
#pragma once

#include <sfmm.hpp>


#define NDIM 3
#define FFT_SIMD_SIZE 4
#define DRIFT_SIMD_SIZE 8
#define FORCE_SIMD_SIZE 8
#define PORDER 8

using fft_float = double;
using fft_simd = sfmm::simd_f64;
using force_float = float;
using force_simd = sfmm::simd_f32;
using pos_fixed = sfmm::fixed32;
using pos_simd = sfmm::simd_fixed32;

