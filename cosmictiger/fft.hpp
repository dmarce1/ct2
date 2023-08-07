/*
 * fft.hpp
 *
 *  Created on: Aug 3, 2023
 *      Author: dmarce1
 */

#ifndef FFT_HPP_
#define FFT_HPP_

#include <cosmictiger/cosmictiger.hpp>
#include <complex>
#include <vector>
#include <sfmm.hpp>

using namespace sfmm;

void transpose_yz(fft_float* X, int N);
void transpose_r2c(fft_float* xin, fft_float* xout, int N);
void transpose_c2r(fft_float* xin, fft_float* xout, int N);
void scramble(fft_float* X, int N1, int N2, int N3);
const std::vector<fft_float>& cos_twiddles(int N);
const std::vector<fft_float>& sin_twiddles(int N);
const std::vector<fft_float>& inv_sin_twiddles(int N);
void fft_real(fft_float* X, int N, int dir = 1);
void fft_real(fft_simd* X, int N, int dir = 1);
void fft_complex(complex<fft_simd>* X, int N, int dir = 1);
void fft_r2c_3d(fft_float* xin, fft_float* xout, int N);
void fft_c2r_3d(fft_float* xin, fft_float* xout, int N);

#endif /* FFT_HPP_ */
