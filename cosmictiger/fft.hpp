/*
 * fft.hpp
 *
 *  Created on: Aug 3, 2023
 *      Author: dmarce1
 */

#ifndef FFT_HPP_
#define FFT_HPP_

#include <complex>
#include <vector>
#include <sfmm.hpp>

void transpose_yz(float* X, int N);
void transpose_r2c(float* xin, float* xout, int N);
void scramble(float* X, int N1, int N2, int N3);
const std::vector<float>& cos_twiddles(int N);
const std::vector<float>& sin_twiddles(int N);
void fft_real(float* X, int N, int N0);
void fft_complex(float* X, float* Y, int N, int N0);
void fft_3d(float* xin, sfmm::complex<float>* xout, int N);

#endif /* FFT_HPP_ */
