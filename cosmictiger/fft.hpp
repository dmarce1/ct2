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

void transpose_yz(float* X, int N);
void transpose_xy(float* X, int N);
void scramble(float* X, int N);
const std::vector<float>& cos_twiddles(int N);
const std::vector<float>& sin_twiddles(int N);
void fft(float* X, float* Y, int N);

#endif /* FFT_HPP_ */
