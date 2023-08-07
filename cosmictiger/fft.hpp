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

void fft_1d_simd(sfmm::complex<fft_simd>* X, int N, int dir = +1, bool scrmbl = true);
const std::vector<fft_float>& cos_twiddles(int N);
const std::vector<fft_float>& sin_twiddles(int N);
const std::vector<fft_float>& inv_sin_twiddles(int N);

#endif /* FFT_HPP_ */
