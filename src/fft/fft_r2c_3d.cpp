#include "cosmictiger/fft.hpp"
#include <hpx/async.hpp>

void fft_r2c_3d(fft_float* xin, fft_float* xout, int N) {
	/*const int N2 = N * N;
	const int No2 = N >> 1;
	const int NNo2 = N * No2;
	const int NNNo2 = N * NNo2;
	const int No2p1 = No2 + 1;
	const int NNo2p1 = N * No2p1;
	const int NNNo2p1 = N * NNo2p1;
	transpose_yz(xin, N);
	scramble(xin, N, N, N);
	for (int n = 0; n < N; n++) {
		fft_real_vector(xin + n * N2, N, N);
	}
	transpose_r2c(xin, xout, N);
	auto* xzero = zout;
	auto* xreal = xout + N2;
	auto* xnyqt = zout + NNNo2;
	auto* ximag = xout + NNNo2p1;
	scramble(xreal, N, N, No2);
	scramble(imag, N, N, No2);
	for (int n = 0; n < N; n++) {
		fft_complex_vector(xreal + n * N2, ximag + n * N2, N, No2);
	}
	fft_real_vector(xzero, N)
	scramble(x, 1, N, N * (No2 + 1));
	scramble(y, 1, N, N * (No2 + 1));
	fft_complex_vector(x, y, N, NNo2p1);*/
}
