#include "cosmictiger/fft.hpp"
#include <hpx/async.hpp>

void fft_3d(float* xin, float* xout, int N) {
	const int N2 = N * N;
	const int No2 = N >> 1;
	transpose_yz(xin, N);
	scramble(xin, N, N, N);
	for (int n = 0; n < N; n++) {
		fft_real(xin + n * N2, N, N);
	}
	transpose_r2c(xin, xout, N);
	auto* x = xout;
	auto* y = xout + N2 * No2;
	scramble(x, N, N, No2);
	scramble(y, N, N, No2);
	for (int n = 0; n < N; n++) {
		fft_complex(x + n * N2, y + n * N2, N, No2);
	}
	scramble(x, 1, N, N * No2);
	scramble(y, 1, N, N * No2);
	fft_complex(x, y, N, N * No2);
}
