#include "cosmictiger/fft.hpp"
#include <hpx/async.hpp>

void fft_c2r_3d(fft_float* xin, fft_float* xout, int N) {
	const int N2 = N * N;
	const int No2 = N >> 1;
	auto* x = xin;
	auto* y = xin + N2 * No2;
	scramble(x, 1, N, N * No2);
	scramble(y, 1, N, N * No2);
	fft_complex(x, y, N, N * No2, -1);
	scramble(x, N, N, No2);
	scramble(y, N, N, No2);
	for (int n = 0; n < N; n++) {
		fft_complex(x + n * N2, y + n * N2, N, No2, -1);
	}
	transpose_c2r(xin, xout, N);
	for (int n = 0; n < N; n++) {
		for (int l = 1; l < N - l; l++) {
			for (int m = 0; m < N; m++) {
				const fft_float tmp = xout[l];
				const int i = (n * N + m) * N;
				xout[l + i] += xout[N - l + i];
				xout[N - l + i] = tmp - xout[N - l + i];
			}
		}
		fft_real(xout + N2 * n, N, N, -1);
		for (int l = 1; l < N - l; l++) {
			for (int m = 0; m < N; m++) {
				const fft_float tmp = xout[l];
				const int i = (n * N + m) * N;
				xout[l + i] += xout[N - l + i];
				xout[N - l + i] = tmp - xout[N - l + i];
			}
		}
	}
	transpose_yz(xin, N);
	scramble(xin, N, N, N);
}
