#include "cosmictiger/fft.hpp"
#include <hpx/async.hpp>

using namespace sfmm;

void fft_r2c_3d(fft_float* Xin, complex<fft_float>* Xout, int N) {
	auto* X = Xin;
	const auto& Wr = cos_twiddles(N);
	const auto& Wi = sin_twiddles(N);
	fft_float* Y = X + N * N * N / 2;
	scramble(X, N);
	for (int k = 0; k < N / 2; k++) {
		fft_complex(X + k * N * N, Y + k * N * N, N, N);
	}
	fft_complex(X, Y, N / 2, N * N);
	transpose_yz(X, N);
	// z x y
	for (int k = 0; k < N / 2; k++) {
		fft_complex(X + k * N * N, Y + k * N * N, N, N);
	}
	const complex<fft_float> J(0, 1);
	for (int k = 0; k <= N / 2; k++) {
		complex<fft_float> W(Wr[k], Wi[k]);
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				const int n0 = j + N * (i + N * (k % (N / 2)));
				const int n1 = ((N - j) % N) + N * (((N - i) % N) + N * ((N / 2 - k) % (N / 2)));
				const complex<fft_float> Z0(X[n0], Y[n0]);
				const complex<fft_float> Z1(X[n1], Y[n1]);
				const auto Xe = (Z0 + Z1) * 0.5;
				const auto Xo = (Z1 - Z0) * J * 0.5;
				Xout[k + (N / 2 + 1) * (j + N * i)] = Xe + Xo * W;
			}
		}
	}
}
