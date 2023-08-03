#include "cosmictiger/fft.hpp"
#include <hpx/async.hpp>

void fft_3d(float* X, float* Y, int N) {
	std::vector<hpx::future<void>> futs(N);
	const int N2 = N * N;
	transpose_yz(X, N);
	transpose_yz(Y, N);
	scramble(X, N);
	scramble(Y, N);
	for (int i = 0; i < N; i++) {
		const int j = i * N2;
		futs[i] = hpx::async(fft, X + j, Y + j, N);
	}
	hpx::wait_all(futs.begin(), futs.end());
	transpose_yz(X, N);
	transpose_yz(Y, N);
	scramble(X, N);
	scramble(Y, N);
	for (int i = 0; i < N; i++) {
		const int j = i * N2;
		futs[i] = hpx::async(fft, X + j, Y + j, N);
	}
	hpx::wait_all(futs.begin(), futs.end());
	transpose_xy(X, N);
	transpose_xy(Y, N);
	scramble(X, N);
	scramble(Y, N);
	for (int i = 0; i < N; i++) {
		const int j = i * N2;
		futs[i] = hpx::async(fft, X + j, Y + j, N);
	}
	hpx::wait_all(futs.begin(), futs.end());
	transpose_xy(X, N);
	transpose_xy(Y, N);
}
