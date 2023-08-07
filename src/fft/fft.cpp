#include <cosmictiger/cosmictiger.hpp>
#include <cosmictiger/fft.hpp>
#include <sfmm.hpp>

using namespace simd;
using namespace sfmm;

static const fft_simd two(2.0);

static void radix4(complex<fft_simd>*, int, int);
static void radix2(complex<fft_simd>*, int, int);
static void scramble(complex<fft_simd>*, int);
static int reverse_increment(int, int);

static int reverse_increment(int i, int N) {
	const int Nm1 = N - 1;
	int j, k;
	j = i;
	j = ~j;
	j = j & Nm1;
	k = 31 - __builtin_clz(j);
	j = 1;
	j = j << k;
	j--;
	k = i;
	k = j & k;
	j++;
	k = j | k;
	i = k;
	return k;
}

static void scramble(complex<fft_simd>* X, int N) {
	int l = 0;
	for (int j = 0; j < N; j++) {
		if (j > l) {
			std::swap(X[j], X[l]);
		}
		l = reverse_increment(l, N);
	}
}

void fft_1d_simd(complex<fft_simd>* X, int N, int dir, bool scrmbl) {
	if (scrmbl) {
		scramble(X, N);
	}
	if ((ilogb(N) % 2 == 1) && (N <= 32)) {
		radix2(X, N, dir);
	} else {
		radix4(X, N, dir);
	}
}

static void radix2(complex<fft_simd>* X, int N, int dir) {
	const auto& Wr = cos_twiddles(N);
	const auto& Wi = (dir > 0) ? sin_twiddles(N) : inv_sin_twiddles(N);
	const int N2 = N >> 1;
	if (N > 2) {
		fft_1d_simd(X, N2, dir, false);
		fft_1d_simd(X + N2, N2, dir, false);
	}
	for (int k2 = 0; k2 < N2; k2++) {
		const complex<fft_simd> W1(Wr[k2], Wi[k2]);
		complex<fft_simd>& U0 = X[k2];
		complex<fft_simd>& U1 = X[N2 + k2];
		U1 = U0 - W1 * U1;
		U0 += U0 - U1;
	}
}

static void radix4(complex<fft_simd>* X, int N, int dir) {
	using namespace simd;
	const auto& Wr = cos_twiddles(N);
	const auto& Wi = dir > 0 ? sin_twiddles(N) : inv_sin_twiddles(N);
	const int N2 = N >> 2;
	if (N > 4) {
		fft_1d_simd(X, N2, dir, false);
		fft_1d_simd(X + N2, N2, dir, false);
		fft_1d_simd(X + 2 * N2, N2, dir, false);
		fft_1d_simd(X + 3 * N2, N2, dir, false);
	}
	for (int k2 = 0; k2 < N2; k2++) {
		const complex<fft_simd> W1(Wr[k2], Wi[k2]);
		const complex<fft_simd> W2 = W1 * W1;
		complex<fft_simd>& U0 = X[k2];
		complex<fft_simd>& U1 = X[k2 + 2 * N2];
		complex<fft_simd>& U2 = X[k2 + N2];
		complex<fft_simd>& U3 = X[k2 + 3 * N2];
		U2 = U0 - W2 * U2;
		U3 = U1 - W2 * U3;
		U0 += U0 - U2;
		U1 += U1 - U3;
		U3.imag() = -U3.imag();
		std::swap(U3.imag(), U3.real());
		U1 = U0 - W1 * U1;
		U3 = U2 - W1 * U3;
		U0 += U0 - U1;
		U2 += U2 - U3;
		std::swap(U2, U3);
	}
}

void fft_2d(complex<fft_float>* X, int N, int dir) {
	std::vector < complex < fft_simd >> Y(N);
	for (int i = 0; i < N; i += FFT_SIMD_SIZE) {
		for (int j = 0; j < FFT_SIMD_SIZE; j++) {
			for (int k = 0; k < N; k++) {
				const int q = N * (i + j) + k;
				Y[k].real()[j] = X[q].real();
				Y[k].imag()[j] = X[q].imag();
			}
		}
		fft_1d_simd(Y.data(), N, dir);
		for (int j = 0; j < FFT_SIMD_SIZE; j++) {
			for (int k = 0; k < N; k++) {
				const int q = N * (i + j) + k;
				X[q].real() = Y[k].real()[j];
				X[q].imag() = Y[k].imag()[j];
			}
		}
	}
	for (int i = 0; i < N; i += FFT_SIMD_SIZE) {
		for (int j = 0; j < FFT_SIMD_SIZE; j++) {
			for (int k = 0; k < N; k++) {
				const int q = i + j + N * k;
				Y[k].real()[j] = X[q].real();
				Y[k].imag()[j] = X[q].imag();
			}
		}
		fft_1d_simd(Y.data(), N, dir);
		for (int j = 0; j < FFT_SIMD_SIZE; j++) {
			for (int k = 0; k < N; k++) {
				const int q = i + j + N * k;
				X[q].real() = Y[k].real()[j];
				X[q].imag() = Y[k].imag()[j];
			}
		}
	}
}

void fft_3d(complex<fft_float>* X, int N, int dir) {
	std::vector < complex < fft_simd >> Y(N * N);
	for (int i = 0; i < N; i++) {
		fft_2d(X + i * N * N, N, dir);
	}
	for (int i = 0; i < N; i += FFT_SIMD_SIZE) {
		for (int k = 0; k < N; k++) {
			for (int j = 0; j < FFT_SIMD_SIZE; j++) {
				for (int l = 0; l < N; l++) {
					const int q = i + j + N * (k + N * l);
					Y[l].real()[j] = X[q].real();
					Y[l].imag()[j] = X[q].imag();
				}
			}
			fft_1d_simd(Y.data(), N, dir);
			for (int j = 0; j < FFT_SIMD_SIZE; j++) {
				for (int l = 0; l < N; l++) {
					const int q = i + j + N * (k + N * l);
					X[q].real() = Y[l].real()[j];
					X[q].imag() = Y[l].imag()[j];
				}
			}
		}
	}
}
