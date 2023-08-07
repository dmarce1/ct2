#include <cosmictiger/cosmictiger.hpp>
#include <cosmictiger/fft.hpp>

static const fft_simd two(2.0);

static void radix4_complex(complex<fft_simd>* X, int N, int dir);
static void radix2_complex(complex<fft_simd>* X, int N, int dir);

void fft_complex(complex<fft_simd>* X, int N, int dir) {
	if ((ilogb(N) % 2 == 1) && (N <= 32)) {
		radix2_complex(X, N, dir);
	} else {
		radix4_complex(X, N, dir);
	}
}

static void radix2_complex(complex<fft_simd>* X, int N, int dir) {
	using namespace simd;
	const auto& Wr = cos_twiddles(N);
	const auto& Wi = (dir > 0) ? sin_twiddles(N) : inv_sin_twiddles(N);
	const int N2 = N >> 1;
	if (N > 2) {
		fft_complex(X, N2, dir);
		fft_complex(X + N2, N2, dir);
	}
	for (int k2 = 0; k2 < N2; k2++) {
		const complex<fft_simd> w1(Wr[k2], Wi[k2]);
		complex<fft_simd>& e0 = X[k2];
		complex<fft_simd>& e1 = X[k2 + N2];
		e1 = e0 - w1 * e1;
		e0 += e0 - e1;
	}
}

static void radix4_complex(complex<fft_simd>* X, int N, int dir) {
	using namespace simd;
	const auto& Wr = cos_twiddles(N);
	const auto& Wi = dir > 0 ? sin_twiddles(N) : inv_sin_twiddles(N);
	const int N2 = N >> 2;
	if (N > 4) {
		fft_complex(X, N2, dir);
		fft_complex(X + N2, N2, dir);
		fft_complex(X + 2 * N2, N2, dir);
		fft_complex(X + 3 * N2, N2, dir);
	}
	for (int k2 = 0; k2 < N2; k2++) {
		const complex<fft_simd> w1(Wr[k2], Wi[k2]);
		const complex<fft_simd> w2 = w1 * w1;
		complex<fft_simd>& e0 = X[k2];
		complex<fft_simd>& e1 = X[k2 + 2 * N2];
		complex<fft_simd>& e2 = X[k2 + N2];
		complex<fft_simd>& e3 = X[k2 + 3 * N2];
		e2 = e0 - w2 * e2;
		e3 = e1 - w2 * e3;
		e0 += e0 - e2;
		e1 += e1 - e3;
		e3.imag() = -e3.imag();
		std::swap(e3.real(), e3.imag());
		e1 = e0 - w1 * e1;
		e3 = e2 - w1 * e3;
		e0 += e0 - e1;
		e2 += e2 - e3;
		std::swap(e2, e3);
	}
}

