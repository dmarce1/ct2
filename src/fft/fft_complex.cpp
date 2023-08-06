#include <cosmictiger/cosmictiger.hpp>
#include <cosmictiger/fft.hpp>

static const fft_simd two(2.0);

static void radix4_complex(fft_float* X, fft_float* Y, int N, int N0, int dir);
static void radix2_complex(fft_float* X, fft_float* Y, int N, int N0, int dir);

void fft_complex(fft_float* X, fft_float* Y, int N, int N0, int dir) {
	if ((ilogb(N) % 2 == 1) && (N <= 32)) {
		radix2_complex(X, Y, N, N0, dir);
	} else {
		radix4_complex(X, Y, N, N0, dir);
	}
}

static void radix2_complex(fft_float* X, fft_float* Y, int N, int N0, int dir) {
	using namespace simd;
	fft_simd tr0, ti0;
	const auto& Wr = cos_twiddles(N);
	const auto& Wi = (dir > 0) ? sin_twiddles(N) : inv_sin_twiddles(N);
	const int N2 = N >> 1;
	const int N0N2 = N0 * N2;
	if (N > 2) {
		fft_complex(X, Y, N2, N0, dir);
		fft_complex(X + N0N2, Y + N0N2, N2, N0, dir);
	}
	for (int k2 = 0; k2 < N2; k2++) {
		const fft_simd wr1 = Wr[k2];
		const fft_simd wi1 = Wi[k2];
		for (int n = 0; n < N0; n += FFT_SIMD_SIZE) {
			const int i0 = N0 * k2 + n;
			const int i1 = i0 + N0N2;
			fft_simd& er0 = *((fft_simd*) (X + i0));
			fft_simd& er1 = *((fft_simd*) (X + i1));
			fft_simd& ei0 = *((fft_simd*) (Y + i0));
			fft_simd& ei1 = *((fft_simd*) (Y + i1));
			tr0 = er0 + wi1 * ei1;
			ti0 = ei0 - wi1 * er1;
			er1 = tr0 - wr1 * er1;
			ei1 = ti0 - wr1 * ei1;
			er0 += er0 - er1;
			ei0 += ei0 - ei1;
		}
	}
}

static void radix4_complex(fft_float* X, fft_float* Y, int N, int N0, int dir) {
	using namespace simd;
	fft_simd tr0, ti0;
	const auto& Wr = cos_twiddles(N);
	const auto& Wi = dir > 0 ? sin_twiddles(N) : inv_sin_twiddles(N);
	const int N2 = N >> 2;
	const int N0N2 = N0 * N2;
	if (N > 4) {
		fft_complex(X, Y, N2, N0, dir);
		fft_complex(X + N0N2, Y + N0N2, N2, N0, dir);
		fft_complex(X + 2 * N0N2, Y + 2 * N0N2, N2, N0, dir);
		fft_complex(X + 3 * N0N2, Y + 3 * N0N2, N2, N0, dir);
	}
	for (int k2 = 0; k2 < N2; k2++) {
		const fft_simd wr1 = Wr[k2];
		const fft_simd wi1 = Wi[k2];
		const fft_simd wr2 = wr1 * wr1 - wi1 * wi1;
		const fft_simd wi2 = fft_simd(2) * wi1 * wr1;
		for (int n = 0; n < N0; n += FFT_SIMD_SIZE) {
			const int i0 = N0 * k2 + n;
			const int i1 = i0 + N0N2;
			const int i2 = i1 + N0N2;
			const int i3 = i2 + N0N2;
			fft_simd& er0 = *((fft_simd*) (X + i0));
			fft_simd& ei0 = *((fft_simd*) (Y + i0));
			fft_simd& er1 = *((fft_simd*) (X + i2));
			fft_simd& ei1 = *((fft_simd*) (Y + i2));
			fft_simd& er2 = *((fft_simd*) (X + i1));
			fft_simd& ei2 = *((fft_simd*) (Y + i1));
			fft_simd& er3 = *((fft_simd*) (X + i3));
			fft_simd& ei3 = *((fft_simd*) (Y + i3));
			tr0 = er0 + wi2 * ei2;
			ti0 = ei0 - wi2 * er2;
			er2 = tr0 - wr2 * er2;
			ei2 = ti0 - wr2 * ei2;
			tr0 = er1 + wi2 * ei3;
			ti0 = ei1 - wi2 * er3;
			er3 = tr0 - wr2 * er3;
			ei3 = ti0 - wr2 * ei3;
			er0 += er0 - er2;
			ei0 += ei0 - ei2;
			er1 += er1 - er3;
			ei1 += ei1 - ei3;
			tr0 = er0 + wi1 * ei1;
			ti0 = ei0 - wi1 * er1;
			er1 = tr0 - wr1 * er1;
			ei1 = ti0 - wr1 * ei1;
			tr0 = er2 + wr1 * ei3;
			ti0 = ei2 - wr1 * er3;
			er3 = tr0 + wi1 * er3;
			ei3 = ti0 + wi1 * ei3;
			er0 += er0 - er1;
			ei0 += ei0 - ei1;
			er2 += er2 - er3;
			ei2 += ei2 - ei3;
			std::swap(er2, er3);
			std::swap(ei2, ei3);
		}
	}
}

