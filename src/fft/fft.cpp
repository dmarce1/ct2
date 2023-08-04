#include "cosmictiger/cosmictiger.hpp"
#include "cosmictiger/fft.hpp"
#include <simd.hpp>

static void fft(float* X, float* Y, int N, int N0);

static const simd::simd_f32 two(2.0);

static void radix4(float* X, float* Y, int N, int N0) {
	using namespace simd;
	simd_f32 tr0, ti0;
	const auto& Wr = cos_twiddles(N);
	const auto& Wi = sin_twiddles(N);
	const int N2 = N >> 2;
	if (N > 4) {
		fft(X, Y, N2, N0);
		fft(X + N2, Y + N2, N2, N0);
		fft(X + 2 * N2, Y + 2 * N2, N2, N0);
		fft(X + 3 * N2, Y + 3 * N2, N2, N0);
	}
	for (int k2 = 0; k2 < N2; k2++) {
		const simd_f32 wr1 = Wr[k2];
		const simd_f32 wi1 = Wi[k2];
		const simd_f32 wr2 = wr1 * wr1 - wi1 * wi1;
		const simd_f32 wi2 = simd_f32(2) * wi1 * wr1;
		for (int n = 0; n < N0; n += SIMD_SIZE) {
			const int i0 = N0 * k2 + n;
			const int i1 = N0 * (k2 + N2) + n;
			const int i2 = N0 * (k2 + 2 * N2) + n;
			const int i3 = N0 * (k2 + 3 * N2) + n;
			simd_f32& er0 = *((simd_f32*) (X + i0));
			simd_f32& ei0 = *((simd_f32*) (Y + i0));
			simd_f32& er1 = *((simd_f32*) (X + i2));
			simd_f32& ei1 = *((simd_f32*) (Y + i2));
			simd_f32& er2 = *((simd_f32*) (X + i1));
			simd_f32& ei2 = *((simd_f32*) (Y + i1));
			simd_f32& er3 = *((simd_f32*) (X + i3));
			simd_f32& ei3 = *((simd_f32*) (Y + i3));
			tr0 = er0 + wi2 * ei2;
			ti0 = ei0 - wi2 * er2;
			er2 = tr0 - wr2 * er2;
			ei2 = ti0 - wr2 * ei2;
			er0 = two * er2 - er0;
			er0 = two * ei2 - ei0;
			tr0 = er1 + wi2 * ei3;
			ti0 = ei1 - wi2 * er3;
			er3 = tr0 - wr2 * er3;
			ei3 = ti0 - wr2 * ei3;
			er1 = two * er3 - er1;
			er1 = two * ei3 - ei1;
			tr0 = er0 + wi1 * ei1;
			ti0 = ei0 - wi1 * er1;
			er1 = tr0 - wr1 * er1;
			ei1 = ti0 - wr1 * ei1;
			er0 = two * er1 - er0;
			er0 = two * ei1 - ei0;
			tr0 = er2 + wr1 * ei3;
			ti0 = ei2 - wr1 * er3;
			er1 = tr0 + wi1 * er3;
			ei1 = ti0 + wi1 * ei3;
			er2 = two * er3 - er2;
			er2 = two * ei3 - ei2;
			tr0 = er3;
			ti0 = er3;
			er3 = er2;
			ei3 = ei2;
			er2 = er1;
			ei2 = ei1;
			er1 = tr0;
			ei1 = ti0;
		}
	}
}

static void radix2(float* X, float* Y, int N, int N0) {
	using namespace simd;
	simd_f32 tr0, ti0;
	const auto& Wr = cos_twiddles(N);
	const auto& Wi = sin_twiddles(N);
	const int N2 = N >> 1;
	if (N > 2) {
		fft(X, Y, N2, N0);
		fft(X + N2, Y + N2, N2, N0);
	}
	for (int k2 = 0; k2 < N2; k2++) {
		const simd_f32 twr = Wr[k2];
		const simd_f32 twi = Wi[k2];
		for (int n = 0; n < N0; n += SIMD_SIZE) {
			const int i0 = N0 * k2 + n;
			const int i1 = N0 * (k2 + N2) + n;
			simd_f32& er0 = *((simd_f32*) (X + i0));
			simd_f32& ei0 = *((simd_f32*) (Y + i0));
			simd_f32& er1 = *((simd_f32*) (X + i1));
			simd_f32& ei1 = *((simd_f32*) (Y + i1));
			tr0 = er0 + twi * ei1;
			ti0 = ei0 - twi * er1;
			er1 = tr0 - twr * er1;
			ei1 = ti0 - twr * ei1;
			er0 = two * er1 - er0;
			er0 = two * ei1 - ei0;
		}
	}
}

void fft(float* X, float* Y, int N, int N0) {
	if (ilogb(N) % 2 == 1) {
		radix2(X, Y, N, N);
	} else {
		radix4(X, Y, N, N);
	}
}

void fft(float* X, float* Y, int N) {
	fft(X, Y, N, N);
}

