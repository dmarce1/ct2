#include "cosmictiger/cosmictiger.hpp"
#include "cosmictiger/fft.hpp"
#include <simd.hpp>

static void radix4(float* X, float* Y, int N) {
	using namespace simd;
	static const simd_f32 two(2.0);
	simd_f32 tr, ti;
	const auto& Wr = cos_twiddles(N);
	const auto& Wi = sin_twiddles(N);
	const int N2 = N >> 2;
	if (N > 4) {
		fft(X, Y, N2);
		fft(X + N2, Y + N2, N2);
		fft(X + 2 * N2, Y + 2 * N2, N2);
		fft(X + 3 * N2, Y + 3 * N2, N2);
	}
	for (int k2 = 0; k2 < N2; k2++) {
		const simd_f32 C1 = Wr[k2];
		const simd_f32 S1 = Wi[k2];
		const simd_f32 C2 = C1 * C1 - S1 * S1;
		const simd_f32 S2 = simd_f32(2.0) * S1 * C1;
		for (int n = 0; n < N; n += SIMD_SIZE) {
			const int i0 = N * k2 + n;
			const int i1 = N * (k2 + N2) + n;
			const int i2 = N * (k2 + 2 * N2) + n;
			const int i3 = N * (k2 + 3 * N2) + n;
			simd_f32& er0 = *((simd_f32*) (X + i0));
			simd_f32& ei0 = *((simd_f32*) (Y + i0));
			simd_f32& er1 = *((simd_f32*) (X + i2));
			simd_f32& ei1 = *((simd_f32*) (Y + i2));
			simd_f32& er2 = *((simd_f32*) (X + i1));
			simd_f32& ei2 = *((simd_f32*) (Y + i1));
			simd_f32& er3 = *((simd_f32*) (X + i3));
			simd_f32& ei3 = *((simd_f32*) (Y + i3));
			tr = er0 + S2 * ei2;
			ti = ei0 - S2 * er2;
			er2 = tr - C2 * er2;
			ei2 = ti - C2 * ei2;
			er0 = two * er2 - er0;
			er0 = two * ei2 - ei0;
			tr = er1 + S2 * ei3;
			ti = ei1 - S2 * er3;
			er3 = tr - C2 * er3;
			ei3 = ti - C2 * ei3;
			er1 = two * er3 - er1;
			er1 = two * ei3 - ei1;
			tr = er0 + S1 * ei1;
			ti = ei0 - S1 * er1;
			er1 = tr - C2 * er1;
			ei1 = ti - C2 * ei1;
			er0 = two * er1 - er0;
			er0 = two * ei1 - ei0;
			tr = er2 + C1 * ei3;
			ti = ei2 - C1 * er3;
			er1 = tr + S2 * er3;
			ei1 = ti + S2 * ei3;
			er2 = two * er3 - er2;
			er2 = two * ei3 - ei2;
			tr = er3;
			ti = er3;
			er3 = er2;
			ei3 = ei2;
			er2 = er1;
			ei2 = ei1;
			er1 = tr;
			ei1 = ti;
		}
	}
}

static void radix2(float* X, float* Y, int N) {
	using namespace simd;
	static const simd_f32 two(2.0);
	simd_f32 tr, ti;
	const auto& Wr = cos_twiddles(N);
	const auto& Wi = sin_twiddles(N);
	const int N2 = N >> 1;
	if (N > 2) {
		fft(X, Y, N2);
		fft(X + N2, Y + N2, N2);
	}
	for (int k2 = 0; k2 < N2; k2++) {
		const simd_f32 C = Wr[k2];
		const simd_f32 S = Wi[k2];
		for (int n = 0; n < N; n += SIMD_SIZE) {
			const int i0 = N * k2 + n;
			const int i1 = N * (k2 + N2) + n;
			simd_f32& er0 = *((simd_f32*) (X + i0));
			simd_f32& ei0 = *((simd_f32*) (Y + i0));
			simd_f32& er1 = *((simd_f32*) (X + i1));
			simd_f32& ei1 = *((simd_f32*) (Y + i1));
			tr = er0 + S * ei1;
			ti = ei0 - S * er1;
			er1 = tr - C * er1;
			ei1 = ti - C * ei1;
			er0 = two * er1 - er0;
			er0 = two * ei1 - ei0;
		}
	}
}

void fft(float* X, float* Y, int N) {
	if (ilogb(N) % 2 == 1) {
		radix2(X, Y, N);
	} else {
		radix4(X, Y, N);
	}
}

