#include <cosmictiger/cosmictiger.hpp>
#include <cosmictiger/fft.hpp>


static const simd::simd_f32 two(2.0);

static void radix4_real(float* X, int N, int N0);
static void radix2_real(float* X, int N, int N0);

void fft_real(float* X, int N, int N0) {
	if (ilogb(N) % 2 == 1) {
		radix2_real(X, N, N0);
	} else {
		radix4_real(X, N, N0);
	}
}

static void radix4_real(float* X, int N, int N0) {
	using namespace simd;
	simd_f32 tr0, ti0, tr1, tr2, tr3;
	const auto& Wr = cos_twiddles(N);
	const auto& Wi = sin_twiddles(N);
	const int N2 = N >> 2;
	const int N2o2 = N >> 3;
	const int N0N2 = N0 * N2;
	if (N > 4) {
		fft_real(X, N2, N0);
		fft_real(X + N0N2, N2, N0);
		fft_real(X + 2 * N0N2, N2, N0);
		fft_real(X + 3 * N0N2, N2, N0);
	}
	for (int n = 0; n < N0; n += SIMD_SIZE) {
		const int& i0 = n;
		const int i1 = i0 + N0N2;
		const int i2 = i1 + N0N2;
		const int i3 = i2 + N0N2;
		simd_f32& er0 = *((simd_f32*) (X + i0));
		simd_f32& er1 = *((simd_f32*) (X + i2));
		simd_f32& er2 = *((simd_f32*) (X + i1));
		simd_f32& er3 = *((simd_f32*) (X + i3));
		tr0 = er0 + er2;
		tr2 = er0 - er2;
		tr1 = er1 + er3;
		tr3 = er1 - er3;
		er0 = tr0 + tr1;
		er2 = tr0 - tr1;
		er1 = tr2;
		er3 = tr3;
		std::swap(er1, er2);
	}
	if (N2o2) {
		for (int n = 0; n < N0; n += SIMD_SIZE) {
			const int& i0 = n;
			const int i1 = i0 + N0N2;
			const int i2 = i1 + N0N2;
			const int i3 = i2 + N0N2;
			simd_f32& er0 = *((simd_f32*) (X + i0));
			simd_f32& er1 = *((simd_f32*) (X + i2));
			simd_f32& er2 = *((simd_f32*) (X + i1));
			simd_f32& er3 = *((simd_f32*) (X + i3));
			tr0 = er1 + er3;
			tr2 = er1 - er3;
			tr1 = M_SQRT1_2 * tr2;
			tr3 = M_SQRT1_2 * tr0;
			tr0 = er0;
			tr2 = er2;
			er0 = tr0 + tr1;
			er3 = -tr2 - tr3;
			er1 = tr0 - tr1;
			er2 = tr2 - tr3;
			std::swap(er1, er2);
		}
	}
	for (int k2 = 1; k2 < N2o2; k2++) {
		const simd_f32 wr1 = Wr[k2];
		const simd_f32 wi1 = Wi[k2];
		const simd_f32 wr2 = wr1 * wr1 - wi1 * wi1;
		const simd_f32 wi2 = simd_f32(2) * wi1 * wr1;
		for (int n = 0; n < N0; n += SIMD_SIZE) {
			const int i0 = N0 * k2 + n;
			const int i1 = i0 + N0N2;
			const int i2 = i1 + N0N2;
			const int i3 = i2 + N0N2;
			const int j0 = N0 * (N2 - k2) + n;
			const int j1 = j0 + N0N2;
			const int j2 = j1 + N0N2;
			const int j3 = j2 + N0N2;
			simd_f32& er0 = *((simd_f32*) (X + i0));
			simd_f32& ei0 = *((simd_f32*) (X + j0));
			simd_f32& er1 = *((simd_f32*) (X + i2));
			simd_f32& ei1 = *((simd_f32*) (X + j2));
			simd_f32& er2 = *((simd_f32*) (X + i1));
			simd_f32& ei2 = *((simd_f32*) (X + j1));
			simd_f32& er3 = *((simd_f32*) (X + i3));
			simd_f32& ei3 = *((simd_f32*) (X + j3));
			tr0 = er0 + wi2 * ei2;
			ti0 = ei0 - wi2 * er2;
			er2 = tr0 - wr2 * er2;
			ei2 = ti0 - wr2 * ei2;
			er0 = two * er2 - er0;
			ei0 = two * ei2 - ei0;
			tr0 = er1 + wi2 * ei3;
			ti0 = ei1 - wi2 * er3;
			er3 = tr0 - wr2 * er3;
			ei3 = ti0 - wr2 * ei3;
			er1 = two * er3 - er1;
			ei1 = two * ei3 - ei1;
			tr0 = er0 + wi1 * ei1;
			ti0 = ei0 - wi1 * er1;
			er1 = tr0 - wr1 * er1;
			ei1 = ti0 - wr1 * ei1;
			er0 = two * er1 - er0;
			ei0 = two * ei1 - ei0;
			tr0 = er2 + wr1 * ei3;
			ti0 = ei2 - wr1 * er3;
			er3 = tr0 + wi1 * er3;
			ei3 = -ti0 - wi1 * ei3;
			er2 = two * er3 - er2;
			ei2 = two * ei3 + ei2;
			ti0 = ei0;
			ei0 = er2;
			er2 = er3;
			er3 = ei2;
			ei2 = er1;
			er1 = ei1;
			ei1 = ei3;
			ei3 = ti0;
		}
	}
}

static void radix2_real(float* X, int N, int N0) {
	using namespace simd;
	simd_f32 tr0, ti0, tr1;
	const auto& Wr = cos_twiddles(N);
	const auto& Wi = sin_twiddles(N);
	const int N2 = N >> 1;
	const int N2o2 = N >> 2;
	const int N0N2 = N0 * N2;
	if (N > 2) {
		fft_real(X, N2, N0);
		fft_real(X + N0N2, N2, N0);
	}
	for (int n = 0; n < N0; n += SIMD_SIZE) {
		const int& i0 = n;
		const int i1 = i0 + N0N2;
		simd_f32& er0 = *((simd_f32*) (X + i0));
		simd_f32& er1 = *((simd_f32*) (X + i1));
		tr0 = er0 + er1;
		tr1 = er0 - er1;
		er0 = tr0;
		er1 = tr1;
	}
	if (N2o2) {
		for (int n = 0; n < N0; n += SIMD_SIZE) {
			const int& i0 = n;
			const int i1 = i0 + N0N2;
			simd_f32& er1 = *((simd_f32*) (X + i1));
			er1 = -er1;
		}
	}
	for (int k2 = 1; k2 < N2o2; k2++) {
		const simd_f32 wr1 = Wr[k2];
		const simd_f32 wi1 = Wi[k2];
		for (int n = 0; n < N0; n += SIMD_SIZE) {
			const int i0 = N0 * k2 + n;
			const int i1 = i0 + N0N2;
			const int j0 = N0 * (N2 - k2) + n;
			const int j1 = j0 + N0N2;
			simd_f32& er0 = *((simd_f32*) (X + i0));
			simd_f32& ei0 = *((simd_f32*) (X + j0));
			simd_f32& er1 = *((simd_f32*) (X + i1));
			simd_f32& ei1 = *((simd_f32*) (X + j1));
			tr0 = er0 + wi1 * ei1;
			ti0 = ei0 - wi1 * er1;
			er1 = tr0 - wr1 * er1;
			ei1 = ti0 - wr1 * ei1;
			er0 = two * er1 - er0;
			ei0 = two * ei1 - ei0;
			std::swap(ei0, er1);
		}
	}
}

