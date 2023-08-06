#include <cosmictiger/cosmictiger.hpp>
#include <cosmictiger/fft.hpp>

static const fft_simd two(2.0);

static void radix4_real(fft_float* X, int N, int N0, int dir);
static void radix2_real(fft_float* X, int N, int N0, int dir);

void fft_real(fft_float* X, int N, int N0, int dir) {
	if ((ilogb(N) % 2 == 1) && (N <= 32)) {
		radix2_real(X, N, N0, dir);
	} else {
		radix4_real(X, N, N0, dir);
	}
}

static void radix2_real(fft_float* X, int N, int N0, int dir) {
	using namespace simd;
	const int N2 = N >> 1;
	const int N2o2 = N >> 2;
	const int N0N2 = N0 * N2;
	if (N > 2) {
		fft_real(X, N2, N0, dir);
		fft_real(X + N0N2, N2, N0, dir);
	}
	{
		fft_simd tr0, tr1;
		for (int n = 0; n < N0; n += FFT_SIMD_SIZE) {
			const int& i0 = n;
			const int i1 = i0 + N0N2;
			fft_simd& er0 = *((fft_simd*) (X + i0));
			fft_simd& er1 = *((fft_simd*) (X + i1));
			tr0 = er0 + er1;
			tr1 = er0 - er1;
			er0 = tr0;
			er1 = tr1;
		}
	}
	if (N2o2) {
		const int _3N0N2o2 = N0N2 + (N0N2 >> 1);
		for (int n = 0; n < N0; n += FFT_SIMD_SIZE) {
			fft_simd& er1 = *((fft_simd*) (X + _3N0N2o2 + n));
			er1 = -er1;
		}
	}
	if (N2o2 > 1) {
		fft_simd tr0, ti0, tr1;
		const auto& Wr = cos_twiddles(N);
		const auto& Wi = dir > 0 ? sin_twiddles(N) : inv_sin_twiddles(N);
		for (int k2 = 1; k2 < N2o2; k2++) {
			const fft_simd wr1 = Wr[k2];
			const fft_simd wi1 = Wi[k2];
			for (int n = 0; n < N0; n += FFT_SIMD_SIZE) {
				const int i0 = N0 * k2 + n;
				const int i1 = i0 + N0N2;
				const int j0 = N0 * (N2 - k2) + n;
				const int j1 = j0 + N0N2;
				fft_simd& er0 = *((fft_simd*) (X + i0));
				fft_simd& ei0 = *((fft_simd*) (X + j0));
				fft_simd& er1 = *((fft_simd*) (X + i1));
				fft_simd& ei1 = *((fft_simd*) (X + j1));
				tr0 = er0 + wi1 * ei1;
				ti0 = ei0 - wi1 * er1;
				er1 = tr0 - wr1 * er1;
				ei1 = wr1 * ei1 - ti0;
				er0 += er0 - er1;
				ei0 += ei0 + ei1;
				ti0 = ei0;
				ei0 = er1;
				er1 = ei1;
				ei1 = ti0;
			}
		}
	}
}

static void radix4_real(fft_float* X, int N, int N0, int dir) {
	using namespace simd;
	const int N2 = N >> 2;
	const int N2o2 = N >> 3;
	const int N0N2 = N0 * N2;
	if (N > 4) {
		fft_real(X, N2, N0, dir);
		fft_real(X + N0N2, N2, N0, dir);
		fft_real(X + 2 * N0N2, N2, N0, dir);
		fft_real(X + 3 * N0N2, N2, N0, dir);
	}
	{
		fft_simd tr0, tr1, tr2, tr3;
		for (int n = 0; n < N0; n += FFT_SIMD_SIZE) {
			const int& i0 = n;
			const int i1 = i0 + N0N2;
			const int i2 = i1 + N0N2;
			const int i3 = i2 + N0N2;
			fft_simd& er0 = *((fft_simd*) (X + i0));
			fft_simd& er1 = *((fft_simd*) (X + i2));
			fft_simd& er2 = *((fft_simd*) (X + i1));
			fft_simd& er3 = *((fft_simd*) (X + i3));
			tr0 = er0 + er2;
			tr2 = er0 - er2;
			tr1 = er1 + er3;
			tr3 = er3 - er1;
			er0 = tr0 + tr1;
			er2 = tr0 - tr1;
			er1 = tr2;
			er3 = tr3;
			std::swap(er1, er2);
		}
	}
	if (N2o2) {
		fft_simd tr0, tr1, tr2, tr3;
		const int N0N2o2 = N0N2 >> 1;
		for (int n = 0; n < N0; n += FFT_SIMD_SIZE) {
			const int i0 = N0N2o2 + n;
			const int i1 = i0 + N0N2;
			const int i2 = i1 + N0N2;
			const int i3 = i2 + N0N2;
			fft_simd& er0 = *((fft_simd*) (X + i0));
			fft_simd& er1 = *((fft_simd*) (X + i2));
			fft_simd& er2 = *((fft_simd*) (X + i1));
			fft_simd& er3 = *((fft_simd*) (X + i3));
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
	if (N2o2 > 1) {
		fft_simd tr0, ti0;
		const auto& Wr = cos_twiddles(N);
		const auto& Wi = (dir > 0) ? sin_twiddles(N) : inv_sin_twiddles(N);
		for (int k2 = 1; k2 < N2o2; k2++) {
			const fft_simd wr1 = Wr[k2];
			const fft_simd wi1 = Wi[k2];
			const fft_simd wr2 = wr1 * wr1 - wi1 * wi1;
			const fft_simd wi2 = fft_simd(2) * wi1 * wr1;
			for (int n = 0; n < N0; n += FFT_SIMD_SIZE) {
				const int i0 = N0 * k2 + n;
				const int i1 = i0 + N0N2;
				const int i2 = i1 + N0N2;
				const int i3 = i2 + N0N2;
				const int j0 = N0 * (N2 - k2) + n;
				const int j1 = j0 + N0N2;
				const int j2 = j1 + N0N2;
				const int j3 = j2 + N0N2;
				fft_simd& er0 = *((fft_simd*) (X + i0));
				fft_simd& ei0 = *((fft_simd*) (X + j0));
				fft_simd& er1 = *((fft_simd*) (X + i2));
				fft_simd& ei1 = *((fft_simd*) (X + j2));
				fft_simd& er2 = *((fft_simd*) (X + i1));
				fft_simd& ei2 = *((fft_simd*) (X + j1));
				fft_simd& er3 = *((fft_simd*) (X + i3));
				fft_simd& ei3 = *((fft_simd*) (X + j3));
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
				ei1 = wr1 * ei1 - ti0;
				tr0 = er2 + wr1 * ei3;
				ti0 = ei2 - wr1 * er3;
				er3 = tr0 + wi1 * er3;
				ei3 = ti0 + wi1 * ei3;
				er0 += er0 - er1;
				ei0 += ei0 + ei1;
				er2 += er2 - er3;
				ei2 = ei3 - two * ei2;
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
}

