#include <cosmictiger/cosmictiger.hpp>
#include <sfmm.hpp>

#define NT 4


static void transpose_yz(float* X, int N, int x, int y, int z, int dx, int depth) {
	if (dx == NT) {
		for (int ix = x; ix < x + NT; ix++) {
			for (int iy = y; iy < y + NT; iy++) {
				for (int iz = z; iz < z + NT; iz++) {
					const int i = N * (N * ix + iy) + iz;
					const int j = N * (N * ix + iz) + iy;
					std::swap(X[i], X[j]);
				}
			}
		}
	} else {
		constexpr int z00 = 0;
		const int do2 = dx >> 1;
		const int dp1 = depth + 1;
		transpose_yz(X, N, x + z00, y + z00, z + z00, do2, dp1);
		transpose_yz(X, N, x + do2, y + z00, z + z00, do2, dp1);
		transpose_yz(X, N, x + z00, y + do2, z + z00, do2, dp1);
		transpose_yz(X, N, x + do2, y + do2, z + z00, do2, dp1);
		transpose_yz(X, N, x + z00, y + z00, z + do2, do2, dp1);
		transpose_yz(X, N, x + do2, y + z00, z + do2, do2, dp1);
		transpose_yz(X, N, x + z00, y + do2, z + do2, do2, dp1);
		transpose_yz(X, N, x + do2, y + do2, z + do2, do2, dp1);
	}
}

static void transpose_r2c(float* xin, float* xout, int N, int x, int y, int z, int dx, int depth) {
	const int No2 = N / 2;
	if (z > No2) {
		return;
	}
	const int N3o2 = N * N * N / 2;
	if (dx == NT) {
		for (int ix = x; ix < x + NT; ix++) {
			for (int iy = 0; iy < std::min(y + NT, No2 + 1); iy++) {
				if (iy == 0 || iy == No2) {
					for (int iz = y; iz < z + NT; iz++) {
						const int ir = N * (N * ix + iy) + iz;
						const int jr = N * (N * ix + iz) + iy;
						const int ji = N * (N * ix + iz) + iy + N3o2;
						xout[jr] = xin[ir];
						xout[ji] = 0.0;
					}
				} else {
					for (int iz = y; iz < z + NT; iz++) {
						const int ir = N * (N * ix + iy) + iz;
						const int ii = N * (N * ix + (N - iy)) + iz;
						const int jr = N * (N * ix + iz) + iy;
						const int ji = N * (N * ix + iz) + iy + N3o2;
						xout[jr] = xin[ir];
						xout[ji] = xin[ii];
					}
				}
			}
		}
	} else {
		constexpr int z00 = 0;
		const int do2 = dx >> 1;
		const int dp1 = depth + 1;
		transpose_r2c(xin, xout, N, x + z00, y + z00, z + z00, do2, dp1);
		transpose_r2c(xin, xout, N, x + do2, y + z00, z + z00, do2, dp1);
		transpose_r2c(xin, xout, N, x + z00, y + do2, z + z00, do2, dp1);
		transpose_r2c(xin, xout, N, x + do2, y + do2, z + z00, do2, dp1);
		transpose_r2c(xin, xout, N, x + z00, y + z00, z + do2, do2, dp1);
		transpose_r2c(xin, xout, N, x + do2, y + z00, z + do2, do2, dp1);
		transpose_r2c(xin, xout, N, x + z00, y + do2, z + do2, do2, dp1);
		transpose_r2c(xin, xout, N, x + do2, y + do2, z + do2, do2, dp1);
	}
}

void transpose_yz(float* X, int N) {
	transpose_yz(X, N, 0, 0, 0, N, 0);
}

void transpose_r2c(float* xin, float* xout, int N) {
	transpose_r2c(xin, (float*) xout, N, 0, 0, 0, N, 0);
}

