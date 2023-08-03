#include <algorithm>
#include <hpx/async.hpp>

#define N0 4
#define MAXPARDEPTH 3
#define NCHILD 8

static void transpose_yz(float* X, int N, int x, int y, int z, int dx, int depth) {
	if (dx == N0) {
		for (int ix = x; ix < x + N0; ix++) {
			for (int iy = y; iy < y + N0; iy++) {
				for (int iz = z; iz < z + N0; iz++) {
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
		if (depth > MAXPARDEPTH) {
			transpose_yz(X, N, x + z00, y + z00, z + z00, do2, dp1);
			transpose_yz(X, N, x + do2, y + z00, z + z00, do2, dp1);
			transpose_yz(X, N, x + z00, y + do2, z + z00, do2, dp1);
			transpose_yz(X, N, x + do2, y + do2, z + z00, do2, dp1);
			transpose_yz(X, N, x + z00, y + z00, z + do2, do2, dp1);
			transpose_yz(X, N, x + do2, y + z00, z + do2, do2, dp1);
			transpose_yz(X, N, x + z00, y + do2, z + do2, do2, dp1);
			transpose_yz(X, N, x + do2, y + do2, z + do2, do2, dp1);
		} else {
			std::array<hpx::future<void>, NCHILD> futs;
			futs[0] = hpx::async(transpose_yz, X, N, x + z00, y + z00, z + z00, do2, dp1);
			futs[1] = hpx::async(transpose_yz, X, N, x + do2, y + z00, z + z00, do2, dp1);
			futs[2] = hpx::async(transpose_yz, X, N, x + z00, y + do2, z + z00, do2, dp1);
			futs[3] = hpx::async(transpose_yz, X, N, x + do2, y + do2, z + z00, do2, dp1);
			futs[4] = hpx::async(transpose_yz, X, N, x + z00, y + z00, z + do2, do2, dp1);
			futs[5] = hpx::async(transpose_yz, X, N, x + do2, y + z00, z + do2, do2, dp1);
			futs[6] = hpx::async(transpose_yz, X, N, x + z00, y + do2, z + do2, do2, dp1);
			futs[7] = hpx::async(transpose_yz, X, N, x + do2, y + do2, z + do2, do2, dp1);
			hpx::wait_all(futs.begin(), futs.end());
		}
	}
}

static void transpose_xy(float* X, int N, int x, int y, int dx, int depth) {
	static thread_local std::vector<float> tmp;
	if (dx == N0) {
		tmp.resize(N);
		const int copy_sz = sizeof(float) * N;
		for (int ix = x; ix < x + N0; ix++) {
			for (int iy = y; iy < y + N0; iy++) {
				const int i = N * (N * ix + iy);
				const int j = N * (N * iy + ix);
				std::memcpy(tmp.data(), X + i, copy_sz);
				std::memcpy(X + i, X + j, copy_sz);
				std::memcpy(X + j, tmp.data(), copy_sz);
			}
		}
	} else {
		constexpr int z00 = 0;
		const int do2 = dx >> 1;
		const int dp1 = depth + 1;
		if (depth > MAXPARDEPTH) {
			transpose_xy(X, N, x + z00, y + z00, do2, dp1);
			transpose_xy(X, N, x + do2, y + z00, do2, dp1);
			transpose_xy(X, N, x + z00, y + do2, do2, dp1);
			transpose_xy(X, N, x + do2, y + do2, do2, dp1);
		} else {
			std::array<hpx::future<void>, NCHILD / 2> futs;
			futs[0] = hpx::async(transpose_xy, X, N, x + z00, y + z00, do2, dp1);
			futs[1] = hpx::async(transpose_xy, X, N, x + do2, y + z00, do2, dp1);
			futs[2] = hpx::async(transpose_xy, X, N, x + z00, y + do2, do2, dp1);
			futs[3] = hpx::async(transpose_xy, X, N, x + do2, y + do2, do2, dp1);
			hpx::wait_all(futs.begin(), futs.end());
		}
	}
}

void transpose_yz(float* X, int N) {
	transpose_yz(X, N, 0, 0, 0, N, 0);
}

void transpose_xy(float* X, int N) {
	transpose_xy(X, N, 0, 0, N, 0);
}
