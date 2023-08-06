#include <cosmictiger/cosmictiger.hpp>
#include <cstring>
#include <unordered_map>
#include <vector>

#include <hpx/async.hpp>

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

void scramble(fft_float* X, int N1, int N2, int N3) {
	static thread_local std::vector<fft_float> tmp;
	const int copy_sz = sizeof(fft_float) * N3;
	for (int i = 0; i < N1; i++) {
		tmp.resize(N3);
		int l = 0;
		for (int j = 0; j < N2; j++) {
			const int n = N3 * (N2 * i + j);
			const int m = N3 * (N2 * i + l);
			if (n > m) {
				std::memcpy(tmp.data(), X + n, copy_sz);
				std::memcpy(X + n, X + m, copy_sz);
				std::memcpy(X + m, tmp.data(), copy_sz);
			}
			l = reverse_increment(l, N2);
		}
	}
}

