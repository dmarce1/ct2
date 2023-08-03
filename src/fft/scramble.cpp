#include <cstring>
#include <unordered_map>
#include <vector>

#include <hpx/async.hpp>

static int reverse_increment(int i, int N) {
	const int Nm1 = N - 1;
	int j, k;
	j = ~i & Nm1;
	k = __builtin_ffs(j) - 1;
	j = (1 << k) - 1;
	k = i & j;
	k &= j;
	k |= (j + 1);
	return k;
}

void scramble(float* X, int N) {
	std::vector<hpx::future<void>> futs(N);
	const int copy_sz = sizeof(float) * N;
	for (int i = 0; i < N; i++) {
		futs[i] = hpx::async([N, i, X, copy_sz]() {
			static thread_local std::vector<float> tmp;
			tmp.resize(N);
			int l = 0;
			for( int j = 0; j < N; j++) {
				const int n = N * (N * i + j);
				const int m = N * (N * i + l);
				std::memcpy(tmp.data(), X + n, copy_sz);
				std::memcpy(X + n, X + m, copy_sz);
				std::memcpy(X + m, tmp.data(), copy_sz);
				l = reverse_increment(l, N);
			}
		});
	}
	hpx::wait_all(futs.begin(), futs.end());
}

