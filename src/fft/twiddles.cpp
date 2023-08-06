#include <cosmictiger/cosmictiger.hpp>
#include <cmath>
#include <unordered_map>
#include <vector>

const std::vector<fft_float>& cos_twiddles(int N) {
	constexpr fft_float twopi = 2.0 * M_PI;
	static std::unordered_map<int, std::vector<fft_float>> values;
	auto i = values.find(N);
	if (i == values.end()) {
		std::vector<fft_float> W(N);
		int j = 0;
		for (int n = 0; n < N; n++) {
			W[n] = cos(-twopi * n / N);
		}
		values[N] = std::move(W);
		i = values.find(N);
	}
	return i->second;
}


const std::vector<fft_float>& sin_twiddles(int N) {
	constexpr fft_float twopi = 2.0 * M_PI;
	static std::unordered_map<int, std::vector<fft_float>> values;
	auto i = values.find(N);
	if (i == values.end()) {
		std::vector<fft_float> W(N);
		int j = 0;
		for (int n = 0; n < N; n++) {
			W[n] = sin(-twopi * n / N);
		}
		values[N] = std::move(W);
		i = values.find(N);
	}
	return i->second;
}


const std::vector<fft_float>& inv_sin_twiddles(int N) {
	constexpr fft_float twopi = 2.0 * M_PI;
	static std::unordered_map<int, std::vector<fft_float>> values;
	auto i = values.find(N);
	if (i == values.end()) {
		std::vector<fft_float> W(N);
		int j = 0;
		for (int n = 0; n < N; n++) {
			W[n] = sin(twopi * n / N);
		}
		values[N] = std::move(W);
		i = values.find(N);
	}
	return i->second;
}


