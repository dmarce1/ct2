#include <cosmictiger/fft.hpp>
#include <cosmictiger/cosmictiger.hpp>
#include <hpx/hpx.hpp>
#include <hpx/hpx_init.hpp>
#include <fftw3.h>
#include <sfmm.hpp>
#include <cosmictiger/util.hpp>

using namespace sfmm;

void fftw_real(const std::vector<fft_float>& xin, std::vector<complex<fft_float>>& xout) {
	const int N = xin.size();
	static std::unordered_map<int, fftw_plan> plans;
	static std::unordered_map<int, fft_float*> in;
	static std::unordered_map<int, fftw_complex*> out;
	if (plans.find(N) == plans.end()) {
		in[N] = (fft_float*) malloc(sizeof(fft_float) * N);
		out[N] = (fftw_complex*) malloc(sizeof(fftw_complex) * (N / 2 + 1));
		plans[N] = fftw_plan_dft_r2c_1d(N, in[N], out[N], FFTW_MEASURE);
	}
	auto* i = in[N];
	auto* o = out[N];
	for (int n = 0; n < N; n++) {
		i[n] = xin[n];
	}
	fftw_execute(plans[N]);
	for (int n = 0; n < N / 2 + 1; n++) {
		xout[n].real() = (o[n][0]);
		xout[n].imag() = (o[n][1]);
	}
}

void fftw_real_3d(const std::vector<fft_float>& xin, std::vector<complex<fft_float>>& xout) {
	const int N = std::lround(std::pow(xin.size(), 1.0 / 3.0));
	static std::unordered_map<int, fftw_plan> plans;
	static std::unordered_map<int, fft_float*> in;
	static std::unordered_map<int, fftw_complex*> out;
	if (plans.find(N) == plans.end()) {
		in[N] = (fft_float*) malloc(sizeof(fft_float) * N * N * N);
		out[N] = (fftw_complex*) malloc(sizeof(fftw_complex) * N * N * (N / 2 + 1));
		plans[N] = fftw_plan_dft_r2c_3d(N, N, N, in[N], out[N], FFTW_MEASURE);
	}
	auto* i = in[N];
	auto* o = out[N];
	for (int n = 0; n < N; n++) {
		i[n] = xin[n];
	}
	fftw_execute(plans[N]);
	for (int n = 0; n < N * N * (N / 2 + 1); n++) {
		xout[n].real() = (o[n][0]);
		xout[n].imag() = (o[n][1]);
	}
}

void fftw_cmplx(std::vector<complex<fft_float>>& x) {
	const int N = x.size();
	static std::unordered_map<int, fftw_plan> plans;
	static std::unordered_map<int, fftw_complex*> in;
	static std::unordered_map<int, fftw_complex*> out;
	if (plans.find(N) == plans.end()) {
		in[N] = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
		out[N] = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
		plans[N] = fftw_plan_dft_1d(N, in[N], out[N], FFTW_FORWARD, FFTW_MEASURE);
	}
	auto* i = in[N];
	auto* o = out[N];
	for (int n = 0; n < N; n++) {
		i[n][0] = x[n].real();
		i[n][1] = x[n].imag();
	}
	fftw_execute(plans[N]);
	for (int n = 0; n < N; n++) {
		x[n].real() = (o[n][0]);
		x[n].imag() = (o[n][1]);
	}
}

double test_complex(int N) {
	std::vector<fft_float> X(N * FFT_SIMD_SIZE), Y(N * FFT_SIMD_SIZE);
	std::vector < complex < fft_float >> Z(N);
	for (int n = 0; n < FFT_SIMD_SIZE * N; n++) {
		if (n % FFT_SIMD_SIZE == 0) {
			Z[n / FFT_SIMD_SIZE].real() = X[n] = rand1();
			Z[n / FFT_SIMD_SIZE].imag() = Y[n] = rand1();
		} else {
			X[n] = X[n - 1];
			Y[n] = Y[n - 1];
		}
	}
	scramble(X.data(), 1, N, FFT_SIMD_SIZE);
	scramble(Y.data(), 1, N, FFT_SIMD_SIZE);
	fft_complex(X.data(), Y.data(), N, FFT_SIMD_SIZE);
	fftw_cmplx (Z);
	double error = 0.0;
	for (int n = 0; n < N; n++) {
		const auto dx = Z[n].real() - X[FFT_SIMD_SIZE * n];
		const auto dy = Z[n].imag() - Y[FFT_SIMD_SIZE * n];
		error += sqr(dx);
		error += sqr(dy);
	}
	error /= 2 * N;
	error = std::sqrt(error);
	return error;
}

double test_real(int N) {
	std::vector<fft_float> X(N * FFT_SIMD_SIZE);
	std::vector<fft_float> Zin(N);
	std::vector < complex < fft_float >> Zout(N / 2 + 1);
	for (int n = 0; n < FFT_SIMD_SIZE * N; n++) {
		if (n % FFT_SIMD_SIZE == 0) {
			Zin[n / FFT_SIMD_SIZE] = X[n] = rand1();
		} else {
			X[n] = X[n - 1];
		}
	}
	scramble(X.data(), 1, N, FFT_SIMD_SIZE);
	fft_real(X.data(), N, FFT_SIMD_SIZE);
	fftw_real(Zin, Zout);
	double error = 0.0;
	for (int n = 0; n < N / 2 + 1; n++) {
		if (n != 0 && n != N / 2) {
			const auto dx = Zout[n].real() - X[FFT_SIMD_SIZE * n];
			const auto dy = Zout[n].imag() - X[FFT_SIMD_SIZE * (N - n)];
			error += sqr(dx);
			error += sqr(dy);
			//	printf("%i | %15e %15e | %15e %15e | %15e %15e\n", n, X[FFT_SIMD_SIZE * n], X[FFT_SIMD_SIZE * (N - n)], Zout[n].real(), Zout[n].imag(), dx, dy);
		} else {
			const auto dx = Zout[n].real() - X[FFT_SIMD_SIZE * n];
			error += sqr(dx);
			//	printf("%i | %15e %15s | %15e %15s | %15e %15s\n", n, X[FFT_SIMD_SIZE * n], "", Zout[n].real(), "", dx, "");
		}
	}
	error /= N;
	error = std::sqrt(error);
	return error;
}

int hpx_main(int argc, char *argv[]) {
	printf("Real\n");
	for (int n = 2; n <= 16384; n *= 2) {
		const double err = test_real(n);
		printf("%i %e\n", n, err);
	}
	printf("\nComplex\n");
	for (int n = 2; n <= 16384; n *= 2) {
		const double err = test_complex(n);
		printf("%i %e\n", n, err);
	}
	return hpx::finalize();
}

int main(int argc, char *argv[]) {
	std::vector<std::string> cfg = {"hpx.commandline.allow_unknown=1"};
	return hpx::init(argc, argv);
}
