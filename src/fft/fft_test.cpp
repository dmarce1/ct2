#include <cosmictiger/fft.hpp>
#include <cosmictiger/cosmictiger.hpp>
#include <hpx/hpx.hpp>
#include <hpx/hpx_init.hpp>
#include <fftw3.h>
#include <sfmm.hpp>
#include <cosmictiger/util.hpp>

#define TEST_MIN 8
#define TEST_MAX 4096

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

void fftw_r2c_3d(const std::vector<fft_float>& xin, std::vector<complex<fft_float>>& xout) {
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
	std::vector < complex < fft_simd >> X(N);
	std::vector < complex < fft_float >> Z(N);
	for (int n = 0; n < N; n++) {
		X[n].real() = Z[n].real() = rand1();
		X[n].imag() = Z[n].imag() = rand1();
	}
	scramble((fft_float*) X.data(), 1, N, 2 * FFT_SIMD_SIZE);
	fft_complex(X.data(), N);
	fftw_cmplx (Z);
	double error = 0.0;
	for (int n = 0; n < N; n++) {
		const auto dx = Z[n].real() - X[n].real()[0];
		const auto dy = Z[n].imag() - X[n].imag()[0];
		error += sqr(dx);
		error += sqr(dy);
	}
	error /= 2 * N;
	error = std::sqrt(error);
	return error;
}

double test_real(int N) {
	std::vector<fft_simd> X(N);
	std::vector<fft_float> Zin(N);
	std::vector < complex < fft_float >> Zout(N / 2 + 1);
	for (int n = 0; n < N; n++) {
		X[n] = Zin[n] = rand1();
	}
	scramble((fft_float*) X.data(), 1, N, FFT_SIMD_SIZE);
	fft_real(X.data(), N);
	fftw_real(Zin, Zout);
	double error = 0.0;
	for (int n = 0; n < N / 2 + 1; n++) {
		if (n != 0 && n != N / 2) {
			const auto dx = Zout[n].real() - X[n][0];
			const auto dy = Zout[n].imag() - X[N - n][0];
			error += sqr(dx);
			error += sqr(dy);
			//	printf("%i | %15e %15e | %15e %15e | %15e %15e\n", n, X[FFT_SIMD_SIZE * n], X[FFT_SIMD_SIZE * (N - n)], Zout[n].real(), Zout[n].imag(), dx, dy);
		} else {
			const auto dx = Zout[n].real() - X[n][0];
			error += sqr(dx);
			//	printf("%i | %15e %15s | %15e %15s | %15e %15s\n", n, X[FFT_SIMD_SIZE * n], "", Zout[n].real(), "", dx, "");
		}
	}
	error /= N;
	error = std::sqrt(error);
	return error;
}

double test_real_scalar(int N) {
	std::vector<fft_float> X(N);
	std::vector<fft_float> Zin(N);
	std::vector < complex < fft_float >> Zout(N / 2 + 1);
	for (int n = 0; n < N; n++) {
		Zin[n] = X[n] = rand1();
	}
	scramble(X.data(), 1, N, 1);
	fft_real(X.data(), N);
	fftw_real(Zin, Zout);
	double error = 0.0;
	for (int n = 0; n < N / 2 + 1; n++) {
		if (n != 0 && n != N / 2) {
			const auto dx = Zout[n].real() - X[n];
			const auto dy = Zout[n].imag() - X[N - n];
			error += sqr(dx);
			error += sqr(dy);
			//	printf("%i | %15e %15e | %15e %15e | %15e %15e\n", n, X[FFT_SIMD_SIZE * n], X[FFT_SIMD_SIZE * (N - n)], Zout[n].real(), Zout[n].imag(), dx, dy);
		} else {
			const auto dx = Zout[n].real() - X[n];
			error += sqr(dx);
			//	printf("%i | %15e %15s | %15e %15s | %15e %15s\n", n, X[FFT_SIMD_SIZE * n], "", Zout[n].real(), "", dx, "");
		}
	}
	error /= N;
	error = std::sqrt(error);
	return error;
}

/*double test_r2c_3d(int N) {
 const int N3 = N * N * N;
 const int N3cmplx = N * N * (N / 2 + 1);
 std::vector<fft_float> Xin(N3 * FFT_SIMD_SIZE);
 std::vector<fft_float> Xout(2 * N3cmplx);
 std::vector<fft_float> Zin(N3);
 std::vector < complex < fft_float >> Zout(N3cmplx);
 for (int n = 0; n < N * N * N; n++) {
 //		Zin[n] = Xin[n] = rand1();
 Zin[n] = Xin[n] = 0.0;
 }
 Xin[0] = Zin[0] = 1.0;
 fft_r2c_3d(Xin.data(), Xout.data(), N);
 fftw_r2c_3d(Zin, Zout);
 double error = 0.0;
 for (int i = 0; i < N; i++) {
 for (int j = 0; j < N; j++) {
 for (int k = 0; k < N / 2 + 1; k++) {
 const int ijk = (N / 2 + 1) * (N * i + j) + k;
 if (ijk != 0 && ijk != N * N * (N / 2 + 1)) {
 const auto dx = Zout[ijk].real() - Xout[ijk];
 const auto dy = Zout[ijk].imag() - Xout[ijk + N3cmplx];
 error += sqr(dx);
 error += sqr(dy);
 printf("%i %i %i | %15e %15e | %15e %15e | %15e %15e\n", i, j, k, Xout[ijk], Xout[ijk + N3cmplx], Zout[ijk].real(), Zout[ijk].imag(), dx, dy);
 } else {
 const auto dx = Zout[ijk].real() - Xout[ijk];
 error += sqr(dx);
 printf("%i %i %i | %15e %15s | %15e %15s | %15e %15s\n", i, j, k, Xout[ijk], "", Zout[ijk].real(), "", dx, "");
 }
 }
 }
 }
 error /= N;
 error = std::sqrt(error);
 return error;
 }*/

int hpx_main(int argc, char *argv[]) {
	/*printf("\n3D\n");
	 for (int n = TEST_MIN; n <= TEST_MAX; n *= 2) {
	 //	const double err = test_r2c_3d(n);
	 //	printf("%i %e\n", n, err);
	 }*/
	printf("\nReal (Scalar)\n");
	for (int n = TEST_MIN; n <= TEST_MAX; n *= 2) {
		const double err = test_real_scalar(n);
		printf("%i %e\n", n, err);
	}
	printf("\nReal (Vector)\n");
	for (int n = TEST_MIN; n <= TEST_MAX; n *= 2) {
		const double err = test_real(n);
		printf("%i %e\n", n, err);
	}
	printf("\nComplex\n");
	for (int n = TEST_MIN; n <= TEST_MAX; n *= 2) {
		const double err = test_complex(n);
		printf("%i %e\n", n, err);
	}
	return hpx::finalize();
}

int main(int argc, char *argv[]) {
	std::vector<std::string> cfg = {"hpx.commandline.allow_unknown=1"};
	return hpx::init(argc, argv);
}
