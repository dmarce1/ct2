#include <cosmictiger/constants.hpp>
#include <cosmictiger/options.hpp>
#include <cosmictiger/util.hpp>
#include <cmath>

template<class T>
inline T sqr(T a) {
	return a * a;
}

void cosmos_update0(double& adotdot, double& adot, double& a, double dt) {
	const double omega_m = get_options().omega_m;
	const double omega_r = get_options().omega_r;
	const double omega_lam = get_options().omega_lam;
	const auto H = constants::H0 * get_options().code_to_s * get_options().hubble;
	const double adot0 = adot;

	adotdot = a * sqr(H) * (-omega_r / (sqr(sqr(a))) - 0.5 * omega_m / (sqr(a) * a) + omega_lam);
	const double da1 = adot * a * dt;
	const double dadot1 = adotdot * a * dt;
	a += da1;
	adot += dadot1;

	adotdot = a * sqr(H) * (-omega_r / (sqr(sqr(a))) - 0.5 * omega_m / (sqr(a) * a) + omega_lam);
	const double da2 = adot * a * dt;
	const double dadot2 = adotdot * a * dt;
	a += -da1 + 0.25 * (da1 + da2);
	adot += -dadot1 + 0.25 * (dadot1 + dadot2);

	adotdot = a * sqr(H) * (-omega_r / (sqr(sqr(a))) - 0.5 * omega_m / (sqr(a) * a) + omega_lam);
	const double da3 = adot * a * dt;
	const double dadot3 = adotdot * a * dt;
	a += -0.25 * (da1 + da2) + (1.0 / 6.0) * (da1 + da2 + 4.0 * da3);
	adot += -0.25 * (dadot1 + dadot2) + (1.0 / 6.0) * (dadot1 + dadot2 + 4.0 * dadot3);

}

void cosmos_update(double& adotdot, double& adot, double& a, double dt0) {
	constexpr int N = 128;
	const double dt = dt0 / N;
	for (int i = 0; i < N; i++) {
		cosmos_update0(adotdot, adot, a, dt);
	}
}

double cosmos_dadt(double a) {
	const auto H = constants::H0 * get_options().code_to_s * get_options().hubble;
	const auto omega_m = get_options().omega_m;
	const auto omega_r = get_options().omega_r;
	const auto omega_lam = get_options().omega_lam;
	return H * a * std::sqrt(omega_r / (sqr(sqr(a))) + omega_m / (a * sqr(a)) + omega_lam);
}

static double cosmos_dadtau(double a) {
	return a * cosmos_dadt(a);
}

double cosmos_conformal_time(double a0, double a1) {
	int N = 2;
	double t = 0.0;
	double tlast;
	do {
		tlast = t;
		t = 0.0;
		double a = a0;
		const double loga0 = log(a0);
		const double loga1 = log(a1);
		const double dloga = (loga1 - loga0) / N;
		double loga = loga0;
		a = exp(loga);
		double dadt = cosmos_dadtau(a);
		for (int i = 0; i < N; i++) {
			t += 0.5 * a / dadt * dloga;
			loga = loga0 + dloga * (i + 1);
			a = exp(loga);
			dadt = cosmos_dadtau(a);
			t += 0.5 * a / dadt * dloga;
		}
		N *= 2;
	} while (fabs(tlast - t) > 1e-9 * t);
	return t;
}

double cosmos_time(double a0, double a1) {
	double a = a0;
	double t = 0.0;
	const int N = 1024 * 1024;
	const double loga0 = log(a0);
	const double loga1 = log(a1);
	const double dloga = (loga1 - loga0) / N;
	double loga = loga0;
	a = exp(loga);
	double dadt = cosmos_dadt(a);
	for (int i = 0; i < N; i++) {
		t += 0.5 * a / dadt * dloga;
		loga = loga0 + dloga * (i + 1);
		a = exp(loga);
		dadt = cosmos_dadt(a);
		t += 0.5 * a / dadt * dloga;
	}
	return t;
}


