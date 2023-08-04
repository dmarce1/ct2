#define PARTICLES_CPP

#include <array>
#include <cstdlib>
#include <limits>
#include <vector>
#include <hpx/async.hpp>
#include "cosmictiger/cosmictiger.hpp"
#include "cosmictiger/particle.hpp"
#include "cosmictiger/fft.hpp"

#define MAXPARDEPTH 10

static std::array<fixed32*, NDIM> particle_pos_ptr;
static size_t particle_cnt = 0;
static size_t rung_begin;
static int Ncell;
static void* particle_ptr;
static std::vector<std::vector<std::vector<std::pair<size_t, size_t>>> >cells;
static std::vector<std::vector<std::vector<sfmm::multipole<float, PORDER>>> >Mn;
static std::vector<std::vector<std::vector<sfmm::expansion<float, PORDER>>> >Ln;

size_t particle_count() {
	return particle_cnt;
}

inline static void particle_swap(size_t a, size_t b) {
	std::swap(particle_x(a), particle_x(b));
	std::swap(particle_y(a), particle_y(b));
	std::swap(particle_z(a), particle_z(b));
	std::swap(particle_vx(a), particle_vx(b));
	std::swap(particle_vy(a), particle_vy(b));
	std::swap(particle_vz(a), particle_vz(b));
	std::swap(particle_rung(a), particle_rung(b));
}

static size_t particle_sort_by_pos(std::pair<size_t, size_t> rng, int xdim, fixed32 xmid) {
	size_t lo = rng.first;
	size_t hi = rng.second;
	const fixed32* xptr = particle_pos_ptr[xdim];
	int flops = 0;
	while (lo < hi) {
		if (xptr[lo] >= xmid) {
			while (lo != hi) {
				hi--;
				if (xptr[hi] < xmid) {
					particle_swap(hi, lo);
					break;
				}
			}
		}
		lo++;
	}
	return hi;
}

static void particle_to_grid(std::pair<size_t, size_t> rng, std::array<size_t, NDIM> I, std::array<size_t, NDIM> N, int depth) {
	size_t vol = N[0] * N[1] * N[2];
	if (vol == 1) {
		cells[I[0]][I[1]][I[1]] = rng;
	} else {
		const int xdim = depth % NDIM;
		const int dp1 = depth + 1;
		auto No2 = N;
		auto Il = I;
		auto Ir = I;
		auto rl = rng;
		auto rr = rng;
		No2[xdim] >>= 1;
		Ir[xdim] += No2[xdim];
		fixed32 xmid = (double) (I[xdim] + No2[xdim]) / (double) Ncell;
		rl.second = rr.first = particle_sort_by_pos(rng, xdim, xmid);
		if (depth > MAXPARDEPTH) {
			particle_to_grid(rr, Ir, No2, dp1);
			particle_to_grid(rl, Il, No2, dp1);
		} else {
			auto fut = hpx::async(particle_to_grid, rr, Ir, No2, dp1);
			particle_to_grid(rl, Il, No2, dp1);
			fut.get();
		}
	}
}

static std::vector<std::vector<std::vector<sfmm::expansion<sfmm::complex_float, PORDER>>> >& particle_greens() {
	using namespace sfmm;
	const size_t vol = Ncell * Ncell * Ncell;
	std::vector<float> X(vol);
	std::vector<float> Y(vol);
	std::unordered_map<int, std::vector < std::vector<std::vector<expansion<complex_float, PORDER>>> >> values;
	auto i = values.find(Ncell);
	if( i == values.end()) {
		static const complex_float J(0.0, 1.0);
		std::vector<std::vector<std::vector<expansion<float, PORDER>>>> Gn(Ncell, std::vector < std::vector<expansion<float, PORDER>>>(Ncell, std::vector<expansion<float, PORDER>>(Ncell)));
		std::vector < std::vector<std::vector<expansion<complex_float, PORDER>>>> Gk(Ncell, std::vector < std::vector<expansion<complex_float, PORDER>>>(Ncell, std::vector<expansion<complex_float, PORDER>>(Ncell / 2 + 1)));
		vec3<float> x;
		for (int j = 0; j < Ncell; j++) {
			x[0] = (float) j / Ncell;
			for (int k = 0; k < Ncell; k++) {
				x[1] = (float) k / Ncell;
				for (int l = 0; l < Ncell; l++) {
					x[2] = (float) l / Ncell;
					expansion<float, PORDER> g;
					expansion<float, PORDER> ge;
					if( j * j + k * k + l * l <= 10) {
						for( int n = 0; n < expansion<float, PORDER>::size(); n++) {
							g[n] = 0.0;
						}
					} else {
						greens(g, x);
					}
					greens_ewald(ge, x);
					g += ge;
					for( int n = 0; n < expansion<float, PORDER>::size(); n++) {
						Gn[j][k][l][n] = g[n];
					}
				}
			}
		}
		for (int i = 0; i < expansion<float, PORDER>::size(); i += 2) {
			for (int j = 0; j < Ncell; j++) {
				for (int k = 0; k < Ncell; k++) {
					for (int l = 0; l < Ncell; l++) {
						const int cell = Ncell * (Ncell * j + k) + l;
						X[cell] = Gn[j][k][l][i];
						if (i + 1 < expansion<float, PORDER>::size()) {
							Y[cell] = Gn[j][k][l][i + 1];
						} else {
							Y[cell] = 0.0;
						}
					}
				}
			}
			fft_3d(X.data(), Y.data(), Ncell);
			for (int j0 = 0; j0 < Ncell; j0++) {
				const int j1 = j0 == 0 ? 0 : Ncell - 0;
				for (int k0 = 0; k0 < Ncell; k0++) {
					const int k1 = k0 == 0 ? 0 : Ncell - k0;
					for (int l0 = 0; l0 < Ncell / 2 + 1; l0++) {
						const int l1 = l0 == 0 ? 0 : Ncell - l0;
						const int n0 = Ncell * (Ncell * j0 + k0) + l0;
						const int n1 = Ncell * (Ncell * j1 + k1) + l1;
						complex<float> Z0;
						complex<float> Z1;
						Z0.real() = X[n0];
						Z0.imag() = Y[n0];
						Z1.real() = X[n1];
						Z1.imag() = -Y[n1];
						Gk[j0][k0][l0][i + 0] = (Z0 + Z1) * 0.5f;
						if (i + 1 < expansion<float, PORDER>::size()) {
							Gk[j0][k0][l0][i + 1] = J * (Z1 - Z0) * 0.5f;
						}
					}
				}
			}
		}
		i = values.find(Ncell);
	}
	return i->second;
}

void particle_compute_multipoles() {
	using namespace sfmm;
	Mn.resize(Ncell, std::vector < std::vector<multipole<float, PORDER>>>(Ncell, std::vector<multipole<float, PORDER>>(Ncell)));
	std::vector<hpx::future<void>> futs;
	for (int j = 0; j < Ncell; j++) {
		for (int k = 0; k < Ncell; k++) {
			for (int l = 0; l < Ncell; l++) {
				futs.push_back(hpx::async([j, k, l]() {
					const auto rng = cells[j][k][l];
					multipole < simd_f32, PORDER > M = 0.0;
					int i;
					vec3<simd_fixed32> y;
					y[0] = fixed32((j + 0.5) / Ncell);
					y[1] = fixed32((k + 0.5) / Ncell);
					y[2] = fixed32((l + 0.5) / Ncell);
					for (i = rng.first; i < rng.second; i += SIMD_SIZE) {
						vec3<simd_fixed32> x;
						const simd_f32 m = simd_f32::mask(std::min((size_t) SIMD_SIZE, rng.second - i));
						x[0] = *((simd_fixed32*) &particle_x(i));
						x[1] = *((simd_fixed32*) &particle_y(i));
						x[2] = *((simd_fixed32*) &particle_z(i));
						P2M(M, m, distance(y, x));
					}
					Mn[j][k][l] = reduce_sum(M);
				}));
			}
		}
	}
	hpx::wait_all(futs.begin(), futs.end());
}

void particle_compute_expansions() {
	using namespace sfmm;
	std::vector < std::vector<std::vector<multipole<complex_float, PORDER>>> >Mk(Ncell, std::vector < std::vector<multipole<complex_float, PORDER>>>(Ncell, std::vector<multipole<complex_float, PORDER>>(Ncell / 2 + 1)));
	std::vector < std::vector<std::vector<expansion<complex_float, PORDER>>> >Lk(Ncell, std::vector < std::vector<expansion<complex_float, PORDER>>>(Ncell, std::vector<expansion<complex_float, PORDER>>(Ncell / 2 + 1)));
	Ln.resize(Ncell, std::vector < std::vector<expansion<float, PORDER>>>(Ncell, std::vector<expansion<float, PORDER>>(Ncell)));
	const size_t vol = Ncell * Ncell * Ncell;
	std::vector<float> X(vol);
	std::vector<float> Y(vol);
	static const complex<float> J(0.0, 1.0);
	for (int i = 0; i < multipole<float, PORDER>::size(); i += 2) {
		for (int j = 0; j < Ncell; j++) {
			for (int k = 0; k < Ncell; k++) {
				for (int l = 0; l < Ncell; l++) {
					const int cell = Ncell * (Ncell * j + k) + l;
					X[cell] = Mn[j][k][l][i];
					if (i + 1 < multipole<float, PORDER>::size()) {
						Y[cell] = Mn[j][k][l][i + 1];
					} else {
						Y[cell] = 0.0;
					}
				}
			}
		}
		fft_3d(X.data(), Y.data(), Ncell);
		for (int j0 = 0; j0 < Ncell; j0++) {
			const int j1 = j0 == 0 ? 0 : Ncell - 0;
			for (int k0 = 0; k0 < Ncell; k0++) {
				const int k1 = k0 == 0 ? 0 : Ncell - k0;
				for (int l0 = 0; l0 < Ncell / 2 + 1; l0++) {
					const int l1 = l0 == 0 ? 0 : Ncell - l0;
					const int n0 = Ncell * (Ncell * j0 + k0) + l0;
					const int n1 = Ncell * (Ncell * j1 + k1) + l1;
					complex<float> Z0;
					complex<float> Z1;
					Z0.real() = X[n0];
					Z0.imag() = Y[n0];
					Z1.real() = X[n1];
					Z1.imag() = -Y[n1];
					Mk[j0][k0][l0][i] = (Z0 + Z1) * 0.5f;
					if (i + 1 < multipole<float, PORDER>::size()) {
						Mk[j0][k0][l0][i + 1] = J * (Z1 - Z0) * 0.5f;
					}
				}
			}
		}
	}
	const auto& Gk = particle_greens();
	for (int j = 0; j < Ncell; j++) {
		for (int k = 0; k < Ncell; k++) {
			for (int l = 0; l < Ncell; l++) {
				MG2L(Lk[j][k][l], Mk[j][k][l], Gk[j][k][l]);
			}
		}
	}
	for (int i = 0; i < expansion<float, PORDER>::size(); i += 2) {
		for (int j = 0; j < Ncell; j++) {
			for (int k = 0; k < Ncell; k++) {
				for (int l = 0; l < Ncell; l++) {
					const int cell = Ncell * (Ncell * j + k) + l;
					X[cell] = Gk[j][k][l][i].real();
					Y[cell] = Gk[j][k][l][i].imag();
					if (i + 1 < expansion<float, PORDER>::size()) {
						X[cell] -= Gk[j][k][l][i + 1].imag();
						Y[cell] += Gk[j][k][l][i + 1].real();
					}
				}
			}
		}
		fft_3d(X.data(), Y.data(), Ncell);
		for (int j0 = 0; j0 < Ncell; j0++) {
			for (int k0 = 0; k0 < Ncell; k0++) {
				for (int l0 = 0; l0 < Ncell / 2 + 1; l0++) {
					const int n0 = Ncell * (Ncell * j0 + k0) + l0;
					Ln[j0][k0][l0][i] = X[n0];
					if (i + 1 < expansion<float, PORDER>::size()) {
						Ln[j0][k0][l0][i + 1] = Y[n0];
					}
				}
			}
		}
	}
}

void particle_allocate(size_t count) {
	particle_cnt = count;
	const size_t nbytes = count * (NDIM * (sizeof(fixed32) + sizeof(float) + sizeof(char)));
	char* ptr = new char[nbytes];
	particle_x_ptr = (fixed32*) ptr;
	particle_y_ptr = (fixed32*) (ptr + NDIM * count);
	particle_z_ptr = (fixed32*) (ptr + 2 * NDIM * count);
	particle_vx_ptr = (float*) (ptr + 3 * NDIM * count);
	particle_vy_ptr = (float*) (ptr + 4 * NDIM * count);
	particle_vz_ptr = (float*) (ptr + 5 * NDIM * count);
	particle_rung_ptr = (char*) (ptr + 6 * NDIM * count);
	particle_pos_ptr[0] = particle_x_ptr;
	particle_pos_ptr[1] = particle_y_ptr;
	particle_pos_ptr[2] = particle_z_ptr;
	particle_ptr = ptr;
}

void particle_sort_by_rung(int minrung) {
	if (minrung == 0) {
		rung_begin = 0;
	} else {
		size_t lo = rung_begin;
		size_t hi = particle_count();
		while (lo < hi) {
			if (particle_rung(lo) >= minrung) {
				while (lo != hi) {
					hi--;
					if (particle_rung(hi) < minrung) {
						particle_swap(hi, lo);
						break;
					}
				}
			}
			lo++;
		}
		rung_begin = hi;
	}
}

void particle_to_grid(int ncell) {
	cells.resize(Ncell, std::vector < std::vector<std::pair<size_t, size_t>>>(Ncell, std::vector<std::pair<size_t, size_t >>(Ncell)));
	Ncell = ncell;
	std::pair < size_t, size_t > rng;
	std::array < size_t, NDIM > I;
	std::array < size_t, NDIM > N;
	rng.first = rung_begin;
	rng.second = particle_cnt;
	for (int n = 0; n < NDIM; n++) {
		I[n] = 0;
		N[n] = Ncell;
	}
	particle_to_grid(rng, I, N, 0);
}
