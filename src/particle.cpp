#define PARTICLES_CPP

#include <array>
#include <cstdlib>
#include <limits>
#include <vector>
#include <hpx/async.hpp>
#include "cosmictiger/cosmictiger.hpp"
#include "cosmictiger/particle.hpp"
#include "cosmictiger/fft.hpp"
#include "cosmictiger/options.hpp"
#include "cosmictiger/util.hpp"
#include <stack>

#define MAXPARDEPTH 10

static std::array<pos_fixed*, NDIM> particle_pos_ptr;
static size_t particle_cnt = 0;
static size_t rung_begin;
static std::stack<size_t> rung_begins;
static int Ncell;
static void* particle_ptr;
static std::vector<std::vector<std::vector<std::pair<size_t, size_t>>> >cells;
static std::vector<std::vector<std::vector<sfmm::multipole<fft_float, PORDER>>> >Mn;
static std::vector<std::vector<std::vector<sfmm::expansion<fft_float, PORDER>>> >Ln;

size_t particle_count() {
	return particle_cnt;
}

void particle_pop_rungs() {
	rung_begin = rung_begins.top();
	rung_begins.pop();
}

void particle_push_rungs() {
	rung_begins.push(rung_begin);
}

std::vector<size_t> particle_rung_counts() {
	std::vector < size_t > counts;
	for (size_t i = 0; i < particle_cnt; i++) {
		const auto rung = particle_rung(i);
		if (counts.size() <= rung) {
			counts.resize(rung + 1, 0);
		}
		counts[rung]++;
	}
	return counts;
}

std::pair<size_t, size_t> particle_current_range() {
	std::pair < size_t, size_t > rc;
	rc.first = rung_begin;
	rc.second = particle_cnt;
	return rc;
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

static size_t particle_sort_by_pos(std::pair<size_t, size_t> rng, int xdim, pos_fixed xmid) {
	size_t lo = rng.first;
	size_t hi = rng.second;
	const pos_fixed* xptr = particle_pos_ptr[xdim];
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
		pos_fixed xmid = (double) (I[xdim] + No2[xdim]) / (double) Ncell;
		rl.second = rr.first = particle_sort_by_pos(rng, xdim, xmid);
		if (depth > MAXPARDEPTH) {
			particle_to_grid(rr, Ir, No2, dp1);
			particle_to_grid(rl, Il, No2, dp1);
		} else {
			auto fut = hpx::async([rr, Ir, No2, dp1]() {particle_to_grid(rr, Ir, No2, dp1);});
			particle_to_grid(rl, Il, No2, dp1);
			fut.get();
		}
	}
}

static std::vector<std::vector<std::vector<sfmm::expansion<sfmm::complex<fft_float>, PORDER>>> >& particle_greens() {
	using namespace sfmm;
	const size_t vol = Ncell * Ncell * Ncell;
	std::vector<fft_float> X(vol);
	std::vector<fft_float> Y(vol + Ncell * Ncell);
	std::unordered_map<int, std::vector < std::vector<std::vector<expansion<complex<fft_float>, PORDER>>> >> values;
	auto i = values.find(Ncell);
	if( i == values.end()) {
		static const complex<fft_float> J(0.0, 1.0);
		std::vector<std::vector<std::vector<expansion<fft_float, PORDER>>>> Gn(Ncell, std::vector < std::vector<expansion<fft_float, PORDER>>>(Ncell, std::vector<expansion<fft_float, PORDER>>(Ncell)));
		std::vector < std::vector<std::vector<expansion<complex<fft_float>, PORDER>>>> Gk(Ncell, std::vector < std::vector<expansion<complex<fft_float>, PORDER>>>(Ncell, std::vector<expansion<complex<fft_float>, PORDER>>(Ncell / 2 + 1)));
		vec3<fft_float> x;
		for (int j = 0; j < Ncell; j++) {
			x[0] = (fft_float) j / Ncell;
			for (int k = 0; k < Ncell; k++) {
				x[1] = (fft_float) k / Ncell;
				for (int l = 0; l < Ncell; l++) {
					x[2] = (fft_float) l / Ncell;
					expansion<fft_float, PORDER> g;
					expansion<fft_float, PORDER> ge;
					if( j * j + k * k + l * l <= 10) {
						for( int n = 0; n < expansion<fft_float, PORDER>::size(); n++) {
							g[n] = 0.0;
						}
					} else {
						greens(g, x);
					}
					greens_ewald(ge, x);
					g += ge;
					for( int n = 0; n < expansion<fft_float, PORDER>::size(); n++) {
						Gn[j][k][l][n] = g[n];
					}
				}
			}
		}
		const int dimag = Ncell * Ncell * (Ncell / 2 + 1);
		const fft_float norm = 1.0 / (Ncell * Ncell * Ncell);
		for (int i = 0; i < expansion<fft_float, PORDER>::size(); i ++) {
			for (int j = 0; j < Ncell; j++) {
				for (int k = 0; k < Ncell; k++) {
					for (int l = 0; l < Ncell; l++) {
						const int cell = Ncell * (Ncell * j + k) + l;
						X[cell] = Gn[j][k][l][i];
					}
				}
			}
			fft_r2c_3d(X.data(), Y.data(), Ncell);
			for (int j0 = 0; j0 < Ncell; j0++) {
				for (int k0 = 0; k0 < Ncell; k0++) {
					for (int l0 = 0; l0 < Ncell / 2 + 1; l0++) {
						const int n0 = Ncell * (Ncell * j0 + k0) + l0;
						complex<fft_float> Z;
						Z.real() = Y[n0];
						Z.imag() = Y[n0 + dimag];
						Gk[j0][k0][l0][i + 0] = Z * norm;
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
	Mn.resize(Ncell, std::vector < std::vector<multipole<fft_float, PORDER>>>(Ncell, std::vector<multipole<fft_float, PORDER>>(Ncell)));
	std::vector<hpx::future<void>> futs;
	for (int j = 0; j < Ncell; j++) {
		for (int k = 0; k < Ncell; k++) {
			for (int l = 0; l < Ncell; l++) {
				futs.push_back(hpx::async([j, k, l]() {
					const auto rng = cells[j][k][l];
					multipole < simd_f32, PORDER > M = 0.0;
					int i;
					vec3<pos_simd> y;
					y[0] = pos_fixed((j + 0.5) / Ncell);
					y[1] = pos_fixed((k + 0.5) / Ncell);
					y[2] = pos_fixed((l + 0.5) / Ncell);
					for (i = rng.first; i < rng.second; i += FFT_SIMD_SIZE) {
						vec3<pos_simd> x;
						const simd_f32 m = simd_f32::mask(std::min((size_t) FFT_SIMD_SIZE, rng.second - i));
						x[0] = *((pos_simd*) &particle_x(i));
						x[1] = *((pos_simd*) &particle_y(i));
						x[2] = *((pos_simd*) &particle_z(i));
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
	std::vector < std::vector<std::vector<multipole<complex<fft_float>, PORDER>>> >Mk(Ncell, std::vector < std::vector<multipole<complex<fft_float>, PORDER>>>(Ncell, std::vector<multipole<complex<fft_float>, PORDER>>(Ncell / 2 + 1)));
	std::vector < std::vector<std::vector<expansion<complex<fft_float>, PORDER>>> >Lk(Ncell, std::vector < std::vector<expansion<complex<fft_float>, PORDER>>>(Ncell, std::vector<expansion<complex<fft_float>, PORDER>>(Ncell / 2 + 1)));
	Ln.resize(Ncell, std::vector < std::vector<expansion<fft_float, PORDER>>>(Ncell, std::vector<expansion<fft_float, PORDER>>(Ncell)));
	const size_t vol = Ncell * Ncell * Ncell;
	std::vector<fft_float> X(vol);
	std::vector<fft_float> Y(vol + Ncell * Ncell);
	static const complex<fft_float> J(0.0, 1.0);
	const int dimag = Ncell * Ncell * (Ncell / 2 + 1);
	for (int i = 0; i < multipole<fft_float, PORDER>::size(); i ++) {
		for (int j = 0; j < Ncell; j++) {
			for (int k = 0; k < Ncell; k++) {
				for (int l = 0; l < Ncell; l++) {
					const int cell = Ncell * (Ncell * j + k) + l;
					X[cell] = Mn[j][k][l][i];
				}
			}
		}
		fft_r2c_3d(X.data(), Y.data(), Ncell);
		for (int j0 = 0; j0 < Ncell; j0++) {
			for (int k0 = 0; k0 < Ncell; k0++) {
				for (int l0 = 0; l0 < Ncell / 2 + 1; l0++) {
					const int n0 = Ncell * (Ncell * j0 + k0) + l0;
					complex<fft_float> Z;
					Z.real() = Y[n0];
					Z.imag() = Y[n0 + dimag];
					Mk[j0][k0][l0][i] = Z;
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
	for (int i = 0; i < expansion<fft_float, PORDER>::size(); i++) {
		for (int j = 0; j < Ncell; j++) {
			for (int k = 0; k < Ncell; k++) {
				for (int l = 0; l < Ncell; l++) {
					const int cell = Ncell * (Ncell * j + k) + l;
					X[cell] = Lk[j][k][l][i].real();
					X[cell + dimag] = Lk[j][k][l][i].imag();
				}
			}
		}
		fft_c2r_3d(X.data(), Y.data(), Ncell);
		for (int j0 = 0; j0 < Ncell; j0++) {
			for (int k0 = 0; k0 < Ncell; k0++) {
				for (int l0 = 0; l0 < Ncell / 2 + 1; l0++) {
					const int n0 = Ncell * (Ncell * j0 + k0) + l0;
					Ln[j0][k0][l0][i] = Y[n0];
				}
			}
		}
	}
}

void particle_allocate(size_t count) {
	particle_cnt = count;
	const size_t nbytes = count * (NDIM * (sizeof(pos_fixed) + sizeof(force_float) + sizeof(char)));
	char* ptr = new char[nbytes];
	particle_x_ptr = (pos_fixed*) ptr;
	particle_y_ptr = (pos_fixed*) (ptr + NDIM * count);
	particle_z_ptr = (pos_fixed*) (ptr + 2 * NDIM * count);
	particle_vx_ptr = (force_float*) (ptr + 3 * NDIM * count);
	particle_vy_ptr = (force_float*) (ptr + 4 * NDIM * count);
	particle_vz_ptr = (force_float*) (ptr + 5 * NDIM * count);
	particle_rung_ptr = (char*) (ptr + 6 * NDIM * count);
	particle_pos_ptr[0] = particle_x_ptr;
	particle_pos_ptr[1] = particle_y_ptr;
	particle_pos_ptr[2] = particle_z_ptr;
	particle_ptr = ptr;
}

int particle_cell_dim() {
	return Ncell;
}

std::pair<size_t, size_t> particle_get_cell(int i, int j, int k) {
	return cells[i][j][k];
}

sfmm::expansion<fft_float, PORDER> particle_get_expansion(int i, int j, int k) {
	return Ln[i][j][k];
}

std::vector<std::pair<size_t, size_t>> particle_cell_neighbors(int i, int j, int k) {
	const auto bwidth = get_options().bwidth;
	std::vector<std::pair<size_t, size_t>> neighbors;
	for (int i0 = -bwidth; i0 <= bwidth; i0++) {
		const int i1 = (i + i0 + Ncell) % Ncell;
		for (int j0 = -bwidth; j0 <= bwidth; j0++) {
			const int j1 = (j + j0 + Ncell) % Ncell;
			for (int k0 = -bwidth; k0 <= bwidth; k0++) {
				const int k1 = (k + k0 + Ncell) % Ncell;
				if (i1 * i1 + j1 * j1 + k1 * k1 <= bwidth * bwidth) {
					neighbors.push_back(cells[i][j][k]);
				}
			}
		}
	}
	return neighbors;
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

void particle_set_minrung(int minrung) {
	for (size_t i = 0; i < particle_cnt; i++) {
		auto& rung = particle_rung(i);
		rung = std::max((int) rung, minrung);
	}
}

void particle_to_grid() {
	const int ppcell = get_options().ppcell;
	const std::pair<size_t, size_t> rng = particle_current_range();
	const size_t count = rng.second - rng.first;
	std::array < size_t, NDIM > I;
	std::array < size_t, NDIM > N;
	Ncell = std::pow(count / ppcell, 1.0 / 3.0) + 4;
	cells.resize(Ncell, std::vector < std::vector<std::pair<size_t, size_t>>>(Ncell, std::vector<std::pair<size_t, size_t>>(Ncell)));
	for (int n = 0; n < NDIM; n++) {
		I[n] = 0;
		N[n] = Ncell;
	}
	particle_to_grid(rng, I, N, 0);
}

void particles_init() {
	for (int n = 0; n < particle_cnt; n++) {
		particle_x(n) = rand1();
		particle_y(n) = rand1();
		particle_z(n) = rand1();
		particle_vx(n) = 0.f;
		particle_vy(n) = 0.f;
		particle_vz(n) = 0.f;
		particle_rung(n) = 0;
	}
}
