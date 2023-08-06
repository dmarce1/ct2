#include <cosmictiger/kick.hpp>

int kick_cell(const std::pair<size_t, size_t>& snk, const sfmm::vec3<double> C0, const std::vector<std::pair<size_t, size_t>>& srcs, const sfmm::expansion<float, PORDER>& L0, const kick_params& params) {
	using namespace sfmm;
	const float h = get_options().hsoft;
	const float GM = get_options().GM;
	const float h2 = h * h;
	const float hinv = 1.0f / h;
	const float hinv3 = sqr(hinv) * hinv;
	int max_rung = 0;
	expansion < force_simd, PORDER > L;
	vec3 < simd_fixed32 > C;
	for (int dim = 0; dim < NDIM; dim++) {
		C[dim] = fixed32(C0[dim]);
	}
	for (int n = 0; n < expansion<float, PORDER>::size(); n++) {
		L[n] = L0[n];
	}
	for (size_t n = snk.first; n < snk.second; n += FORCE_SIMD_SIZE) {
		force_type < force_simd > f;
		f.init();
		vec3 < force_simd > dx;
		const simd_fixed32 x = *((simd_fixed32*) (&particle_x(n)));
		const simd_fixed32 y = *((simd_fixed32*) (&particle_y(n)));
		const simd_fixed32 z = *((simd_fixed32*) (&particle_z(n)));
		dx[0] = distance(x, C[0]);
		dx[1] = distance(y, C[1]);
		dx[2] = distance(z, C[2]);
		L2P(f, L, dx);
		for (size_t i = n; i < std::min(n + FORCE_SIMD_SIZE, snk.second); i++) {
			force_type < force_simd > df;
			df.init();
			const simd_fixed32 x0 = particle_x(i);
			const simd_fixed32 y0 = particle_y(i);
			const simd_fixed32 z0 = particle_z(i);
			for (int l = 0; l < srcs.size(); l++) {
				const auto& src = srcs[l];
				for (size_t j = src.first; j < src.second; j++) {
					const simd_fixed32 x1 = *((simd_fixed32*) (&particle_x(j)));
					const simd_fixed32 y1 = *((simd_fixed32*) (&particle_y(j)));
					const simd_fixed32 z1 = *((simd_fixed32*) (&particle_z(j)));
					vec3<force_simd> fn, ff;
					const force_simd m = -force_simd::mask(std::min((size_t) FORCE_SIMD_SIZE, src.second - j));
					dx[0] = distance(x1, x0);
					dx[1] = distance(y1, y0);
					dx[2] = distance(z1, z0);
					const force_simd r2 = sqr(dx[0]) + sqr(dx[1]) + sqr(dx[2]);
					const force_simd rzero = r2 < force_simd(1.17549435082228750797e-37);
					const force_simd rinv = rsqrt(r2 + rzero);
					const force_simd rinv3 = sqr(rinv) * rinv;
					const force_simd pf = rinv;
					const force_simd pn = (force_simd(1.5) * hinv - force_simd(0.5) * r2 * hinv3);
					force_simd wn = r2 < h2;
					force_simd wf = force_simd(1) - wn;
					ff[0] = dx[0] * rinv3;
					ff[1] = dx[1] * rinv3;
					ff[2] = dx[2] * rinv3;
					fn[0] = dx[0] * hinv3;
					fn[1] = dx[1] * hinv3;
					fn[2] = dx[2] * hinv3;
					wn *= m;
					wf *= m;
					df.potential += pn * wn + pf * wf;
					for (int dim = 0; dim < NDIM; dim++) {
						df.force[dim] += fn[dim] * wn + ff[dim] * wf;
					}
				}
			}
			f.potential += reduce_sum(df.potential);
			f.potential *= GM;
			for (int dim = 0; dim < NDIM; dim++) {
				f.force[dim] += reduce_sum(f.force[dim]);
				f.force[dim] *= GM;
			}
			auto& vx = particle_vx(i);
			auto& vy = particle_vy(i);
			auto& vz = particle_vz(i);
			const float sgn = params.top ? 1.f : -1.f;
			if (params.ascending) {
				const float dt = 0.5f * params.t0 / (1 << params.rung);
				if (!params.first) {
					vx += sgn * f.force[0][i] * dt;
					vy += sgn * f.force[1][i] * dt;
					vz += sgn * f.force[2][i] * dt;
				}
			}
			if (params.descending) {
				const float g2 = sqr(f.force[0][i]) + sqr(f.force[1][i]) + sqr(f.force[2][i]);
				float dt = fminf(params.eta * sqrt(params.scale * h / sqrtf(g2)), params.t0);
				dt = fminf(params.max_dt, dt);
				int rung = particle_rung(i);
				rung = std::max(params.rung + int((int) ceilf(log2f(params.t0 / dt)) > params.rung), rung - 1);
				particle_rung(i) = rung;
				max_rung = std::max(rung, max_rung);
				dt = 0.5f * params.t0 / (1 << params.rung);
				vx += sgn * f.force[0][i] * dt;
				vy += sgn * f.force[1][i] * dt;
				vz += sgn * f.force[2][i] * dt;
			}
		}
	}
	return max_rung;
}

int kick(const kick_params& params) {
	using namespace sfmm;
	const int N = particle_cell_dim();
	int maxrung = 0;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N; k++) {
				const auto snk = particle_get_cell(i, j, k);
				vec3<double> C0;
				C0[0] = (i + 0.5) / N;
				C0[1] = (j + 0.5) / N;
				C0[2] = (k + 0.5) / N;
				const auto srcs = particle_cell_neighbors(i, j, k);
				const auto L = particle_get_expansion(i, j, k);
				const int rung = kick_cell(snk, C0, srcs, L, params);
				maxrung = std::max(maxrung, rung);
			}
		}
	}
	return maxrung;
}
