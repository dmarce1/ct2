#include <cosmictiger/cosmictiger.hpp>
#include <cosmictiger/options.hpp>
#include <cosmictiger/particle.hpp>
#include <sfmm.hpp>
#include <cmath>

struct kick_params {
	bool ascending;
	bool descending;
	bool sign;
	bool top;
	bool first;
	int rung;
	float t0;
	float eta;
	float max_dt;
	float scale;
};

int kick(const std::pair<size_t, size_t>& snk, const sfmm::vec3<float> C0, const std::vector<std::pair<size_t, size_t>>& srcs, const sfmm::expansion<float, PORDER>& L0, const kick_params& params) {
	using namespace sfmm;
	const float h = get_options().hsoft;
	const float GM = get_options().GM;
	const float h2 = h * h;
	const float hinv = 1.0f / h;
	const float hinv3 = sqr(hinv) * hinv;
	int max_rung = 0;
	expansion < simd_f32, PORDER > L;
	vec3 < simd_fixed32 > C;
	for (int dim = 0; dim < NDIM; dim++) {
		C[dim] = fixed32(C0[dim]);
	}
	for (int n = 0; n < expansion<float, PORDER>::size(); n++) {
		L[n] = L0[n];
	}
	for (size_t n = snk.first; n < snk.second; n += SIMD_SIZE) {
		force_type < simd_f32 > f;
		f.init();
		vec3 < simd_f32 > dx;
		const simd_fixed32 x = *((simd_fixed32*) (&particle_x(n)));
		const simd_fixed32 y = *((simd_fixed32*) (&particle_y(n)));
		const simd_fixed32 z = *((simd_fixed32*) (&particle_z(n)));
		dx[0] = distance(x, C[0]);
		dx[1] = distance(y, C[1]);
		dx[2] = distance(z, C[2]);
		L2P(f, L, dx);
		for (size_t i = n; i < std::min(n + SIMD_SIZE, snk.second); i++) {
			force_type < simd_f32 > df;
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
					vec3<simd_f32> fn, ff;
					simd_f32 m = simd_f32::mask(std::min((size_t) SIMD_SIZE, src.second - j));
					dx[0] = distance(x1, x0);
					dx[1] = distance(y1, y0);
					dx[2] = distance(z1, z0);
					const simd_f32 r2 = sqr(dx[0]) + sqr(dx[1]) + sqr(dx[2]);
					const simd_f32 rzero = r2 < simd_f32(1.17549435082228750797e-37);
					const simd_f32 rinv = rsqrt(r2 + rzero);
					const simd_f32 rinv3 = sqr(rinv) * rinv;
					const simd_f32 pf = rinv;
					const simd_f32 pn = (simd_f32(1.5) * hinv - simd_f32(0.5) * r2 * hinv3);
					simd_f32 wn = r2 < h2;
					simd_f32 wf = simd_f32(1) - wn;
					ff[0] = dx[0] * rinv3;
					ff[1] = dx[1] * rinv3;
					ff[2] = dx[2] * rinv3;
					fn[0] = dx[0] * hinv3;
					fn[1] = dx[1] * hinv3;
					fn[2] = dx[2] * hinv3;
					m = -m;
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
			float dt;
			const float sgn = params.sign * (params.top ? 1.f : -1.f);
			if (params.ascending) {
				dt = 0.5f * params.t0 / (1 << params.rung);
				if (!params.first) {
					vx += sgn * f.force[0][i] * dt;
					vy += sgn * f.force[1][i] * dt;
					vz += sgn * f.force[2][i] * dt;
				}
			}
			if (params.descending) {
				const float g2 = sqr(f.force[0][i]) + sqr(f.force[1][i]) + sqr(f.force[2][i]);
				dt = fminf(params.eta * sqrt(params.scale * h / sqrtf(g2)), params.t0);
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
