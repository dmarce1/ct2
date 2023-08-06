#include <cosmictiger/cosmictiger.hpp>
#include <cosmictiger/particle.hpp>

void drift(float dt, float scale) {
	using namespace sfmm;
	const float lambda = dt / scale;
	for (size_t i = 0; i < particle_count(); i += DRIFT_SIMD_SIZE) {
		const simd_f32 vx = *((simd_f32*) (&particle_vx(i)));
		const simd_f32 vy = *((simd_f32*) (&particle_vy(i)));
		const simd_f32 vz = *((simd_f32*) (&particle_vz(i)));
		simd_fixed32& x = *((simd_fixed32*) (&particle_x(i)));
		simd_fixed32& y = *((simd_fixed32*) (&particle_y(i)));
		simd_fixed32& z = *((simd_fixed32*) (&particle_z(i)));
		x = x + simd_fixed32(vx * lambda);
		y = y + simd_fixed32(vy * lambda);
		z = z + simd_fixed32(vz * lambda);
	}
}
