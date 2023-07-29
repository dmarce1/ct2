#define PARTICLES_CPP

#include <stdlib.h>
#include "cosmictiger/cosmictiger.hpp"
#include "cosmictiger/particles.hpp"

size_t particle_count = 0;


void particles_allocate(size_t count) {
 	particle_count = count;
	const size_t nbytes = count * (NDIM * (sizeof(fixed32) + sizeof(float) + sizeof(char)));
	char* ptr = new char[nbytes];
	particles_x_ptr = (fixed32*) ptr;
	particles_y_ptr = (fixed32*) (ptr + NDIM * count);
	particles_z_ptr = (fixed32*) (ptr + 2 * NDIM * count);
	particles_vx_ptr = (float*) (ptr + 3 * NDIM * count);
	particles_vy_ptr = (float*) (ptr + 4 * NDIM * count);
	particles_vz_ptr = (float*) (ptr + 5 * NDIM * count);
	particles_rung_ptr = (char*) (ptr + 5 * NDIM * count);
}
