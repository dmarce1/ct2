#define PARTICLES_CPP

#include <array>
#include <cstdlib>
#include <limits>
#include <vector>
#include <hpx/async.hpp>
#include "cosmictiger/cosmictiger.hpp"
#include "cosmictiger/particles.hpp"

size_t particle_count = 0;
size_t particle_end;
void* particle_ptr;

std::vector<std::pair<int, int>> cells;
int Ncell;

void particles_allocate(size_t count) {
	particle_count = count;
	const size_t nbytes = count * (NDIM * (sizeof(fixed32) + sizeof(float) + sizeof(char)));
	char* ptr = new char[nbytes];
	particle_x_ptr = (fixed32*) ptr;
	particle_y_ptr = (fixed32*) (ptr + NDIM * count);
	particle_z_ptr = (fixed32*) (ptr + 2 * NDIM * count);
	particle_vx_ptr = (float*) (ptr + 3 * NDIM * count);
	particle_vy_ptr = (float*) (ptr + 4 * NDIM * count);
	particle_vz_ptr = (float*) (ptr + 5 * NDIM * count);
	particle_rung_ptr = (char*) (ptr + 5 * NDIM * count);
	particle_ptr = ptr;
}

