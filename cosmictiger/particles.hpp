/*
 * particle.hpp
 *
 *  Created on: Jul 26, 2023
 *      Author: dmarce1
 */

#ifndef PARTICLES_HPP_
#define PARTICLES_HPP_

#include <cstdint>

using fixed32 = std::uint32_t;

#ifdef PARTICLES_CPP
#define PARTICLES_EXTERN
#else
#define PARTICLES_EXTERN extern
#endif

PARTICLES_EXTERN fixed32* particle_x_ptr;
PARTICLES_EXTERN fixed32* particle_y_ptr;
PARTICLES_EXTERN fixed32* particle_z_ptr;
PARTICLES_EXTERN float* particle_vx_ptr;
PARTICLES_EXTERN float* particle_vy_ptr;
PARTICLES_EXTERN float* particle_vz_ptr;
PARTICLES_EXTERN char* particle_rung_ptr;

fixed32& particle_x(int i) {
	return particle_x_ptr[i];
}

fixed32& particle_y(int i) {
	return particle_y_ptr[i];
}

fixed32& particle_z(int i) {
	return particle_z_ptr[i];
}

float& particle_vx(int i) {
	return particle_vx_ptr[i];
}

float& particle_vy(int i) {
	return particle_vy_ptr[i];
}

float& particle_vz(int i) {
	return particle_vz_ptr[i];
}

char& particle_rung(int i) {
	return particle_rung_ptr[i];
}

void particle_set_x(fixed32 x, int n) {
	particle_x_ptr[n] = x;
}

void particle_set_y(fixed32 x, int n) {
	particle_y_ptr[n] = x;
}

void particle_set_z(fixed32 x, int n) {
	particle_z_ptr[n] = x;
}


#endif /* PARTICLES_HPP_ */
