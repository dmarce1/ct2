/*
 * particles.hpp
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

PARTICLES_EXTERN fixed32* particles_x_ptr;
PARTICLES_EXTERN fixed32* particles_y_ptr;
PARTICLES_EXTERN fixed32* particles_z_ptr;
PARTICLES_EXTERN float* particles_vx_ptr;
PARTICLES_EXTERN float* particles_vy_ptr;
PARTICLES_EXTERN float* particles_vz_ptr;
PARTICLES_EXTERN char* particles_rung_ptr;

fixed32 particles_x(int i) {
	return particles_x_ptr[i];
}

fixed32 particles_y(int i) {
	return particles_y_ptr[i];
}

fixed32 particles_z(int i) {
	return particles_z_ptr[i];
}

float particles_vx(int i) {
	return particles_vx_ptr[i];
}

float particles_vy(int i) {
	return particles_vy_ptr[i];
}

float particles_vz(int i) {
	return particles_vz_ptr[i];
}

int particles_rung(int i) {
	return particles_rung_ptr[i];
}

void particles_set_x(fixed32 x, int n) {
	particles_x_ptr[n] = x;
}

void particles_set_y(fixed32 x, int n) {
	particles_y_ptr[n] = x;
}

void particles_set_z(fixed32 x, int n) {
	particles_z_ptr[n] = x;
}


#endif /* PARTICLES_HPP_ */
