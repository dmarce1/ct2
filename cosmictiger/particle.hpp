/*
 * particle.hpp
 *
 *  Created on: Jul 26, 2023
 *      Author: dmarce1
 */

#ifndef PARTICLES_HPP_
#define PARTICLES_HPP_

#include <cstdint>
#include <sfmm.hpp>

using fixed32 = sfmm::fixed32;

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

inline fixed32& particle_x(int i) {
	return particle_x_ptr[i];
}

inline fixed32& particle_y(int i) {
	return particle_y_ptr[i];
}

inline fixed32& particle_z(int i) {
	return particle_z_ptr[i];
}

inline float& particle_vx(int i) {
	return particle_vx_ptr[i];
}

inline float& particle_vy(int i) {
	return particle_vy_ptr[i];
}

inline float& particle_vz(int i) {
	return particle_vz_ptr[i];
}

inline char& particle_rung(int i) {
	return particle_rung_ptr[i];
}


#endif /* PARTICLES_HPP_ */
