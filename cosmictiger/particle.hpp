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
#include <cosmictiger/cosmictiger.hpp>


#ifdef PARTICLES_CPP
#define PARTICLES_EXTERN
#else
#define PARTICLES_EXTERN extern
#endif

PARTICLES_EXTERN pos_fixed* particle_x_ptr;
PARTICLES_EXTERN pos_fixed* particle_y_ptr;
PARTICLES_EXTERN pos_fixed* particle_z_ptr;
PARTICLES_EXTERN force_float* particle_vx_ptr;
PARTICLES_EXTERN force_float* particle_vy_ptr;
PARTICLES_EXTERN force_float* particle_vz_ptr;
PARTICLES_EXTERN char* particle_rung_ptr;

inline pos_fixed& particle_x(int i) {
	return particle_x_ptr[i];
}

inline pos_fixed& particle_y(int i) {
	return particle_y_ptr[i];
}

inline pos_fixed& particle_z(int i) {
	return particle_z_ptr[i];
}

inline force_float& particle_vx(int i) {
	return particle_vx_ptr[i];
}

inline force_float& particle_vy(int i) {
	return particle_vy_ptr[i];
}

inline force_float& particle_vz(int i) {
	return particle_vz_ptr[i];
}

inline char& particle_rung(int i) {
	return particle_rung_ptr[i];
}

size_t particle_count();
void particles_init();
void particle_sort_by_rung(int minrung);
std::vector<size_t> particle_rung_counts();
void particle_to_grid();
std::pair<size_t, size_t> particle_current_range();
void particle_to_grid();
void particle_pop_rungs();
void particle_push_rungs();
void particle_set_minrung(int minrung);
std::vector<std::pair<size_t, size_t>> particle_cell_neighbors(int i, int j, int k);
int particle_cell_dim();
std::pair<size_t, size_t> particle_get_cell(int i, int j, int k);
sfmm::expansion<fft_float, PORDER> particle_get_expansion(int i, int j, int k);

#endif /* PARTICLES_HPP_ */
