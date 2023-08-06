/*
 * options.hpp
 *
 *  Created on: Aug 4, 2023
 *      Author: dmarce1
 */

#ifndef OPTIONS_HPP_
#define OPTIONS_HPP_

#include <cstdlib>
#include <string>

struct options {
	float bwidth;
	std::string config_file;
	float code_to_cm;
	float code_to_g;
	float code_to_s;
	float eta;
	float GM;
	float hsoft;
	float hubble;
	float omega_lam;
	float omega_m;
	float omega_r;
	float mass_res;
	size_t Ngrid;
	size_t nparts;
	int ppcell;
	float Theta;
};

const options& get_options();
bool process_options(int argc, char *argv[]);



#endif /* OPTIONS_HPP_ */
