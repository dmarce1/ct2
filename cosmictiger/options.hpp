/*
 * options.hpp
 *
 *  Created on: Aug 4, 2023
 *      Author: dmarce1
 */

#ifndef OPTIONS_HPP_
#define OPTIONS_HPP_

struct options {
	std::string config_file;
	float hsoft;
	float bwidth;
	float GM;
};

const options& get_options();
bool process_options(int argc, char *argv[]);



#endif /* OPTIONS_HPP_ */
