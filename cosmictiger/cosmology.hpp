/*
 * cosmology.hpp
 *
 *  Created on: Aug 6, 2023
 *      Author: dmarce1
 */

#ifndef COSMOLOGY_HPP_
#define COSMOLOGY_HPP_



double cosmos_conformal_time(double a0, double a1);
double cosmos_dadt(double a);
double cosmos_time(double a0, double a1);
void cosmos_update(double& adotdot, double& adot, double& a, double dt0);


#endif /* COSMOLOGY_HPP_ */
