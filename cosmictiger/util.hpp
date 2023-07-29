/*
 * util.hpp
 *
 *  Created on: Jul 26, 2023
 *      Author: dmarce1
 */

#ifndef UTIL_HPP_
#define UTIL_HPP_

inline int round_up(int i, int r) {
	return r * ((i - 1) / r + 1);
}

#endif /* UTIL_HPP_ */
