/*
 * util.hpp
 *
 *  Created on: Jul 26, 2023
 *      Author: dmarce1
 */

#ifndef UTIL_HPP_
#define UTIL_HPP_

#include <cstdint>

using time_type = std::uint64_t;


inline int round_up(int i, int r) {
	return r * ((i - 1) / r + 1);
}

inline time_type inc(time_type t, int max_rung) {
	if (max_rung > 0) {
		t += (time_type(1) << time_type(64 - max_rung));
	}
	return t;
}

inline int min_rung(time_type t) {
	int min_rung = 64;
	while (((t & 1) == 0) && (min_rung != 0)) {
		min_rung--;
		t >>= 1;
	}
	return min_rung;
}


double rand1();


#endif /* UTIL_HPP_ */
