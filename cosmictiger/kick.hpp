#pragma once


#include <cosmictiger/cosmictiger.hpp>
#include <cosmictiger/options.hpp>
#include <cosmictiger/particle.hpp>
#include <sfmm.hpp>
#include <cmath>

struct kick_params {
	bool ascending;
	bool descending;
	bool top;
	bool first;
	int rung;
	float t0;
	float eta;
	float max_dt;
	float scale;
};


int kick_step(int& minrung, int max_rung, double scale, double adot, double tau, double t0, int minrung0, bool nocrop);
int kick(const kick_params& params);

