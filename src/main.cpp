#include <hpx/hpx.hpp>
#include <hpx/hpx_init.hpp>
#include <cosmictiger/options.hpp>
#include <cosmictiger/particle.hpp>
#include <cosmictiger/kick.hpp>
#include <cosmictiger/drift.hpp>

int kick_step(int& minrung, int max_rung, double scale, double adot, double tau, double t0, int minrung0, bool nocrop) {
	using namespace sfmm;
	const size_t Ngrid = get_options().Ngrid;
	const size_t nparts = Ngrid * Ngrid * Ngrid;
	std::vector<int> levels(std::max(max_rung - minrung, 0) + 1);
	int k = 0;
	for (int i = max_rung; i > minrung; i--) {
		levels[k++] = i;
	}
	levels[k++] = minrung;
	bool ascending = true;
	bool clip_top = false;
	bool top;
	for (int li = 0; li < levels.size(); li++) {
		if (levels[li] == minrung) {
			ascending = false;
			top = true;
		} else {
			top = false;
		}
		if (!ascending && !top) {
			particle_push_rungs();
		}
		if (!ascending) {
			particle_sort_by_rung(levels[li]);
		}
		auto rng = particle_current_range();
		if (top && minrung == minrung0) {
			auto counts = particle_rung_counts();
			if (counts.size() > minrung0 + 1) {
				const auto total = nparts;
				size_t fast = 0;
				size_t slow = counts[minrung0];
				for (int i = minrung0 + 1; i < counts.size(); i++) {
					fast += counts[i];
				}
				if (3 * fast > slow) {
					clip_top = !nocrop;
				}
			}
		}
		particle_to_grid();
		kick_params kparams;
		if (clip_top && levels[li] == minrung0 + 1) {
			kparams.top = true;
		} else {
			kparams.top = top;
		}
		kparams.ascending = ascending || top;
		if (clip_top && top) {
			kparams.descending = false;
		} else {
			kparams.descending = !ascending || top;
		}
		kparams.eta = get_options().eta;
		kparams.scale = scale;
		kparams.first = tau == 0.0;
		kparams.rung = levels[li];
		kparams.t0 = t0;
		kparams.max_dt = 0.01f / adot;
		int max_rung = kick(kparams);
		if (clip_top && top) {
			particle_set_minrung(minrung0 + 1);
		}
		if (clip_top && top) {
			levels.push_back(minrung0 + 1);
		} else if (!ascending || top) {
			if (max_rung > levels[li]) {
				levels.push_back(levels[li] + 1);
			}
		}
		if ((!ascending || top) && !(clip_top && top)) {
			const double dt = t0 / (1 << levels[li]);
			drift(dt, scale);
		}
		if (ascending && !top) {
			particle_pop_rungs();
		}
	}
	if (clip_top) {
		minrung++;
	}
	return max_rung;

}

int hpx_main(int argc, char *argv[]) {
	process_options(argc, argv);
	return hpx::finalize();
}

int main(int argc, char *argv[]) {
	std::vector<std::string> cfg = {"hpx.commandline.allow_unknown=1"};
	return hpx::init(argc, argv);
}
