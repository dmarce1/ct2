#include <hpx/hpx.hpp>
#include <cosmictiger/constants.hpp>
#include <cosmictiger/util.hpp>
#include <cosmictiger/options.hpp>
#include <cmath>

template<class T>
inline T sqr(T a) {
	return a * a;
}

#define SHOW( opt ) show(#opt, opts.opt)

options global_opts;

static void show(const char* name, bool opt) {
	printf("%-20s: %c\n", name, opt ? 'T' : 'F');
}

static void show(const char* name, int opt) {
	printf("%-20s: %i\n", name, opt);
}

static void show(const char* name, size_t opt) {
	printf("%-20s: %li\n", name, opt);
}

static void show(const char* name, float opt) {
	printf("%-20s: %e\n", name, opt);
}

static void show(const char* name, std::string opt) {
	printf("%-20s: %s\n", name, opt.c_str());
}

static void set_options(const options& opts) {
	global_opts = opts;
}

const options& get_options() {
	return global_opts;
}

bool process_options(int argc, char *argv[]) {
	options opts;
	namespace po = hpx::program_options;
	bool rc;
	po::options_description command_opts("options");

	command_opts.add_options()                                                                       //
	("help", "produce help message")                                                                       //
	("config_file", po::value < std::string > (&(opts.config_file))->default_value(""), "configuration file")                                                  //
	("hsoft", po::value<float>(&(opts.hsoft))->default_value(1.0 / 50.0), "dark matter softening ") //
	("bwidth", po::value<float>(&(opts.bwidth))->default_value(3.0), "width of direct window") //
	("omega_m", po::value<float>(&(opts.omega_m))->default_value(0.25), "") //
	("mass_res", po::value<float>(&(opts.mass_res))->default_value(1e9), "mass of a single particle in solar masses") //
	("hubble", po::value<float>(&(opts.hubble))->default_value(0.7), "") //
	("eta", po::value<float>(&(opts.eta))->default_value(0.15), "accuracy parameter") //
	("Theta", po::value<float>(&(opts.Theta))->default_value(1.0), "CMB temp normalized to 2.73 K") //
	("Ngrid", po::value < size_t > (&(opts.Ngrid))->default_value(32), "grid dimension") //
	("ppcell", po::value<int>(&(opts.ppcell))->default_value(256), "particles per cell") //
			;
	printf("Processing options\n");
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, command_opts), vm);
	po::notify(vm);
	if (vm.count("help")) {
		std::cout << command_opts << "\n";
		rc = false;
	} else {
		if (!opts.config_file.empty()) {
			std::ifstream cfg_fs { vm["config_file"].as<std::string>() };
			if (cfg_fs) {
				po::store(po::parse_config_file(cfg_fs, command_opts), vm);
				rc = true;
			} else {
				printf("Configuration file %s not found!\n", opts.config_file.c_str());
				return false;
			}
		} else {
			rc = true;
		}
	}

	if (rc) {
		po::notify(vm);
	}
	opts.GM = 1.0 / (opts.Ngrid * opts.Ngrid * opts.Ngrid);
	opts.nparts = opts.Ngrid * opts.Ngrid * opts.Ngrid;
	opts.code_to_g = constants::M0 * opts.mass_res;
	opts.code_to_cm = std::pow(opts.code_to_g * (8.0 * M_PI) * opts.nparts * constants::G / (3.0 * opts.omega_m * sqr(constants::H0 * opts.hubble)), 1.0 / 3.0);
	opts.code_to_s = opts.code_to_cm / constants::c;
	constexpr double Neff = 3.086;
	opts.omega_r = 32.0 * M_PI / 3.0 * constants::G * constants::sigma * (1 + Neff * (7. / 8.0) * std::pow(4. / 11., 4. / 3.)) * std::pow(constants::H0, -2) * std::pow(constants::c, -3) * std::pow(2.73 * opts.Theta, 4) * std::pow(opts.hubble, -2);
	opts.omega_lam = 1.0 - opts.omega_r - opts.omega_m;
	SHOW(bwidth);
	SHOW(config_file);
	SHOW(code_to_cm);
	SHOW(code_to_g);
	SHOW(code_to_s);
	SHOW(eta);
	SHOW(GM);
	SHOW(hsoft);
	SHOW(hubble);
	SHOW(omega_lam);
	SHOW(omega_m);
	SHOW(omega_r);
	SHOW(mass_res);
	SHOW(Ngrid);
	SHOW(nparts);
	SHOW(ppcell);
	SHOW(Theta);
	set_options(opts);
	return rc;
}

