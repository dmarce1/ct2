#include <hpx/hpx.hpp>
#include <cosmictiger/options.hpp>

#define SHOW( opt ) show(#opt, opts.opt)

options global_opts;

static void show(const char* name, bool opt) {
	printf("%-20s: %c\n", name, opt ? 'T' : 'F');
}

static void show(const char* name, int opt) {
	printf("%-20s: %i\n", name, opt);
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
	("GM", po::value<float>(&(opts.GM))->default_value(1.0), "gravitational constant") //
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
	SHOW(bwidth);
	SHOW(hsoft);
	set_options(opts);
	return rc;
}

