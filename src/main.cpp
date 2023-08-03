#include <hpx/hpx.hpp>
#include <hpx/hpx_init.hpp>
#include <simd.hpp>
#include <sfmm.hpp>

#include "cosmictiger/cosmictiger.hpp"

int hpx_main(int argc, char *argv[]) {
	return hpx::finalize();
}

int main(int argc, char *argv[]) {
	std::vector<std::string> cfg = {"hpx.commandline.allow_unknown=1"};
	return hpx::init(argc, argv);
}
