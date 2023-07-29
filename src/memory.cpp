#include <stdio.h>
#include <stdlib.h>
#include "cosmictiger/util.hpp"

void* operator new(std::size_t n) {
	void* memptr;
	if (posix_memalign(&memptr, 32, round_up(n, 32)) != 0) {
		printf("operator new failed!\n");
		abort();
	}
	return memptr;
}

void operator delete(void * p) {
	free(p);
}

void *operator new[](std::size_t n) {
	void* memptr;
	if (posix_memalign(&memptr, 32, round_up(n, 32)) != 0) {
		printf("operator new[] failed!\n");
		abort();
	}
	return memptr;
}

void operator delete[](void *p) {
	free(p);
}
