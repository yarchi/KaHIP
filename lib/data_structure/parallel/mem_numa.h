#include <cstdint>

#ifdef __gnu_linux__
#include <numa.h>
#endif

#include <string>
#include "tools/macros_assertions.h"

namespace parallel {

void numa_memory_interleave(uint32_t num_threads, uint32_t threads_per_socket) {
#ifdef __gnu_linux__
        if (numa_available() < 0) {
		printf("No NUMA support available on this system.\n");
		exit(1);
	}

	uint32_t num_sockets = std::ceil((num_threads + 0.0) / threads_per_socket);
	ALWAYS_ASSERT(num_sockets > 0);
	std::string mask_string = "0-" + std::to_string(num_sockets - 1);
	auto* mask = numa_parse_nodestring(mask_string.c_str());
        ALWAYS_ASSERT(mask != 0);

        numa_set_interleave_mask(mask);
#endif
}
}
