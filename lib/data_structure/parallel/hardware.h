#pragma once

#include "tools/macros_assertions.h"

#include <cstdint>

namespace parallel {
static constexpr uint64_t g_l1_cache_size = 32 * 1024;
static constexpr uint64_t g_l2_cache_size = 256 * 1024;
static constexpr uint64_t g_l3_cache_size = 20480 * 1024;

static constexpr uint32_t g_threads_per_socket = 8;

static size_t get_mem_for_thread(uint32_t proc_id, uint32_t num_threads) {
        ALWAYS_ASSERT(proc_id < num_threads);
        uint32_t proc_in_full_socket = (num_threads / g_threads_per_socket) * g_threads_per_socket;
        uint32_t proc_per_socket;
        if (proc_id >= proc_in_full_socket) {
                proc_per_socket = num_threads - proc_in_full_socket;
        } else {
                proc_per_socket = g_threads_per_socket;
        }

        // two tables of size 1024 with uint32_t elements
        size_t tabular_hash_mem = 8 * 1024;
        ALWAYS_ASSERT(g_l3_cache_size > tabular_hash_mem * num_threads);
        return std::max(g_l3_cache_size / proc_per_socket, g_l2_cache_size - tabular_hash_mem);
}

}