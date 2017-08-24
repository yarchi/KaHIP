#pragma once

#include <cstdint>

namespace parallel {
static constexpr uint64_t g_l1_cache_size = 32 * 1024 * sizeof(char);
static constexpr uint64_t g_half_l2_cache_size = 256 * 1024 * sizeof(char) / 2;
static constexpr uint64_t g_l2_cache_size = 256 * 1024 * sizeof(char);
static constexpr uint64_t g_l3_cache_size = 20480 * 1024 * sizeof(char);

static constexpr uint32_t g_threads_per_socket = 8;

}