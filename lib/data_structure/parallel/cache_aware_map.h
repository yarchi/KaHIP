#pragma once

#include <type_traits>

#include "data_structure/parallel/hash_table.h"

namespace parallel {

template <
        typename small_map_type,
        typename large_map_type,
        uint64_t l1_cache_size = 32 * 1024 * sizeof(char),
        uint64_t l2_cache_size = 256 * 1024 * sizeof(char)
>
class cache_aware_map_impl {
public:
        static_assert(std::is_same<typename small_map_type::mapped_type, typename large_map_type::value_type>::value, "value_type shoud be the same");

        using key_type = typename small_map_type::key_type;
        using value_type = typename small_map_type::mapped_type;

        static_assert(std::is_integral<key_type>::value, "key_type shoud be integral");

        cache_aware_map_impl(uint64_t _max_size)
#ifndef ONLY_ARRAY
                :       small_map(get_max_size_to_fit_l1())
                ,       max_size(_max_size)
                ,       cur_container(cur_container_type::SMALL)
#else
                :       max_size(_max_size)
                ,       cur_container(cur_container_type::LARGE)
#endif
        {
                try_swap_containers();
        }

        inline bool contains(key_type key, value_type& value) {
#ifndef ONLY_ARRAY
                if (cur_container == cur_container_type::SMALL) {
                        return small_map.contains(key, value);
                } else {
                        value = large_map[key];
                        return value != non_initialized;
                }
#else
                value = large_map[key];
                return value != non_initialized;
#endif
        }

        inline void clear() {
#ifndef ONLY_ARRAY
                small_map.clear();
                small_map.reserve(get_max_size_to_fit_l1());
                cur_container = cur_container_type::SMALL;
#endif
                for (auto key : accessed) {
                        large_map[key] = non_initialized;
                }
                accessed.clear();
        }

        inline value_type& operator[](key_type key) {
#ifndef ONLY_ARRAY
                if (cur_container == cur_container_type::SMALL) {
                        value_type& ref = small_map[key];

                        // check if growth happened and if yes then check
                        // hash table fits into L2 cache
                        if (!try_swap_containers())
                                return ref;
                }
#endif
                if (large_map[key] == non_initialized) {
                        accessed.push_back(key);
                }
                return large_map[key];
        }

        size_t memory_size() const {
                if (cur_container == cur_container_type::SMALL) {
                        return small_map.ht_memory_size();
                } else {
                        return max_size * sizeof(value_type);
                }
        }
private:
        enum class cur_container_type {
                SMALL,
                LARGE
        };

        static constexpr size_t size_factor = small_map_type::size_factor;
        static constexpr value_type non_initialized = std::numeric_limits<value_type>::max();

        small_map_type small_map;
        large_map_type large_map;
        std::vector<key_type> accessed;
        uint64_t max_size;
        cur_container_type cur_container;

        inline bool try_swap_containers() {
                if (cur_container == cur_container_type::SMALL && small_map.ht_memory_size() > l2_cache_size) {
                        large_map.assign(max_size, non_initialized);
                        accessed.reserve(128);

                        for (const auto& rec : small_map) {
                                large_map[rec.first] = rec.second;
                        }

                        small_map.clear();
                        cur_container = cur_container_type::LARGE;
                        return true;
                }
                return false;
        }

        static constexpr size_t get_max_size_to_fit_l1() {
                // (SizeFactor + 1.1) * max_size * sizeof(Element) + sizeof(Position) * max_size Bytes = 16 * 1024 Bytes,
                // where 16 * 1024 Bytes is half of L1 cache and (2 + 1.1) * max_size * 8 + 4 * max_size Bytes
                // is the size of a hash table. We calculate that max_size ~ 560.
                using position_type = typename small_map_type::Position;
                using element_type = typename small_map_type::Element;

                return round_up_to_previous_power_2(16 * 1024 / (sizeof(element_type) * (size_factor + 1.1)));
        }
};

template <
        typename small_map_type,
        typename large_map_type,
        uint64_t l1_cache_size,
        uint64_t l2_cache_size
>
constexpr typename cache_aware_map_impl<small_map_type, large_map_type, l1_cache_size, l2_cache_size>::value_type
        cache_aware_map_impl<small_map_type, large_map_type, l1_cache_size, l2_cache_size>::non_initialized;

template <typename key_type, typename value_type>
using cache_aware_map = cache_aware_map_impl<hash_map<key_type, value_type>, std::vector<value_type>>;

}
