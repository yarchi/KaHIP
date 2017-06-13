#pragma once

#include <type_traits>

#include "data_structure/parallel/hash_table.h"

namespace parallel {

template <typename _key_type, typename _value_type>
class hash_table_map {
public:
        static_assert(std::is_integral<_key_type>::value, "key_type shoud be integral");
        using key_type = _key_type;
        using value_type = _value_type;
        using map_type = hash_map<key_type, value_type>;

        explicit hash_table_map(uint64_t)
                :       map(get_max_size_to_fit_l1())
        {}

        inline bool contains(key_type key, value_type& value) {
                return map.contains(key, value);
        }

        inline void clear() {
                map.clear();
                map.reserve(get_max_size_to_fit_l1());
        }

        inline value_type& operator[](key_type key) {
                return map[key];
        }

        inline size_t ht_memory_size() const {
                return map.ht_memory_size();
        }

        inline typename map_type::Iterator begin() {
                return map.begin();
        }

        inline typename map_type::Iterator end() {
                return map.end();
        }

private:
        map_type map;

        static constexpr size_t size_factor = map_type::size_factor;

        static constexpr size_t get_max_size_to_fit_l1() {
                // (SizeFactor + 1.1) * max_size * sizeof(Element) + sizeof(Position) * max_size Bytes = 16 * 1024 Bytes,
                // where 16 * 1024 Bytes is half of L1 cache and (2 + 1.1) * max_size * 8 + 4 * max_size Bytes
                // is the size of a hash table. We calculate that max_size ~ 560.
                using element_type = typename map_type::Element;
                return round_up_to_previous_power_2( 32 * 1024 / (sizeof(element_type) * (size_factor)) );
        }
};

template <typename _key_type, typename _value_type>
class array_map {
public:
        static_assert(std::is_integral<_key_type>::value, "key_type shoud be integral");
        using key_type = _key_type;
        using value_type = _value_type;

        explicit array_map(uint64_t _max_size = 0)
                :       max_size(_max_size)
        {
                init(max_size);
        }

        inline void init(uint64_t _max_size) {
                map.assign(_max_size, non_initialized);
                accessed.reserve(128);
        }

        inline bool contains(key_type key, value_type& value) {
                value = map[key];
                return value != non_initialized;
        }

        inline void clear() {
                for (auto key : accessed) {
                        map[key] = non_initialized;
                }
                accessed.clear();
        }

        inline value_type& operator[](key_type key) {
                if (map[key] == non_initialized) {
                        accessed.push_back(key);
                }
                return map[key];
        }

        inline typename std::vector<value_type>::iterator begin() {
                return map.begin();
        }

        inline typename std::vector<value_type>::iterator end() {
                return map.end();
        }

private:
        static constexpr value_type non_initialized = std::numeric_limits<value_type>::max();

        std::vector<value_type> map;
        std::vector<key_type> accessed;
        uint64_t max_size;
};

template <typename _key_type, typename _value_type>
constexpr typename array_map<_key_type, _value_type>::value_type array_map<_key_type, _value_type>::non_initialized;

template <
        typename _key_type,
        typename _value_type,
        uint64_t l1_cache_size = 32 * 1024 * sizeof(char),
        uint64_t l2_cache_size = 256 * 1024 * sizeof(char)
>
class cache_aware_map_impl {
public:
        using key_type = _key_type;
        using value_type = _value_type;
        using small_map_type = hash_table_map<key_type, value_type>;
        using large_map_type = array_map<key_type, value_type>;

        static_assert(std::is_same<typename small_map_type::value_type, typename large_map_type::value_type>::value,
                      "value_type shoud be the same");

        static_assert(std::is_integral<key_type>::value, "key_type shoud be integral");

        cache_aware_map_impl(uint64_t _max_size)
                :       small_map(_max_size)
                ,       max_size(_max_size)
                ,       cur_container(cur_container_type::SMALL)
        {
                try_swap_containers();
        }

        inline bool contains(key_type key, value_type& value) {
                if (cur_container == cur_container_type::SMALL) {
                        return small_map.contains(key, value);
                } else {
                        return large_map.contains(key, value);
                }
        }

        inline void clear() {
                small_map.clear();
                large_map.clear();
        }

        inline value_type& operator[](key_type key) {
                if (cur_container == cur_container_type::SMALL) {
                        value_type& ref = small_map[key];

                        // check if growth happened and if yes then check
                        // hash table fits into L2 cache
                        if (!try_swap_containers())
                                return ref;
                }
                return large_map[key];
        }

        size_t memory_size() const {
                if (cur_container == cur_container_type::SMALL) {
                        return small_map.ht_memory_size();
                } else {
                        return large_map.size() * sizeof(value_type);
                }
        }
private:
        enum class cur_container_type {
                SMALL,
                LARGE
        };

        small_map_type small_map;
        large_map_type large_map;
        std::vector<key_type> accessed;
        uint64_t max_size;
        cur_container_type cur_container;

        inline bool try_swap_containers() {
                if (cur_container == cur_container_type::SMALL && small_map.ht_memory_size() > l2_cache_size) {
                        large_map.init(max_size);

                        for (const auto& rec : small_map) {
                                large_map[rec.first] = rec.second;
                        }

                        small_map.clear();
                        cur_container = cur_container_type::LARGE;
                        return true;
                }
                return false;
        }
};

template <typename key_type, typename value_type>
using cache_aware_map = cache_aware_map_impl<key_type, value_type>;

}
