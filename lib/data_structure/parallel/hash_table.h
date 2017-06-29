#pragma once

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#define XXH_PRIVATE_API
#include "data_structure/parallel/lib/xxhash.h"
#include "tools/macros_assertions.h"

namespace parallel {
template <typename T>
struct simple_hash {
        using hash_type = uint64_t;
        inline hash_type operator()(const T& x) const {
                return x;
        }
};

template <typename T>
struct xxhash {
        static const size_t significant_digits = 64;
        using hash_type = uint64_t;
        inline hash_type operator()(const T& x) const {
                return XXH64(&x, sizeof(x), 0);
        }
};

template <typename T>
struct xxhash<std::pair<T, T>> {
        using hash_type = uint64_t;

        xxhash<T> h;

        uint64_t operator() (const std::pair<T, T>& x) const {
                return h(x.first) ^ h(x.second);
        }
};

template <typename Key>
class MurmurHash {
public:
        static const size_t significant_digits = 64;
        using hash_type = std::uint64_t;

        explicit MurmurHash(uint32_t seed = 0)
                :   _seed(seed)
        {}

        void reset(uint32_t seed) {
                _seed = seed;
        }

        inline hash_type operator() (const Key& key) const {
                return hash(reinterpret_cast<const void*>(&key), sizeof(key), _seed);
        }

private:
        uint32_t _seed;

        inline hash_type hash(const void* key, uint32_t len, uint32_t seed) const {
                const uint64_t m = 0xc6a4a7935bd1e995;
                const int r = 47;

                uint64_t h = seed ^ (len * m);

                const uint64_t * data = (const uint64_t *)key;
                const uint64_t * end = data + (len/8);

                while(data != end)
                {
                        uint64_t k = *data++;

                        k *= m;
                        k ^= k >> r;
                        k *= m;

                        h ^= k;
                        h *= m;
                }

                const unsigned char * data2 = (const unsigned char*)data;

                switch(len & 7)
                {
                        case 7: h ^= uint64_t(data2[6]) << 48;
                        case 6: h ^= uint64_t(data2[5]) << 40;
                        case 5: h ^= uint64_t(data2[4]) << 32;
                        case 4: h ^= uint64_t(data2[3]) << 24;
                        case 3: h ^= uint64_t(data2[2]) << 16;
                        case 2: h ^= uint64_t(data2[1]) << 8;
                        case 1: h ^= uint64_t(data2[0]);
                                h *= m;
                };

                h ^= h >> r;
                h *= m;
                h ^= h >> r;

                return h;
        }
};

template <typename T>
struct MurmurHash<std::pair<T, T>> {
        using hash_type = uint64_t;

        MurmurHash<T> h;

        uint64_t operator() (const std::pair<T, T>& x) const {
                return h(x.first) ^ h(x.second);
        }
};

constexpr static uint32_t round_up_to_next_power_2(uint32_t v) {
        v--;
        v |= v >> 1;
        v |= v >> 2;
        v |= v >> 4;
        v |= v >> 8;
        v |= v >> 16;
        v++;
        return v;
}

constexpr static uint32_t round_up_to_previous_power_2(uint32_t x) {
        x = x | (x >> 1);
        x = x | (x >> 2);
        x = x | (x >> 4);
        x = x | (x >> 8);
        x = x | (x >> 16);
        return x - (x >> 1);
}

constexpr static bool is_power_2(uint64_t x) {
        return ((x != 0) && !(x & (x - 1)));
}

template <typename hash_table_type>
static constexpr size_t get_max_size_to_fit_l1() {
        // (SizeFactor + 1.1) * max_size * sizeof(Element) Bytes = 16 * 1024 Bytes,
        // where 16 * 1024 Bytes is half of L1 cache and (2 + 1.1) * max_size * 8
        // is the size of a hash table. We calculate that max_size ~ 512.
        return round_up_to_previous_power_2(16 * 1024 / (sizeof(typename hash_table_type::Element) * (hash_table_type::size_factor + 1.1)));
}

template <typename HashTable>
class HashTableIterator {
private:
        using HashMap = HashTable;
        using Element = typename HashMap::Element;
        using Position = typename HashMap::Position;

public:
        HashTableIterator(const HashMap& hm, const Position offset) :
                _hm(hm),
                _offset(offset)
        {}

        const Element& operator* () {
                return _hm._ht[_hm._poses[_offset]];
        }

        HashTableIterator& operator++ () {
                ++_offset;
                return *this;
        }

        bool operator== (const HashTableIterator& it) {
                return _offset == it._offset;
        }

        bool operator!= (const HashTableIterator& it) {
                return !(*this == it);
        }

private:
        const HashMap& _hm;
        Position _offset;
};

#undef PERFORMANCE_STATISTICS

template <typename Key, typename Value, typename Hash = xxhash<Key>,
        bool TGrowable = false, bool Cache = true, size_t SizeFactor = 2>
class HashMap {
public:
        using Element = std::pair<Key, Value>;
        using key_type = Key;
        using mapped_type = Value;
        using hash_type = typename Hash::hash_type;
        using Position = uint32_t;

private:
        using TSelf = HashMap<Key, Value, Hash, TGrowable, Cache, SizeFactor>;

        static constexpr hash_type max_hash_value = std::numeric_limits<hash_type>::max();

#ifdef PERFORMANCE_STATISTICS
        #pragma message("PERFORMANCE_STATISTICS FOR HASH TABLE IS ON")
        static size_t num_access;
        static size_t num_contain;
        static size_t overall_max_size;
        static size_t num_probes;
        static size_t num_find_pos;
#endif
public:
        using Iterator = HashTableIterator<TSelf>;

        friend Iterator;

        static constexpr size_t size_factor = SizeFactor;

        explicit HashMap(const uint64_t max_size = 0) :
                _empty_element(std::numeric_limits<Key>::max()),
                _max_size(std::max(round_up_to_next_power_2(max_size), 16u)),
                _ht_size(_max_size * SizeFactor),
                //_ht(_ht_size + _max_size * 1.1, std::make_pair(_empty_element, Value())),
                _ht(_ht_size, std::make_pair(_empty_element, Value())),
                _poses(),
                _hash(),
                _last_key(_empty_element),
                _last_position(0) {
                ALWAYS_ASSERT(is_power_2(_ht_size));
                _poses.reserve(_max_size);
        }

        inline constexpr size_t ht_memory_size() const {
                return _ht.size() * sizeof(Element);
        }

        HashMap(const TSelf&) = default;
        HashMap(TSelf&&) = default;

        TSelf& operator=(TSelf& other) = default;
        TSelf& operator=(TSelf&& other) = default;

#ifdef PERFORMANCE_STATISTICS
        ~HashMap() {
                overall_max_size = std::max<uint32_t>(overall_max_size, size());
        }
#endif

        Iterator begin() const {
                return Iterator(*this, 0);
        }

        Iterator end() const {
                return Iterator(*this, size());
        }

        void reserve(const uint32_t max_size) {
                _ht_size = max_size * SizeFactor;
                ALWAYS_ASSERT(is_power_2(_ht_size));
                _max_size = max_size;
                //_ht.resize(_ht_size + _max_size * 1.1, std::make_pair(_empty_element, Value()));
                _ht.resize(_ht_size, std::make_pair(_empty_element, Value()));
                _poses.reserve(_max_size);
        }

        inline uint64_t size() const {
                return _poses.size();
        }

        inline bool empty() const {
                return size() == 0;
        }

        inline Value& operator[](const Key& key) {
#ifdef PERFORMANCE_STATISTICS
                ++num_access;
#endif
                if (TGrowable && size() == _max_size)
                        resize();

                const Position pos = findPosition(key);
                if (_ht[pos].first == _empty_element) {
                        _ht[pos].first = key;
                        _ht[pos].second = Value();
                        _poses.push_back(pos);
                }

                return _ht[pos].second;
        }

        inline bool contains(const Key& key, Value& value) {
#ifdef PERFORMANCE_STATISTICS
                ++num_contain;
#endif
                size_t pos = findPosition(key);
                Element& elem = _ht[pos];
                if (elem.first != _empty_element) {
                        value = elem.second;
                        return true;
                } else {
                        return false;
                }
        }

        inline bool contains(const Key& key) {
#ifdef PERFORMANCE_STATISTICS
                ++num_contain;
#endif
                return _ht[findPosition(key)].first != _empty_element;
        }

        inline void insert(const Element& elem) {
                insertImpl(elem.first, elem.second);
        }

        inline void insert(const Key& key, const Value& value) {
                insertImpl(key, value);
        }

        inline void clear() {
#ifdef PERFORMANCE_STATISTICS
                overall_max_size = std::max<uint32_t>(overall_max_size, size());
#endif
                for (auto pos : _poses) {
                        _ht[pos].first = _empty_element;
                }

                _poses.clear();

                _last_key = _empty_element;
                _last_position = 0;
        }

        inline void swap(TSelf& hash_map) {
                std::swap(_ht_size, hash_map._ht_size);
                std::swap(_max_size, hash_map._max_size);
                _ht.swap(hash_map._ht);
                _poses.swap(hash_map._poses);
                std::swap(_last_key, hash_map._last_key);
                std::swap(_last_position, hash_map._last_position);
                std::swap(_empty_element, hash_map._empty_element);
        }

#ifdef PERFORMANCE_STATISTICS
        static void reset_statistics() {
                num_access = 0;
                num_contain = 0;
                overall_max_size = 0;
                num_find_pos = 0;
                num_probes = 0;
        }

        static void print_statistics() {
                std::cout << "Num access\t" << num_access << std::endl;
                std::cout << "Num contain\t" << num_contain << std::endl;
                std::cout << "Max overall max size\t" << overall_max_size << std::endl;
                std::cout << "Num find pos\t" << num_find_pos << std::endl;
                std::cout << "Num probes\t" << num_probes << std::endl;
                std::cout << "Average num prob per find pos\t" << (num_probes + 0.0) / num_find_pos << std::endl;
        }

        void print_full_statistics() const {
                std::cout << "Num access\t" << num_access << std::endl;
                std::cout << "Num contain\t" << num_contain << std::endl;
                std::cout << "Max overall max size\t" << overall_max_size << std::endl;
                std::cout << "Size\t" << size() << std::endl;
                std::cout << "Mem (table onle)\t" << _ht.size() * sizeof(Element) << std::endl;
                std::cout << "Num find pos\t" << num_find_pos << std::endl;
                std::cout << "Num probes\t" << num_probes << std::endl;
                std::cout << "Average num prob per find pos\t" << (num_probes + 0.0) / num_find_pos << std::endl;
        }
#endif
private:
        void resize() {
                TSelf new_hash_map(2 * _max_size);

                for (auto pos : _poses)
                        new_hash_map.insert(_ht[pos].first, _ht[pos].second);

                swap(new_hash_map);
        }

        inline void insertImpl(const Key& key, const Value& value) {
                if (TGrowable && size() == _max_size)
                        resize();

                const Position pos = findPosition(key);
                if (_ht[pos].first == _empty_element) {
                        _ht[pos].first = key;
                        _ht[pos].second = value;
                        _poses.push_back(pos);
                }
        }

        inline Position findPosition(const Key& key) {
#ifdef PERFORMANCE_STATISTICS
                ++num_find_pos;
#endif
                if (Cache && key == _last_key) {
                        return _last_position;
                }

                const Position startPosition = _hash(key) & (_ht_size - 1);
                for (Position pos = startPosition; pos < _ht.size(); ++pos) {
#ifdef PERFORMANCE_STATISTICS
                        ++num_probes;
#endif
                        if (_ht[pos].first == _empty_element || _ht[pos].first == key) {
                                if (Cache) {
                                        _last_key = key;
                                        _last_position = pos;
                                }
                                return pos;
                        }
                }

                for (Position pos = 0; pos < startPosition; ++pos) {
                        if (_ht[pos].first == _empty_element || _ht[pos].first == key) {
                                if (Cache) {
                                        _last_key = key;
                                        _last_position = pos;
                                }
                                return pos;
                        }
                }

                throw std::runtime_error("Hash table overflowed");
        }

        Key _empty_element;
        uint64_t _max_size;
        uint64_t _ht_size;
        std::vector<Element> _ht;
        std::vector<Position> _poses;
        Hash _hash;
        Key _last_key;
        Position _last_position;
};

#ifdef PERFORMANCE_STATISTICS
template <typename Key, typename Value, typename Hash, bool TGrowable, bool Cache, size_t SizeFactor>
size_t HashMap<Key, Value, Hash, TGrowable, Cache, SizeFactor>::num_access(0);

template <typename Key, typename Value, typename Hash, bool TGrowable, bool Cache, size_t SizeFactor>
size_t HashMap<Key, Value, Hash, TGrowable, Cache, SizeFactor>::num_contain(0);

template <typename Key, typename Value, typename Hash, bool TGrowable, bool Cache, size_t SizeFactor>
size_t HashMap<Key, Value, Hash, TGrowable, Cache, SizeFactor>::overall_max_size(0);

template <typename Key, typename Value, typename Hash, bool TGrowable, bool Cache, size_t SizeFactor>
size_t HashMap<Key, Value, Hash, TGrowable, Cache, SizeFactor>::num_probes(0);

template <typename Key, typename Value, typename Hash, bool TGrowable, bool Cache, size_t SizeFactor>
size_t HashMap<Key, Value, Hash, TGrowable, Cache, SizeFactor>::num_find_pos(0);
#endif

template <typename Key, typename Value, typename Hash = xxhash<Key>,
        bool TGrowable = false, bool Cache = true, size_t SizeFactor = 2>
class HashMapWithErase {
public:
        using Element = std::pair<Key, Value>;
        using key_type = Key;
        using mapped_type = Value;

private:
        using TSelf = HashMapWithErase<Key, Value, Hash, TGrowable, Cache, SizeFactor>;
        using Position = uint32_t;

public:

        explicit HashMapWithErase(const uint64_t max_size = 0) :
                _empty_element(std::numeric_limits<Key>::max()),
                _deleted_element(_empty_element - 1),
                _ht_size(max_size * SizeFactor),
                _max_size(max_size),
                _ht(_ht_size + _max_size * 1.1, std::make_pair(_empty_element, Value())),
                _poses(),
                _hash(),
                _last_key(_empty_element),
                _last_position(0) {
                ALWAYS_ASSERT(is_power_2(_ht_size));
                _poses.reserve(max_size);
        }

        static constexpr size_t size_factor = SizeFactor;

        HashMapWithErase(const TSelf&) = default;
        HashMapWithErase(TSelf&&) = default;

        TSelf& operator=(TSelf& other) = default;
        TSelf& operator=(TSelf&& other) = default;

        void reserve(const uint32_t max_size) {
                _ht_size = max_size * SizeFactor;
                _max_size = max_size;
                _ht.resize(_ht_size + _max_size * 1.1, std::make_pair(_empty_element.first, Value()));
                _poses.reserve(_max_size);
        }

        inline size_t size() const {
                return _poses.size();
        }

        inline bool empty() const {
                return size() == 0;
        }

        inline void erase(const Key& key) {
                const Position pos = findPosition(key);
                if (_ht[pos].first != _empty_element && _ht[pos].first != _deleted_element) {
                        _ht[pos].first = _deleted_element;
                }
        }

        inline Value& operator[] (const Key& key) {
                if (TGrowable && size() == _max_size / 2)
                        resize();

                const Position pos = findPosition(key);
                if (_ht[pos].first == _empty_element || _ht[pos].first == _deleted_element) {
                        _ht[pos].first = key;
                        _ht[pos].second = Value();
                        _poses.push_back(pos);
                }

                return _ht[pos].second;
        }

        inline bool contains(const Key& key) {
                // Sometimes findPosition can return _delete_element
                // which means that in searched until the end of _ht and did not
                // find element but found _deleted_element. This means that
                // element is NOT in hash table
                const Position pos = findPosition(key);
                if (_ht[pos].first == _empty_element || _ht[pos].first == _deleted_element) {
                        return false;
                }

                return true;
        }

        inline void insert(const Element& elem) {
                insertImpl(elem.first, elem.second);
        }

        inline void insert(const Key& key, const Value& value) {
                insertImpl(key, value);
        }

        inline void clear() {
                for (auto pos : _poses) {
                        _ht[pos].first = _empty_element;
                }
                _poses.clear();

                _last_key = _empty_element;
                _last_position = 0;
        }

        inline void swap(TSelf& hash_map) {
                std::swap(_ht_size, hash_map._ht_size);
                std::swap(_max_size, hash_map._max_size);
                _ht.swap(hash_map._ht);
                _poses.swap(hash_map._poses);
                std::swap(_last_key, hash_map._last_key);
                std::swap(_last_position, hash_map._last_position);
                std::swap(_empty_element, hash_map._empty_element);
                std::swap(_deleted_element, hash_map._deleted_element);
        }

private:
        void resize() {
                TSelf new_hash_map(2 * _max_size);

                for (auto pos : _poses)
                        new_hash_map.insert(_ht[pos]);

                swap(new_hash_map);
        }

        inline void insertImpl(const Key& key, const Value& value) {
                if (TGrowable && size() == _max_size / 2)
                        resize();

                const Position pos = findPosition(key);
                if (_ht[pos].first == _empty_element || _ht[pos].first == _deleted_element) {
                        _ht[pos].first = key;
                        _ht[pos].second = value;
                        _poses.push_back(pos);
                }
        }

        inline Position findPosition(const Key& key) {
                if (Cache && key == _last_key) {
                        return _last_position;
                }

                const Position startPosition = _hash(key) & (_ht_size - 1);
                Position firstDeleted = std::numeric_limits<Position>::max();
                for (Position pos = startPosition; pos < _ht.size(); ++pos) {
                        if (_ht[pos].first == _empty_element || _ht[pos].first == key) {
                                if (firstDeleted != std::numeric_limits<Position>::max() && _ht[pos].first == key) {
                                        std::swap(_ht[firstDeleted], _ht[pos]);

                                        pos = firstDeleted;
                                }

                                if (Cache) {
                                        _last_key = key;
                                        _last_position = pos;
                                }

                                return pos;
                        }
                        if (firstDeleted == std::numeric_limits<Position>::max() &&
                            _ht[pos].first == _deleted_element) {
                                firstDeleted = pos;
                        }
                }

                if (firstDeleted != std::numeric_limits<Position>::max()) {
                        if (Cache) {
                                _last_key = key;
                                _last_position = firstDeleted;
                        }
                        return firstDeleted;
                }

                std::cerr << "hash table overflowed" << std::endl;
                std::exit(-1);
        }

        key_type _empty_element;
        key_type _deleted_element;
        uint64_t _ht_size;
        uint64_t _max_size;
        std::vector<Element> _ht;
        std::vector<Position> _poses;
        Hash _hash;
        Key _last_key;
        Position _last_position;
};

template <typename Key, typename Hash = xxhash<Key>, bool TGrowable = false,
        bool Cache = true, size_t SizeFactor = 2>
class HashSet {
public:
        using Element = Key;
        using hash_type = typename Hash::hash_type;
private:
        using TSelf = HashSet<Key, Hash, TGrowable, Cache, SizeFactor>;
        using Position = uint32_t;

        static constexpr hash_type max_hash_value = std::numeric_limits<hash_type>::max();
public:
        static constexpr size_t size_factor = SizeFactor;

        explicit HashSet(const uint64_t max_size = 0) :
                _empty_element(std::numeric_limits<Key>::max()),
                _ht_size(max_size * SizeFactor),
                _max_size(max_size),
                _ht(_ht_size + _max_size * 1.1, _empty_element),
                _poses(),
                _hash(),
                _last_key(_empty_element),
                _last_position(0) {
                ALWAYS_ASSERT(is_power_2(_ht_size));
                _poses.reserve(max_size);
        }

        HashSet(const TSelf&) = default;
        HashSet(TSelf&&) = default;

        TSelf& operator= (TSelf& other) = default;
        TSelf& operator= (TSelf&& other) = default;

        void reserve(const uint32_t max_size) {
                _ht_size = max_size * SizeFactor;
                _max_size = max_size;
                _ht.resize(_ht_size + _max_size * 1.1, _empty_element);
                _poses.reserve(_max_size);
        }

        inline uint64_t size() const {
                return _poses.size();
        }

        inline bool empty() const {
                return size() == 0;
        }

        inline bool contains(const Key& key) {
                if (_ht[findPosition(key)] == _empty_element) {
                        return false;
                }
                return true;
        }

        inline void insert(const Key& key) {
                insertImpl(key);
        }

        inline void clear() {
                for (auto pos : _poses) {
                        _ht[pos] = _empty_element;
                }

                _poses.clear();

                _last_key = _empty_element;
                _last_position = 0;
        }

        inline void swap(TSelf& hash_set) {
                std::swap(_ht_size, hash_set._ht_size);
                std::swap(_max_size, hash_set._max_size);
                _ht.swap(hash_set._ht);
                _poses.swap(hash_set._poses);
                std::swap(_last_key, hash_set._last_key);
                std::swap(_last_position, hash_set._last_position);
        }

private:
        void resize() {
                TSelf new_hash_map(2 * _max_size);

                for (auto pos : _poses)
                        new_hash_map.insert(_ht[pos]);

                swap(new_hash_map);
        }

        inline void insertImpl(const Key& key) {
                if (TGrowable && size() == _max_size / 2)
                        resize();

                const Position pos = findPosition(key);
                if (_ht[pos] == _empty_element) {
                        _ht[pos] = key;
                        _poses.push_back(pos);
                }
        }

        inline Position findPosition(const Key& key) {
                if (Cache && key == _last_key) {
                        return _last_position;
                }

                //const Position startPosition = _hash(key) & (_ht_size - 1);
                const Position startPosition = _hash(key) % _ht_size;
                for (Position pos = startPosition; pos < _ht.size(); ++pos) {
                        if (_ht[pos] == _empty_element || _ht[pos] == key) {
                                if (Cache) {
                                        _last_key = key;
                                        _last_position = pos;
                                }
                                return pos;
                        }
                }

                throw std::runtime_error("Hash table overflowed");
        }

        const Key _empty_element;
        uint64_t _ht_size;
        uint64_t _max_size;
        std::vector<Element> _ht;
        std::vector<Position> _poses;
        Hash _hash;
        Key _last_key;
        Position _last_position;
};

//template <typename key_type, typename value_type>
//using hash_map = HashMap<key_type, value_type, simple_hash<key_type>, true>;

//template <typename key_type, typename value_type>
//using hash_map_with_erase = HashMapWithErase<key_type, value_type, simple_hash<key_type>, true>;
//
//template <typename key_type>
//using hash_set = HashSet<key_type, simple_hash<key_type>, true>;

//template <typename key_type, typename value_type>
//using hash_map = HashMap<key_type, value_type, xxhash<key_type>, true>;
//
//template <typename key_type, typename value_type>
//using hash_map_with_erase = HashMapWithErase<key_type, value_type, xxhash<key_type>, true>;
//
//template <typename key_type>
//using hash_set = HashSet<key_type, xxhash<key_type>, true>;

template <typename key_type, typename value_type>
using hash_map = HashMap<key_type, value_type, MurmurHash<key_type>, true, false>;

template <typename key_type, typename value_type>
using hash_map_with_erase = HashMapWithErase<key_type, value_type, MurmurHash<key_type>, true>;

template <typename key_type>
using hash_set = HashSet<key_type, MurmurHash<key_type>, true>;

}
