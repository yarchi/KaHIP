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
        using hash_type = uint64_t;
        inline hash_type operator()(const T& x) const {
                return XXH64(&x, sizeof(x), 0);
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

constexpr static bool is_power_2(uint64_t x) {
        return ((x != 0) && !(x & (x - 1)));
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


template <typename Key, typename Value, typename Hash = xxhash<Key>,
        bool TGrowable = false, bool Cache = true, size_t SizeFactor = 2>
class HashMap {
public:
        using Element = std::pair<Key, Value>;
        using key_type = Key;
        using mapped_type = Value;
        using hash_type = typename Hash::hash_type;

private:
        using TSelf = HashMap<Key, Value, Hash, TGrowable, Cache, SizeFactor>;
        using Position = uint32_t;

        static constexpr hash_type max_hash_value = std::numeric_limits<hash_type>::max();
public:
        using Iterator = HashTableIterator<TSelf>;

        friend Iterator;

        explicit HashMap(const uint64_t max_size = 0) :
                _empty_element(std::numeric_limits<Key>::max()),
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

        HashMap(const TSelf&) = default;
        HashMap(TSelf&&) = default;

        TSelf& operator=(TSelf& other) = default;
        TSelf& operator=(TSelf&& other) = default;

        static constexpr size_t get_max_size_to_fit_l1() {
                // (SizeFactor + 1.1) * max_size * sizeof(Element) + sizeof(Position) * max_size Bytes = 16 * 1024 Bytes,
                // where 16 * 1024 Bytes is half of L1 cache and (2 + 1.1) * max_size * 8 + 4 * max_size Bytes
                // is the size of a hash table. We calculate that max_size ~ 560.
                return round_up_to_next_power_2(16 * 1024 / (sizeof(Element) * (SizeFactor + 1.1) + sizeof(Position)) - sizeof(TSelf)) / SizeFactor;
        }

        Iterator begin() const {
                return Iterator(*this, 0);
        }

        Iterator end() const {
                return Iterator(*this, size());
        }

        void reserve(const uint32_t max_size) {
                _ht_size = max_size * SizeFactor;
                _max_size = max_size;
                _ht.resize(_ht_size + _max_size * 1.1, std::make_pair(_empty_element, Value()));
                _poses.reserve(_max_size);
        }

        inline uint64_t size() const {
                return _poses.size();
        }

        inline bool empty() const {
                return size() == 0;
        }

        inline Value& operator[](const Key& key) {
                if (TGrowable && size() == _max_size / 2)
                        resize();

                const Position pos = findPosition(key);
                if (_ht[pos].first == _empty_element) {
                        _ht[pos].first = key;
                        _ht[pos].second = Value();
                        _poses.push_back(pos);
                }

                return _ht[pos].second;
        }

        inline bool contains(const Key& key) {
                return _ht[findPosition(key)].first != _empty_element;
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
        }

private:
        void resize() {
                TSelf new_hash_map(2 * _max_size);

                for (auto pos : _poses)
                        new_hash_map.insert(_ht[pos].first, _ht[pos].second);

                swap(new_hash_map);
        }

        inline void insertImpl(const Key& key, const Value& value) {
                if (TGrowable && size() == _max_size / 2)
                        resize();

                const Position pos = findPosition(key);
                if (_ht[pos].first == _empty_element) {
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
                for (Position pos = startPosition; pos < _ht.size(); ++pos) {
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
        uint64_t _ht_size;
        uint64_t _max_size;
        std::vector<Element> _ht;
        std::vector<Position> _poses;
        Hash _hash;
        Key _last_key;
        Position _last_position;
};

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
                _pos_in_position(_ht_size + _max_size * 1.1),
                _hash(),
                _last_key(_empty_element),
                _last_position(0) {
                ALWAYS_ASSERT(is_power_2(_ht_size));
                _poses.reserve(max_size);
        }

        HashMapWithErase(const TSelf&) = default;
        HashMapWithErase(TSelf&&) = default;

        TSelf& operator=(TSelf& other) = default;
        TSelf& operator=(TSelf&& other) = default;

        static constexpr size_t get_max_size_to_fit_l1() {
                // (SizeFactor + 1.1) * max_size * sizeof(Element) + (sizeof(Position) + sizeof(uint32_t)) * max_size Bytes = 16 * 1024 Bytes,
                // where 16 * 1024 Bytes is half of L1 cache and (2 + 1.1) * max_size * 8 + 4 * max_size Bytes
                // is the size of a hash table. We calculate that max_size ~ 560.
                return round_up_to_next_power_2(16 * 1024 / (sizeof(Element) * (SizeFactor + 1.1) + sizeof(Position) +
                        sizeof(uint32_t)) - sizeof(TSelf)) / SizeFactor;
        }

        void reserve(const uint32_t max_size) {
                _ht_size = max_size * SizeFactor;
                _max_size = max_size;
                _ht.resize(_ht_size + _max_size * 1.1, std::make_pair(_empty_element.first, Value()));
                _pos_in_position.resize(_ht_size + _max_size * 1.1);
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
                        const uint32_t pos_of_deleted = _pos_in_position[pos];
                        _pos_in_position[_poses.back()] = pos_of_deleted;
                        std::swap(_poses[pos_of_deleted], _poses.back());
                        _poses.resize(_poses.size() - 1);
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
                        _pos_in_position[pos] = _poses.size() - 1;
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
                _pos_in_position.swap(hash_map._pos_in_position);
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
                        _pos_in_position[pos] = _poses.size() - 1;
                }
        }

        inline Position findPosition(const Key& key) {
                if (Cache && key == _last_key) {
                        return _last_position;
                }

                const Position startPosition = _hash(key) & (_ht_size - 1);
                Position lastDeleted = std::numeric_limits<Position>::max();
                for (Position pos = startPosition; pos < _ht.size(); ++pos) {
                        if (_ht[pos].first == _empty_element || _ht[pos].first == key) {
                                if (Cache) {
                                        _last_key = key;
                                        _last_position = pos;
                                }
                                return pos;
                        }
                        if (_ht[pos].first == _deleted_element) {
                                lastDeleted = pos;
                        }
                }

                if (lastDeleted != std::numeric_limits<Position>::max()) {
                        if (Cache) {
                                _last_key = key;
                                _last_position = lastDeleted;
                        }
                        return lastDeleted;
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
        std::vector<uint32_t> _pos_in_position;
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

        static constexpr size_t get_max_size_to_fit_l1() {
                // (SizeFactor + 1.1) * max_size * sizeof(Element) + sizeof(Position) * max_size Bytes = 16 * 1024 Bytes,
                // where 16 * 1024 Bytes is half of L1 cache and (2 + 1.1) * max_size * 8 + 4 * max_size Bytes
                // is the size of a hash table. We calculate that max_size ~ 560.
                return round_up_to_next_power_2(16 * 1024 / (sizeof(Element) * (SizeFactor + 1.1) + sizeof(Position)) - sizeof(TSelf)) / SizeFactor;
        }

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

                const Position startPosition = _hash(key) & (_ht_size - 1);
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

template <typename key_type, typename value_type>
using hash_map = HashMap<key_type, value_type, simple_hash<key_type>, true>;

template <typename key_type, typename value_type>
using hash_map_with_erase = HashMapWithErase<key_type, value_type, simple_hash<key_type>, true>;

template <typename key_type>
using hash_set = HashSet<key_type, simple_hash<key_type>, true>;

}
