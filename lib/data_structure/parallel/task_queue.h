#pragma once

#include "data_structure/parallel/algorithm.h"
#include "data_structure/parallel/atomics.h"

namespace parallel {

template <typename T>
class thread_container {
public:
        explicit thread_container(size_t size = 100) {
                m_elems.reserve(size);
                m_counter.store(0, std::memory_order_relaxed);
        }

        inline bool try_pop(T& elem) {
                size_t offset = m_counter.fetch_add(1, std::memory_order_relaxed);

                if (offset < m_elems.size()) {
                        elem = m_elems[offset];
                        return true;
                } else {
                        return false;
                }
        }

        inline void push_back(const T& elem) {
                m_elems.push_back(elem);
        }

        inline void clear() {
                m_elems.clear();
        }

private:
        std::vector<T> m_elems;
        AtomicWrapper<T> m_counter;
};

template <typename T>
class task_queue {
public:
        explicit task_queue(size_t size) {
                m_thread_containers.resize(size);
                m_counter.store(0, std::memory_order_relaxed);
        }

        inline bool try_pop(T& elem) {
                size_t offset = m_counter.load(std::memory_order_relaxed);

                bool res = false;
                do {
                        offset = m_counter.load(std::memory_order_relaxed);

                        if ((res = m_thread_containers[offset].get().try_pop(elem)) || empty()) {
                                break;
                        }
                } while (!m_counter.compare_exchange_weak(offset, offset + 1, std::memory_order_relaxed));

                return res;
        }

        inline bool empty() const {
                return m_counter.load(std::memory_order_relaxed) < m_thread_containers.size();
        }

        inline thread_container<T>& operator[](size_t thread_id) {
                return m_thread_containers[thread_id].get();
        }

        inline const thread_container<T>& operator[](size_t thread_id) const {
                return m_thread_containers[thread_id].get();
        }

        inline void push(const T& elem) {
                m_thread_containers[0].get().push_back(elem);
        }

        inline void clear() {
                for (auto& container : m_thread_containers) {
                        container.get().clear();
                }
        }
private:
        Cvector<thread_container<T>> m_thread_containers;
        std::atomic<size_t> m_counter;
};

};