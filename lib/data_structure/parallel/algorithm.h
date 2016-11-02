#pragma once

#include <algorithm>
#include <cmath>
#include <iterator>
#include <future>
#include <vector>

namespace parallel {

template <typename Functor>
static Functor range(size_t first, size_t last, Functor func) {
        for (auto elem = first; elem != last; ++elem) {
                func(elem);
        }
        return func;
}

template <typename Functor, typename TaskConsumer>
static Functor range(size_t first, size_t last, TaskConsumer& task_consumer, Functor func, bool sync) {
        size_t size = last - first;
        size_t proc = task_consumer.NumThreads() + 1;
        if (size <= 5 * proc) {
                return range(first, last, func);
        }

        size_t work_per_thread = size / proc;

        std::vector<std::future<Functor>> futures;
        if (sync) {
                futures.reserve(proc - 1);
        }

        for (size_t i = 0; i + 1 < proc; ++i) {
                if (sync) {
                        futures.push_back(task_consumer.Submit(range<Functor>, first,
                                                               first + work_per_thread, func));
                } else {
                        task_consumer.Submit(range<Functor>, first, first + work_per_thread, func);
                }
                first += work_per_thread;
        }

        range(first, last, func);

        if (sync) {
                std::for_each(futures.begin(), futures.end(), [](auto& future) {
                        future.get();
                });
        }
        return func;
}

template <typename Functor, typename TaskConsumer>
Functor range_sync(size_t first, size_t last, TaskConsumer& task_consumer, Functor func)
{
        return range(first, last, task_consumer, func, true);
}

template <typename Functor, typename TaskConsumer>
Functor range_async(size_t first, size_t last, TaskConsumer& task_consumer, Functor func)
{
        return range(first, last, task_consumer, func, false);
}

template <typename Iterator, typename Functor, typename TaskConsumer>
static Functor for_each(Iterator first, Iterator last, TaskConsumer& task_consumer, Functor func,
                 std::random_access_iterator_tag, bool sync) {
        size_t size = std::distance(first, last);
        size_t proc = task_consumer.NumThreads() + 1;
        if (size <= 5 * proc) {
                return std::for_each(first, last, func);
        }

        size_t work_per_thread = size / proc;

        std::vector<std::future<Functor>> futures;
        if (sync) {
                futures.reserve(proc - 1);
        }

        for (size_t i = 0; i + 1 < proc; ++i) {
                if (sync) {
                        futures.push_back(task_consumer.Submit(std::for_each<Iterator, Functor>, first,
                                                               first + work_per_thread, func));
                } else {
                        task_consumer.Submit(std::for_each<Iterator, Functor>, first, first + work_per_thread, func);
                }
                first += work_per_thread;
        }

        std::for_each(first, last, func);

        if (sync) {
                std::for_each(futures.begin(), futures.end(), [](auto& future) {
                        future.get();
                });
        }
        return func;
}

template <typename Iterator, typename Functor, typename TaskConsumer>
Functor for_each_sync(Iterator first, Iterator last, TaskConsumer& task_consumer, Functor func)
{
        return for_each(first, last, task_consumer, func,
                        typename std::iterator_traits<Iterator>::iterator_category(), true);
}

template <typename Iterator, typename Functor, typename TaskConsumer>
Functor for_each_async(Iterator first, Iterator last, TaskConsumer& task_consumer, Functor func)
{
        return for_each(first, last, task_consumer, func,
                        typename std::iterator_traits<Iterator>::iterator_category(), false);
}

const uint32_t g_cache_line_size = 64 * sizeof(char);

template <typename T>
class alignas(g_cache_line_size) CacheAlignedData {
public:
        using Type = T;

        inline operator T() const {
                return m_elem;
        }

        CacheAlignedData(T elem = T())
                :       m_elem(elem)
        {
                static_assert(sizeof(Self) % g_cache_line_size == 0, "No cache line alignment");
        }

        CacheAlignedData(const CacheAlignedData& other)
                :       m_elem(other.m_elem)
        {}

        CacheAlignedData(CacheAlignedData&& other)
                :       m_elem(other.m_elem)
        {}

        T& get() {
                return m_elem;
        }

        const T& get() const {
                return m_elem;
        }

        bool operator! () const {
                return !m_elem;
        }

        CacheAlignedData& operator= (const T& elem) {
                m_elem = elem;
                return *this;
        }

        CacheAlignedData& operator= (const CacheAlignedData& other) {
                m_elem = other.m_elem;
                return *this;
        }

        CacheAlignedData& operator= (CacheAlignedData&& other) {
                m_elem = other.m_elem;
                return *this;
        }

        CacheAlignedData& operator+= (const T& elem) {
                m_elem += elem;
                return *this;
        }

        CacheAlignedData& operator-= (const T& elem) {
                m_elem -= elem;
                return *this;
        }

        CacheAlignedData& operator*= (const T& elem) {
                m_elem *= elem;
                return *this;
        }

        CacheAlignedData& operator/= (const T& elem) {
                m_elem /= elem;
                return *this;
        }

        T& operator++ () {
                ++m_elem;
                return *this;
        }

        T operator++ (int) {
                T tmp = m_elem;
                ++m_elem;
                return tmp;
        }
private:
        using Self = CacheAlignedData<T>;
        T m_elem;

};

template <typename T>
using Cvector = std::vector<CacheAlignedData<T>>;

template <typename T, typename TaskConsumer>
class ParallelVector {
public:
#ifdef CACHE_LINE_ALIGNMENT
        static const bool m_cache_line_alignment = true;
#else
        static const bool m_cache_line_alignment = false;
#endif
        using Self = ParallelVector<T, TaskConsumer>;
        using Type = typename std::conditional<m_cache_line_alignment, CacheAlignedData<T>, T>::type;

        const Type& operator[] (size_t index) const {
                return m_ptr[index];
        }

        Type& operator[] (size_t index) {
                return m_ptr[index];
        }

        size_t size() const {
                return m_size;
        }

        void swap(Self& other) {
                std::swap(m_ptr, other.m_ptr);
                std::swap(m_size, other.m_size);
        }

        ~ParallelVector() {
                parallel_destruction(m_task_consumer, m_ptr, m_size);
        }

        const Type* begin() const {
                return m_ptr;
        }

        Type* begin() {
                return m_ptr;
        }

        const Type* end() const {
                return m_ptr + m_size;
        }

        Type* end() {
                return m_ptr + m_size;
        }
private:
        explicit ParallelVector(TaskConsumer& task_consumer, size_t size, T init_value = 0)
                :       m_ptr(nullptr)
                ,       m_task_consumer(task_consumer)
                ,       m_size(size)
        {
                m_ptr = parallel_construction<Type>(m_task_consumer, size, init_value);
        }

        template <typename V>
        V* parallel_construction(TaskConsumer& task_consumer, size_t size, T init_value) {
                V* ptr = reinterpret_cast<V*>(::operator new(size * sizeof(V)));
                for_each_sync(ptr, ptr + size, task_consumer, [&init_value](V& p) {
                        new (&p) V(init_value);
                });
                return ptr;
        }

        template <typename V>
        void parallel_destruction(TaskConsumer& task_consumer, V* ptr, size_t size) {
                for_each_sync(ptr, ptr + size, task_consumer, [](V& p) {
                        p.~V();
                });
        }

        template <typename G, typename V>
        friend ParallelVector<G, V> Get_parallel_vector(V&, size_t, G);

        Type* m_ptr;
        TaskConsumer& m_task_consumer;
        size_t m_size;
};

template <typename T, typename TaskConsumer>
static ParallelVector<T, TaskConsumer> Get_parallel_vector(TaskConsumer& task_consumer, size_t size, T init_value = 0) {
        return ParallelVector<T, TaskConsumer>(task_consumer, size, init_value);
};

}