#pragma once

#include <algorithm>
#include <cmath>
#include <iterator>
#include <future>
#include <vector>

namespace parallel {
template<typename _value_type>
class integer_iterator : public std::iterator<std::random_access_iterator_tag, _value_type> {
private:
        using base = std::iterator<std::random_access_iterator_tag, _value_type>;
public:
        typedef typename base::iterator_category iterator_category;
        typedef typename base::value_type value_type;
        typedef typename base::difference_type difference_type;
        typedef typename base::pointer pointer;
        typedef typename base::reference reference;

        typedef integer_iterator<value_type> self;

        inline integer_iterator(value_type value)
                : m_value(value) {}

        inline value_type operator*() const {
                return m_value;
        }

        inline difference_type operator-(self& other) const {
                return m_value - other.m_value;
        }

        inline self operator+(difference_type val) const {
                return self(m_value + val);
        }

        inline self& operator+=(difference_type val) {
                m_value += val;
                return *this;
        }

        inline self& operator++() {
                ++m_value;
                return *this;
        }

        inline self operator++(int) {
                value_type old = m_value;
                ++m_value;
                return self(old);
        }

        inline bool operator==(const self& other) const {
                return m_value == other.m_value;
        }

        inline bool operator!=(const self& other) const {
                return m_value != other.m_value;
        }

private:
        static_assert(std::is_integral<value_type>::value, "Value should be of integral type");

        value_type m_value;
};

template<typename Iterator, typename Function>
Function apply_to_range(Iterator first, Iterator last, Function f)
{
        f(first, last);
        return f;
}

template<bool Sync, template <typename, typename> typename AggregateFuncWrapper,
        typename Iterator, typename Functor, typename TaskConsumer>
static Functor for_each(Iterator first, Iterator last, TaskConsumer& task_consumer, Functor func,
                        std::random_access_iterator_tag) {
        size_t size = std::distance(first, last);
        size_t proc = task_consumer.NumThreads() + 1;
        AggregateFuncWrapper<Iterator, Functor> wrapper;

        if (size <= 5 * proc) {
                return wrapper(first, last, func);
        }

        size_t work_per_thread = size / proc;

        std::vector<std::future<Functor>> futures;
        if (Sync) {
                futures.reserve(proc - 1);
        }

        for (size_t i = 0; i + 1 < proc; ++i) {
                if (Sync) {
                        futures.push_back(task_consumer.Submit(i, wrapper, first,
                                                               first + work_per_thread, func));
//                        futures.push_back(task_consumer.Submit(wrapper, first,
//                                                               first + work_per_thread, func));
                } else {
                        task_consumer.Submit(i, wrapper,
                                             first, first + work_per_thread, func);
//                        task_consumer.Submit(wrapper,
//                                             first, first + work_per_thread, func);
                }
                first += work_per_thread;
        }

        wrapper(first, last, func);

        if (Sync) {
                std::for_each(futures.begin(), futures.end(), [](auto& future) {
                        future.get();
                });
        }
        return func;
}

template <typename iterator_type, bool is_integral = std::is_integral<iterator_type>::value>
struct iterator
{
        using type = iterator_type;
        static type& get(type& iter) {
                return iter;
        }
};

template <typename value_type>
struct iterator<value_type, true>
{
        using type = integer_iterator<value_type>;

        static type get(value_type val) {
                return integer_iterator<value_type>(val);
        }
};

template <typename Iterator, typename Functor>
struct for_each_wrapper {
        inline Functor operator() (Iterator begin, Iterator end, Functor f) const {
                return std::for_each(begin, end, f);
        }
};

template <typename Iterator, typename Functor>
struct apply_to_range_wrapper {
        inline Functor operator() (Iterator begin, Iterator end, Functor f) const {
                return parallel::apply_to_range(begin, end, f);
        }
};

template<typename Iterator, typename Functor, typename TaskConsumer>
Functor for_each_sync(Iterator first, Iterator last, TaskConsumer& task_consumer, Functor func) {
        using iter_type = typename iterator<Iterator>::type;
        auto getter = iterator<Iterator>::get;

        return for_each<true, for_each_wrapper>(getter(first), getter(last), task_consumer, func,
                                                typename std::iterator_traits<iter_type>::iterator_category());
}

template<typename Iterator, typename Functor, typename TaskConsumer>
Functor for_each_async(Iterator first, Iterator last, TaskConsumer& task_consumer, Functor func) {
        using iter_type = typename iterator<Iterator>::type;
        auto getter = iterator<Iterator>::get;
        return for_each<false, for_each_wrapper>(getter(first), getter(last), task_consumer, func,
                                                 typename std::iterator_traits<iter_type>::iterator_category());
}

template<typename Iterator, typename Functor, typename TaskConsumer>
Functor apply_to_range_sync(Iterator first, Iterator last, TaskConsumer& task_consumer, Functor func) {
        using iter_type = typename iterator<Iterator>::type;
        auto getter = iterator<Iterator>::get;

        return for_each<true, apply_to_range_wrapper>(getter(first), getter(last), task_consumer, func,
                                                      typename std::iterator_traits<iter_type>::iterator_category());
}

template<typename Iterator, typename Functor, typename TaskConsumer>
Functor apply_to_range_async(Iterator first, Iterator last, TaskConsumer& task_consumer, Functor func) {
        using iter_type = typename iterator<Iterator>::type;
        auto getter = iterator<Iterator>::get;
        return for_each<false, apply_to_range_wrapper>(getter(first), getter(last), task_consumer, func,
                                                       typename std::iterator_traits<iter_type>::iterator_category());
}

constexpr uint32_t g_cache_line_size = 64 * sizeof(char);

template <typename T>
class alignas(g_cache_line_size) CacheAlignedData {
public:
        using Type = T;

        CacheAlignedData()
                :       m_elem()
        {
                static_assert(sizeof(Self) % g_cache_line_size == 0, "No cache line alignment");
        }

        CacheAlignedData(T elem)
                :       m_elem(elem)
        {
                static_assert(sizeof(Self) % g_cache_line_size == 0, "No cache line alignment");
        }

        template <typename... Args>
        CacheAlignedData(Args&&... args)
                :       m_elem(std::forward<Args>(args)...)
        {}


        CacheAlignedData(const CacheAlignedData& other) = default;
        CacheAlignedData(CacheAlignedData&& other) = default;

        inline T& get() {
                return m_elem;
        }

        inline const T& get() const {
                return m_elem;
        }

        inline bool operator! () const {
                return !m_elem;
        }

        inline CacheAlignedData& operator= (const T& elem) {
                m_elem = elem;
                return *this;
        }

        inline CacheAlignedData& operator= (const CacheAlignedData& other) {
                m_elem = other.m_elem;
                return *this;
        }

        inline CacheAlignedData& operator= (CacheAlignedData&& other) {
                m_elem = other.m_elem;
                return *this;
        }

        inline CacheAlignedData& operator+= (const T& elem) {
                m_elem += elem;
                return *this;
        }

        inline CacheAlignedData& operator-= (const T& elem) {
                m_elem -= elem;
                return *this;
        }

        inline CacheAlignedData& operator*= (const T& elem) {
                m_elem *= elem;
                return *this;
        }

        inline CacheAlignedData& operator/= (const T& elem) {
                m_elem /= elem;
                return *this;
        }

        inline T& operator++ () {
                ++m_elem;
                return *this;
        }

        inline T operator++ (int) {
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