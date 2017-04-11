#pragma once

#include <random>

namespace parallel {
class random {
public:
        explicit random(uint32_t seed)
                :       m_rnd_bit(0, 1)
                ,       m_mt(seed)
        {}

        inline bool bit() {
                return m_rnd_bit(m_mt);
        }

        // returns random number in [a, b]
        template <typename T = int>
        inline T random_number(T a = std::numeric_limits<T>::min(), T b = std::numeric_limits<T>::max()) {
                std::uniform_int_distribution<T> rnd(a, b);
                return rnd(m_mt);
        }

        void set_seed(uint32_t seed) {
                m_mt.seed(seed);
        }

        template <typename T>
        inline void shuffle(std::vector<T>& vec) {
                shuffle(vec.begin(), vec.end());
        }

        template <typename iterator_type>
        inline void shuffle(iterator_type begin, iterator_type end) {
                size_t size = end - begin;
                if (size < 2) {
                        return;
                }

                for (size_t i = 0; i < size; ++i) {
                        std::uniform_int_distribution<uint32_t> rnd(i, size - 1);
                        size_t rnd_ind = rnd(m_mt);
                        std::swap(*(begin + i), *(begin + rnd_ind));
                }
        }

private:
        std::uniform_int_distribution<uint32_t> m_rnd_bit;
        std::mt19937 m_mt;
};
}