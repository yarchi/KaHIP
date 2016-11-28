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

        void set_seed(uint32_t seed) {
                m_mt.seed(seed);
        }
private:
        std::uniform_int_distribution<uint32_t> m_rnd_bit;
        std::mt19937 m_mt;
};
}