#include <cstdint>
#include <vector>

namespace parallel {
/**
 * [key_2, key_4, key_5, key_3, key_1           , val_1, pos_of_key_1, val_2, pos_of_key_2, val_3, pos_of_key_3, val_4, pos_of_key_4, val_5, pos_of_key_5, -1, -1]
 * [******************m_max_size**************, ******************m_max_size**************, ******************m_max_size**************]
 **/
template <typename type>
class fast_set {
public:
        using iterator_type = typename std::vector<type>::iterator;

        fast_set(PartitionID max_size)
                :       m_max_size(max_size)
                ,       m_size(0)
                ,       m_data(3 * max_size)
        {
                for (size_t i = m_max_size; i < m_data.size(); ++i) {
                        m_data[i] = m_empty;
                }
        }

        void insert(type key, type val) {
                type pos = m_size++;
                m_data[pos] = key;
                m_data[m_max_size + 2 * key] = val;
                m_data[m_max_size + 2 * key + 1] = pos;
        }

        void remove(type key) {
                if (m_size == 0) {
                        return;
                }
                type pos = m_data[m_max_size + 2 * key + 1];

                m_data[m_max_size + 2 * key] = m_empty;
                m_data[m_max_size + 2 * key + 1] = m_empty;

                std::swap(m_data[m_size - 1], m_data[pos]);
                --m_size;

                key = m_data[pos];
                m_data[m_max_size + 2 * key + 1] = pos;
        }

        iterator_type begin() {
                return m_data.begin();
        }

        iterator_type end() {
                return m_data.begin() + m_size;
        }

        type& operator[](type key) {
                if (m_data[m_max_size + 2 * key] == m_empty) {
                        insert(key, 0);
                }
                return m_data[m_max_size + 2 * key];
        }

        const type& operator[](type key) const {
                return m_data[m_max_size + 2 * key];
        }
private:
        static constexpr type m_empty = std::numeric_limits<type>::max();

        uint32_t m_max_size;
        uint32_t m_size;
        std::vector<type> m_data;
};

}
