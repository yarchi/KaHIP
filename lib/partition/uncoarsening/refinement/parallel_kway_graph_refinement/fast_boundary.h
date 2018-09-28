#pragma once
#include "data_structure/graph_access.h"
#include "data_structure/parallel/graph_algorithm.h"
#include "data_structure/parallel/hash_table.h"
#include "data_structure/parallel/task_queue.h"
#include "data_structure/parallel/time.h"
#include "definitions.h"

#include <functional>
#include <vector>
#include <utility>

#include <tbb/concurrent_queue.h>

namespace parallel {


class fast_boundary {
public:
        fast_boundary(graph_access& G, const PartitionConfig& config)
                :       m_G(G)
                ,       m_config(config)
                ,       m_blocks_info(m_G.get_partition_count())
        {}

        inline NodeWeight get_block_weight(PartitionID partition) {
                return m_blocks_info[partition].block_weight;
        }
        inline NodeID get_block_size(PartitionID partition) {
                return m_blocks_info[partition].block_size;
        }

        inline void set_block_weight(PartitionID partition, NodeWeight weight) {
                m_blocks_info[partition].block_weight = weight;
        }
        inline void set_block_size(PartitionID partition, NodeID size) {
                m_blocks_info[partition].block_size = size;
        }

        void balance_singletons() {
                for(size_t i = 0; i < m_singletons.size(); ++i) {
                        NodeWeight min = m_blocks_info[0].block_weight;
                        PartitionID p = 0;
                        for(size_t j = 0; j < m_blocks_info.size(); ++j) {
                                if (m_blocks_info[j].block_weight < min) {
                                        min = m_blocks_info[j].block_weight;
                                        p = j;
                                }
                        }

                        NodeID node = m_singletons[i];
                        if (m_blocks_info[p].block_weight + m_G.getNodeWeight(node) <= m_config.upper_bound_partition) {
                                m_blocks_info[m_G.getPartitionIndex(node)].block_weight -= m_G.getNodeWeight(node);
                                m_blocks_info[p].block_weight += m_G.getNodeWeight(node);
                                m_G.setPartitionIndex(node, p);
                        }
                }
        }

protected:
        using hash_map_with_erase_type = parallel::HashMapWithErase<NodeID, int32_t, parallel::xxhash<NodeID>, true, false>;

        struct block_data_type {
                NodeWeight block_weight = 0;
                NodeID block_size = 0;
        };

        graph_access& m_G;
        const PartitionConfig& m_config;
        std::vector<block_data_type> m_blocks_info;
        std::vector<NodeID> m_singletons;
};

class fast_sequential_boundary : public fast_boundary {
public:
        using iterator_type = hash_map_with_erase_type::Iterator;

        fast_sequential_boundary(graph_access& G, const PartitionConfig& config)
                :       fast_boundary(G, config)
                ,       m_boundary(1)
        {}

        iterator_type begin() const {
                return m_boundary.begin();
        }

        iterator_type end() const {
                return m_boundary.end();
        }

        void move(NodeID vertex, PartitionID from, PartitionID to) {
                ALWAYS_ASSERT(m_G.getPartitionIndex(vertex) == to);

                int32_t external_neighbors = 0;
                forall_out_edges(m_G, e, vertex) {
                        NodeID target = m_G.getEdgeTarget(e);
                        PartitionID target_block = m_G.getPartitionIndex(target);

                        if (target_block == to) {
                                int32_t new_count = --m_boundary[target];
                                if (new_count == 0) {
                                        m_boundary.erase(target);
                                }
                        }

                        if (target_block == from) {
                                ++m_boundary[target];
                        }

                        if (target_block != to) {
                                ++external_neighbors;
                        }
                } endfor
                if (external_neighbors == 0) {
                        m_boundary.erase(vertex);
                } else {
                        m_boundary[vertex] = external_neighbors;
                }
        }

        void construct_boundary() {
                std::vector<std::pair<NodeID, int32_t>> preliminary_boundary;
                preliminary_boundary.reserve(1000);
                m_singletons.reserve(100);

                forall_nodes(m_G, n) {
                        int32_t num_external_neighbors = 0;

                        if (m_G.getNodeDegree(n) == 0) {
                                m_singletons.push_back(n);
                        }

                        PartitionID cur_block = m_G.getPartitionIndex(n);
                        forall_out_edges(m_G, e, n) {
                                NodeID target = m_G.getEdgeTarget(e);
                                PartitionID target_block = m_G.getPartitionIndex(target);

                                if (cur_block != target_block) {
                                        ++num_external_neighbors;
                                }
                        } endfor
                        if (num_external_neighbors > 0) {
                                preliminary_boundary.emplace_back(n, num_external_neighbors);
                        }

                        ++m_blocks_info[cur_block].block_size;
                        m_blocks_info[cur_block].block_weight += m_G.getNodeWeight(n);
                } endfor

                m_boundary.reserve(2 * preliminary_boundary.size());
                for (const auto& elem : preliminary_boundary) {
                        m_boundary[elem.first] = elem.second;
                }
        }

        void check_boundary() {
                std::vector<NodeID> this_boundary;
                this_boundary.reserve(m_boundary.size());

                for (const auto& elem : m_boundary) {
                        this_boundary.push_back(elem.first);
                }

                std::sort(this_boundary.begin(), this_boundary.end());

                std::vector<NodeID> expected_boundary;
                expected_boundary.reserve(m_boundary.size());
                std::cout << "Boundary diff: " << std::endl;
                forall_nodes(m_G, n) {
                        PartitionID cur_block = m_G.getPartitionIndex(n);
                        bool boundary = false;
                        forall_out_edges(m_G, e, n) {
                                NodeID target = m_G.getEdgeTarget(e);
                                PartitionID target_block = m_G.getPartitionIndex(target);

                                if (cur_block != target_block) {
                                        boundary = true;
                                        break;
                                }
                        } endfor
                        if (boundary) {
                                expected_boundary.push_back(n);
                                if (!m_boundary.contains(n)) {
                                        std::cout << "(" << n << ", FN) ";
                                }
                        }
                        if (!boundary) {
                                if (m_boundary.contains(n)) {
                                        std::cout << "(" << n << ", FP) ";
                                }
                        }
                } endfor
                std::cout << std::endl;
                std::cout << "check size = " << m_boundary.size() << std::endl;
                if (expected_boundary != this_boundary) {
                        std::cout << "expected size = " << expected_boundary.size() << std::endl;
                        std::cout << "this size = " << this_boundary.size() << std::endl;
                        ALWAYS_ASSERT(expected_boundary == this_boundary);
                }
        }
private:
        hash_map_with_erase_type m_boundary;

};

class fast_parallel_boundary : public fast_boundary {
private:
public:
        fast_parallel_boundary(graph_access& G, const PartitionConfig& config)
                :       fast_boundary(G, config)
                ,       m_primary_hash((uint32_t) m_config.seed)
                ,       m_secondary_hash((uint32_t) m_config.seed + 1)
                ,       m_second_level_size(m_config.num_threads * m_config.num_threads)
        {}

        void move(NodeID vertex, PartitionID from, PartitionID to) {
                ALWAYS_ASSERT(m_G.getPartitionIndex(vertex) == to);

                int32_t external_neighbors = 0;
                forall_out_edges(m_G, e, vertex) {
                        NodeID target = m_G.getEdgeTarget(e);
                        PartitionID target_block = m_G.getPartitionIndex(target);

                        if (target_block == to) {
                                auto& target_boundary_ht = m_boundaries_per_thread[num_hash_table(target)].get();
                                int32_t new_count = --target_boundary_ht[target];
                                if (new_count == 0) {
                                        target_boundary_ht.erase(target);
                                }
                        }

                        if (target_block == from) {
                                auto& target_boundary_ht = m_boundaries_per_thread[num_hash_table(target)].get();
                                ++target_boundary_ht[target];
                        }

                        if (target_block != to) {
                                ++external_neighbors;
                        }
                } endfor

                auto& boundary_ht = m_boundaries_per_thread[num_hash_table(vertex)].get();
                if (external_neighbors == 0) {
                        boundary_ht.erase(vertex);
                } else {
                        boundary_ht[vertex] = external_neighbors;
                }
        }

        void construct_boundary() {
                using std::placeholders::_1;

                container_collection_type containers(m_config.num_threads, Cvector<thread_container_type>(m_config.num_threads));

                if (parallel::g_thread_pool.NumThreads() == 0) {
                        distribute_boundary_vertices(m_G, containers, [this](uint32_t vertex, uint32_t) {
                                return external_neighbors(vertex);
                        });
                } else {
                        if (m_config.use_numa_aware_graph && parallel::graph_is_large(m_G, m_config)) {
                                parallel::numa_aware_graph na_graph(m_G, m_config.num_threads, m_config.threads_per_socket);
                                na_graph.construct();
                                distribute_boundary_vertices(na_graph.get_handle(0), containers, [this, &na_graph](uint32_t vertex, uint32_t thread_id) {
                                        return is_boundary(na_graph.get_handle(thread_id), vertex);
                                });
                        } else {
                                distribute_boundary_vertices(m_G, containers, [this](uint32_t vertex, uint32_t) {
                                        return is_boundary(m_G, vertex);
                                });
                        }
                }

                // ------------------
                CLOCK_START;
                m_boundaries_per_thread.resize(m_config.num_threads);

                auto task_copy_to_hash_tables = [this, &containers] (uint32_t thread_id) {
                        NodeID total_size = 0;
                        auto& container = containers[thread_id];
                        for (const auto& sub_container : container) {
                                total_size += sub_container.get().size();
                        }

                        m_boundaries_per_thread[thread_id].get().reserve(std::max(total_size, 16u));

                        for (const auto& sub_container : container) {
                                for (const auto& elem : sub_container.get()) {
                                        m_boundaries_per_thread[thread_id].get().insert(elem.first, elem.second);
                                }
                        }
                        container.clear();
                        return total_size;
                };

                std::vector<std::future<NodeID>> futures_other;
                futures_other.reserve(parallel::g_thread_pool.NumThreads());
                for (size_t i = 0; i < parallel::g_thread_pool.NumThreads(); ++i) {
                        futures_other.push_back(parallel::g_thread_pool.Submit(i, task_copy_to_hash_tables, i + 1));
                }
                std::cout << task_copy_to_hash_tables(0) << " ";

                std::for_each(futures_other.begin(), futures_other.end(), [&](auto& future){
                        std::cout << future.get() << " ";
                });
                std::cout << std::endl;
                CLOCK_END("Copy to hash table");
        }

        const hash_map_with_erase_type& operator[] (uint32_t thread_id) const {
                return m_boundaries_per_thread[thread_id].get();
        }

        hash_map_with_erase_type& operator[] (uint32_t thread_id) {
                return m_boundaries_per_thread[thread_id].get();
        }

        void check_boundary() {
                std::vector<NodeID> this_boundary;

                size_t size = 0;
                for (uint32_t thread_id = 0; thread_id < m_config.num_threads; ++thread_id) {
                        size += m_boundaries_per_thread[thread_id].get().size();
                }
                this_boundary.reserve(size);

                for (uint32_t thread_id = 0; thread_id < m_config.num_threads; ++thread_id) {
                        for (const auto& elem : m_boundaries_per_thread[thread_id].get()) {
                                this_boundary.push_back(elem.first);
                        }
                }

                std::sort(this_boundary.begin(), this_boundary.end());

                std::vector<NodeID> expected_boundary;
                expected_boundary.reserve(size);
                std::cout << "Boundary diff: " << std::endl;
                forall_nodes(m_G, n) {
                        PartitionID cur_block = m_G.getPartitionIndex(n);
                        bool boundary = false;
                        forall_out_edges(m_G, e, n) {
                                NodeID target = m_G.getEdgeTarget(e);
                                PartitionID target_block = m_G.getPartitionIndex(target);

                                if (cur_block != target_block) {
                                        boundary = true;
                                        break;
                                }
                        } endfor

                        if (boundary) {
                                expected_boundary.push_back(n);
//                                if (!m_boundaries_per_thread[num_hash_table(n)].get().contains(n)) {
//                                     std::cout << "! " << n << " ";
//                                }
                        } else {
//                                if (m_boundaries_per_thread[num_hash_table(n)].get().contains(n)) {
//                                        std::cout << n << " ";
//                                }
                        }
                } endfor
                std::cout << std::endl;
                std::cout << "check size = " << size << std::endl;
                if (expected_boundary == this_boundary) {
                        std::cout << "expected size = " << expected_boundary.size() << std::endl;
                        std::cout << "this size = " << this_boundary.size() << std::endl;
                        ALWAYS_ASSERT(expected_boundary == this_boundary);
                }
        }

        ~fast_parallel_boundary() {
                auto task = [this](uint32_t thread_id) {
                        m_boundaries_per_thread[thread_id].get().clear();
                        // to resize hash tables in parallel
                        m_boundaries_per_thread[thread_id].get().reserve(1);
                };
                std::vector<std::future<void>> futures;
                futures.reserve(parallel::g_thread_pool.NumThreads());
                for (size_t i = 0; i < parallel::g_thread_pool.NumThreads(); ++i) {
                        futures.push_back(parallel::g_thread_pool.Submit(i, task, i + 1));
                }
                task(0);

                std::for_each(futures.begin(), futures.end(), [&](auto& future) {
                        future.get();
                });
        }

        void begin_movements() {
                if (m_containers.size() != m_config.num_threads) {
                        m_containers.resize(m_config.num_threads, Cvector<thread_container<NodeID>>(m_second_level_size));
                }
        }

        void add_vertex_to_check(NodeID vertex) {
                uint32_t thread_id = num_hash_table(vertex);
                uint32_t bucket_id = (uint32_t) m_secondary_hash(vertex) % m_second_level_size;
                m_containers[thread_id][bucket_id].get().concurrent_emplace_back(vertex);
        }

        void finish_movements() {
                auto task = [this](uint32_t thread_id) {
                        auto& container = m_containers[thread_id];
                        parallel::HashSet<NodeID, parallel::xxhash<NodeID>, true> processed_vertices(1024u);
                        for (auto& sub_container : container) {
                                for (const auto& vertex : sub_container.get()) {
                                        if (processed_vertices.contains(vertex)) {
                                                continue;
                                        }
                                        processed_vertices.insert(vertex);
                                        int32_t num_external_neighbors = is_boundary(m_G, vertex);
                                        if (num_external_neighbors > 0) {
                                                m_boundaries_per_thread[thread_id].get().insert(vertex, num_external_neighbors);
                                        } else {
                                                m_boundaries_per_thread[thread_id].get().erase(vertex);
                                        }
                                }
                                sub_container.get().clear();
                        }
                        return processed_vertices.size();
                };

                std::vector<std::future<size_t>> futures;
                futures.reserve(parallel::g_thread_pool.NumThreads());
                for (size_t i = 0; i < parallel::g_thread_pool.NumThreads(); ++i) {
                        futures.push_back(parallel::g_thread_pool.Submit(i, task, i + 1));
                }
                NodeID size = task(0);

                std::for_each(futures.begin(), futures.end(), [&](auto& future){
                        NodeID sz = future.get();
                        size += sz;
                });
                std::cout << "TOTAL SIZE = " << size << std::endl;
        }
private:
        using thread_container_type = thread_container<std::pair<NodeID, int32_t>>;
        using container_collection_type = std::vector<Cvector<thread_container_type>>;


        Cvector<hash_map_with_erase_type> m_boundaries_per_thread;
        std::vector<parallel::Cvector<parallel::thread_container<NodeID>>> m_containers;
        parallel::xxhash<NodeID> m_primary_hash;
        parallel::xxhash<NodeID> m_secondary_hash;
        uint32_t m_second_level_size;

        template <typename TGraph, typename TFunctor>
        void distribute_boundary_vertices(TGraph& graph, container_collection_type& containers, TFunctor&& bnd_func) {
                CLOCK_START;
                std::vector<std::future<std::vector<block_data_type>>> futures;
                futures.reserve(parallel::g_thread_pool.NumThreads());


                NodeID block_size = (NodeID) sqrt(graph.number_of_nodes());
                block_size = std::max(10 * block_size, 1000u);

                std::atomic<NodeID> offset(0);
                auto task_distribute = [this, &containers, &offset, block_size, &bnd_func, &graph] (uint32_t thread_id) {
                        std::vector<block_data_type> blocks_info(m_G.get_partition_count());
                        while (true) {
                                size_t cur_index = offset.load(std::memory_order_relaxed);

                                if (cur_index >= graph.number_of_nodes()) {
                                        break;
                                }

                                NodeID begin = offset.fetch_add(block_size, std::memory_order_relaxed);
                                NodeID end = begin + block_size;
                                end = end <= graph.number_of_nodes() ? end : graph.number_of_nodes();

                                if (begin >= graph.number_of_nodes()) {
                                        break;
                                }

                                for (NodeID node = begin; node != end; ++node) {
                                        PartitionID cur_part = graph.getPartitionIndex(node);
                                        int32_t num_external_neighbors = bnd_func(node, thread_id);

                                        ++blocks_info[cur_part].block_size;
                                        blocks_info[cur_part].block_weight += graph.getNodeWeight(node);
                                        if (num_external_neighbors > 0) {
                                                auto& container = containers[num_hash_table(node)][m_secondary_hash(node) % m_config.num_threads].get();
                                                container.concurrent_emplace_back(node, num_external_neighbors);
                                        }
                                }
                        }
                        return blocks_info;
                };

                for (size_t i = 0; i < parallel::g_thread_pool.NumThreads(); ++i) {
                        futures.push_back(parallel::g_thread_pool.Submit(i, task_distribute, i + 1));
                }

                m_blocks_info = task_distribute(0);
                std::for_each(futures.begin(), futures.end(), [&](auto& future){
                        std::vector<block_data_type> tmp_blocks_info(future.get());

                        for (size_t i = 0; i < m_blocks_info.size(); ++i) {
                                m_blocks_info[i].block_size += tmp_blocks_info[i].block_size;
                                m_blocks_info[i].block_weight += tmp_blocks_info[i].block_weight;
                        }
                });
                std::cout << std::endl;
                CLOCK_END("Distribute boundary vertices");
        }

        inline NodeID num_hash_table(NodeID vertex) const {
                return (NodeID) m_primary_hash(vertex) % m_config.num_threads;
        }

        inline int32_t external_neighbors(NodeID vertex) const {
                int32_t num_external_neighbors = 0;
                PartitionID cur_part = m_G.getPartitionIndex(vertex);
                forall_out_edges(m_G, e, vertex){
                        NodeID target = m_G.getEdgeTarget(e);
                        PartitionID part = m_G.getPartitionIndex(target);
                        if (cur_part != part) {
                                ++num_external_neighbors;
                        }
                } endfor
                return num_external_neighbors;
        }

        template <typename TGraph>
        inline int32_t is_boundary(TGraph& graph, NodeID vertex) const {
                int32_t num_external_neighbors = 0;
                PartitionID cur_part = graph.getPartitionIndex(vertex);
                forall_out_edges(graph, e, vertex) {
                        NodeID target = graph.getEdgeTarget(e);
                        PartitionID part = graph.getPartitionIndex(target);

                        if (cur_part == part) {
                                continue;
                        } else {
                                ++num_external_neighbors;
                                break;
                        }
                } endfor
                return num_external_neighbors;
        }
};

//using boundary_type = parallel::fast_sequential_boundary;
using boundary_type = parallel::fast_parallel_boundary;

}