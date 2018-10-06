#include "coarsening/matching/local_max.h"
#include "data_structure/parallel/thread_pool.h"
#include "data_structure/parallel/time.h"

#include <tbb/concurrent_queue.h>

#include <parallel/algorithm>

namespace parallel {

void local_max_matching::match(const PartitionConfig& partition_config,
                               graph_access& G,
                               Matching& edge_matching,
                               CoarseMapping& mapping,
                               NodeID& no_of_coarse_vertices,
                               NodePermutationMap& permutation) {
        switch (partition_config.matching_type) {
                case MATCHING_SEQUENTIAL_LOCAL_MAX:
                        sequential_match(partition_config, G, edge_matching, mapping, no_of_coarse_vertices,
                                         permutation);
                        break;
                case MATCHING_PARALLEL_LOCAL_MAX:
                        // with queue is slower since we do not process vertices in increasing order of their degree
                        // and shuffle them
                        parallel_match(partition_config, G, edge_matching, mapping, no_of_coarse_vertices, permutation);
                        break;
                default:
                        std::cout << "Incorrect matching type expected sequential local max or parallel local max"
                                  << std::endl;
                        abort();
        }
}

void local_max_matching::sequential_match(const PartitionConfig& partition_config,
                                          graph_access& G,
                                          Matching& edge_matching,
                                          CoarseMapping& mapping,
                                          NodeID& no_of_coarse_vertices,
                                          NodePermutationMap& permutation) {
        std::cout << "MAX VERTEX WEIGHT = " << partition_config.max_vertex_weight << std::endl;
        CLOCK_START;
        std::vector<std::pair<NodeID, uint32_t>> max_neighbours;

        permutation.clear();
        permutation.reserve(G.number_of_nodes());
        edge_matching.reserve(G.number_of_nodes());
        max_neighbours.reserve(G.number_of_nodes());

        for (size_t i = 0; i < G.number_of_nodes(); ++i) {
                permutation.push_back(i);
                edge_matching.push_back(i);
                max_neighbours.emplace_back(i, 0);
        }
        std::random_shuffle(permutation.begin(), permutation.end());

        // need two queues to save max_neighbour in array
        std::unique_ptr<std::queue<NodeID>> node_queue = std::make_unique<std::queue<NodeID>>();
        std::unique_ptr<std::queue<NodeID>> node_queue_next = std::make_unique<std::queue<NodeID>>();
        for (NodeID node : permutation) {
                node_queue->push(node);
        }
        random rnd(partition_config.seed);
        CLOCK_END("Coarsening: Matching: Init");

        CLOCK_START_N;
        uint32_t round = 1;
        uint32_t threshold = (uint32_t) (2.0 * G.number_of_nodes() / 3.0);
        uint32_t remaining_vertices = G.number_of_nodes();
        while (!node_queue->empty() && round < m_max_round + 1 && remaining_vertices > threshold) {
                std::cout << round - 1 << std::endl;
                std::cout << remaining_vertices << std::endl;
                CLOCK_START;
                while (!node_queue->empty()) {
                        NodeID node = node_queue->front();
                        node_queue->pop();

                        // already mathed
                        if (edge_matching[node] != node) {
                                continue;
                        }

                        NodeID max_neighbour = find_max_neighbour_sequential(node, round, G, partition_config,
                                                                             max_neighbours, edge_matching, rnd);

                        if (max_neighbour == m_none) {
                                continue;
                        }


                        NodeID max_neighbour_neighbour = find_max_neighbour_sequential(max_neighbour, round, G,
                                                                                       partition_config, max_neighbours,
                                                                                       edge_matching, rnd);

                        if (max_neighbour_neighbour == node) {
                                // match
                                edge_matching[node] = max_neighbour;
                                edge_matching[max_neighbour] = node;
                                --remaining_vertices;
                        } else {
                                node_queue_next->push(node);
                        }
                }
                ++round;
                std::swap(node_queue, node_queue_next);
                CLOCK_END("Round time");
        }
        std::cout << remaining_vertices << std::endl;
        CLOCK_END("Coarsening: Matching: Main");

        CLOCK_START_N;
        no_of_coarse_vertices = 0;
        mapping.resize(G.number_of_edges());
        for (NodeID i = 0; i < permutation.size(); ++i) {
                NodeID node = permutation[i];
                if (node <= edge_matching[node]) {
                        mapping[edge_matching[node]] = no_of_coarse_vertices;
                        mapping[node] = no_of_coarse_vertices;
                        ++no_of_coarse_vertices;
                }
        }
        std::cout << no_of_coarse_vertices << std::endl;
        CLOCK_END("Coarsening: Matching: Remap");
}

void local_max_matching::parallel_match_with_queue(const PartitionConfig& partition_config,
                                        graph_access& G,
                                        Matching& edge_matching,
                                        CoarseMapping& mapping,
                                        NodeID& no_of_coarse_vertices,
                                        NodePermutationMap& permutation) {

        // init
        CLOCK_START;
        uint32_t num_threads = partition_config.num_threads;

        parallel::ParallelVector<AtomicWrapper<int>> vertex_mark(G.number_of_nodes());
        parallel::ParallelVector<atomic_pair_type> max_neighbours(G.number_of_nodes());

        parallel::parallel_for_index(0u, G.number_of_nodes(), [&](NodeID node) {
                vertex_mark[node] = MatchingPhases::NOT_STARTED;
                new (max_neighbours.begin() + node) atomic_pair_type(0, 0);
        });

        permutation.clear();
        permutation.reserve(G.number_of_nodes());
        edge_matching.reserve(G.number_of_nodes());
        for (size_t i = 0; i < G.number_of_nodes(); ++i) {
                permutation.push_back(i);
                edge_matching.push_back(i);
        }

        {
                CLOCK_START;
                parallel::random_shuffle(permutation.begin(), permutation.end(), partition_config.num_threads);
                CLOCK_END("Shuffle");
        }

        Cvector<random> randoms;
        randoms.reserve(num_threads);
        for (uint32_t id = 0; id < num_threads; ++id) {
                randoms.emplace_back(partition_config.seed + id);
        }

        const size_t block_size_edges = std::max((uint32_t) sqrt(G.number_of_edges()), 1000u);
        std::unique_ptr<tbb::concurrent_queue<block_type>> node_queue = std::make_unique<tbb::concurrent_queue<block_type>>();
        std::unique_ptr<tbb::concurrent_queue<block_type>> node_queue_next = std::make_unique<tbb::concurrent_queue<block_type>>();

        block_type block;
        block.reserve(block_size_edges);
        size_t cur_block_size = 0;
        for (auto node : permutation) {
                block.push_back(node);
                cur_block_size += G.getNodeDegree(node);
                if (cur_block_size >= block_size_edges) {
                        node_queue->push(std::move(block));
                        block.clear();
                        block.reserve(block_size_edges);
                        cur_block_size = 0;
                }
        }
        if (!block.empty()) {
                node_queue->push(std::move(block));
        }

        CLOCK_END("Coarsening: Matching: Init");
        uint32_t round = 1;

        auto task = [&](uint32_t id) {
                size_t this_thread_block_size = 0;
                block_type cur_block;
                block_type next_block;
                next_block.reserve(block_size_edges);

                auto& rnd = randoms[id].get();
                NodeID matched_vertices = 0;

                while (node_queue->try_pop(cur_block)) {
                        for (NodeID node : cur_block) {
                                if (vertex_mark[node].load(std::memory_order_relaxed) == MatchingPhases::MATCHED) {
                                        continue;
                                }

                                NodeID max_neighbour = find_max_neighbour_parallel(node, round, G, partition_config,
                                                                                   max_neighbours, vertex_mark, rnd);
                                if (max_neighbour == m_none) {
                                        continue;
                                }

                                NodeID max_neighbour_neighbour = find_max_neighbour_parallel(max_neighbour, round, G,
                                                                                             partition_config,
                                                                                             max_neighbours,
                                                                                             vertex_mark, rnd);

                                if (max_neighbour_neighbour == node) {
                                        // match edge
                                        if (node < max_neighbour) {
                                                edge_matching[node] = max_neighbour;
                                                edge_matching[max_neighbour] = node;
                                                ++matched_vertices;
                                                vertex_mark[node].store(MatchingPhases::MATCHED,
                                                                        std::memory_order_release);
                                                vertex_mark[max_neighbour].store(MatchingPhases::MATCHED,
                                                                                 std::memory_order_release);
                                        }
                                } else {
                                        // put in the other queue vertex
                                        next_block.push_back(node);
                                        this_thread_block_size += G.getNodeDegree(node);
                                        if (this_thread_block_size >= block_size_edges) {
                                                node_queue_next->push(std::move(next_block));
                                                next_block.clear();
                                                next_block.reserve(block_size_edges);
                                                this_thread_block_size = 0;
                                        }
                                }
                        }
                }
                if (!next_block.empty()) {
                        node_queue_next->push(std::move(next_block));
                }
                return matched_vertices;
        };

        CLOCK_START_N;
        // add loop for rounds with two queues
        std::vector<std::future<NodeID>> futures;
        futures.reserve(num_threads);
        NodeID threshold = (NodeID) (2.0 * G.number_of_nodes() / 3.0);
        NodeID coarse_vertices = G.number_of_nodes();
        while (!node_queue->empty() && round < m_max_round + 1 &&  coarse_vertices > threshold) {
                CLOCK_START;
                std::cout << round - 1 << std::endl;
                std::cout << coarse_vertices << std::endl;
                futures.clear();
                for (uint32_t id = 1; id < num_threads; ++id) {
                        futures.push_back(g_thread_pool.Submit(id - 1, task, id));
                }

                coarse_vertices -= task(0);
                std::for_each(futures.begin(), futures.end(), [&coarse_vertices](auto& future) {
                        coarse_vertices -= future.get();
                });
                node_queue.swap(node_queue_next);
                ++round;
                CLOCK_END("Round time");
        }
        std::cout << coarse_vertices << std::endl;

        CLOCK_END("Coarsening: Matching: Main");

        CLOCK_START_N;
        no_of_coarse_vertices = coarse_vertices;
        remap_matching(partition_config, G, edge_matching, mapping, no_of_coarse_vertices, permutation);
        CLOCK_END("Coarsening: Matching: Remap");
}

void local_max_matching::parallel_match(const PartitionConfig& partition_config,
                                        graph_access& G,
                                        Matching& edge_matching,
                                        CoarseMapping& mapping,
                                        NodeID& no_of_coarse_vertices,
                                        NodePermutationMap& permutation) {

        // init
        CLOCK_START;
        uint32_t num_threads = partition_config.num_threads;

        parallel::ParallelVector<AtomicWrapper<int>> vertex_mark(G.number_of_nodes());
        parallel::ParallelVector<atomic_pair_type> max_neighbours(G.number_of_nodes());

        parallel::parallel_for_index(0u, G.number_of_nodes(), [&](NodeID node) {
                vertex_mark[node] = MatchingPhases::NOT_STARTED;
                new (max_neighbours.begin() + node) atomic_pair_type(0, 0);
        });

        permutation.clear();
        permutation.reserve(G.number_of_nodes());
        edge_matching.reserve(G.number_of_nodes());
        for (size_t i = 0; i < G.number_of_nodes(); ++i) {
                permutation.push_back(i);
                edge_matching.push_back(i);
        }

        {
                CLOCK_START;
                parallel::random_shuffle(permutation.begin(), permutation.end(), partition_config.num_threads);
                CLOCK_END("Shuffle");
        }

        Cvector<random> randoms;
        randoms.reserve(num_threads);
        for (uint32_t id = 0; id < num_threads; ++id) {
                randoms.emplace_back(partition_config.seed + id);
        }
        CLOCK_END("Coarsening: Matching: Init");

        uint32_t round = 1;
        std::atomic<NodeID> offset(0);
        const NodeID block_size = std::max<NodeID>(sqrt(G.number_of_nodes()), 1000);
        auto task = [&](uint32_t id) {
                auto& rnd = randoms[id].get();
                NodeID matched_vertices = 0;

                while (true) {
                        NodeID begin = offset.fetch_add(block_size, std::memory_order_relaxed);
                        NodeID end = begin + block_size;
                        end = end <= G.number_of_nodes() ? end : G.number_of_nodes();

                        if (begin >= G.number_of_nodes()) {
                                break;
                        }

                        for (NodeID node = begin; node != end; ++node) {
                                if (vertex_mark[node].load(std::memory_order_relaxed) == MatchingPhases::MATCHED) {
                                        continue;
                                }

                                NodeID max_neighbour = find_max_neighbour_parallel(node, round, G, partition_config,
                                                                                   max_neighbours,
                                                                                   vertex_mark, rnd);
                                if (max_neighbour == m_none) {
                                        continue;
                                }

                                NodeID max_neighbour_neighbour = find_max_neighbour_parallel(max_neighbour, round, G,
                                                                                             partition_config,
                                                                                             max_neighbours,
                                                                                             vertex_mark, rnd);

                                if (max_neighbour_neighbour == node) {
                                        // match edge
                                        if (node < max_neighbour) {
                                                ALWAYS_ASSERT(edge_matching[node] == node);
                                                ALWAYS_ASSERT(edge_matching[max_neighbour] == max_neighbour);
                                                edge_matching[node] = max_neighbour;
                                                edge_matching[max_neighbour] = node;
                                                ++matched_vertices;
                                                vertex_mark[node].store(MatchingPhases::MATCHED,
                                                                        std::memory_order_release);
                                                vertex_mark[max_neighbour].store(MatchingPhases::MATCHED,
                                                                                 std::memory_order_release);
                                        }
                                }
                        }
                }
                return matched_vertices;
        };

        CLOCK_START_N;
        std::vector<std::future<NodeID>> futures;
        futures.reserve(num_threads);
        NodeID threshold = (NodeID) (2.0 * G.number_of_nodes() / 3.0);
        NodeID coarse_vertices = G.number_of_nodes();
        while (round < m_max_round + 1 &&  coarse_vertices > threshold) {
                CLOCK_START;
                offset.store(0, std::memory_order_acquire);
                std::cout << round - 1 << std::endl;
                std::cout << coarse_vertices << std::endl;
                futures.clear();
                for (uint32_t id = 1; id < num_threads; ++id) {
                        futures.push_back(g_thread_pool.Submit(id - 1, task, id));
                }

                coarse_vertices -= task(0);
                std::for_each(futures.begin(), futures.end(), [&coarse_vertices](auto& future) {
                        coarse_vertices -= future.get();
                });
                ++round;
                CLOCK_END("Round time");
        }
        std::cout << coarse_vertices << std::endl;

        CLOCK_END("Coarsening: Matching: Main");

        CLOCK_START_N;
        no_of_coarse_vertices = coarse_vertices;
        remap_matching(partition_config, G, edge_matching, mapping, no_of_coarse_vertices, permutation);
        CLOCK_END("Coarsening: Matching: Remap");
}

void local_max_matching::remap_matching(const PartitionConfig& partition_config, graph_access& G,
                                        Matching& edge_matching, CoarseMapping& mapping, NodeID& no_of_coarse_vertices,
                                        NodePermutationMap& permutation) {
        parallel::ParallelVector<NodeID> aux_edge_matching(G.number_of_nodes());
        parallel::parallel_for_index(NodeID(0), G.number_of_nodes(), [&](NodeID node) {
                NodeID matched = edge_matching[node];
                ALWAYS_ASSERT(node == edge_matching[matched]);
                if (node <= matched) {
                        aux_edge_matching[node] = 1;
                } else {
                        aux_edge_matching[node] = 0;
                }
        });

        std::vector<NodeID> n(G.number_of_nodes());
        parallel::partial_sum(aux_edge_matching.begin(), aux_edge_matching.end(), n.begin(),
                              partition_config.num_threads);

        ALWAYS_ASSERT(no_of_coarse_vertices == n.back());

        mapping.resize(G.number_of_nodes());
        parallel::parallel_for_index(NodeID(0), G.number_of_nodes(), [&](NodeID node) {
                if (node <= edge_matching[node]) {
                        NodeID coarse_node = n[node] - 1;
                        mapping[edge_matching[node]] = coarse_node;
                        mapping[node] = coarse_node;
                }
        });
}

NodeID local_max_matching::find_max_neighbour_parallel(NodeID node, const uint32_t round, graph_access& G,
                                                       const PartitionConfig& partition_config,
                                                       ParallelVector<atomic_pair_type>& max_neighbours,
                                                       ParallelVector<AtomicWrapper<int>>& vertex_mark,
                                                       random& rnd) const {
        if (max_neighbours[node].first.load(std::memory_order_relaxed) == m_none) {
                return m_none;
        }

        int mark = MatchingPhases::NOT_STARTED;

        if (max_neighbours[node].second.load(std::memory_order_acquire) == round) {
                while (vertex_mark[node].load(std::memory_order_acquire) == MatchingPhases::STARTED);
                return max_neighbours[node].first.load(std::memory_order_relaxed);
        }

        if (vertex_mark[node].compare_exchange_strong(mark, MatchingPhases::STARTED, std::memory_order_acquire)) {

                // region open -- this region is evaluated ONLY when vertex_mark[node] == MatchingPhases::STARTED by ONLY one thread
                uint32_t elem_round = round - 1;
                if (max_neighbours[node].second.compare_exchange_strong(elem_round, round, std::memory_order_relaxed)) {
                        // find best neighbour for target
                        NodeID max_neighbour = find_max_neighbour_parallel(node, G, partition_config, vertex_mark, rnd);

                        max_neighbours[node].first.store(max_neighbour, std::memory_order_relaxed);
                }
                // region close

                int mark = MatchingPhases::STARTED;
                vertex_mark[node].compare_exchange_strong(mark, MatchingPhases::NOT_STARTED, std::memory_order_release);
        } else {
                while (vertex_mark[node].load(std::memory_order_acquire) == MatchingPhases::STARTED);
        }
        return max_neighbours[node].first.load(std::memory_order_relaxed);
}

NodeID local_max_matching::find_max_neighbour_sequential(NodeID node, const uint32_t round, graph_access& G,
                                                  const PartitionConfig& partition_config,
                                                  std::vector<std::pair<NodeID, uint32_t>>& max_neighbours,
                                                  Matching& edge_matching, random& rnd) const {
        NodeID max_neighbour = max_neighbours[node].first;

        // recaluclate max_neighbour even if max_neighbour from previous iteration is NOT matched
        // because if there are a lot of edges with the same rating then only few of them will match
        if (max_neighbour == m_none || round == max_neighbours[node].second) {
                return max_neighbour;
        }
        max_neighbours[node].second = round;

        NodeWeight node_weight = G.getNodeWeight(node);
        EdgeRatingType max_rating = 0;
        max_neighbour = m_none;

        forall_out_edges(G, e, node) {
                NodeID target = G.getEdgeTarget(e);
                EdgeRatingType edge_rating = G.getEdgeRating(e);
                NodeWeight coarser_weight = G.getNodeWeight(target) + node_weight;

                if ((edge_rating > max_rating || (edge_rating == max_rating && rnd.bit())) &&
                    edge_matching[target] == target &&
                    coarser_weight <= partition_config.max_vertex_weight) {

                        if (partition_config.graph_allready_partitioned &&
                            G.getPartitionIndex(node) != G.getPartitionIndex(target)) {
                                continue;
                        }

                        if (partition_config.combine &&
                            G.getSecondPartitionIndex(node) != G.getSecondPartitionIndex(target)) {
                                continue;
                        }

                        max_neighbour = target;
                        max_rating = edge_rating;
                }
        } endfor

        max_neighbours[node].first = max_neighbour;
        return max_neighbour;
}

NodeID
local_max_matching::find_max_neighbour_parallel(NodeID node, graph_access& G, const PartitionConfig& partition_config,
                                                ParallelVector<AtomicWrapper<int>>& vertex_mark, random& rnd) const {
        NodeWeight node_weight = G.getNodeWeight(node);
        EdgeRatingType max_rating = 0;
        NodeID max_target = m_none;
        forall_out_edges(G, e, node) {
                NodeID target = G.getEdgeTarget(e);
                EdgeRatingType edge_rating = G.getEdgeRating(e);
                NodeWeight coarser_weight = G.getNodeWeight(target) + node_weight;

                if ((edge_rating > max_rating || (edge_rating == max_rating && rnd.bit())) &&
                    vertex_mark[target].load(std::memory_order_relaxed) != MatchingPhases::MATCHED &&
                    coarser_weight <= partition_config.max_vertex_weight) {

                        if (partition_config.graph_allready_partitioned &&
                            G.getPartitionIndex(node) != G.getPartitionIndex(target)) {
                                continue;
                        }

                        if (partition_config.combine &&
                            G.getSecondPartitionIndex(node) != G.getSecondPartitionIndex(target)) {
                                continue;
                        }

                        max_target = target;
                        max_rating = edge_rating;
                }
        } endfor
        return max_target;
}

}