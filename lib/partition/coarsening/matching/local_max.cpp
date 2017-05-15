#include "coarsening/matching/local_max.h"
#include "data_structure/parallel/thread_pool.h"
#include "data_structure/parallel/time.h"

#include <tbb/concurrent_queue.h>

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
        CLOCK_START;
        std::vector<NodeID> vertex_permutation;
        std::vector<std::pair<NodeID, uint32_t>> max_neighbours;

        vertex_permutation.reserve(G.number_of_nodes());
        edge_matching.reserve(G.number_of_nodes());
        max_neighbours.reserve(G.number_of_nodes());

        for (size_t i = 0; i < G.number_of_nodes(); ++i) {
                vertex_permutation.push_back(i);
                edge_matching.push_back(i);
                max_neighbours.emplace_back(i, 0);
        }
        std::random_shuffle(vertex_permutation.begin(), vertex_permutation.end());

        // need two queues to save max_neighbour in array
        std::unique_ptr<std::queue<NodeID>> node_queue = std::make_unique<std::queue<NodeID>>();
        std::unique_ptr<std::queue<NodeID>> node_queue_next = std::make_unique<std::queue<NodeID>>();
        for (NodeID node : vertex_permutation) {
                node_queue->push(node);
        }
        random rnd(partition_config.seed);
        CLOCK_END("Coarsening: Matching: Init");

        CLOCK_START_N;
        uint32_t round = 0;
        while (!node_queue->empty()) {
                while (!node_queue->empty()) {
                        NodeID node = node_queue->front();
                        node_queue->pop();

                        // already mathed
                        if (edge_matching[node] != node) {
                                continue;
                        }

                        if (max_neighbours[node].second != round) {
                                max_neighbours[node].first = node;
                                max_neighbours[node].second = round;
                        }
                        NodeID max_neighbour = max_neighbours[node].first;
                        if (max_neighbour == node) {
                                max_neighbour = find_max_neighbour_sequential(node, G, partition_config, edge_matching,
                                                                              rnd);
                                max_neighbours[node].first = max_neighbour;
                        }

                        if (max_neighbour == m_none) {
                                continue;
                        }

                        if (max_neighbours[max_neighbour].second != round) {
                                max_neighbours[max_neighbour].first = max_neighbour;
                                max_neighbours[max_neighbour].second = round;
                        }
                        NodeID max_neighbour_neighbour = max_neighbours[max_neighbour].first;
                        if (max_neighbour_neighbour == max_neighbour) {
                                max_neighbour_neighbour = find_max_neighbour_sequential(max_neighbour, G,
                                                                                        partition_config,
                                                                                        edge_matching,
                                                                                        rnd);
                                max_neighbours[max_neighbour].first = max_neighbour_neighbour;
                        }


                        if (max_neighbour_neighbour == node) {
                                // match
                                edge_matching[node] = max_neighbour;
                                edge_matching[max_neighbour] = node;
                        } else {
                                node_queue_next->push(node);
                        }
                }
                ++round;
                std::swap(node_queue, node_queue_next);
        }
        CLOCK_END("Coarsening: Matching: Main");
}

void local_max_matching::parallel_match(const PartitionConfig& partition_config,
                                        graph_access& G,
                                        Matching& edge_matching,
                                        CoarseMapping& mapping,
                                        NodeID& no_of_coarse_vertices,
                                        NodePermutationMap& permutation) {

        // init
        //auto vertex_mark = Get_parallel_vector<AtomicWrapper<bool>>(parallel::g_thread_pool, false);
        CLOCK_START;
        uint32_t num_threads = partition_config.num_threads;
        std::vector<AtomicWrapper<int>> vertex_mark(G.number_of_nodes(), MatchingPhases::NOT_STARTED);
        std::vector<std::pair<AtomicWrapper<NodeID>, AtomicWrapper<uint32_t>>> max_neighbours(G.number_of_nodes());

        Cvector<random> randoms;
        randoms.reserve(num_threads);
        for (uint32_t id = 0; id < num_threads; ++id) {
                randoms.emplace_back(partition_config.seed + id);
        }

        std::vector<NodeID> vertex_permutation;
        vertex_permutation.reserve(G.number_of_nodes());
        edge_matching.reserve(G.number_of_nodes());
        for (size_t i = 0; i < G.number_of_nodes(); ++i) {
                vertex_permutation.push_back(i);
                edge_matching.push_back(i);
        }
        std::random_shuffle(vertex_permutation.begin(), vertex_permutation.end());

        const size_t block_size = 100;
        tbb::concurrent_queue<block_type> node_queue;
        block_type block;
        block.reserve(block_size);
        for (auto node : vertex_permutation) {
                block.push_back(node);
                if (block.size() >= block_size) {
                        node_queue.push(std::move(block));
                        block.clear();
                        block.reserve(block_size);
                }
        }

        CLOCK_END("Coarsening: Matching: Init");
        uint32_t round = 1;

        auto task = [&](uint32_t id) {
                block_type cur_block;
                block_type next_block;
                next_block.reserve(block_size);

                auto& rnd = randoms[id].get();

                while (node_queue.try_pop(block)) {
                        for (NodeID node : block) {
                                if (vertex_mark[node].load(std::memory_order_relaxed) == MatchingPhases::MATCHED) {
                                        continue;
                                }

                                int mark = MatchingPhases::NOT_STARTED;
                                if (vertex_mark[node].compare_exchange_strong(mark, MatchingPhases::STARTED,
                                                                              std::memory_order_acquire)) {

                                        uint32_t elem_round = round - 1;
                                        if (max_neighbours[node].second.compare_exchange_strong(elem_round, round,
                                                                                                std::memory_order_relaxed)) {
                                                // find best neighbour for target
                                                NodeID max_neighbour = find_max_neighbour_parallel(node, G,
                                                                                                   partition_config,
                                                                                                   vertex_mark, rnd);

                                                max_neighbours[node].first.store(max_neighbour,
                                                                                 std::memory_order_relaxed);
                                        }

                                        vertex_mark[node].store(MatchingPhases::FOUND_LOCAL_MAX,
                                                                std::memory_order_release);
                                } else {
                                        while (vertex_mark[node].load(std::memory_order_acquire) ==
                                               MatchingPhases::STARTED);
                                }

                                NodeID max_neighbour = max_neighbours[node].first.load(std::memory_order_relaxed);
                                if (max_neighbour == m_none) {
                                        continue;
                                }

                                mark = MatchingPhases::NOT_STARTED;
                                if (vertex_mark[max_neighbour].compare_exchange_strong(mark, MatchingPhases::STARTED,
                                                                                       std::memory_order_acq_rel)) {

                                        uint32_t elem_round = round - 1;
                                        if (max_neighbours[max_neighbour].second.compare_exchange_strong(elem_round,
                                                                                                         round,
                                                                                                         std::memory_order_relaxed)) {
                                                // find best neighbour of best neighbour
                                                NodeID max_target_neighbour = find_max_neighbour_parallel(max_neighbour, G,
                                                                                                          partition_config,
                                                                                                          vertex_mark, rnd);

                                                max_neighbours[max_neighbour].first.store(max_target_neighbour,
                                                                                          std::memory_order_relaxed);
                                        }

                                        vertex_mark[max_neighbour].store(MatchingPhases::FOUND_LOCAL_MAX,
                                                                         std::memory_order_release);
                                } else {
                                        while (vertex_mark[max_neighbour].load(std::memory_order_acquire) ==
                                               MatchingPhases::STARTED);
                                }

                                if (max_neighbours[max_neighbour].first.load(std::memory_order_relaxed) == node) {
                                        // match edge
                                        if (node < max_neighbour) {
                                                edge_matching[node] = max_neighbour;
                                                edge_matching[max_neighbour] = node;
                                                vertex_mark[node].store(MatchingPhases::MATCHED,
                                                                        std::memory_order_release);
                                                vertex_mark[max_neighbour].store(MatchingPhases::MATCHED,
                                                                                 std::memory_order_release);
                                        }
                                } else {
                                        // put in the other queue vertex !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                        vertex_mark[node].store(MatchingPhases::NOT_STARTED, std::memory_order_release);
                                        next_block.push_back(node);
                                        if (next_block.size() >= block_size) {
                                                node_queue.push(std::move(next_block));
                                                next_block.clear();
                                                next_block.reserve(block_size);
                                        }
                                }
                        }
                }
        };

        CLOCK_START_N;
        // add loop for rounds with two queues
        std::vector<std::future<void>> futures;
        futures.reserve(num_threads);

        for (uint32_t id = 1; id < num_threads; ++id) {
                futures.push_back(g_thread_pool.Submit(task, id));
        }

        task(0);
        std::for_each(futures.begin(), futures.end(), [](auto& future) {
                future.get();
        });
        ++round;

        CLOCK_END("Coarsening: Matching: Main");
}

NodeID
local_max_matching::find_max_neighbour_sequential(NodeID node, graph_access& G, const PartitionConfig& partition_config,
                                                  Matching& edge_matching, random& rnd) const {
        NodeWeight node_weight = G.getNodeWeight(node);
        EdgeWeight max_weight = 0;
        NodeID max_target = m_none;

        forall_out_edges(G, e, node) {
                NodeID target = G.getEdgeTarget(e);
                EdgeWeight edge_weight = G.getEdgeWeight(e);
                NodeWeight coarser_weight = G.getNodeWeight(target) + node_weight;

                if ((edge_weight > max_weight || (edge_weight == max_weight && rnd.bit())) &&
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

                        max_target = target;
                        max_weight = edge_weight;
                }
        } endfor
        return max_target;
}

NodeID
local_max_matching::find_max_neighbour_parallel(NodeID node, graph_access& G, const PartitionConfig& partition_config,
                                                std::vector<AtomicWrapper<int>>& vertex_mark, random& rnd) const {
        NodeWeight node_weight = G.getNodeWeight(node);
        EdgeWeight max_weight = 0;
        NodeID max_target = m_none;
        forall_out_edges(G, e, node) {
                NodeID target = G.getEdgeTarget(e);
                EdgeWeight edge_weight = G.getEdgeWeight(e);
                NodeWeight coarser_weight = G.getNodeWeight(target) + node_weight;

                if ((edge_weight > max_weight || (edge_weight == max_weight && rnd.bit())) &&
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
                        max_weight = edge_weight;
                }
        } endfor
        return max_target;
}

}