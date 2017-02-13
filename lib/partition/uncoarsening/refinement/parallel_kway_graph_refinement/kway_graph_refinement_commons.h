#pragma once

#include <vector>

#include "data_structure/graph_access.h"
#include "data_structure/parallel/algorithm.h"
#include "data_structure/parallel/atomics.h"
#include "data_structure/parallel/hash_table.h"
#include "data_structure/parallel/spin_lock.h"
#include "data_structure/parallel/thread_config.h"
#include "definitions.h"
#include "partition/partition_config.h"
#include "partition/uncoarsening/refinement/quotient_graph_refinement/complete_boundary.h"

namespace parallel {
class thread_data_refinement_core : public parallel::thread_config {
public:
        using nodes_partitions_hash_table = parallel::hash_map<NodeID, PartitionID>;
        //using nodes_partitions_hash_table = std::unordered_map<NodeID, PartitionID>;

        // global data
        PartitionConfig& config;
        graph_access& G;
        complete_boundary& boundary;
        int step_limit;
        std::vector <AtomicWrapper<bool>>& moved_idx;
        Cvector <AtomicWrapper<NodeWeight>>& parts_weights;
        Cvector <AtomicWrapper<NodeWeight>>& parts_sizes;
        Cvector <AtomicWrapper<int>>& moved_count;
        int upper_bound_gain_improvement;
        AtomicWrapper<uint32_t>& time_stamp;

        // local thread data
        boundary_starting_nodes start_nodes;
        nodes_partitions_hash_table nodes_partitions;
        std::vector<std::pair<int, int>> min_cut_indices;
        std::vector<NodeID> transpositions;
        std::vector<PartitionID> from_partitions;
        std::vector<PartitionID> to_partitions;
        std::vector<EdgeWeight> gains;
        std::vector<NodeID> moved;
        std::vector<uint32_t> time_stamps;

        // local statistics about time in all iterations
        double total_thread_time;
        uint32_t tried_movements;
        uint32_t accepted_movements;
        double total_thread_try_move_time;
        double total_thread_accepted_move_time;
        double total_thread_unroll_move_time;
        uint32_t scaned_neighbours;
        double time_move_nodes;
        uint32_t transpositions_size;
        int unperformed_gain;

        thread_data_refinement_core(uint32_t _id,
                                    uint32_t _seed,
                                    PartitionConfig& _config,
                                    graph_access& _G,
                                    complete_boundary& _boundary,
                                    std::vector <AtomicWrapper<bool>>& _moved_idx,
                                    Cvector <AtomicWrapper<NodeWeight>>& _parts_weights,
                                    Cvector <AtomicWrapper<NodeWeight>>& _parts_sizes,
                                    Cvector <AtomicWrapper<int>>& _moved_count,
                                    AtomicWrapper<uint32_t>& _moved_idx_counter,
                                    AtomicWrapper<uint32_t>& _time_stamp)
                :       parallel::thread_config(_id, _seed)
                ,       config(_config)
                ,       G(_G)
                ,       boundary(_boundary)
                ,       step_limit(0)
                ,       moved_idx(_moved_idx)
                ,       parts_weights(_parts_weights)
                ,       parts_sizes(_parts_sizes)
                ,       moved_count(_moved_count)
                ,       upper_bound_gain_improvement(0)
                ,       time_stamp(_time_stamp)
                ,       nodes_partitions(131072)
                ,       total_thread_time(0.0)
                ,       tried_movements(0)
                ,       accepted_movements(0)
                ,       total_thread_try_move_time(0.0)
                ,       total_thread_unroll_move_time(0.0)
                ,       scaned_neighbours(0)
                ,       time_move_nodes(0.0)
                ,       transpositions_size(0)
                ,       unperformed_gain(0)
                ,       m_moved_idx_counter(_moved_idx_counter)
        {
                m_local_degrees.resize(config.k);

                // needed for the computation of internal and external degrees
                m_round = 0;

                min_cut_indices.reserve(100);
                transpositions.reserve(100);
                from_partitions.reserve(100);
                to_partitions.reserve(100);
                gains.reserve(100);
                moved.reserve(100);
                start_nodes.reserve(100);
        }

        thread_data_refinement_core(const thread_data_refinement_core& td) = default;
        thread_data_refinement_core(thread_data_refinement_core&& td) = default;

        thread_data_refinement_core& operator=(const thread_data_refinement_core&) = delete;
        thread_data_refinement_core& operator=(thread_data_refinement_core&&) = delete;

        inline PartitionID get_local_partition(NodeID node) {
                return nodes_partitions.contains(node) ? nodes_partitions[node] : G.getPartitionIndex(node);
                //return nodes_partitions.find(node) != nodes_partitions.end() ? nodes_partitions[node] : G.getPartitionIndex(node);
        }

        void reset_thread_data() {
                for (auto node : moved) {
                        moved_idx[node].store(false, std::memory_order_relaxed);
                }
                m_moved_idx_counter.fetch_add(1, std::memory_order_release);

                moved.clear();

                nodes_partitions.clear();
                //nodes_partitions.reserve(nodes_partitions_hash_table::get_max_size_to_fit_l1());
                nodes_partitions.reserve(131072);

                min_cut_indices.clear();
                transpositions.clear();
                from_partitions.clear();
                to_partitions.clear();
                gains.clear();
                start_nodes.clear();
                time_stamps.clear();

                while (!is_all_data_reseted());
        }

        inline Gain compute_gain(NodeID node, PartitionID from, PartitionID& to, EdgeWeight& ext_degree) {
                //for all incident partitions compute gain
                //return max gain and "to" partition
                EdgeWeight max_degree = 0;
                to = INVALID_PARTITION;

                m_round++;//can become zero again
                forall_out_edges(G, e, node)
                {
                        NodeID target = G.getEdgeTarget(e);
                        PartitionID target_partition = get_local_partition(target);

                        if (m_local_degrees[target_partition].round == m_round) {
                                m_local_degrees[target_partition].local_degree += G.getEdgeWeight(e);
                        } else {
                                m_local_degrees[target_partition].local_degree = G.getEdgeWeight(e);
                                m_local_degrees[target_partition].round = m_round;
                        }


                        if (target_partition != from && m_local_degrees[target_partition].local_degree >= max_degree) {
                                if (m_local_degrees[target_partition].local_degree > max_degree) {
                                        max_degree = m_local_degrees[target_partition].local_degree;
                                        to = target_partition;
                                } else {
                                        //break ties randomly
                                        bool accept = rnd.bit();
                                        if (accept) {
                                                max_degree = m_local_degrees[target_partition].local_degree;
                                                to = target_partition;
                                        }
                                }
                        }
                }
                endfor

                //std::cout << "Node " << node << ", degree " << G.getNodeDegree(node) << ", cur part " << from << ", best part: " << to << std::endl;

//                if (to == 4294967295) {
//                        forall_out_edges(G, e, node) {
//                                NodeID target = G.getEdgeTarget(e);
//                                PartitionID target_partition = get_local_partition(target);
//                                std::cout << "target " << target << ", part: " << target_partition << std::endl;
//                        } endfor
//                }

                if (to != INVALID_PARTITION) {
                        ext_degree = max_degree;
                } else {
                        ext_degree = 0;
                }

                if (m_local_degrees[from].round != m_round) {
                        m_local_degrees[from].local_degree = 0;
                }

                return max_degree - m_local_degrees[from].local_degree;
        }

private:
        AtomicWrapper<uint32_t>& m_moved_idx_counter;

        inline bool is_all_data_reseted() const {
                return m_moved_idx_counter.load(std::memory_order_acquire) == config.num_threads;
        }

        //for efficient computation of internal and external degrees
        struct round_struct {
                round_struct()
                        :       round(0)
                        ,       local_degree(0)
                {}

                uint32_t round;
                EdgeWeight local_degree;
        };

        std::vector<round_struct> m_local_degrees;
        uint32_t m_round;
};
}