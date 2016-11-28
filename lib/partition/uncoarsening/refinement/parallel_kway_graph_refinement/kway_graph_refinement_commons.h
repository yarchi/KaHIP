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
        using nodes_partitions_hash_table = parallel::HashMap<NodeID, PartitionID, parallel::xxhash<NodeID>, true>;

        PartitionConfig& config;
        graph_access& G;
        complete_boundary& boundary;
        boundary_starting_nodes& start_nodes;
        int step_limit;
        std::vector <AtomicWrapper<bool>>& moved_idx;
        bool compute_touched_partitions;
        std::unordered_map <PartitionID, PartitionID>& touched_blocks;
        Cvector <AtomicWrapper<NodeWeight>>& parts_weights;
        Cvector <AtomicWrapper<NodeWeight>>& parts_sizes;
        int upper_bound_gain_improvement;

        nodes_partitions_hash_table nodes_partitions;
        std::vector<std::pair<PartitionID, uint32_t>> min_cut_indices;
        std::vector <NodeID> transpositions;
        std::vector <PartitionID> from_partitions;
        std::vector <PartitionID> to_partitions;
        std::vector <EdgeWeight> gains;
        std::vector <NodeID> moved;


        thread_data_refinement_core(uint32_t _id,
                                    uint32_t _seed,
                                    PartitionConfig& _config,
                                    graph_access& _G,
                                    complete_boundary& _boundary,
                                    boundary_starting_nodes& _start_nodes,
                                    int _step_limit,
                                    std::vector <AtomicWrapper<bool>>& _moved_idx,
                                    bool _compute_touched_partitions,
                                    std::unordered_map <PartitionID, PartitionID>& _touched_blocks,
                                    Cvector <AtomicWrapper<NodeWeight>>& _parts_weights,
                                    Cvector <AtomicWrapper<NodeWeight>>& _parts_sizes)
                :       parallel::thread_config(_id, _seed)
                ,       config(_config)
                ,       G(_G)
                ,       boundary(_boundary)
                ,       start_nodes(_start_nodes)
                ,       step_limit(_step_limit)
                ,       moved_idx(_moved_idx)
                ,       compute_touched_partitions(_compute_touched_partitions)
                ,       touched_blocks(_touched_blocks)
                ,       parts_weights(_parts_weights)
                ,       parts_sizes(_parts_sizes)
                ,       upper_bound_gain_improvement(0)
                ,       nodes_partitions(nodes_partitions_hash_table::get_max_size_to_fit_l1())
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
        }

        thread_data_refinement_core(const thread_data_refinement_core& td) = default;
//                :       parallel::thread_config(td)
//                ,       config(td.config)
//                ,       G(td.G)
//                ,       boundary(td.boundary)
//                ,       start_nodes(td.start_nodes)
//                ,       step_limit(td.step_limit)
//                ,       moved_idx(td.moved_idx)
//                ,       compute_touched_partitions(td.compute_touched_partitions)
//                ,       touched_blocks(td.touched_blocks)
//                ,       parts_weights(td.parts_weights)
//                ,       parts_sizes(td.parts_sizes)
//                ,       upper_bound_gain_improvement(td.upper_bound_gain_improvement)
//                ,       nodes_partitions(td.nodes_partitions)
//                ,       transpositions(td.transpositions)
//                ,       from_partitions(td.from_partitions)
//                ,       to_partitions(td.to_partitions)
//                ,       gains(td.gains)
//                ,       min_cut_indices(td.min_cut_indices)
//                ,       m_local_degrees(td.m_local_degrees)
//                ,       m_round(td.m_round)
//        {
//                m_local_degrees.resize(config.k);
//
//                // needed for the computation of internal and external degrees
//                m_round = 0;
//        }

        thread_data_refinement_core(thread_data_refinement_core&& td) = default;
//                :       parallel::thread_config(td)
//                ,       config(td.config)
//                ,       G(td.G)
//                ,       boundary(td.boundary)
//                ,       start_nodes(td.start_nodes)
//                ,       step_limit(td.step_limit)
//                ,       moved_idx(td.moved_idx)
//                ,       compute_touched_partitions(td.compute_touched_partitions)
//                ,       touched_blocks(td.touched_blocks)
//                ,       parts_weights(td.parts_weights)
//                ,       parts_sizes(td.parts_sizes)
//                ,       upper_bound_gain_improvement(td.upper_bound_gain_improvement)
//                ,       nodes_partitions(std::move(td.nodes_partitions))
//                ,       transpositions(std::move(td.transpositions))
//                ,       from_partitions(td.from_partitions)
//                ,       to_partitions(td.to_partitions)
//                ,       gains(td.gains)
//                ,       min_cut_indices(td.min_cut_indices)
//        {}

        thread_data_refinement_core& operator=(const thread_data_refinement_core&) = delete;
        thread_data_refinement_core& operator=(thread_data_refinement_core&&) = delete;

        inline PartitionID get_local_partition(NodeID node) const {
                return nodes_partitions.contains(node) ? nodes_partitions[node] : G.getPartitionIndex(node);
        }

        void clear() {
                nodes_partitions.clear();
                min_cut_indices.clear();
                transpositions.clear();
                from_partitions.clear();
                to_partitions.clear();
                gains.clear();

                for (auto node : moved) {
                        moved_idx[node].store(false, std::memory_order_relaxed);
                }
                moved.clear();

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