/******************************************************************************
 * kway_graph_refinement_core.cpp
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 ******************************************************************************
 * Copyright (C) 2013-2015 Christian Schulz <christian.schulz@kit.edu>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 2 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#include <algorithm>

#include "data_structure/priority_queues/bucket_pq.h"
#include "data_structure/priority_queues/maxNodeHeap.h"
#include "tools/random_functions.h"
#include "uncoarsening/refinement/kway_graph_refinement/kway_stop_rule.h"
//#include "uncoarsening/refinement/kway_graph_refinement/quality_metrics.h"
#include "uncoarsening/refinement/parallel_kway_graph_refinement/kway_graph_refinement_core.h"

namespace parallel {
constexpr unsigned int kway_graph_refinement_core::sentinel;
constexpr int kway_graph_refinement_core::signed_sentinel;

kway_graph_refinement_core::kway_graph_refinement_core() {
}

kway_graph_refinement_core::~kway_graph_refinement_core() {

}

std::tuple<EdgeWeight, int, int>
kway_graph_refinement_core::single_kway_refinement_round_par(thread_data_refinement_core& td) const {
        std::unordered_map <PartitionID, PartitionID> touched_blocks;
        td.compute_touched_partitions = false;
        return single_kway_refinement_round_internal_par(td, touched_blocks);
}

std::tuple<EdgeWeight, int, int>
kway_graph_refinement_core::single_kway_refinement_round_par(thread_data_refinement_core& td,
                                                             std::unordered_map <PartitionID, PartitionID>& touched_blocks_ref) const {
        td.compute_touched_partitions = true;
        return single_kway_refinement_round_internal_par(td, touched_blocks_ref);
}

std::tuple<EdgeWeight, int, int>
kway_graph_refinement_core::single_kway_refinement_round_internal_par(thread_data_refinement_core& td,
                                                                      std::unordered_map <PartitionID, PartitionID>& touched_blocks) const {

        std::unique_ptr <refinement_pq> queue;
        if (td.config.use_bucket_queues) {
                EdgeWeight max_degree = td.G.getMaxDegree();
                queue = std::make_unique<bucket_pq>(max_degree);
        } else {
                queue = std::make_unique<maxNodeHeap>();
        }

        init_queue_with_boundary(td, queue);

        if (queue->empty()) {
                td.transpositions.push_back(sentinel);
                td.from_partitions.push_back(sentinel);
                td.to_partitions.push_back(sentinel);
                td.gains.push_back(signed_sentinel);
                return std::make_tuple(0, 0, -1);
        }

        int max_number_of_swaps = (int) (td.G.number_of_nodes());

        EdgeWeight cut = std::numeric_limits<int>::max() / 2; // so we dont need to compute the edge cut
        EdgeWeight initial_cut = cut;

        //roll forwards
        EdgeWeight best_cut = cut;
        int number_of_swaps = 0;
        int movements = 0;
        int overall_moved = 0;

        std::unique_ptr <kway_stop_rule> stopping_rule;
        switch (td.config.kway_stop_rule) {
                case KWAY_SIMPLE_STOP_RULE:
                        stopping_rule = std::make_unique<kway_simple_stop_rule>(td.config);
                        break;
                case KWAY_ADAPTIVE_STOP_RULE:
                        stopping_rule = std::make_unique<kway_adaptive_stop_rule>(td.config);
                        break;

        }

        int previously_moved = td.transpositions.size();
        int min_cut_index = previously_moved - 1;
        for (number_of_swaps = 0, movements = 0; movements < max_number_of_swaps; movements++, number_of_swaps++) {
                if (queue->empty()) {
                        break;
                }

                if (stopping_rule->search_should_stop(min_cut_index - previously_moved,
                                                      number_of_swaps, td.step_limit)) {
                        break;
                }

                Gain gain = queue->maxValue();
                NodeID node = queue->deleteMax();

                PartitionID from = td.get_local_partition(node);
#ifndef NDEBUG
                PartitionID maxgainer;
                EdgeWeight ext_degree;
                ASSERT_EQ(gain, td.compute_gain(node, from, maxgainer, ext_degree));
                ASSERT_TRUE(ext_degree > 0);
#endif

                bool successfull;
                int moved;
                PartitionID to;
                std::tie(successfull, moved) = local_move_node(td, node, from, to, queue);
                overall_moved += moved;

                if (successfull) {
                        cut -= gain;
                        stopping_rule->push_statistics(gain);

                        bool accept_equal = td.rnd.bit();
                        if (cut < best_cut || (cut == best_cut && accept_equal)) {
                                if (cut < best_cut)
                                        stopping_rule->reset_statistics();
                                best_cut = cut;
                                min_cut_index = previously_moved + number_of_swaps;
                        }
                        td.from_partitions.push_back(from);
                        td.to_partitions.push_back(to);
                        td.transpositions.push_back(node);
                        td.gains.push_back(gain);
                        ALWAYS_ASSERT(min_cut_index < (int64_t) td.transpositions.size());
                } else {
                        number_of_swaps--; //because it wasnt swaps
                }
        }

        overall_moved -= unroll_moves(td, min_cut_index);

        td.transpositions.push_back(sentinel);
        td.from_partitions.push_back(sentinel);
        td.to_partitions.push_back(sentinel);
        td.gains.push_back(signed_sentinel);
        return std::make_tuple(initial_cut - best_cut, overall_moved, min_cut_index);
}

uint32_t kway_graph_refinement_core::unroll_moves(thread_data_refinement_core& td, int min_cut_index) const {
        uint32_t unrolled_moves = 0;
        // unrolled_moves <  td.transpositions.size() - (min_cut_index + 1)
        while (unrolled_moves + min_cut_index + 1 < td.transpositions.size()) {
                size_t index = td.transpositions.size() - 1 - unrolled_moves;
                NodeID node = td.transpositions[index];
                PartitionID from = td.from_partitions[index];
                PartitionID to = td.to_partitions[index];
                local_move_back_node(td, node, from, to);

                ++unrolled_moves;
        }
        return unrolled_moves;
}

std::pair<EdgeWeight, uint32_t> kway_graph_refinement_core::apply_moves(std::vector<thread_data_refinement_core>& threads_data) const {

        uint32_t overall_moved = 0;
        EdgeWeight overall_gain = 0;

        moved_nodes_hash_set moved_nodes(moved_nodes_hash_set::get_max_size_to_fit_l1());
        std::vector <NodeID> moved_nodes_vec;
        for (size_t id = 0; id < threads_data.size(); ++id) {
                overall_gain += apply_moves(threads_data[id], moved_nodes, moved_nodes_vec);
                overall_moved += moved_nodes_vec.size();
        }
        return std::make_pair(overall_gain, overall_moved);
}

EdgeWeight kway_graph_refinement_core::apply_moves(thread_data_refinement_core& td,
                                                 moved_nodes_hash_set& moved_nodes,
                                                 std::vector <NodeID>& moved_nodes_vec) const {
        ALWAYS_ASSERT(td.transpositions.size() == td.from_partitions.size());
        ALWAYS_ASSERT(td.transpositions.size() == td.to_partitions.size());
        ALWAYS_ASSERT(td.transpositions.size() == td.gains.size());
        auto min_cut_iter = td.min_cut_indices.begin();

        moved_nodes_vec.clear();
        moved_nodes_vec.reserve(td.transpositions.size());
        EdgeWeight cut_improvement = 0;

        for (int index = 0; index < (int) td.transpositions.size(); ++index) {
                int min_cut_index = min_cut_iter->first;
                uint32_t next_index = min_cut_iter->second;
                ++min_cut_iter;

                if (min_cut_index == -1) {
                        index = next_index;
                        continue;
                }
                int cut = 0;
                int start_index = index;
                while (index <= min_cut_index) {
                        NodeID node = td.transpositions[index];
                        PartitionID from = td.from_partitions[index];
                        PartitionID to = td.to_partitions[index];
                        EdgeWeight gain = td.gains[index];


                        // check if any nodes where moved by other threads,
                        // if yes then stop moving
                        bool move_node = true;
                        forall_out_edges(td.G, e, node) {
                                NodeID target = td.G.getEdgeTarget(e);
                                if (moved_nodes.contains(target)) {
                                        move_node = false;
                                        unroll_relaxed_moves(td, start_index, index, cut_improvement);
                                        return cut_improvement;
                                }
                        }
                        endfor

                        // move node
                        if (move_node) {
                                bool success = relaxed_move_node(td, node, from, to);
                                if (success) {
                                        moved_nodes_vec.push_back(node);
                                        cut_improvement += gain;
                                        cut += gain;
                                } else {
                                        unroll_relaxed_moves(td, start_index, index, cut_improvement);
                                        return cut_improvement;
                                }
                        }
                        ++index;
                }

                index = next_index;
        }

        for (auto node : moved_nodes_vec) {
                moved_nodes.insert(node);
        }
        return cut_improvement;
}

void kway_graph_refinement_core::init_queue_with_boundary(thread_data_refinement_core& td,
                                                          std::unique_ptr <refinement_pq>& queue) const {
        if (td.config.permutation_during_refinement == PERMUTATION_QUALITY_FAST) {
                random_functions::permutate_vector_fast(td.start_nodes, false);
        } else if (td.config.permutation_during_refinement == PERMUTATION_QUALITY_GOOD) {
                random_functions::permutate_vector_good(td.start_nodes, false);
        }

        for (unsigned int i = 0; i < td.start_nodes.size(); i++) {
                NodeID node = td.start_nodes[i];

                bool expected = false;
                if (td.moved_idx[node].compare_exchange_strong(expected, true, std::memory_order_relaxed)) {
                        PartitionID max_gainer;
                        EdgeWeight ext_degree;
                        //compute gain
                        PartitionID from = td.get_local_partition(node);
                        Gain gain = td.compute_gain(node, from, max_gainer, ext_degree);
                        queue->insert(node, gain);
                }
        }
}

inline bool kway_graph_refinement_core::kway_graph_refinement_core::relaxed_move_node(thread_data_refinement_core& td,
                                                                                      NodeID node,
                                                                                      PartitionID from,
                                                                                      PartitionID to) const {
        ASSERT_TRUE(td.boundary.assert_bnodes_in_boundaries());
        ASSERT_TRUE(td.boundary.assert_boundaries_are_bnodes());

        NodeWeight this_nodes_weight = td.G.getNodeWeight(node);

        if (td.boundary.getBlockWeight(to) + this_nodes_weight >= td.config.upper_bound_partition)
                return false;

        if (td.boundary.getBlockNoNodes(from) == 1) // assure that no block gets accidentally empty
                return false;

        td.G.setPartitionIndex(node, to);

        boundary_pair pair;
        pair.k = td.config.k;
        pair.lhs = from;
        pair.rhs = to;

        td.boundary.postMovedBoundaryNodeUpdates(node, &pair, true, true);

        td.boundary.setBlockNoNodes(from, td.boundary.getBlockNoNodes(from) - 1);
        td.boundary.setBlockNoNodes(to, td.boundary.getBlockNoNodes(to) + 1);
        td.boundary.setBlockWeight(from, td.boundary.getBlockWeight(from) - this_nodes_weight);
        td.boundary.setBlockWeight(to, td.boundary.getBlockWeight(to) + this_nodes_weight);

        ASSERT_TRUE(td.boundary.assert_bnodes_in_boundaries());
        ASSERT_TRUE(td.boundary.assert_boundaries_are_bnodes());

        return true;
}

void kway_graph_refinement_core::unroll_relaxed_moves(thread_data_refinement_core& td, int start, int end,
                                                      int& cut_improvement) const {
        for (int index = end - 1; index >= start; --index) {
                NodeID node = td.transpositions[index];
                PartitionID from = td.from_partitions[index];
                PartitionID to = td.to_partitions[index];
                cut_improvement -= td.gains[index];
                relaxed_move_node_back(td, node, from, to);
        }
}

void kway_graph_refinement_core::relaxed_move_node_back(thread_data_refinement_core& td, NodeID node, PartitionID from,
                                                        PartitionID to) const {
        ALWAYS_ASSERT(td.G.getPartitionIndex(node) == to);

        td.G.setPartitionIndex(node, from);

        boundary_pair pair;
        pair.k   = td.config.k;
        pair.lhs = from;
        pair.rhs = to;

        //update all boundaries
        td.boundary.postMovedBoundaryNodeUpdates(node, &pair, true, true);

        NodeWeight this_nodes_weight = td.G.getNodeWeight(node);
        td.boundary.setBlockNoNodes(from, td.boundary.getBlockNoNodes(from) + 1);
        td.boundary.setBlockNoNodes(to, td.boundary.getBlockNoNodes(to) - 1);
        td.boundary.setBlockWeight(from, td.boundary.getBlockWeight(from) + this_nodes_weight);
        td.boundary.setBlockWeight(to, td.boundary.getBlockWeight(to) - this_nodes_weight);
}

inline bool kway_graph_refinement_core::local_move_back_node(thread_data_refinement_core& td,
                                                                        NodeID node,
                                                                        PartitionID from,
                                                                        PartitionID to) const {
        td.nodes_partitions[node] = from;
        NodeWeight this_nodes_weight = td.G.getNodeWeight(node);

        td.parts_weights[from].get().fetch_add(this_nodes_weight, std::memory_order_relaxed);
        td.parts_weights[to].get().fetch_sub(this_nodes_weight, std::memory_order_relaxed);
        td.parts_sizes[to].get().fetch_sub(1, std::memory_order_relaxed);
        td.parts_sizes[from].get().fetch_add(1, std::memory_order_relaxed);

        return true;
};

inline std::pair<bool, int>
kway_graph_refinement_core::kway_graph_refinement_core::local_move_node(thread_data_refinement_core& td,
                                                                        NodeID node,
                                                                        PartitionID from,
                                                                        PartitionID& to,
                                                                        std::unique_ptr <refinement_pq>& queue) const {

        EdgeWeight node_ext_deg;

        td.compute_gain(node, from, to, node_ext_deg);

        NodeWeight this_nodes_weight = td.G.getNodeWeight(node);
        NodeWeight part_weight = td.parts_weights[to].get().load(std::memory_order_relaxed);

        if (td.parts_sizes[from].get().load(std::memory_order_relaxed) == 1) {
                return std::make_pair(false, 0);
        }

        do {
                if (part_weight + this_nodes_weight >= td.config.upper_bound_partition) {
                        return std::make_pair(false, 0);
                }
        } while (!td.parts_weights[to].get().compare_exchange_weak(part_weight,
                                                                   part_weight + this_nodes_weight,
                                                                   std::memory_order_relaxed));

        td.nodes_partitions[node] = to;
        td.parts_weights[from].get().fetch_sub(this_nodes_weight, std::memory_order_relaxed);
        td.parts_sizes[to].get().fetch_add(1, std::memory_order_relaxed);
        td.parts_sizes[from].get().fetch_sub(1, std::memory_order_relaxed);

        int moved = 0;
        //update gain of neighbors / the boundaries have allready been updated
        forall_out_edges(td.G, e, node)
        {
                NodeID target = td.G.getEdgeTarget(e);
                PartitionID targets_to;
                EdgeWeight ext_degree; // the local external degree
                PartitionID target_from = td.get_local_partition(target);

                Gain gain = td.compute_gain(target, target_from, targets_to, ext_degree);

                if (queue->contains(target)) {
                        assert(td.moved_idx[target].load(std::memory_order_relaxed));
                        if (ext_degree > 0) {
                                queue->changeKey(target, gain);
                        } else {
                                queue->deleteNode(target);
                        }
                } else {
                        if (ext_degree > 0) {
                                bool expected = false;
                                if (td.moved_idx[target].compare_exchange_strong(expected, true,
                                                                                 std::memory_order_relaxed)) {
                                        queue->insert(target, gain);
                                        ++moved;
                                }
                        }
                }
        }
        endfor

        return std::make_pair(true, moved);
}
}
