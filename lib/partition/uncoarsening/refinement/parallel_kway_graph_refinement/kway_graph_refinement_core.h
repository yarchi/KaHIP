#pragma once

#include <unordered_map>
#include <vector>
#include <memory>

#include "definitions.h"
#include "data_structure/priority_queues/priority_queue_interface.h"
#include "uncoarsening/refinement/kway_graph_refinement/kway_graph_refinement_commons.h"
#include "uncoarsening/refinement/parallel_kway_graph_refinement/kway_graph_refinement_commons.h"

#include <tbb/concurrent_queue.h>

namespace parallel {

class kway_graph_refinement_core {
public:
        kway_graph_refinement_core();

        virtual ~kway_graph_refinement_core();

        std::tuple<EdgeWeight, int, uint32_t> single_kway_refinement_round(thread_data_refinement_core& td);

        std::pair<EdgeWeight, uint32_t> apply_moves(uint32_t num_threads,
                                                    Cvector <thread_data_refinement_core>& threads_data,
                                                    bool compute_touched_partitions,
                                                    std::unordered_map<PartitionID, PartitionID>& touched_blocks,
                                                    std::vector<NodeID>& reactivated_vertices) const;

        std::pair<EdgeWeight, uint32_t> apply_moves(Cvector <thread_data_refinement_core>& threads_data,
                                                    bool compute_touched_partitions,
                                                    std::unordered_map<PartitionID, PartitionID>& touched_blocks,
                                                    std::vector<NodeID>& reactivated_vertices,
                                                    tbb::concurrent_queue<uint32_t>& finished_threads,
                                                    std::vector<std::future<bool>>& futures,
                                                    bool& is_more_that_5percent_moved) const;

private:
        using moved_nodes_hash_map = parallel::hash_map_with_erase<NodeID, std::pair<uint32_t, PartitionID>>;
        using moved_hash_set = parallel::hash_set<NodeID>; // moved by local search

        static constexpr unsigned int sentinel = std::numeric_limits<unsigned int>::max();
        static constexpr int signed_sentinel = std::numeric_limits<int>::max();

        std::tuple<EdgeWeight, int, uint32_t> single_kway_refinement_round_internal(thread_data_refinement_core& td);

        std::pair<EdgeWeight, uint32_t> local_search_from_one_node(thread_data_refinement_core& td,
                                                                   moved_nodes_hash_map& moved_nodes,
                                                                   NodeID start_node,
                                                                   uint32_t max_number_of_swaps,
                                                                   bool compute_touched_partitions,
                                                                   std::unordered_map<PartitionID, PartitionID>& touched_blocks) const;

        std::pair<EdgeWeight, uint32_t> gain_recalculation(thread_data_refinement_core& td,
                                                           moved_nodes_hash_map& moved_nodes, int start, int end,
                                                           bool compute_touched_partitions,
                                                           std::unordered_map<PartitionID, PartitionID>& touched_blocks) const;

        void init_queue_with_boundary(thread_data_refinement_core& config,
                                      std::unique_ptr<refinement_pq>& queue);

        inline bool local_move_node(thread_data_refinement_core& config, NodeID node, PartitionID from, PartitionID& to,
                                    std::unique_ptr<refinement_pq>& queue, Gain gain);

        void unroll_relaxed_moves(thread_data_refinement_core& td, moved_nodes_hash_map& moved_nodes,
                                  int start, int end, int& cut_improvement) const;

        void relaxed_move_node_back(thread_data_refinement_core& td, NodeID node, PartitionID from,
                                    PartitionID to) const;

        inline bool relaxed_move_node(thread_data_refinement_core& td, NodeID node, PartitionID from,
                                      PartitionID to) const;

        void unroll_moves(thread_data_refinement_core& td, int min_cut_index) const;

        inline bool local_move_back_node(thread_data_refinement_core& td, NodeID node, PartitionID from,
                                         PartitionID to) const;

        EdgeWeight apply_moves(thread_data_refinement_core& config, moved_nodes_hash_map& moved_nodes,
                               bool compute_touched_partitions,
                               std::unordered_map<PartitionID, PartitionID>& touched_blocks,
                               std::vector<NodeID>& reactivated_vertices) const;

        bool is_moved(moved_nodes_hash_map& moved_nodes, NodeID node, uint32_t thread_id) const;

        std::pair<EdgeWeight, uint32_t>
        apply_moves_with_time_stamp(Cvector<thread_data_refinement_core>& threads_data,
                                    bool compute_touched_partitions,
                                    std::unordered_map<PartitionID, PartitionID>& touched_blocks,
                                    std::vector<std::future<uint32_t>>& futures) const;

        inline bool move_node(thread_data_refinement_core& td,
                              moved_hash_set& moved,
                              NodeID node,
                              std::unique_ptr<refinement_pq>& queue,
                              kway_graph_refinement_commons* commons) const;
};
}
