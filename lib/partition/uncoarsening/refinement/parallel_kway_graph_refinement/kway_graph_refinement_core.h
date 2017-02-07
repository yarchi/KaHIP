#pragma once

#include <unordered_map>
#include <vector>
#include <memory>

#include "definitions.h"
#include "data_structure/priority_queues/priority_queue_interface.h"
#include "uncoarsening/refinement/parallel_kway_graph_refinement/kway_graph_refinement_commons.h"

namespace parallel {

class kway_graph_refinement_core {
public:
        kway_graph_refinement_core();

        virtual ~kway_graph_refinement_core();

        std::tuple<EdgeWeight, int, int, uint32_t> single_kway_refinement_round(thread_data_refinement_core& td) const;

        std::pair <EdgeWeight, uint32_t> apply_moves(Cvector<thread_data_refinement_core>& threads_data,
                                                     bool compute_touched_partitions,
                                                     std::unordered_map <PartitionID, PartitionID>& touched_blocks,
                                                     std::vector<std::future<uint32_t>>& futures) const;

private:
        using moved_nodes_hash_set = parallel::hash_set<NodeID>;

        static constexpr unsigned int sentinel = std::numeric_limits < unsigned int>::max();
        static constexpr int signed_sentinel = std::numeric_limits<int>::max();

        std::tuple<EdgeWeight, int, int, uint32_t> single_kway_refinement_round_internal(thread_data_refinement_core& td) const;

        void init_queue_with_boundary(thread_data_refinement_core& config,
                                      std::unique_ptr <refinement_pq>& queue) const;

        inline std::pair<bool, int> local_move_node(thread_data_refinement_core& config,
                                                    NodeID node,
                                                    PartitionID from,
                                                    PartitionID& to,
                                                    std::unique_ptr <refinement_pq>& queue) const;

        void unroll_relaxed_moves(thread_data_refinement_core& td, int start, int end, int& cut_improvement) const;
        void relaxed_move_node_back(thread_data_refinement_core& td, NodeID node, PartitionID from,
                                    PartitionID to) const;

        inline bool relaxed_move_node(thread_data_refinement_core& td, NodeID node, PartitionID from,
                                      PartitionID to) const;

        uint32_t unroll_moves(thread_data_refinement_core& td, int min_cut_index) const;

        inline bool local_move_back_node(thread_data_refinement_core& td, NodeID node, PartitionID from,
                                         PartitionID to) const;

        EdgeWeight apply_moves(thread_data_refinement_core& config, moved_nodes_hash_set& moved_nodes,
                               std::vector <NodeID>& moved_nodes_vec,
                               bool compute_touched_partitions,
                               std::unordered_map <PartitionID, PartitionID>& touched_blocks) const;
};
}
