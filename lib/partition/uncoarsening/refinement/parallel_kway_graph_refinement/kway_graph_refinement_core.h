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

        EdgeWeight apply_moves(uint32_t num_threads, Cvector <thread_data_refinement_core>& threads_data,
                               std::vector<NodeID>& reactivated_vertices) const;

private:
        static constexpr unsigned int sentinel = std::numeric_limits<unsigned int>::max();
        static constexpr int signed_sentinel = std::numeric_limits<int>::max();

        std::tuple<EdgeWeight, int, uint32_t> single_kway_refinement_round_internal(thread_data_refinement_core& td);

        void init_queue_with_boundary(thread_data_refinement_core& config,
                                      std::unique_ptr<refinement_pq>& queue);

        inline bool local_move_node(thread_data_refinement_core& config, NodeID node, PartitionID from, PartitionID& to,
                                    std::unique_ptr<refinement_pq>& queue, Gain gain);

        void unroll_relaxed_moves(thread_data_refinement_core& td,
                                  std::vector<NodeID>& transpositions,
                                  std::vector<PartitionID>& from_partitions,
                                  std::vector<Gain>& gains,
                                  int& cut_improvement) const;

        void relaxed_move_node_back(thread_data_refinement_core& td, NodeID node, PartitionID from,
                                    PartitionID to) const;

        inline bool relaxed_move_node(thread_data_refinement_core& td, NodeID node, PartitionID from,
                                      PartitionID to) const;

        uint32_t unroll_moves(thread_data_refinement_core& td, int min_cut_index) const;

        inline bool local_move_back_node(thread_data_refinement_core& td, NodeID node, PartitionID from,
                                         PartitionID to) const;

        EdgeWeight apply_moves(thread_data_refinement_core& td, std::vector<NodeID>& reacticated_vertices) const;
};
}
