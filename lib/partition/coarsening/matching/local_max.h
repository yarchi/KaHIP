#pragma once

#include "data_structure/parallel/algorithm.h"
#include "data_structure/parallel/atomics.h"
#include "data_structure/parallel/random.h"

#include "matching.h"

namespace parallel {
class local_max_matching : public matching {
public:
        void match(const PartitionConfig& partition_config,
                   graph_access& G,
                   Matching& edge_matching,
                   CoarseMapping& mapping,
                   NodeID& no_of_coarse_vertices,
                   NodePermutationMap& permutation) override;

private:
        static constexpr NodeID m_none = std::numeric_limits<NodeID>::max();
        static constexpr uint32_t m_max_round = 2;

        using block_type = std::vector<NodeID>;


        void parallel_match(const PartitionConfig& partition_config,
                            graph_access& G,
                            Matching& edge_matching,
                            CoarseMapping& mapping,
                            NodeID& no_of_coarse_vertices,
                            NodePermutationMap& permutation);

        void sequential_match(const PartitionConfig& partition_config,
                              graph_access& G,
                              Matching& edge_matching,
                              CoarseMapping& mapping,
                              NodeID& no_of_coarse_vertices,
                              NodePermutationMap& permutation);

        NodeID find_max_neighbour_sequential(NodeID node, graph_access& G, const PartitionConfig& partition_config,
                                             Matching& vertex_mark, random& rnd) const;

        NodeID find_max_neighbour_parallel(NodeID node, const uint32_t round, graph_access& G,
                                           const PartitionConfig& partition_config,
                                           std::vector<std::pair<AtomicWrapper<NodeID>, AtomicWrapper<uint32_t>>>& max_neighbours,
                                           std::vector<AtomicWrapper<int>>& vertex_mark, random& rnd) const;

        NodeID find_max_neighbour_parallel(NodeID node, graph_access& G, const PartitionConfig& partition_config,
                                           std::vector<AtomicWrapper<int>>& vertex_mark, random& rnd) const;

        enum MatchingPhases {
                NOT_STARTED = 0,
                STARTED,
                FOUND_LOCAL_MAX,
                MATCHED
        };
};
}