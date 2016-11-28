#include "partition/uncoarsening/parallel_kway_graph_refinemet/kway_graph_refinement/kway_graph_refinement_commons.h"

namespace parallel {

class thread_data_factory {
public:
        using thread_data_refinement_core = parallel::thread_data_refinement_core;

        thread_data_factory(PartitionConfig& config,
                            graph_access& G,
                            complete_boundary& boundary,
                            boundary_starting_nodes& start_nodes,
                            int step_limit,
                            bool compute_touched_partitions,
                            std::unordered_map <PartitionID, PartitionID>& touched_blocks)
                :       m_G(G)
                ,       m_moved_idx(G.number_of_nodes())
                ,       m_parts_weights(config.k)
                ,       m_parts_sizes(config.k)
        {
                for (PartitionID block = 0; block < G.get_partition_count(); ++block) {
                        m_parts_weights[block].get().store(boundary.getBlockWeight(block), std::memory_order_relaxed);
                        m_parts_sizes[block].get().store(boundary.getBlockNoNodes(block), std::memory_order_relaxed);
                }

                thread_data.reserve(config.num_threads);

                for (uint32_t id = 0; id < config.num_threads; ++id) {
                        thread_data.emplace_back(id,
                                                 id,
                                                 config,
                                                 G,
                                                 boundary,
                                                 start_nodes,
                                                 step_limit,
                                                 m_moved_idx,
                                                 compute_touched_partitions,
                                                 touched_blocks,
                                                 parts_weights,
                                                 parts_size);
                }

        }


        void update_partitions_weights() {
                for (PartitionID block = 0; block < m_G.get_partition_count(); ++block) {
                        m_parts_weights[block].get().store(boundary.getBlockWeight(block), std::memory_order_relaxed);
                        m_parts_sizes[block].get().store(boundary.getBlockNoNodes(block), std::memory_order_relaxed);
                }
        }

        thread_data_refinement_core& get_thread_data(uint32_t id) {
                return thread_data[id];
        }

        const thread_data_refinement_core& get_thread_data(uint32_t id) const {
                return thread_data[id];
        }

private:
        Cvector<thread_data_refinement_core> thread_data;

        graph_access& m_G;

        std::vector <AtomicWrapper<bool>> m_moved_idx;
        Cvector <AtomicWrapper<NodeWeight>> m_parts_weights;
        Cvector <AtomicWrapper<NodeWeight>> m_parts_size;

};

class multitry_kway_fm {
public:
        using thread_data_refinement_core = parallel::thread_data_refinement_core;

        multitry_kway_fm(PartitionConfig& config,
                         graph_access& G,
                         complete_boundary& boundary,
                         boundary_starting_nodes& start_nodes,
                         int step_limit,
                         bool compute_touched_partitions)
                :       m_factory(config,
                                  graph,
                                  boundary,
                                  start_nodes,
                                  step_limit,
                                  compute_touched_partitions,
                                  touched_blocks)
        {}

        virtual ~multitry_kway_fm();

        int perform_refinement(PartitionConfig& config, graph_access& G,
                               complete_boundary& boundary, unsigned rounds,
                               bool init_neighbors, unsigned alpha);

        int perform_refinement_around_parts(PartitionConfig& config,
                                            graph_access& G,
                                            complete_boundary& boundary,
                                            bool init_neighbors,
                                            unsigned alpha,
                                            PartitionID& lhs, PartitionID& rhs,
                                            std::unordered_map <PartitionID, PartitionID>& touched_blocks);

private:
        thread_data_factory m_factory;

        int start_more_locallized_search(PartitionConfig& config, graph_access& G,
                                             complete_boundary& boundary,
                                             bool init_neighbors,
                                             bool compute_touched_blocks,
                                             std::unordered_map<PartitionID, PartitionID>& touched_blocks,
                                             std::vector<NodeID>& todolist);
};

}