/******************************************************************************
 * label_propagation_refinement.h  
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


#ifndef LABEL_PROPAGATION_REFINEMENT_R4XW141Y
#define LABEL_PROPAGATION_REFINEMENT_R4XW141Y

#include "data_structure/parallel/algorithm.h"
#include "data_structure/parallel/atomics.h"

#include "data_structure/parallel/pool_allocator.h"
#include "data_structure/parallel/thread_pool.h"

#include "definitions.h"
#include "../refinement.h"

#include <tbb/scalable_allocator.h>
#include <tbb/concurrent_queue.h>
#include <tbb/memory_pool.h>

#include <vector>

class label_propagation_refinement : public refinement {
public:
        label_propagation_refinement();
        virtual ~label_propagation_refinement();

        virtual EdgeWeight perform_refinement(PartitionConfig & config, 
                                              graph_access & G, 
                                              complete_boundary&);

private:
        //using Allocator = growt::PoolAllocator<NodeID>;
        using Allocator = tbb::scalable_allocator<NodeID>;
        using Block = std::vector<NodeID, Allocator>;
        using ConcurrentQueue = tbb::concurrent_queue<Block>;
        using Pair = std::pair<NodeID, NodeID>;

        Allocator m_block_allocator;

        inline uint32_t get_block_size(graph_access& G, PartitionConfig& config) const {
                uint32_t block_size = 1000;
                if (config.block_size_unit == BlockSizeUnit::NODES) {
                        block_size = std::max((uint32_t) sqrt(G.number_of_nodes()), block_size);
                }
                if (config.block_size_unit == BlockSizeUnit::EDGES) {
                        block_size = std::max((uint32_t) sqrt(G.number_of_edges()), block_size);
                }
                return block_size;
        }

        std::chrono::system_clock::time_point begin, end;
        EdgeWeight sequential_label_propagation(PartitionConfig & config,
                                                graph_access & G);

        EdgeWeight parallel_label_propagation(PartitionConfig & config,
                                              graph_access & G);

        EdgeWeight parallel_label_propagation_with_queue(graph_access& G,
                                                         PartitionConfig& config,
                                                         parallel::Cvector<parallel::AtomicWrapper<NodeWeight>>& cluster_sizes,
                                                         std::vector<std::vector<PartitionID>>& hash_maps,
                                                         std::vector<Pair>& permutation);

        EdgeWeight parallel_label_propagation(graph_access& G,
                                              PartitionConfig& config,
                                              parallel::Cvector<parallel::AtomicWrapper<NodeWeight>>& cluster_sizes,
                                              std::vector<std::vector<PartitionID>>& hash_maps,
                                              std::vector<Pair>& permutation);

        void seq_init_for_edge_unit(graph_access& G, const uint32_t block_size,
                                std::vector<Pair>& permutation,
                                std::vector<parallel::AtomicWrapper<NodeWeight>>& cluster_sizes,
                                std::unique_ptr<ConcurrentQueue>& queue);

        void par_init_for_edge_unit(graph_access& G, const uint32_t block_size,
                                    std::vector<Pair>& permutation,
                                    std::vector<parallel::AtomicWrapper<NodeWeight>>& cluster_sizes,
                                    std::unique_ptr<ConcurrentQueue>& queue);

        void init_for_node_unit(graph_access& G, const uint32_t block_size,
                                std::vector<Pair>& permutation,
                                std::vector<parallel::AtomicWrapper<NodeWeight>>& cluster_sizes,
                                std::unique_ptr<ConcurrentQueue>& queue);
};


#endif /* end of include guard: LABEL_PROPAGATION_REFINEMENT_R4XW141Y */
