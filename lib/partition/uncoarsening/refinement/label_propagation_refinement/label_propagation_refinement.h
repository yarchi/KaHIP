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
#include "data_structure/parallel/thread_pool.h"

#include "definitions.h"
#include "../refinement.h"

class label_propagation_refinement : public refinement {
public:
        label_propagation_refinement();
        virtual ~label_propagation_refinement();

        virtual EdgeWeight perform_refinement(PartitionConfig & config, 
                                              graph_access & G, 
                                              complete_boundary & boundary);

private:
        std::chrono::system_clock::time_point begin, end;
        enum class parallel_lp_type {
                parallel_queue, parallel
        };

        EdgeWeight sequential_label_propagation(PartitionConfig & config,
                                                graph_access & G,
                                                complete_boundary & boundary);

        EdgeWeight parallel_label_propagation(PartitionConfig & config,
                                              graph_access & G,
                                              complete_boundary & boundary);

        EdgeWeight parallel_label_propagation_with_queue(graph_access& G,
                                                         PartitionConfig& config,
                                                         parallel::TThreadPool& pool,
                                                         parallel::Cvector<parallel::AtomicWrapper<NodeWeight>>& cluster_sizes,
                                                         std::vector<std::vector<PartitionID>>& hash_maps,
                                                         std::vector<NodeID>& permutation);

        EdgeWeight parallel_label_propagation(graph_access& G,
                                              PartitionConfig& config,
                                              parallel::TThreadPool& pool,
                                              parallel::Cvector<parallel::AtomicWrapper<NodeWeight>>& cluster_sizes,
                                              std::vector<std::vector<PartitionID>>& hash_maps,
                                              std::vector<NodeID>& permutation);
};


#endif /* end of include guard: LABEL_PROPAGATION_REFINEMENT_R4XW141Y */
