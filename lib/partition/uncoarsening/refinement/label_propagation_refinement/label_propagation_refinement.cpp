/******************************************************************************
 * label_propagation_refinement.cpp 
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

#include "label_propagation_refinement.h"
#include "partition/coarsening/clustering/node_ordering.h"
#include "tools/random_functions.h"

//#define TBB_PREVIEW_GLOBAL_CONTROL 1
//#include <tbb/global_control.h>

#include <parallel/algorithm>

#include <omp.h>

#include <chrono>
#include <random>

#include "ittnotify.h"

using namespace parallel;

label_propagation_refinement::label_propagation_refinement()
        :       m_block_allocator(2)
{}

label_propagation_refinement::~label_propagation_refinement() {

}

EdgeWeight label_propagation_refinement::perform_refinement(PartitionConfig& config, graph_access& G,
                                                            complete_boundary& boundary) {

        //if (!config.parallel_local_search || G.number_of_nodes() < 10000000) {
        if (!config.parallel_local_search) {
                auto begin = std::chrono::high_resolution_clock::now();
                auto res = sequential_label_propagation(config, G, boundary);
                auto end = std::chrono::high_resolution_clock::now();

                std::cout << "Sequential lp:\t" << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()
                          << std::endl;
                return res;
        } else {
                // init memory pool for blocks
                m_block_allocator.reset();
                double size = config.block_size_unit == BlockSizeUnit::EDGES ? G.number_of_edges() : G.number_of_nodes();
                m_block_allocator.init((2 * std::ceil(size / config.block_size) + config.num_threads) * config.block_size * sizeof(NodeID));

                auto begin = std::chrono::high_resolution_clock::now();
                auto res = parallel_label_propagation(config, G, boundary);
                auto end = std::chrono::high_resolution_clock::now();

                std::cout << "Parallel lp:\t" << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()
                          << std::endl;

                return res;
        }
}

EdgeWeight label_propagation_refinement::sequential_label_propagation(PartitionConfig & partition_config,
                                                                      graph_access & G,
                                                                      complete_boundary & boundary) {
        auto begin = std::chrono::high_resolution_clock::now();
        __itt_resume();
        NodeWeight block_upperbound = partition_config.upper_bound_partition;

        // in this case the _matching paramter is not used 
        // coarse_mappng stores cluster id and the mapping (it is identical)
        std::vector<PartitionID> hash_map(partition_config.k,0);
        std::vector<NodeID> permutation(G.number_of_nodes());
        std::vector<NodeWeight> cluster_sizes(partition_config.k, 0);

        auto b = std::chrono::high_resolution_clock::now();
        node_ordering n_ordering;
        n_ordering.order_nodes(partition_config, G, permutation);
        auto e = std::chrono::high_resolution_clock::now();
        std::cout << "Sequential init of permutations lp:\t" << std::chrono::duration_cast<std::chrono::microseconds>(e - b).count()
                  << std::endl;

        std::queue< NodeID > * Q             = new std::queue< NodeID >();
        std::queue< NodeID > * next_Q        = new std::queue< NodeID >();
        std::vector<bool> * Q_contained      = new std::vector<bool>(G.number_of_nodes(), false);
        std::vector<bool> * next_Q_contained = new std::vector<bool> (G.number_of_nodes(), false);
        forall_nodes(G, node) {
                cluster_sizes[G.getPartitionIndex(node)] += G.getNodeWeight(node);
                Q->push(permutation[node]);
        } endfor
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "Init sequential lp:\t" << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()
                  << std::endl;

        begin = std::chrono::high_resolution_clock::now();

        int parts_size = partition_config.k * 2;
        std::vector<PartitionID> partIDs;
        partIDs.reserve(parts_size);
        unsigned int change_counter = 0;

        for( int j = 0; j < partition_config.label_iterations_refinement; j++) {
                while( !Q->empty() ) {
                        NodeID node = Q->front();
                        Q->pop();
                        (*Q_contained)[node] = false;

                        if (G.getNodeDegree(node) <= parts_size) {
                                //now move the node to the cluster that is most common in the neighborhood
                                forall_out_edges(G, e, node) {
                                        NodeID target = G.getEdgeTarget(e);
                                        PartitionID part = G.getPartitionIndex(target);
                                        hash_map[part] += G.getEdgeWeight(e);
                                        // ONLY FOR UNIT WEIGHTS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                                        //++hash_map[part];
                                        // ONLY FOR UNIT WEIGHTS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                                        partIDs.push_back(part);
                                } endfor
                        } else {
                                //now move the node to the cluster that is most common in the neighborhood
                                forall_out_edges(G, e, node) {
                                        NodeID target = G.getEdgeTarget(e);
                                        hash_map[G.getPartitionIndex(target)] += G.getEdgeWeight(e);
                                } endfor
                        }

                        //second sweep for finding max and resetting array
                        PartitionID max_block = G.getPartitionIndex(node);
                        PartitionID my_block  = G.getPartitionIndex(node);
                        PartitionID max_value = 0;

                        if (G.getNodeDegree(node) <= parts_size) {
                                for (auto cur_block : partIDs) {
                                        if (hash_map[cur_block] == 0) {
                                                continue;
                                        }

                                        PartitionID cur_value     = hash_map[cur_block];

                                        if((cur_value > max_value  || (cur_value == max_value && random_functions::nextBool()))
                                           && (cluster_sizes[cur_block] + G.getNodeWeight(node) < block_upperbound || (cur_block == my_block && cluster_sizes[my_block] <= partition_config.upper_bound_partition)))
                                        {
                                                max_value = cur_value;
                                                max_block = cur_block;
                                        }

                                        hash_map[cur_block] = 0;
                                }
                                partIDs.clear();
                        } else {
                                forall_out_edges(G, e, node) {
                                        NodeID target             = G.getEdgeTarget(e);
                                        PartitionID cur_block     = G.getPartitionIndex(target);
                                        if (hash_map[cur_block] == 0) {
                                                continue;
                                        }

                                        PartitionID cur_value     = hash_map[cur_block];

                                        if((cur_value > max_value  || (cur_value == max_value && random_functions::nextBool()))
                                           && (cluster_sizes[cur_block] + G.getNodeWeight(node) < block_upperbound || (cur_block == my_block && cluster_sizes[my_block] <= partition_config.upper_bound_partition)))
                                        {
                                                max_value = cur_value;
                                                max_block = cur_block;
                                        }

                                        hash_map[cur_block] = 0;
                                } endfor
                        }

                        cluster_sizes[G.getPartitionIndex(node)]  -= G.getNodeWeight(node);
                        cluster_sizes[max_block]         += G.getNodeWeight(node);
                        bool changed_label                = G.getPartitionIndex(node) != max_block; 
                        change_counter                   += changed_label;
                        G.setPartitionIndex(node, max_block);
                        //std::cout <<  "maxblock " <<  max_block  << std::endl;

                        if(changed_label) {
                                forall_out_edges(G, e, node) {
                                        NodeID target = G.getEdgeTarget(e);
                                        if(!(*next_Q_contained)[target]) {
                                                next_Q->push(target);
                                                (*next_Q_contained)[target] = true;
                                        } 
                                } endfor
                        }
                } 

                std::swap( Q, next_Q);
                std::swap( Q_contained, next_Q_contained);

        }

        end = std::chrono::high_resolution_clock::now();
        std::cout << "Main sequential lp:\t" << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()
                  << std::endl;
        std::cout << "Improved:\t" << change_counter << std::endl;

        delete Q;
        delete next_Q;
        delete Q_contained;
        delete next_Q_contained;
        __itt_pause();


        // in this case the _matching paramter is not used 
        // coarse_mappng stores cluster id and the mapping (it is identical)
        //std::vector<PartitionID> hash_map(G.number_of_nodes(),0);
        //std::vector<NodeID> permutation(G.number_of_nodes());
        //std::vector<NodeWeight> cluster_sizes(partition_config.k,0);

        //forall_nodes(G, node) {
                //cluster_sizes[G.getPartitionIndex(node)] += G.getNodeWeight(node);
        //} endfor
        
        //random_functions::permutate_vector_fast(permutation, true);
        //NodeWeight block_upperbound = partition_config.upper_bound_partition;

        //for( int j = 0; j < partition_config.label_iterations; j++) {
                //forall_nodes(G, i) {
                        //NodeID node = permutation[i];
                        ////move the node to the cluster that is most common in the neighborhood

                        //forall_out_edges(G, e, node) {
                                //NodeID target = G.getEdgeTarget(e);
                                //hash_map[G.getPartitionIndex(target)]+=G.getEdgeWeight(e);
                        //} endfor

                        ////second sweep for finding max and resetting array
                        //PartitionID max_block = G.getPartitionIndex(node);
                        //PartitionID my_block  = G.getPartitionIndex(node);

                        //PartitionID max_value = 0;
                        //forall_out_edges(G, e, node) {
                                //NodeID target             = G.getEdgeTarget(e);
                                //PartitionID cur_block     = G.getPartitionIndex(target);
                                //PartitionID cur_value     = hash_map[cur_block];
                                //if((cur_value > max_value  || (cur_value == max_value && random_functions::nextBool())) 
                                //&& (cluster_sizes[cur_block] + G.getNodeWeight(node) < block_upperbound || cur_block == my_block))
                                //{
                                        //max_value = cur_value;
                                        //max_block = cur_block;
                                //}

                                //hash_map[cur_block] = 0;
                        //} endfor
                        //cluster_sizes[G.getPartitionIndex(node)] -= G.getNodeWeight(node);
                        //cluster_sizes[max_block] += G.getNodeWeight(node);
                        //G.setPartitionIndex(node,max_block);
                //} endfor
        //}
        
        return 0;

}

void label_propagation_refinement::init_for_node_unit(graph_access& G, const size_t block_size, TThreadPool& pool,
                                                      std::vector<Pair>& permutation,
                                                      std::vector<AtomicWrapper<NodeWeight>>& cluster_sizes,
                                                      std::unique_ptr<ConcurrentQueue>& queue) {
        Block block(m_block_allocator);
        block.reserve(block_size);
        forall_nodes(G, node) {
                cluster_sizes[G.getPartitionIndex(node)].fetch_add(G.getNodeWeight(node),
                                                                         std::memory_order_relaxed);
                block.push_back(permutation[node].first);
                if (block.size() == block.capacity()) {
                        // block is full
                        queue->push(std::move(block));
                        block.clear();
                        block.reserve(block_size);
                }
        } endfor

        if (!block.empty()) {
                queue->push(std::move(block));
        }
}

void label_propagation_refinement::seq_init_for_edge_unit(graph_access& G, const size_t block_size, TThreadPool& pool,
                                                      std::vector<Pair>& permutation,
                                                      std::vector<AtomicWrapper<NodeWeight>>& cluster_sizes,
                                                      std::unique_ptr<ConcurrentQueue>& queue) {
        Block block(m_block_allocator);
        block.reserve(block_size);
        size_t cur_block_size = 0;
        forall_nodes(G, node) {
                cluster_sizes[G.getPartitionIndex(node)].fetch_add(G.getNodeWeight(node),
                                                                   std::memory_order_relaxed);
                block.push_back(permutation[node].first);
                cur_block_size += permutation[node].second;
                if (cur_block_size >= block_size) {
                        // block is full
                        queue->push(std::move(block));
                        block.clear();
                        block.reserve(block_size);
                        cur_block_size = 0;
                }
        } endfor
        if (!block.empty()) {
                queue->push(std::move(block));
        }
}

void label_propagation_refinement::par_init_for_edge_unit(graph_access& G, const size_t block_size, TThreadPool& pool,
                                                      std::vector<Pair>& permutation,
                                                      std::vector<AtomicWrapper<NodeWeight>>& cluster_sizes,
                                                      std::unique_ptr<ConcurrentQueue>& queue) {

        forall_nodes(G, node) {
                cluster_sizes[G.getPartitionIndex(node)].fetch_add(G.getNodeWeight(node),
                                                                   std::memory_order_relaxed);
        } endfor

        apply_to_range_sync(NodeID(0), G.number_of_nodes(), pool, [&](auto begin, auto end) {
                Block block(m_block_allocator);
                block.reserve(100);
                size_t cur_block_size = 0;

                for (auto it = begin; it != end; ++it) {
                        NodeID node = *it;
                        block.push_back(permutation[node].first);
                        cur_block_size += permutation[node].second;
                        if (cur_block_size >= block_size) {
                                // block is full
                                queue->push(std::move(block));
                                block.clear();
                                block.reserve(100);
                                cur_block_size = 0;
                        }
                }
                if (!block.empty()) {
                        queue->push(std::move(block));
                }
        });
}

EdgeWeight label_propagation_refinement::parallel_label_propagation_with_queue(graph_access& G,
                                                                               PartitionConfig& config,
                                                                               TThreadPool& pool,
                                                                               Cvector<AtomicWrapper<NodeWeight>>& cluster_sizes1,
                                                                               std::vector<std::vector<PartitionID>>& hash_maps,
                                                                               std::vector<Pair>& permutation) {
        std::cout << "Num threads:\t" << config.num_threads << std::endl;
        __itt_resume();
        const NodeWeight block_upperbound = config.upper_bound_partition;
        auto queue = std::make_unique<ConcurrentQueue>();
        auto next_queue = std::make_unique<ConcurrentQueue>();

        std::vector<AtomicWrapper<bool>> queue_contains(G.number_of_nodes());
        std::vector<AtomicWrapper<bool>> next_queue_contains(G.number_of_nodes());

        size_t block_size = 1;
        const bool use_edge_unit = config.block_size_unit == BlockSizeUnit::EDGES;
        if (config.block_size_unit == BlockSizeUnit::NODES) {
                block_size = config.block_size;
        } else if (config.block_size_unit == BlockSizeUnit::EDGES) {
                //block_size = G.number_of_edges() / (config.num_threads * 16);
                block_size = config.block_size;
        }

        std::cout << "Block size:\t" << block_size << std::endl;

        auto b = std::chrono::high_resolution_clock::now();

        std::vector<AtomicWrapper<NodeWeight>> cluster_sizes(config.k);

        if (use_edge_unit) {
                if (G.number_of_nodes() < 10000000) {
                        seq_init_for_edge_unit(G, block_size, pool, permutation, cluster_sizes, queue);
                } else {
                        par_init_for_edge_unit(G, block_size, pool, permutation, cluster_sizes, queue);
                }
        } else {
                init_for_node_unit(G, block_size, pool, permutation, cluster_sizes, queue);
        }

        auto e = std::chrono::high_resolution_clock::now();
        std::cout << "Init queue lp:\t" << std::chrono::duration_cast<std::chrono::microseconds>(e - b).count()
                  << std::endl;

        std::vector<std::future<NodeWeight>> futures;
        futures.reserve(pool.NumThreads());
        NodeWeight num_changed_label = 0;

        end = std::chrono::high_resolution_clock::now();
        std::cout << "Init parallel lp:\t" << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()
                  << std::endl;


        begin = std::chrono::high_resolution_clock::now();
        //std::atomic<uint32_t> total(0);
        //thread_local std::vector<PartitionID> hash_map(config.k);
        std::cout << "Num blocks\t" << queue->unsafe_size() << std::endl;
        for (int j = 0; j < config.label_iterations_refinement; j++) {
                if (queue->empty()) {
                        break;
                }

                auto process = [&](const size_t id) {
//                        if (hash_map.size() != config.k)
//                                hash_map.resize(config.k);

                        NodeWeight num_changed_label = 0;
                        Block cur_block(m_block_allocator);
                        Block new_block(m_block_allocator);
                        size_t new_block_size = 0;
                        new_block.reserve(block_size);

                        std::uniform_int_distribution<unsigned int> rnd(0,1);
                        std::mt19937 mt;
                        mt.seed(id);

                        EdgeWeight parts_size = config.k * 2;
                        std::vector<PartitionID> partIDs;
                        partIDs.reserve(parts_size);

                        while (queue->try_pop(cur_block)) {
                                for (auto node : cur_block) {
                                        // REMOVE ME <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                                        //total.fetch_add(G.getNodeDegree(node), std::memory_order_relaxed);
                                        // REMOVE ME <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

                                        queue_contains[node].store(false, std::memory_order_relaxed);
                                        auto& hash_map = hash_maps[id];
                                        if (G.getNodeDegree(node) <= parts_size) {
                                                //now move the node to the cluster that is most common in the neighborhood
                                                forall_out_edges(G, e, node) {
                                                        NodeID target = G.getEdgeTarget(e);
                                                        PartitionID part = G.getPartitionIndex(target);
                                                        hash_map[part] += G.getEdgeWeight(e);
                                                        // ONLY FOR UNIT WEIGHTS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                                                        //++hash_map[part];
                                                        // ONLY FOR UNIT WEIGHTS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                                                        partIDs.push_back(part);
                                                } endfor
                                        } else {
                                                //now move the node to the cluster that is most common in the neighborhood
                                                forall_out_edges(G, e, node) {
                                                        NodeID target = G.getEdgeTarget(e);
                                                        hash_map[G.getPartitionIndex(target)] += G.getEdgeWeight(e);

                                                        // ONLY FOR UNIT WEIGHTS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                                                        //++hash_map[G.getPartitionIndex(target)];
                                                        // ONLY FOR UNIT WEIGHTS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                                                } endfor
                                        }

                                        //second sweep for finding max and resetting array
                                        PartitionID my_block = G.getPartitionIndex(node);
                                        PartitionID max_block  = my_block;
                                        NodeWeight max_cluster_size;
                                        PartitionID max_value = 0;

                                        if (G.getNodeDegree(node) <= parts_size) {
                                                for (auto cur_part : partIDs) {
                                                        if (hash_map[cur_part] == 0) {
                                                                continue;
                                                        }

                                                        PartitionID cur_value = hash_map[cur_part];
                                                        NodeWeight cur_cluster_size = cluster_sizes[cur_part].load(
                                                                std::memory_order_acquire);

                                                        if ((cur_value > max_value || (cur_value == max_value
                                                                                       && (bool) rnd(mt)))
                                                            &&
                                                            (cur_cluster_size + G.getNodeWeight(node) < block_upperbound
                                                             || (cur_part == my_block &&
                                                                 cur_cluster_size <= block_upperbound))) {
                                                                max_value = cur_value;
                                                                max_block = cur_part;
                                                                max_cluster_size = cur_cluster_size;
                                                        }

                                                        hash_map[cur_part] = 0;
                                                }
                                                partIDs.clear();
                                        } else {
                                                forall_out_edges(G, e, node) {
                                                        NodeID target = G.getEdgeTarget(e);
                                                        PartitionID cur_part = G.getPartitionIndex(target);
                                                        if (hash_map[cur_part] == 0) {
                                                                continue;
                                                        }

                                                        PartitionID cur_value = hash_map[cur_part];
                                                        NodeWeight cur_cluster_size = cluster_sizes[cur_part].load(
                                                                std::memory_order_acquire);

                                                        if ((cur_value > max_value || (cur_value == max_value
                                                                                       && (bool) rnd(mt)))
                                                            &&
                                                            (cur_cluster_size + G.getNodeWeight(node) < block_upperbound
                                                             || (cur_part == my_block &&
                                                                 cur_cluster_size <= block_upperbound))) {
                                                                max_value = cur_value;
                                                                max_block = cur_part;
                                                                max_cluster_size = cur_cluster_size;
                                                        }

                                                        hash_map[cur_part] = 0;
                                                }
                                                endfor
                                        }


                                        bool changed_label = my_block != max_block;
                                        if (changed_label) {
                                                // try update size of the cluster
                                                bool perform_move = true;
                                                auto& atomic_val = cluster_sizes[max_block];
                                                while (!atomic_val.compare_exchange_weak(
                                                        max_cluster_size,
                                                        max_cluster_size + G.getNodeWeight(node),
                                                        std::memory_order_acq_rel)) {
                                                        if (max_cluster_size + G.getNodeWeight(node) >
                                                            block_upperbound) {
                                                                perform_move = false;
                                                                break;
                                                        }
                                                }

                                                if (perform_move) {
                                                        cluster_sizes[G.getPartitionIndex(node)].
                                                                fetch_add(-G.getNodeWeight(node),
                                                                          std::memory_order_acq_rel);

                                                        G.setPartitionIndex(node, max_block);

                                                        ++num_changed_label;

                                                        forall_out_edges(G, e, node)
                                                        {
                                                                NodeID target = G.getEdgeTarget(e);
                                                                if (!next_queue_contains[target].exchange(
                                                                        true, std::memory_order_acq_rel)) {
                                                                        new_block.push_back(target);

                                                                        new_block_size += use_edge_unit ? G.getNodeDegree(target) : 1;
                                                                        if (new_block_size >= block_size) {
                                                                                next_queue->push(std::move(new_block));
                                                                                new_block.clear();
                                                                                new_block.reserve(block_size);
                                                                                new_block_size = 0;
                                                                        }
                                                                }
                                                        }
                                                        endfor
                                                }
                                        }
                                }
                        }
                        if (!new_block.empty()) {
                                next_queue->push(new_block);
                                new_block_size = 0;
                        }
                        return num_changed_label;
                };


                for (size_t i = 0; i < pool.NumThreads(); ++i) {
                        futures.push_back(pool.Submit(process, i + 1));
                }

                num_changed_label += process(0);
                std::for_each(futures.begin(), futures.end(), [&](auto& future){
                        num_changed_label += future.get();
                });
                //std::cout << "Queue size\t" << total << std::endl;
                //total = 0;
                std::swap(queue, next_queue);
                std::swap(queue_contains, next_queue_contains);
                futures.clear();
        }
        __itt_pause();
        end = std::chrono::high_resolution_clock::now();
        std::cout << "Main parallel lp:\t" << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()
                  << std::endl;
        std::cout << "Improved:\t" << num_changed_label << std::endl;
        return num_changed_label;
}

EdgeWeight label_propagation_refinement::parallel_label_propagation(graph_access& G,
                                                                    PartitionConfig& config,
                                                                    TThreadPool& pool,
                                                                    Cvector<AtomicWrapper<NodeWeight>>& cluster_sizes,
                                                                    std::vector<std::vector<PartitionID>>& hash_maps,
                                                                    std::vector<Pair>& permutation) {
        NodeWeight block_upperbound = config.upper_bound_partition;
        std::vector<std::future<NodeWeight>> futures;
        futures.reserve(pool.NumThreads());
        NodeWeight num_changed_label = 0;

        end = std::chrono::high_resolution_clock::now();
        std::cout << "Init parallel lp:\t" << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()
                  << std::endl;

        begin = std::chrono::high_resolution_clock::now();
        for (int j = 0; j < config.label_iterations_refinement; j++) {
                auto process = [&](const size_t id, NodeID begin, NodeID end) {
                        auto& hash_map = hash_maps[id];
                        NodeWeight num_changed_label = 0;
                        std::uniform_int_distribution<unsigned int> rnd(0,1);
                        std::mt19937 mt;
                        mt.seed(id);
                        for (NodeID node = begin; node != end; ++node) {
                                //now move the node to the cluster that is most common in the neighborhood
                                forall_out_edges(G, e, node) {
                                        NodeID target = G.getEdgeTarget(e);
                                        hash_map[G.getPartitionIndex(target)] += G.getEdgeWeight(e);
                                } endfor

                                //second sweep for finding max and resetting array
                                PartitionID max_block = G.getPartitionIndex(node);
                                PartitionID my_block  = G.getPartitionIndex(node);
                                NodeWeight max_cluster_size = cluster_sizes[max_block].get().load(std::memory_order_acquire);

                                PartitionID max_value = 0;
                                forall_out_edges(G, e, node) {
                                        NodeID target             = G.getEdgeTarget(e);
                                        PartitionID cur_block     = G.getPartitionIndex(target);
                                        PartitionID cur_value     = hash_map[cur_block];
                                        NodeWeight cur_cluster_size = cluster_sizes[cur_block].get().load(std::memory_order_acquire);

                                        if((cur_value > max_value  || (cur_value == max_value
                                                                       && (bool) rnd(mt)))
                                           && (cur_cluster_size + G.getNodeWeight(node) < block_upperbound
                                               || (cur_block == my_block && cur_cluster_size <= block_upperbound)))
                                        {
                                                max_value = cur_value;
                                                max_block = cur_block;
                                                max_cluster_size = cur_cluster_size;
                                        }

                                        hash_map[cur_block] = 0;
                                } endfor

                                bool changed_label = my_block != max_block;
                                if (changed_label) {
                                        // try update size of the cluster
                                        bool perform_move = true;
                                        auto& atomic_val = cluster_sizes[max_block].get();
                                        while (!atomic_val.compare_exchange_weak(
                                                max_cluster_size,
                                                max_cluster_size + G.getNodeWeight(node),
                                                std::memory_order_acq_rel)) {
                                                if (max_cluster_size + G.getNodeWeight(node) >
                                                    block_upperbound) {
                                                        perform_move = false;
                                                        break;
                                                }
                                        }

                                        if (perform_move) {
                                                cluster_sizes[G.getPartitionIndex(node)].get().fetch_add(
                                                        -G.getNodeWeight(node), std::memory_order_acq_rel);
                                                G.setPartitionIndex(node, max_block);
                                                ++num_changed_label;
                                        }
                                }
                        }
                        return num_changed_label;
                };

                size_t work_per_thread = G.number_of_nodes() / (pool.NumThreads() + 1);
                NodeID first = 0;
                for (size_t i = 0; i < pool.NumThreads(); ++i) {
                        futures.push_back(pool.Submit(process, i + 1, first, first + work_per_thread));
                        first += work_per_thread;
                }

                num_changed_label += process(0, first, G.number_of_nodes());
                std::for_each(futures.begin(), futures.end(), [&](auto& future){
                        num_changed_label += future.get();
                });
                futures.clear();
        }
        end = std::chrono::high_resolution_clock::now();
        std::cout << "Main parallel (no queue) lp:\t" << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()
                  << std::endl;
        std::cout << "Improved:\t" << num_changed_label << std::endl;
        return num_changed_label;
}

EdgeWeight label_propagation_refinement::parallel_label_propagation(PartitionConfig& config, graph_access& G,
                                                                    complete_boundary& boundary) {
        // TODO: pool should be passed to all functions
        begin = std::chrono::high_resolution_clock::now();
        auto b = std::chrono::high_resolution_clock::now();
        std::vector <std::vector<PartitionID>> hash_maps(config.num_threads, std::vector<PartitionID>(config.k));
        Cvector <AtomicWrapper<NodeWeight>> cluster_sizes(config.k);
        auto e = std::chrono::high_resolution_clock::now();
        std::cout << "Init other vectors lp:\t" << std::chrono::duration_cast<std::chrono::microseconds>(e - b).count()
                  << std::endl;

        b = std::chrono::high_resolution_clock::now();
        std::vector<Pair> permutation;
        permutation.reserve(G.number_of_nodes());
//        node_ordering n_ordering;
//        n_ordering.order_nodes(config, G, permutation);
        forall_nodes(G, node) {
                permutation.emplace_back(node, G.getNodeDegree(node));
        } endfor
//        tbb::global_control c(tbb::global_control::max_allowed_parallelism, config.num_threads + 1);
//        tbb::parallel_sort(permutation.begin(), permutation.end(), [&](Pair lhs, Pair rhs) {
//                return G.getNodeDegree(lhs) < G.getNodeDegree(rhs);
//        });
        omp_set_dynamic(false);
        omp_set_num_threads(config.num_threads);
        __gnu_parallel::sort(permutation.begin(), permutation.end(), [&](const Pair& lhs, const Pair& rhs) {
                return lhs.second < rhs.second;
        });
        e = std::chrono::high_resolution_clock::now();
        std::cout << "Parallel init of permutations lp:\t"
                  << std::chrono::duration_cast<std::chrono::microseconds>(e - b).count()
                  << std::endl;

        TThreadPool pool(config.num_threads - 1);
        if (config.parallel_lp_type == ParallelLPType::NO_QUEUE) {
                return parallel_label_propagation(G, config, pool, cluster_sizes, hash_maps, permutation);
        } else if (config.parallel_lp_type == ParallelLPType::QUEUE) {
                return parallel_label_propagation_with_queue(G, config, pool, cluster_sizes, hash_maps, permutation);
        } else {
                return 0;
        }
}

