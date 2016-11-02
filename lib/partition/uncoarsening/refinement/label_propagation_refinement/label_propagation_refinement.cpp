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

#include <tbb/concurrent_queue.h>

#define TBB_PREVIEW_GLOBAL_CONTROL 1
#include <tbb/global_control.h>

#include <tbb/parallel_sort.h>

#include <chrono>
#include <random>

#include "ittnotify.h"

using namespace parallel;

label_propagation_refinement::label_propagation_refinement() {
                
}

label_propagation_refinement::~label_propagation_refinement() {

}

EdgeWeight label_propagation_refinement::perform_refinement(PartitionConfig& config, graph_access& G,
                                                            complete_boundary& boundary) {

        if (!config.parallel_local_search || G.number_of_nodes() < 100000) {
                auto begin = std::chrono::high_resolution_clock::now();
                auto res = sequential_label_propagation(config, G, boundary);
                auto end = std::chrono::high_resolution_clock::now();

                std::cout << "Sequential lp:\t" << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()
                          << std::endl;
                return res;
        } else {
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
        unsigned int change_counter = 0;
        for( int j = 0; j < partition_config.label_iterations_refinement; j++) {
                while( !Q->empty() ) {
                        NodeID node = Q->front();
                        Q->pop();
                        (*Q_contained)[node] = false;

                        //now move the node to the cluster that is most common in the neighborhood
                        forall_out_edges(G, e, node) {
                                NodeID target = G.getEdgeTarget(e);
                                hash_map[G.getPartitionIndex(target)]+=G.getEdgeWeight(e);
                                //std::cout <<  "curblock " <<  G.getPartitionIndex(target)  << std::endl;
                        } endfor

                        //second sweep for finding max and resetting array
                        PartitionID max_block = G.getPartitionIndex(node);
                        PartitionID my_block  = G.getPartitionIndex(node);
                        PartitionID max_value = 0;

                        forall_out_edges(G, e, node) {
                                NodeID target             = G.getEdgeTarget(e);
                                PartitionID cur_block     = G.getPartitionIndex(target);
                                PartitionID cur_value     = hash_map[cur_block];

                                if((cur_value > max_value  || (cur_value == max_value && random_functions::nextBool()))
                                && (cluster_sizes[cur_block] + G.getNodeWeight(node) < block_upperbound || (cur_block == my_block && cluster_sizes[my_block] <= partition_config.upper_bound_partition)))
                                {
                                        max_value = cur_value;
                                        max_block = cur_block;
                                }

                                hash_map[cur_block] = 0;
                        } endfor

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


EdgeWeight label_propagation_refinement::parallel_label_propagation_with_queue(graph_access& G,
                                                                               PartitionConfig& config,
                                                                               TThreadPool& pool,
                                                                               Cvector<AtomicWrapper<NodeWeight>>& cluster_sizes1,
                                                                               std::vector<std::vector<PartitionID>>& hash_maps,
                                                                               std::vector<NodeID>& permutation) {
        using Block = std::vector<NodeID>;
        __itt_resume();
        std::cout << "Num threads:\t" << config.num_threads << std::endl;

        const NodeWeight block_upperbound = config.upper_bound_partition;
        auto queue = std::make_unique<tbb::concurrent_queue<Block>>();
        auto next_queue = std::make_unique<tbb::concurrent_queue<Block>>();

        std::vector<AtomicWrapper<bool>> queue_contains(G.number_of_nodes());
        std::vector<AtomicWrapper<bool>> next_queue_contains(G.number_of_nodes());

        const size_t block_size = config.block_size;
        std::cout << "Block size:\t" << block_size << std::endl;

        auto b = std::chrono::high_resolution_clock::now();
        Block block;
        block.reserve(block_size);

        std::vector<AtomicWrapper<NodeWeight>> cluster_sizes(config.k);

        forall_nodes(G, node) {
                cluster_sizes[G.getPartitionIndex(node)].fetch_add(G.getNodeWeight(node),
                                                                         std::memory_order_relaxed);
                block.push_back(permutation[node]);
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
        //thread_local std::vector<PartitionID> hash_map(config.k);
        for (int j = 0; j < config.label_iterations_refinement; j++) {
                if (queue->empty()) {
                        break;
                }

                auto process = [&](const size_t id) {
//                        if (hash_map.size() != config.k)
//                                hash_map.resize(config.k);

                        NodeWeight num_changed_label = 0;
                        Block cur_block;
                        Block new_block;
                        new_block.reserve(block_size);

                        std::uniform_int_distribution<unsigned int> rnd(0,1);
                        std::mt19937 mt;
                        mt.seed(id);
                        while (queue->try_pop(cur_block)) {
                                for (auto node : cur_block) {
                                        queue_contains[node].store(false, std::memory_order_relaxed);
                                        auto& hash_map = hash_maps[id];
                                        //now move the node to the cluster that is most common in the neighborhood
                                        forall_out_edges(G, e, node) {
                                                NodeID target = G.getEdgeTarget(e);
                                                hash_map[G.getPartitionIndex(target)] += G.getEdgeWeight(e);
                                        } endfor

                                        //second sweep for finding max and resetting array
                                        PartitionID max_block = G.getPartitionIndex(node);
                                        PartitionID my_block  = G.getPartitionIndex(node);
                                        NodeWeight max_cluster_size;

                                        PartitionID max_value = 0;

                                        forall_out_edges(G, e, node) {
                                                NodeID target             = G.getEdgeTarget(e);
                                                PartitionID cur_block     = G.getPartitionIndex(target);
                                                PartitionID cur_value     = hash_map[cur_block];
                                                NodeWeight cur_cluster_size = cluster_sizes[cur_block].load(std::memory_order_acquire);

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
                                                                        if (new_block.size() == new_block.capacity()) {
                                                                                next_queue->push(std::move(new_block));
                                                                                new_block.clear();
                                                                                new_block.reserve(block_size);
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
                std::swap(queue, next_queue);
                std::swap(queue_contains, next_queue_contains);
                futures.clear();
        }
        end = std::chrono::high_resolution_clock::now();
        __itt_pause();
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
                                                                    std::vector<NodeID>& permutation) {
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
        TThreadPool pool(config.num_threads - 1);

        begin = std::chrono::high_resolution_clock::now();
        auto b = std::chrono::high_resolution_clock::now();
        std::vector<std::vector<PartitionID>> hash_maps(config.num_threads, std::vector<PartitionID>(config.k));
        Cvector<AtomicWrapper<NodeWeight>> cluster_sizes(config.k);
        auto e = std::chrono::high_resolution_clock::now();
        std::cout << "Init other vectors lp:\t" << std::chrono::duration_cast<std::chrono::microseconds>(e - b).count()
                  << std::endl;

        b = std::chrono::high_resolution_clock::now();
        std::vector<NodeID> permutation(G.number_of_nodes());
//        node_ordering n_ordering;
//        n_ordering.order_nodes(config, G, permutation);
        forall_nodes(G, node) {
                permutation[node] = node;
        } endfor
        tbb::global_control c(tbb::global_control::max_allowed_parallelism, config.num_threads + 1);
        tbb::parallel_sort(permutation.begin(), permutation.end(), [&](NodeID lhs, NodeID rhs) {
                return G.getNodeDegree(lhs) < G.getNodeDegree(rhs);
        });
        e = std::chrono::high_resolution_clock::now();
        std::cout << "Parallel init of permutations lp:\t" << std::chrono::duration_cast<std::chrono::microseconds>(e - b).count()
                  << std::endl;

        //const parallel_lp_type lp_type = parallel_lp_type::parallel;
        const parallel_lp_type lp_type = parallel_lp_type::parallel_queue;

        if (lp_type == parallel_lp_type::parallel) {
                return parallel_label_propagation(G, config, pool, cluster_sizes, hash_maps, permutation);
        } else if (lp_type == parallel_lp_type::parallel_queue) {
                return parallel_label_propagation_with_queue(G, config, pool, cluster_sizes, hash_maps, permutation);
        } else {
                throw std::logic_error("Unknown parallel lp type");
        }
}

