/******************************************************************************
 * contraction.cpp 
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

#include "contraction.h"
#include "data_structure/parallel/hash_table.h"
#include "data_structure/parallel/time.h"
#include "../uncoarsening/refinement/quotient_graph_refinement/complete_boundary.h"
#include "macros_assertions.h"

#include "ittnotify.h"

contraction::contraction() {

}

contraction::~contraction() {

}

// for documentation see technical reports of christian schulz  
void contraction::contract(const PartitionConfig & partition_config, 
                           graph_access & G, 
                           graph_access & coarser, 
                           const Matching & edge_matching,
                           const CoarseMapping & coarse_mapping,
                           const NodeID & no_of_coarse_vertices,
                           const NodePermutationMap & permutation) const {

        if(partition_config.matching_type == CLUSTER_COARSENING) {
                if (!partition_config.fast_contract_clustering) {
                        return contract_clustering(partition_config, G, coarser, edge_matching, coarse_mapping, no_of_coarse_vertices, permutation);
                } else {
                        return fast_contract_clustering(partition_config, G, coarser, edge_matching, coarse_mapping, no_of_coarse_vertices, permutation);
                }
        }

        if(partition_config.combine) {
                coarser.resizeSecondPartitionIndex(no_of_coarse_vertices);
        }

        std::vector<NodeID> new_edge_targets(G.number_of_edges());
        forall_edges(G, e) {
                new_edge_targets[e] = coarse_mapping[G.getEdgeTarget(e)];
        } endfor

        std::vector<EdgeID> edge_positions(no_of_coarse_vertices, UNDEFINED_EDGE);

        //we dont know the number of edges jet, so we use the old number for 
        //construction of the coarser graph and then resize the field according
        //to the number of edges we really got
        coarser.start_construction(no_of_coarse_vertices, G.number_of_edges());

        NodeID cur_no_vertices = 0;

        forall_nodes(G, n) {
                NodeID node = permutation[n];
                //we look only at the coarser nodes
                if(coarse_mapping[node] != cur_no_vertices) 
                        continue;
                
                NodeID coarseNode = coarser.new_node();
                coarser.setNodeWeight(coarseNode, G.getNodeWeight(node));

                if(partition_config.combine) {
                        coarser.setSecondPartitionIndex(coarseNode, G.getSecondPartitionIndex(node));
                }

                // do something with all outgoing edges (in auxillary graph)
                forall_out_edges(G, e, node) {
                        visit_edge(G, coarser, edge_positions, coarseNode, e, new_edge_targets);                        
                } endfor

                //this node was really matched
                NodeID matched_neighbor = edge_matching[node];
                if(node != matched_neighbor) {
                        //update weight of coarser node
                        NodeWeight new_coarse_weight = G.getNodeWeight(node) + G.getNodeWeight(matched_neighbor);
                        coarser.setNodeWeight(coarseNode, new_coarse_weight);

                        forall_out_edges(G, e, matched_neighbor) {
                                visit_edge(G, coarser, edge_positions, coarseNode, e, new_edge_targets);
                        } endfor
                }
                forall_out_edges(coarser, e, coarseNode) {
                       edge_positions[coarser.getEdgeTarget(e)] = UNDEFINED_EDGE;
                } endfor
                
                cur_no_vertices++;
        } endfor

        ASSERT_RANGE_EQ(edge_positions, 0, edge_positions.size(), UNDEFINED_EDGE); 
        ASSERT_EQ(no_of_coarse_vertices, cur_no_vertices);
        
        //this also resizes the edge fields ... 
        coarser.finish_construction();
}

void contraction::contract_clustering(const PartitionConfig & partition_config, 
                              graph_access & G, 
                              graph_access & coarser, 
                              const Matching & edge_matching,
                              const CoarseMapping & coarse_mapping,
                              const NodeID & no_of_coarse_vertices,
                              const NodePermutationMap & permutation) const {

        if(partition_config.combine) {
                coarser.resizeSecondPartitionIndex(no_of_coarse_vertices);
        }

        //save partition map -- important if the graph is allready partitioned
        std::vector< int > partition_map(G.number_of_nodes());
        int k = G.get_partition_count();
        forall_nodes(G, node) {
                partition_map[node] = G.getPartitionIndex(node);
                G.setPartitionIndex(node, coarse_mapping[node]);
        } endfor

        G.set_partition_count(no_of_coarse_vertices);

        complete_boundary bnd(&G);
        bnd.build();
        bnd.getUnderlyingQuotientGraph(coarser);

        G.set_partition_count(k);
        forall_nodes(G, node) {
                G.setPartitionIndex(node, partition_map[node]);
                coarser.setPartitionIndex(coarse_mapping[node], G.getPartitionIndex(node));

                if(partition_config.combine) {
                        coarser.setSecondPartitionIndex(coarse_mapping[node], G.getSecondPartitionIndex(node));
                }

        } endfor

}

void contraction::fast_contract_clustering(const PartitionConfig& partition_config,
                                      graph_access& G,
                                      graph_access& coarser,
                                      const Matching&,
                                      const CoarseMapping& coarse_mapping,
                                      const NodeID& no_of_coarse_vertices,
                                      const NodePermutationMap&) const {
        __itt_resume();
        if (partition_config.combine) {
                coarser.resizeSecondPartitionIndex(no_of_coarse_vertices);
        }

        block_informations empty;
        std::vector<block_informations> block_infos(no_of_coarse_vertices, empty);

        // build set of new edges
        auto get_pair = [](NodeID cluster_a, NodeID cluster_b) {
                if (cluster_a <= cluster_b) {
                        return std::make_pair(cluster_a, cluster_b);
                } else {
                        return std::make_pair(cluster_b, cluster_a);
                }
        };

        double avg_degree = (G.number_of_edges() + 0.0) / G.number_of_nodes();
        size_t num_cut_edges = std::min<size_t>(avg_degree * no_of_coarse_vertices, G.number_of_edges() / 2);
        using key_type = std::pair<NodeID, NodeID>;
        parallel::HashMap<key_type, EdgeWeight, parallel::MurmurHash<key_type>, true> new_edges(num_cut_edges);

        forall_nodes(G, n) {
                PartitionID source_cluster = coarse_mapping[n];
                block_infos[source_cluster].block_weight += G.getNodeWeight(n);
                block_infos[source_cluster].block_no_nodes += 1;


                forall_out_edges(G, e, n) {
                        NodeID targetID = G.getEdgeTarget(e);
                        PartitionID target_cluster = coarse_mapping[targetID];
                        bool is_cut_edge = source_cluster != target_cluster;

                        if (is_cut_edge) {
                                new_edges[get_pair(source_cluster, target_cluster)] += G.getEdgeWeight(e);
                        }
                } endfor
        } endfor

        // construct graph
        std::vector<std::vector<std::pair<PartitionID, EdgeWeight>>> building_tool(no_of_coarse_vertices);
        for (auto& data : building_tool) {
                data.reserve(avg_degree);
        }

        for (const auto& edge_data : new_edges) {
                auto& edge = edge_data.first;
                auto edge_weight = edge_data.second / 2;
                building_tool[edge.first].emplace_back(edge.second, edge_weight);
                building_tool[edge.second].emplace_back(edge.first, edge_weight);
        }

        coarser.start_construction(building_tool.size(), 2 * new_edges.size());

        for (size_t p = 0; p < building_tool.size(); p++) {
                NodeID node = coarser.new_node();
                coarser.setNodeWeight(node, block_infos[p].block_weight);

                for(size_t j = 0; j < building_tool[p].size(); j++) {
                        EdgeID e = coarser.new_edge(node, building_tool[p][j].first);
                        coarser.setEdgeWeight(e, building_tool[p][j].second);
                }
        }
        coarser.finish_construction();

        forall_nodes(G, node) {
                coarser.setPartitionIndex(coarse_mapping[node], G.getPartitionIndex(node));

                if(partition_config.combine) {
                        coarser.setSecondPartitionIndex(coarse_mapping[node], G.getSecondPartitionIndex(node));
                }
        } endfor
        __itt_pause();
}

// for documentation see technical reports of christian schulz  
void contraction::contract_partitioned(const PartitionConfig & partition_config, 
                                       graph_access & G, 
                                       graph_access & coarser, 
                                       const Matching & edge_matching,
                                       const CoarseMapping & coarse_mapping,
                                       const NodeID & no_of_coarse_vertices,
                                       const NodePermutationMap & permutation) const {

        if(partition_config.matching_type == CLUSTER_COARSENING) {
                return contract_clustering(partition_config, G, coarser, edge_matching, coarse_mapping, no_of_coarse_vertices, permutation);
        }


        std::vector<NodeID> new_edge_targets(G.number_of_edges());
        forall_edges(G, e) {
                new_edge_targets[e] = coarse_mapping[G.getEdgeTarget(e)];
        } endfor

        std::vector<EdgeID> edge_positions(no_of_coarse_vertices, UNDEFINED_EDGE);

        //we dont know the number of edges jet, so we use the old number for 
        //construction of the coarser graph and then resize the field according
        //to the number of edges we really got
        coarser.set_partition_count(G.get_partition_count());
        coarser.start_construction(no_of_coarse_vertices, G.number_of_edges());

        if(partition_config.combine) {
                coarser.resizeSecondPartitionIndex(no_of_coarse_vertices);
        }

        NodeID cur_no_vertices = 0;

        PRINT(std::cout <<  "contracting a partitioned graph"  << std::endl;)
        forall_nodes(G, n) {
                NodeID node = permutation[n];
                //we look only at the coarser nodes
                if(coarse_mapping[node] != cur_no_vertices) 
                        continue;
                
                NodeID coarseNode = coarser.new_node();
                coarser.setNodeWeight(coarseNode, G.getNodeWeight(node));
                coarser.setPartitionIndex(coarseNode, G.getPartitionIndex(node));

                if(partition_config.combine) {
                        coarser.setSecondPartitionIndex(coarseNode, G.getSecondPartitionIndex(node));
                }
                // do something with all outgoing edges (in auxillary graph)
                forall_out_edges(G, e, node) {
                                visit_edge(G, coarser, edge_positions, coarseNode, e, new_edge_targets);                        
                } endfor

                //this node was really matched
                NodeID matched_neighbor = edge_matching[node];
                if(node != matched_neighbor) {
                        //update weight of coarser node
                        NodeWeight new_coarse_weight = G.getNodeWeight(node) + G.getNodeWeight(matched_neighbor);
                        coarser.setNodeWeight(coarseNode, new_coarse_weight);

                        forall_out_edges(G, e, matched_neighbor) {
                                visit_edge(G, coarser, edge_positions, coarseNode, e, new_edge_targets);
                        } endfor
                }
                forall_out_edges(coarser, e, coarseNode) {
                       edge_positions[coarser.getEdgeTarget(e)] = UNDEFINED_EDGE;
                } endfor
                
                cur_no_vertices++;
        } endfor

        ASSERT_RANGE_EQ(edge_positions, 0, edge_positions.size(), UNDEFINED_EDGE); 
        ASSERT_EQ(no_of_coarse_vertices, cur_no_vertices);
        
        //this also resizes the edge fields ... 
        coarser.finish_construction();
}



