#pragma once

#include "data_structure/parallel/bits.h"
#include "data_structure/graph_access.h"

#include <algorithm>
#include <vector>
#include <unordered_map>

namespace parallel {

static void shuffle_graph(graph_access& graph, graph_access& shuffled_graph) {
        std::vector<NodeID> nodes_perm;
        nodes_perm.reserve(graph.number_of_nodes());

        for (size_t i = 0; i < graph.number_of_nodes(); ++i) {
                nodes_perm.push_back(i);
        }

        std::random_shuffle(nodes_perm.begin(), nodes_perm.end());

        std::vector<std::vector<std::pair<NodeID, EdgeWeight>>> new_edges(graph.number_of_nodes());
        std::vector<NodeWeight> new_weights(graph.number_of_nodes());

        for (NodeID node = 0; node < nodes_perm.size(); ++node) {
                NodeID node_perm = nodes_perm[node];
                new_edges[node_perm].reserve(graph.getNodeDegree(node));
                new_weights[node_perm] = graph.getNodeWeight(node);
                forall_out_edges(graph, e, node){
                                        NodeID target = graph.getEdgeTarget(e);
                                        NodeID target_perm = nodes_perm[target];
                                        EdgeWeight weight = graph.getEdgeWeight(e);

                                        new_edges[node_perm].emplace_back(target_perm, weight);
                                }endfor
        }

        shuffled_graph.start_construction(graph.number_of_nodes(), graph.number_of_edges());

        for (size_t i = 0; i < new_edges.size(); ++i) {
                NodeID node = shuffled_graph.new_node();
                shuffled_graph.setNodeWeight(node, new_weights[i]);

                for (size_t j = 0; j < new_edges[i].size(); ++j) {
                        EdgeID e = shuffled_graph.new_edge(node, new_edges[i][j].first);
                        shuffled_graph.setEdgeWeight(e, new_edges[i][j].second);
                }
        }

        shuffled_graph.finish_construction();
}

static void update_hist(const std::vector<uint32_t>& input, std::unordered_map<uint32_t, uint32_t>& hist,
                        uint32_t count) {
        for (size_t i = 0; i + 1 < input.size(); ++i) {
                uint32_t res = parallel::least_significant_bit(input[i] ^ input[i + 1]);
                hist[res] += count;
        }
}

static std::unordered_map<uint32_t, uint32_t> get_bit_diff_hist(graph_access& graph,
                                                                const std::vector<uint32_t>* evaluate) {
        std::vector<NodeID> neighbors;
        std::unordered_map<uint32_t, uint32_t> hist;
        for (NodeID node = 0; node < graph.number_of_nodes(); ++node) {
                uint32_t count = 1;
                if (evaluate != nullptr) {
                        count = !(*evaluate)[node];

                        if (count == 0) {
                                continue;
                        }
                }

                neighbors.clear();
                neighbors.reserve(graph.getNodeDegree(node));

                forall_out_edges(graph, e, node){
                                        NodeID target = graph.getEdgeTarget(e);
                                        neighbors.push_back(target);
                                }endfor

                update_hist(neighbors, hist, count);
        }

        return hist;
};

static std::unordered_map<uint32_t, uint32_t> get_bit_diff_hist(graph_access& graph) {
        return get_bit_diff_hist(graph, nullptr);
}

static std::unordered_map<uint32_t, uint32_t> get_bit_diff_hist(graph_access& graph,
                                                                const std::vector<uint32_t>& evaluate) {
        return get_bit_diff_hist(graph, &evaluate);
}

static void sort_edges(graph_access& graph, graph_access& sorted_graph) {
        std::vector<std::vector<std::pair<NodeID, EdgeWeight>>> new_edges(graph.number_of_nodes());

        for (NodeID node = 0; node < graph.number_of_nodes(); ++node) {
                new_edges[node].reserve(graph.getNodeDegree(node));
                forall_out_edges(graph, e, node){
                                        NodeID target = graph.getEdgeTarget(e);
                                        EdgeWeight weight = graph.getEdgeWeight(e);

                                        new_edges[node].emplace_back(target, weight);
                                }endfor
                std::sort(new_edges[node].begin(), new_edges[node].end());
        }

        sorted_graph.start_construction(graph.number_of_nodes(), graph.number_of_edges());

        for (size_t i = 0; i < new_edges.size(); ++i) {
                NodeID node = sorted_graph.new_node();
                sorted_graph.setNodeWeight(node, graph.getNodeWeight(node));

                for (size_t j = 0; j < new_edges[i].size(); ++j) {
                        EdgeID e = sorted_graph.new_edge(node, new_edges[i][j].first);
                        sorted_graph.setEdgeWeight(e, new_edges[i][j].second);
                }
        }

        sorted_graph.finish_construction();
}

}
