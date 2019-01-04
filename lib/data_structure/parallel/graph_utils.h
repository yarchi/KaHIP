#pragma once

#include "data_structure/graph_access.h"

#include <algorithm>
#include <vector>
#include <fstream>

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
                forall_out_edges(graph, e, node) {
                        NodeID target = graph.getEdgeTarget(e);
                        NodeID target_perm = nodes_perm[target];
                        EdgeWeight weight = graph.getEdgeWeight(e);

                        new_edges[node_perm].emplace_back(target_perm, weight);
                } endfor
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

template <typename T>
static void get_deltas(const std::vector<T>& input, std::vector<T>& deltas) {
//        for (size_t i = 0; i + 1 < input.size(); ++i) {
//                deltas.push_back(input[i + 1] - input[i]);
//        }

        size_t res = 0;
        // 8 elements of size 8 bytes
        size_t cache_line_size = 8;
        size_t i = 0;
        while (i < input.size()) {
                T start = input[i];
                T cur = start;
                while (cur < start + cache_line_size) {
                        ++res;
                        ++i;
                        cur = input[i];
                }
        }
}

double get_average_chain_length(graph_access& graph, NodeID delta) {
        std::vector<NodeID> neighbors;

        NodeID total_chain_length = 0;
        NodeID num_chains = 0;

        for (NodeID node = 0; node < graph.number_of_nodes(); ++node) {
                neighbors.clear();

                neighbors.reserve(graph.getNodeDegree(node));

                forall_out_edges(graph, e, node) {
                        NodeID target = graph.getEdgeTarget(e);
                        neighbors.push_back(target);
                } endfor

                if (neighbors.empty()) {
                        continue;
                }
                std::sort(neighbors.begin(), neighbors.end());
                NodeID chain_begin = 0;
                for (size_t i = 0; i + 1 < neighbors.size(); ++i) {
                        if (neighbors[i + 1] - neighbors[i] > delta) {
                                ++num_chains;
                                total_chain_length += neighbors[i] - neighbors[chain_begin];
                                chain_begin = i + 1;
                        }
                }
        }

        return (total_chain_length + 0.0) / num_chains;
}

static std::pair<double, double> average_chain_length(graph_access& graph) {
        std::vector<NodeID> neighbors;
        std::vector<NodeID> deltas;
        std::vector<NodeID> medians;
        medians.reserve(graph.number_of_nodes());
        size_t total_median = 0;
        for (NodeID node = 0; node < graph.number_of_nodes(); ++node) {
                neighbors.clear();
                deltas.clear();

                neighbors.reserve(graph.getNodeDegree(node));
                deltas.reserve(graph.getNodeDegree(node));

                forall_out_edges(graph, e, node) {
                        NodeID target = graph.getEdgeTarget(e);
                        neighbors.push_back(target);
                } endfor

                std::sort(neighbors.begin(), neighbors.end());
                get_deltas(neighbors, deltas);
                if (!deltas.empty()) {
                        std::nth_element(deltas.begin(), deltas.begin() + deltas.size() / 2, deltas.end());
                        NodeID median = deltas[deltas.size() / 2];
                        medians.push_back(median);
                        total_median += median;
                }
        }

        std::nth_element(medians.begin(), medians.begin() + medians.size() / 2, medians.end());


        auto res1 = get_average_chain_length(graph, (total_median + 0.0) / graph.number_of_nodes());
        auto res2 = get_average_chain_length(graph, medians[medians.size() / 2]);

        return {res1, res2};
}

static void sort_edges(graph_access& graph, graph_access& sorted_graph) {
        std::vector<std::vector<std::pair<NodeID, EdgeWeight>>> new_edges(graph.number_of_nodes());

        for (NodeID node = 0; node < graph.number_of_nodes(); ++node) {
                new_edges[node].reserve(graph.getNodeDegree(node));
                forall_out_edges(graph, e, node) {
                        NodeID target = graph.getEdgeTarget(e);
                        EdgeWeight weight = graph.getEdgeWeight(e);

                        new_edges[node].emplace_back(target, weight);
                } endfor
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

static void print_graph_stat(graph_access& graph) {
        uint64_t edge_weight = 0;
        uint64_t node_weight = 0;
        double avg_degree = (graph.number_of_edges() + 0.0) / graph.number_of_nodes();

        std::vector<std::pair<uint64_t, uint32_t>> topEdgeDegree;
        topEdgeDegree.reserve(graph.number_of_nodes());
        for (NodeID node = 0; node < graph.number_of_nodes(); ++node) {
                node_weight += graph.getNodeWeight(node);
                topEdgeDegree.emplace_back(graph.getNodeDegree(node), node);
                forall_out_edges(graph, e, node) {
                        NodeID target = graph.getEdgeTarget(e);
                        EdgeWeight weight = graph.getEdgeWeight(e);

                        edge_weight += weight;
                } endfor

        }



        std::sort(topEdgeDegree.begin(), topEdgeDegree.end(), std::greater<std::pair<uint64_t, uint32_t>>());
        uint64_t cur_deg = 0;

        size_t i = 0;
        double quantile = 0.05;
        uint64_t total_deg = 2 * graph.number_of_edges();

        std::cout << "Total degree\t" << total_deg << std::endl;

        std::vector<uint32_t> to_remove;
        while (cur_deg <= total_deg * 0.5 && i < topEdgeDegree.size())
        {
                cur_deg += topEdgeDegree[i].first;
                std::cout << "Vertex " << topEdgeDegree[i].second << " has degree " << topEdgeDegree[i].first << std::endl;
                to_remove.push_back(topEdgeDegree[i].second);
                if (cur_deg >= total_deg * quantile) {
                        std::cout << "Top " << i + 1 << " vertices have total degree of " <<  cur_deg
                                  << ". This is " << (cur_deg + 0.0) / total_deg << "of all total degree" << std::endl;
                        quantile += 0.05;
                        break;
                }
                i++;
        }

        std::cout << "Remove vertices" << std::endl;
        std::vector<uint32_t> not_remove(graph.number_of_nodes(), 1);
        std::vector<uint32_t> new_nums(graph.number_of_nodes() + 1, 1);
        new_nums.back() = 0;
        for (auto n : to_remove) {
                not_remove[n] = 0;
                new_nums[n] = 0;
        }

        for (uint32_t node = 0; node < graph.number_of_nodes(); ++node)
                new_nums[node + 1] += new_nums[node];

        std::vector<std::vector<std::pair<uint64_t, uint64_t>>> adj_list(graph.number_of_nodes());
        size_t edgs = 0;
        for (uint32_t node = 0; node < graph.number_of_nodes(); ++node) {
                if (not_remove[node] == 0) {
                        continue;
                }
                uint32_t new_num = new_nums[node] - 1;
                forall_out_edges(graph, e, node) {
                        NodeID target = graph.getEdgeTarget(e);
                        EdgeWeight weight = graph.getEdgeWeight(e);

                        if (not_remove[target] == 1) {
                                ALWAYS_ASSERT(new_num != new_nums[target] - 1);
                                adj_list[new_num].emplace_back(new_nums[target] - 1, weight);
                                ++edgs;
                        }

                } endfor
        }

        std::cout << "Write graph" << std::endl;
        std::ofstream out("coarse_graph_without_high_degree");
        out << new_nums.back() << ' ' << edgs / 2 << " 11\n";
        for (uint32_t node = 0; node < graph.number_of_nodes(); ++node) {
                if (not_remove[node] == 0) {
                        continue;
                }
                out << graph.getNodeWeight(node) << ' ';

                uint32_t new_num = new_nums[node] - 1;
                std::sort(adj_list[new_num].begin(), adj_list[new_num].end());
                for (const auto& edge : adj_list[new_num]) {
                        out << edge.first + 1 << ' ' << edge.second << ' ';
                }
                out << '\n';
        }

        out.close();
}