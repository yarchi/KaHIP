/******************************************************************************
 * kaffpa.cpp
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
#include <iostream>
#include <argtable2.h>
#include <fstream>
#include <queue>
#include <string.h>
#include <regex.h>
#include <stack>

#include "data_structure/graph_access.h"
#include "graph_io.h"
#include "macros_assertions.h"
#include "parse_parameters.h"
#include "partition/partition_config.h"
#include "timer.h"
#include "partition/uncoarsening/refinement/parallel_kway_graph_refinement/kway_graph_refinement_commons.h"

void bfs_traversal_to_file(graph_access& G, const std::string& output) {
        using vertex_type = NodeID;
        ALWAYS_ASSERT(G.number_of_nodes() > 0);

        std::vector<parallel::ht_query> ops;
        ops.reserve(G.number_of_edges());

        std::vector<uint8_t> visited(G.number_of_nodes());

        size_t num_comp = 0;
        for (size_t node = 0; node < G.number_of_nodes(); ++node) {
                if (visited[node]) {
                        continue;
                }
                ++num_comp;
                std::queue<vertex_type> queue;
                queue.push(node);
                while (!queue.empty()) {
                        vertex_type v = queue.front();
                        queue.pop();
                        ops.push_back({v, parallel::query_type::FIND});
                        if (visited[v]) {
                                continue;
                        }
                        visited[v] = true;

                        forall_out_edges(G, e, v)
                                        {
                                                NodeID target = G.getEdgeTarget(e);
                                                if (!visited[target]) {
                                                        queue.push(target);
                                                }
                                        }
                        endfor
                }
        }
        for (size_t node = 0; node < G.number_of_nodes(); ++node) {
                if (!visited[node]) {
                        std::cout << "Not visited: " << node << std::endl;
                }
        }
        std::cout << "Done BFS" << std::endl;
        std::cout << "Num comp " << num_comp << std::endl;
        std::cout << "Size " << ops.size() << std::endl;
        std::cout << "Correct " << (ops.size() - num_comp == G.number_of_edges() / 2) << std::endl;
        std::ofstream out(output);
        for (const auto& rec : ops) {
                out << rec.key << ' ' << rec.query << std::endl;
        }
        //out.close();
        std::cout << "Done writing to file" << std::endl;
}

void dfs_traversal_to_file(graph_access& G, const std::string& output) {
        using vertex_type = NodeID;

        std::stack<vertex_type> st_neighbors;
        ALWAYS_ASSERT(G.number_of_nodes() > 0);
        std::cout << "Start DFS" << std::endl;

        std::vector<parallel::ht_query> ops;
        ops.reserve(G.number_of_edges());
        std::vector<uint8_t> visited(G.number_of_nodes());

        size_t num_comp = 0;
        for (size_t node = 0; node < G.number_of_nodes(); ++node) {
                if (visited[node]) {
                        continue;
                }
                ++num_comp;
                std::stack<vertex_type> st;
                st.push(node);
                while (!st.empty()) {
                        vertex_type v = st.top();
                        st.pop();
                        ops.push_back({v, parallel::query_type::FIND});
                        if (visited[v]) {
                                continue;
                        }
                        visited[v] = true;

                        ALWAYS_ASSERT(st_neighbors.empty());
                        forall_out_edges(G, e, v)
                                        {
                                                NodeID target = G.getEdgeTarget(e);
                                                if (!visited[target]) {
                                                        st_neighbors.push(target);
                                                }
                                        }
                        endfor
                        while (!st_neighbors.empty()) {
                                NodeID target = st_neighbors.top();
                                st_neighbors.pop();
                                st.push(target);
                        }
                }
        }
        std::cout << "Num comp " << num_comp << std::endl;
        std::cout << "Size " << ops.size() << std::endl;
        std::cout << "Correct " << (ops.size() - num_comp == G.number_of_edges() / 2) << std::endl;
        std::cout << "Done DFS" << std::endl;
        std::ofstream out(output);
        for (const auto& rec : ops) {
                out << rec.key << ' ' << rec.query << std::endl;
        }
        out.close();
        std::cout << "Done writing to file" << std::endl;
}

void two_hop_traversal(graph_access& G, const std::string& output) {
        ALWAYS_ASSERT(G.number_of_nodes() > 0);

        std::vector<parallel::ht_query> ops;
        ops.reserve(G.number_of_edges());

        for (size_t node = 0; node < G.number_of_nodes(); ++node) {
                forall_out_edges(G, e, node) {
                        NodeID neighbor = G.getEdgeTarget(e);
                        forall_out_edges(G, e, neighbor) {
                                NodeID neighbor_neighbor = G.getEdgeTarget(e);
                                ops.push_back({neighbor_neighbor, parallel::query_type::FIND});
                                if (ops.size() > 500'000'000) {
                                        break;
                                }
                        }endfor
                        if (ops.size() > 500'000'000) {
                                break;
                        }
                } endfor
                if (ops.size() > 500'000'000) {
                        break;
                }
        }
        std::cout << "Done all neighbors traversal" << std::endl;
        std::cout << "Size " << ops.size() << std::endl;
        std::ofstream out(output);
        for (const auto& rec : ops) {
                out << rec.key << ' ' << rec.query << std::endl;
        }
        out.close();
        std::cout << "Done writing to file" << std::endl;
}

int main(int argn, char** argv) {
        std::cout << "Git revision\t" << GIT_DESC << std::endl;
        PartitionConfig partition_config;
        std::string graph_filename;

        bool is_graph_weighted = false;
        bool suppress_output = false;
        bool recursive = false;

        int ret_code = parse_parameters(argn, argv,
                                        partition_config,
                                        graph_filename,
                                        is_graph_weighted,
                                        suppress_output, recursive);

        if (ret_code) {
                return 0;
        }

        partition_config.LogDump(stdout);
        graph_access G;
        ALWAYS_ASSERT(partition_config.main_core == 0);

        timer t;
        graph_io::readGraphWeighted(G, graph_filename);
        std::cout << "io time: " << t.elapsed() << std::endl;
        std::cout << "N = " << G.number_of_nodes() << std::endl;
        std::cout << "M = " << G.number_of_edges() << std::endl;

        if (partition_config.graph_algo_type == "bfs") {
                bfs_traversal_to_file(G, partition_config.filename_output);
        } else if (partition_config.graph_algo_type == "dfs") {
                dfs_traversal_to_file(G, partition_config.filename_output);
        } else if (partition_config.graph_algo_type == "2_hop") {
                two_hop_traversal(G, partition_config.filename_output);
        } else {
                std::cout << "incorrect algo type" << std::endl;
                abort();
        }


        return 0;
}
