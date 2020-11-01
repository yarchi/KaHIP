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

void convert_to_mtx(graph_access& G, const std::string& output) {
        ALWAYS_ASSERT(G.number_of_nodes() > 0);


        std::ofstream out(output);
        out << "%symmetric\n";
        out << G.number_of_nodes() << ' ' << G.number_of_nodes() << ' ' << G.number_of_edges() / 2 << '\n';
        for (size_t node = 0; node < G.number_of_nodes(); ++node) {
                forall_out_edges(G, e, node) {
                        NodeID neighbor = G.getEdgeTarget(e);
                        out << node + 1 << ' ' << neighbor +  1 << '\n';
                } endfor
        }
        std::cout << "Done conversion" << std::endl;
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

        convert_to_mtx(G, partition_config.filename_output);

        return 0;
}
