#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

#include "data_structure/parallel/thread_pool.h"

void parse(std::string& line, size_t& a, size_t& b) {
        size_t* cur = &a;
        for (size_t i = 0; i < line.size(); ++i) {
                if (line[i] != ' ') {
                        *cur = *cur * 10 + (line[i] - '0');
                } else {
                        cur = &b;
                }
        }
}

std::vector<std::vector<uint32_t>> read_edge_list(const std::string& path) {
        std::cout << "Read edge list" << std::endl;
        std::ifstream f(path);

        if (!f)
        {
                std::cout << "File not found!" << std::endl;
                exit(0);
        }

        std::string line;
        std::getline(f, line);

        size_t num_nodes = std::stoull(line);

        size_t parsed = 0;
        size_t gb = 0;
        std::vector<std::vector<uint32_t>> adj_list(num_nodes);
        while (std::getline(f, line)) {
                if (parsed / (1024 * 1024 * 1024ull) >= gb + 5) {
                        gb = parsed / (1024 * 1024 * 1024ull);
                        std::cout << gb << " GB" << std::endl;
                }

                parsed += line.size() + 1;

                size_t v = 0, u = 0;
                parse(line, v, u);
                if (v == 0 || u == 0) {
                        std::cout << v << ' ' << u << std::endl;
                        std::cout << "node equals zero" << std::endl;
                        exit(0);
                }
                v--;
                u--;
                if (v != u) {
                        adj_list[v].push_back(u);
                        adj_list[u].push_back(v);
                }
                // since two lines are the same in input file (a bug to fix)
                std::getline(f, line);
                parsed += line.size() + 1;
        }

        return adj_list;
}

void sort_edges(std::vector<std::vector<uint32_t>>& adj_list) {
        std::cout << "Sort adjacency list" << std::endl;
        std::atomic<size_t> offset(0);
        auto sort = [&]() {
                size_t num = offset.fetch_add(1000);
                while (num < adj_list.size()) {
                        for (size_t i = num; i < std::min(num + 1000, adj_list.size()); ++i) {
                                std::sort(adj_list[i].begin(), adj_list[i].end());
                                auto it = std::unique(adj_list[i].begin(), adj_list[i].end());
                                adj_list[i].resize(it - adj_list[i].begin());
                        }
                        num = offset.fetch_add(1000);
                }
        };
        parallel::submit_for_all(sort);
}

void write_metis(const std::string& out_path, std::vector<std::vector<uint32_t>>& adj_list) {
        std::cout << "Write adjacency list" << std::endl;
        std::ofstream out(out_path);

        size_t num_edges = 0;
        for (const auto& edges : adj_list) {
                num_edges += edges.size();
        }

        out << adj_list.size() << ' ' << num_edges / 2 << std::endl;
        for (const auto& edges : adj_list) {
                for (auto edge : edges) {
                        out << edge + 1 << ' ';
                }
                out << std::endl;
        }

}

int main(int argc, const char* argv[]) {
        parallel::PinToCore(0);
        parallel::g_thread_pool.Resize(std::stoi(argv[3]) - 1);
        auto adj_list = read_edge_list(argv[1]);
        sort_edges(adj_list);
        write_metis(argv[2], adj_list);
        return 0;
}
