#pragma once

#include "partition/uncoarsening/refinement/kway_graph_refinement/multitry_kway_fm.h"

#include "partition/uncoarsening/refinement/parallel_kway_graph_refinement/kway_graph_refinement_commons.h"

#include <tbb/concurrent_queue.h>

class multitry_kway_fm;

namespace parallel {

class thread_data_factory {
public:
        using thread_data_refinement_core = parallel::thread_data_refinement_core;

        thread_data_factory(PartitionConfig& config,
                            graph_access& G,
                            complete_boundary& boundary)
                :       time_setup_start_nodes(0.0)
                ,       time_local_search(0.0)
                ,       time_init(0.0)
                ,       time_generate_moves(0.0)
                ,       time_move_nodes(0.0)
                ,       m_config(config)
                ,       m_G(G)
                ,       m_boundary(boundary)
                ,       m_moved_idx(G.number_of_nodes())
                ,       m_parts_weights(config.k)
                ,       m_parts_sizes(config.k)
                ,       m_moved_count(config.num_threads)
                ,       m_moved_idx_counter(0)
        {
                for (PartitionID block = 0; block < G.get_partition_count(); ++block) {
                        m_parts_weights[block].get().store(boundary.getBlockWeight(block), std::memory_order_relaxed);
                        m_parts_sizes[block].get().store(boundary.getBlockNoNodes(block), std::memory_order_relaxed);
                }

                m_thread_data.reserve(config.num_threads);

                for (uint32_t id = 0; id < config.num_threads; ++id) {
                        m_thread_data.emplace_back(id,
                                                 id,
                                                 config,
                                                 G,
                                                 boundary,
                                                 m_moved_idx,
                                                 m_parts_weights,
                                                 m_parts_sizes,
                                                 m_moved_count,
                                                 m_moved_idx_counter);
                }

        }


        void reset_global_data() {
                for (PartitionID block = 0; block < m_G.get_partition_count(); ++block) {
                        m_parts_weights[block].get().store(m_boundary.getBlockWeight(block), std::memory_order_relaxed);
                        m_parts_sizes[block].get().store(m_boundary.getBlockNoNodes(block), std::memory_order_relaxed);
                }

                for (uint32_t id = 0; id < m_config.num_threads; ++id) {
                         m_moved_count[id].get().store(0, std::memory_order_relaxed);
                }

                m_moved_idx_counter.store(0, std::memory_order_relaxed);
                queue.clear();
        }

        thread_data_refinement_core& get_thread_data(uint32_t id) {
                return m_thread_data[id].get();
        }

        const thread_data_refinement_core& get_thread_data(uint32_t id) const {
                return m_thread_data[id].get();
        }

        Cvector<thread_data_refinement_core>& get_all_threads_data() {
                return m_thread_data;
        }

        const Cvector<thread_data_refinement_core>& get_all_threads_data() const {
                return m_thread_data;
        }

        void print_iteration_statistics() {
                statistics_type stat;
                
                std::cout << "Time setup start nodes\t" << time_setup_start_nodes << " s" << std::endl;
                std::cout << "Time local search\t" << time_local_search << " s" << std::endl;

                std::cout << "Time init\t" << time_init << " s" << std::endl;
                std::cout << "Time generate moves\t" << time_generate_moves << " s" << std::endl;
                std::cout << "Time move nodes\t" << time_move_nodes << " s" << std::endl;

                stat.time_setup_start_nodes = time_setup_start_nodes;
                stat.time_local_search = time_local_search;
                stat.time_init = time_init;
                stat.time_generate_moves = time_generate_moves;
                stat.time_move_nodes = time_move_nodes;


                uint32_t total_tried_movements = 0;
                uint32_t total_accepted_movements = 0;
                uint32_t total_scaned_neighbours = 0;

                double total = 0.0;
                double total_tried = 0.0;
                double total_accepted = 0.0;
                double total_unroll = 0.0;
                for (uint32_t id = 0; id < m_config.num_threads; ++id) {
                        std::cout << "proc_id\t" << id << " | "
                                  << "time\t" << m_thread_data[id].get().total_thread_time << " s | "
                                  << "tried moves\t" << m_thread_data[id].get().tried_movements << " | "
                                  << "accepted moves\t" << m_thread_data[id].get().accepted_movements << " | "
                                  << "scanned neighbours\t" << m_thread_data[id].get().scaned_neighbours << " | "
                                  << "try moves time\t" << m_thread_data[id].get().total_thread_try_move_time << " s | "
                                  << "accepted moves time\t" << m_thread_data[id].get().total_thread_accepted_move_time << " s | "
                                  << "unroll moves time\t" << m_thread_data[id].get().total_thread_unroll_move_time << " s | "
                                  << "move nodes time\t" << m_thread_data[id].get().time_move_nodes << " s | "
                                  << "transpositions size\t" << m_thread_data[id].get().transpositions_size
                                  << std::endl;

                        total_tried_movements += m_thread_data[id].get().tried_movements;
                        total_accepted_movements += m_thread_data[id].get().accepted_movements;
                        total_scaned_neighbours += m_thread_data[id].get().scaned_neighbours;

                        statistics_type::proc_stat proc_stat;
                        proc_stat.proc_id = id;
                        proc_stat.total_thread_time = m_thread_data[id].get().total_thread_time;
                        proc_stat.tried_movements = m_thread_data[id].get().tried_movements;
                        proc_stat.accepted_movements = m_thread_data[id].get().accepted_movements;
                        proc_stat.scaned_neighbours = m_thread_data[id].get().scaned_neighbours;
                        proc_stat.total_thread_try_move_time = m_thread_data[id].get().total_thread_try_move_time;
                        proc_stat.total_thread_accepted_move_time = m_thread_data[id].get().total_thread_accepted_move_time;
                        proc_stat.total_thread_unroll_move_time = m_thread_data[id].get().total_thread_unroll_move_time;
                        
                        stat.proc_stats.push_back(proc_stat);
                        
                        total += m_thread_data[id].get().total_thread_time;
                        total_tried += m_thread_data[id].get().total_thread_try_move_time;
                        total_accepted += m_thread_data[id].get().total_thread_accepted_move_time;
                        total_unroll += m_thread_data[id].get().total_thread_unroll_move_time;
                }

                stat.total_tried_movements = total_tried_movements;
                stat.total_accepted_movements = total_accepted_movements;
                stat.total_scanned_neighbours = total_scaned_neighbours;

                std::cout << "Total tried moves\t" << total_tried_movements << std::endl;
                std::cout << "Total accepted movse\t" << total_accepted_movements << std::endl;
                std::cout << "Total scanned neighbours\t" << total_scaned_neighbours << std::endl;

                std::cout << "Average TIME per thread\t" << total / m_config.num_threads << " s" << std::endl;
                std::cout << "Average TIME tried moves per thread\t" << total_tried / m_config.num_threads << " s" << std::endl;
                std::cout << "Average TIME accepted moves per thread\t" << total_accepted / m_config.num_threads << " s" << std::endl;
                std::cout << "Average TIME unroll per thread\t" << total_unroll / m_config.num_threads << " s" << std::endl;

                stat.avg_thread_time = total / m_config.num_threads;
                stat.avg_tried = total_tried / m_config.num_threads;
                stat.avg_accepted = total_accepted / m_config.num_threads;
                stat.avg_unroll = total_unroll / m_config.num_threads;

                m_statistics.push_back(stat);
        }

        static void print_full_statistics() {
                if (m_statistics.empty())
                        return;

                uint32_t num_threads = m_statistics.front().proc_stats.size();
                statistics_type stat;
                stat.proc_stats.resize(num_threads);

                for (auto& st : m_statistics) {
                        stat += st;
                }

                std::cout << "Time setup start nodes\t" << stat.time_setup_start_nodes << " s" << std::endl;
                std::cout << "Time local search\t" << stat.time_local_search << " s" << std::endl;

                std::cout << "Time init\t" << stat.time_init << " s" << std::endl;
                std::cout << "Time generate moves\t" << stat.time_generate_moves << " s" << std::endl;
                std::cout << "Time move nodes\t" << stat.time_move_nodes << " s" << std::endl;

                for (auto& pr : stat.proc_stats) {
                        std::cout << "proc_id\t" << pr.proc_id << " | "
                                  << "time\t" << pr.total_thread_time << " s | "
                                  << "tried moves\t" << pr.tried_movements << " | "
                                  << "accepted moves\t" << pr.accepted_movements << " | "
                                  << "scanned neighbours\t" << pr.scaned_neighbours << " | "
                                  << "try moves time\t" << pr.total_thread_try_move_time << " s | "
                                  << "accepted moves time\t" << pr.total_thread_accepted_move_time << " s | "
                                  << "unroll moves time\t" << pr.total_thread_unroll_move_time << " s | "
                                  << std::endl;
                }

                std::cout << "Total tried moves\t" << stat.total_tried_movements << std::endl;
                std::cout << "Total accepted movse\t" << stat.total_accepted_movements << std::endl;
                std::cout << "Total scanned neighbours\t" << stat.total_scanned_neighbours << std::endl;

                std::cout << "Average TIME per thread\t" << stat.avg_thread_time << " s" << std::endl;
                std::cout << "Average TIME tried moves per thread\t" << stat.avg_tried << " s" << std::endl;
                std::cout << "Average TIME accepted moves per thread\t" << stat.avg_accepted << " s" << std::endl;
                std::cout << "Average TIME unroll per thread\t" << stat.avg_unroll << " s" << std::endl;
        }

        virtual ~thread_data_factory() {
                print_iteration_statistics();
        }

        tbb::concurrent_queue <NodeID> queue;

        double time_setup_start_nodes;
        double time_local_search;

        double time_init;
        double time_generate_moves;
        double time_move_nodes;
private:
        struct statistics_type {
                double time_setup_start_nodes = 0.0;
                double time_local_search = 0.0;
                double time_init = 0.0;
                double time_generate_moves = 0.0;
                double time_move_nodes = 0.0;

                double avg_thread_time = 0.0;
                double avg_tried = 0.0;
                double avg_accepted = 0.0;
                double avg_unroll = 0.0;

                uint32_t total_tried_movements = 0;
                uint32_t total_accepted_movements = 0;
                uint32_t total_scanned_neighbours = 0;

                struct proc_stat {
                        uint32_t proc_id = 0;
                        double total_thread_time = 0.0;
                        uint32_t tried_movements = 0;
                        uint32_t accepted_movements = 0;
                        uint32_t scaned_neighbours = 0;
                        double total_thread_try_move_time = 0.0;
                        double total_thread_accepted_move_time = 0.0;
                        double total_thread_unroll_move_time = 0.0;

                        proc_stat& operator+= (const proc_stat& ps) {
                                proc_id = ps.proc_id;

                                total_thread_time += ps.total_thread_time;
                                tried_movements += ps.tried_movements;
                                accepted_movements += ps.accepted_movements;
                                scaned_neighbours += ps.scaned_neighbours;
                                total_thread_try_move_time += ps.total_thread_try_move_time;
                                total_thread_accepted_move_time += ps.total_thread_accepted_move_time;
                                total_thread_unroll_move_time += ps.total_thread_unroll_move_time;

                                return *this;
                        }
                };
                
                std::vector<proc_stat> proc_stats;

                statistics_type operator+= (const statistics_type& stat) {
                        time_setup_start_nodes += stat.time_setup_start_nodes;
                        time_local_search += stat.time_local_search;
                        time_init += stat.time_init;
                        time_generate_moves += stat.time_generate_moves;
                        time_move_nodes += stat.time_move_nodes;

                        total_tried_movements += stat.total_tried_movements;
                        total_accepted_movements += stat.total_accepted_movements;
                        total_scanned_neighbours += stat.total_scanned_neighbours;

                        avg_thread_time += stat.avg_thread_time;
                        avg_tried += stat.avg_tried;
                        avg_accepted += stat.avg_accepted;
                        avg_unroll += stat.avg_unroll;

                        ALWAYS_ASSERT(proc_stats.size() == stat.proc_stats.size());
                        for (size_t i = 0; i < proc_stats.size(); ++i) {
                                proc_stats[i] += stat.proc_stats[i];
                        }
                        return *this;
                }
        };
        
        static std::vector<statistics_type> m_statistics;
        Cvector<thread_data_refinement_core> m_thread_data;

        PartitionConfig& m_config;
        graph_access& m_G;
        complete_boundary& m_boundary;

        // global data
        std::vector <AtomicWrapper<bool>> m_moved_idx;
        Cvector <AtomicWrapper<NodeWeight>> m_parts_weights;
        Cvector <AtomicWrapper<NodeWeight>> m_parts_sizes;
        Cvector <AtomicWrapper<int>> m_moved_count;
        AtomicWrapper<uint32_t> m_moved_idx_counter;
};

class multitry_kway_fm : public ::multitry_kway_fm {
public:
        using thread_data_refinement_core = parallel::thread_data_refinement_core;

        multitry_kway_fm(PartitionConfig& config, graph_access& G, complete_boundary& boundary)
                :       m_factory(config, G, boundary)
        {}

        static void print_full_statistics() {
                thread_data_factory::print_full_statistics();
        }

        virtual ~multitry_kway_fm() {

        }

        virtual int perform_refinement(PartitionConfig& config, graph_access& G,
                                       complete_boundary& boundary, unsigned rounds,
                                       bool init_neighbors, unsigned alpha) override;

        virtual int perform_refinement_around_parts(PartitionConfig& config,
                                            graph_access& G,
                                            complete_boundary& boundary,
                                            bool init_neighbors,
                                            unsigned alpha,
                                            PartitionID& lhs, PartitionID& rhs,
                                            std::unordered_map <PartitionID, PartitionID>& touched_blocks) override;

private:
        thread_data_factory m_factory;

        int start_more_locallized_search(PartitionConfig& config, graph_access& G,
                                             complete_boundary& boundary,
                                             bool init_neighbors,
                                             bool compute_touched_blocks,
                                             std::unordered_map<PartitionID, PartitionID>& touched_blocks,
                                             std::vector<NodeID>& todolist);
};

}