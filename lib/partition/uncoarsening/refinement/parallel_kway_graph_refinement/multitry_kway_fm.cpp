#include "data_structure/parallel/thread_pool.h"
#include "data_structure/parallel/time.h"

#include "uncoarsening/refinement/parallel_kway_graph_refinement/kway_graph_refinement_core.h"
#include "uncoarsening/refinement/parallel_kway_graph_refinement/multitry_kway_fm.h"

#include <iomanip>
#include <sstream>

#ifdef __gnu_linux__
//#include "ittnotify.h"
#endif

namespace parallel {

std::vector<thread_data_factory::statistics_type> thread_data_factory::m_statistics;

int multitry_kway_fm::perform_refinement(PartitionConfig& config, graph_access& G,
                                         complete_boundary& boundary, unsigned rounds,
                                         bool init_neighbors, unsigned alpha) {
        unsigned tmp_alpha = config.kway_adaptive_limits_alpha;
        KWayStopRule tmp_stop = config.kway_stop_rule;
        config.kway_adaptive_limits_alpha = alpha;
        config.kway_stop_rule = KWAY_ADAPTIVE_STOP_RULE;
        int overall_improvement = 0;

        while (true) {
                CLOCK_START;

                setup_start_nodes_all(G, boundary);

                if (m_factory.queue.size() == 0) {
                        break;
                }

                m_factory.time_setup_start_nodes += CLOCK_END_TIME;


                CLOCK_START_N;
                std::unordered_map<PartitionID, PartitionID> touched_blocks;
                EdgeWeight improvement = start_more_locallized_search(config, G, boundary, init_neighbors,
                                                                      false,
                                                                      touched_blocks);

                m_factory.time_local_search += CLOCK_END_TIME;
                if (improvement == 0) {
                        break;
                }

                overall_improvement += improvement;
        }

        config.kway_adaptive_limits_alpha = tmp_alpha;
        config.kway_stop_rule = tmp_stop;
        ASSERT_TRUE(overall_improvement >= 0);
        return (int) overall_improvement;
}

int multitry_kway_fm::perform_refinement_around_parts(PartitionConfig& config, graph_access& G,
                                                      complete_boundary& boundary, bool init_neighbors,
                                                      unsigned alpha,
                                                      PartitionID& lhs, PartitionID& rhs,
                                                      std::unordered_map<PartitionID, PartitionID>& touched_blocks) {
        unsigned tmp_alpha = config.kway_adaptive_limits_alpha;
        KWayStopRule tmp_stop = config.kway_stop_rule;
        config.kway_adaptive_limits_alpha = alpha;
        config.kway_stop_rule = KWAY_ADAPTIVE_STOP_RULE;
        int overall_improvement = 0;

        for (unsigned i = 0; i < config.local_multitry_rounds; i++) {
                CLOCK_START;

                setup_start_nodes_around_blocks(G, boundary, lhs, rhs);

                if (m_factory.queue.size() == 0) {
                        break;
                }
                m_factory.time_setup_start_nodes += CLOCK_END_TIME;


                CLOCK_START_N;
                EdgeWeight improvement = start_more_locallized_search(config, G, boundary, init_neighbors, true,
                                                                      touched_blocks);

                m_factory.time_local_search += CLOCK_END_TIME;
                if (improvement == 0) {
                        break;
                }

                overall_improvement += improvement;
        }
        config.kway_adaptive_limits_alpha = tmp_alpha;
        config.kway_stop_rule = tmp_stop;
        ASSERT_TRUE(overall_improvement >= 0);
        return (int) overall_improvement;
}

int multitry_kway_fm::start_more_locallized_search(PartitionConfig& config, graph_access& G,
                                                   complete_boundary& boundary,
                                                   bool init_neighbors,
                                                   bool compute_touched_blocks,
                                                   std::unordered_map<PartitionID, PartitionID>& touched_blocks) {
        uint32_t num_threads = config.num_threads;
        parallel::kway_graph_refinement_core refinement_core;
        int local_step_limit = 50;

        CLOCK_START;
#ifdef COMPARE_WITH_SEQUENTIAL_KAHIP
        m_factory.m_config.num_threads = 1;
        num_threads = 1;

        std::vector<NodeID> todolist;
        todolist.reserve(m_factory.queue.size());
        NodeID node;
        while (m_factory.queue.try_pop(node, 0)) {
                todolist.push_back(node);
        }
        std::sort(todolist.begin(), todolist.end());
        auto it = std::unique(todolist.begin(), todolist.end());
        todolist.resize(it - todolist.begin());
        random_functions::permutate_vector_good(todolist, false);
#endif
        int total_gain_improvement = 0;

#ifdef COMPARE_WITH_SEQUENTIAL_KAHIP
        while (!todolist.empty()) {
#else
        // we need the external loop for move strategy when conflicted nodes are reactivated for the next
        // parallel phase
        while (!m_factory.queue.empty()) {
#endif
                auto task = [&](uint32_t id) {
                        CLOCK_START;
                        NodeID node;

                        auto& td = m_factory.get_thread_data(id);
                        td.reset_thread_data();

                        td.step_limit = local_step_limit;
                        uint32_t nodes_processed = 0;
                        bool res = false;
#ifdef COMPARE_WITH_SEQUENTIAL_KAHIP
                        while (!todolist.empty()) {
                                int random_idx = random_functions::nextInt(0, todolist.size() - 1);
                                NodeID node = todolist[random_idx];
#else
                        //__itt_resume();
                        while (m_factory.queue.try_pop(node, id)) {
#endif
                                // this change changes num part accesses since it changes random source
#ifdef COMPARE_WITH_SEQUENTIAL_KAHIP
#else
                                if (!td.moved_idx[node].load(std::memory_order_relaxed)) {
#endif
                                        PartitionID maxgainer;
                                        EdgeWeight extdeg = 0;
                                        PartitionID from = td.get_local_partition(node);
                                        td.compute_gain(node, from, maxgainer, extdeg);
#ifdef COMPARE_WITH_SEQUENTIAL_KAHIP
                                        if (!td.moved_idx[node].load(std::memory_order_relaxed) && extdeg > 0) {
#else
                                        if (extdeg > 0) {
#endif
                                                td.start_nodes.clear();
                                                td.start_nodes.reserve(G.getNodeDegree(node) + 1);
                                                td.start_nodes.push_back(node);

                                                if (init_neighbors) {
                                                        forall_out_edges(G, e, node) {
                                                                ++td.scaned_neighbours;
                                                                NodeID target = G.getEdgeTarget(e);
                                                                if (!td.moved_idx[target].load(std::memory_order_relaxed)) {
                                                                        extdeg = 0;
                                                                        PartitionID from = td.get_local_partition(target);
                                                                        td.compute_gain(target, from, maxgainer, extdeg);
                                                                        if (extdeg > 0) {
                                                                                td.start_nodes.push_back(target);
                                                                        }
                                                                }
                                                        } endfor
                                                }

                                                nodes_processed += td.start_nodes.size();

                                                int improvement = 0;
                                                int min_cut_index = 0;
                                                uint32_t tried_movements = 0;
                                                size_t moved_before = td.moved.size();

                                                std::tie(improvement, min_cut_index, tried_movements) =
                                                        refinement_core.single_kway_refinement_round(td);

                                                if (improvement < 0) {
                                                        std::cout << "buf error improvement < 0" << std::endl;
                                                        abort();
                                                }

                                                td.upper_bound_gain_improvement += improvement;

                                                ALWAYS_ASSERT(td.transpositions.size() > 0);
                                                td.min_cut_indices.emplace_back(min_cut_index, td.transpositions.size() - 1);
                                                if (!td.config.kway_all_boundary_nodes_refinement) {
                                                        td.moved_count[id].get().fetch_add(td.moved.size() - moved_before,
                                                                std::memory_order_relaxed);
                                                }
                                                td.tried_movements += tried_movements;
                                        }
#ifndef COMPARE_WITH_SEQUENTIAL_KAHIP
                                }
#endif

                                if (!td.config.kway_all_boundary_nodes_refinement) {
                                        int overall_movement = 0;
                                        for (uint32_t id = 0; id < num_threads; ++id) {
                                                int moved = td.moved_count[id].get().load(std::memory_order_relaxed);
                                                overall_movement += moved;
                                        }

                                        if (overall_movement > 0.05 * G.number_of_nodes()) {
                                                ++td.stop_faction_of_nodes_moved;
                                                res = true;
                                                break;
                                        }
                                }
#ifdef COMPARE_WITH_SEQUENTIAL_KAHIP
                                std::swap(todolist[random_idx], todolist[todolist.size() - 1]);
                                todolist.pop_back();
#endif
                        }
                        //__itt_pause();

                        td.num_threads_finished.fetch_add(1, std::memory_order_acq_rel);
                        td.total_thread_time += CLOCK_END_TIME;
//                        if (id > 0) {
//                                finished_threads.push(id);
//                        }
                        return res;
                };

                CLOCK_START_N;
                std::vector<std::future<bool>> futures;
                futures.reserve(num_threads - 1);

                for (uint32_t id = 1; id < num_threads; ++id) {
                        futures.push_back(parallel::g_thread_pool.Submit(id - 1, task, id));
                        //futures.push_back(parallel::g_thread_pool.Submit(task, id));
                }

                bool is_more_that_5percent_moved = task(0);

                {
                        CLOCK_START;
                        std::for_each(futures.begin(), futures.end(), [&](auto& future) {
                                if (future.get()) {
                                        is_more_that_5percent_moved = true;
                                }
                        });
                        m_factory.time_wait += CLOCK_END_TIME;
                }
                m_factory.time_generate_moves += CLOCK_END_TIME;

                std::vector<NodeID> reactivated_vertices;
                reactivated_vertices.reserve(100);

                int real_gain_improvement = 0;
                uint32_t real_nodes_movement = 0;
                CLOCK_START_N;
                std::tie(real_gain_improvement, real_nodes_movement) = refinement_core.apply_moves(num_threads,
                        m_factory.get_all_threads_data(), compute_touched_blocks, touched_blocks, reactivated_vertices);

                total_gain_improvement += real_gain_improvement;

                m_factory.partial_reset_global_data();

                if (config.kway_all_boundary_nodes_refinement && real_gain_improvement == 0) {
                        break;
                }

                m_factory.get_thread_data(0).rnd.shuffle(reactivated_vertices);
                for (auto vertex : reactivated_vertices) {
                        m_factory.queue.push(vertex);
                }

                m_factory.time_move_nodes += CLOCK_END_TIME;

                ALWAYS_ASSERT(real_gain_improvement >= 0);

                if (!config.kway_all_boundary_nodes_refinement && is_more_that_5percent_moved) {
                        break;
                }
        }
        m_factory.reset_global_data();
        ALWAYS_ASSERT(total_gain_improvement >= 0);
        return total_gain_improvement;
}

int multitry_kway_fm::start_more_locallized_search_experimental(PartitionConfig& config, graph_access& G,
                                                                complete_boundary& boundary,
                                                                bool init_neighbors,
                                                                bool compute_touched_blocks,
                                                                std::unordered_map<PartitionID, PartitionID>& touched_blocks,
                                                                std::vector<NodeID>& todolist) {
        uint32_t num_threads = config.num_threads;
        parallel::kway_graph_refinement_core refinement_core;
        int local_step_limit = 50;

        ALWAYS_ASSERT(config.apply_move_strategy != ApplyMoveStrategy::REACTIVE_VERTICES);

        m_factory.reset_global_data();

//        uint32_t not_moved_start = todolist.size();
//        std::atomic<uint32_t> not_moved_cur = not_moved_start;

        std::atomic_bool first_execution;
        first_execution.store(true, std::memory_order_relaxed);

        int total_gain_improvement = 0;
        while (!todolist.empty()) {
                //all experiments with this parameter:
                //const uint32_t batch_multiplicative_const = 10;
                const uint32_t batch_multiplicative_const = 1;
                uint32_t batch_size = batch_multiplicative_const * num_threads;
                //uint32_t batch_size = todolist.size();

                CLOCK_START;
                for (size_t i = 0; i < batch_size && !todolist.empty();) {
                        size_t random_idx = random_functions::nextInt(0, todolist.size() - 1);
                        NodeID node = todolist[random_idx];
                        if (first_execution.load(std::memory_order_relaxed) || !m_factory.is_moved(node)) {
                                m_factory.queue.push(node);
                                ++i;
                        }

                        std::swap(todolist[random_idx], todolist.back());
                        todolist.pop_back();
                }
                m_factory.time_init += CLOCK_END_TIME;


                auto task = [&](uint32_t id) {
                        CLOCK_START;
                        NodeID node;

                        auto& td = m_factory.get_thread_data(id);
                        if (first_execution.load(std::memory_order_relaxed)) {
                                td.reset_thread_data();
                        } else {
                                td.partial_reset_thread_data();
                        }

                        td.step_limit = local_step_limit;
                        uint32_t nodes_processed = 0;
                        while (m_factory.queue.try_pop(node, id)) {
                                PartitionID maxgainer;
                                EdgeWeight extdeg = 0;
                                PartitionID from = td.get_local_partition(node);
                                td.compute_gain(node, from, maxgainer, extdeg);

                                if (!td.moved_idx[node].load(std::memory_order_relaxed) && extdeg > 0) {
                                        td.start_nodes.clear();
                                        td.start_nodes.reserve(G.getNodeDegree(node) + 1);
                                        td.start_nodes.push_back(node);

                                        if (init_neighbors) {
                                                forall_out_edges(G, e, node) {
                                                        NodeID target = G.getEdgeTarget(e);
                                                        if (!td.moved_idx[target].load(std::memory_order_relaxed)) {
                                                                extdeg = 0;
                                                                PartitionID from = td.get_local_partition(target);
                                                                td.compute_gain(target, from, maxgainer, extdeg);
                                                                if (extdeg > 0) {
                                                                        td.start_nodes.push_back(target);
                                                                }
                                                        }
                                                } endfor
                                        }

                                        nodes_processed += td.start_nodes.size();

                                        int improvement = 0;
                                        int min_cut_index = 0;
                                        uint32_t tried_movements = 0;
                                        size_t moved_before = td.moved.size();

                                        std::tie(improvement, min_cut_index, tried_movements) =
                                                refinement_core.single_kway_refinement_round(td);
                                        if (improvement < 0) {
                                                std::cout << "buf error improvement < 0" << std::endl;
                                        }

                                        td.upper_bound_gain_improvement += improvement;

                                        ALWAYS_ASSERT(td.transpositions.size() > 0);
                                        td.min_cut_indices.emplace_back(min_cut_index, td.transpositions.size() - 1);
                                        td.moved_count[id].get().fetch_add(td.moved.size() - moved_before,
                                                                           std::memory_order_relaxed);

                                        td.tried_movements += tried_movements;
                                }

                                int overall_movement = 0;
                                for (uint32_t id = 0; id < num_threads; ++id) {
                                        int moved = td.moved_count[id].get().load(std::memory_order_relaxed);
                                        overall_movement += moved;
                                }

                                if (overall_movement > 0.05 * G.number_of_nodes()) {
                                        td.total_thread_time += CLOCK_END_TIME;
                                        ++td.stop_faction_of_nodes_moved;
                                        return true;
                                }
                        }
                        td.total_thread_time += CLOCK_END_TIME;
                        return false;
                };

                CLOCK_START_N;
                std::vector<std::future<bool>> futures;
                futures.reserve(num_threads - 1);

                for (uint32_t id = 1; id < num_threads; ++id) {
                        futures.push_back(parallel::g_thread_pool.Submit(id - 1, task, id));
                        //futures.push_back(parallel::g_thread_pool.Submit(task, id));
                }

                bool is_more_that_5percent_moved = task(0);
                std::for_each(futures.begin(), futures.end(), [&](auto& future) {
                        if (future.get()) {
                                is_more_that_5percent_moved = true;
                        }
                });

                m_factory.time_generate_moves += CLOCK_END_TIME;

                std::vector<NodeID> reactivated_vertices;

                int real_gain_improvement = 0;
                uint32_t real_nodes_movement = 0;
                CLOCK_START_N;
                std::tie(real_gain_improvement, real_nodes_movement) = refinement_core.apply_moves(num_threads,
                        m_factory.get_all_threads_data(), compute_touched_blocks, touched_blocks, reactivated_vertices);

                total_gain_improvement += real_gain_improvement;
                first_execution.store(false, std::memory_order_release);

                m_factory.partial_reset_global_data();

                m_factory.time_move_nodes += CLOCK_END_TIME;

                ALWAYS_ASSERT(real_gain_improvement >= 0);

                if (is_more_that_5percent_moved) {
                        break;
                }
        }
        ALWAYS_ASSERT(total_gain_improvement >= 0);
        return total_gain_improvement;
}

void multitry_kway_fm::setup_start_nodes_around_blocks(graph_access& G, complete_boundary& boundary, PartitionID & lhs,
                                                       PartitionID & rhs) {
        std::vector<PartitionID> lhs_neighbors;
        boundary.getNeighbors(lhs, lhs_neighbors);

        std::vector<PartitionID> rhs_neighbors;
        boundary.getNeighbors(rhs, rhs_neighbors);

        auto sub_task = [this, &G, &boundary](uint32_t thread_id, std::atomic<size_t>& counter, PartitionID part,
                                              const std::vector<PartitionID>& neighbour_parts) {
                auto& thread_container = m_factory.queue[thread_id];

                size_t offset = counter.fetch_add(1, std::memory_order_relaxed);

                while (offset < neighbour_parts.size()) {
                        PartitionID neighbor_part = neighbour_parts[offset];

                        const PartialBoundary& partial_boundary_part =
                                boundary.getDirectedBoundaryThreadSafe(part, part, neighbor_part);

                        for (const auto& elem : partial_boundary_part.internal_boundary) {
                                NodeID cur_bnd_node = elem.first;
                                ALWAYS_ASSERT(G.getPartitionIndex(cur_bnd_node) == part);
                                thread_container.push_back(cur_bnd_node);
                        }

                        const PartialBoundary& partial_boundary_neighbor_part =
                                boundary.getDirectedBoundaryThreadSafe(neighbor_part, part, neighbor_part);

                        for (const auto& elem : partial_boundary_neighbor_part.internal_boundary) {
                                NodeID cur_bnd_node = elem.first;
                                ALWAYS_ASSERT(G.getPartitionIndex(cur_bnd_node) == neighbor_part);
                                thread_container.push_back(cur_bnd_node);
                        }
                        offset = counter.fetch_add(1, std::memory_order_relaxed);
                }

                auto& td = m_factory.get_thread_data(thread_id);
                td.rnd.shuffle(thread_container.begin(), thread_container.end());
        };

        std::atomic<size_t> lhs_counter(0);
        std::atomic<size_t> rhs_counter(0);
        auto task = [&](uint32_t thread_id) {
                sub_task(thread_id, lhs_counter, lhs, lhs_neighbors);
                sub_task(thread_id, rhs_counter, rhs, rhs_neighbors);
        };

        std::vector<std::future<void>> futures;
        futures.reserve(parallel::g_thread_pool.NumThreads());

        for (uint32_t id = 0; id < parallel::g_thread_pool.NumThreads(); ++id) {
                futures.push_back(parallel::g_thread_pool.Submit(id, task, id + 1));
                //futures.push_back(parallel::g_thread_pool.Submit(task, id + 1));
        }

        task(0);
        std::for_each(futures.begin(), futures.end(), [](auto& future) {
                future.get();
        });

        shuffle_task_queue();
}

void multitry_kway_fm::setup_start_nodes_all(graph_access& G, complete_boundary& boundary) {
        QuotientGraphEdges quotient_graph_edges;
        boundary.getQuotientGraphEdges(quotient_graph_edges);

        auto sub_task = [this, &G, &boundary, &quotient_graph_edges](uint32_t thread_id, std::atomic<size_t>& counter) {
                auto& thread_container = m_factory.queue[thread_id];

#undef NO_DUPLICATES

#ifdef NO_DUPLICATES
                parallel::hash_set<NodeID> added_nodes(get_max_size_to_fit_l1<parallel::hash_set<NodeID>>());
#endif
                size_t offset = counter.fetch_add(1, std::memory_order_relaxed);

                while (offset < quotient_graph_edges.size()) {
                        boundary_pair& ret_value = quotient_graph_edges[offset];
                        PartitionID lhs = ret_value.lhs;
                        PartitionID rhs = ret_value.rhs;

                        auto& partial_boundary_lhs = boundary.getDirectedBoundaryThreadSafe(lhs, lhs, rhs);
                        for (const auto& elem : partial_boundary_lhs.internal_boundary) {
                                NodeID cur_bnd_node = elem.first;
                                ALWAYS_ASSERT(G.getPartitionIndex(cur_bnd_node) == lhs);
#ifdef NO_DUPLICATES
                                if (!added_nodes.contains(cur_bnd_node)) {
                                        thread_container.push_back(cur_bnd_node);
                                        added_nodes.insert(cur_bnd_node);
                                }
#else
                                thread_container.push_back(cur_bnd_node);
#endif
                        }

                        auto& partial_boundary_rhs = boundary.getDirectedBoundaryThreadSafe(rhs, lhs, rhs);
                        for (const auto& elem : partial_boundary_rhs.internal_boundary) {
                                NodeID cur_bnd_node = elem.first;
                                ALWAYS_ASSERT(G.getPartitionIndex(cur_bnd_node) == rhs);
#ifdef NO_DUPLICATES
                                if (!added_nodes.contains(cur_bnd_node)) {
                                        thread_container.push_back(cur_bnd_node);
                                        added_nodes.insert(cur_bnd_node);
                                }
#else
                                thread_container.push_back(cur_bnd_node);
#endif
                        }

                        offset = counter.fetch_add(1, std::memory_order_relaxed);
                }

                auto& td = m_factory.get_thread_data(thread_id);
                td.rnd.shuffle(thread_container.begin(), thread_container.end());
        };

        std::atomic<size_t> counter(0);
        auto task = [&](uint32_t thread_id) {
                sub_task(thread_id, counter);
        };

        std::vector<std::future<void>> futures;
        futures.reserve(parallel::g_thread_pool.NumThreads());

        for (uint32_t id = 0; id < parallel::g_thread_pool.NumThreads(); ++id) {
                futures.push_back(parallel::g_thread_pool.Submit(id, task, id + 1));
                //futures.push_back(parallel::g_thread_pool.Submit(task, id + 1));
        }

        task(0);
        std::for_each(futures.begin(), futures.end(), [](auto& future) {
                future.get();
        });

        shuffle_task_queue();
}

void  multitry_kway_fm::shuffle_task_queue() {
        auto& td = m_factory.get_thread_data(0);

        td.rnd.shuffle(m_factory.queue.begin(), m_factory.queue.end());
}

}
