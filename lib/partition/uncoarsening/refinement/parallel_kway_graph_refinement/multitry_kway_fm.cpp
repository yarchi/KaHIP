namespace parallel {

int multitry_kway_fm::perform_refinement(PartitionConfig& config, graph_access& G,
                                             complete_boundary& boundary, unsigned rounds,
                                             bool init_neighbors, unsigned alpha) {

        commons = kway_graph_refinement_commons::getInstance(config);

        unsigned tmp_alpha = config.kway_adaptive_limits_alpha;
        KWayStopRule tmp_stop = config.kway_stop_rule;
        config.kway_adaptive_limits_alpha = alpha;
        config.kway_stop_rule = KWAY_ADAPTIVE_STOP_RULE;

        int overall_improvement = 0;
        for (unsigned i = 0; i < rounds; i++) {
                boundary_starting_nodes start_nodes;
                boundary.setup_start_nodes_all(G, start_nodes);
                if (start_nodes.size() == 0) {
                        return 0;
                }// nothing to refine

                std::unordered_map <PartitionID, PartitionID> touched_blocks;
                EdgeWeight improvement = start_more_locallized_search(config, G, boundary,
                                                                          init_neighbors, false, touched_blocks,
                                                                          start_nodes);
                if (improvement == 0) break;
                overall_improvement += improvement;

        }

        ASSERT_TRUE(overall_improvement >= 0);

        config.kway_adaptive_limits_alpha = tmp_alpha;
        config.kway_stop_rule = tmp_stop;

        return (int) overall_improvement;

}

int multitry_kway_fm::perform_refinement_around_parts(PartitionConfig& config, graph_access& G,
                                                          complete_boundary& boundary, bool init_neighbors,
                                                          unsigned alpha,
                                                          PartitionID& lhs, PartitionID& rhs,
                                                          std::unordered_map <PartitionID, PartitionID>& touched_blocks) {
        commons = kway_graph_refinement_commons::getInstance(config);

        unsigned tmp_alpha = config.kway_adaptive_limits_alpha;
        KWayStopRule tmp_stop = config.kway_stop_rule;
        config.kway_adaptive_limits_alpha = alpha;
        config.kway_stop_rule = KWAY_ADAPTIVE_STOP_RULE;
        int overall_improvement = 0;

        for (unsigned i = 0; i < config.local_multitry_rounds; i++) {
                boundary_starting_nodes start_nodes;
                boundary.setup_start_nodes_around_blocks(G, lhs, rhs, start_nodes);

                if (start_nodes.size() == 0) { return 0; }// nothing to refine

                EdgeWeight improvement = start_more_locallized_search(config, G, boundary, init_neighbors, true,
                                                                          touched_blocks, start_nodes);
                if (improvement == 0) break;

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
                                                       std::unordered_map <PartitionID, PartitionID>& touched_blocks,
                                                       std::vector <NodeID>& todolist) {
        uint32_t num_threads = config.num_threads;
        commons = kway_graph_refinement_commons::getInstance(config);
        parallel::kway_graph_refinement_core refinement_core;
        int local_step_limit = 0;
        int upper_bound_gain_improvement = 0;

        tbb::concurrent_queue <NodeID> queue;
        while (!todolist.empty()) {
                size_t random_idx = random_functions::nextInt(0, todolist.size() - 1);
                NodeID node = todolist[random_idx];
                queue.push(node);

                std::swap(todolist[random_idx], todolist.back());
                todolist.pop_back();
        }

        m_factory.update_partitions_weights();

        //CLOCK_START;
        std::vector <AtomicWrapper<bool>> moved_idx(G.number_of_nodes());
        //CLOCK_END("Init move_nodes");

        //CLOCK_START_N;
        Cvector <AtomicWrapper<int>> moved_count(num_threads);
        Cvector <std::vector<NodeID>> transpositions(num_threads);
        Cvector <std::vector<PartitionID>> from_partitions(num_threads);
        Cvector <std::vector<PartitionID>> to_partitions(num_threads);
        Cvector <std::vector<EdgeWeight>> gains(num_threads);
        Cvector <std::vector<std::pair<PartitionID, uint32_t>>> min_cut_indices(num_threads);
        Cvector <AtomicWrapper<NodeWeight>> parts_weights(config.k);
        Cvector <AtomicWrapper<NodeWeight>> parts_sizes(config.k);

        for (PartitionID block = 0; block < G.get_partition_count(); ++block) {
                parts_weights[block].get().store(boundary.getBlockWeight(block), std::memory_order_relaxed);
                parts_sizes[block].get().store(boundary.getBlockNoNodes(block), std::memory_order_relaxed);
        }

        std::vector<parallel::thread_data_refinement_core> threads_data;
        threads_data.reserve(config.k);
        //CLOCK_END_N("Init auxilary vectors");

        auto task = [&](uint32_t id) {
                NodeID node;

                boundary_starting_nodes real_start_nodes;
                parallel::thread_data_refinement_core::nodes_partitions_hash_table nodes_partitions(
                        parallel::thread_data_refinement_core::nodes_partitions_hash_table::get_max_size_to_fit_l1());

                parallel::thread_data_refinement_core td(
                        id,
                        id,
                        config,
                        G,
                        boundary,
                        real_start_nodes,
                        local_step_limit,
                        moved_idx,
                        compute_touched_blocks,
                        touched_blocks,
                        nodes_partitions,
                        transpositions[id].get(),
                        from_partitions[id].get(),
                        to_partitions[id].get(),
                        gains[id].get(),
                        parts_weights,
                        parts_sizes,
                        min_cut_indices[id].get()
                );


                td.transpositions.reserve(100);
                td.from_partitions.reserve(100);
                td.to_partitions.reserve(100);
                td.gains.reserve(100);
                td.min_cut_indices.reserve(100);

                while (queue.try_pop(node)) {
                        PartitionID maxgainer;
                        EdgeWeight extdeg = 0;
                        PartitionID from = td.get_local_partition(node);
                        td.compute_gain(node, from, maxgainer, extdeg);

                        if (!moved_idx[node].load(std::memory_order_relaxed) && extdeg > 0) {
                                real_start_nodes.clear();
                                real_start_nodes.reserve(G.getNodeDegree(node) + 1);
                                real_start_nodes.push_back(node);

                                if (init_neighbors) {
                                        forall_out_edges(G, e, node)
                                        {
                                                NodeID target = G.getEdgeTarget(e);
                                                if (!moved_idx[target].load(std::memory_order_relaxed)) {
                                                        extdeg = 0;
                                                        td.compute_gain(target, from, maxgainer, extdeg);
                                                        if (extdeg > 0) {
                                                                real_start_nodes.push_back(target);
                                                        }
                                                }
                                        }
                                        endfor
                                }

                                int improvement = 0;
                                int movement = 0;
                                int min_cut_index = 0;
                                if (compute_touched_blocks) {
                                        std::tie(improvement, movement, min_cut_index) =
                                                refinement_core.single_kway_refinement_round_par(td, touched_blocks);
                                        if (improvement < 0) {
                                                std::cout << "buf error improvement < 0" << std::endl;
                                        }
                                } else {
                                        std::tie(improvement, movement, min_cut_index) =
                                                refinement_core.single_kway_refinement_round_par(td);
                                        if (improvement < 0) {
                                                std::cout << "buf error improvement < 0" << std::endl;
                                        }
                                }

                                td.upper_bound_gain_improvement += improvement;

                                ALWAYS_ASSERT(transpositions[id].get().size() > 0);
                                td.min_cut_indices.emplace_back(min_cut_index,
                                                                transpositions[id].get().size() - 1);
                                moved_count[id].get().fetch_add(movement, std::memory_order_relaxed);
                        }

                        int overall_movement = 0;
                        for (uint32_t id = 0; id < num_threads; ++id) {
                                int moved = moved_count[id].get().load(std::memory_order_relaxed);
                                overall_movement += moved;
                        }

                        if (overall_movement > 0.05 * G.number_of_nodes())
                                break;
                }
                return td;
        };

        //CLOCK_START_N;
        std::vector <std::future<parallel::thread_data_refinement_core>> futures;
        futures.reserve(num_threads - 1);

        for (uint32_t id = 1; id < num_threads; ++id) {
                futures.push_back(parallel::g_thread_pool.Submit(task, id));
        }

        threads_data.push_back(task(0));
        upper_bound_gain_improvement += threads_data.back().upper_bound_gain_improvement;
        for (uint32_t id = 1; id < num_threads; ++id) {
                threads_data.push_back(futures[id - 1].get());
                upper_bound_gain_improvement += threads_data.back().upper_bound_gain_improvement;

        }
        //CLOCK_END_N("Generate moves time");

        //std::cout << "Upperbound for gain improvement:\t" << upper_bound_gain_improvement << std::endl;

        int real_gain_improvement = 0;
        uint32_t real_nodes_movement = 0;
        //CLOCK_START_N;
        std::tie(real_gain_improvement, real_nodes_movement) = refinement_core.apply_moves(threads_data);
        //CLOCK_END_N("Move nodes time");

        //std::cout << "Real gain improvement\t" << real_gain_improvement << std::endl;
        //std::cout << "Number of nodes moved\t" << real_nodes_movement << std::endl;

        ALWAYS_ASSERT(real_gain_improvement >= 0);

        return real_gain_improvement;
}

}