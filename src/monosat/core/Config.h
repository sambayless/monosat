/**************************************************************************************************
 The MIT License (MIT)

 Copyright (c) 2014, Sam Bayless

 Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
 associated documentation files (the "Software"), to deal in the Software without restriction,
 including without limitation the rights to use, copy, modify, merge, publish, distribute,
 sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all copies or
 substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
 NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
 OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 **************************************************************************************************/
#ifndef CONFIG_H_
#define CONFIG_H_

#include "monosat/utils/Options.h"
#include "monosat/utils/System.h"

namespace Monosat {
extern BoolOption opt_show_version_and_quit;
extern IntOption opt_verb;
extern IntOption opt_verb_optimize;
extern BoolOption opt_pre;
extern DoubleOption opt_var_decay;
extern DoubleOption opt_clause_decay;
extern DoubleOption opt_theory_decay;
extern DoubleOption opt_random_var_freq;
extern DoubleOption opt_random_seed;
extern IntOption opt_ccmin_mode;
extern IntOption opt_phase_saving;
extern BoolOption opt_rnd_init_act;
extern BoolOption opt_luby_restart;
extern IntOption opt_restart_first;
extern DoubleOption opt_restart_inc;
extern DoubleOption opt_garbage_frac;
extern BoolOption opt_restarts;
extern BoolOption opt_rnd_restart;
extern BoolOption opt_rnd_theory_restart;
extern StringOption opt_record_file;
extern IntOption opt_limit_optimization_conflicts;
extern IntOption opt_limit_optimization_time;
extern BoolOption opt_limit_optimization_time_per_arg;
extern BoolOption opt_pb_theory;
extern bool opt_record;
extern DoubleOption opt_rnd_optimization_freq;
extern DoubleOption opt_rnd_optimization_restart_freq;
extern IntOption opt_theory_conflict_max;
extern DoubleOption opt_random_theory_freq;
extern BoolOption opt_theory_internal_vsids_fsm;
extern DoubleOption opt_random_theory_order_freq;
extern DoubleOption opt_randomize_theory_order_restart_freq;
extern BoolOption opt_theory_order_initial_sort;
extern BoolOption opt_randomize_theory_order_all;
extern BoolOption opt_theory_decision_round_robin;

extern BoolOption opt_optimization_init_solve;
extern BoolOption opt_decide_objectives_first;
extern BoolOption opt_strict_search_optimization;
extern BoolOption opt_interpolate;
extern IntOption opt_eager_prop;
extern IntOption opt_subsearch;

extern BoolOption opt_parser_immediate_mode;
extern BoolOption opt_remap_vars;
extern BoolOption opt_decide_optimization_lits;
extern IntOption opt_optimization_search_type;

extern IntOption opt_clausify_amo;
extern BoolOption opt_amo_eager_prop;

extern StringOption opt_debug_learnt_clauses;
extern BoolOption opt_debug_model;
extern int64_t opt_n_learnts;
extern BoolOption opt_write_bv_bounds;
extern BoolOption opt_write_bv_analysis;
extern FILE* opt_write_learnt_clauses;

//extern StringOption opt_fsm_model;


extern BoolOption opt_graph;
extern BoolOption opt_inc_graph;
extern IntOption opt_dec_graph;
extern BoolOption opt_conflict_shortest_path;
extern BoolOption opt_conflict_min_cut;

extern StringOption opt_maxflow_alg;
extern StringOption opt_reach_alg;
extern StringOption opt_con_alg;
extern StringOption opt_allpairs_alg;
extern StringOption opt_undir_allpairs_alg;
extern StringOption opt_dist_alg;
extern StringOption opt_mst_alg;

extern StringOption opt_components_alg;
extern StringOption opt_cycle_alg;
extern IntOption opt_adaptive_history_clear;
extern BoolOption opt_disable_history_clears;
extern IntOption opt_dynamic_history_clear;
extern BoolOption opt_lazy_backtrack;
extern BoolOption opt_lazy_backtrack_decisions;
extern IntOption opt_lazy_conflicts;
extern BoolOption opt_keep_lazy_conflicts;
extern BoolOption opt_lazy_backtrack_redecide;
extern BoolOption opt_theory_order_vsids;
extern BoolOption opt_theory_order_swapping;
extern BoolOption opt_theory_order_swapping_first_on_unit;
extern BoolOption opt_theory_order_conflict_sort_vsids;
extern BoolOption opt_theory_order_swapping_preserve_order;
extern BoolOption opt_theory_order_conflict_clear_on_restart;
extern BoolOption opt_theory_order_conflict_conservative;
extern BoolOption opt_theory_order_conflict_on_unit;
extern IntOption opt_theory_order_conflict_restart;
extern BoolOption opt_theory_order_conflict_sort_counter;
extern IntOption opt_theory_order_swapping_max_invovled;
extern BoolOption opt_theory_order_swapping_reset_counts_new_conflict;
extern BoolOption opt_theory_order_swapping_luby;
extern BoolOption opt_theory_order_conflict_luby;
extern IntOption opt_theory_order_conflict_count_preserve;
extern IntOption opt_theory_order_conflict_count_analysis;
extern BoolOption opt_theory_order_swapping_prioritize_last_decision;
extern BoolOption opt_decide_theories_only_prop_decision;
extern BoolOption opt_theory_order_conflict_skip_middle;

extern BoolOption opt_monolothic_theory_decisions;
extern BoolOption opt_vsids_both;
extern DoubleOption opt_theory_vsids_balance;
extern BoolOption opt_use_var_decay_for_theory_vsids;
extern BoolOption opt_vsids_solver_as_theory;
extern BoolOption opt_theory_internal_vsids;
extern BoolOption opt_theory_prioritize_conflicts;
extern BoolOption opt_theory_priority_clear;
extern BoolOption opt_theory_propagate_assumptions;

extern BoolOption opt_check_solution;
extern BoolOption opt_print_reach;
extern BoolOption opt_print_graph;
extern BoolOption opt_print_theory_debug;

extern IntOption opt_learn_reaches;
extern StringOption opt_priority;

extern BoolOption opt_print_conflicts;
extern BoolOption opt_rnd_phase;
extern BoolOption opt_init_rnd_phase;

extern BoolOption opt_reach_prop;
extern BoolOption opt_decide_theories;
extern BoolOption opt_decide_graph_distance;
extern BoolOption opt_decide_graph_bv;
extern BoolOption opt_cmp_lits_decidable;
extern BoolOption opt_decide_bv_bitwise;
extern BoolOption opt_decide_bv_intrinsic;
extern BoolOption opt_decide_theories_reverse;
extern BoolOption opt_use_random_path_for_decisions;
extern BoolOption opt_use_optimal_path_for_decisions;
extern DoubleOption opt_decide_graph_re_rnd;
extern BoolOption opt_print_decision_path;
extern BoolOption opt_force_distance_solver;
extern DoubleOption opt_allpairs_percentage;
extern BoolOption opt_decide_graph_neg;
extern BoolOption opt_decide_graph_pos;
extern BoolOption opt_ignore_theories;
extern BoolOption opt_check_pure_theory_lits;
extern BoolOption opt_mst_min_cut;
extern BoolOption opt_connected_components_min_cut;
extern BoolOption opt_optimize_mst;
extern BoolOption opt_skip_deletions;
extern BoolOption opt_skip_additions;
extern BoolOption opt_permanent_theory_conflicts;
extern IntOption opt_temporary_theory_conflicts;
extern IntOption opt_temporary_theory_reasons;
extern BoolOption opt_force_directed;
extern BoolOption opt_decide_graph_chokepoints;
extern IntOption opt_sort_graph_decisions;
extern IntOption opt_flow_router_heuristic;
extern IntOption opt_flow_router_policy;
extern BoolOption opt_ignore_assign_edge_weights;

extern BoolOption opt_decide_fsm_neg;
extern BoolOption opt_decide_fsm_pos;

extern BoolOption opt_compute_max_distance;
extern BoolOption opt_detect_pure_theory_lits;
extern BoolOption opt_detect_pure_lits;
extern IntOption opt_detect_satisfied_predicates;
extern BoolOption opt_propagate_theories_during_simplification;
extern BoolOption opt_propagate_theories_during_fast_simplification;
extern BoolOption opt_shrink_theory_conflicts;

extern IntOption opt_width;
extern IntOption opt_height;
extern IntOption opt_bits;
extern BoolOption opt_encode_reach_underapprox_as_sat;
extern IntOption opt_encode_dist_underapprox_as_sat;
extern BoolOption opt_sat_distance_encoding_unconstrained_default;
extern BoolOption opt_csv;
extern BoolOption opt_rnd_shuffle;
extern DoubleOption opt_rnd_shortest_path;
extern DoubleOption opt_rnd_shortest_edge;
extern BoolOption opt_components_learn_connect;
extern BoolOption opt_learn_unreachable_component;
extern BoolOption opt_dinics_recursive;


extern BoolOption opt_conflict_dfs;
extern BoolOption opt_conflict_from_source;
extern BoolOption opt_allow_reach_decisions;
extern BoolOption opt_conflict_1uip;
extern BoolOption opt_use_kt_for_conflicts;
//extern BoolOption opt_maxflow_backward;
extern BoolOption opt_conflict_min_cut_maxflow;
extern IntOption opt_history_clear;
extern BoolOption opt_kt_preserve_order;

extern IntOption opt_maxflow_decisions_type;

extern BoolOption opt_lazy_maxflow_decisions;
extern BoolOption opt_maxflow_allow_cycles;
extern BoolOption opt_old_lazy_maxflow_decisions;
extern IntOption opt_maxflow_decisions_q;
extern IntOption opt_maxflow_decision_paths;
extern BoolOption opt_reach_detector_combined_maxflow;
extern IntOption opt_adaptive_conflict_mincut;
extern BoolOption opt_shortest_path_prune_dist;
extern BoolOption opt_graph_bv_prop;

extern IntOption opt_graph_prop_skip;
extern IntOption opt_bv_prop_skip;
extern IntOption opt_fsm_prop_skip;

extern BoolOption opt_fsm_negate_underapprox;
extern BoolOption opt_fsm_edge_prop;
extern BoolOption opt_fsm_forced_edge_prop;
extern BoolOption opt_fsm_chokepoint_prop;

extern BoolOption opt_fsm_as_graph;
extern IntOption opt_fsm_symmetry_breaking;
extern BoolOption opt_fsm_track_used_transitions;

extern BoolOption opt_learn_acyclic_flows;


extern IntOption opt_min_edgeset;
extern BoolOption opt_only_prop_edgeset;

extern BoolOption opt_graph_cache_propagation;
extern IntOption opt_graph_use_cache_for_decisions;
extern OptionSet opt_route;
extern OptionSet opt_route2;


enum class ReachAlg {
    ALG_SAT,
    ALG_DFS,
    ALG_DIJKSTRA,
    ALG_DISTANCE,
    ALG_BFS,
    ALG_RAMAL_REPS,
    ALG_RAMAL_REPS_BATCHED,
    ALG_RAMAL_REPS_BATCHED2
};

//For undirected reachability
enum class ConnectivityAlg {
    ALG_SAT, ALG_DFS, ALG_DIJKSTRA, ALG_DISTANCE, ALG_BFS, ALG_THORUP
};
extern ConnectivityAlg undirectedalg;
extern ReachAlg reachalg;

enum class AllPairsAlg {
    ALG_FLOYDWARSHALL, ALG_DIJKSTRA_ALLPAIRS
};
extern AllPairsAlg allpairsalg;
enum class AllPairsConnectivityAlg {
    ALG_FLOYDWARSHALL, ALG_DIJKSTRA_ALLPAIRS, ALG_THORUP
};
extern AllPairsConnectivityAlg undirected_allpairsalg;
enum class MinCutAlg {
    ALG_EDMONSKARP, ALG_EDKARP_ADJ,
    // ALG_IBFS, //omitted for licensing reasons
            ALG_EDKARP_DYN,
    ALG_DINITZ,
    ALG_DINITZ_LINKCUT,
    ALG_KOHLI_TORR
};
extern MinCutAlg mincutalg;
enum class MinSpanAlg {
    ALG_KRUSKAL, ALG_PRIM, ALG_SPIRA_PAN
};
extern MinSpanAlg mstalg;
enum class ComponentsAlg {
    ALG_DISJOINT_SETS

};
extern ComponentsAlg componentsalg;

enum class CycleAlg {
    ALG_DFS_CYCLE,
    ALG_PK_CYCLE
};
extern CycleAlg cyclealg;


enum class ConvexHullAlg {
    ALG_MONOTONE_HULL, ALG_QUICKHULL

};
extern ConvexHullAlg hullAlg;

enum class DistAlg {
    ALG_SAT, ALG_DIJKSTRA, ALG_DISTANCE, ALG_RAMAL_REPS, ALG_RAMAL_REPS_BATCHED, ALG_RAMAL_REPS_BATCHED2
};

extern DistAlg distalg;

extern IntOption opt_time;

static inline double rtime(int level = 1){
    if(level <= opt_time){
        return fastTime();
    }else{
        return 0;
    }
}

}

#endif /* CONFIG_H_ */
