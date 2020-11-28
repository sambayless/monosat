/****************************************************************************************
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

#include "monosat/core/Config.h"

using namespace Monosat;

#ifdef DEBUG_CORE
int dbg_total_iterations = 0;
#endif

static const char* _cat = "CORE";
static const char* _cat_opt = "OPT";
static const char* _cat_theory = "THEORY";

static const char* _cat_sms = "SMS";
static const char* _cat_graph = "GRAPH";
static const char* _cat_bv = "BV";
static const char* _cat_amo = "AMO";
static const char* _cat_fsm = "FSM";


BoolOption Monosat::opt_show_version_and_quit("MAIN", "version", "Display version information and quit", false);
IntOption Monosat::opt_verb("MAIN", "verb", "Verbosity level (0=silent, 1=some, 2=more).", 0, IntRange(0, 5));
IntOption Monosat::opt_verb_optimize("MAIN", "verb-opt", "Verbosity level for optimization (0=silent, 1=some, 2=more).",
                                     0, IntRange(0, 5));

DoubleOption Monosat::opt_var_decay(_cat, "var-decay", "The variable activity decay factor", 0.95,
                                    DoubleRange(0, false, 1, false));
DoubleOption Monosat::opt_clause_decay(_cat, "cla-decay", "The clause activity decay factor", 0.999,
                                       DoubleRange(0, false, 1, false));
DoubleOption Monosat::opt_theory_decay(_cat, "theory-decay", "The theory activity decay factor", 0.95,
                                       DoubleRange(0, false, 1, false));
DoubleOption Monosat::opt_random_var_freq(_cat, "rnd-freq",
                                          "The frequency with which the decision heuristic tries to choose a random variable",
                                          0,
                                          DoubleRange(0, true, 1, true));
DoubleOption Monosat::opt_random_seed(_cat, "rnd-seed", "Used by the random variable selection", 91648253,
                                      DoubleRange(0, false, HUGE_VAL, false));
IntOption Monosat::opt_ccmin_mode(_cat, "ccmin-mode", "Controls conflict clause minimization (0=none, 1=basic, 2=deep)",
                                  2, IntRange(0, 2));
IntOption Monosat::opt_phase_saving(_cat, "phase-saving",
                                    "Controls the level of phase saving (0=none, 1=limited, 2=full)", 2,
                                    IntRange(0, 2));
BoolOption Monosat::opt_rnd_init_act(_cat, "rnd-init", "Randomize the initial activity", false);
BoolOption Monosat::opt_restarts(_cat, "restarts", "If false, then conflict-based restarts are disabled.", true);
BoolOption Monosat::opt_luby_restart(_cat, "luby", "Use the Luby restart sequence", true);
IntOption Monosat::opt_restart_first(_cat, "rfirst", "The base restart interval", 100, IntRange(1, INT32_MAX));
DoubleOption Monosat::opt_restart_inc(_cat, "rinc", "Restart interval increase factor", 2,
                                      DoubleRange(1, false, HUGE_VAL, false));
DoubleOption Monosat::opt_garbage_frac(_cat, "gc-frac",
                                       "The fraction of wasted memory allowed before a garbage collection is triggered",
                                       0.20,
                                       DoubleRange(0, false, HUGE_VAL, false));
BoolOption Monosat::opt_pre("MAIN", "pre", "Completely turn on/off any preprocessing.", true);
IntOption Monosat::opt_time(_cat, "verb-time", "Detail level of timing benchmarks (these add some overhead)", 0,
                            IntRange(0, 5));


BoolOption  Monosat::opt_pb_theory("PB", "pb-theory",
                                   "Use experimental lazy pseudo-boolean constraint solver (not recommended!) instead of minisat+",
                                   false);

IntOption Monosat::opt_limit_optimization_conflicts(_cat_opt, "opt-conflict-limit",
                                                    "Limit number of conflicts during optimization rounds before giving up (-1 to allow infinite conflicts)",
                                                    -1, IntRange(-1, INT32_MAX));
IntOption Monosat::opt_limit_optimization_time(_cat_opt, "opt-time-limit",
                                               "Limit time spent in optimization (after finding an initial satisfying solution) before giving up (-1 to allow infinite time during optimization)",
                                               -1, IntRange(-1, INT32_MAX));
BoolOption Monosat::opt_limit_optimization_time_per_arg(_cat_opt, "opt-time-limit-per-arg",
                                                        "Apply optimization time limit separately to each individual optimization argument, instead of globally to the whole optimization process",
                                                        false);


StringOption Monosat::opt_record_file(_cat, "debug-log",
                                      "Log (very expensive) debugging info at extensions of the following path (empty string (recommended) disables)",
                                      "");
bool Monosat::opt_record = false;

int64_t Monosat::opt_n_learnts = 0;
BoolOption Monosat::opt_debug_model(_cat, "debug-model", "", false);


StringOption Monosat::opt_debug_learnt_clauses(_cat, "debug-learnts",
                                               "Write all learned clauses to the following file (empty string (recommended) disables)", "");
FILE* Monosat::opt_write_learnt_clauses = nullptr;

//StringOption Monosat::StringOption opt_fsm_model(_cat_fsm,"File to write fsm model, if fsm theory is used","","");

BoolOption Monosat::opt_write_bv_analysis(_cat, "debug-analysis", "", false);
BoolOption Monosat::opt_write_bv_bounds(_cat, "debug-bounds", "", false);

IntOption Monosat::opt_theory_conflict_max(_cat, "theory-conflict-limit",
                                           "The maximum number of consecutive times the theory solver's can conflict before theory decisions are temporarily disabled (0 to all infinite theory conflicts)",
                                           0,
                                           IntRange(0, INT32_MAX));

DoubleOption Monosat::opt_random_theory_freq(_cat, "rnd-theory-freq",
                                             "The frequency with which the decision theory solvers are selected to make decisions",
                                             1,
                                             DoubleRange(0, true, 1, true));

DoubleOption Monosat::opt_random_theory_order_freq(_cat, "rnd-theory-rnd-freq",
                                                   "The frequency with which the decision theory uses vsids to make decisions, if theory-order-vsids is enabled",
                                                   1,
                                                   DoubleRange(0, true, 1, true));

DoubleOption Monosat::opt_randomize_theory_order_restart_freq(_cat_theory, "theory-order-rnd-restart-freq",
                                                              "Frequency with which to randomize the decision heuristic order (on restarts)",
                                                              0,
                                                              DoubleRange(0, true, 1, true));

BoolOption Monosat::opt_optimization_init_solve(_cat, "optimize-init-solve",
                                                "Perform an initial, free solving round before applying optimization (this is potentially expensive, but also gives the solver a chance to find good initial bounds)",
                                                true);
BoolOption Monosat::opt_decide_objectives_first(_cat, "decide-objectives-first",
                                                "Decide any free variables in objective functions first, when performing optimization",
                                                true);
BoolOption Monosat::opt_theory_order_initial_sort(_cat_theory, "theory-order-initial-sort", "", true);

BoolOption Monosat::opt_randomize_theory_order_all(_cat, "rnd-theory-order-all",
                                                   "If theory decisions are used, randomize the order that theories are decided at every decision (not just every restart)",
                                                   false);
BoolOption Monosat::opt_theory_decision_round_robin(_cat, "round-robin-theory-order",
                                                    "If theory decisions are used, use a round-robin policy for which theory to decide first",
                                                    false);

BoolOption Monosat::opt_parser_immediate_mode(_cat, "immediate-parser", "", false);
BoolOption Monosat::opt_remap_vars(_cat, "remap-vars",
                                   "Remap variables in the GNF internally in the solver, to minimize space required",
                                   true);
BoolOption Monosat::opt_decide_optimization_lits(_cat_opt, "decide-opt-lits",
                                                 "Allow decisions on literals introduced by optimization constraints",
                                                 true);
BoolOption  Monosat::opt_strict_search_optimization(_cat_opt, "strict-search",
                                                    "Use strict (rather than >=/<=) constraints during optimization (be warned, this is not always safe!)",
                                                    false);
IntOption  Monosat::opt_optimization_search_type(_cat_opt, "search",
                                                 "0: Use binary search for optimization. 1: Use linear search. 2: Use binary search, but switch to linear search for hard values.",
                                                 0, IntRange(0, 2));
DoubleOption Monosat::opt_rnd_optimization_freq(_cat_graph, "search-freq", "", 0.5, DoubleRange(0, true, 1, true));
DoubleOption Monosat::opt_rnd_optimization_restart_freq(_cat_graph, "search-restart-freq", "", 0.05,
                                                        DoubleRange(0, true, 1, true));

IntOption Monosat::opt_clausify_amo(_cat_amo, "clausify-amo",
                                    "Eagerly convert at-most-one predicates to clauses, if they have this many or fewer arguments (after removing level 0 constants). 0 to never clausify amo predicates.",
                                    30, IntRange(0, INT32_MAX));

BoolOption Monosat::opt_amo_eager_prop(_cat_amo, "amo-eager-prop",
                                       "Propagate a-m-o literals as soon as they are implied, instead of waiting for theory propagation",
                                       true);

BoolOption Monosat::opt_interpolate(_cat_sms, "interpolate",
                                    "Store learnt interface clauses to form interpolants between modules", false);
IntOption Monosat::opt_eager_prop(_cat_sms, "eager-prop",
                                  "Controls whether unit propagation is allowed to cross subsolver boundaries. 0= Disable. 1= Enable. 2=Enable, but don't cross the last interpolant. 3= Enable, but don't cross the last interpolant, or any earlier solver. 4= Enable, but dont cross the last interpolant, or any earlier solver, unless they have already had their interpolants strengthened",
                                  1, IntRange(0, 5));
IntOption Monosat::opt_subsearch(_cat_sms, "subsearch",
                                 "Control how the solver performs search on the subsolvers: 0=abort as soon as a conflict backtracks past the supersolvers decisionlevel. 1=Abort only once a conflict on the super-interface variables is found, allowing backtracks past those variables in the process. 2=Abort only when the the super-solvers assignment is proven to be in conflict. 3=Don't continue subsearach if the subsolver has backtracked past super-solver decisions. 4=Don't continue past the last interpolant level if any solver has backtracked past a super solver's decisions",
                                 2, IntRange(0, 4));
BoolOption Monosat::opt_permanent_theory_conflicts(_cat_sms, "permanent-theory-conflicts",
                                                   "True if conflict clauses from theory solvers are treated as permanent clauses; false if they are treated as learnt clauses (that is, allowed to be deleted by the solver)",
                                                   false);
IntOption Monosat::opt_temporary_theory_conflicts(_cat_sms, "temporary-theory-conflicts",
                                                  "True if theory conflict clauses larger than this size should be discarded immediately.",
                                                  INT32_MAX,
                                                  IntRange(0, INT32_MAX));
IntOption Monosat::opt_temporary_theory_reasons(_cat_sms, "temporary-theory-reasons",
                                                "True if theory reason clauses larger than this size should be discarded immediately.",
                                                INT32_MAX,
                                                IntRange(0, INT32_MAX));
BoolOption Monosat::opt_graph(_cat_graph, "graph", "Use graph theory solver", true);
BoolOption Monosat::opt_inc_graph(_cat_graph, "inc", "Use incremental graph reachability", false);
IntOption Monosat::opt_dec_graph(_cat_graph, "dec", "Use decremental graph reachability", 0, IntRange(0, 2));
StringOption Monosat::opt_maxflow_alg(_cat_graph, "maxflow",
                                      "Select max s-t-flow algorithm (edmondskarp, edmondskarp-adj, edmondskarp-dynamic,dinitz,dinitz-linkcut, kohli-torr)",
                                      "kohli-torr"); //ibfs
StringOption Monosat::opt_reach_alg(_cat_graph, "reach",
                                    "Select reachability algorithm (bfs,dfs, dijkstra,ramal-reps,cnf)", "ramal-reps");
StringOption Monosat::opt_dist_alg(_cat_graph, "dist",
                                   "Select reachability algorithm (bfs,dfs, dijkstra,ramal-reps,cnf)", "ramal-reps");

StringOption Monosat::opt_con_alg(_cat_graph, "connect",
                                  "Select undirected reachability algorithm (bfs,dfs, dijkstra, thorup,cnf)", "bfs");
StringOption Monosat::opt_undir_allpairs_alg(_cat_graph, "connect-allpairs",
                                             "Select allpairs reachability algorithm (floyd-warshall,dijkstra, thorup)",
                                             "floyd-warshall");
StringOption Monosat::opt_mst_alg(_cat_graph, "mst", "Select minimum spanning tree algorithm (kruskal,prim,spira-pan)",
                                  "spira-pan");

StringOption Monosat::opt_allpairs_alg(_cat_graph, "allpairs",
                                       "Select allpairs reachability algorithm (floyd-warshall,dijkstra)",
                                       "floyd-warshall");
StringOption Monosat::opt_components_alg(_cat_graph, "components",
                                         "Select connected-components algorithm (disjoint-sets, link-cut)",
                                         "disjoint-sets");
StringOption Monosat::opt_cycle_alg(_cat_graph, "cycles",
                                    "Select cycle detection algorithm (dfs, pk)", "pk");

BoolOption Monosat::opt_conflict_shortest_path(_cat_graph, "conflict-shortest-path",
                                               "Use shortest path (instead of arbitrary path) for conflict resolution (in theories that support this)",
                                               true);
BoolOption Monosat::opt_conflict_min_cut(_cat_graph, "conflict-min-cut",
                                         "Use min-cut (instead of arbitrary cut) for conflict resolution (in theories that support this)",
                                         false);
BoolOption Monosat::opt_conflict_1uip(_cat_graph, "conflict-1uip",
                                      "Use 1 uip (instead of arbitrary or min cut) for conflict resolution (in theories that support this)",
                                      false);
BoolOption Monosat::opt_conflict_min_cut_maxflow(_cat_graph, "conflict-min-cut-maxflow",
                                                 "Use min-cut (instead of arbitrary cut) for conflict resolution for maximum flow properties",
                                                 false);

IntOption Monosat::opt_history_clear(_cat_graph, "history-clear",
                                     "Rate at which the history of dynamic graphs is cleared", 1000,
                                     IntRange(1, INT32_MAX));
IntOption Monosat::opt_adaptive_history_clear(_cat_graph, "adaptive-history-clear",
                                              "If >0, ignore the history clear option, and instead set the history clear rate to be this value multiplied by the number of edges in the graph",
                                              0, IntRange(0, INT32_MAX));
BoolOption Monosat::opt_disable_history_clears(_cat_graph, "disable-history-clear", "", false);
IntOption Monosat::opt_dynamic_history_clear(_cat_graph, "dynamic-history-clear",
                                             "0=dont use dynamic history clears,1=use opportunistic dynamic history clears (falling back on normal history clears if that fails), 2=force dynamic history clears",
                                             0, IntRange(0, 2));

BoolOption Monosat::opt_lazy_backtrack(_cat_theory, "lazy-backtrack", "", false);
BoolOption Monosat::opt_lazy_backtrack_decisions(_cat_theory, "lazy-backtrack-decisions", "", false);
IntOption Monosat::opt_lazy_conflicts(_cat_theory, "lazy-conflicts",
                                      "0= unassign all lazy lits and reprop, 1=unassign all lazy lits in the clause, reprop, 2=unassign one lit, reprop, 3=skip lazy conflict analysis",
                                      0, IntRange(0, 3));
BoolOption Monosat::opt_keep_lazy_conflicts(_cat_theory, "keep-lazy-conflicts",
                                            "Keep clauses from lazy conflicts (only relevant if lazy-backtracking is enabled)",
                                            true);
BoolOption Monosat::opt_lazy_backtrack_redecide(_cat_theory, "lazy-backtrack-redecide", "", false);
BoolOption Monosat::opt_vsids_both(_cat_theory, "vsids-both",
                                   "Use vsids decision heuristic for both theories and main solver, in combination",
                                   false);
DoubleOption Monosat::opt_theory_vsids_balance(_cat_theory, "vsids-balance", "", 1,
                                               DoubleRange(0, false, 10000000, true));
BoolOption Monosat::opt_vsids_solver_as_theory(_cat_theory, "vsids-solver-as-theory",
                                               "Use vsids decision heuristic for both theories and main solver, treating the main solver as a theory",
                                               false);
BoolOption Monosat::opt_use_var_decay_for_theory_vsids(_cat_theory, "use-var-decay-for-theory-vsids-both", "", true);
BoolOption Monosat::opt_theory_order_vsids(_cat_theory, "theory-order-vsids",
                                           "Use vsids decision heuristic outside of theory solvers, to pick which theory solver to make decisions next",
                                           true);
BoolOption  Monosat::opt_theory_order_swapping(_cat_theory, "theory-order-swapping",
                                               "Use RUC-style theory order swapping decision heuristic to pick which theory solver to make decisions next",
                                               false);
BoolOption Monosat::opt_theory_order_swapping_preserve_order(_cat_theory, "theory-order-swapping-preserve", "", false);
IntOption  Monosat::opt_theory_order_conflict_restart(_cat_theory, "theory-order-conflict-restart",
                                                      "Use RUC-style theory order conflict counter-driven theory order restarts; restart the theory order after this many conflicts have occured in a theory (0 to disable)",
                                                      0, IntRange(0, INT32_MAX));
BoolOption Monosat::opt_theory_order_conflict_sort_counter(_cat_theory, "theory-order-restart-sorting",
                                                           "Sort decision heuristics by conflict counter on restarts",
                                                           false);
BoolOption Monosat::opt_theory_order_conflict_sort_vsids(_cat_theory, "theory-order-conflict-sort-vsids",
                                                         "Sort decision heuristics by vsids on conflict counter restarts, while putting the highest conflict counter first",
                                                         false);
BoolOption Monosat::opt_decide_theories_only_prop_decision(_cat_theory, "only-prop-theories-on-heuristic-change", "",
                                                           false);
BoolOption Monosat::opt_theory_order_swapping_first_on_unit(_cat_theory, "theory-order-swapping-first-on-unit",
                                                            "Move a decision heuristic to the begining of the heuristics if it is involved in a lone conflict",
                                                            false);
BoolOption Monosat::opt_theory_order_conflict_clear_on_restart(_cat_theory, "theory-order-conflict-restart-clear", "",
                                                               false);
BoolOption Monosat::opt_theory_order_conflict_conservative(_cat_theory, "theory-order-conflict-restart-conservative",
                                                           "", true);
BoolOption Monosat::opt_theory_order_conflict_on_unit(_cat_theory, "theory-order-conflict-on-unit", "", true);
BoolOption Monosat::opt_theory_order_swapping_reset_counts_new_conflict(_cat_theory,
                                                                        "theory-order-swapping-reset-counts-new-conflict",
                                                                        "", false);
BoolOption Monosat::opt_theory_order_conflict_skip_middle(_cat_theory, "theory-order-conflict-skip-middle", "", false);

IntOption Monosat::opt_theory_order_conflict_count_preserve(_cat_theory, "theory-order-conflict-count-preserve", "", 0,
                                                            IntRange(0, 4));

IntOption Monosat::opt_theory_order_conflict_count_analysis(_cat_theory, "theory-order-conflict-count-analysis", "", 0,
                                                            IntRange(0, 2));
BoolOption Monosat::opt_monolothic_theory_decisions(_cat_theory, "theory-combined-decisions",
                                                    "Make heuristic decisions within theories all at once, rather than individually per predicate.",
                                                    false);
BoolOption Monosat::opt_theory_order_swapping_prioritize_last_decision(_cat_theory,
                                                                       "theory-order-swapping-last-decision",
                                                                       "Prioritize the last decision heuristic when analyzing conflicts",
                                                                       false);
IntOption Monosat::opt_theory_order_swapping_max_invovled(_cat_theory, "theory-order-swapping-max",
                                                          "Maximum number of heuristics to treat as involved in a heuristic conflict.",
                                                          2, IntRange(2, INT32_MAX));
BoolOption Monosat::opt_theory_order_swapping_luby(_cat_theory, "theory-order-swapping-luby",
                                                   "Use Luby series to set the maximum number of heuristics involved in a heuristic conflict.",
                                                   false);
BoolOption Monosat::opt_theory_internal_vsids(_cat_graph, "theory-internal-vsids",
                                              "Use vsids decision heuristic within theory solvers", false);
BoolOption Monosat::opt_theory_internal_vsids_fsm(_cat_fsm, "theory-internal-vsids-fsm",
                                                  "Use vsids decision heuristic within fsm theory solvers", true);
BoolOption Monosat::opt_theory_prioritize_conflicts(_cat_graph, "theory-prioritize-conflicts", "", false);
BoolOption Monosat::opt_theory_priority_clear(_cat_graph, "theory-prioritize-clear", "", false);
BoolOption Monosat::opt_theory_order_conflict_luby(_cat_theory, "theory-order-conflict-luby",
                                                   "Use Luby series to set the maximum number of conflicts in any one heuristic before performing a restart.",
                                                   false);


BoolOption Monosat::opt_check_solution(_cat_graph, "check-solution", "Double check solution", true);
BoolOption Monosat::opt_print_reach(_cat_graph, "print-reach", "Print reachability graphs", false);
BoolOption Monosat::opt_print_graph(_cat_graph, "print-graph", "Print digraph", false);
BoolOption Monosat::opt_print_theory_debug(_cat_graph, "print-theory-debug",
                                           "Print internal theory solution assignments (intended for debugging only; output format is not stable)",
                                           false);

BoolOption Monosat::opt_compute_max_distance(_cat_graph, "max-distance",
                                             "Optimize shortest path computation to abort when a path is longer than the longest path in the constraints.",
                                             true);
BoolOption Monosat::opt_learn_unreachable_component(_cat_graph, "learn-component", "", false);
BoolOption Monosat::opt_force_directed(_cat_graph, "force-directed",
                                       "Use directed reachability algorithms for undirected reachability queries (by duplicating edges as needed)\n",
                                       false);

BoolOption Monosat::opt_rnd_restart(_cat, "rnd-restart", "Randomize activity on restart", false);
BoolOption Monosat::opt_rnd_theory_restart(_cat, "rnd-theory-restart", "Randomize theory activity on restart", false);

IntOption Monosat::opt_learn_reaches(_cat_graph, "learn-reach",
                                     "Learn using reach variables: 0 = Never, 1=Paths, 2=Cuts,3=Always", 0,
                                     IntRange(0, 3));

StringOption Monosat::opt_priority(_cat_graph, "priority", "Decision variable priority list", "");

BoolOption Monosat::opt_print_conflicts(_cat, "print-conflicts", "", false);

BoolOption Monosat::opt_rnd_phase(_cat, "rnd-phase", "", false);
BoolOption Monosat::opt_init_rnd_phase(_cat, "init-rnd-phase", "", false);

BoolOption Monosat::opt_encode_reach_underapprox_as_sat(_cat_graph, "reach-underapprox-cnf",
                                                        "Compute the under-approximate side of reachability constraints using CNF (only requires linear number of constraints), instead of the chosen algorithm",
                                                        true);

IntOption Monosat::opt_encode_dist_underapprox_as_sat(_cat_graph, "dist-underapprox-cnf",
                                                      "Compute the under-approximate side of distance constraints using CNF, instead of the chosen algorithm (0=don't use CNF encoding)",
                                                      0, IntRange(0, 2));
BoolOption Monosat::opt_sat_distance_encoding_unconstrained_default(_cat_graph,
                                                                    "dist-underapprox-cnf-dst-unconstrained", "", true);
BoolOption Monosat::opt_reach_prop(_cat_graph, "prop-reach", "", false);

BoolOption Monosat::opt_decide_theories(_cat_graph, "decide-theories", "", false);
BoolOption Monosat::opt_decide_graph_distance(_cat_graph, "decide-graph-dist", "", false);
BoolOption Monosat::opt_decide_graph_bv(_cat_graph, "decide-graph-bv", "", false);
BoolOption Monosat::opt_cmp_lits_decidable(_cat_graph, "decide-cmp-lits",
                                           "Controls whether or not comparison lits introduced by the bv solver (but not in the original formula) can be chosen as decisions by the SAT solver",
                                           false);
BoolOption Monosat::opt_decide_bv_intrinsic(_cat_graph, "decide-bv-intrinsic", "", false);
BoolOption Monosat::opt_decide_bv_bitwise(_cat_graph, "decide-bv-bitwise", "", false);
BoolOption Monosat::opt_decide_theories_reverse(_cat_graph, "decide-theories-reverse",
                                                "Decide theories in reverse order, if theory decisions are enabled",
                                                false);
BoolOption Monosat::opt_use_random_path_for_decisions(_cat_graph, "decide-graph-rnd", "", false);
BoolOption Monosat::opt_use_optimal_path_for_decisions(_cat_graph, "decide-graph-opt",
                                                       "When selecting a path during decisions, find the shortest path, excluding the weight of already assigned edges.",
                                                       false);

DoubleOption Monosat::opt_decide_graph_re_rnd(_cat_graph, "decide-graph-re-rnd",
                                              "Randomly make new random graphs for graph decisions instead of sticking with just one",
                                              0.01);

BoolOption Monosat::opt_print_decision_path(_cat_graph, "decide-graph-print", "", false);
BoolOption Monosat::opt_force_distance_solver(_cat_graph, "force-distance",
                                              "Force the graph distance solver to be used, instead of the optimized reachability solver",
                                              false);

DoubleOption Monosat::opt_allpairs_percentage(_cat_graph, "allpairs-frac",
                                              "Fraction of nodes as source reach querries at which to trigger using allpairs solver instead of separate theory solvers. 0 to force all querries to go through allpairs solver. 1 to disable.",
                                              1, DoubleRange(0, true, 1, true));

BoolOption Monosat::opt_decide_graph_neg(_cat_graph, "decide-graph-neg", "", false);

BoolOption Monosat::opt_decide_graph_pos(_cat_graph, "decide-graph-pos", "", true);

BoolOption Monosat::opt_ignore_theories(_cat_graph, "ignore-theories", "", false);

//BoolOption Minisat::opt_check_pure_theory_lits(_cat_graph,"pure-theory-lits","",false);

BoolOption Monosat::opt_decide_graph_chokepoints(_cat_graph, "decide-graph-chokepoints", "", false);
IntOption Monosat::opt_sort_graph_decisions(_cat_graph, "decide-graph-sort",
                                            "0=dont sort, 1=sort by shortest, 2=sort by longest", 0, IntRange(0, 2));

BoolOption Monosat::opt_decide_fsm_neg(_cat_fsm, "decide-fsm-neg",
                                       "", false);

BoolOption Monosat::opt_decide_fsm_pos(_cat_fsm, "decide-fsm-pos",
                                       "", true);
BoolOption Monosat::opt_detect_pure_lits(_cat, "detect-pure-lits",
                                         "Detect pure literals in the main solver (only during simplification, not during search)",
                                         true);

IntOption Monosat::opt_detect_satisfied_predicates(_cat, "detect-sat-predicates",
                                                   "Take extra effort to detect satisfied"
                                                   " theory atoms, so that they don't need to be checked subsequently. Check every Nth round (0 to never check, 1 to always check)",
                                                   0, IntRange(0, INT32_MAX));

IntOption Monosat::opt_flow_router_heuristic(_cat, "flow-router-heuristic",
                                             "Enabled if >0. If opt_flow_router_heuristic=2, then suppress flow decisions when only 1 path is left in the flow.",
                                             0, IntRange(0, 5));
IntOption Monosat::opt_flow_router_policy(_cat, "flow-router-policy", "Flow router propagation enabled if >0.", 1,
                                          IntRange(0, 2));

BoolOption Monosat::opt_ignore_assign_edge_weights(_cat_graph, "ignore-edge-weight-heuristic",
                                                   "Ignore any edge weight heuristic directives", false);


BoolOption Monosat::opt_detect_pure_theory_lits(_cat, "detect-pure-theory-lits",
                                                "Detect pure literals in the theory solvers", true);

BoolOption Monosat::opt_mst_min_cut(_cat_graph, "mst-min-cut",
                                    "Search for a min-cut during conflict resolution of disconnected minimum spanning trees",
                                    true);
BoolOption Monosat::opt_connected_components_min_cut(_cat_graph, "cc-mincut",
                                                     "Search for a min-cut during conflict resolution of connected components",
                                                     true);
BoolOption Monosat::opt_optimize_mst(_cat_graph, "opt-mst",
                                     "Find the solution that minimizes the spanning tree by making repeated sat checks.",
                                     false);
BoolOption Monosat::opt_skip_deletions(_cat_graph, "skip-deletions", "", false);
BoolOption Monosat::opt_skip_additions(_cat_graph, "skip-additions", "", false);

BoolOption Monosat::opt_allow_reach_decisions(_cat_graph, "allow-reach-decision", "", true);

BoolOption Monosat::opt_conflict_dfs(_cat_graph, "conflict-dfs",
                                     "Use a DFS (instead of a BFS) to find the conflict cut", true);
BoolOption Monosat::opt_conflict_from_source(_cat_graph, "conflict-from-source",
                                             "Search from the source (instead of back from the target) to find the conflict cut",
                                             false);
BoolOption Monosat::opt_use_kt_for_conflicts(_cat_graph, "use-kt-for-conflicts",
                                             "When kohli-torr is selected as the maximum flow algorithm, also use kohli-torr for conflcits and decisions (instead of edmonds karp)",
                                             true);
//BoolOption Monosat::opt_maxflow_backward(_cat_graph,"maxflow-backward","Reverse source,sink, and all edge capacities when computing maxflow (this can force different flow assignments by individual algorithms, but doesn't effect the set of valid solutions)",false);
BoolOption Monosat::opt_kt_preserve_order(_cat_graph, "kt-preserve-order",
                                          "Attempt to preserve the order of flow assigned by the kohli-torr maxflow algorithm",
                                          false);


IntOption Monosat::opt_maxflow_decisions_type(_cat_graph, "maxflow-decisions", "0=None, 1=flow-based, 2=cut-based", 1,
                                              IntRange(0, 2));

BoolOption Monosat::opt_lazy_maxflow_decisions(_cat_graph, "lazy-maxflow-decisions", "", true);
BoolOption Monosat::opt_maxflow_allow_cycles(_cat_graph, "allow-maxflow-cycles",
                                             "Allow (superfluous) cycles in the maxflow solution", false);

BoolOption Monosat::opt_old_lazy_maxflow_decisions(_cat_graph, "old-lazy-maxflow-decisions", "", false);
IntOption Monosat::opt_maxflow_decisions_q(_cat_graph, "maxflow-decisions-q",
                                           "Use a FIFO (instead of LIFO) decision order in the maxflow theory", 1,
                                           IntRange(0, 4));
IntOption Monosat::opt_maxflow_decision_paths(_cat_graph, "maxflow-decision-paths",
                                              "Make maxflow decisions path by path (0 = never, 1 = forward, 2= backward)",
                                              0, IntRange(0, 3));

BoolOption Monosat::opt_reach_detector_combined_maxflow(_cat_graph, "reach-combined-maxflow", "", true);
IntOption Monosat::opt_adaptive_conflict_mincut(_cat_graph, "adaptive-conflict-mincut",
                                                "First try applying conflict detection without mincut analysis (which is faster), then try again with mincut analysis if the learnt clause is >= this length (0 to disable, 1 to always use mincut analysis)",
                                                0, IntRange(0, INT32_MAX));

BoolOption Monosat::opt_shortest_path_prune_dist(_cat_graph, "shortest-paths-prune-dist",
                                                 "Prune edges based on distances from learnt clauses for the shortest paths theory",
                                                 false);

ConvexHullAlg Monosat::hullAlg = ConvexHullAlg::ALG_MONOTONE_HULL;

BoolOption Monosat::opt_propagate_theories_during_simplification(_cat, "theory-prop-during-simp",
                                                                 "Apply propagation to theory solvers during simplification. Can be very expensive (depending on the theory).",
                                                                 true);


BoolOption Monosat::opt_propagate_theories_during_fast_simplification(_cat, "theory-prop-during-fast-simp",
                                                                      "Apply propagation to theory solvers during (fast) simplification. Can be very expensive (depending on the theory).",
                                                                      true);


BoolOption Monosat::opt_theory_propagate_assumptions(_cat, "theory-prop-during-assumps",
                                                     "Apply theory propagation while assigning assumptions (can be slow, if there are many assumptions)",
                                                     false);

BoolOption Monosat::opt_shrink_theory_conflicts(_cat, "shrink-theory-conflicts", "", false);

BoolOption Monosat::opt_rnd_shuffle(_cat_graph, "rnd-shuffle",
                                    "Inject randomness into the solver by shuffling the order of propagation of graph constraints.",
                                    true);

DoubleOption Monosat::opt_rnd_shortest_path(_cat_graph, "rnd-shortest-path-freq",
                                            "Frequency with which to select a random (but still shortest) path in path selection heuristics.",
                                            0, DoubleRange(0, true, 1, true));

DoubleOption Monosat::opt_rnd_shortest_edge(_cat_graph, "rnd-shortest-edge-freq",
                                            "Frequency with which to select a random (but still shortest) edge in path selection heuristics.",
                                            0, DoubleRange(0, true, 1, true));


BoolOption Monosat::opt_components_learn_connect(_cat_graph, "components-learn-connect", "", false);

BoolOption Monosat::opt_dinics_recursive(_cat_graph, "dinitz-recursive",
                                         "Use the recursive (default: iterative) Dinic's Maximum-flow implementation",
                                         false);

IntOption  Monosat::opt_graph_prop_skip(_cat_graph, "graph-theory-skip",
                                        "Only process every nth graph theory propagation ('1' skips no propagations)",
                                        1, IntRange(1, INT32_MAX));

IntOption  Monosat::opt_bv_prop_skip(_cat_bv, "bv-theory-skip",
                                     "Only process every nth bv theory propagation ('1' skips no propagations)", 1,
                                     IntRange(1, INT32_MAX));

IntOption  Monosat::opt_fsm_prop_skip(_cat_fsm, "fsm-theory-skip",
                                      "Only process every nth fsm theory propagation ('1' skips no propagations)", 1,
                                      IntRange(1, INT32_MAX));

BoolOption  Monosat::opt_graph_bv_prop(_cat_graph, "graph-bv-prop", "", false);

BoolOption Monosat::opt_fsm_track_used_transitions(_cat_fsm, "fsm-track-used", "", true);
IntOption Monosat::opt_fsm_symmetry_breaking(_cat_fsm, "fsm-symmetry-breaking",
                                             "", 0, IntRange(0, 2));
BoolOption Monosat::opt_fsm_negate_underapprox(_cat_fsm, "fsm-negate-under",
                                               "", false);//this option appears to be broken
BoolOption Monosat::opt_fsm_forced_edge_prop(_cat_fsm, "fsm-forced-edge-prop", "", false);
BoolOption Monosat::opt_fsm_edge_prop(_cat_fsm, "fsm-edge-prop", "", false);
BoolOption Monosat::opt_fsm_as_graph(_cat_fsm, "fsm-as-graph", "", true);
BoolOption Monosat::opt_fsm_chokepoint_prop(_cat_fsm, "fsm-chokepoint", "", false);

BoolOption Monosat::opt_learn_acyclic_flows(_cat_graph, "learn-acyclic-flows", "", false);

IntOption Monosat::opt_min_edgeset(_cat_graph, "min-edge-set",
                                   "set to 0 to accept all edge sets, or -1 to disable all edge sets", 10,
                                   IntRange(-1, INT32_MAX));
BoolOption Monosat::opt_only_prop_edgeset(_cat_graph, "only-prop-edgesets", "", false);


IntOption Monosat::opt_width("GRAPH", "width", "Width of graph.\n", 0, IntRange(0, INT32_MAX));
IntOption Monosat::opt_height("GRAPH", "height", "Height of graph.\n", 0, IntRange(0, INT32_MAX));
IntOption Monosat::opt_bits("GRAPH", "bits", "Bits per position in graph.\n", 1, IntRange(0, INT32_MAX));

BoolOption Monosat::opt_graph_cache_propagation(_cat_graph, "cache-graph-propagation", "", true);
IntOption Monosat::opt_graph_use_cache_for_decisions(_cat_graph, "cache-graph-decisions",
                                                     "(when using cache-graph-propagation, 2== use cache for all decisions, 1 = clear cache before decisions, 0=never use cache for decisions",
                                                     2, IntRange(0, 2));

OptionSet Monosat::opt_route(_cat_graph, "ruc", "Use RUC-syle heuristics for routing circuits",
                             "-cache-graph-propagation -cache-graph-decisions=2 -only-prop-theories-on-heuristic-change -theory-order-swapping-preserve -no-theory-order-conflict-sort-vsids -no-theory-order-initial-sort -no-theory-order-swapping-first-on-unit -theory-order-conflict-restart=10 -theory-order-swapping-luby -no-theory-order-restart-sorting -no-theory-order-swapping-last-decision -decide-theories -theory-order-swapping -no-theory-order-vsids",
                             false);
OptionSet Monosat::opt_route2(_cat_graph, "route",
                              "Use RUC-syle heuristics for routing circuits, with alternate configuration options",
                              "-theory-order-conflict-count-analysis=0 -theory-order-conflict-luby -theory-order-conflict-on-unit -theory-order-conflict-restart=50 -theory-order-conflict-restart-clear -no-theory-order-conflict-restart-conservative -no-theory-order-conflict-skip-middle -no-theory-order-conflict-sort-vsids -theory-order-initial-sort  -no-theory-order-restart-sorting -theory-order-swapping -no-theory-order-swapping-first-on-unit -no-theory-order-swapping-luby -no-theory-order-swapping-reset-counts-new-conflict  -no-theory-order-swapping-preserve  -theory-order-swapping-max=3  -no-theory-order-swapping-last-decision   -theory-order-rnd-restart-freq=0 -theory-order-conflict-count-preserve=1  -no-restarts -only-prop-theories-on-heuristic-change -decide-theories -cache-graph-propagation -cache-graph-decisions=1 -no-theory-order-vsids -rnd-shortest-path-freq=0.05 ",
                              false);

BoolOption Monosat::opt_csv("GRAPH", "csv", "Output in CSV format", false);

MinCutAlg Monosat::mincutalg = MinCutAlg::ALG_EDMONSKARP;
ReachAlg Monosat::reachalg = ReachAlg::ALG_RAMAL_REPS_BATCHED;
ConnectivityAlg Monosat::undirectedalg = ConnectivityAlg::ALG_DFS;
DistAlg Monosat::distalg = DistAlg::ALG_RAMAL_REPS_BATCHED;
AllPairsAlg Monosat::allpairsalg = AllPairsAlg::ALG_DIJKSTRA_ALLPAIRS;
AllPairsConnectivityAlg Monosat::undirected_allpairsalg = AllPairsConnectivityAlg::ALG_DIJKSTRA_ALLPAIRS;
ComponentsAlg Monosat::componentsalg = ComponentsAlg::ALG_DISJOINT_SETS;
MinSpanAlg Monosat::mstalg = MinSpanAlg::ALG_KRUSKAL;
CycleAlg Monosat::cyclealg = CycleAlg::ALG_PK_CYCLE;

