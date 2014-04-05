/*
 * Config.h
 *
 *  Created on: 2012-10-20
 *      Author: sam
 */

#ifndef CONFIG_H_
#define CONFIG_H_

#include "utils/Options.h"
namespace Minisat {
extern IntOption opt_verb;
extern  DoubleOption  opt_var_decay;
extern  DoubleOption  opt_clause_decay ;
extern  DoubleOption  opt_random_var_freq ;
extern  DoubleOption  opt_random_seed;
extern  IntOption     opt_ccmin_mode  ;
extern  IntOption     opt_phase_saving ;
extern  BoolOption    opt_rnd_init_act ;
extern  BoolOption    opt_luby_restart ;
extern  IntOption     opt_restart_first ;
extern  DoubleOption  opt_restart_inc   ;
extern  DoubleOption  opt_garbage_frac  ;
extern  BoolOption    opt_rnd_restart ;

extern  BoolOption opt_interpolate;
extern IntOption opt_eager_prop;
extern IntOption opt_subsearch;

extern BoolOption opt_graph;
extern BoolOption opt_inc_graph;
extern IntOption opt_dec_graph;
extern BoolOption opt_conflict_shortest_path;
extern BoolOption opt_conflict_min_cut;
extern IntOption opt_restrict_decisions;

extern StringOption opt_min_cut_alg;
extern StringOption opt_reach_alg;
extern StringOption opt_con_alg;
extern StringOption opt_allpairs_alg;
extern StringOption opt_undir_allpairs_alg;

extern StringOption opt_components_alg;
extern BoolOption opt_check_solution;
extern BoolOption opt_print_reach;

extern IntOption opt_learn_reaches;
extern StringOption opt_priority;

extern BoolOption opt_print_conflicts;
extern BoolOption opt_rnd_phase;
extern BoolOption opt_init_rnd_phase;

extern BoolOption opt_reach_prop;
extern BoolOption opt_decide_graph;
extern BoolOption opt_decide_graph_distance;
extern BoolOption opt_use_random_path_for_decisions;
extern BoolOption opt_use_optimal_path_for_decisions;
extern DoubleOption opt_decide_graph_re_rnd;
extern BoolOption opt_print_decision_path;
extern BoolOption opt_force_distance_solver;
extern DoubleOption opt_allpairs_percentage;
extern BoolOption opt_decide_graph_neg;
extern BoolOption opt_decide_graph_pos;
extern BoolOption opt_ignore_graph;
extern BoolOption opt_check_pure_theory_lits;
extern BoolOption opt_mst_min_cut;
extern BoolOption opt_connected_components_min_cut;
extern BoolOption opt_optimize_mst;
extern BoolOption opt_skip_deletions;
extern BoolOption opt_skip_additions;
extern BoolOption opt_permanent_theory_conflicts;
enum class ReachAlg{
	 ALG_DFS,
	 ALG_DIJKSTRA,
	 ALG_DISTANCE,
	 ALG_BFS
};

//For undirected reachability
enum class ConnectivityAlg{
	 ALG_DFS,
	 ALG_DIJKSTRA,
	 ALG_DISTANCE,
	 ALG_BFS,
	 ALG_THORUP
};
extern ConnectivityAlg undirectedalg;
extern ReachAlg reachalg;
extern ReachAlg distalg;
enum class AllPairsAlg{
	 ALG_FLOYDWARSHALL,
	 ALG_DIJKSTRA_ALLPAIRS
};
extern AllPairsAlg allpairsalg;
enum class AllPairsConnectivityAlg{
	 ALG_FLOYDWARSHALL,
	 ALG_DIJKSTRA_ALLPAIRS,
	 ALG_THORUP
};
extern AllPairsConnectivityAlg undirected_allpairsalg;
enum class MinCutAlg{
	 ALG_EDMONSKARP,
	 ALG_EDKARP_ADJ,
	 ALG_IBFS
};
extern MinCutAlg mincutalg;
enum class ComponentsAlg{
	ALG_DISJOINT_SETS

};
extern ComponentsAlg componentsalg;
}

#endif /* CONFIG_H_ */
