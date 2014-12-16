/****************************************************************************************[Solver.h]
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
#ifndef FSM_REACHDETECTOR_H_
#define FSM_REACHDETECTOR_H_
#include "utils/System.h"

#include "dgl/DynamicGraph.h"

#include "DynamicFSM.h"

#include "core/SolverTypes.h"
#include "mtl/Map.h"


#include "utils/System.h"
#include "FSMDetector.h"

using namespace dgl;
namespace Monosat {

class FSMTheorySolver;

class FSMReachDetector: public FSMDetector {
public:
	FSMTheorySolver* outer;
	DynamicFSM &f_under;
	DynamicFSM &f_over;

	int within;
	int source;
	double rnd_seed;


	CRef underprop_marker;
	CRef overprop_marker;
	CRef forced_edge_marker;


	struct Change {
		Var v;
		int u;
	};
	vec<bool> is_changed;
	vec<Change> changed;

s
	//stats
	
	int stats_full_updates = 0;
	int stats_fast_updates = 0;
	int stats_fast_failed_updates = 0;
	int stats_skip_deletes = 0;
	int stats_skipped_updates = 0;
	int stats_num_skipable_deletions = 0;
	int stats_learnt_components = 0;
	int stats_learnt_components_sz = 0;
	double mod_percentage = 0.2;
	int stats_pure_skipped = 0;
	int stats_shrink_removed = 0;
	double stats_full_update_time = 0;
	double stats_fast_update_time = 0;

	void printStats() {
		//printf("Reach detector\n");
		FSMDetector::printStats();
		if (opt_detect_pure_theory_lits)
			printf("\tPropagations skipped by pure literal detection: %d\n", stats_pure_skipped);
		if (opt_shrink_theory_conflicts) {
			printf("\t%d lits removed by shrinking conflicts\n", stats_shrink_removed);
		}
		if (opt_learn_unreachable_component) {
			printf("\t%d components learned, average component size: %f\n", stats_learnt_components,
					stats_learnt_components_sz / (float) stats_learnt_components);
		}
	}
	

	void unassign(Lit l) {
		FSMDetector::unassign(l);
		int index = var(l) - first_reach_var;
		if (index >= 0 && index < reach_lit_map.size() && reach_lit_map[index] != -1) {
			int node = reach_lit_map[index];
			if (!is_changed[node]) {
				changed.push( { var(l), node });
				is_changed[node] = true;
			}
		}
	}
	
	int getNode(Var reachVar) {
		assert(reachVar >= first_reach_var);
		int index = reachVar - first_reach_var;
		assert(index < reach_lit_map.size());
		assert(reach_lit_map[index] >= 0);
		return reach_lit_map[index];
	}
	void backtrack(int level) {
		to_decide.clear();
		last_decision_status = -1;
		
	}


	void buildSATConstraints(bool onlyUnderApprox = false, int within_steps = -1);
	bool propagate(vec<Lit> & conflict);
	void buildReachReason(int node, vec<Lit> & conflict);
	void buildNonReachReason(int node, vec<Lit> & conflict, bool force_maxflow = false);
	void buildForcedEdgeReason(int reach_node, int forced_edge_id, vec<Lit> & conflict);
	void buildReason(Lit p, vec<Lit> & reason, CRef marker);
	bool checkSatisfied();
	void printSolution(std::ostream& write_to);

	void addLit(int from, int to, Var reach_var);
	Lit decide(int level);
	void preprocess();
	void dbg_sync_reachability();

	FSMReachDetector(int _detectorID, FSMTheorySolver * _outer, DynamicGraph &_g, DynamicGraph &_antig,
			int _source, double seed = 1);
	virtual ~FSMReachDetector() {
		
	}
	
	const char* getName() {
		return "Reachability Detector";
	}
	
	bool dbg_cut(std::vector<MaxFlowEdge> & cut, DynamicGraph & graph, int source, int node) {
#ifndef NDEBUG
		
		DynamicGraph t;
		for (int i = 0; i < graph.nodes(); i++)
			t.addNode();
		std::vector<int> capacity;
		for (int id = 0; id < graph.edges(); id++) {
			t.addEdge(graph.getEdge(id).from, graph.getEdge(id).to, id);
			if (id % 2 == 0) {
				bool incut = false;
				for (int i = 0; i < cut.size(); i++) {
					if (cut[i].id == id) {
						incut = true;
						break;
					}
				}
				if (incut) {
					capacity.push_back(0);
					t.disableEdge(id);
				} else
					capacity.push_back(1);
			} else {
				capacity.push_back(0xFFFF);
			}
			
		}
		EdmondsKarpAdj<std::vector<int>, int> check(t, capacity, source, node);
		std::vector<MaxFlowEdge> check_cut;
		int flow = check.minCut(check_cut);
		assert(flow < 0xFFFF);
		if (flow < 0xFFFF) {
			exit(4);
		}

#endif
		return true;
	}
	
};
}
;
#endif /* REACHDETECTOR_H_ */
