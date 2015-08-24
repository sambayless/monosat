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
#ifndef REACHDETECTOR_H_
#define REACHDETECTOR_H_
#include "utils/System.h"
#include "GraphTheoryTypes.h"
#include "dgl/DynamicGraph.h"
#include "dgl/Reach.h"
#include "dgl/Dijkstra.h"
#include "dgl/BFS.h"
#include "dgl/DFS.h"

#include "core/SolverTypes.h"
#include "mtl/Map.h"
#include "dgl/MaxFlow.h"

#include "dgl/EdmondsKarp.h"
#include "dgl/EdmondsKarpAdj.h"
#include "dgl/Chokepoint.h"
#include "WeightedDijkstra.h"

#include "utils/System.h"
#include "Detector.h"
using namespace dgl;
namespace Monosat {
template<typename Weight>
class GraphTheorySolver;
template<typename Weight>
class ReachDetector: public Detector {
public:
	GraphTheorySolver<Weight> * outer;
	DynamicGraph<Weight> &g_under;
	DynamicGraph<Weight> &g_over;
	DynamicGraph<Weight> cutgraph;
	int within;
	int source;
	double rnd_seed;
	int constraintsBuiltOver = -1;
	int constraintsBuiltUnder = -1;

	CRef underprop_marker;
	CRef overprop_marker;
	CRef forced_edge_marker;

	Reach * underapprox_detector = nullptr;
	Reach * overapprox_reach_detector = nullptr;
	Reach * underapprox_path_detector = nullptr;
	Reach * overapprox_path_detector = nullptr;
	Reach * cutgraph_detector = nullptr;
	Reach * underapprox_fast_detector = nullptr;
	Distance<int> * negative_distance_detector = nullptr;
	vec<bool> original_reach_lits;
	vec<Lit> reach_lits;
	vec<Lit> cnf_reach_lits;
	Var first_reach_var;
	vec<int> order_vec;
	vec<int> reach_lit_map;
	vec<int> force_reason;

	vec<Lit> to_decide;
	int last_decision_status = -1;
	/*
	 struct DistLit{
	 Lit l;
	 int min_distance;

	 };*/
	vec<vec<Lit> > dist_lits;

	std::vector<ForceReason> forced_edges;
	struct Change {
		Var v;
		int u;
	};
	vec<bool> is_changed;
	vec<Change> changed;

	vec<Lit> extra_conflict;
	vec<int> removed_edges;
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
		Detector::printStats();
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
	
	struct ReachStatus {
		ReachDetector & detector;
		bool polarity;
		void setReachable(int u, bool reachable);
		bool isReachable(int u) const {
			return false;
		}
		
		void setMininumDistance(int u, bool reachable, Weight distance);

		ReachStatus(ReachDetector & _outer, bool _polarity) :
				detector(_outer), polarity(_polarity) {
		}
	};
	ReachStatus *positiveReachStatus = nullptr;
	ReachStatus *negativeReachStatus = nullptr;
	MaxFlow<long> * conflict_flow = nullptr;
	std::vector<MaxFlow<long> *> conflict_flows;

	WeightedDijkstra<Weight,double> * rnd_path = nullptr;
	std::vector<double> rnd_weight;
	/*struct OptimalWeightEdgeStatus{
	 ReachDetector & detector;
	 int operator [] (int edge) const ;
	 int size()const;
	 OptimalWeightEdgeStatus(ReachDetector & _outer):detector(_outer){}

	 };

	 OptimalWeightEdgeStatus opt_weight;
	 WeightedDijkstra<OptimalWeightEdgeStatus> * opt_path;*/
	Reach * chokepoint_detector = nullptr;
/*	struct CutStatus {
		long one = 1;
		long inf = 0xFFFF;
		ReachDetector & outer;

		const long &operator [](int id) const {
			if (id % 2 == 0) {
				return one;
			} else {
				return inf;
			}
		}
		int size() const {
			return outer.g_under.edges() * 2;
		}
		CutStatus(ReachDetector & _outer) :
				outer(_outer) {
		}
		
	} cutStatus;*/
	std::vector<MaxFlowEdge> cut;
	/*		struct ChokepointStatus{
	 ReachDetector & detector;
	 bool mustReach(int node);
	 bool operator() (int edge_id);
	 ChokepointStatus(ReachDetector & _outer):detector(_outer){

	 }
	 }chokepoint_status;

	 Chokepoint<ChokepointStatus > chokepoint;*/

	//This class can stand in as a reach algorithm if we have encoded reachability directly into the cnf.
	class CNFReachability: public Reach {
		ReachDetector & detector;
		bool over_approx;
		int num_updates = 0;

	public:
		CNFReachability(ReachDetector & _outer, bool over_approx) :
				detector(_outer), over_approx(over_approx) {
			
		}
		int numUpdates() const {
			return num_updates;
		}
		
		void update() {
			num_updates++;
		}
		void setSource(int s) {
			assert(s == detector.source);
		}
		int getSource() {
			return detector.source;
		}
		
		bool connected_unsafe(int t) {
			return connected(t);
		}
		bool connected_unchecked(int t) {
			return connected(t);
		}
		bool connected(int t) {
			assert(detector.reach_lits[t] != lit_Undef);
			if (over_approx)
				return detector.outer->value(detector.reach_lits[t]) != l_False;
			else
				return detector.outer->value(detector.reach_lits[t]) == l_True;
		}
		int previous(int node) {
			assert(false);
			exit(5); //not supported
		}
		int incomingEdge(int node) {
			assert(false);
			exit(5); //not supported
		}
		
	};
	void unassign(Lit l) {
		Detector::unassign(l);
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
	/*	Lit getLit(int node){

	 return reach_lits[node];

	 }*/

	void buildSATConstraints(bool onlyUnderApprox = false, int within_steps = -1);
	bool propagate(vec<Lit> & conflict);
	void buildReachReason(int node, vec<Lit> & conflict);
	void buildNonReachReason(int node, vec<Lit> & conflict, bool force_maxflow = false);
	void buildForcedEdgeReason(int reach_node, int forced_edge_id, vec<Lit> & conflict);
	void buildReason(Lit p, vec<Lit> & reason, CRef marker);
	bool checkSatisfied();
	void printSolution(std::ostream& write_to);

	void addLit(int from, int to, Var reach_var);
	Lit decide();
	void preprocess();
	void dbg_sync_reachability();

	ReachDetector(int _detectorID, GraphTheorySolver<Weight> * _outer, DynamicGraph<Weight>  &_g, DynamicGraph<Weight>  &_antig,
			int _source, double seed = 1); //:Detector(_detectorID),outer(_outer),within(-1),source(_source),rnd_seed(seed),positive_reach_detector(NULL),negative_reach_detector(NULL),positive_path_detector(NULL),positiveReachStatus(NULL),negativeReachStatus(NULL),chokepoint_status(*this),chokepoint(chokepoint_status, _antig,source){}
	virtual ~ReachDetector() {
		if (rnd_path)
			delete rnd_path;
		if (chokepoint_detector)
			delete chokepoint_detector;
		if (cutgraph_detector)
			delete cutgraph_detector;

		if (positiveReachStatus)
			delete positiveReachStatus;
		if (negativeReachStatus)
			delete negativeReachStatus;

		if(underapprox_path_detector && underapprox_path_detector != underapprox_detector)
			delete underapprox_path_detector;
		
		if(overapprox_path_detector && overapprox_path_detector != overapprox_reach_detector)
			delete overapprox_path_detector;

		if (underapprox_fast_detector && underapprox_fast_detector != underapprox_path_detector){
			delete underapprox_fast_detector;
		}

		if (negative_distance_detector && negative_distance_detector != overapprox_path_detector)
			delete negative_distance_detector;

		if (underapprox_detector)
			delete underapprox_detector;

		if (overapprox_reach_detector)
			delete overapprox_reach_detector;




		if(conflict_flow)
			delete conflict_flow;

		for (auto * c: conflict_flows){
			if (c)
				delete c;
		}
	}
	
	const char* getName() {
		return "Reachability Detector";
	}
	
	bool dbg_cut(std::vector<MaxFlowEdge> & cut, DynamicGraph<Weight>  & graph, int source, int node) {
#ifndef NDEBUG
		
		DynamicGraph<int>  t;
		for (int i = 0; i < graph.nodes(); i++)
			t.addNode();

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
					t.setEdgeWeight(id, 0);
					t.disableEdge(id);
				} else
					t.setEdgeWeight(id, 1);
			} else {
				t.setEdgeWeight(id, 0xFFFF);
			}
			
		}
		EdmondsKarpAdj<int> check(t,  source, node);
		std::vector<MaxFlowEdge> check_cut;
		int flow = check.minCut(check_cut);
		assert(flow < 0xFFFF);
		if (flow < 0xFFFF) {
			exit(4);
		}

#endif
		return true;
	}

	//Return the path (in terms of nodes)
	bool getModel_Path(int node, std::vector<int> & store_path);
	bool getModel_PathByEdgeLit(int node, std::vector<Lit> & store_path);
};
}
;
#endif /* REACHDETECTOR_H_ */
