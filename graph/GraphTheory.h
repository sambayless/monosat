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

#ifndef DGRAPH_H_
#define DGRAPH_H_

#include "utils/System.h"
#include "core/Theory.h"

#include "dgl/Reach.h"
#include "dgl/Dijkstra.h"
#include "dgl/BFS.h"

#include "core/SolverTypes.h"
#include "mtl/Map.h"
#include "dgl/MaxFlow.h"
#include "dgl/DynamicGraph.h"
#include "dgl/EdmondsKarp.h"
#include "dgl/EdmondsKarpAdj.h"
#include "dgl/EdmondsKarpDynamic.h"
#include "dgl/KohliTorr.h"
#include "dgl/Dinics.h"
#include "dgl/DinicsLinkCut.h"

#include "dgl/Chokepoint.h"
#include "WeightedDijkstra.h"
#include "GraphTheoryTypes.h"
#include "utils/System.h"
#include "core/Solver.h"

#include "AllPairsDetector.h"
#include "ReachDetector.h"
#include "comparison/ComparisonBVTheory.h"
#include "comparison/BitVector.h"
#include "DistanceDetector.h"
#include "MSTDetector.h"
#include "MaxflowDetector.h"
#include "ConnectedComponentsDetector.h"
#include "CycleDetector.h"
#include "SteinerDetector.h"
#include <vector>
#include <gmpxx.h>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <sstream>

using namespace dgl;
namespace Monosat {
template<typename Weight>
class GraphTheorySolver;

template<typename Weight>
class GraphTheorySolver: public Theory {
public:
	
	double rnd_seed;

private:
	Solver * S;
	int local_q = 0;
public:
	int id;
	bool all_edges_unit = true;
	bool all_edges_positive=true;
	vec<lbool> assigns;
	MSTDetector<Weight> * mstDetector = nullptr;
	vec<ReachabilityConstraint> unimplemented_reachability_constraints;

	DynamicGraph<Weight> g_under;
	DynamicGraph<Weight> g_over;
	/**
	 * The cutgraph is (optionally) used for conflict analysis by some graph theories.
	 * It has two edges for every edge in the real graph (with indices edgeID*2 and edgeID*2+1).
	 * If edge ID is assigned to FALSE, then edge ID*2 is assigned to enabled in the cutgraph, and
	 * edge ID*2+1 is disabled.
	 * Otherwise, if edge ID is unassigned or true, then edge ID*2 is disabled in the cutgraph, and
	 * edge ID*2+1 is enabled.
	 */
	DynamicGraph<long> cutGraph;

	struct ComparisonStatus;
	//if bitvectors weights are supplied, then this manages the resulting weights.
	ComparisonBVTheorySolver<long, ComparisonStatus> * comparator=nullptr;

	struct ComparisonStatus{
		GraphTheorySolver & outer;
		ComparisonBVTheorySolver<long, ComparisonStatus> & comp;
		void comparisonAltered(int bvID, int comparisonID){

		}
		void bvAltered(int bvID){
		/*	if(!outer.bv_needs_update[bvID]){
				outer.bv_needs_update[bvID]=true;
				outer.bvs_to_update.push(bvID);
			}*/
			int edgeID = outer.getBVEdge(bvID);
			outer.g_under.setEdgeWeight(edgeID,outer.edge_bv_weights[bvID].getUnder());
			outer.g_over.setEdgeWeight(edgeID, outer.edge_bv_weights[bvID].getOver());
		}

	};

	using BitVector = typename ComparisonBVTheorySolver<Weight, ComparisonStatus>::BitVector;

	vec<Assignment> trail;
	vec<int> trail_lim;

	struct ReachInfo {
		int source;
		bool distance = 0;
		Detector * detector;

		ReachInfo() :
			source(-1), detector(nullptr) {
		}
	};

	vec<ReachInfo> dist_info;
	vec<ReachInfo> reach_info;
	vec<ReachInfo> connect_info;
public:
	vec<Detector*> detectors;
	vec<ReachDetector<Weight>*> reach_detectors;
	vec<DistanceDetector<Weight>*> distance_detectors;
	vec<MaxflowDetector<Weight>*> flow_detectors;
	ConnectedComponentsDetector<Weight>* component_detector = nullptr;
	CycleDetector<Weight> * cycle_detector = nullptr;
	vec<SteinerDetector<Weight>*> steiner_detectors;

	vec<int> marker_map;

	std::vector<MaxFlowEdge> cut;

	//Full matrix
	//vec<vec<Edge> > edges;
	
	//Just a list of the edges
	vec<Edge> edge_list;
	vec<vec<Edge> > undirected_adj;
	vec<vec<Edge> > inv_adj;

	//vector of the weights for each edge
	std::vector<Weight> edge_weights;
	std::vector<BitVector> edge_bv_weights;

	struct BVInfo{
		int edgeID;
	};

	std::vector<BVInfo> bitvectors;
	vec<bool> bv_needs_update;
	vec<int> bvs_to_update;

	bool requiresPropagation = true;

	vec<char> seen;
	vec<int> to_visit;
	vec<Lit> tmp_clause;
	//Data about local theory variables, and how they connect to the sat solver's variables
	struct VarData {
		int isEdge :1;
		int isBV:1;
		int occursPositive :1;
		int occursNegative :1;
		int detector_edge :28;	//the detector this variable belongs to, or its edge number, if it is an edge variable
		Var solverVar;
	};

	vec<VarData> vars;
	int theory_index = 0;
public:
	
	double mctime = 0;
	double reachtime = 0;
	double unreachtime = 0;
	double pathtime = 0;
	double propagationtime = 0;
	long stats_propagations = 0;
	long stats_num_conflicts = 0;
	long stats_decisions = 0;
	long stats_num_reasons = 0;

	double reachupdatetime = 0;
	double unreachupdatetime = 0;
	double stats_initial_propagation_time = 0;
	double stats_decision_time = 0;
	double stats_reason_initial_time = 0;
	double stats_reason_time = 0;
	long num_learnt_paths = 0;
	long learnt_path_clause_length = 0;
	long num_learnt_cuts = 0;
	long learnt_cut_clause_length = 0;
	long stats_pure_skipped = 0;
	long stats_mc_calls = 0;
	long stats_propagations_skipped = 0;
	vec<Lit> reach_cut;

	struct CutStatus {
		int one = 1;
		int inf = 0xFFFF;
		GraphTheorySolver & outer;

		const int &operator [](int id) const {
			return one;
			/*if(outer.value(outer.edge_list[id].v) ==l_False){
			 return one;
			 }else{
			 return inf;
			 }*/
		}
		int size() const {
			return outer.edge_list.size();
		}
		CutStatus(GraphTheorySolver & _outer) :
				outer(_outer) {
		}
		
	} cutStatus;

	struct PropCutStatus {
		GraphTheorySolver & outer;
		int operator ()(int id) const {
			
			if (outer.value(outer.edge_list[id].v) == l_Undef) {
				return 1;
			} else {
				assert(outer.value(outer.edge_list[id].v)==l_True);
				return 0xF0F0F0;
			}
		}
		PropCutStatus(GraphTheorySolver & _outer) :
				outer(_outer) {
		}
		
	} propCutStatus;

	GraphTheorySolver(Solver * S_, int _id = -1) :
			S(S_), id(_id), cutStatus(*this), propCutStatus(*this) {
#ifdef RECORD
		{
			char t[30];
			sprintf(t, "TEST_GRAPH%d", id);
			g_under.outfile = fopen(t, "w");
		}
		{
			char t[30];
			sprintf(t, "TEST_ANTI_GRAPH%d", id);
			g_over.outfile = fopen(t, "w");
		}
		{
			char t[30];
			sprintf(t, "TEST_CUT_GRAPH%d", id);
			cutGraph.outfile = fopen(t, "w");
		}
#endif
		
		if (opt_adaptive_history_clear > 0) {
			g_under.adaptive_history_clear = true;
			g_over.adaptive_history_clear = true;
			cutGraph.adaptive_history_clear = true;
			g_under.historyClearInterval = opt_adaptive_history_clear;
			g_over.historyClearInterval = opt_adaptive_history_clear;
			cutGraph.historyClearInterval = opt_adaptive_history_clear;
		} else {
			g_under.historyClearInterval = opt_history_clear;
			g_over.historyClearInterval = opt_history_clear;
			cutGraph.historyClearInterval = opt_history_clear;
		}
		
		rnd_seed = opt_random_seed;
	}
	
	void printStats(int detailLevel) {
		if (detailLevel > 0) {
			for (Detector * d : detectors)
				d->printStats();
		}
		
		printf("Graph %d stats:\n", getGraphID());
		/*		printf("Decision Time: %f\n", stats_decision_time);
		 printf("Prop Time: %f (initial: %f)\n", propagationtime,stats_initial_propagation_time);
		 printf("Conflict Time: %f (initial: %f)\n", stats_reason_time,stats_reason_initial_time);*/

		/*	printf("Reach Time: %f (update Time: %f)\n", reachtime,reachupdatetime);
		 printf("Unreach Time: %f (Update Time: %f)\n", unreachtime,unreachupdatetime);
		 printf("Path Time: %f (#Paths: %d, AvgLength %f, total: %d)\n", pathtime, num_learnt_paths, (learnt_path_clause_length /  ((float) num_learnt_paths+1)),learnt_path_clause_length);
		 printf("Min-cut Time: %f (%d calls, %f average, #Cuts: %d, AvgLength %f, total: %d)\n", mctime, stats_mc_calls,(mctime/(stats_mc_calls ? stats_mc_calls:1)),  num_learnt_cuts, (learnt_cut_clause_length /  ((float) num_learnt_cuts+1)),learnt_cut_clause_length);
		 */
		printf("%d nodes, %d edges\n", g_under.nodes(), g_under.edges());
		printf("History Clears: over_approx %ld, under_approx %ld, cut_graph %ld\n", g_over.historyclears,
				g_under.historyclears, cutGraph.historyclears);
		printf("Propagations: %ld (%f s, avg: %f s, %ld skipped)\n", stats_propagations, propagationtime,
				(propagationtime) / ((double) stats_propagations + 1), stats_propagations_skipped);
		printf("Decisions: %ld (%f s, avg: %f s)\n", stats_decisions, stats_decision_time,
				(stats_decision_time) / ((double) stats_decisions + 1));
		printf("Conflicts: %ld\n", stats_num_conflicts);
		printf("Reasons: %ld (%f s, avg: %f s)\n", stats_num_reasons, stats_reason_time,
				(stats_reason_time) / ((double) stats_num_reasons + 1));
		
		fflush(stdout);
	}
	
	void writeTheoryWitness(std::ostream& write_to) {
		
		for (Detector * d : detectors) {
			write_to << "Graph " << this->getGraphID() << ", detector " << d->getID() << ":\n";
			d->printSolution(write_to);
		}
		
	}
	
	inline int getTheoryIndex() {
		return theory_index;
	}
	inline void setTheoryIndex(int id) {
		theory_index = id;
	}
	inline int getGraphID() {
		return id;
	}
	inline bool isEdgeVar(Var v) {
		assert(v < vars.size());
		return vars[v].isEdge;
	}
	inline int getEdgeID(Var v) {
		assert(isEdgeVar(v));
		return vars[v].detector_edge;
	}
	inline int getDetector(Var v) {
		assert(!isEdgeVar(v));
		return vars[v].detector_edge;
	}
	
	inline Var getEdgeVar(int edgeID) {
		Var v = edge_list[edgeID].v;
		assert(v < vars.size());
		assert(vars[v].isEdge);
		return v;
	}
	
	void makeEqual(Lit l1, Lit l2) {
		Lit o1 = toSolver(l1);
		Lit o2 = toSolver(l2);
		S->addClause(~o1, o2);
		S->addClause(o1, ~o2);
	}
	void makeEqualInSolver(Lit l1, Lit l2) {
		S->addClause(~l1, l2);
		S->addClause(l1, ~l2);
	}
	void addClause(Lit l1) {
		Lit o1 = toSolver(l1);
		S->addClause(o1);
	}
	void addClause(Lit l1, Lit l2) {
		Lit o1 = toSolver(l1);
		Lit o2 = toSolver(l2);
		S->addClause(o1, o2);
	}
	void addClause(Lit l1, Lit l2, Lit l3) {
		Lit o1 = toSolver(l1);
		Lit o2 = toSolver(l2);
		Lit o3 = toSolver(l3);
		S->addClause(o1, o2, o3);
	}
	void addClause(vec<Lit> & c) {
		tmp_clause.clear();
		c.copyTo(tmp_clause);
		toSolver(tmp_clause);
		S->addClause(tmp_clause);
	}
	void addClauseSafely(vec<Lit> & c) {
		tmp_clause.clear();
		c.copyTo(tmp_clause);
		toSolver(tmp_clause);
		
		S->addClauseSafely(tmp_clause);
	}
	
	Var newVar(int forDetector = -1, bool connectToTheory = false) {
		Var s = S->newVar();
		return newVar(s, forDetector, false, connectToTheory);
	}
	Var newBVVar(Var solverVar, int bvID, int edgeID){
		while (S->nVars() <= solverVar)
			S->newVar();
		Var v = vars.size();
		vars.push();
		vars[v].isEdge = false;
		vars[v].isBV=true;
		vars[v].detector_edge = bvID;
		vars[v].solverVar = solverVar;
		assigns.push(l_Undef);

		S->setTheoryVar(solverVar, getTheoryIndex(), v);
		assert(toSolver(v) == solverVar);

		return v;
	}

	Var newVar(Var solverVar, int detector, bool isEdge = false, bool connectToTheory = true) {
		while (S->nVars() <= solverVar)
			S->newVar();
		Var v = vars.size();
		vars.push();
		vars[v].isEdge = isEdge;
		vars[v].isBV=false;
		vars[v].detector_edge = detector;
		vars[v].solverVar = solverVar;
		assigns.push(l_Undef);
		if (connectToTheory) {
			S->setTheoryVar(solverVar, getTheoryIndex(), v);
			assert(toSolver(v) == solverVar);
		}
		if (!isEdge && detector >= 0)
			detectors[detector]->addVar(v);
		return v;
	}
	inline int level(Var v) {
		return S->level(toSolver(v));
	}
	inline int decisionLevel() {
		return trail_lim.size(); //S->decisionLevel();
	}
	inline int nVars() const {
		return vars.size(); //S->nVars();
	}
	inline Var toSolver(Var v) {
		//return v;
		assert(v < vars.size());
		//assert(S->hasTheory(vars[v].solverVar));
		//assert(S->getTheoryVar(vars[v].solverVar)==v);
		return vars[v].solverVar;
	}
	
	inline Lit toSolver(Lit l) {
		//assert(S->hasTheory(vars[var(l)].solverVar));
		//assert(S->getTheoryVar(vars[var(l)].solverVar)==var(l));
		return mkLit(vars[var(l)].solverVar, sign(l));
	}
	
	void toSolver(vec<Lit> & c) {
		for (int i = 0; i < c.size(); i++) {
			c[i] = toSolver(c[i]);
		}
	}
	
	inline lbool value(Var v) {
		if (assigns[v] != l_Undef)
			assert(S->value(toSolver(v)) == assigns[v]);
		
		return assigns[v]; //S->value(toSolver(v));
	}
	inline lbool value(Lit l) {
		if (assigns[var(l)] != l_Undef) {
			assert(S->value(toSolver(l)) == (assigns[var(l)] ^ sign(l)));
		}
		return assigns[var(l)] ^ sign(l);
	}
	inline lbool dbg_value(Var v) {
		return S->value(toSolver(v));
	}
	inline lbool dbg_value(Lit l) {
		return S->value(toSolver(l));
	}
	inline bool enqueue(Lit l, CRef reason) {
		assert(assigns[var(l)]==l_Undef);
		
		Lit sl = toSolver(l);
		if (S->enqueue(sl, reason)) {
			enqueueTheory(l);
			return true;
		} else {
			return false;
		}
	}
	
	~GraphTheorySolver() {
	}
	;
	int newNode() {
		
		inv_adj.push();
		undirected_adj.push();
		reach_info.push();
		connect_info.push();
		dist_info.push();
		g_over.addNode();
		cutGraph.addNode();
		seen.growTo(nNodes());
		
		return g_under.addNode();
	}
	void newNodes(int n) {
		for (int i = 0; i < n; i++)
			newNode();
	}
	int nNodes() {
		return g_under.nodes();
	}
	bool isNode(int n) {
		return n >= 0 && n < nNodes();
	}
	
	bool hasBitVector(int edgeID){
		return edge_list[edgeID].bvID>=0;
	}

	bool dbg_propgation(Lit l) {
#ifndef NDEBUG
		static vec<Lit> c;
		c.clear();
		for (int i = 0; i < S->trail.size(); i++) {
			if (!S->hasTheory(S->trail[i]) || S->getTheoryID(S->trail[i]) != getTheoryIndex())
				continue;
			Lit l = S->getTheoryLit(S->trail[i]);
			Var v = var(l);
			if (isEdgeVar(v)) {
				Edge & e = edge_list[getEdgeID(v)];
				c.push(l);
			}
			
			/*		if(v>=min_edge_var && v<min_edge_var+num_edges){
			 if(edge_list[v-min_edge_var].v<0)
			 continue;
			 Edge e = edge_list[v-min_edge_var];
			 c.push(l);
			 }*/
		}
		c.push(~l);
		//bool res = dbg->solve(c);
		//assert(~res);
#endif
		return true;
	}
	
	void dbg_sync_reachability() {
		
#ifdef DEBUG_GRAPH
		
		for(int i = 0;i<reach_detectors.size();i++) {
			ReachDetector* d = reach_detectors[i];
			d->dbg_sync_reachability();
		}
#endif
		
	}
	void dbg_sync() {
#ifndef NDEBUG
		if (opt_conflict_min_cut) {
			for (int i = 0; i < edge_list.size(); i++) {
				if (value(edge_list[i].v) == l_False) {
					assert(cutGraph.edgeEnabled(i * 2));
					assert(!cutGraph.edgeEnabled(i * 2 + 1));
				} else {
					assert(!cutGraph.edgeEnabled(i * 2));
					assert(cutGraph.edgeEnabled(i * 2 + 1));
				}
			}
		}
#endif
#ifdef DEBUG_DIJKSTRA
		
		for(int i = 0;i<assigns.size();i++) {
			lbool val = assigns[i];

			if(val!=l_Undef) {
				bool found=false;
				for(int j = 0;j<trail.size();j++) {
					if(trail[j].var==i) {
						assert(!found);
						assert(trail[j].assign== (val==l_True));
						found=true;
					}
				}
				assert(found);
			} else {
				for(int j = 0;j<trail.size();j++) {
					if(trail[j].var==i) {
						assert(false);
					}
				}
			}
			if(val!=l_Undef)
			assert(val==S->value(toSolver(i)));
		}
		/*static vec<lbool> assigned;
		 assigned.clear();*.
		 */

		/*for(int i = 0;i<edge_list.size();i++)
		 assigned.push(l_Undef);
		 int lev = 0;
		 int j =0;
		 for(int i = 0;i<trail.size();i++){
		 while(j<trail_lim.size() && i>=trail_lim[j]){
		 lev++;
		 j++;
		 }

		 Assignment e = trail[i];
		 if(e.isEdge){

		 Lit l = mkLit( e.var,!e.assign);
		 assert(value(l)==l_True);
		 int expected_level = level(var(l));
		 assert(level(var(l))==lev);
		 int edge_num = getEdgeID(e.var); //e.var-min_edge_var;
		 if(edge_list[edge_num].v<0)
		 continue;
		 assert(assigned[edge_num]==l_Undef);

		 assigned[edge_num] = sign(l)?l_False:l_True;
		 }else{
		 //Reachability assignment


		 }
		 }

		 for(int i = 0;i<assigned.size();i++){
		 assert(ass[i]== assigned[i]);
		 }

		 for(int i = 0;i<S->trail.size();i++){

		 Lit l = S->trail[i];
		 if(S->getTheoryID(l)==getTheoryIndex()){
		 l = S->getTheoryLit(l);
		 Var v = var(l);

		 int lev = level(v);
		 if(isEdgeVar(v)){
		 int edge_num = getEdgeID(v);
		 if(edge_list[edge_num].v<0)
		 continue;
		 lbool assigned_val=assigned[edge_num];
		 assert(assigned_val== (sign(l)?l_False:l_True));
		 }
		 if(v>= min_edge_var && v<min_edge_var+num_edges){
		 int edge_num = v-min_edge_var;
		 if(edge_list[edge_num].v<0)
		 continue;
		 lbool assigned_val=assigned[edge_num];
		 assert(assigned_val== (sign(l)?l_False:l_True));
		 }

		 }
		 }*/

#endif
	}
	void dbg_full_sync() {
#ifdef DEBUG_GRAPH
		dbg_sync();
		for(int i =0;i<edge_list.size();i++) {

			if(edge_list[i].edgeID>=0 && g_under.edgeEnabled(i)) {
				assert(value(edge_list[i].v)==l_True);
			} else if (edge_list[i].edgeID>=0 && !g_over.edgeEnabled(i)) {
				assert(value(edge_list[i].v)==l_False);
			}
		}

#endif
	}
	
	void backtrackUntil(int level) {
		static int it = 0;
		
		bool changed = false;
		//need to remove and add edges in the two graphs accordingly.
		if (trail_lim.size() > level) {
			
			int stop = trail_lim[level];
			for (int i = trail.size() - 1; i >= trail_lim[level]; i--) {
				
				Assignment & e = trail[i];
				assert(assigns[e.var]!=l_Undef);
				if (e.isEdge) {
					assert(dbg_value(e.var)==l_Undef);
					int edge_num = getEdgeID(e.var); //e.var-min_edge_var;
							
					if (e.assign) {
						g_under.disableEdge(e.from, e.to, edge_num);
						assert(!cutGraph.edgeEnabled(edge_num * 2));
					} else {
						g_over.enableEdge(e.from, e.to, edge_num);
						assert(g_over.hasEdge(e.from, e.to));
						if (opt_conflict_min_cut) {
							assert(cutGraph.edgeEnabled(edge_num * 2));
							cutGraph.disableEdge(e.from, e.to, edge_num * 2);
							assert(!cutGraph.edgeEnabled(edge_num * 2 + 1));
							cutGraph.enableEdge(e.from, e.to, edge_num * 2 + 1);
						}
					}
				} else {
					//This is a reachability literal				  
					detectors[getDetector(e.var)]->unassign(mkLit(e.var, !e.assign));
				}
				assigns[e.var] = l_Undef;
				changed = true;
			}
			trail.shrink(trail.size() - stop);
			trail_lim.shrink(trail_lim.size() - level);
			assert(trail_lim.size() == level);
			
			if (changed) {
				requiresPropagation = true;
				/*				g.markChanged();
				 antig.markChanged();
				 cutGraph.markChanged();*/
			}
			
			for (Detector * d : detectors) {
				d->backtrack(level);
			}
		}
		/*		if(local_q>S->qhead)
		 local_q=S->qhead;*/
		assert(dbg_graphsUpToDate());
		/*for(int i = 0;i<reach_detectors.size();i++){
		 if(reach_detectors[i]->positive_reach_detector)
		 reach_detectors[i]->positive_reach_detector->update();
		 if(reach_detectors[i]->negative_reach_detector)
		 reach_detectors[i]->negative_reach_detector->update();
		 }*/

		dbg_sync();
		
	}
	;

	Lit decideTheory() {
		if (!opt_decide_theories)
			return lit_Undef;
		double start = rtime(1);
		
		dbg_full_sync();
		for (int i = 0; i < detectors.size(); i++) {
			Detector * r = detectors[i];
			Lit l = r->decide(decisionLevel());
			if (l != lit_Undef) {
				stats_decisions++;
				r->stats_decisions++;
				stats_decision_time += rtime(1) - start;
				return toSolver(l);
			}
		}
		stats_decision_time += rtime(1) - start;
		return lit_Undef;
	}
	
	void backtrackUntil(Lit p) {
		//need to remove and add edges in the two graphs accordingly.
		
		int i = trail.size() - 1;
		for (; i >= 0; i--) {
			Assignment e = trail[i];
			if (var(p) == e.var) {
					assert(sign(p) != e.assign);
					break;
				}
			if (e.isEdge) {
				int edge_num = getEdgeID(e.var); //e.var-min_edge_var;
				assert(assigns[e.var]!=l_Undef);
				assigns[e.var] = l_Undef;
				if (e.assign) {
					g_under.disableEdge(e.from, e.to, edge_num);
					
					assert(!cutGraph.edgeEnabled(edge_num * 2));
				} else {
					g_over.enableEdge(e.from, e.to, edge_num);
					assert(g_over.hasEdge(e.from, e.to));
					if (opt_conflict_min_cut) {
						assert(cutGraph.edgeEnabled(edge_num * 2));
						cutGraph.disableEdge(e.from, e.to, edge_num * 2);
						assert(!cutGraph.edgeEnabled(edge_num * 2 + 1));
						cutGraph.enableEdge(e.from, e.to, edge_num * 2 + 1);
					}
				}
			} else {

				assigns[e.var] = l_Undef;
				detectors[getDetector(e.var)]->unassign(mkLit(e.var, !e.assign));
			}
		}
		
		trail.shrink(trail.size() - (i + 1));
		//if(i>0){
		requiresPropagation = true;
		/*			g.markChanged();
		 antig.markChanged();
		 cutGraph.markChanged();*/
		//}
		for (Detector * d : detectors) {
			d->backtrack(this->decisionLevel());
		}
		//while(trail_lim.size() && trail_lim.last()>=trail.size())
		//	trail_lim.pop();
		dbg_sync();
		/*		for(int i = 0;i<reach_detectors.size();i++){
		 if(reach_detectors[i]->positive_reach_detector)
		 reach_detectors[i]->positive_reach_detector->update();
		 if(reach_detectors[i]->negative_reach_detector)
		 reach_detectors[i]->negative_reach_detector->update();
		 }*/
	}
	;

	void newDecisionLevel() {
		trail_lim.push(trail.size());
	}
	;

	void buildReason(Lit p, vec<Lit> & reason) {
		CRef marker = S->reason(var(toSolver(p)));
		assert(marker != CRef_Undef);
		int pos = CRef_Undef - marker;
		int d = marker_map[pos];
		//double initial_start = rtime(1);
		double start = rtime(1);
		backtrackUntil(p);
		
		assert(d < detectors.size());
		detectors[d]->buildReason(p, reason, marker);
		toSolver(reason);
		double finish = rtime(1);
		stats_reason_time += finish - start;
		stats_num_reasons++;
		//stats_reason_initial_time+=start-initial_start;
		
	}
	
	bool dbg_reachable(int from, int to, bool undirected = false) {
#ifdef DEBUG_DIJKSTRA
		
		if(undirected) {
			UnweightedDijkstra<Reach::NullStatus, true> d(from,g_under);
			d.update();
			return d.connected(to);
		} else {
			UnweightedDijkstra<> d(from,g_under);
			d.update();
			return d.connected(to);
		}

#else
		return true;
#endif
	}
	
	bool dbg_notreachable(int from, int to, bool undirected = false) {
		
#ifndef NDEBUG
		//drawFull(from,to);
		
		DynamicGraph<Weight> g;
		for (int i = 0; i < nNodes(); i++) {
			g.addNode();
		}
		
		for (int i = 0; i < edge_list.size(); i++) {
			if (edge_list[i].v < 0)
				continue;
			Edge e = edge_list[i];
			if (value(e.v) != l_False) {
				Weight w = edge_weights[e.edgeID];
				g.addEdge(e.from, e.to,e.edgeID,w);
			}
		}
		if (undirected) {
			UnweightedDijkstra<Weight, Reach::NullStatus, true> d(from, g);
			
			return !d.connected(to);
		} else {
			UnweightedDijkstra<Weight,Reach::NullStatus, false> d(from, g);
			
			return !d.connected(to);
		}
#endif
		return true;
		
	}
	
	bool dbg_graphsUpToDate() {
#ifdef DEBUG_GRAPH
		for(int i = 0;i<edge_list.size();i++) {
			if(edge_list[i].v<0)
			continue;
			Edge e = edge_list[i];
			lbool val = value(e.v);

			if(val==l_True || val==l_Undef) {
				assert(g_over.edgeEnabled(i));
				//assert(antig.hasEdge(e.from,e.to));
			} else {
				assert(!g_over.edgeEnabled(i));
				//assert(!antig.hasEdge(e.from,e.to));
			}
			if(val==l_True) {
				assert(g_under.edgeEnabled(i));
				//assert(g.hasEdge(e.from,e.to));
			} else {
				assert(!g_under.edgeEnabled(i));
				//assert(!g.hasEdge(e.from,e.to));
			}
		}

#endif
		return true;
	}
	
	/*	int getEdgeID(Var v){
	 assert(v>= min_edge_var && v<min_edge_var+edge_list.size());

	 //this is an edge assignment
	 int edge_num = v-min_edge_var;
	 return edge_num;
	 }*/

	void preprocess() {
		for (int i = 0; i < detectors.size(); i++) {
			detectors[i]->preprocess();
		}
	}
	void setLiteralOccurs(Lit l, bool occurs) {
		if (isEdgeVar(var(l))) {
			//don't do anything
		} else {
			//this is a graph property detector var
			if (!sign(l) && vars[var(l)].occursPositive != occurs)
				detectors[getDetector(var(l))]->setOccurs(l, occurs);
			else if (sign(l) && vars[var(l)].occursNegative != occurs)
				detectors[getDetector(var(l))]->setOccurs(l, occurs);
		}
		
	}
	
	void enqueueTheory(Lit l) {
		Var v = var(l);
		
		int lev = level(v);
		
		assert(decisionLevel() <= lev);
		
		while (lev > trail_lim.size()) {
			newDecisionLevel();
		}
		
		if (assigns[var(l)] != l_Undef) {
			return;			//this is already enqueued.
		}
		assert(assigns[var(l)]==l_Undef);
		assigns[var(l)] = sign(l) ? l_False : l_True;
		requiresPropagation = true;
		//printf("enqueue %d\n", dimacs(l));
		
#ifndef NDEBUG
		{
			for (int i = 0; i < trail.size(); i++) {
				assert(trail[i].var != v);
			}
		}
#endif
		
#ifdef RECORD
		if (g_under.outfile) {
			fprintf(g_under.outfile, "enqueue %d\n", dimacs(l));
			
			fprintf(g_under.outfile, "\n");
			fflush(g_under.outfile);
		}
		if (g_over.outfile) {
			fprintf(g_over.outfile, "enqueue %d\n", dimacs(l));
			fprintf(g_over.outfile, "\n");
			fflush(g_over.outfile);
		}
#endif
		
		if(isBVVar(var(l))){
			int bvID = getBV(var(l));
			if(isBVEdge(bvID)){
				int edgeID = getBVEdge(bvID);
			}
			if(!bv_needs_update[bvID]){
				bv_needs_update[bvID]=true;
				bvs_to_update.push(bvID);
			}

		}else if (isEdgeVar(var(l))) {
			
			//this is an edge assignment
			int edge_num = getEdgeID(var(l)); //v-min_edge_var;
			assert(edge_list[edge_num].v == var(l));
			
			int from = edge_list[edge_num].from;
			int to = edge_list[edge_num].to;
			trail.push( { true, !sign(l), from, to, v });
			
			Assignment e = trail.last();
			assert(e.from == from);
			assert(e.to == to);
			
			if (!sign(l)) {
				g_under.enableEdge(from, to, edge_num);
			} else {
				g_over.disableEdge(from, to, edge_num);
				if (opt_conflict_min_cut) {
					assert(cutGraph.edgeEnabled(edge_num * 2 + 1));
					assert(!cutGraph.edgeEnabled(edge_num * 2));
					cutGraph.enableEdge(e.from, e.to, edge_num * 2);
					cutGraph.disableEdge(e.from, e.to, edge_num * 2 + 1);
				}
			}
			
		} else {
			
			trail.push( { false, !sign(l), 0, 0, v });
			//this is an assignment to a non-edge atom. (eg, a reachability assertion)
			detectors[getDetector(var(l))]->assign(l);
		}
		
	}
	;
	bool propagateTheory(vec<Lit> & conflict) {
		static int itp = 0;
		if (++itp == 62279) {
			int a = 1;
		}
		stats_propagations++;
		dbg_sync();
		if (!requiresPropagation) {
			stats_propagations_skipped++;
			assert(dbg_graphsUpToDate());
			return true;
		}

		bool any_change = false;
		double startproptime = rtime(1);
		
		conflict.clear();
		//Can probably speed this up alot by a) constant propagating reaches that I care about at level 0, and b) Removing all detectors for nodes that appear only in the opposite polarity (or not at all) in the cnf.
		//That second one especially.
		
		//At level 0, need to propagate constant reaches/source nodes/edges...

		for (int bvID:bvs_to_update){
			assert(bv_needs_update[bvID]);
			if(isBVEdge(bvID)){
				int edgeID = getBVEdge(bvID);
				g_under.setEdgeWeight(edgeID,edge_bv_weights[bvID].getUnder());
				g_over.setEdgeWeight(edgeID, edge_bv_weights[bvID].getOver());
				bv_needs_update[bvID]=false;
			}
		}
		bvs_to_update.clear();

		dbg_sync();
		assert(dbg_graphsUpToDate());
		
		for (int d = 0; d < detectors.size(); d++) {
			assert(conflict.size() == 0);
			bool r = detectors[d]->propagate(conflict);
			if (!r) {
				stats_num_conflicts++;
				toSolver(conflict);
				propagationtime += rtime(1) - startproptime;
				return false;
			}
		}
		
		dbg_full_sync();
		
		requiresPropagation = false;
		g_under.clearChanged();
		g_over.clearChanged();
		cutGraph.clearChanged();
		
		g_under.clearHistory();
		g_over.clearHistory();
		cutGraph.clearHistory();
		//detectors_to_check.clear();
		
		double elapsed = rtime(1) - startproptime;
		propagationtime += elapsed;
		dbg_sync();
		dbg_sync_reachability();
		return true;
	}
	;

	bool solveTheory(vec<Lit> & conflict) {
		requiresPropagation = true;		//Just to be on the safe side... but this shouldn't really be required.
		bool ret = propagateTheory(conflict);
		//Under normal conditions, this should _always_ hold (as propagateTheory should have been called and checked by the parent solver before getting to this point).
		assert(ret);
		return ret;
	}
	;

	void drawFull(int from, int to) {
		printf("digraph{\n");
		for (int i = 0; i < nNodes(); i++) {
			if (i == from) {
				printf("n%d [label=\"From\", style=filled, fillcolor=blue]\n", i);
			} else if (i == to) {
				printf("n%d [label=\"To\", style=filled, fillcolor=red]\n", i);
			} else
				printf("n%d\n", i);
		}
		
		for (int i = 0; i < edge_list.size(); i++) {
			if (edge_list[i].v < 0)
				continue;
			Edge & e = edge_list[i];
			const char * s = "black";
			if (value(e.v) == l_True)
				s = "blue";
			else if (value(e.v) == l_False)
				s = "red";
			printf("n%d -> n%d [label=\"v%d\",color=\"%s\"]\n", e.from, e.to, e.v, s);
		}
		
		printf("}\n");
	}
	void drawFull() {
		printf("digraph{\n");
		for (int i = 0; i < nNodes(); i++) {
			printf("n%d\n", i);
		}
		
		for (int i = 0; i < edge_list.size(); i++) {
			Edge & e = edge_list[i];
			const char * s = "black";
			if (value(e.v) == l_True)
				s = "blue";
			else if (value(e.v) == l_False)
				s = "red";
			else {
				int a = 1;
			}
			printf("n%d -> n%d [label=\"v%d\",color=\"%s\"]\n", e.from, e.to, e.v, s);
		}
		
		printf("}\n");
	}
	
	bool check_solved() {
		if (opt_print_graph) {
			drawFull();
		}
		for (int i = 0; i < edge_list.size(); i++) {
			if (edge_list[i].v < 0)
				continue;
			Edge & e = edge_list[i];
			lbool val = value(e.v);
			if (val == l_Undef) {
				return false;
			}
			
			if (val == l_True) {
				/*	if(!g.hasEdge(e.from,e.to)){
				 return false;
				 }
				 if(!antig.hasEdge(e.from,e.to)){
				 return false;
				 }*/
				if (!g_under.edgeEnabled(e.edgeID)) {
					return false;
				}
				if (!g_over.edgeEnabled(e.edgeID)) {
					return false;
				}
			} else {
				/*if(g.hasEdge(e.from,e.to)){
				 return false;
				 }*/
				if (g_under.edgeEnabled(e.edgeID)) {
					return false;
				}
				if (g_over.edgeEnabled(e.edgeID)) {
					return false;
				}
				/*if(antig.hasEdge(e.from,e.to)){
				 return false;
				 }*/

			}
			
		}
		for (int i = 0; i < detectors.size(); i++) {
			if (!detectors[i]->checkSatisfied()) {
				return false;
			}
		}
		return true;
	}
	
	bool dbg_solved() {
#ifdef DEBUG_GRAPH
		for(int i = 0;i<edge_list.size();i++) {
			if(edge_list[i].v<0)
			continue;
			Edge & e = edge_list[i];
			lbool val = value(e.v);
			assert(val!=l_Undef);

			if(val==l_True) {
				assert(g_under.hasEdge(e.from,e.to));
				assert(g_over.hasEdge(e.from,e.to));
			} else {
				assert(!g_under.hasEdge(e.from,e.to));
				assert(!g_over.hasEdge(e.from,e.to));
			}

		}
		assert(check_solved());
		/*for(int i = 0;i<reach_detectors.size();i++){
		 ReachDetector* d  = reach_detectors[i];
		 for(int j = 0;j< d->reach_lits.size();j++){
		 Lit l = d->reach_lits[j];
		 if(l!=lit_Undef){
		 int node = d->getNode(var(l));

		 if(value(l)==l_True){
		 assert(d->positive_reach_detector->connected(node));
		 }else if (value(l)==l_False){
		 assert(! d->negative_reach_detector->connected(node));
		 }else{
		 assert(!d->positive_reach_detector->connected(node));
		 assert(d->negative_reach_detector->connected(node));
		 }
		 }
		 }
		 }*/
#endif
		return true;
	}
	
	void drawCurrent() {
		int from = -1;
		int to = -1;
		printf("digraph{\n");
		for (int i = 0; i < nNodes(); i++) {
			if (i == from) {
				printf("n%d [label=\"From\", style=filled, fillcolor=blue]\n", i);
			} else if (i == to) {
				printf("n%d [label=\"To\", style=filled, fillcolor=red]\n", i);
			} else
				printf("n%d\n", i);
		}

		for (int i = 0; i < edge_list.size(); i++) {
			if (edge_list[i].v < 0)
				continue;
			Edge & e = edge_list[i];
			const char * s = "black";
			if (value(e.v) == l_True)
				s = "blue";
			else if (value(e.v) == l_False)
				continue;
				//s = "red";
			printf("n%d -> n%d [label=\"v%d\",color=\"%s\"]\n", e.from, e.to, e.v, s);
		}
		
		printf("}\n");
	}
	int nEdges() {
		return edge_list.size();
	}
	CRef newReasonMarker(int detectorID) {
		CRef reasonMarker = S->newReasonMarker(this);
		int mnum = CRef_Undef - reasonMarker;
		marker_map.growTo(mnum + 1);
		marker_map[mnum] = detectorID;
		return reasonMarker;
	}

	bool isBVEdge(int bvID){
		return bitvectors[bvID].edgeID>=0;
	}
	int getBVEdge(int bvID){
		assert(isBVEdge(bvID));
		return bitvectors[bvID].edgeID;
	}

	bool isBVVar(Var v){
		return vars[v].isBV;
	}

	int getBV(Var v){
		assert(isBVVar(v));
		return vars[v].detector_edge;
	}

	Lit newEdge(int from, int to, Var outerVar,vec<Var> & bitVector) {
			assert(outerVar!=var_Undef);
			assert(edge_weights.size()==0);
			/*	if(outerVar==var_Undef)
			 outerVar = S->newVar();*/
			vec<Var> internalBV;
			int bvID = bitvectors.size();
			int index = edge_list.size();
			for (Var v:bitVector){
				internalBV.push(newBVVar(v,bvID,index));
			}

		/*	if(!comparator){
				comparator = new ComparisonBVTheorySolver<Weight,ComparisonStatus>(ComparisonStatus(*this, *comparator));
			}*/
			BitVector bv = comparator->newBitvector(bvID,bitVector);
			bitvectors.push(bvID);

			all_edges_unit &= (bv.getUnder()== 1 && bv.getOver()==1);
			all_edges_positive &= bv.getUnder()>0;
			edge_list.push();
			Var v = newVar(outerVar, index, true);

			undirected_adj[to].push( { v, outerVar, from, to, index });
			undirected_adj[from].push( { v, outerVar, to, from, index });
			inv_adj[to].push( { v, outerVar, from, to, index });

			//num_edges++;
			edge_list[index].v = v;
			edge_list[index].outerVar = outerVar;
			edge_list[index].from = from;
			edge_list[index].to = to;
			edge_list[index].edgeID = index;
			edge_list[index].bvID = bv.getID();

			if (edge_bv_weights.size() <= index) {
				edge_bv_weights.resize(index + 1);
			}
			edge_bv_weights[index] = bv;

			//edges[from][to]= {v,outerVar,from,to,index};
			g_under.addEdge(from, to, index,bv.getUnder());
			g_under.disableEdge(from, to, index);
			g_over.addEdge(from, to, index,bv.getOver());
			cutGraph.addEdge(from, to, index * 2,1);
			cutGraph.addEdge(from, to, index * 2 + 1,0xFFFF);
			cutGraph.disableEdge(from, to, index * 2);

	#ifdef RECORD
			if (g_under.outfile) {

				fprintf(g_under.outfile, "edge_bv_weight %d", index);
				for(Var v:bitVector)
					fprintf(g_under.outfile," %d", v+1);
				fprintf(g_under.outfile,"\n");
				fflush(g_under.outfile);
			}
			if (g_over.outfile) {
				fprintf(g_over.outfile, "edge_bv_weight %d", index);
				for(Var v:bitVector)
					fprintf(g_over.outfile," %d", v+1);
				fprintf(g_over.outfile,"\n");
				fflush(g_over.outfile);
			}
			if (cutGraph.outfile) {

				fprintf(cutGraph.outfile, "edge_weight %d %d\n", index * 2, 1);
				fprintf(cutGraph.outfile, "edge_weight %d %d\n", index * 2 + 1, 0xFFFF);
				fflush(cutGraph.outfile);
			}
	#endif
			return mkLit(v, false);
		}
	Lit newEdge(int from, int to, Var outerVar = var_Undef, Weight weight = 1) {
		assert(outerVar!=var_Undef);
		/*	if(outerVar==var_Undef)
		 outerVar = S->newVar();*/
		assert(bitvectors.size()==0);
		all_edges_unit &= (weight == 1);
		all_edges_positive &= weight>0;
		int index = edge_list.size();
		edge_list.push();
		Var v = newVar(outerVar, index, true);
		
		/*
		 if(num_edges>0){
		 }else
		 min_edge_var=v;

		 int index = v-min_edge_var;*/
		/*
		 while(edge_list.size()<=index){
		 edge_list.push({-1,-1,-1,-1,-1,1});
		 assigns.push(l_Undef);
		 }*/
		undirected_adj[to].push( { v, outerVar, from, to, index });
		undirected_adj[from].push( { v, outerVar, to, from, index });
		inv_adj[to].push( { v, outerVar, from, to, index });
		
		//num_edges++;
		edge_list[index].v = v;
		edge_list[index].outerVar = outerVar;
		edge_list[index].from = from;
		edge_list[index].to = to;
		edge_list[index].edgeID = index;
		//edge_list[index].weight=weight;
		if (edge_weights.size() <= index) {
			edge_weights.resize(index + 1);
		}
		edge_weights[index] = weight;
		
		//edges[from][to]= {v,outerVar,from,to,index};
		g_under.addEdge(from, to, index,weight);
		g_under.disableEdge(from, to, index);
		g_over.addEdge(from, to, index,weight);
		cutGraph.addEdge(from, to, index * 2,1);
		cutGraph.addEdge(from, to, index * 2 + 1,0xFFFF);
		cutGraph.disableEdge(from, to, index * 2);
		
#ifdef RECORD
		if (g_under.outfile) {
			std::stringstream wt;
			wt << weight;
			fprintf(g_under.outfile, "edge_weight %d %s\n", index, wt.str().c_str());
			fflush(g_under.outfile);
		}
		if (g_over.outfile) {
			std::stringstream wt;
			wt << weight;
			fprintf(g_over.outfile, "edge_weight %d %s\n", index, wt.str().c_str());
			fflush(g_over.outfile);
		}
		if (cutGraph.outfile) {
			
			fprintf(cutGraph.outfile, "edge_weight %d %d\n", index * 2, 1);
			fprintf(cutGraph.outfile, "edge_weight %d %d\n", index * 2 + 1, 0xFFFF);
			fflush(cutGraph.outfile);
		}
#endif
		return mkLit(v, false);
	}
	/*	int getEdgeID(int from, int to){
	 assert(edges[from][to].edgeID>=0);
	 return edges[from][to].edgeID;
	 }*/
/*	Weight getWeight(int edgeID) {
		return edge_weights[edgeID];
	}*/
	void reachesWithinSteps(int from, int to, Var reach_var, int within_steps) {
		
		assert(from < g_under.nodes());
		if (within_steps <= -1)
			within_steps = g_under.nodes();
		
		if (dist_info[from].source < 0) {
			DistanceDetector<Weight> * d = new DistanceDetector<Weight>(detectors.size(), this,  g_under, g_over,
					from, drand(rnd_seed));
			detectors.push(d);
			distance_detectors.push(d);
			assert(detectors.last()->getID() == detectors.size() - 1);
			dist_info[from].source = from;
			dist_info[from].detector = detectors.last();
		}
		
		DistanceDetector<Weight> * d = (DistanceDetector<Weight>*) dist_info[from].detector;
		assert(d);
		
		d->addUnweightedShortestPathLit(from, to, reach_var, within_steps);
		
	}
	void reachesWithinDistance(int from, int to, Var reach_var, Weight distance) {
		
		assert(from < g_under.nodes());
		
		if (dist_info[from].source < 0) {
			DistanceDetector<Weight> * d = new DistanceDetector<Weight>(detectors.size(), this,  g_under, g_over,
					from, drand(rnd_seed));
			detectors.push(d);
			distance_detectors.push(d);
			assert(detectors.last()->getID() == detectors.size() - 1);
			dist_info[from].source = from;
			dist_info[from].detector = detectors.last();
		}
		
		DistanceDetector<Weight> * d = (DistanceDetector<Weight>*) dist_info[from].detector;
		assert(d);
		
		d->addWeightedShortestPathLit(from, to, reach_var, distance);
		
	}
	void implementConstraints() {
		if (!S->okay())
			return;
		if (opt_allpairs_percentage >= 1) {
			for (int i = 0; i < unimplemented_reachability_constraints.size(); i++) {
				ReachabilityConstraint c = unimplemented_reachability_constraints[i];
				reaches_private(c.from, c.to, c.reach_var, c.distance);
			}
			
		} else if (opt_allpairs_percentage == 0) {
			for (int i = 0; i < unimplemented_reachability_constraints.size(); i++) {
				ReachabilityConstraint c = unimplemented_reachability_constraints[i];
				allpairs(c.from, c.to, c.reach_var, c.distance);
			}
		} else {
			{
				vec<bool> seen;
				int count = 0;
				seen.growTo(nNodes());
				for (int i = 0; i < unimplemented_reachability_constraints.size(); i++) {
					ReachabilityConstraint c = unimplemented_reachability_constraints[i];
					if (!seen[c.from]) {
						seen[c.from] = true;
						count++;
					}
				}
				double frac = ((double) count) / ((double) nNodes());
				
				if (opt_verb > 0 && frac >= opt_allpairs_percentage) {
					printf("Allpairs solver triggered for graph %d by percentage of source nodes: %d/%d=%f>%f\n",
							getGraphID(), count, nNodes(), frac, (double) opt_allpairs_percentage);
				}
				
				for (int i = 0; i < unimplemented_reachability_constraints.size(); i++) {
					ReachabilityConstraint c = unimplemented_reachability_constraints[i];
					if (frac >= opt_allpairs_percentage)
						allpairs(c.from, c.to, c.reach_var, c.distance);
					else
						reaches_private(c.from, c.to, c.reach_var, c.distance);
				}
			}
		}
		unimplemented_reachability_constraints.clear();
		
	}
	void allpairs_undirected(int from, int to, Var reach_var, int within_steps = -1) {
		
	}
	
	void allpairs(int from, int to, Var reach_var, int within_steps = -1) {
		//for now, reachesWithinSteps to be called instead
		
		assert(from < g_under.nodes());
		if (within_steps > g_under.nodes())
			within_steps = -1;
		
		if (reach_info[from].source < 0) {
			
			detectors.push(new AllPairsDetector<Weight>(detectors.size(), this, g_under, g_over, drand(rnd_seed)));
			//reach_detectors.push(reach_detectors.last());
			
			assert(detectors.last()->getID() == detectors.size() - 1);
			
			//reach_detectors.last()->negative_dist_detector = new Dijkstra(from,antig);
			//reach_detectors.last()->source=from;
			
			reach_info[from].source = from;
			reach_info[from].detector = detectors.last();
			
			//reach_detectors.last()->within=within_steps;
			
		}
		
		AllPairsDetector<Weight> * d = (AllPairsDetector<Weight>*) reach_info[from].detector;
		assert(d);
		
		d->addLit(from, to, reach_var, within_steps);
		
	}
	
	void reaches_private(int from, int to, Var reach_var, int within_steps = -1) {
		//for now, reachesWithinSteps to be called instead
		if (within_steps >= 0 || opt_force_distance_solver) {
			reachesWithinSteps(from, to, reach_var, within_steps);
			return;
		}
		
		assert(from < g_under.nodes());
		if (within_steps > g_under.nodes())
			within_steps = -1;
		
		if (reach_info[from].source < 0) {
			
			ReachDetector<Weight>*rd = new ReachDetector<Weight>(detectors.size(), this, g_under, g_over, from,
					drand(rnd_seed));
			detectors.push(rd);
			reach_detectors.push(rd);
			
			assert(detectors.last()->getID() == detectors.size() - 1);
			
			//reach_detectors.last()->negative_dist_detector = new Dijkstra(from,antig);
			//reach_detectors.last()->source=from;
			
			reach_info[from].source = from;
			reach_info[from].detector = detectors.last();
			
			//reach_detectors.last()->within=within_steps;
			
		}
		
		ReachDetector<Weight> * d = (ReachDetector<Weight>*) reach_info[from].detector;
		assert(d);
		assert(within_steps == -1);
		d->addLit(from, to, reach_var);
		
	}
	
	void reaches(int from, int to, Var reach_var, int within_steps = -1) {
		unimplemented_reachability_constraints.push( { from, to, within_steps, reach_var });
		//to allow us to alter the solving algorithm based on the number and type of constraints, we aren't implementing them here directly any more - instead,
		//we just store the constraints in this vector, then implement them later when 'implementConstraints' is called.
	}
	
	void reachesAny(int from, Var firstVar, int within_steps = -1) {
		for (int i = 0; i < g_under.nodes(); i++) {
			reaches(from, i, firstVar + i, within_steps);
		}
	}
	
	void reachesAny(int from, vec<Lit> & reachlits_out, int within_steps = -1) {
		for (int i = 0; i < g_under.nodes(); i++) {
			Var reachVar = S->newVar();
			//reaches(from,i,reachVar,within_steps);
			reaches(from, i, reachVar, within_steps);
			reachlits_out.push(mkLit(reachVar, false));
		}
	}


	//v will be true if the minimum weight is <= the specified value
	void minimumSpanningTree(Var v, Weight minimum_weight, bool inclusive) {
		if (!mstDetector) {
			mstDetector = new MSTDetector<Weight>(detectors.size(), this, g_under, g_over,  drand(rnd_seed));
			detectors.push(mstDetector);
		}
		mstDetector->addWeightLit(v, minimum_weight,inclusive);
	}
	void edgeInMinimumSpanningTree(Var edgeVar, Var var) {
		if (!mstDetector) {
			mstDetector = new MSTDetector<Weight>(detectors.size(), this, g_under, g_over,  drand(rnd_seed));
			detectors.push(mstDetector);
		}
		if (!S->hasTheory(edgeVar) || (S->getTheoryID(edgeVar) != getTheoryIndex())
				|| !isEdgeVar(S->getTheoryVar(edgeVar))) {
			fprintf(stderr, "%d is not an edge variable for theory %d! Aborting\n", edgeVar + 1, getTheoryIndex());
			exit(1);
		}
		edgeVar = S->getTheoryVar(edgeVar);
		int edgeid = getEdgeID(edgeVar);
		assert(edgeid >= 0);
		if (edge_list[edgeid].v == var_Undef) {
			printf("MST edge constraint for undefined edge %d with variable %d, aborting.\n", edgeid, edgeVar + 1);
			exit(1);
		}
		mstDetector->addTreeEdgeLit(edgeid, var);
	}
	void maxFlow(int from, int to, Weight  max_flow, Var v, bool inclusive=true) {
		
		for (int i = 0; i < flow_detectors.size(); i++) {
			if (flow_detectors[i]->source == from && flow_detectors[i]->target == to) {
				flow_detectors[i]->addFlowLit(max_flow, v,inclusive);
				return;
			}
		}
		MaxflowDetector<Weight> *f = new MaxflowDetector<Weight>(detectors.size(), this,  g_under, g_over, from,
				to, drand(rnd_seed));
		flow_detectors.push(f);
		detectors.push(f);
		f->addFlowLit(max_flow, v,inclusive);
	}

	void minConnectedComponents(int min_components, Var v) {
		if (!component_detector) {
			component_detector = new ConnectedComponentsDetector<Weight>(detectors.size(), this, g_under, g_over,
					drand(rnd_seed));
			detectors.push(component_detector);
		}
		component_detector->addConnectedComponentsLit(v, min_components);
	}
	void acyclic(Var v,bool directed) {
		if (!cycle_detector) {
			cycle_detector = new CycleDetector<Weight>(detectors.size(), this, g_under, g_over, true, drand(rnd_seed));
			detectors.push(cycle_detector);
		}
		cycle_detector->addAcyclicLit(directed, v);
	}
	
	void steinerTree(const vec<std::pair<int, Var> > & terminals, int steinerTreeID) {
		steiner_detectors.growTo(steinerTreeID + 1);
		assert(!steiner_detectors[steinerTreeID]);
		steiner_detectors[steinerTreeID] = new SteinerDetector<Weight>(detectors.size(), this,  g_under, g_over,
				drand(rnd_seed));
		detectors.push(steiner_detectors[steinerTreeID]);
		for (int i = 0; i < terminals.size(); i++) {
			steiner_detectors[steinerTreeID]->addTerminalNode(terminals[i].first, terminals[i].second);
		}
	}
	
	void addSteinerWeightConstraint(int steinerTreeID, Weight weight, Var outerVar) {
		if (steinerTreeID >= steiner_detectors.size()) {
			fprintf(stderr, "invalid steinerTreeID %d\n", steinerTreeID);
			exit(1);
		}
		steiner_detectors[steinerTreeID]->addWeightLit(weight, outerVar);
	}
	
	/*	void inTerminalSet(int node, int terminalSet, Var outerVar){
	 terminalSets.growTo(terminalSet+1);
	 while(terminalSets[terminalSet].nodes()<node){
	 terminalSets[terminalSet].addNode();
	 }
	 Var v= newVar(outerVar,node,true);
	 Lit l = mkLit(v,false);
	 if(terminalSets[terminalSet].getNodeVar(node)<0){
	 terminalSets[terminalSet].setNodeVar(v);
	 }
	 }*/

	void printSolution() {
		
		for (auto * d : detectors) {
			assert(d);
			d->printSolution();
		}
	}
	
};

}
;

#endif /* DGRAPH_H_ */
