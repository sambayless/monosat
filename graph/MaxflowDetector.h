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
#ifndef MAXFLOWDETECTOR_H_
#define MAXFLOWDETECTOR_H_
#include "utils/System.h"
#include "dgl/KohliTorr.h"
#include "GraphTheoryTypes.h"
#include "dgl/DynamicGraph.h"
#include "dgl/MaxFlow.h"
#include "dgl/EdmondsKarp.h"
#include "core/SolverTypes.h"
#include "mtl/Map.h"
#include "mtl/Deque.h"
#include "utils/System.h"
#include "Detector.h"
#include "bv/BVTheorySolver.h"
using namespace dgl;
namespace Monosat {
template<typename Weight>
class GraphTheorySolver;
template<typename Weight = int>
class MaxflowDetector: public Detector, public DynamicGraphAlgorithm  {
public:
	GraphTheorySolver<Weight> * outer;
	std::vector<Weight> capacities;

	DynamicGraph<Weight> &g_under;
	DynamicGraph<Weight> &g_over;

	int source;
	int target;
	double rnd_seed;
	CRef underprop_marker;
	CRef overprop_marker;

	MaxFlow<Weight>* underapprox_detector = nullptr;
	MaxFlow<Weight> * overapprox_detector = nullptr;
	MaxFlow<Weight> * underapprox_conflict_detector = nullptr;
	MaxFlow<Weight> * overapprox_conflict_detector = nullptr;
	int last_decision_status = -1;
	int last_decision_q_pos = 0;
	int alg_id=-1;

	long stats_decision_calculations = 0;
	double stats_flow_calc_time = 0;
	double stats_flow_recalc_time = 0;
	double stats_redecide_time = 0;

	Lit last_decision_lit = lit_Undef;

	//vec<Lit> decisions;
	vec<bool> is_potential_decision;
	vec<bool> in_decision_q;
	vec<int> potential_decisions;
	struct DecisionS{
		int path_decision:1;
		int edgeID:31;
		DecisionS():path_decision(0), edgeID(-1){}
		DecisionS(const int edgeID):path_decision(0), edgeID(edgeID){}
		DecisionS(int edgeID, bool pathDecision):path_decision(pathDecision), edgeID(edgeID){}
		operator int() const { return edgeID; }
	};
	Deque<DecisionS> potential_decisions_q;
	vec<Lit> to_decide;
	std::vector<int> q;

	DynamicGraph<long> learn_graph;
	vec<int> back_edges;
	int learngraph_history_qhead = 0;
	int learngraph_history_clears = -1;
	MaxFlow<long> * learn_cut = nullptr;
	//int current_decision_edge=-1;
	//vec<Lit>  reach_lits;
	Var first_reach_var;
	vec<int> reach_lit_map;
	vec<int> force_reason;

	struct DistLit {
		Lit l;
		Weight max_flow=-1;
		BitVector<Weight>  bv;
		bool inclusive;//If inclusive true, l is true iff the maximum flow is >= max_flow; else, l is true iff the maximum flow is > max_flow.
	};
	vec<DistLit> flow_lits;

	std::vector<MaxFlowEdge> cut;

	vec<MaxFlowEdge> tmp_cut;
	vec<int> visit;
	vec<bool> seen;
	vec<bool> seen_path;


	void backtrack(int level) {
		to_decide.clear();
		last_decision_status = -1;
		//LevelDetector::backtrack(level);
	}
	void collectChangedEdges();
	void collectDisabledEdges();
	void updateHistory(){
		collectDisabledEdges();
	}
	bool propagate(vec<Lit> & conflict){
		return propagate(conflict,false);
	}
	bool propagate(vec<Lit> & conflict, bool backtrackOnly);
	void buildMaxFlowTooHighReason(Weight flow, vec<Lit> & conflict);
	void buildMaxFlowTooLowReason(Weight flow, vec<Lit> & conflict, bool force_maxflow = false);
	void buildForcedEdgeReason(int reach_node, int forced_edge_id, vec<Lit> & conflict);
	void buildReason(Lit p, vec<Lit> & reason, CRef marker);
	bool checkSatisfied();
	void undecide(Lit l);
	Lit decide();
	void printStats() {
		Detector::printStats();
		if (mincutalg == MinCutAlg::ALG_KOHLI_TORR) {
			KohliTorr<Weight> * kt = (KohliTorr<Weight> *) overapprox_detector;
			printf(
					"\tDecision flow calculations: %ld, (redecide: %f s) flow_calc %f s, flow_discovery %f s, (%ld) (maxflow %f,flow assignment %f)\n",
					stats_decision_calculations, stats_redecide_time, stats_flow_calc_time, stats_flow_recalc_time,
					kt->stats_flow_calcs, kt->stats_flow_time, kt->stats_calc_time);
		} else
			printf("\tDecision flow calculations: %ld\n", stats_decision_calculations);
		
	}
	//Lit decideByPath(int level);
	void dbg_decisions();
	void printSolution(std::ostream & write_to);
	void addFlowLit(Weight max_flow, Var reach_var, bool inclusive);
	void addFlowBVLessThan(const BitVector<Weight>  &bv, Var v, bool inclusive);
	MaxflowDetector(int _detectorID, GraphTheorySolver<Weight> * _outer,
			DynamicGraph<Weight>  &_g, DynamicGraph<Weight>  &_antig, int _source, int _target, double seed = 1); //:Detector(_detectorID),outer(_outer),within(-1),source(_source),rnd_seed(seed),positive_reach_detector(NULL),negative_reach_detector(NULL),positive_path_detector(NULL),positiveReachStatus(NULL),negativeReachStatus(NULL){}
	~MaxflowDetector() {
		
	}
	const char* getName() {
		return "Max-flow Detector";
	}
	
/*	void decideEdge(int edgeID,  bool assign = true) {
		assert(decisions.size() >= decisionLevel());

		newDecisionLevel(outer->decisionLevel()+1);
		
		Lit l = mkLit(outer->getEdgeVar(edgeID), !assign);

		assert(!decisions.contains(l));
		assert(!decisions.contains(~l));
		
		//decisions.push(l);
		
		assert(decisions.size() == decisionLevel());
		dbg_decisions();
	}*/
	/*void localBacktrack() {
		dbg_decisions();
		//undo a decision edge, and return it to the set of potential decisions
		assert(decisions.size() == decisionLevel() + 1);
		
		Lit l = decisions.last();
		decisions.pop();
		
		assert(outer->isEdgeVar(var(l)));
		int edgeID = outer->getEdgeID(var(l));

		//assert(!is_potential_decision[edgeID]);
		assert(!potential_decisions.contains(edgeID));
		assert(!potential_decisions_q.contains(edgeID));

		//this check is optional, but if used, must use the _conflict_ detector (to match the decisions, which are also using the conflict detector).
		if (overapprox_conflict_detector->getEdgeFlow(edgeID) > 0) {			//this check is optional

			is_potential_decision[edgeID] = true;
			if (opt_maxflow_decisions_q == 0)
				potential_decisions_q.insertBack(edgeID);
			else if (opt_maxflow_decisions_q == 1) {
				potential_decisions_q.insertBack(edgeID);//insert in LIFO order, no FIFO, because we are unwinding the decisions
			} else if (opt_maxflow_decisions_q == 2) {
				potential_decisions_q.insert(edgeID);			//insert in FIFO order instead
			} else if (opt_maxflow_decisions_q == 3) {
				potential_decisions_q.insertBack(edgeID);//insert in LIFO order, no FIFO, because we are unwinding the decisions
			} else if (opt_maxflow_decisions_q == 4) {
				potential_decisions_q.insert(edgeID);//insert in LIFO order, no FIFO, because we are unwinding the decisions
			}
		} else {
			is_potential_decision[edgeID] = false;			//discard this edge from the set of potential decisions
		}
		dbg_decisions();
	}*/

	void preprocess(){

		is_potential_decision.growTo(g_under.edges(), false);

		in_decision_q.growTo(g_under.edges(), false);
	}
	
private:
	void buildDinitzLinkCut();
};
template<>
void MaxflowDetector<int>::buildDinitzLinkCut();

}
;

#endif /* DistanceDetector_H_ */
