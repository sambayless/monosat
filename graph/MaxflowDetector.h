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

#include "GraphTheoryTypes.h"
#include "dgl/graph/DynamicGraph.h"
#include "dgl/MaxFlow.h"
#include "dgl/EdmondsKarp.h"
#include "core/SolverTypes.h"
#include "mtl/Map.h"

#include "utils/System.h"
#include "Detector.h"
using namespace dgl;
namespace Monosat{
template<typename Weight>
class GraphTheorySolver;
template<typename Weight=int>
class MaxflowDetector:public Detector{
public:
		GraphTheorySolver<Weight> * outer;
		std::vector<Weight> capacities;
		DynamicGraph & over_graph;
		 DynamicGraph &g;
		 DynamicGraph &antig;
		//int within;
		int source;
		int target;
		double rnd_seed;
		CRef reach_marker;
		CRef non_reach_marker;
		CRef forced_reach_marker;
		MaxFlow<Weight>* positive_detector;
		MaxFlow<Weight> * negative_detector;
		MaxFlow<Weight> * positive_conflict_detector;
		MaxFlow<Weight> * negative_conflict_detector;
		int last_decision_status=-1;
		Lit last_decision_lit = lit_Undef;
		vec<Lit> to_decide;
		std::vector<int> q;

		DynamicGraph learn_graph;
		vec<int> back_edges;
		std::vector<int> learn_caps;
		MaxFlow<int> * learn_cut=nullptr;

		//vec<Lit>  reach_lits;
		Var first_reach_var;
		vec<int> reach_lit_map;
		vec<int> force_reason;

		struct DistLit{
			Lit l;
			Weight max_flow;

		};
		vec<DistLit> flow_lits;


		vec<MaxFlowEdge> tmp_cut;
		vec<int> visit;
		vec<bool> seen;
		vec<int> prev;
		vec<int> dist;
/*		int getNode(Var reachVar){
			assert(reachVar>=first_reach_var);
			int index = reachVar-first_reach_var;
			assert(index< reach_lit_map.size());
			assert(reach_lit_map[index]>=0);
			return reach_lit_map[index];
		}*/
		void backtrack(int level){
			to_decide.clear();
			last_decision_status=-1;
		}
		bool propagate(vec<Lit> & conflict);
		void buildMaxFlowTooHighReason(Weight flow,vec<Lit> & conflict);
		void buildMaxFlowTooLowReason(Weight flow,vec<Lit> & conflict);
		void buildForcedEdgeReason(int reach_node, int forced_edge_id,vec<Lit> & conflict);
		void buildReason(Lit p, vec<Lit> & reason, CRef marker);
		bool checkSatisfied();
		Lit decide();
		void printSolution();
		void addFlowLit(Weight max_flow,Var reach_var);
		MaxflowDetector(int _detectorID, GraphTheorySolver<Weight> * _outer,std::vector<Weight> & capacities, DynamicGraph &_g, DynamicGraph &_antig, int _source, int _target,double seed=1);//:Detector(_detectorID),outer(_outer),within(-1),source(_source),rnd_seed(seed),positive_reach_detector(NULL),negative_reach_detector(NULL),positive_path_detector(NULL),positiveReachStatus(NULL),negativeReachStatus(NULL){}
		virtual ~MaxflowDetector(){

		}
		const char* getName(){
			return "Max-flow Detector";
		}
private:
		void buildDinitzLinkCut();
};
template<>
void MaxflowDetector<int>::buildDinitzLinkCut();

};


#endif /* DistanceDetector_H_ */
