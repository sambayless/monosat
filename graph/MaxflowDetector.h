/*
 * DistanceDetector.h
 *
 *  Created on: 2014-01-08
 *      Author: sam
 */

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
namespace Minisat{
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
/*		int getNode(Var reachVar){
			assert(reachVar>=first_reach_var);
			int index = reachVar-first_reach_var;
			assert(index< reach_lit_map.size());
			assert(reach_lit_map[index]>=0);
			return reach_lit_map[index];
		}*/

		bool propagate(vec<Lit> & conflict);
		void buildMaxFlowTooHighReason(Weight flow,vec<Lit> & conflict);
		void buildMaxFlowTooLowReason(Weight flow,vec<Lit> & conflict);
		void buildForcedEdgeReason(int reach_node, int forced_edge_id,vec<Lit> & conflict);
		void buildReason(Lit p, vec<Lit> & reason, CRef marker);
		bool checkSatisfied();
		Lit decide();
		void addFlowLit(Weight max_flow,Var reach_var);
		MaxflowDetector(int _detectorID, GraphTheorySolver<Weight> * _outer,std::vector<Weight> & capacities, DynamicGraph &_g, DynamicGraph &_antig, int _source, int _target,double seed=1);//:Detector(_detectorID),outer(_outer),within(-1),source(_source),rnd_seed(seed),positive_reach_detector(NULL),negative_reach_detector(NULL),positive_path_detector(NULL),positiveReachStatus(NULL),negativeReachStatus(NULL){}
		virtual ~MaxflowDetector(){

		}
		const char* getName(){
			return "Max-flow Detector";
		}
};
};


#endif /* DistanceDetector_H_ */
