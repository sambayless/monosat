/*
 * DistanceDetector.h
 *
 *  Created on: 2014-01-08
 *      Author: sam
 */

#ifndef MAXFLOWDETECTOR_H_
#define MAXFLOWDETECTOR_H_
#include "utils/System.h"

#include "Graph.h"
#include "MaxFlow.h"
#include "EdmondsKarp.h"
#include "core/SolverTypes.h"
#include "mtl/Map.h"


#include "utils/System.h"
#include "Detector.h"
namespace Minisat{
class GraphTheorySolver;
class MaxflowDetector:public Detector{
public:
		GraphTheorySolver * outer;
		DynamicGraph<PositiveEdgeStatus> & over_graph;
		 DynamicGraph<PositiveEdgeStatus> &g;
			 DynamicGraph<NegativeEdgeStatus> &antig;
		//int within;
		int source;
		int target;
		double rnd_seed;
		CRef reach_marker;
		CRef non_reach_marker;
		CRef forced_reach_marker;
		MaxFlow * positive_detector;
		MaxFlow * negative_detector;


		//vec<Lit>  reach_lits;
		Var first_reach_var;
		vec<int> reach_lit_map;
		vec<int> force_reason;

		struct DistLit{
			Lit l;
			int max_flow;

		};
		vec<DistLit> flow_lits;


		vec<MaxFlow::Edge> tmp_cut;
		vec<int> visit;
		vec<bool> seen;
/*		int getNode(Var reachVar){
			assert(reachVar>=first_reach_var);
			int index = reachVar-first_reach_var;
			assert(index< reach_lit_map.size());
			assert(reach_lit_map[index]>=0);
			return reach_lit_map[index];
		}*/

		bool propagate(vec<Assignment> & trail,vec<Lit> & conflict);
		void buildMaxFlowTooHighReason(int node,vec<Lit> & conflict);
		void buildMaxFlowTooLowReason(int node,vec<Lit> & conflict);
		void buildForcedEdgeReason(int reach_node, int forced_edge_id,vec<Lit> & conflict);
		void buildReason(Lit p, vec<Lit> & reason, CRef marker);
		bool checkSatisfied();
		Lit decide();
		void addFlowLit(int max_flow,Var reach_var);
		MaxflowDetector(int _detectorID, GraphTheorySolver * _outer, DynamicGraph<PositiveEdgeStatus> &_g, DynamicGraph<NegativeEdgeStatus> &_antig, int _source, int _target,double seed=1);//:Detector(_detectorID),outer(_outer),within(-1),source(_source),rnd_seed(seed),positive_reach_detector(NULL),negative_reach_detector(NULL),positive_path_detector(NULL),positiveReachStatus(NULL),negativeReachStatus(NULL){}
		virtual ~MaxflowDetector(){

		}
};
};
#endif /* DistanceDetector_H_ */
