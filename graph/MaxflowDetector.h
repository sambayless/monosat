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
		//int within;
		int source;
		double rnd_seed;

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
		vec<vec<DistLit> > flow_lits;
		struct Change{
				Lit l;
				int u;
			};
		vec<Change> changed;

		vec<MaxFlow::Edge> tmp_cut;
		vec<int> visit;
		vec<bool> seen;
		int getNode(Var reachVar){
			assert(reachVar>=first_reach_var);
			int index = reachVar-first_reach_var;
			assert(index< reach_lit_map.size());
			assert(reach_lit_map[index]>=0);
			return reach_lit_map[index];
		}

		bool propagate(vec<Assignment> & trail,vec<Lit> & conflict);
		void buildReachReason(int node,vec<Lit> & conflict);
		void buildNonReachReason(int node,vec<Lit> & conflict);
		void buildForcedEdgeReason(int reach_node, int forced_edge_id,vec<Lit> & conflict);
		void buildReason(Lit p, vec<Lit> & reason, CRef marker);
		bool checkSatisfied();
		Lit decide();
		void addLit(int from, int to, Var reach_var,int within_steps=-1);
		MaxflowDetector(int _detectorID, GraphTheorySolver * _outer, DynamicGraph<PositiveEdgeStatus> &_g, DynamicGraph<NegativeEdgeStatus> &_antig, int _source,double seed=1);//:Detector(_detectorID),outer(_outer),within(-1),source(_source),rnd_seed(seed),positive_reach_detector(NULL),negative_reach_detector(NULL),positive_path_detector(NULL),positiveReachStatus(NULL),negativeReachStatus(NULL){}
		virtual ~MaxflowDetector(){

		}
};
};
#endif /* DistanceDetector_H_ */
