/*
 * TreeReachDetector.h
 *
 *  Created on: 2014-01-08
 *      Author: sam
 */

#ifndef TREEREACHDETECTOR_H_
#define TREEREACHDETECTOR_H_
#include "utils/System.h"

#include "Graph.h"
#include "Reach.h"
#include "Dijkstra.h"
#include "BFSReachability.h"
#include "DFSReachability.h"
#include "UnweightedDistance.h"
#include "core/SolverTypes.h"
#include "mtl/Map.h"
#include "MaxFlow.h"
#include "IBFS.h"
#include "LinkCutForest.h"
#include "EdmondsKarp.h"
#include "EdmondsKarpAdj.h"
#include "Chokepoint.h"
#include "WeightedDijkstra.h"
#include "utils/System.h"
#include "DynamicConnectivity.h"
#include "Detector.h"
namespace Minisat{
class GraphTheorySolver;
class TreeReachDetector:public Detector{
public:
		GraphTheorySolver * outer;
		 DynamicGraph<PositiveEdgeStatus> &g;
		 DynamicGraph<NegativeEdgeStatus> &antig;
		int within;
		int source;
		double rnd_seed;
		CRef reach_marker;
		CRef non_reach_marker;


		LinkCutForest<PositiveEdgeStatus>  * positive_reach_detector;

		DynamicConnectivity<NegativeEdgeStatus> * negative_reach_detector;
		//Reach *  positive_path_detector;

		vec<Lit>  reach_lits;
		Var first_reach_var;
		vec<int> reach_lit_map;

		struct Change{
				Lit l;
				int u;
			};
		vec<Change> changed;

		vec<Change> & getChanged(){
			return changed;
		}

		WeightedDijkstra<NegativeEdgeStatus,vec<double>> * rnd_path;
		vec<double> rnd_weight;
		struct OptimalWeightEdgeStatus{
			TreeReachDetector & detector;
			int operator [] (int edge) const ;
			int size()const;
			OptimalWeightEdgeStatus(TreeReachDetector & _outer):detector(_outer){}

		};
		OptimalWeightEdgeStatus opt_weight;
		WeightedDijkstra<NegativeEdgeStatus,OptimalWeightEdgeStatus> * opt_path;
		bool check_positive;
		bool check_negative;
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
		void addLit(int from, int to, Var reach_var);
		Lit decide();
		void preprocess();
		void dbg_sync_reachability();

		TreeReachDetector(int _detectorID, GraphTheorySolver * _outer, DynamicGraph<PositiveEdgeStatus> &_g, DynamicGraph<NegativeEdgeStatus> &_antig, int _source,double seed=1);//:Detector(_detectorID),outer(_outer),within(-1),source(_source),rnd_seed(seed),positive_reach_detector(NULL),negative_reach_detector(NULL),positive_path_detector(NULL),positiveReachStatus(NULL),negativeReachStatus(NULL),chokepoint_status(*this),chokepoint(chokepoint_status, _antig,source){}
		virtual ~TreeReachDetector(){

		}
};
};
#endif /* REACHDETECTOR_H_ */
