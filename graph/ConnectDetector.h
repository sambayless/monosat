/*
 * ReachDetector.h
 *
 *  Created on: 2014-01-08
 *      Author: sam
 */

#ifndef CONNECTDETECTOR_H_
#define CONNECTDETECTOR_H_
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
#include "EdmondsKarp.h"
#include "EdmondsKarpAdj.h"
#include "Chokepoint.h"
#include "WeightedDijkstra.h"

#include "utils/System.h"
#include "Detector.h"
namespace Minisat{
class GraphTheorySolver;
class ConnectDetector:public Detector{
public:
		GraphTheorySolver * outer;
		 DynamicGraph<PositiveEdgeStatus> &g;
		 DynamicGraph<NegativeEdgeStatus> &antig;
		int within;
		int source;
		int constraintsBuilt;
		double rnd_seed;
		CRef reach_marker;
		CRef non_reach_marker;
		CRef forced_reach_marker;

		Reach * positive_reach_detector;
		Reach * negative_reach_detector;
		Reach *  positive_path_detector;

		vec<vec<Lit>>  dist_lits;

		vec<Lit>  reach_lits;
		Var first_reach_var;
		vec<int> reach_lit_map;
		vec<int> force_reason;


		vec<ForceReason> forced_edges;
		struct Change{
				Lit l;
				int u;
			};
		vec<Change> changed;

		vec<Change> & getChanged(){
			return changed;
		}

		//stats

		int stats_full_updates;
		int stats_fast_updates;
		int stats_fast_failed_updates;
		int stats_skip_deletes;
		int stats_skipped_updates;
		int stats_num_skipable_deletions;
		double mod_percentage;

		double stats_full_update_time;
		double stats_fast_update_time;


		struct ReachStatus{
			ConnectDetector & detector;
			bool polarity;
			void setReachable(int u, bool reachable);
			bool isReachable(int u) const{
				return false;
			}

			void setMininumDistance(int u, bool reachable, int distance);


			ReachStatus(ConnectDetector & _outer, bool _polarity):detector(_outer), polarity(_polarity){}
		};
		ReachStatus *positiveReachStatus;
		ReachStatus *negativeReachStatus;
		WeightedDijkstra<NegativeEdgeStatus,vec<double>> * rnd_path;
		vec<double> rnd_weight;
		struct OptimalWeightEdgeStatus{
			ConnectDetector & detector;
			int operator [] (int edge) const ;
			int size()const;
			OptimalWeightEdgeStatus(ConnectDetector & _outer):detector(_outer){}

		};
		OptimalWeightEdgeStatus opt_weight;
		WeightedDijkstra<NegativeEdgeStatus,OptimalWeightEdgeStatus> * opt_path;
		bool check_positive;
		bool check_negative;

		struct ChokepointStatus{
			ConnectDetector & detector;
			bool mustReach(int node);
			bool operator() (int edge_id);
			ChokepointStatus(ConnectDetector & _outer):detector(_outer){

			}
		}chokepoint_status;

		Chokepoint<ChokepointStatus ,NegativeEdgeStatus> chokepoint;

		int getNode(Var reachVar){
			assert(reachVar>=first_reach_var);
			int index = reachVar-first_reach_var;
			assert(index< reach_lit_map.size());
			assert(reach_lit_map[index]>=0);
			return reach_lit_map[index];
		}

	/*	Lit getLit(int node){

			return reach_lits[node];

		}*/
		void buildSATConstraints(int distance=-1);
		bool propagate(vec<Lit> & conflict);
		void buildReachReason(int node,vec<Lit> & conflict);
		void buildNonReachReason(int node,vec<Lit> & conflict);
		void buildForcedEdgeReason(int reach_node, int forced_edge_id,vec<Lit> & conflict);
		void buildReason(Lit p, vec<Lit> & reason, CRef marker);
		bool checkSatisfied();
		void addLit(int from, int to, Var reach_var);
		Lit decide();
		void preprocess();
		void dbg_sync_reachability();

		ConnectDetector(int _detectorID, GraphTheorySolver * _outer, DynamicGraph<PositiveEdgeStatus> &_g, DynamicGraph<NegativeEdgeStatus> &_antig, int _source,double seed=1);//:Detector(_detectorID),outer(_outer),within(-1),source(_source),rnd_seed(seed),positive_reach_detector(NULL),negative_reach_detector(NULL),positive_path_detector(NULL),positiveReachStatus(NULL),negativeReachStatus(NULL),chokepoint_status(*this),chokepoint(chokepoint_status, _antig,source){}
		virtual ~ConnectDetector(){

		}
};
};
#endif /* REACHDETECTOR_H_ */
