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
#include "dgl/Reach.h"
#include "dgl/Dijkstra.h"
#include "dgl/BFS.h"
#include "dgl/DFS.h"
#include "dgl/UnweightedDistance.h"
#include "core/SolverTypes.h"
#include "mtl/Map.h"
#include "dgl/MaxFlow.h"
#include "dgl/IBFS.h"
#include "dgl/EdmondsKarp.h"
#include "dgl/EdmondsKarpAdj.h"
#include "dgl/Chokepoint.h"
#include "dgl/WeightedDijkstra.h"

#include "utils/System.h"
#include "Detector.h"
namespace Minisat{
class GraphTheorySolver;
class ConnectDetector:public Detector{
public:
		GraphTheorySolver * outer;
		 DynamicGraph &g;
		 DynamicGraph &antig;
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


		std::vector<ForceReason> forced_edges;
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
		WeightedDijkstra<vec<double>> * rnd_path;
		vec<double> rnd_weight;
		struct OptimalWeightEdgeStatus{
			ConnectDetector & detector;
			int operator [] (int edge) const ;
			int size()const;
			OptimalWeightEdgeStatus(ConnectDetector & _outer):detector(_outer){}

		};
		OptimalWeightEdgeStatus opt_weight;
		WeightedDijkstra<OptimalWeightEdgeStatus> * opt_path;
		bool check_positive;
		bool check_negative;

		struct ChokepointStatus{
			ConnectDetector & detector;
			bool mustReach(int node);
			bool operator() (int edge_id);
			ChokepointStatus(ConnectDetector & _outer):detector(_outer){

			}
		}chokepoint_status;

		Chokepoint<ChokepointStatus> chokepoint;

		int getNode(Var reachVar){
			assert(reachVar>=first_reach_var);
			int index = reachVar-first_reach_var;
			assert(index< reach_lit_map.size());
			assert(reach_lit_map[index]>=0);
			return reach_lit_map[index];
		}

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

		ConnectDetector(int _detectorID, GraphTheorySolver * _outer, DynamicGraph &_g, DynamicGraph &_antig, int _source,double seed=1);//:Detector(_detectorID),outer(_outer),within(-1),source(_source),rnd_seed(seed),positive_reach_detector(NULL),negative_reach_detector(NULL),positive_path_detector(NULL),positiveReachStatus(NULL),negativeReachStatus(NULL),chokepoint_status(*this),chokepoint(chokepoint_status, _antig,source){}
		virtual ~ConnectDetector(){

		}
		const char* getName(){
			return "Connected Detector";
		}
};
};
#endif /* REACHDETECTOR_H_ */
