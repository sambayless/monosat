/*
 * ReachDetector.h
 *
 *  Created on: 2014-01-08
 *      Author: sam
 */

#ifndef REACHDETECTOR_H_
#define REACHDETECTOR_H_
#include "utils/System.h"
#include "GraphTheoryTypes.h"
#include "dgl/graph/DynamicGraph.h"
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
using namespace dgl;
namespace Minisat{
template<typename Weight>
class GraphTheorySolver;
template<typename Weight>
class ReachDetector:public Detector{
public:
		GraphTheorySolver<Weight> * outer;
		 DynamicGraph &g;
		 DynamicGraph &antig;
		 DynamicGraph cutgraph;
		int within;
		int source;
		double rnd_seed;
		int constraintsBuilt;
		CRef reach_marker;
		CRef non_reach_marker;
		CRef forced_reach_marker;

		Reach * positive_reach_detector;
		Reach * negative_reach_detector;
		Reach *  positive_path_detector;
		Reach *  cutgraph_reach_detector;
		vec<Lit>  reach_lits;
		Var first_reach_var;
		vec<int> order_vec;
		vec<int> reach_lit_map;
		vec<int> force_reason;
/*
		struct DistLit{
			Lit l;
			int min_distance;

		};*/
		vec<vec<Lit> > dist_lits;


		std::vector<ForceReason> forced_edges;
		struct Change{
				Lit l;
				int u;
			};
		vec<Change> changed;

		vec<Change> & getChanged(){
			return changed;
		}
		vec<Lit> extra_conflict;
		vec<int> removed_edges;
		//stats

		int stats_full_updates=0;
		int stats_fast_updates=0;
		int stats_fast_failed_updates=0;
		int stats_skip_deletes=0;
		int stats_skipped_updates=0;
		int stats_num_skipable_deletions=0;
		int stats_learnt_components = 0;
		int stats_learnt_components_sz =0;
		double mod_percentage=0.2;
		int stats_pure_skipped=0;
		int stats_shrink_removed=0;
		double stats_full_update_time=0;
		double stats_fast_update_time=0;

		void printStats(){
			//printf("Reach detector\n");
			Detector::printStats();
			if(opt_detect_pure_theory_lits)
				printf("\tPropagations skipped by pure literal detection: %d\n", stats_pure_skipped);
			if(opt_shrink_theory_conflicts){
				printf("\t%d lits removed by shrinking conflicts\n",stats_shrink_removed);
			}
			if(opt_learn_unreachable_component){
				printf("\t%d components learned, average component size: %f\n",stats_learnt_components,stats_learnt_components_sz / (float)stats_learnt_components);
			}
		}
		struct ReachStatus{
			ReachDetector & detector;
			bool polarity;
			void setReachable(int u, bool reachable);
			bool isReachable(int u) const{
				return false;
			}

			void setMininumDistance(int u, bool reachable, int distance);

			ReachStatus(ReachDetector & _outer, bool _polarity):detector(_outer), polarity(_polarity){}
		};
		ReachStatus *positiveReachStatus;
		ReachStatus *negativeReachStatus;

		WeightedDijkstra<double> * rnd_path;
		std::vector<double> rnd_weight;
		/*struct OptimalWeightEdgeStatus{
			ReachDetector & detector;
			int operator [] (int edge) const ;
			int size()const;
			OptimalWeightEdgeStatus(ReachDetector & _outer):detector(_outer){}

		};

		OptimalWeightEdgeStatus opt_weight;
		WeightedDijkstra<OptimalWeightEdgeStatus> * opt_path;*/
		Reach * chokepoint_detector;

		struct ChokepointStatus{
			ReachDetector & detector;
			bool mustReach(int node);
			bool operator() (int edge_id);
			ChokepointStatus(ReachDetector & _outer):detector(_outer){

			}
		}chokepoint_status;

		Chokepoint<ChokepointStatus > chokepoint;

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
		void buildSATConstraints(bool onlyUnderApprox=false,int within_steps=-1);
		bool propagate(vec<Lit> & conflict);
		void buildReachReason(int node,vec<Lit> & conflict);
		void buildNonReachReason(int node,vec<Lit> & conflict);
		void buildForcedEdgeReason(int reach_node, int forced_edge_id,vec<Lit> & conflict);
		void buildReason(Lit p, vec<Lit> & reason, CRef marker);
		bool checkSatisfied();
		void printSolution();
		void addLit(int from, int to, Var reach_var);
		Lit decide();
		void preprocess();
		void dbg_sync_reachability();

		ReachDetector(int _detectorID, GraphTheorySolver<Weight> * _outer, DynamicGraph &_g, DynamicGraph &_antig, int _source,double seed=1);//:Detector(_detectorID),outer(_outer),within(-1),source(_source),rnd_seed(seed),positive_reach_detector(NULL),negative_reach_detector(NULL),positive_path_detector(NULL),positiveReachStatus(NULL),negativeReachStatus(NULL),chokepoint_status(*this),chokepoint(chokepoint_status, _antig,source){}
		virtual ~ReachDetector(){

		}

		const char* getName(){
			return "Reachability Detector";
		}
};
};
#endif /* REACHDETECTOR_H_ */
