/*
 * DistanceDetector.h
 *
 *  Created on: 2014-01-08
 *      Author: sam
 */

#ifndef DISTANCEETECTOR_H_
#define DISTANCEETECTOR_H_
#include "utils/System.h"

#include "GraphTheoryTypes.h"
#include "dgl/graph/DynamicGraph.h"
#include "dgl/Reach.h"
#include "dgl/Distance.h"
#include "dgl/Dijkstra.h"
#include "dgl/BFS.h"
#include "dgl/UnweightedDistance.h"
#include "core/SolverTypes.h"
#include "mtl/Map.h"
#include "dgl/WeightedDijkstra.h"
#include <gmpxx.h>
#include "utils/System.h"
#include "Detector.h"
#include <vector>
using namespace dgl;
namespace Minisat{
template<typename Weight>
class GraphTheorySolver;



template<typename Weight>
class DistanceDetector:public Detector{
public:
		GraphTheorySolver<Weight> * outer;
		 DynamicGraph &g;
		 DynamicGraph &antig;
		 std::vector<Weight> & weights;
		//int within;
		int source;
		double rnd_seed;
		int constraintsBuilt;
		CRef reach_marker;
		CRef non_reach_marker;
		CRef weighted_reach_marker;
		CRef weighted_non_reach_marker;
		CRef forced_reach_marker;
		Distance<int> * positive_reach_detector;
		Distance<int> * negative_reach_detector;

		Distance<Weight> * positive_weighted_distance_detector;
		Distance<Weight> * negative_weighted_distance_detector;

		Reach *  positive_path_detector;

		//vec<Lit>  reach_lits;
		Var first_reach_var;
		vec<int> reach_lit_map;
		vec<int> force_reason;
		int max_unweighted_distance;

		int stats_pure_skipped;
		vec<vec<Lit> > unweighted_sat_lits;

		struct UnweightedDistLit{
			Lit l;
			int min_unweighted_distance;

		};
		vec<vec<UnweightedDistLit> > unweighted_dist_lits;


		struct WeightedDistLit{
			Lit l;
			int u;
			Weight min_distance;
			bool operator <(WeightedDistLit & b)const{
				return min_distance<b.min_distance;
			}
		};
		vec<WeightedDistLit> weighted_dist_lits;

		struct Change{
				Lit l;
				int u;
			};
		vec<Change> changed;
		vec<Var> tmp_nodes;
		vec<Change> & getChanged(){
			return changed;
		}
		std::vector<double> rnd_weight;

		WeightedDijkstra<double> * rnd_path;
		/*struct OptimalWeightEdgeStatus{
			DistanceDetector & detector;
			int operator [] (int edge) const ;
			int size() const;
			OptimalWeightEdgeStatus(DistanceDetector & _outer):detector(_outer){}

		};*/
		//OptimalWeightEdgeStatus opt_weight;
		//WeightedDijkstra<OptimalWeightEdgeStatus> * opt_path;
		struct ReachStatus{
			DistanceDetector & detector;
			bool polarity;
			void setReachable(int u, bool reachable);
			bool isReachable(int u) const{
				return false;
			}

			void setMininumDistance(int u, bool reachable, int distance);

			ReachStatus(DistanceDetector & _outer, bool _polarity):detector(_outer), polarity(_polarity){}
		};
		struct DistanceStatus{
			DistanceDetector & detector;
			bool polarity;
			void setReachable(int u, bool reachable);
			bool isReachable(int u) const{
				return false;
			}

			void setMininumDistance(int u, bool reachable, Weight& distance);

			DistanceStatus(DistanceDetector & _outer, bool _polarity):detector(_outer), polarity(_polarity){}
		};
		ReachStatus *positiveReachStatus;
		ReachStatus *negativeReachStatus;
		DistanceStatus *positiveDistanceStatus;
		DistanceStatus *negativeDistanceStatus;
		int getNode(Var reachVar){
			assert(reachVar>=first_reach_var);
			int index = reachVar-first_reach_var;
			assert(index< reach_lit_map.size());
			assert(reach_lit_map[index]>=0);
			return reach_lit_map[index];
		}

		void printStats(){
			//printf("Distance detector\n");
			Detector::printStats();
			if(opt_detect_pure_theory_lits)
				printf("\tPropagations skipped by pure literal detection: %d\n", stats_pure_skipped);
		}

	/*	Lit getLit(int node){

			return reach_lits[node];

		}*/

		bool propagate(vec<Lit> & conflict);
		void buildReachReason(int node,vec<Lit> & conflict);
		void buildNonReachReason(int node,vec<Lit> & conflict);
		void buildDistanceLEQReason(int to,Weight & min_distance,vec<Lit> & conflict);
		void buildDistanceGTReason(int to,Weight & min_distance,vec<Lit> & conflict);
		void buildForcedEdgeReason(int reach_node, int forced_edge_id,vec<Lit> & conflict);
		void buildReason(Lit p, vec<Lit> & reason, CRef marker);
		bool checkSatisfied();
		Lit decide();
		void addUnweightedShortestPathLit(int from, int to, Var reach_var,int within_steps=-1);
		void addWeightedShortestPathLit(int from, int to, Var reach_var,Weight within_distance);
		DistanceDetector(int _detectorID, GraphTheorySolver<Weight> * _outer, std::vector<Weight> & weights, DynamicGraph &_g, DynamicGraph &_antig, int _source,double seed=1);//:Detector(_detectorID),outer(_outer),within(-1),source(_source),rnd_seed(seed),positive_reach_detector(NULL),negative_reach_detector(NULL),positive_path_detector(NULL),positiveReachStatus(NULL),negativeReachStatus(NULL){}
		virtual ~DistanceDetector(){

		}
		const char* getName(){
			return "Shortest Path Detector";
		}
private:
		void buildSATConstraints(int distance=-1);
};

};
#endif /* DistanceDetector_H_ */
