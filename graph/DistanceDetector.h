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

#ifndef DISTANCEETECTOR_H_
#define DISTANCEETECTOR_H_
#include "utils/System.h"

#include "GraphTheoryTypes.h"
#include "dgl/DynamicGraph.h"
#include "dgl/Reach.h"
#include "dgl/Distance.h"
#include "dgl/Dijkstra.h"
#include "dgl/BFS.h"
#include "dgl/MaxFlow.h"
#include "core/SolverTypes.h"
#include "mtl/Map.h"
#include "WeightedDijkstra.h"
#include <gmpxx.h>
#include "utils/System.h"
#include "Detector.h"
#include "bv/BVTheorySolver.h"
#include <vector>
using namespace dgl;
namespace Monosat {
template<typename Weight>
class GraphTheorySolver;

template<typename Weight>
class DistanceDetector: public Detector {
public:
	GraphTheorySolver<Weight> * outer;

	DynamicGraph<Weight>  &g_under;
	DynamicGraph<Weight>  &g_over;
	//int within;
	int source;
	double rnd_seed;
	int constraintsBuilt;
	CRef unweighted_underprop_marker;
	CRef unweighted_overprop_marker;
	CRef weighted_underprop_marker;
	CRef weighted_overprop_marker;

	Distance<int> * underapprox_unweighted_distance_detector = nullptr;
	Distance<int> * overapprox_unweighted_distance_detector = nullptr;

	Distance<Weight> * underapprox_weighted_distance_detector = nullptr;
	Distance<Weight> * overapprox_weighted_distance_detector = nullptr;
	Distance<Weight> * underapprox_weighted_path_detector = nullptr;
	Reach * underapprox_path_detector = nullptr;

	//vec<Lit>  reach_lits;
	Var first_reach_var;
	struct ReachLit{
		int to;
		int within;
	};
	vec<ReachLit> reach_lit_map;
	vec<int> force_reason;
	bool has_unweighted_shortest_paths_overapprox = false;
	bool has_weighted_shortest_paths_overapprox = false;
	vec<int> unweighted_over_approx_shortest_paths;
	vec<Weight> over_approx_shortest_paths;
	MaxFlow<long> * conflict_flow = nullptr;

	int max_unweighted_distance;

	long stats_pure_skipped = 0;
	long stats_distance_gt_reasons = 0;
	long stats_distance_leq_reasons = 0;
	long stats_unweighted_gt_reasons = 0;
	long stats_unweighted_leq_reasons = 0;
	long stats_gt_unweighted_edges_skipped = 0;
	long stats_gt_weighted_edges_skipped = 0;

	vec<vec<Lit> > unweighted_sat_lits;

	struct UnweightedDistLit {
		Lit l;
		int min_unweighted_distance;
		
	};

	vec<vec<UnweightedDistLit> > unweighted_dist_lits;

	struct WeightedDistLit {
		Lit l;
		int u;
		Weight min_distance;
		bool strictComparison;
		bool operator <(WeightedDistLit & b) const {
			return min_distance < b.min_distance;
		}
	};
	vec<WeightedDistLit> weighted_dist_lits;

	struct WeightedDistBVLit {
		Lit l;
		int u;
		BitVector<Weight> bv;
		bool strictComparison;
	};
	vec<WeightedDistBVLit> weighted_dist_bv_lits;

	struct Change {
		//Var v;
		int u;
		//int min_distance;
	};
	vec<Change> changed;
	vec<bool> is_changed;
	vec<Var> tmp_nodes;

	std::vector<double> rnd_weight;

	WeightedDijkstra<Weight,double> * rnd_path;
	/*struct OptimalWeightEdgeStatus{
	 DistanceDetector & detector;
	 int operator [] (int edge) const ;
	 int size() const;
	 OptimalWeightEdgeStatus(DistanceDetector & _outer):detector(_outer){}

	 };*/
	//OptimalWeightEdgeStatus opt_weight;
	//WeightedDijkstra<OptimalWeightEdgeStatus> * opt_path;
	struct ReachStatus {
		DistanceDetector & detector;
		bool polarity;
		void setReachable(int u, bool reachable);
		bool isReachable(int u) const {
			return false;
		}
		
		void setMininumDistance(int u, bool reachable, int distance);

		ReachStatus(DistanceDetector & _outer, bool _polarity) :
				detector(_outer), polarity(_polarity) {
		}
	};
	struct DistanceStatus {
		DistanceDetector & detector;
		bool polarity;
		void setReachable(int u, bool reachable);
		bool isReachable(int u) const {
			return false;
		}
		
		void setMininumDistance(int u, bool reachable, Weight& distance);

		DistanceStatus(DistanceDetector & _outer, bool _polarity) :
				detector(_outer), polarity(_polarity) {
		}
	};

	/*struct CutStatus {
		long one = 1;
		long inf = 0xFFFF;
		DistanceDetector & outer;

		const long &operator [](int id) const {
			if (id % 2 == 0) {
				return one;
			} else {
				return inf;
			}
		}
		int size() const {
			return outer.g_under.edges() * 2;
		}
		CutStatus(DistanceDetector & _outer) :
				outer(_outer) {
		}
		
	} cutStatus;*/
	std::vector<MaxFlowEdge> cut;

	ReachStatus *positiveReachStatus;
	ReachStatus *negativeReachStatus;
	DistanceStatus *positiveDistanceStatus;
	DistanceStatus *negativeDistanceStatus;
	int getNode(Var reachVar) {
		assert(reachVar >= first_reach_var);
		int index = reachVar - first_reach_var;
		assert(index < reach_lit_map.size());
		assert(reach_lit_map[index].to >= 0);
		return reach_lit_map[index].to;
	}
	int getMaximumDistance(Var reachVar) {
			assert(reachVar >= first_reach_var);
			int index = reachVar - first_reach_var;
			assert(index < reach_lit_map.size());
			assert(reach_lit_map[index].within >= 0);
			return reach_lit_map[index].within;
		}
	void printStats() {
		//printf("Distance detector\n");
		Detector::printStats();
		if (opt_verb > 0) {
			if (opt_detect_pure_theory_lits)
				printf("\tPropagations skipped by pure literal detection: %ld\n", stats_pure_skipped);
			printf("\tUnweighted Reasons (leq,gt): %ld,%ld\n", stats_unweighted_leq_reasons, stats_unweighted_gt_reasons);
			printf("\tWeighted Reasons (leq,gt): %ld,%ld\n", stats_distance_leq_reasons, stats_distance_gt_reasons);
			printf("\tConflict Edges Skipped (unweighted %ld, weighted %ld)\n", stats_gt_unweighted_edges_skipped,
					stats_gt_weighted_edges_skipped);
		}
	}
	
	void printSolution(std::ostream& write_to);

	/*	Lit getLit(int node){

	 return reach_lits[node];

	 }*/

	void unassign(Lit l) {
		Detector::unassign(l);
		int index = var(l) - first_reach_var;
		
		if (index >= 0 && index < reach_lit_map.size() && reach_lit_map[index].to != -1) {
			int node = reach_lit_map[index].to;
			if (!is_changed[node]) {
				changed.push( { node });
				is_changed[node] = true;
			}
		}
	}
	void preprocess();
	bool propagate(vec<Lit> & conflict);
	void buildUnweightedDistanceLEQReason(int node, vec<Lit> & conflict);
	void buildUnweightedDistanceGTReason(int node, int within_steps, vec<Lit> & conflict);
	void buildDistanceLEQReason(int to, Weight & min_distance, vec<Lit> & conflict, bool strictComparison=false);
	void buildDistanceGTReason(int to, Weight & min_distance, vec<Lit> & conflict, bool strictComparison=true);

	void buildReason(Lit p, vec<Lit> & reason, CRef marker);
	bool checkSatisfied();
	Lit decide();
	void updateShortestPaths(bool unweighted);
	void addUnweightedShortestPathLit(int from, int to, Var reach_var, int within_steps = -1);
	void addWeightedShortestPathLit(int from, int to, Var reach_var, Weight within_distance, bool strictComparison);
	void addWeightedShortestPathBVLit(int from, int to, Var reach_var, const BitVector<Weight> & bv, bool strictComparison);

	DistanceDetector(int _detectorID, GraphTheorySolver<Weight> * _outer,
			DynamicGraph<Weight>  &_g, DynamicGraph<Weight>  &_antig, int _source, double seed = 1);//:Detector(_detectorID),outer(_outer),within(-1),source(_source),rnd_seed(seed),positive_reach_detector(NULL),negative_reach_detector(NULL),positive_path_detector(NULL),positiveReachStatus(NULL),negativeReachStatus(NULL){}
	virtual ~DistanceDetector() {
		if (positiveDistanceStatus )
			delete positiveDistanceStatus;
		if (negativeDistanceStatus )
			delete negativeDistanceStatus;
		if (overapprox_unweighted_distance_detector)
			delete overapprox_unweighted_distance_detector;

		if (underapprox_path_detector && underapprox_path_detector != underapprox_unweighted_distance_detector)
			delete underapprox_path_detector;

		if (underapprox_weighted_path_detector && underapprox_weighted_path_detector != underapprox_weighted_distance_detector)
			delete underapprox_weighted_path_detector;

		if(positiveReachStatus)
			delete positiveReachStatus;
		if (negativeReachStatus)
			delete negativeReachStatus;

		if (underapprox_unweighted_distance_detector)
			delete underapprox_unweighted_distance_detector;



		if (underapprox_weighted_distance_detector)
			delete underapprox_weighted_distance_detector;

		if (overapprox_weighted_distance_detector)
			delete overapprox_weighted_distance_detector;



		if (rnd_path)
			delete rnd_path;

		if (conflict_flow)
			delete conflict_flow;
	}
	const char* getName() {
		return "Shortest Path Detector";
	}
private:
	void buildUnweightedSATConstraints(bool onlyUnderApprox,int distance = -1);
};

}
;
#endif /* DistanceDetector_H_ */
