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

#ifndef ALLPAIRSDETECTOR_H_
#define ALLPAIRSDETECTOR_H_
#include "utils/System.h"
#include "GraphTheoryTypes.h"
#include "dgl/DynamicGraph.h"
#include "dgl/AllPairs.h"
#include <vector>
#include "core/SolverTypes.h"
#include "mtl/Map.h"
#include "utils/System.h"
#include "Detector.h"

//#define DEBUG_ALLPAIRS

using namespace dgl;
namespace Monosat {
template<typename Weight>
class GraphTheorySolver;
template<typename Weight>
class AllPairsDetector: public Detector {
public:
	GraphTheorySolver<Weight> * outer;
	DynamicGraph<Weight>  &g_under;
	DynamicGraph<Weight>  &g_over;
	int within = 0;

	double rnd_seed;
#ifdef DEBUG_ALLPAIRS
	AllPairs * dbg_positive_reach_detector;
	AllPairs * dbg_negative_reach_detector;
#endif
	
	CRef underprop_marker;
	CRef overprop_marker;

	AllPairs * underapprox_reach_detector=nullptr;
	AllPairs * overapprox_reach_detector=nullptr;
	AllPairs * underapprox_path_detector=nullptr;

	//vec<Lit>  reach_lits;
	Var first_reach_var;

	vec<int> force_reason;

	struct DistLit {
		Lit l;
		int min_distance;
		int source;
		int node;
	};
	vec<bool> installed_sources;
	vec<int> sources;
	vec<vec<vec<DistLit> >> dist_lits;
	vec<DistLit> reach_lit_map;

	std::vector<int> tmp_path;

	struct Change {
		Lit l;
		int u;
		int source;
	};
	vec<Change> changed;

	vec<Change> & getChanged() {
		return changed;
	}
	struct IgnoreStatus {
		
		void setReachable(int from, int u, bool reachable) {
			
		}
		bool isReachable(int from, int u) const {
			return false;
		}
		
		void setMininumDistance(int from, int u, bool reachable, int distance) {
			
		}
		IgnoreStatus() {
		}
	} ignoreStatus;
	struct ReachStatus {
		AllPairsDetector & detector;
		bool polarity;
		void setReachable(int from, bool reachable) {
			//for compatability with reach algs
		}
		void setReachable(int from, int u, bool reachable);
		bool isReachable(int from, int u) const {
			return false;
		}
		
		void setMininumDistance(int from, int u, bool reachable, int distance);

		ReachStatus(AllPairsDetector & _outer, bool _polarity) :
				detector(_outer), polarity(_polarity) {
		}
	};
	ReachStatus *positiveReachStatus;
	ReachStatus *negativeReachStatus;

	int getNode(Var reachVar) {
		assert(reachVar >= first_reach_var);
		int index = reachVar - first_reach_var;
		assert(index < reach_lit_map.size());
		assert(reach_lit_map[index].node >= 0);
		return reach_lit_map[index].node;
	}
	int getSource(Var reachVar) {
		assert(reachVar >= first_reach_var);
		int index = reachVar - first_reach_var;
		assert(index < reach_lit_map.size());
		assert(reach_lit_map[index].source >= 0);
		return reach_lit_map[index].source;
	}
	/*	Lit getLit(int node){

	 return reach_lits[node];

	 }*/
	bool propagate(vec<Lit> & conflict);
	void buildReachReason(int from, int to, vec<Lit> & conflict);
	void buildNonReachReason(int from, int to, vec<Lit> & conflict);

	void buildReason(Lit p, vec<Lit> & reason, CRef marker);
	bool checkSatisfied();
	Lit decide();
	void addLit(int from, int to, Var reach_var, int within_steps = -1);
	AllPairsDetector(int _detectorID, GraphTheorySolver<Weight> * _outer, DynamicGraph<Weight>  &_g, DynamicGraph<Weight>  &_antig,
			double seed = 1); //:Detector(_detectorID),outer(_outer),within(-1),source(_source),rnd_seed(seed),positive_reach_detector(NULL),negative_reach_detector(NULL),positive_path_detector(NULL),positiveReachStatus(NULL),negativeReachStatus(NULL){}
	virtual ~AllPairsDetector() {
		if (underapprox_path_detector && underapprox_path_detector!=overapprox_reach_detector){
			delete underapprox_path_detector;
		}


		if(underapprox_reach_detector){
			delete underapprox_reach_detector;

		}
		if (overapprox_reach_detector){
			delete overapprox_reach_detector;

		}

#ifdef DEBUG_ALLPAIRS
	{
		delete dbg_positive_reach_detector;
		delete dbg_negative_reach_detector;
	}
#endif
	}
	const char* getName() {
		return "All-pairs Reachability Detector";
	}
};
}
;
#endif /* AllPairsDetector_H_ */
