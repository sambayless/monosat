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
#ifndef MST_DETECTOR_H_
#define MST_DETECTOR_H_
#include "utils/System.h"

#include "GraphTheoryTypes.h"
#include "dgl/DynamicGraph.h"
#include "dgl/MinimumSpanningTree.h"

#include "core/SolverTypes.h"
#include "mtl/Map.h"

#include "dgl/alg/DisjointSets.h"
#include "utils/System.h"
#include "Detector.h"

using namespace dgl;
namespace Monosat {
template<typename Weight>
class GraphTheorySolver;
template<typename Weight = int>
class MSTDetector: public Detector {
public:
	GraphTheorySolver<Weight> * outer;

	DynamicGraph<Weight> & g_under;
	DynamicGraph<Weight> & g_over;

	double rnd_seed;
	CRef underprop_marker;
	CRef overprop_marker;
	CRef underprop_edge_marker;
	CRef overprop_edge_marker;

	MinimumSpanningTree<Weight> * underapprox_detector = nullptr;
	MinimumSpanningTree<Weight> * overapprox_detector = nullptr;
	MinimumSpanningTree<Weight> * underapprox_conflict_detector = nullptr;
	MinimumSpanningTree<Weight> * overapprox_conflict_detector = nullptr;

	Weight lowest_weight_lit = -1;
	Weight highest_weight_lit = -1;

	vec<int> force_reason;

	struct MSTWeightLit {
		Lit l;
		bool inclusive;
		Weight min_weight;
		MSTWeightLit() :
				l(lit_Undef),inclusive(false), min_weight(-1) {
		}
	};
	bool checked_unique;
	bool all_unique;
	vec<MSTWeightLit> weight_lits;
	struct MSTEdgeLit {
		Lit l;
		int edgeID;
		MSTEdgeLit() :
				l(lit_Undef), edgeID(-1) {
		}
	};
	vec<int> tree_edge_lits_map;
	vec<MSTEdgeLit> tree_edge_lits;
	Var first_reach_var;

	struct ChangedEdge {
		Var v;
		int edgeID;
	};
	vec<bool> is_edge_changed;
	vec<ChangedEdge> changed_edges;

	vec<Var> tmp_nodes;
	vec<bool> seen;
	vec<bool> black;
	vec<int> ancestors;

	vec<Lit> tmp_conflict;
	vec<int> visit;
	DisjointSets sets;

	struct MSTStatus {
		MSTDetector & detector;
		bool polarity;

		void setMinimumSpanningTree(Weight& weight, bool connected);
		void inMinimumSpanningTree(int edgeID, bool in_tree);

		MSTStatus(MSTDetector & _outer, bool _polarity) :
				detector(_outer), polarity(_polarity) {
		}
	};
	MSTStatus *positiveReachStatus;
	MSTStatus *negativeReachStatus;

	void unassign(Lit l) {
		Detector::unassign(l);
		int index = var(l) - first_reach_var;
		if (index >= 0 && index < tree_edge_lits_map.size() && tree_edge_lits_map[index] != -1) {
			int edgeID = tree_edge_lits_map[index];
			if (!is_edge_changed[edgeID]) {
				changed_edges.push( { var(l), edgeID });
				is_edge_changed[edgeID] = true;
			}
		}
	}
	void preprocess();
	bool propagate(vec<Lit> & conflict);
	void buildMinWeightTooSmallReason(Weight & weight, vec<Lit> & conflict);
	void buildMinWeightTooLargeReason(Weight & weight, vec<Lit> & conflict);
	void buildEdgeInTreeReason(int edge, vec<Lit> & conflict);
	void buildEdgeNotInTreeReason(int edge, vec<Lit> & conflict);

	void buildReason(Lit p, vec<Lit> & reason, CRef marker);
	bool checkSatisfied();
	Lit decide();
	void addTreeEdgeLit(int edge_id, Var reach_var);
	void addWeightLit(Var weight_var, Weight & min_weight, bool inclusive);
	void printSolution(std::ostream & write_to);
	MSTDetector(int _detectorID, GraphTheorySolver<Weight> * _outer, DynamicGraph<Weight>  &_g, DynamicGraph<Weight>  &_antig,
			double seed = 1); //:Detector(_detectorID),outer(_outer),within(-1),source(_source),rnd_seed(seed),positive_reach_detector(NULL),negative_reach_detector(NULL),positive_path_detector(NULL),positiveReachStatus(NULL),negativeReachStatus(NULL){}
	virtual ~MSTDetector() {
		if (positiveReachStatus)
			delete positiveReachStatus;
		if (negativeReachStatus)
			delete negativeReachStatus;

		if (underapprox_conflict_detector && underapprox_conflict_detector != underapprox_detector)
			delete underapprox_conflict_detector;

		if (overapprox_conflict_detector && overapprox_conflict_detector != overapprox_detector)
			delete overapprox_conflict_detector;

		if (underapprox_detector)
			delete underapprox_detector;
		if (overapprox_detector)
			delete overapprox_detector;


	}
	const char* getName() {
		return "MST Detector";
	}
	Weight getModel_SpanningTreeWeight(){
		return underapprox_detector->weight();
	}
private:
	void TarjanOLCA(int node, vec<Lit> & conflict);
	bool walkback(Weight & weight, int from, int to);
	void TarjanOLCA_edge(int node, int edgeid, int lowest_endpoint, vec<Lit> & conflict);
	Weight walkback_edge(Weight &weight, int edgeid, int from, int to, bool &found);
};
}
;

#endif /* DistanceDetector_H_ */
