/*
 * DistanceDetector.h
 *
 *  Created on: 2014-01-08
 *      Author: sam
 */

#ifndef MST_DETECTOR_H_
#define MST_DETECTOR_H_
#include "utils/System.h"

#include "GraphTheoryTypes.h"
#include "dgl/graph/DynamicGraph.h"
#include "dgl/MinimumSpanningTree.h"

#include "core/SolverTypes.h"
#include "mtl/Map.h"

#include "dgl/alg/DisjointSets.h"
#include "utils/System.h"
#include "Detector.h"

using namespace dgl;
namespace Monosat{
template<typename Weight>
class GraphTheorySolver;
template<typename Weight=int>
class MSTDetector:public Detector{
public:
		GraphTheorySolver<Weight> * outer;
		//int within;
		DynamicGraph & g;
		DynamicGraph & antig;

		double rnd_seed;
		CRef reach_marker;
		CRef non_reach_marker;
		CRef reach_edge_marker;
		CRef non_reach_edge_marker;

		MinimumSpanningTree<Weight> * positive_reach_detector;
		MinimumSpanningTree<Weight> * negative_reach_detector;
		MinimumSpanningTree<Weight> *  positive_conflict_detector;
		MinimumSpanningTree<Weight> * negative_conflict_detector;
		//Reach *  positive_path_detector;

		//vec<Lit>  reach_lits;

		//Map<float,int> weight_lit_map;
		vec<int> force_reason;
		std::vector<Weight> & edge_weights;
		struct MSTWeightLit{
			Lit l;
			Weight min_weight;
			MSTWeightLit():l(lit_Undef),min_weight(-1){}
		};
		bool checked_unique;
		bool all_unique;
		vec<MSTWeightLit>  weight_lits;
		struct MSTEdgeLit{
			Lit l;
			int edgeID;
			MSTEdgeLit():l(lit_Undef),edgeID(-1){}
		};
		vec<int> tree_edge_lits_map;
		vec<MSTEdgeLit>  tree_edge_lits;
		Var first_reach_var;
		struct ChangedWeight{
			Lit l;
			Weight weight;
		};
		vec<ChangedWeight> changed_weights;

		struct ChangedEdge{
			Lit l;
			int edgeID;
		};
		vec<ChangedEdge> changed_edges;

		vec<Var> tmp_nodes;
		vec<bool> seen;
		vec<bool> black;
		vec<int> ancestors;

		vec<Lit> tmp_conflict;
		vec<int> visit;
		DisjointSets sets;

		struct MSTStatus{
			MSTDetector & detector;
			bool polarity;

			void setMinimumSpanningTree(Weight& weight);
			void inMinimumSpanningTree(int edgeID, bool in_tree);

			MSTStatus(MSTDetector & _outer, bool _polarity):detector(_outer), polarity(_polarity){}
		};
		MSTStatus *positiveReachStatus;
		MSTStatus *negativeReachStatus;



		bool propagate(vec<Lit> & conflict);
		void buildMinWeightTooSmallReason(Weight & weight,vec<Lit> & conflict);
		void buildMinWeightTooLargeReason(Weight & weight,vec<Lit> & conflict);
		void buildEdgeInTreeReason(int edge,vec<Lit> & conflict);
		void buildEdgeNotInTreeReason(int edge,vec<Lit> & conflict);
		//void buildForcedMinWeightReason(int reach_node, int forced_edge_id,vec<Lit> & conflict);
		void buildReason(Lit p, vec<Lit> & reason, CRef marker);
		bool checkSatisfied();
		Lit decide();
		void addTreeEdgeLit(int edge_id, Var reach_var);
		void addWeightLit(Var weight_var,Weight & min_weight);

		MSTDetector(int _detectorID, GraphTheorySolver<Weight> * _outer, DynamicGraph &_g, DynamicGraph &_antig, std::vector<Weight> & _edge_weights,  double seed=1);//:Detector(_detectorID),outer(_outer),within(-1),source(_source),rnd_seed(seed),positive_reach_detector(NULL),negative_reach_detector(NULL),positive_path_detector(NULL),positiveReachStatus(NULL),negativeReachStatus(NULL){}
		virtual ~MSTDetector(){

		}
		const char* getName(){
			return "MST Detector";
		}
private:
		void TarjanOLCA(int node, vec<Lit> & conflict);
		bool walkback(Weight & weight, int from, int to);
		void TarjanOLCA_edge(int node, int edgeid,int lowest_endpoint, vec<Lit> & conflict);
		Weight walkback_edge(Weight &weight,int edgeid, int from, int to, bool &found);
};
};



#endif /* DistanceDetector_H_ */
