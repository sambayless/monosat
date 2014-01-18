/*
 * DistanceDetector.h
 *
 *  Created on: 2014-01-08
 *      Author: sam
 */

#ifndef MST_DETECTOR_H_
#define MST_DETECTOR_H_
#include "utils/System.h"

#include "Graph.h"
#include "MinimumSpanningTree.h"

#include "core/SolverTypes.h"
#include "mtl/Map.h"
#include "Kruskal.h"
#include "DisjointSets.h"
#include "utils/System.h"
#include "Detector.h"
namespace Minisat{
class GraphTheorySolver;
class MSTDetector:public Detector{
public:
		GraphTheorySolver * outer;
		//int within;
		DynamicGraph<PositiveEdgeStatus> & g;
		DynamicGraph<NegativeEdgeStatus> & antig;

		double rnd_seed;

		MinimumSpanningTree * positive_reach_detector;
		MinimumSpanningTree * negative_reach_detector;
		//Reach *  positive_path_detector;

		//vec<Lit>  reach_lits;

		Map<float,int> weight_lit_map;
		vec<int> force_reason;
		vec<int> edge_weights;
		struct MSTWeightLit{
			Lit l;
			int min_weight;
			MSTWeightLit():l(lit_Undef),min_weight(-1){}
		};
		vec<MSTWeightLit>  weight_lits;
		struct MSTEdgeLit{
			Lit l;
			int edgeID;
			MSTEdgeLit():l(lit_Undef),edgeID(-1){}
		};
		vec<MSTEdgeLit>  tree_edge_lits;
		struct ChangedWeight{
			Lit l;
			int weight;
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
		DisjointSets sets;

		struct MSTStatus{
			MSTDetector & detector;
			bool polarity;

			void setMinimumSpanningTree(int weight);
			void inMinimumSpanningTree(int edgeID, bool in_tree);

			MSTStatus(MSTDetector & _outer, bool _polarity):detector(_outer), polarity(_polarity){}
		};
		MSTStatus *positiveReachStatus;
		MSTStatus *negativeReachStatus;



		bool propagate(vec<Assignment> & trail,vec<Lit> & conflict);
		void buildMinWeightReason(int weight,vec<Lit> & conflict);
		void buildNonMinWeightReason(int weight,vec<Lit> & conflict);
		void buildEdgeReason(int edge,vec<Lit> & conflict);
		void buildNonEdgeReason(int edge,vec<Lit> & conflict);
		//void buildForcedMinWeightReason(int reach_node, int forced_edge_id,vec<Lit> & conflict);
		void buildReason(Lit p, vec<Lit> & reason, CRef marker);
		bool checkSatisfied();
		Lit decide();
		void addTreeEdgeLit(int edge_id, Var reach_var);
		void addWeightLit(Var weight_var,int min_weight);

		MSTDetector(int _detectorID, GraphTheorySolver * _outer, DynamicGraph<PositiveEdgeStatus> &_g, DynamicGraph<NegativeEdgeStatus> &_antig,  double seed=1);//:Detector(_detectorID),outer(_outer),within(-1),source(_source),rnd_seed(seed),positive_reach_detector(NULL),negative_reach_detector(NULL),positive_path_detector(NULL),positiveReachStatus(NULL),negativeReachStatus(NULL){}
		virtual ~MSTDetector(){

		}
private:
		void TarjanOLCA(int node, vec<Lit> & conflict);
		bool walkback(int weight, int from, int to);
		void TarjanOLCA_edge(int node, int edgeid,int lowest_endpoint, vec<Lit> & conflict);
		int walkback_edge(int weight,int edgeid, int from, int to, bool &found);
};
};
#endif /* DistanceDetector_H_ */
