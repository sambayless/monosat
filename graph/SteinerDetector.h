/*
 * DistanceDetector.h
 *
 *  Created on: 2014-01-08
 *      Author: sam
 */

#ifndef STEINER_DETECTOR_H_
#define STEINER_DETECTOR_H_
#include "utils/System.h"

#include "Graph.h"
#include "dgl/MinimumSpanningTree.h"
#include "dgl/SteinerTree.h"
#include "core/SolverTypes.h"
#include "mtl/Map.h"

#include "dgl/alg/DisjointSets.h"
#include "dgl/graph/DynamicNodes.h"
#include "utils/System.h"
#include "Detector.h"

using namespace dgl;
namespace Minisat{
class GraphTheorySolver;
class SteinerDetector:public Detector{
public:
		GraphTheorySolver * outer;

		DynamicGraph & g;
		DynamicGraph & antig;

		DynamicNodes underTerminalSet;
		DynamicNodes overTerminalSet;
		double rnd_seed;
		CRef reach_marker;
		CRef non_reach_marker;
		CRef reach_edge_marker;
		CRef non_reach_edge_marker;

;


		SteinerTree * positive_reach_detector;
		SteinerTree * negative_reach_detector;
		SteinerTree *  positive_conflict_detector;
		SteinerTree * negative_conflict_detector;

		vec<Var> terminal_map;
		vec<int> terminal_var_map;

		Map<float,int> weight_lit_map;

		vec<int> & edge_weights;
		struct WeightLit{
			Lit l;
			int min_weight;
			WeightLit():l(lit_Undef),min_weight(-1){}
		};
		bool checked_unique;
		bool all_unique;
		vec<WeightLit>  weight_lits;

		Var first_reach_var;
		struct ChangedWeight{
			Lit l;
			int weight;
		};
		vec<ChangedWeight> changed_weights;
		vec<Var> tmp_nodes;
		vec<bool> seen;
		vec<bool> black;
		vec<int> ancestors;

		vec<Lit> tmp_conflict;
		vec<int> visit;
		DisjointSets sets;
		struct SteinerStatus{
			SteinerDetector & detector;
			bool polarity;

			void setMinimumSteinerTree(int weight);

			SteinerStatus(SteinerDetector & _outer, bool _polarity):detector(_outer), polarity(_polarity){}
		};

		SteinerStatus *negativeStatus;
		SteinerStatus *positiveStatus;

		bool propagate(vec<Lit> & conflict);
		void buildMinWeightTooSmallReason(int weight,vec<Lit> & conflict);
		void buildMinWeightTooLargeReason(int weight,vec<Lit> & conflict);

		void buildReason(Lit p, vec<Lit> & reason, CRef marker);
		bool checkSatisfied();
		Lit decide();

		void addWeightLit(Var weight_var,int min_weight);
		void addTerminalNode(int node,Var theoryVar);
		SteinerDetector(int _detectorID, GraphTheorySolver * _outer, DynamicGraph &_g, DynamicGraph &_antig,  double seed=1);//:Detector(_detectorID),outer(_outer),within(-1),source(_source),rnd_seed(seed),positive_reach_detector(NULL),negative_reach_detector(NULL),positive_path_detector(NULL),positiveReachStatus(NULL),negativeReachStatus(NULL){}

		virtual void assign(Lit l){
			 Detector::assign(l);
			 if(var(l)<terminal_var_map.size()){
			 int node =terminal_var_map[var(l)];
			 if(node>=0){
				 if(sign(l)){
					 underTerminalSet.setNodeEnabled(node,false);
					 overTerminalSet.setNodeEnabled(node,false);
				 }else{
					 underTerminalSet.setNodeEnabled(node,true);
					 overTerminalSet.setNodeEnabled(node,true);
				 }
			 }
			 }
		}
		virtual	 void unassign(Lit l){
			 Detector::unassign(l);
			 if(var(l)<terminal_var_map.size()){
			 int node =terminal_var_map[var(l)];
			 if(node>=0){
				 underTerminalSet.setNodeEnabled(node,false);
				 overTerminalSet.setNodeEnabled(node,true);
			 }
			 }
		}

		const char* getName(){
			return "Steiner Detector";
		}

};
};
#endif /* DistanceDetector_H_ */
