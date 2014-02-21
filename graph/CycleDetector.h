
#ifndef CYCLE_DETECTOR_H_
#define CYCLE_DETECTOR_H_
#include "utils/System.h"

#include "Graph.h"
#include "ConnectedComponents.h"
#include "Cycle.h"
#include "core/SolverTypes.h"
#include "mtl/Map.h"
#include "DFSCycle.h"
#include "utils/System.h"
#include "Detector.h"
namespace Minisat{
class GraphTheorySolver;
class CycleDetector:public Detector{
public:
		GraphTheorySolver * outer;
		//int within;
		DynamicGraph<PositiveEdgeStatus> & g;
		DynamicGraph<NegativeEdgeStatus> & antig;

		double rnd_seed;
		CRef directed_cycle_marker;
		CRef no_directed_cycle_marker;

		CRef undirected_cycle_marker;
		CRef no_undirected_cycle_marker;
		Cycle * positive_reach_detector;
		Cycle * negative_reach_detector;

		//Reach *  positive_path_detector;

		//vec<Lit>  reach_lits;


		vec<Var> tmp_nodes;
		vec<bool> seen;
		vec<bool> black;
		vec<int> ancestors;

		vec<bool> edge_in_clause;
		vec<int> visit;

		Lit undirected_cycle_lit;
		Lit directed_cycle_lit;




		bool propagate(vec<Assignment> & trail,vec<Lit> & conflict);

		void buildNoUndirectedCycleReason(vec<Lit> & conflict);
			void buildNoDirectedCycleReason(vec<Lit> & conflict);
			void buildUndirectedCycleReason(vec<Lit> & conflict);
			void buildDirectedCycleReason(vec<Lit> & conflict);
		//void buildForcedMinWeightReason(int reach_node, int forced_edge_id,vec<Lit> & conflict);
		void buildReason(Lit p, vec<Lit> & reason, CRef marker);
		bool checkSatisfied();
		Lit decide();
		void addCycleDetectorLit(bool undirected, Var v);

		CycleDetector(int _detectorID, GraphTheorySolver * _outer, DynamicGraph<PositiveEdgeStatus> &_g, DynamicGraph<NegativeEdgeStatus> &_antig, bool detect_directed_cycles=true,   double seed=1);//:Detector(_detectorID),outer(_outer),within(-1),source(_source),rnd_seed(seed),positive_reach_detector(NULL),negative_reach_detector(NULL),positive_path_detector(NULL),positiveReachStatus(NULL),negativeReachStatus(NULL){}
		virtual ~CycleDetector(){

		}

};
};
#endif /* DistanceDetector_H_ */
