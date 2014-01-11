/*
 * AllPairsDetector.h
 *
 *  Created on: 2014-01-08
 *      Author: sam
 */

#ifndef ALLPAIRSDETECTOR_H_
#define ALLPAIRSDETECTOR_H_
#include "utils/System.h"

#include "Graph.h"
#include "AllPairs.h"
#include "FloydWarshall.h"
#include "DijkstraAllPairs.h"
#include "core/SolverTypes.h"
#include "mtl/Map.h"

#define DEBUG_ALLPAIRS
#include "utils/System.h"
#include "Detector.h"
namespace Minisat{
class GraphTheorySolver;
class AllPairsDetector:public Detector{
public:
		GraphTheorySolver * outer;
		int within;

		double rnd_seed;
#ifdef DEBUG_ALLPAIRS
		AllPairs * dbg_positive_reach_detector;
		AllPairs * dbg_negative_reach_detector;
#endif
		AllPairs * positive_reach_detector;
		AllPairs * negative_reach_detector;
		AllPairs *  positive_path_detector;

		//vec<Lit>  reach_lits;
		Var first_reach_var;

		vec<int> force_reason;

		struct DistLit{
			Lit l;
			int min_distance;
			int source;
			int node;
		};
		vec<bool> installed_sources;
		vec<int> sources;
		vec<vec<vec<DistLit> >> dist_lits;
		vec<DistLit> reach_lit_map;

		vec<int> tmp_path;

		struct Change{
				Lit l;
				int u;
				int source;
			};
		vec<Change> changed;

		vec<Change> & getChanged(){
			return changed;
		}
		struct IgnoreStatus{

					void setReachable(int from,int u, bool reachable){

					}
					bool isReachable(int from,int u) const{
						return false;
					}

					void setMininumDistance(int from,int u, bool reachable, int distance){

					}
					IgnoreStatus(){}
				}ignoreStatus;
		struct ReachStatus{
			AllPairsDetector & detector;
			bool polarity;
			void setReachable(int from,int u, bool reachable);
			bool isReachable(int from,int u) const{
				return false;
			}

			void setMininumDistance(int from,int u, bool reachable, int distance);


			ReachStatus(AllPairsDetector & _outer, bool _polarity):detector(_outer), polarity(_polarity){}
		};
		ReachStatus *positiveReachStatus;
		ReachStatus *negativeReachStatus;

		int getNode(Var reachVar){
			assert(reachVar>=first_reach_var);
			int index = reachVar-first_reach_var;
			assert(index< reach_lit_map.size());
			assert(reach_lit_map[index].node>=0);
			return reach_lit_map[index].node;
		}
		int getSource(Var reachVar){
			assert(reachVar>=first_reach_var);
			int index = reachVar-first_reach_var;
			assert(index< reach_lit_map.size());
			assert(reach_lit_map[index].source>=0);
			return reach_lit_map[index].source;
		}
	/*	Lit getLit(int node){

			return reach_lits[node];

		}*/
		bool propagate(vec<Assignment> & trail,vec<Lit> & conflict);
		void buildReachReason(int from, int to,vec<Lit> & conflict);
		void buildNonReachReason(int from, int to,vec<Lit> & conflict);
		void buildForcedEdgeReason(int from, int to, int forced_edge_id,vec<Lit> & conflict);
		void buildReason(Lit p, vec<Lit> & reason, CRef marker);
		bool checkSatisfied();
		Lit decide();
		void addLit(int from, int to, Var reach_var,int within_steps=-1);
		AllPairsDetector(int _detectorID, GraphTheorySolver * _outer, DynamicGraph<PositiveEdgeStatus> &_g, DynamicGraph<NegativeEdgeStatus> &_antig,  int within_steps,double seed=1);//:Detector(_detectorID),outer(_outer),within(-1),source(_source),rnd_seed(seed),positive_reach_detector(NULL),negative_reach_detector(NULL),positive_path_detector(NULL),positiveReachStatus(NULL),negativeReachStatus(NULL){}
		virtual ~AllPairsDetector(){

		}
};
};
#endif /* AllPairsDetector_H_ */
