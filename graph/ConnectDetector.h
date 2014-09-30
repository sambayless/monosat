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
#ifndef CONNECTDETECTOR_H_
#define CONNECTDETECTOR_H_
#include "utils/System.h"

#include "GraphTheoryTypes.h"
#include "dgl/graph/DynamicGraph.h"
#include "dgl/Reach.h"
#include "dgl/Dijkstra.h"
#include "dgl/BFS.h"
#include "dgl/DFS.h"

#include "core/SolverTypes.h"
#include "mtl/Map.h"
#include "dgl/MaxFlow.h"

#include "dgl/EdmondsKarp.h"
#include "dgl/EdmondsKarpAdj.h"
#include "dgl/Chokepoint.h"
#include "WeightedDijkstra.h"

#include "utils/System.h"
#include "Detector.h"
using namespace dgl;
namespace Monosat{
template<typename Weight>
class GraphTheorySolver;
template<typename Weight>
class ConnectDetector:public Detector{
public:
		GraphTheorySolver<Weight> * outer;
		 DynamicGraph &g;
		 DynamicGraph &antig;
		int within;
		int source;
		int constraintsBuilt;
		double rnd_seed;
		CRef reach_marker;
		CRef non_reach_marker;
		CRef forced_reach_marker;

		Reach * positive_reach_detector;
		Reach * negative_reach_detector;
		Reach *  positive_path_detector;

		vec<vec<Lit>>  dist_lits;

		vec<Lit>  reach_lits;
		Var first_reach_var;
		vec<int> reach_lit_map;
		vec<int> force_reason;


		std::vector<ForceReason> forced_edges;
		struct Change{
				Lit l;
				int u;
			};
		vec<Change> changed;

		vec<Change> & getChanged(){
			return changed;
		}

		//stats

		int stats_full_updates;
		int stats_fast_updates;
		int stats_fast_failed_updates;
		int stats_skip_deletes;
		int stats_skipped_updates;
		int stats_num_skipable_deletions;
		double mod_percentage;

		double stats_full_update_time;
		double stats_fast_update_time;


		struct ReachStatus{
			ConnectDetector & detector;
			bool polarity;
			void setReachable(int u, bool reachable);
			bool isReachable(int u) const{
				return false;
			}

			void setMininumDistance(int u, bool reachable, int distance);


			ReachStatus(ConnectDetector & _outer, bool _polarity):detector(_outer), polarity(_polarity){}
		};
		ReachStatus *positiveReachStatus;
		ReachStatus *negativeReachStatus;
		WeightedDijkstra<double> * rnd_path;
		std::vector<double> rnd_weight;

		//WeightedDijkstra<double> * opt_path;
		bool check_positive;
		bool check_negative;

		struct ChokepointStatus{
			ConnectDetector & detector;
			bool mustReach(int node);
			bool operator() (int edge_id);
			ChokepointStatus(ConnectDetector & _outer):detector(_outer){

			}
		}chokepoint_status;

		Chokepoint<ChokepointStatus> chokepoint;

		int getNode(Var reachVar){
			assert(reachVar>=first_reach_var);
			int index = reachVar-first_reach_var;
			assert(index< reach_lit_map.size());
			assert(reach_lit_map[index]>=0);
			return reach_lit_map[index];
		}

		void buildSATConstraints(int distance=-1);
		bool propagate(vec<Lit> & conflict);
		void buildReachReason(int node,vec<Lit> & conflict);
		void buildNonReachReason(int node,vec<Lit> & conflict);
		void buildForcedEdgeReason(int reach_node, int forced_edge_id,vec<Lit> & conflict);
		void buildReason(Lit p, vec<Lit> & reason, CRef marker);
		bool checkSatisfied();
		void addLit(int from, int to, Var reach_var);
		Lit decide();
		void preprocess();
		void dbg_sync_reachability();

		ConnectDetector(int _detectorID, GraphTheorySolver<Weight> * _outer, DynamicGraph &_g, DynamicGraph &_antig, int _source,double seed=1);//:Detector(_detectorID),outer(_outer),within(-1),source(_source),rnd_seed(seed),positive_reach_detector(NULL),negative_reach_detector(NULL),positive_path_detector(NULL),positiveReachStatus(NULL),negativeReachStatus(NULL),chokepoint_status(*this),chokepoint(chokepoint_status, _antig,source){}
		virtual ~ConnectDetector(){

		}
		const char* getName(){
			return "Connected Detector";
		}
};
};
#endif /* REACHDETECTOR_H_ */
