#ifndef LINN_CUT_CONNECTEDCOMPONENTS_H_
#define LINN_CUT_CONNECTEDCOMPONENTS_H_



#include "mtl/Vec.h"
#include "mtl/Heap.h"
#include "DynamicGraph.h"
#include "core/Config.h"
#include "ConnectedComponents.h"
#include "DisjointSets.h"
#include "LinkCut.h"
#include <limits>
using namespace Minisat;



template<class EdgeStatus=DefaultEdgeStatus>
class LinkCutForest:public ConnectedComponents{
public:

	DynamicGraph<EdgeStatus> & g;

	int last_modification;

	int last_addition;
	int last_deletion;
	int history_qhead;

	int last_history_clear;
	bool hasParents;
	int INF;
	LinkCut sets;
	int cycleID;

#ifndef NDEBUG
	DisjointSets dbg_sets;
#endif
	vec<int> q;
	vec<int> check;
	const int reportPolarity;

public:
	//stats

	int stats_full_updates;
	int stats_fast_updates;
	int stats_fast_failed_updates;
	int stats_skip_deletes;
	int stats_skipped_updates;
	int stats_num_skipable_deletions;
	double stats_full_update_time;
	double stats_fast_update_time;
	double mod_percentage;

	LinkCutForest(DynamicGraph<EdgeStatus> & graph, int _reportPolarity=0 ):g(graph), last_modification(-1),last_addition(-1),last_deletion(-1),history_qhead(0),last_history_clear(0),INF(0),reportPolarity(_reportPolarity){

		mod_percentage=0.2;
		stats_full_updates=0;
		stats_fast_updates=0;
		stats_skip_deletes=0;
		stats_skipped_updates=0;
		stats_full_update_time=0;
		stats_fast_update_time=0;
		stats_num_skipable_deletions=0;
		stats_fast_failed_updates=0;
		hasParents=false;
		cycleID=-1;
	}

	void setNodes(int n){
		q.capacity(n);
		check.capacity(n);

		INF=std::numeric_limits<int>::max();
		//sets.AddElements(n);
		while(sets.nNodes()<n)
			sets.addNode();
	}

	void update( ){
		static int iteration = 0;
		int local_it = ++iteration ;

		if(last_modification>0 && g.modifications==last_modification){
			stats_skipped_updates++;
			return;
		}
		stats_full_updates++;
		double startdupdatetime = cpuTime();
		if(last_deletion==g.deletions){
			stats_num_skipable_deletions++;
		}
		hasParents=false;

		setNodes(g.nodes);


#ifndef NDEBUG
		dbg_sets.Reset();
		for(int i = 0;i<g.edges;i++){
			if(g.edgeEnabled(i)){
				int u = g.all_edges[i].from;
				int v = g.all_edges[i].to;
				dbg_sets.UnionElements(u,v);
			}
		}

#endif
		cycleID = -1;
		if(g.historyclears!=last_history_clear){
			last_history_clear=g.historyclears;
			history_qhead=0;
			//start from scratch
		}else{
			//incremental/decremental update
			for(;history_qhead<g.history.size();history_qhead++){
				int edgeid = g.history[history_qhead].id;
				bool add = g.history[history_qhead].addition;
				int u =  g.history[history_qhead].u;
				int v =  g.history[history_qhead].v;
				if(add){
					if(sets.connected(u,v)){
						//then adding this edge would produce an (undirected) cycle.
						cycleID= edgeid;
						break;
					}
					sets.link(u,v);
				}else{
					if(sets.connected(u,v)){
						sets.cut(u,v);
					}
				}
			}
		}

		if(cycleID>-1){


		}else{
			assert(dbg_sets.NumSets()== sets.numRoots());
		}

		last_modification=g.modifications;
		last_deletion = g.deletions;
		last_addition=g.additions;

		history_qhead=g.history.size();
		last_history_clear=g.historyclears;

		stats_full_update_time+=cpuTime()-startdupdatetime;;
	}

	bool connected(int from, int to){
		update();
		return sets.findRoot(from)==sets.findRoot(to);
	}

	bool hasCycle(){
		update();
		return cycleID>=0;
	}

	int getCycleEdge(){
		assert(hasCycle());
		update();
		return cycleID;
	}

	 int numComponents(){
		 update();
		 return sets.numRoots();
	 }
	int getComponent(int node){
		update();
		return sets.findRoot(node);
	}


	bool dbg_uptodate(){
		return true;
	};


};


#endif /* DISJOINTSETCONNECTEDCOMPONENTS_H_ */
