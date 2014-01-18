
#ifndef PRIMS_H_
#define PRIMS_H_

#include "mtl/Vec.h"
#include "mtl/Heap.h"
#include "DynamicGraph.h"
#include "core/Config.h"
#include "MinimumSpanningTree.h"
#include "DisjointSets.h"
using namespace Minisat;



template<class Status,class EdgeStatus=DefaultEdgeStatus>
class Kruskal:public MinimumSpanningTree{
public:

	DynamicGraph<EdgeStatus> & g;
	Status &  status;
	int last_modification;
	int min_weight;
	int last_addition;
	int last_deletion;
	int history_qhead;

	int last_history_clear;

	int INF;
	DisjointSets sets;
	vec<int> mst;
	vec<int> q;
	vec<int> check;
	const int reportPolarity;

	//vec<char> old_seen;
	vec<char> seen;;
//	vec<int> changed;
	vec<int> edge_weights;
    struct EdgeLt {
        const vec<int>&  edge_weights;

        bool operator () (int x, int y) const {
        	return edge_weights[x]<edge_weights[y];
        }
        EdgeLt(const vec<int>&  _edge_weights) : edge_weights(_edge_weights) { }
    };

	Heap<EdgeLt> edge_heap;


	vec<int> prev;

	struct DefaultReachStatus{
			vec<bool> stat;
				void setReachable(int u, bool reachable){
					stat.growTo(u+1);
					stat[u]=reachable;
				}
				bool isReachable(int u) const{
					return stat[u];
				}
				DefaultReachStatus(){}
			};

public:


	Kruskal(DynamicGraph<EdgeStatus> & graph, Status & _status, int _reportPolarity=0 ):g(graph), status(_status), last_modification(-1),last_addition(-1),last_deletion(-1),history_qhead(0),last_history_clear(0),INF(0),reportPolarity(_reportPolarity),edge_heap(EdgeLt(edge_weights)){
		marked=false;
		mod_percentage=0.2;
		stats_full_updates=0;
		stats_fast_updates=0;
		stats_skip_deletes=0;
		stats_skipped_updates=0;
		stats_full_update_time=0;
		stats_fast_update_time=0;
		stats_num_skipable_deletions=0;
		stats_fast_failed_updates=0;
		min_weight=-1;
	}

	void setNodes(int n){
		q.capacity(n);
		check.capacity(n);
		seen.growTo(n);
		prev.growTo(n);
		INF=g.nodes+1;
		sets.AddElements(n);

	}

	void update( ){
		static int iteration = 0;
		int local_it = ++iteration ;
		stats_full_updates++;
		double startdupdatetime = cpuTime();
		if(last_modification>0 && g.modifications==last_modification){
			stats_skipped_updates++;
			return;
		}

		if(last_deletion==g.deletions){
			stats_num_skipable_deletions++;
		}

		sets.Reset();
		setNodes(g.nodes);

		min_weight=0;

		mst.clear();


		for(int i = 0;i<g.edges;i++){
			if(g.edgeEnabled(i)){
				edge_heap.insert(i);
			}
		}

		while(edge_heap.size()){
			int edge_id = edge_heap.removeMin();
			int u = g.getEdge(edge_id).from;
			int v = g.getEdge(edge_id).to;
			int set1 = sets.FindSet(u);
			int set2 = sets.FindSet(v);
			if(set1!=set2){
				assert(g.edgeEnabled(edge_id));
				mst.push(edge_id);
				if(reportPolarity>-1)
					status.inMinimumSpanningTree(edge_id,true);
				min_weight+=getEdgeWeight(edge_id);
				sets.Union(set1,set2);
			}
		}

		status.setMinimumSpanningTree(min_weight);

		assert(dbg_uptodate());

		last_modification=g.modifications;
		last_deletion = g.deletions;
		last_addition=g.additions;

		history_qhead=g.history.size();
		last_history_clear=g.historyclears;



		stats_full_update_time+=cpuTime()-startdupdatetime;;
	}
	vec<int> & getSpanningTree(){
		return mst;
	 }
	bool dbg_mst(){

		return true;
	}
	int getEdgeWeight(int edgeID){
		return edge_weights[edgeID];
	}

	int weight(){
		if(last_modification!=g.modifications)
			update();

		assert(dbg_uptodate());

		return min_weight;
	}

	bool dbg_uptodate(){
		return true;
	};

};

#endif
