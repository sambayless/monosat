
#ifndef PRIMS_H_
#define PRIMS_H_

#include "mtl/Vec.h"
#include "mtl/Heap.h"
#include "DynamicGraph.h"
#include "core/Config.h"
#include "MinimumSpanningTree.h"
#include <limits>
using namespace Minisat;



template<class Status,class EdgeStatus=DefaultEdgeStatus>
class Prim:public MinimumSpanningTree{
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
	vec<int> keys;
	vec<int> parents;
	vec<int> parent_edges;
    struct VertLt {
        const vec<int>&  keys;

        bool operator () (int x, int y) const {
        	return keys[x]<keys[y];
        }
        VertLt(const vec<int>&  _key) : keys(_key) { }
    };

	Heap<VertLt> Q;

	vec<int> mst;
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


	Prim(DynamicGraph<EdgeStatus> & graph, Status & _status, int _reportPolarity=0 ):g(graph), status(_status), last_modification(-1),last_addition(-1),last_deletion(-1),history_qhead(0),last_history_clear(0),INF(0),reportPolarity(_reportPolarity),Q(VertLt(keys)){
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
		INF=std::numeric_limits<int>::max();
		sets.AddElements(n);
		parents.growTo(n);
		keys.growTo(n);
		parent_edges.growTo(n);
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

		sets.Reset();
		setNodes(g.nodes);

		min_weight=0;

		mst.clear();


		for(int i = 0;i<g.nodes;i++){
			parents[i]=-1;
			int key = INF;
			keys[i]=key;
		}
		Q.insert(0);//arbitrary first node to examine

		while(Q.size()){
			int u = Q.removeMin();
			seen[u]=true;
			if(u!=0){
				int parent = parents[u];
				assert(parent!=-1);
				int edgeid = g.adjacency[parent][u].id;
				mst.push(edgeid);
			}
			for(int j = 0;j< g.adjacency[u].size();j++){
				int edgeid = g.adjacency[u][j].id;
				int v = g.adjacency[u][j].node;
				if(!seen[v]){
					int w = edge_weights[edgeid];
					if(w<keys[v]){
						parents[v]=u;
						parent_edges[v]=edgeid;
						keys[v]= w;

						Q.update(v);
					}
				}
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
		update();
		return mst;
	 }
	 int getParent(int node){
		 update();
		 return parents[node];
	 }
	 int getParentEdge(int node){
		 if(getParent(node)!=-1)
			 return parent_edges[node];
		 else
			 return -1;
	 }
	bool dbg_mst(){

		return true;
	}
	int getEdgeWeight(int edgeID){
		return edge_weights[edgeID];
	}
	bool edgeInTree(int edgeid){
		update();
		int u = g.all_edges[edgeid].from;
		int v = g.all_edges[edgeid].to;
		return parents[u]==v || parents[v]==u;
	}
	int weight(){
		update();

		assert(dbg_uptodate());

		return min_weight;
	}

	bool dbg_uptodate(){
		return true;
	};


};

#endif
