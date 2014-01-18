
#ifndef PRIMS_H_
#define PRIMS_H_

#include "mtl/Vec.h"
#include "mtl/Heap.h"
#include "DynamicGraph.h"
#include "core/Config.h"
#include "MinimumSpanningTree.h"
#include "DisjointSets.h"
#include <limits>
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
	bool hasParents;
	int INF;
	DisjointSets sets;
	vec<int> mst;
	vec<int> q;
	vec<int> check;
	const int reportPolarity;

	//vec<char> old_seen;
	vec<bool> in_tree;;
	vec<int> parents;
	vec<int> parent_edges;
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
		hasParents=false;
	}

	void setNodes(int n){
		q.capacity(n);
		check.capacity(n);
		in_tree.growTo(g.edges);

		INF=std::numeric_limits<int>::max();
		sets.AddElements(n);
		parents.growTo(n);
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
		hasParents=false;
		sets.Reset();
		setNodes(g.nodes);

		min_weight=0;

		mst.clear();


		for(int i = 0;i<g.edges;i++){
			if(g.edgeEnabled(i)){
				edge_heap.insert(i);
			}
			in_tree[i]=false;
		}

		while(edge_heap.size()){
			int edge_id = edge_heap.removeMin();
			int u = g.getEdge(edge_id).from;
			int v = g.getEdge(edge_id).to;
			int set1 = sets.FindSet(u);
			int set2 = sets.FindSet(v);
			if(set1!=set2){
				assert(g.edgeEnabled(edge_id));
				assert(parents[v]==-1);
				in_tree[edge_id]=true;
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
		update();
		return mst;
	 }

	int getParent(int node){
		update();
		//kruskals doesn't actually give us the parents, normally. need to construct it on demand, here.
		if (!hasParents)
			buildParents();
		return parents[node];
	}
	 int getParentEdge(int node){
		 if(getParent(node)!=-1)
			 return parent_edges[node];
		 else
			 return -1;
	 }
	bool edgeInTree(int edgeid){
		update();
		return in_tree[edgeid];
	}
	bool dbg_mst(){

		return true;
	}
	int getEdgeWeight(int edgeID){
		return edge_weights[edgeID];
	}

	int weight(){

		update();

		assert(dbg_uptodate());

		return min_weight;
	}

	bool dbg_uptodate(){
		return true;
	};
private:
	void buildParents(){
		hasParents=true;
		for(int i = 0;i<parents.size();i++){
			parents[i]=-1;
			parent_edges[i]=-1;
		}
		q.clear();
		q.push(0);
		while(q.size()){
			int u = q.last();
			assert(parents[u]==-1);
			q.pop();
			for (int j = 0;j<g.adjacency[u].size();j++){
				int edge = g.adjacency[u][j].id;
				int to = g.adjacency[u][j].node;
				if(in_tree[edge]){
					if(parents[to]==-1){
						parents[to]=u;
						parent_edges[to] = edge;
						q.push(to);
					}
				}
			}
		}
	}
};

#endif
