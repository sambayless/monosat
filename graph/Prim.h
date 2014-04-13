
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
	int next_component;
	vec<int> components;
	vec<int> roots;
    struct VertLt {
        const vec<int>&  keys;

        bool operator () (int x, int y) const {
        	return keys[x]<keys[y];
        }
        VertLt(const vec<int>&  _key) : keys(_key) { }
    };
    bool hasComponents;
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

	int stats_full_updates;
	int stats_fast_updates;
	int stats_fast_failed_updates;
	int stats_skip_deletes;
	int stats_skipped_updates;
	int stats_num_skipable_deletions;
	double mod_percentage;

	double stats_full_update_time;
	double stats_fast_update_time;

	Prim(DynamicGraph<EdgeStatus> & graph, Status & _status, int _reportPolarity=0 ):g(graph), status(_status), last_modification(-1),last_addition(-1),last_deletion(-1),history_qhead(0),last_history_clear(0),INF(0),reportPolarity(_reportPolarity),Q(VertLt(keys)){

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
		hasComponents=false;
	}

	void setNodes(int n){
		q.capacity(n);
		check.capacity(n);
		seen.growTo(n);
		prev.growTo(n);
		INF=std::numeric_limits<int>::max();

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
		double startdupdatetime = rtime(2);
		if(last_deletion==g.deletions){
			stats_num_skipable_deletions++;
		}
		hasComponents=false;

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
				int edgeid = g.adjacency_undirected[parent][u].id;
				mst.push(edgeid);
			}
			for(int j = 0;j< g.adjacency_undirected[u].size();j++){
				int edgeid = g.adjacency_undirected[u][j].id;
				int v = g.adjacency_undirected[u][j].node;
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
		int numsets = 0;

		for(int i = 0;i<g.nodes;i++){
			if(parents[i] ==-1) {
				numsets++;
			}

		}
		if(numsets>1){
			min_weight=INF;
		}
		if(numsets==0){
			assert(min_weight==0);
		}
		status.setMinimumSpanningTree(min_weight);

		assert(dbg_uptodate());

		last_modification=g.modifications;
		last_deletion = g.deletions;
		last_addition=g.additions;

		history_qhead=g.history.size();
		last_history_clear=g.historyclears;



		stats_full_update_time+=rtime(2)-startdupdatetime;;
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
	int numComponents(){
		update();
		buildComponents();
		return next_component;
	}
	int getComponent(int node){
		update();
		buildComponents();
		assert(components[node]!=-1);
		return components[node];
	}
	int getRoot(int component=0){
		update();
		if(component>0){
			buildComponents();
			assert(getParent(roots[component])==-1);
			return roots[component];
		}
		assert(getParent(0)==-1);//because we always build the tree from 0
		return 0;
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
	void buildComponents(){
		if(!hasComponents){
					hasComponents=true;
					components.clear();
					components.growTo(g.nodes,-1);
					next_component = 0;
					roots.clear();
					//root_list.clear();

					//identify each connected component.
					//use disjoint sets for this, later


					assert(roots.size()>0);


				}
	}

};

#endif
