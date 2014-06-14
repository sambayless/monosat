
#ifndef PRIMS_H_
#define PRIMS_H_

#include <vector>
#include "mtl/Heap.h"
#include "DynamicGraph.h"
#include "core/Config.h"
#include "MinimumSpanningTree.h"
#include <limits>



namespace dgl{
template<class Status>
class Prim:public MinimumSpanningTree{
public:

	DynamicGraph & g;
	Status &  status;
	int last_modification;
	int min_weight;
	int last_addition;
	int last_deletion;
	int history_qhead;

	int last_history_clear;

	int INF;


	std::vector<int> q;
	std::vector<int> check;
	const int reportPolarity;

	//std::vector<char> old_seen;
	std::vector<char> seen;;
//	std::vector<int> changed;
	//std::vector<int> edge_weights;
	std::vector<int> keys;
	std::vector<int> parents;
	std::vector<int> parent_edges;
	//int next_component;
	std::vector<int> components;
	std::vector<int> roots;
	std::vector<bool> in_tree;
	int numsets;
    struct VertLt {
        const std::vector<int>&  keys;

        bool operator () (int x, int y) const {
        	return keys[x]<keys[y];
        }
        VertLt(const std::vector<int>&  _key) : keys(_key) { }
    };
    //bool hasComponents;
	Heap<VertLt> Q;

	std::vector<int> mst;
	std::vector<int> prev;

	struct DefaultReachStatus{
		std::vector<bool> stat;
		void setReachable(int u, bool reachable){
			stat.resize(u+1);
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

	Prim(DynamicGraph & graph, Status & _status, int _reportPolarity=0 ):g(graph), status(_status), last_modification(-1),last_addition(-1),last_deletion(-1),history_qhead(0),last_history_clear(0),INF(0),reportPolarity(_reportPolarity),Q(VertLt(keys)){

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
		numsets=0;
	}

	void setNodes(int n){
		q.reserve(n);
		check.reserve(n);
		seen.resize(n);
		prev.resize(n);
		INF=std::numeric_limits<int>::max();
		components.resize(g.nodes);
		parents.resize(n);
		keys.resize(n);
		parent_edges.resize(n);
	}

	void update( ){
		static int iteration = 0;
		int local_it = ++iteration ;
#ifdef RECORD
		if(g.outfile){
			fprintf(g.outfile,"m\n");
			fflush(g.outfile);
		}
#endif
		if(g.modifications==89){
			int a=1;
		}
		if(last_modification>0 && g.modifications==last_modification){
			stats_skipped_updates++;
			return;
		}

		stats_full_updates++;
		double startdupdatetime = rtime(2);
		if(last_deletion==g.deletions){
			stats_num_skipable_deletions++;
		}

		setNodes(g.nodes);

		min_weight=0;

		mst.clear();

		in_tree.clear();in_tree.resize(g.edges);

		for(int i = 0;i<g.nodes;i++){
			parents[i]=-1;
			parent_edges[i]=-1;//not really needed
			components[i]=i;
			int key = INF;
			keys[i]=key;
		}
		int n_seen=0;
		int root=0;
		numsets=0;
		seen.clear();
		seen.resize(g.nodes);
		min_weight=0;
		//This outer loop is to get Prim's to compute the minimum spanning _forest_, in the case that the graph is disconnected.
		while(n_seen<g.nodes){
			for(;root<g.nodes && seen[root];root++);//find the first unseen node
			components[root]=root;
			Q.insert(root);//arbitrary first node to examine
			numsets++;
			while(Q.size()){
				int u = Q.removeMin();
				seen[u]=true;
				n_seen++;
				if(u!=root){
					int parent = parents[u];
					assert(parent!=-1);
					int edgeid = parent_edges[u];
					assert(edgeid!=-1);
					min_weight+=g.getWeight(edgeid);
					mst.push_back(edgeid);
					in_tree[edgeid]=true;
					components[u]=root;
				}
				for(int j = 0;j< g.adjacency_undirected[u].size();j++){
					int edgeid = g.adjacency_undirected[u][j].id;
					if(!g.edgeEnabled(edgeid))
						continue;
					int v = g.adjacency_undirected[u][j].node;
					if(!seen[v]){
						int w = g.getWeight(edgeid);
						if(w<keys[v]){
							parents[v]=u;
							parent_edges[v]=edgeid;
							keys[v]= w;

							Q.update(v);
						}
					}
				}
			}
		}


/*		for(int i = 0;i<g.nodes;i++){
			if(parents[i] ==-1) {
				numsets++;
			}

		}*/

		if(numsets==0){
			assert(min_weight==0);
		}
		status.setMinimumSpanningTree(numsets>1 ? INF: min_weight);
		for(int i = 0;i<in_tree.size();i++){
			//Note: for the tree edge detector, polarity is effectively reversed.
			if(reportPolarity<1 && (!g.edgeEnabled(i) || in_tree[i]) ){
				status.inMinimumSpanningTree(i,true);
			}else if(reportPolarity>-1 && (g.edgeEnabled(i) && ! in_tree[i]) ){
				status.inMinimumSpanningTree(i,false);
			}
		}
		assert(dbg_uptodate());

		last_modification=g.modifications;
		last_deletion = g.deletions;
		last_addition=g.additions;

		history_qhead=g.history.size();
		last_history_clear=g.historyclears;



		stats_full_update_time+=rtime(2)-startdupdatetime;;
	}
	std::vector<int> & getSpanningTree(){
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

	bool edgeInTree(int edgeid){
		update();
		//int u = g.all_edges[edgeid].from;
		//int v = g.all_edges[edgeid].to;
		//return parents[u]==v || parents[v]==u;
		return in_tree[edgeid];
	}
	int numComponents(){
		update();
		//buildComponents();
		return numsets;
	}
	int getComponent(int node){
		update();
		//buildComponents();
		assert(components[node]!=-1);
		return components[node];
	}
	int getRoot(int component=0){
		update();
		return components[component];
	/*	if(component>0){
			buildComponents();
			assert(getParent(roots[component])==-1);
			return roots[component];
		}
		assert(getParent(0)==-1);//because we always build the tree from 0
		return 0;*/
	}

	int weight(){
		update();

		assert(dbg_uptodate());
		if(numsets>1){
			return INF;
		}
		return min_weight;
	}

	int forestWeight(){
		update();
		assert(dbg_uptodate());
		return min_weight;
	}

	bool dbg_uptodate(){
#ifndef NDEBUG
		int sumweight = 0;
		in_tree.resize(g.nEdgeIDs());
		for(int i = 0;i<g.edges;i++){
			if(in_tree[i]){
				sumweight+= g.getWeight(i);
			}
		}
		assert(sumweight ==min_weight || min_weight==INF);

#endif
		return true;
	};
private:
/*	void buildComponents(){
		if(!hasComponents){
					hasComponents=true;
					components.clear();
					components.resize(g.nodes,-1);
					next_component = 0;
					roots.clear();
					//root_list.clear();

					//identify each connected component.
					//use disjoint sets for this, later


					assert(roots.size()>0);


				}
	}*/

};
};
#endif
