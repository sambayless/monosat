
#ifndef SPIRA_PAN_H_
#define SPIRA_PAN_H_

#include "mtl/Vec.h"
#include "mtl/Heap.h"
#include "mtl/Sort.h"
#include "DynamicGraph.h"
#include "core/Config.h"
#include "MinimumSpanningTree.h"
#include <algorithm>
#include <limits>
using namespace Minisat;

template<class Status,class EdgeStatus=DefaultEdgeStatus>
class SpiraPan:public MinimumSpanningTree{
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

	vec<int> mst;
	vec<int> q;
	vec<int> check;
	const int reportPolarity;

	//vec<char> old_seen;
	vec<bool> in_tree;
	vec<bool> keep_in_tree;
	vec<int> parents;
	vec<int> parent_edges;
	vec<int> components_to_visit;
	vec<int> component_weight;
    struct VertLt {
        const vec<int>&  keys;

        bool operator () (int x, int y) const {
        	return keys[x]<keys[y];
        }
        VertLt(const vec<int>&  _key) : keys(_key) { }
    };

	Heap<VertLt> Q;

	vec<bool> seen;

	int num_sets=0;


	vec<int> components;

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

	SpiraPan(DynamicGraph<EdgeStatus> & graph, Status & status, int reportPolarity=0 ):g(graph), status(status), last_modification(-1),last_addition(-1),last_deletion(-1),history_qhead(0),last_history_clear(0),INF(0),reportPolarity(reportPolarity),Q(VertLt(component_weight)){

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
		in_tree.growTo(g.nEdgeIDs());
		seen.growTo(n);
		INF=std::numeric_limits<int>::max();

		parents.growTo(n);
		parent_edges.growTo(n);
	}
/*
 * 	//chin houck insertion
	void insert(int newNode){
		if(g.adjacency_undirected[newNode].size()==0)
			return;
		marked.clear();
		marked.growTo(g.nodes);
		incident_edges.clear();
		incident_edges.growTo(g.nodes,-1);
		for(auto & edge:g.adjacency_undirected[newNode]){
			incident_edges[edge.node]=edge.id;
		}
		insert(g.adjacency_undirected[newNode][0]);
	}

	//Insert a node into an _existing_ minimum spanning tree
	void insert(int r, int z, int & t){
		int m = incident_edges[r];
		marked[r]=true;

		for(auto edge:g.adjacency_undirected[r]){//IF this is a dense graph, then this loop is highly sub-optimal!
			//for each incident edge in the old MST:
			if(in_tree[edge.id] && !marked[ edge.node]){
				insert(edge.node,z,t);
				assert(dbg_is_largest_edge_on_path(m,r,z));
				assert(dbg_is_largest_edge_on_path(t,edge.node,z));
				int k = t;
				if(t<0 || g.weights[edge.id]>g.weights[t])
					k=edge.id;
				int h = t;
				if(t<0 || g.weights[edge.id]<g.weights[t])
					h=edge.id;
				if(h>=0){
					keep_in_tree[h]=true;
				}
				if(m<0 || (k>=0 && g.weights[k]<g.weights[m])){
					m=k;
				}
			}
		}
		t=m;
	}*/

	bool dbg_is_largest_edge_on_path(int edge, int from, int to){
#ifndef NDEBUG

#endif
		return true;
	}

	void dbg_parents(){
#ifndef NDEBUG

#endif
	}

	void update( ){
		static int iteration = 0;
		int local_it = ++iteration ;
#ifdef RECORD
		if(g.outfile && mstalg==MinSpanAlg::ALG_SPIRA_PAN){
			fprintf(g.outfile,"m\n");
			fflush(g.outfile);
		}
#endif
		if(last_modification>0 && g.modifications==last_modification)
					return;
		if(last_modification<=0 || g.changed()){
			INF=g.nodes+1;
			hasParents=false;

			setNodes(g.nodes);
			seen.clear();
			seen.growTo(g.nodes);
			min_weight=0;
			num_sets = g.nodes;
			components.clear();
			for(int i = 0;i<g.nodes;i++)
				components.push(i);
			mst.clear();
			parents.growTo(g.nodes,-1);
			for(int i = 0;i<in_tree.size();i++)
				in_tree[i]=false;
		}

		if(last_history_clear!=g.historyclears){
			history_qhead=0;
			last_history_clear=g.historyclears;
		}

		double startdupdatetime = rtime(2);
		assert(components_to_visit.size()==0);
		for (int i = history_qhead;i<g.history.size();i++){

			int edgeid = g.history[i].id;
			if(g.history[i].addition && g.edgeEnabled(edgeid) ){

				int u = g.all_edges[edgeid].from;
				int v = g.all_edges[edgeid].to;
				int w = g.all_edges[edgeid].weight;
				if(components[u] != components[v]){
					//If u,v are in separate components, then this edge must be in the mst (and we need to fix the component markings)
					in_tree[edgeid]=true;
					num_sets--;
					int higher_component = u;
					int lower_component = v;
					int new_c = components[u];
					int old_c = components[v];
					if(components[v]>components[u]){
						std::swap(higher_component,lower_component);
						std::swap(new_c,old_c);
					}
					//ok, now set every node in the higher component to be in the lower component with a simple dfs.
					//fix the parents at the same time.
					min_weight+=g.weights[edgeid];
					parents[higher_component]=lower_component;
					parent_edges[higher_component]=edgeid;
					q.clear();
					q.push(higher_component);
					while(q.size()){
						int n = q.last(); q.pop();
						for(auto & edge:g.adjacency_undirected[n]){
							if(in_tree[edge.id]){
								assert(g.edgeEnabled(edge.id));
								int t = edge.node;
								if(components[t]==old_c){
									components[t]==new_c;
									q.push(t);
								}
								parents[t]=n;
								parent_edges[t]=edge.id;
							}
						}
					}
					dbg_parents();
				}else{
					assert(components[u]==components[v]);
					if(parents[u]==v || parents[v]==u){
						//If there is already another edge (u,v) that is in the tree, then we at most need to swap that edge out for this one.
						//note that this can only be the case if u is the parent of v or vice versa
						int p_edge = parent_edges[u];
						if(parents[v]==u){
							p_edge = parent_edges[v];
						}
						if(g.getWeight(p_edge)> g.getWeight(edgeid)){
							//then swap these edges without changing anything else.
							in_tree[p_edge]=false;
							in_tree[edgeid]=true;
							int delta = g.getWeight(p_edge)- g.getWeight(edgeid);
							min_weight-=delta;
							if(parents[v]==u){
								assert(parent_edges[v]==p_edge);
								parent_edges[v] = edgeid;
							}else{
								assert(parent_edges[u]==p_edge);
								parent_edges[u] = edgeid;
							}
						}
					}else{

						//otherwise, find the cycle induced by adding this edge into the MST (by walking up the tree to find the LCA - if we are doing many insertions, could we swap this out for tarjan's OLCA?).
						int p = u;

						while(p>-1){
							seen[p]=true;
							p=parents[p];
						}
						int min_edge_weight = INF;
						int min_edge = -1;
						bool edge_on_left=false;//records which branch of the tree rooted at p the edge we are replacing is
						p = v;
						while(true){
							assert(p>-1);//u and v must share a parent, because u and v are in the same connected component and we have already computed the mst.
							if(seen[p]){
								break;
							}else{
								if (g.getWeight( parent_edges[p]) < min_edge_weight){
									min_edge_weight=g.weights[ parent_edges[p]];
									min_edge = parent_edges[p];
								}
								p = parents[p];
							}
						}
						assert(seen[p]);
						 p = u;
						while(p>-1){
							assert(seen[p]);
							seen[p]=false;
							if (g.getWeight( parent_edges[p]) < min_edge_weight){
								min_edge_weight=g.weights[ parent_edges[p]];
								min_edge = parent_edges[p];
								edge_on_left=true;
							}
							p=parents[p];
						}
						if(min_edge_weight<g.getWeight(edgeid)){
							//then swap out that edge with the new edge.
							//this will also require us to repair parent edges in the cycle
							in_tree[edgeid]=true;
							in_tree[min_edge]=false;
							int p;
							int last_p;
							if(edge_on_left){
								p = u;
								last_p = v;
							}else{
								p=v;
								last_p = u;
							}

							int last_edge = edgeid;
							while(parent_edges[p]!=min_edge){
								parents[p]=last_p;
								std::swap(parent_edges[v],last_edge);//re-orient the parents.
								last_p = p;
								p = parents[p];
							}
						}
					}

					dbg_parents();
				}
			}else if (!g.history[i].addition &&  !g.edgeEnabled(edgeid)){
				if(!in_tree[edgeid]){
					//If an edge is disabled that was NOT in the MST, then no update is required.
				}else{
					//this is the 'tricky' case for mst.
					//following Spira & Pan, each removed edge splits the spanning tree into separate components that are MST's for those components.
					//we then basically run Prim's to stitch those components back together, if they can be stitched.

					int u = g.all_edges[edgeid].from;
					int v = g.all_edges[edgeid].to;

					in_tree[edgeid]=false;
					min_weight-=g.getWeight(edgeid);
					num_sets++;

					assert(components[u]==components[v]);
					assert(parents[u]==v || parents[v]==u);
					if(parents[u]==v){
						parents[u]=-1;
						parent_edges[u]=-1;
					}else{
						parents[v]=-1;
						parent_edges[v]=-1;
					}
					//if we want to maintain the guarantee components are always assigned the lowest node number that they contain, we'd need to modify the code below a bit.
					int new_c = u;

					if(new_c == components[v]){
						new_c=v;
					}
					int old_c = components[u];
					components_to_visit.push(new_c);
					assert(new_c!= components[u]);
					assert(q.size()==0);
					//relabel the components of the tree that has been split off.
					q.clear();
					q.push(new_c);
					while(q.size()){
						int n = q.last(); q.pop();
						for(auto & edge:g.adjacency_undirected[n]){
							if(in_tree[edge.id]){
								assert(g.edgeEnabled(edge.id));
								int t = edge.node;
								if(components[t]==old_c){
									components[t]=new_c;
									q.push(t);
								}
							}
						}
					}

				}
			}
		}
		if(components_to_visit.size()){
			//execute Prim's (or Boruvka's?) on the connected components, using a heap.


				for(int i = 0;i<components_to_visit.size();i++){

					int c = components_to_visit[i];
					assert(c>=0);
					assert(components[c]==c);
					//ok, try to connect u's component to v
					//ideally, we'd use the smaller of these two components...
					int smallest_edge=-1;
					int smallest_weight = INF;
					Q.insert(c);

						q.clear();
						q.push(c);
						while(q.size()){
							int n = q.last(); q.pop();
							for(auto & edge:g.adjacency_undirected[n]){
								if(g.edgeEnabled(edge.id)){
									int t = edge.node;
									if(!in_tree[edge.id]){
										assert(g.edgeEnabled(edge.id));
										if(components[t]!= c){
											if(g.getWeight(edge.id)<=smallest_weight){
												smallest_weight = g.getWeight(edge.id);
												smallest_edge = edge.id;
											}
										}
									}else{
										assert(components[t]==components[n]);
										q.push(t);//no need to mark t as visited, as this is a tree, so it is acyclic
									}
								}
							}
						}



					if(smallest_edge>-1){
						num_sets--;
						min_weight+=g.getWeight(smallest_edge);
						assert(!in_tree[smallest_edge]);
						in_tree[smallest_edge]=true;

						//connect these two components together using this edge.
						//this requires us to fix the parent edges.

						int f = g.all_edges[smallest_edge].from;
						int t = g.all_edges[smallest_edge].to;
						assert(components[f]==c|| components[t]==c);

						if(components[f]==c){
							std::swap(f,t);
						}
						int new_c = c;
						int old_c = components[f];
						parents[f]=t;
						components[f]=new_c;
						parent_edges[f]=smallest_edge;
						q.clear();
						q.push(f);
						while(q.size()){
							int n = q.last(); q.pop();
							assert(components[n]==new_c);
							for(auto & edge:g.adjacency_undirected[n]){
								if(in_tree[edge.id]){
									assert(g.edgeEnabled(edge.id));
									int t = edge.node;
									if(components[t]==old_c){
										components[t]==new_c;
										q.push(t);
									}
									parents[t]=n;
									parent_edges[t]=edge.id;
								}
							}
						}

					}

				//}
			}/*else{

//replace the edges one by one, as described in Spira Pan 1975

			}*/
		}
		components_to_visit.clear();
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
	bool edgeInTree(int edgeid){
		update();
		return in_tree[edgeid];
	}
	bool dbg_mst(){

		return true;
	}


	int weight(){
		update();
		assert(dbg_uptodate());
		if(num_sets==1)
			return min_weight;
		else
			return INF;
	}
	 int numComponents(){
		 update();
		 return num_sets;
	 }
	int getComponent(int node){
		update();
		return components[node];
	}
	int getRoot(int component=0){
		update();
		return parents[component];
	}

	bool dbg_uptodate(){
#ifndef NDEBUG
		int sumweight = 0;
		in_tree.growTo(g.nEdgeIDs());
		for(int i = 0;i<g.nEdgeIDs();i++){
			if(in_tree[i]){
				sumweight+= g.getWeight(i);
			}
		}
		assert(sumweight ==min_weight || min_weight==INF);


#endif
		return true;
	};



/*
#ifndef NDEBUG
		int rootcount =0;
		for(int i = 0;i<parents.size();i++){
			if(parents[i]==-1)
				rootcount++;
		}
		assert(rootcount==1);
#endif
*/

};

#endif
