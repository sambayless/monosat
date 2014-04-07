
#ifndef DYNAMIC_CONNECTIVITY_H_
#define DYNAMIC_CONNECTIVITY_H_

#include "mtl/Vec.h"
#include "mtl/Heap.h"
#include "DynamicGraph.h"
#include "EulerTree.h"
#include "TreapCustom.h"
#include "core/Config.h"
#include "Reach.h"
#include "mtl/Sort.h"
#include "ThorupDynamicConnectivity.h"
using namespace Minisat;



template<class Status,class EdgeStatus=DefaultEdgeStatus>
class DynamicConnectivity:public Reach, public AllPairs,public ConnectedComponents{
public:

	DynamicGraph<EdgeStatus> & g;
	Status & status;
	int last_modification;
	int last_addition;
	int last_deletion;
	int history_qhead;


	int last_history_clear;
	vec<int> sources;
	//the transitive closure of the graph (only for components connected to one of the sources)
	struct ClosureData{
		char reachable:1;
		char changed:1;
		ClosureData():reachable(false),changed(false){}
	};
	vec<vec<ClosureData>> transitive_closure;
	struct Change{
		int sourceNum;
		int node;
	};
	vec<Change> changes;
	vec<int> prev;
	vec<int> path;
	vec<int> u_component;
	vec<int> v_component;
	vec<int> component;
	int INF;

	ThorupDynamicConnectivity t;
	int default_source;
	const int reportPolarity;

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

	DynamicConnectivity(DynamicGraph<EdgeStatus> & graph,Status & status, int _reportPolarity=0 ):g(graph),status(status), last_modification(-1),last_addition(-1),last_deletion(-1),history_qhead(0),last_history_clear(0),INF(0),reportPolarity(_reportPolarity){
		stats_skipped_updates=0;
		stats_full_updates=0;
		stats_full_update_time=0;
		stats_fast_update_time=0;
		default_source=-1;
	}

	void addSource(int s){
		int sourceNum = sources.size();
		sources.push(s);

		transitive_closure.push();
		transitive_closure[sourceNum].growTo(g.nodes);
		transitive_closure[sourceNum][s].reachable=true;

		last_modification=-1;
		last_addition=-1;
		last_deletion=-1;
	}

	void setSource(int s){
		if(!sources.contains(s)){
			addSource(s);
		}
		default_source = s;

	}
	int getSource(){
		return default_source;
	}

	void setNodes(int n){

		INF=g.nodes+1;

		while(t.nNodes()< g.nodes){
			t.addNode();
		}

	}

	void updateEdge(int u,int v,int edgeid, bool add){

							if(add){
								bool already_connected=false;
								if(!t.connected(u,v)){
									already_connected=false;
									u_component.clear();
									v_component.clear();
									for(int i = 0;i<sources.size();i++){

											//then it may be the case that new nodes are connected to one of the sources we are tracking
											if(transitive_closure[i][u].reachable && !transitive_closure[i][v].reachable){
												//ok, explore the new connected component and enclose it
												if(!v_component.size())
													t.getConnectedComponent(v,v_component);
												for(int n:v_component){
													transitive_closure[i][n].reachable=true;
													if(!transitive_closure[i][n].changed){
														transitive_closure[i][n].changed=true;
														changes.push({i,n});
													}
												}
											}else if (transitive_closure[i][v].reachable &&! transitive_closure[i][u].reachable){
												//ok, explore the new connected component and enclose it
												if(!u_component.size())
													t.getConnectedComponent(u,u_component);
												for(int n:u_component){
													transitive_closure[i][n].reachable=true;
													if(!transitive_closure[i][n].changed){
														transitive_closure[i][n].changed=true;
														changes.push({i,n});
													}
												}
											}
									}
								}else{
									already_connected=true;
								}
								//avoiding an extra connected check here by using the unchecked variatn...
								t.setEdgeEnabledUnchecked(u,v,edgeid,already_connected);
							}else{
								if(t.setEdgeEnabled(u,v,edgeid,add)){
									assert(!t.connected(u,v));
									u_component.clear();
									v_component.clear();
									for(int i = 0;i<sources.size();i++){
										int source = sources[i];
											//then it may be the case that new nodes are connected to one of the sources we are tracking
											if(transitive_closure[i][u].reachable){
												assert(transitive_closure[i][v].reachable);
												//ok, explore the new connected component and enclose it
												if(t.connected(source,u)){
													assert(!t.connected(source,v));
													if(!v_component.size())
														t.getConnectedComponent(v,v_component);
													for(int n:v_component){
														transitive_closure[i][n].reachable=false;
														if(!transitive_closure[i][n].changed){
															transitive_closure[i][n].changed=true;
															changes.push({i,n});
														}
													}
												}else{
													assert(t.connected(source,v));
													assert(!t.connected(source,u));
													if(!u_component.size())
														t.getConnectedComponent(u,u_component);
													for(int n:u_component){
														transitive_closure[i][n].reachable=false;
														if(!transitive_closure[i][n].changed){
															transitive_closure[i][n].changed=true;
															changes.push({i,n});
														}
													}
												}
											}
									}
								}
							}
	}

#ifndef NDEBUG
	DisjointSets dbg_sets;
#endif
	void update(){


		static int iteration = 0;
			int local_it = ++iteration ;


			if(last_modification>0 && g.modifications==last_modification){
				stats_skipped_updates++;
				return;
			}
			stats_full_updates++;

			setNodes(g.nodes);


	#ifndef NDEBUG
			dbg_sets.Reset();

			dbg_sets.AddElements(g.nodes);

			for(int i = 0;i<g.all_edges.size();i++){
				if(g.edgeEnabled(i) && g.all_edges[i].id>=0){
					int u = g.all_edges[i].from;
					int v = g.all_edges[i].to;
					dbg_sets.UnionElements(u,v);
				}
			}

	#endif
			changes.clear();
			if(g.historyclears!=last_history_clear){
				last_history_clear=g.historyclears;
				history_qhead=0;
				//start from scratch
				for(int i = 0;i<g.all_edges.size();i++){
					if(g.all_edges[i].id>=0){
						bool add = g.edgeEnabled(i);
						int u =  g.all_edges[i].from;
						int v =  g.all_edges[i].to;
						updateEdge(u,v,i,add);
					}
				}
			}else{
				//incremental/decremental update
				for(;history_qhead<g.history.size();history_qhead++){
					int edgeid = g.history[history_qhead].id;
					bool add = g.history[history_qhead].addition;
					int u =  g.all_edges[edgeid].from;
					int v =  g.all_edges[edgeid].to;
					updateEdge(u,v,edgeid,add);
				}

			}
#ifndef NDEBUG
				for(int i = 0;i<g.edges;i++){
					if(g.all_edges[i].id>=0){
					assert(t.edges[i].edgeID==g.all_edges[i].id);
					assert(t.edges[i].from==g.all_edges[i].from);
					assert(t.edges[i].to==g.all_edges[i].to);
					assert(t.edgeEnabled(i)==g.edgeEnabled(i));
					}
				}
				for(int i = 0;i<g.nodes;i++){
					for(int j = 0;j<g.nodes;j++){
						bool dbg_connected = dbg_sets.FindSet(i)==dbg_sets.FindSet(j);
						assert(dbg_connected == t.connected(i,j));
					}
				}
#endif
			for(Change c:changes){
				int source =sources[ c.sourceNum];
				assert(transitive_closure[c.sourceNum][c.node].changed);
				transitive_closure[c.sourceNum][c.node].changed=false;
				if(source==default_source){
					bool reachable = transitive_closure[c.sourceNum][c.node].reachable;
					if(reachable && reportPolarity>=0){
						assert(t.connected(source,c.node));
						status.setReachable(c.node,true);
					}else if(!reachable && reportPolarity<=0 ){
						assert(!t.connected(source,c.node));
						status.setReachable(c.node,false);
					}
				}
			}
			changes.clear();
#ifndef NDEBUG
			{
				for(vec<ClosureData> & v:transitive_closure)
					for(ClosureData & c:v)
						assert(!c.changed);
			}
#endif
/*

			if (default_source>=0){
				prev.clear();
				prev.growTo(g.nodes,-1);
				//ok, traverse the nodes connected to this component
				component.clear();
				static int iter=0;
				//this is NOT the right way to do this.
				//need to only see check from t!
				t.getConnectedComponent(default_source,component);

				if(reportPolarity>=0){
					for(int v:component){
						status.setReachable(v,true);
					}
				}
				if(reportPolarity<=0){
					sort(component);
					int nextpos = 0;
					for(int v = 0;v<g.nodes;v++){
						if(nextpos<component.size() && v==component[nextpos]){
							nextpos++;
							assert(dbg_sets.FindSet(v)==dbg_sets.FindSet(default_source));
							//status.setReachable(v,true);

						}else{
							assert(dbg_sets.FindSet(v)!=dbg_sets.FindSet(default_source));
							status.setReachable(v,false);

						}
					}
				}
			}
*/

			assert(dbg_sets.NumSets()== t.numComponents());


			last_modification=g.modifications;
			last_deletion = g.deletions;
			last_addition=g.additions;

			history_qhead=g.history.size();
			last_history_clear=g.historyclears;

			//stats_full_update_time+=cpuTime()-startdupdatetime;;
	}

	void dbg_path(int from,int to, vec<int> & path){
		assert(path.size());
		assert(path[0]==from);
		assert(path.last()==to);
		for(int i = 1;i<path.size();i++){
			int v = path[i];
			int u = path[i-1];
			assert(g.hasEdgeUndirected(u,v));

		}
	}

	bool connected(int from, int to){
			update();
			return t.connected(from,to);
		}


		 int numComponents(){
			 update();
			 return t.numComponents();
		 }
		int getComponent(int node){
			update();
			return t.findRoot(node);
		}

		bool connected_unsafe(int from,int to){
			return t.connected(from,to);
		}
		bool connected_unchecked(int from,int to){
			return t.connected(from,to);
		}

		int distance(int from,int to){
			update();
			return t.connected(from,to) ? 0:INF ;
		}
		int distance_unsafe(int from,int to){
			return t.connected(from,to) ? 0:INF ;
		}
		void getPath(int source, int to, vec<int> & path_store){
			return t.getPath(source,to, path_store);
		}


		 bool connected_unsafe(int t){
			 return connected_unsafe(getSource(),t);
		 }
		 bool connected_unchecked(int t){
			 return connected_unchecked(getSource(),t);
		 }
		 bool connected( int t){
			 return connected(getSource(),t);
		 }
		 int distance( int t){
			 return distance(getSource(),t);
		 }
		 int distance_unsafe(int t){
			 return distance_unsafe(getSource(),t);
		 }
		 int previous( int to){
			 update();
			 if(prev[to]<0){

				t.getPath(getSource(),to,path);
			/*	for(int i:path){
					printf("%d->",i);
				}
				printf("\n");*/

				dbg_path(getSource(),to,path);

				prev.clear();
				prev.growTo(g.nodes,-1);
				assert(path[0]==getSource());
				for(int i = 1;i<path.size();i++){
					prev[path[i]]=path[i-1];
				}

			 }
			 assert(prev[to]!=to);
			 assert(prev[to]>=0);
			 return prev[to];
		 }
		 void getPath( int t, vec<int> & path_store){
			 getPath(getSource(),t,  path_store);
		 }
};

#endif
