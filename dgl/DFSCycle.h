
#ifndef DFS_CYCLE_H_
#define DFS_CYCLE_H_

#include <vector>
#include "mtl/Heap.h"
#include "DynamicGraph.h"
#include "core/Config.h"
#include "Reach.h"


template<class EdgeStatus=DefaultEdgeStatus>
class DFSCycle:public Cycle{
public:

	DynamicGraph<EdgeStatus> & g;

	bool directed;
	int last_modification;
	int last_addition;
	int last_deletion;
	int history_qhead;

	int last_history_clear;


	int INF;


	vec<int> q;
	vec<int> check;
	vec<int> path;
	vec<int> cycle;

	const int reportPolarity;

	//vec<char> old_seen;
	vec<bool> seen;
	vec<bool> in_q;

//	vec<int> changed;

	bool undirected_cycle;
	bool directed_cycle;

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


	DFSCycle(DynamicGraph<EdgeStatus> & graph,bool _directed=true, int _reportPolarity=0 ):g(graph),directed(_directed), last_modification(-1),last_addition(-1),last_deletion(-1),history_qhead(0),last_history_clear(0),INF(0),reportPolarity(_reportPolarity){
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
		undirected_cycle=false;
		directed_cycle=false;
	}


	void setNodes(int n){
		q.capacity(n);
		check.capacity(n);
		seen.growTo(n);
		prev.growTo(n);
		INF=g.nodes+1;
	}

	void update( ){
		static int iteration = 0;
		int local_it = ++iteration ;

		if(last_modification>0 && g.modifications==last_modification){
			stats_skipped_updates++;
			//reportPolarity(undirected_cycle,directed_cycle);
			return;
		}

		if(last_modification>0 &&  last_deletion==g.deletions){
			stats_num_skipable_deletions++;
			//reportPolarity(undirected_cycle,directed_cycle);
			return;
		}

		setNodes(g.nodes);

		stats_full_updates++;
		double startdupdatetime = rtime(2);

		q.clear();
		for(int i = 0;i<g.nodes;i++){
			seen[i]=0;
			in_q[i]=0;
			prev[i]=-1;
		}

		cycle.clear();
		path.clear();
		for(int k = 0;k<g.nodes;k++){//to handle disconnected graph
			if(seen[k])
				continue;
			q.push_(k);
			while(q.size()){//dfs
				int u = q.last();
				if(in_q[u]){
					q.pop();
					path.pop();
					in_q[u]=false;
					continue;
				}else{
					in_q[u]=true;

				}

				assert(seen[u]);


				for(int i = 0;i<g.adjacency[u].size();i++){
					if(!g.edgeEnabled( g.adjacency[u][i].id))
						continue;
					int v = g.adjacency[u][i].node;

					if(!seen[v]){
						seen[v]=1;
						prev[v]=u;
						q.push_(v);
						path.push(g.adjacency[u][i].id);
					}else{
						if(!undirected_cycle){
							//path.copyTo(cycle);

							assert(path.size()==q.size()-1);
							for(int j = 1;j<q.size();j++){
								if(in_q[j]){
									cycle.push(path[j-1]);
								}
							}
							cycle.push(g.adjacency[u][i].id);

						}
						undirected_cycle=true;
						if(!directed){
							break;
						}else if(in_q[i]){
							directed_cycle=true;
							cycle.clear();
							assert(path.size()==q.size()-1);
							for(int j = 1;j<q.size();j++){
								if(in_q[j]){
									cycle.push(path[j-1]);
								}
							}
							cycle.push(g.adjacency[u][i].id);

							break;
						}
					}
				}
			}
		}
		last_modification=g.modifications;
		last_deletion = g.deletions;
		last_addition=g.additions;

		history_qhead=g.history.size();
		last_history_clear=g.historyclears;
		stats_full_update_time+=rtime(2)-startdupdatetime;;
	}


	bool hasDirectedCycle(){
		update();
		return directed_cycle;
	}
	bool hasUndirectedCycle(){
		update();
		return undirected_cycle;
	}

	vec<int> & getUndirectedCycle(){
		update();
		return cycle;
	}
	vec<int> & getDirectedCycle(){
		update();
		return cycle;
	}
};

#endif
