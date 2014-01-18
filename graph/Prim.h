
#ifndef PRIMS_H_
#define PRIMS_H_

#include "mtl/Vec.h"
#include "mtl/Heap.h"
#include "DynamicGraph.h"
#include "core/Config.h"
#include "MinimumSpanningTree.h"
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

	int source;
	int INF;


	vec<int> q;
	vec<int> check;
	const int reportPolarity;

	//vec<char> old_seen;
	vec<char> seen;;
//	vec<int> changed;


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


	Prim(int s,DynamicGraph<EdgeStatus> & graph, Status & _status, int _reportPolarity=0 ):g(graph), status(_status), last_modification(-1),last_addition(-1),last_deletion(-1),history_qhead(0),last_history_clear(0),source(s),INF(0),reportPolarity(_reportPolarity){
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

		setNodes(g.nodes);

		min_weight=0;


		q.clear();
		for(int i = 0;i<g.nodes;i++){
			seen[i]=0;
			prev[i]=-1;
		}
		seen[source]=1;
		q.push_(source);
		for (int i = 0;i<q.size();i++){
			int u = q[i];
			assert(seen[u]);
			if(reportPolarity==1)
				status.setReachable(u,true);

			for(int i = 0;i<g.adjacency[u].size();i++){
				if(!g.edgeEnabled( g.adjacency[u][i].id))
					continue;
				int v = g.adjacency[u][i].node;

				if(!seen[v]){
					seen[v]=1;
					prev[v]=u;
					q.push_(v);
				}
			}
		}

		if(reportPolarity<1){
			for(int u = 0;u<g.nodes;u++){
				if(!seen[u]){
					status.setReachable(u,false);
				}else if(reportPolarity==0){
					status.setReachable(u,true);
				}
			}
		}
		assert(dbg_uptodate());

		last_modification=g.modifications;
		last_deletion = g.deletions;
		last_addition=g.additions;

		history_qhead=g.history.size();
		last_history_clear=g.historyclears;



		stats_full_update_time+=cpuTime()-startdupdatetime;;
	}

	bool dbg_mst(){

		return true;
	}
	int weight(){
		if(last_modification!=g.modifications)
			update();

		assert(dbg_uptodate());

		return min_weight;
	}

};

#endif
