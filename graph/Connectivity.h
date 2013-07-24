
#ifndef CONNECTIVITY_H_
#define CONNECTIVITY_H_

#include "mtl/Vec.h"
#include "mtl/Heap.h"
#include "DynamicGraph.h"
using namespace Minisat;

/*
class GraphListener{
	void addEdge(int u, int v, int mod);
	void removeEdge(int u, int v, int mod);
};*/


class Connectivity{
public:
	DynamicGraph & g;
	int last_modification;
	int last_addition;
	int last_deletion;
	int addition_qhead;
	int deletion_qhead;
	int lastaddlist;
	int lastdellist;
	int source;
	int INF;

	vec<int> q;


	bool marked;
	vec<char> old_seen;
	vec<char> seen;;
	vec<int> changed;


	vec<int> prev;



public:
	int stats_full_updates;
	int stats_fast_updates;
	int stats_skip_deletes;
	int stats_skipped_updates;

	double stats_full_update_time;
	double stats_fast_update_time;

	Connectivity(int s,DynamicGraph & graph):g(graph), last_modification(-1),last_addition(-1),last_deletion(-1),addition_qhead(0),deletion_qhead(0),lastaddlist(0),lastdellist(0),source(s),INF(0),marked(false),stats_full_updates(0),stats_fast_updates(0),stats_skip_deletes(0),stats_skipped_updates(0),stats_full_update_time(0),stats_fast_update_time(0){	}
	Connectivity(const Connectivity& d):g(d.g), last_modification(-1),last_addition(-1),last_deletion(-1),addition_qhead(0),deletion_qhead(0),lastaddlist(0),lastdellist(0),source(d.source),INF(0),marked(false),stats_full_updates(0),stats_fast_updates(0),stats_skip_deletes(0),stats_skipped_updates(0),stats_full_update_time(0),stats_fast_update_time(0){};


	void setSource(int s){
		source = s;
		last_modification=-1;
		last_addition=-1;
		last_deletion=-1;
	}
	int getSource(){
		return source;
	}

	void updateFast(){
		stats_fast_updates++;
		double start_time = cpuTime();

		assert(last_deletion==g.deletions);
		last_modification=g.modifications;
		last_addition=g.additions;
		INF=g.nodes+1;
		seen.growTo(g.nodes);
		prev.growTo(g.nodes);

		if(lastaddlist!=g.addlistclears){
			addition_qhead=0;
			lastaddlist=g.addlistclears;
		}
		int start = q.size();
		//ok, now check if any of the added edges allow for new connectivity
		for (int i = addition_qhead;i<g.addition_list.size();i++){
			int u=g.addition_list[i].u;
			int v=g.addition_list[i].v;

			if(!seen[v]){
				q.push(v);
				seen[v]=1;
				prev[v]=u;
			}
		}
		addition_qhead=g.addition_list.size();

		for(int i = start;i<q.size();i++){
			int u = q[i];
			assert(seen[u]);
			for(int i = 0;i<g.adjacency[u].size();i++){
				int v = g.adjacency[u][i];

				if(!seen[v]){
					//this was changed
					changed.push(v);
					seen[v]=1;
					prev[v]=u;
					q.push(v);
				}
			}
		}
		stats_fast_update_time+=cpuTime()-start_time;
	}
	vec<int> & getChanged(){
		return changed;
	}
	void clearChanged(){
		changed.clear();
	}
	void update( ){
		static int iteration = 0;
		int local_it = ++iteration ;
		if(local_it==17513){
			int a =1;
		}



		stats_full_updates++;
		double startdupdatetime = cpuTime();

		INF=g.nodes+1;
		seen.growTo(g.nodes);
		prev.growTo(g.nodes);
		old_seen.growTo(g.nodes);
		q.clear();
		for(int i = 0;i<g.nodes;i++){
			old_seen[i]=last_modification > 0 ? seen[i]:0;//this won't work properly if we added nodes...
			seen[i]=0;
			prev[i]=-1;
		}
		seen[source]=1;
		q.push(source);
		for (int i = 0;i<q.size();i++){
			int u = q[i];
			assert(seen[u]);
			if(!old_seen[u]){
				changed.push(u);
			}

			for(int i = 0;i<g.adjacency[u].size();i++){
				int v = g.adjacency[u][i];
				if(!seen[v]){
					seen[v]=1;
					prev[v]=u;
					q.push(v);
				}
			}
		}

		for(int u = 0;u<g.nodes;u++){
			if(last_modification <=0  || (old_seen[u]==1 && seen[u]==0)){
				changed.push(u);
			}
		}

		assert(dbg_uptodate());

		last_modification=g.modifications;
		last_deletion = g.deletions;
		last_addition=g.additions;

		addition_qhead=g.addition_list.size();
		lastaddlist=g.addlistclears;

		deletion_qhead=g.deletion_list.size();
		lastdellist=g.dellistclears;

		stats_full_update_time+=cpuTime()-startdupdatetime;;
	}

	bool dbg_uptodate(){
#ifdef DEBUG_DIJKSTRA
		if(last_modification<=0)
			return true;
		Dijkstra d(source,g);
		d.update();
		for(int i = 0;i<g.nodes;i++){

			int dbgdist = d.dist[i];
			if(!seen[i])
				assert(dbgdist==d.INF  );
			else
				assert(dbgdist<d.INF);
		}
#endif
		return true;
	}

	bool connected_unsafe(int t){
		return t<seen.size() && seen[t];
	}
	bool connected(int t){
		if(last_modification!=g.modifications)
			update();

		assert(dbg_uptodate());

		return seen[t];
	}

	int previous(int t){
		return prev[t];
	}

};

#endif
