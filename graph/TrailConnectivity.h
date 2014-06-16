
#ifndef TRAIL_CONNECTIVITY_H_
#define TRAIL_CONNECTIVITY_H_

#include "mtl/Vec.h"
#include "mtl/Heap.h"
#include "DynamicGraph.h"
#include "core/Config.h"
#include "Reach.h"
using namespace Minisat;

template<class Status,class EdgeStatus=DefaultEdgeStatus>
class TrailConnectivity:public Reach{
public:

	DynamicGraph<EdgeStatus> & g;
	Status &  status;
	int last_modification;
	int last_addition;
	int last_deletion;
	int history_qhead;

	int last_history_clear;

	int source;
	int INF;


	vec<int> q;

	const int reportPolarity;

	//vec<char> old_seen;
	vec<char> seen;
	vec<int> prev;

	struct Diff{
		int node;
		bool reachable;
		int prev;
	};

	vec<Reach> trail;
	vec<int> trail_lim;


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


	TrailConnectivity(int s,DynamicGraph<EdgeStatus> & graph, Status & _status, int _reportPolarity=0 ):g(graph), status(_status), last_modification(-1),last_addition(-1),last_deletion(-1),history_qhead(0),last_history_clear(0),source(s),INF(0),reportPolarity(0){
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
	}
	//Connectivity(const Connectivity& d):g(d.g), last_modification(-1),last_addition(-1),last_deletion(-1),history_qhead(0),last_history_clear(0),source(d.source),INF(0),mod_percentage(0.2),stats_full_updates(0),stats_fast_updates(0),stats_skip_deletes(0),stats_skipped_updates(0),stats_full_update_time(0),stats_fast_update_time(0){marked=false;};

	void newDecisionLevel(){
		 trail_lim.push(trail.size());
	}
	void backtrack(int level){
		if(trail_lim.size()>level){
			int stop = trail_lim[level];
			for(int i = trail.size()-1;i>=trail_lim[level];i--){
				Diff r = trail[i];
				if(r.reachable){
					seen[r.node]=false;
					prev[r.node]=-1;
				}else{
					seen[r.node]=true;
					prev[r.node]=r.prev;
				}
			}
			trail.shrink(trail.size()-stop);
			trail_lim.shrink(trail_lim.size()-level);
			assert(trail_lim.size()==level);
		}
	}


	void setSource(int s){
		source = s;
		last_modification=-1;
		last_addition=-1;
		last_deletion=-1;
	}
	int getSource(){
		return source;
	}

	/*void updateFast(){
		stats_fast_updates++;
		double start_time = rtime(2);

		assert(last_deletion==g.deletions);
		last_modification=g.modifications;
		last_addition=g.additions;
		INF=g.nodes()+1;
		seen.growTo(g.nodes());
		prev.growTo(g.nodes());

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
		stats_fast_update_time+=rtime(2)-start_time;
	}*/
/*	vec<int> & getChanged(){
		return changed;
	}
	void clearChanged(){
		changed.clear();
	}*/


	/*
	 * WARNING: THIS FUNDAMENTALLY WONT WORK if there are any cycles in the graph!
	 * inline void delete_update(int to){
		q.clear();
		q.push(to);
		seen[to]=0;
		//Is this really safe? check it very carefully, it could easily be wrong
		while(q.size()){
			int u = q.last();
			q.pop();
			assert(!seen[u]);
			for(int i = 0;i<g.inverted_adjacency[u].size();i++){
				int v = g.inverted_adjacency[u][i];
				if(seen[v]){
					seen[v]=1;
					//Then since to is still seen, we are up to date
					break;
				}
			}
			if(!seen[u]){
				for(int i = 0;i<g.adjacency[u].size();i++){
					int v = g.adjacency[u][i];
					if(seen[v] && prev[v]==to){
						seen[v]=0;
					}
				}
			}else{
#ifdef GRAPH_DEBUG
				for(int i = 0;i<g.adjacency[u].size();i++){
						int v = g.adjacency[u][i];
						assert(seen[v]);
				}
#endif
			}
		}
	}*/

	void setNodes(int n){
		q.capacity(n);
		seen.growTo(n);
		prev.growTo(n);
		INF=g.nodes()+1;
	}

	inline void add_update(int to){
		q.clear();
		q.push_(to);
		//while(q.size()){
		for(int i = 0;i<q.size();i++){
			int u = q[i];
			assert(seen[u]);
			//status.setReachable(u,true);
			//if(!old_seen[u]){
			//	changed.push(u);
			//}
			for(int i = 0;i<g.adjacency[u].size();i++){
				if(!g.edgeEnabled( g.adjacency[u][i].id))
					continue;
				int v = g.adjacency[u][i].to;
				if(!seen[v]){
					seen[v]=1;
					prev[v]=u;
					q.push_(v);
				}
			}
		}
	}
	bool update_additions(){
		double startdupdatetime = rtime(2);

		if(g.historyclears!=last_history_clear){
				last_history_clear=g.historyclears;
				history_qhead=0;
			}

			assert(INF>g.nodes());
			assert(seen.size()>=g.nodes());
			//old_seen.growTo(g.nodes());
			q.clear();

			for(int i = history_qhead;i<g.history.size();i++){
				int from = g.history[i].u;
				int to = g.history[i].v;
				if(g.history[i].addition){
					//incrementally add edge


					if(seen[from] && !seen[to]){
						seen[to]=1;
						prev[to]=from;
						add_update(to);
					}

				}else if (!seen[from] || (seen[to] && seen[prev[to]] &&  prev[to]!=from)){
					//then deleting this edge has no impact on connectivity, so don't need to do anything
				}else{
					stats_fast_failed_updates++;
					stats_fast_update_time+=rtime(2)-startdupdatetime;;
					return false;
				}

			}

			stats_fast_updates++;
			history_qhead = g.history.size();

			assert(dbg_uptodate());

			last_modification=g.modifications;
			last_deletion = g.deletions;
			last_addition=g.additions;

			history_qhead=g.history.size();
			last_history_clear=g.historyclears;

			stats_fast_update_time+=rtime(2)-startdupdatetime;;
			return true;
	}

	//Attempt an incremental update
	/*bool incrementalUpdate(){
		if(g.historyclears!=last_history_clear){
			last_history_clear=g.historyclears;
			history_qhead=0;
		}


		double startdupdatetime = rtime(2);

		q.clear();

		for(int i = history_qhead;i<g.history.size();i++){
			if(g.history[i].addition){
				//incrementally add edge
				int from = g.history[i].u;
				int to = g.history[i].v;

				if(seen[from] && !seen[to]){
					seen[to]=1;
					prev[to]=from;
					add_update(to);
				}

			}else{
				//incrementally delete edge
				int from = g.history[i].u;
				int to = g.history[i].v;
				if(to==0){
					int a =1;
				}
				if(!seen[from] || !seen[to] || prev[to]!=from){
					//then deleting this edge has no impact on connectivity, so don't need to do anything
				}else{

					//IF no other incoming.edges() are seen, then this might be a safe deletion (but we'd need to update any outgoing.edges() that have to as their previous...)

					//Incremental update failed.
					stats_fast_update_time+=rtime(2)-startdupdatetime;;
					return false;
				}

			}

		}
		stats_fast_updates++;
		history_qhead = g.history.size();

		for(int u = 0;u<g.nodes();u++){
			if(last_modification <=0  || (old_seen[u]==1 && seen[u]==0)){
				changed.push(u);
			}
		}

		assert(dbg_uptodate());

		last_modification=g.modifications;
		last_deletion = g.deletions;
		last_addition=g.additions;

		history_qhead=g.history.size();
		last_history_clear=g.historyclears;



		stats_fast_update_time+=rtime(2)-startdupdatetime;;
		return true;
	}*/

	void update( ){
		static int iteration = 0;
		int local_it = ++iteration ;
		stats_full_updates++;
		double startdupdatetime = rtime(2);
		if(last_modification>0 && g.modifications==last_modification){
			stats_skipped_updates++;
			return;
		}

		if(last_deletion==g.deletions){
			stats_num_skipable_deletions++;
		}

		setNodes(g.nodes());

		if(g.historyclears!=last_history_clear){
			last_history_clear=g.historyclears;
			history_qhead=0;
		}else if(opt_inc_graph && last_modification>0 && (g.historyclears <= (last_history_clear+1))){// && (g.history.size()-history_qhead < g.edges()*mod_percentage)){
			if(opt_dec_graph && last_deletion < g.deletions){
				double startddecupdatetime = rtime(2);
				//scan through the deletions and check if any of them matter..
				bool safe=true;
				for(int i = history_qhead;i<g.history.size();i++){
					int from = g.history[i].u;
					int to = g.history[i].v;
					if(g.history[i].addition){
						//safe
					}else if (!seen[from] || (seen[to] && seen[prev[to]] &&  prev[to]!=from)){
						//then deleting this edge has no impact on connectivity, so don't need to do anything
					}else{
						safe= false;
						break;
					}
				}
				if(safe){
					last_deletion=g.deletions;
				}
				stats_fast_update_time+=rtime(2)-startddecupdatetime;;
			}

			if(last_deletion==g.deletions){
				if(update_additions())
					return;
			}
			/*if(incrementalUpdate())
				return;*/
		}



		q.clear();
		for(int i = 0;i<g.nodes();i++){
			seen[i]=0;
			prev[i]=-1;
		}
		seen[source]=1;
		q.push_(source);
		for (int i = 0;i<q.size();i++){
			int u = q[i];
			assert(seen[u]);
			if(reportPolarity>-1)
				status.setReachable(u,true);

			for(int i = 0;i<g.adjacency[u].size();i++){
				if(!g.edgeEnabled( g.adjacency[u][i].id))
					continue;
				int v = g.adjacency[u][i].to;

				if(!seen[v]){
					seen[v]=1;
					prev[v]=u;
					q.push_(v);
				}
			}
		}

		if(reportPolarity<1){
			for(int u = 0;u<g.nodes();u++){
				if(!seen[u]){
					status.setReachable(u,false);
				}
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

	bool dbg_path(int to){
#ifdef DEBUG_DIJKSTRA
		assert(connected(to));
		if(to == source){
			return true;
		}
		int p = prev[to];

		if(p<0){
			return false;
		}
		if(p==to){
			return false;
		}

		return dbg_path(p);


#endif
		return true;
	}

	bool dbg_uptodate(){
#ifdef DEBUG_DIJKSTRA
		if(last_modification<=0)
			return true;
		Dijkstra<EdgeStatus> d(source,g);
		d.update();
		for(int i = 0;i<g.nodes();i++){

			int dbgdist = d.dist[i];
			if(!seen[i])
				assert(dbgdist==d.INF  );
			else
				assert(dbgdist<d.INF);
		}
#endif
		return true;
	}

	bool connected_unsafe(int t)const{
		return t<seen.size() && seen[t];
	}
	bool connected_unchecked(int t)const{
		assert(last_modification==g.modifications);
		return connected_unsafe(t);
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
