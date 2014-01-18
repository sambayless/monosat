
#ifndef CONNECTIVITY_H_
#define CONNECTIVITY_H_

#include "mtl/Vec.h"
#include "mtl/Heap.h"
#include "DynamicGraph.h"
#include "core/Config.h"
#include "Reach.h"
using namespace Minisat;



template<class Status,class EdgeStatus=DefaultEdgeStatus>
class Connectivity:public Reach{
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


	Connectivity(int s,DynamicGraph<EdgeStatus> & graph, Status & _status, int _reportPolarity=0 ):g(graph), status(_status), last_modification(-1),last_addition(-1),last_deletion(-1),history_qhead(0),last_history_clear(0),source(s),INF(0),reportPolarity(_reportPolarity){
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
		check.capacity(n);
		seen.growTo(n);
		prev.growTo(n);
		INF=g.nodes+1;
	}

	inline void add_update(int to, bool update){
		q.clear();
		q.push_(to);
		//while(q.size()){
		for(int i = 0;i<q.size();i++){
			int u = q[i];
			assert(seen[u]);
			if(update)
				status.setReachable(u,seen[u]);
			//status.setReachable(u,true);
			//if(!old_seen[u]){
			//	changed.push(u);
			//}
			for(int i = 0;i<g.adjacency[u].size();i++){
				if(!g.edgeEnabled( g.adjacency[u][i].id))
					continue;
				int v = g.adjacency[u][i].node;
				if(!seen[v]){
					seen[v]=1;
					prev[v]=u;
					if(v==53 && u == 54){
						int a=1;
					}
					q.push_(v);
				}
			}
		}
	}


	inline void delete_update(int to){
		q.clear();
		q.push_(to);
		seen[to]=false;
		prev[to]=-1;
		check.clear();
		check.push_(to);
		//while(q.size()){
		for(int i = 0;i<q.size();i++){
			int u = q[i];
			assert(!seen[u]);

			for(int i = 0;i<g.adjacency[u].size();i++){
/*				if(!g.edgeEnabled( g.adjacency[u][i].id))
					continue;*/
				int v = g.adjacency[u][i].node;
				if(seen[v] && prev[v]==u){
					seen[v]=0;
					prev[v]=-1;
					check.push_(v);
					q.push_(v);

				}
			}
		}

		q.clear();
		for(int i = 0;i<check.size();i++){
			int u = check[i];
			if(!seen[u]){
				for(int i = 0;i<g.inverted_adjacency[u].size();i++){

					if(g.edgeEnabled(g.inverted_adjacency[u][i].id)){
						int from = g.inverted_adjacency[u][i].node;
						int to = u;
						if(seen[from]){
							if(to==231){
								int a =1;
							}
							seen[to]=1;
							prev[to]=from;
							add_update(to,false);
							if(to==53 && from == 54){
												int a=1;
											}
							break;
						}
					}
				}

			}
		}

	}
	bool update_additions(){
		double startdupdatetime = cpuTime();

		if(g.historyclears!=last_history_clear){
				last_history_clear=g.historyclears;
				history_qhead=0;
			}

			assert(INF>g.nodes);
			assert(seen.size()>=g.nodes);
			//old_seen.growTo(g.nodes);
			q.clear();

			for(int i = history_qhead;i<g.history.size();i++){
				int from = g.history[i].u;
				int to = g.history[i].v;
				if(g.history[i].addition){
					//incrementally add edge


					if(seen[from] && !seen[to]){
						seen[to]=1;
						prev[to]=from;
						add_update(to,true);
					}

				}else if (to==source || !seen[from] || (seen[to] && seen[prev[to]] &&  prev[to]!=from)){
					//then deleting this edge has no impact on connectivity, so don't need to do anything
				}else{
					stats_fast_failed_updates++;
					stats_fast_update_time+=cpuTime()-startdupdatetime;;
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

			stats_fast_update_time+=cpuTime()-startdupdatetime;;
			return true;
	}

	bool incrementalUpdate(){
			double startdupdatetime = cpuTime();

			if(g.historyclears!=last_history_clear){
					last_history_clear=g.historyclears;
					history_qhead=0;
				}

				assert(INF>g.nodes);
				assert(seen.size()>=g.nodes);
				//old_seen.growTo(g.nodes);
				q.clear();

				for(int i = history_qhead;i<g.history.size();i++){
					int from = g.history[i].u;
					int to = g.history[i].v;
					if(g.history[i].id==814){
												int a=1;
											}
					if(g.history[i].addition && g.edgeEnabled(g.history[i].id)){
						//incrementally add edge


						if(seen[from] && !seen[to]){
							seen[to]=1;
							prev[to]=from;
							add_update(to,false);
						}

					}else if (!g.history[i].addition && !g.edgeEnabled(g.history[i].id)){

						if (to==source || !seen[from] || (seen[to] && seen[prev[to]] &&  prev[to]!=from)){
							//then deleting this edge has no impact on connectivity, so don't need to do anything
						}else{
							delete_update(to);
						}
					}
				}


				for(int u = 0;u<g.nodes;u++){
					status.setReachable(u,seen[u]);
				}



				stats_fast_updates++;
				history_qhead = g.history.size();

				assert(dbg_uptodate());

				last_modification=g.modifications;
				last_deletion = g.deletions;
				last_addition=g.additions;

				history_qhead=g.history.size();
				last_history_clear=g.historyclears;

				stats_fast_update_time+=cpuTime()-startdupdatetime;;
				return true;
		}
	//Attempt an incremental update
	/*bool incrementalUpdate(){
		if(g.historyclears!=last_history_clear){
			last_history_clear=g.historyclears;
			history_qhead=0;
		}


		double startdupdatetime = cpuTime();

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

					//IF no other incoming edges are seen, then this might be a safe deletion (but we'd need to update any outgoing edges that have to as their previous...)

					//Incremental update failed.
					stats_fast_update_time+=cpuTime()-startdupdatetime;;
					return false;
				}

			}

		}
		stats_fast_updates++;
		history_qhead = g.history.size();

		for(int u = 0;u<g.nodes;u++){
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



		stats_fast_update_time+=cpuTime()-startdupdatetime;;
		return true;
	}*/

	void update( ){
		static int iteration = 0;
		int local_it = ++iteration ;

		if(last_modification>0 && g.modifications==last_modification){
			stats_skipped_updates++;
			return;
		}

		if(last_deletion==g.deletions){
			stats_num_skipable_deletions++;
		}

		setNodes(g.nodes);

		if(g.historyclears!=last_history_clear){
			last_history_clear=g.historyclears;
			history_qhead=0;
		}else if(opt_inc_graph && last_modification>0 && (g.historyclears <= (last_history_clear+1))){// && (g.history.size()-history_qhead < g.edges*mod_percentage)){
			if(opt_dec_graph==2){
				if(incrementalUpdate())
					return;
			}else{
				if(opt_dec_graph==1 && last_deletion < g.deletions){
					double startddecupdatetime = cpuTime();
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
					stats_fast_update_time+=cpuTime()-startddecupdatetime;;
				}

				if(last_deletion==g.deletions){
					if(update_additions())
						return;
				}
			}
			/**/
		}

		stats_full_updates++;
		double startdupdatetime = cpuTime();

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
	void drawFull(){
				printf("digraph{\n");
				for(int i = 0;i< g.nodes;i++){

					if(seen[i]){
						printf("n%d [fillcolor=blue style=filled]\n", i);
					}else{
						printf("n%d \n", i);
					}


				}

				for(int i = 0;i< g.adjacency.size();i++){
					for(int j =0;j<g.adjacency[i].size();j++){
					int id  =g.adjacency[i][j].id;
					int u =  g.adjacency[i][j].node;
					const char * s = "black";
					if( g.edgeEnabled(id))
						s="blue";
					else
						s="red";



					printf("n%d -> n%d [label=\"v%d\",color=\"%s\"]\n", i,u, id, s);
					}
				}

				printf("}\n");
			}
	bool dbg_uptodate(){
#ifdef DEBUG_DIJKSTRA
		if(last_modification<=0)
			return true;
		Dijkstra<EdgeStatus> d(source,g);
		d.update();
		//drawFull();
		for(int i = 0;i<g.nodes;i++){

			int dbgdist = d.dist[i];
			if(!seen[i])
				assert(dbgdist==d.INF  );
			else{

				if(!(dbgdist<d.INF)){
					drawFull();
				}
				assert(dbgdist<d.INF);
			}
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
	int distance(int t){
		if(connected(t))
			return 1;
		else
			return INF;
	}
	int distance_unsafe(int t){
		if(connected_unsafe(t))
			return 1;
		else
			return INF;
	}
	int previous(int t){
		assert(t>=0 && t<prev.size());
		assert(prev[t]>=-1 && prev[t]<prev.size());
		return prev[t];
	}

};

#endif
