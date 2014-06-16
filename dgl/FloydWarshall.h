
#ifndef FLOYD_WARSHALL_H_
#define FLOYD_WARSHALL_H_

#include <vector>
#include "alg/Heap.h"
#include "graph/DynamicGraph.h"
#include "core/Config.h"
#include "AllPairs.h"
#include "mtl/Sort.h"

namespace dgl{

template<class Status>
class FloydWarshall:public AllPairs{
public:

	DynamicGraph & g;
	Status & status;
	int last_modification;
	int last_addition;
	int last_deletion;
	int history_qhead;

	int last_history_clear;


	std::vector<int> sources;
	int INF;

	std::vector<int> order;
//	std::vector<int> q;
//	std::vector<int> check;
	const int reportPolarity;

	std::vector<std::vector<int> > dist;
	std::vector<std::vector<int> >  next;

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

	FloydWarshall(DynamicGraph & graph,Status & _status,  int _reportPolarity=0 ):g(graph),status(_status), last_modification(-1),last_addition(-1),last_deletion(-1),history_qhead(0),last_history_clear(0),INF(0),reportPolarity(0){

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


	void addSource(int s){
		assert(!std::count( sources.begin(),sources.end(), s));
		sources.push_back(s);

		last_modification=-1;
		last_addition=-1;
		last_deletion=-1;
	}



	void setNodes(int n){
		//q.reserve(n);
		//check.reserve(n);
		order.clear();
		for(int i =0;i<n;i++)
			order.push_back(i);
		INF=g.nodes()+1;
		if(dist.size()<n){
			dist.resize(n);

			for(int i =0;i<dist.size();i++)
				dist[i].resize(n);

			next.resize(n);
			for(int i =0;i<dist.size();i++)
				next[i].resize(n);
		}
	}

	struct lt_key{
		std::vector<int> & _dist;
		 bool operator()(int a, int b)const{
			return _dist[a]<_dist[b];
		}
		 lt_key(std::vector<int> & d):_dist(d){};
	};
	void update( ){

		stats_full_updates++;

		if(last_modification>0 && g.modifications==last_modification){
			stats_skipped_updates++;
			return;
		}

		if(last_deletion==g.deletions){
			stats_num_skipable_deletions++;
		}

		setNodes(g.nodes());


		for(int i = 0;i<g.nodes();i++){
			for(int j = 0;j<g.nodes();j++){
				next[i][j]=-1;
				dist[i][j]=INF;
			}
			dist[i][i]=0;
		}

		for(int i = 0;i<g.all_edges.size();i++){
			if(g.edgeEnabled(g.all_edges[i].id)){
				 int u =g.all_edges[i].from;
				 int v =g.all_edges[i].to;
				 if(u!=v)
					 dist[u][v]= 1;
			}
		}


		//for(int l = 0;l<sources.size();l++){
		//	int k = sources[l];
		for(int k = 0;k<g.nodes();k++){
			for(int i = 0;i<g.nodes();i++){
				for (int j = 0;j<g.nodes();j++){
					int d = dist[i][k] + dist[k][j];
					if(d < dist[i][j]){
						dist[i][j] = d;
						next[i][j]=k;
					}
				}
			}
		}
		for(int i = 0;i<sources.size();i++){
			int s = sources[i];
			//sort(order,lt_key(dist[s]));//disabled because it is NOT required
			for(int j = 0;j<order.size();j++){
				int u =j;// order[j];
		/*		if(j<order.size()-1){
					assert(dist[s][u]<=dist[s][order[j+1]]);
				}*/
				//Wrong. This is only required if we are returning learnt clauses that include other reachability lits.
				//it is crucial to return the nodes in order of distance, so that they are enqueued in the correct order in the solver.
				if(dist[s][u]>=INF && reportPolarity<1){
					status.setReachable(s,u,false);
					status.setMininumDistance(s,u,false,INF);
				}else if(dist[s][u]<INF && reportPolarity>-1){
					status.setReachable(s,u,true);
					status.setMininumDistance(s,u,true,dist[s][u]);
				}
			}
		/*	for(int u = 0;u<g.nodes();u++){
				if(dist[s][u]>=INF && reportPolarity<1){
					//status.setReachable(s,u,false);
					//status.setMininumDistance(s,u,false,INF);
				}else if(dist[s][u]<INF && reportPolarity>-1){
					//status.setReachable(s,u,true);
					//status.setMininumDistance(s,u,true,dist[s][u]);
				}
			}*/
		}
		assert(dbg_uptodate());

		last_modification=g.modifications;
		last_deletion = g.deletions;
		last_addition=g.additions;

		history_qhead=g.history.size();
		last_history_clear=g.historyclears;


	}


	void getPath(int from, int to, std::vector<int> & path){
		update();
		path.push_back(from);
		getPath_private(from,to,path);
		assert(path.back()!=to);
		path.push_back(to);
	}
	void getPath_private(int from, int to, std::vector<int> & path){
		assert(dist[from][to]<INF);
		int intermediate = next[from][to];
		if(intermediate>-1){
			getPath_private(from, intermediate, path);
			path.push_back(intermediate);
			getPath_private(intermediate,to,path);

		}

	}
	bool dbg_path(int from,int to){

		return true;
	}
	void drawFull(){
				/*printf("digraph{\n");
				for(int i = 0;i< g.nodes();i++){

					if(seen[i]){
						printf("n%d [fillcolor=blue style=filled]\n", i);
					}else{
						printf("n%d \n", i);
					}


				}

				for(int i = 0;i< g.nodes();i++){
					for(int j =0;j<g.nIncident(u);j++){
					int id  =g.incident(i,j).id;
					int u =  g.incident(i,j).node;
					const char * s = "black";
					if( g.edgeEnabled(id))
						s="blue";
					else
						s="red";



					printf("n%d -> n%d [label=\"v%d\",color=\"%s\"]\n", i,u, id, s);
					}
				}

				printf("}\n");*/
			}
	bool dbg_uptodate(){

		return true;
	}

	bool connected_unsafe(int from,int t){
		return dist[from][t]<INF;
	}
	bool connected_unchecked(int from,int t){
		assert(last_modification==g.modifications);
		return connected_unsafe(from,t);
	}
	bool connected(int from,int t){
		if(last_modification!=g.modifications)
			update();

		assert(dbg_uptodate());

		return dist[from][t]<INF;
	}
	int distance(int from,int t){
		if(connected(from,t))
			return dist[from][t];
		else
			return INF;
	}
	int distance_unsafe(int from,int t){
		if(connected_unsafe(from,t))
			return  dist[from][t];
		else
			return INF;
	}
/*	int previous(int t){
		assert(t>=0 && t<next.size());
		assert(next[t]>=-1 && next[t]<next.size());
		return next[t];
	}*/

};
};
#endif
