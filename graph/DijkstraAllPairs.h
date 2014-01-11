
#ifndef DIJKSTRA_ALLPAIRS_H_
#define DIJKSTRA_ALLPAIRS_H_

#include "mtl/Vec.h"
#include "mtl/Heap.h"
#include "DynamicGraph.h"
#include "core/Config.h"
#include "AllPairs.h"
using namespace Minisat;



template<class Status,class EdgeStatus=DefaultEdgeStatus>
class DijkstraAllPairs:public AllPairs{
public:

	DynamicGraph<EdgeStatus> & g;
	Status & status;
	int last_modification;
	int last_addition;
	int last_deletion;
	int history_qhead;

	int last_history_clear;


	vec<int> sources;
	int INF;



	vec<int> check;
	const int reportPolarity;

	vec<vec<int> > dist;
	vec<vec<int> >  prev;

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
	vec<int> * dist_ptr;
	struct DistCmp{
		vec<int> ** _dist;
		 bool operator()(int a, int b)const{
			return (**_dist)[a]<(**_dist)[b];
		}
		 DistCmp(vec<int> **d):_dist(d){};
	};

	Heap<DistCmp> q;

public:


	DijkstraAllPairs(DynamicGraph<EdgeStatus> & graph,Status & _status,  int _reportPolarity=0 ):g(graph),status(_status), last_modification(-1),last_addition(-1),last_deletion(-1),history_qhead(0),last_history_clear(0),INF(0),reportPolarity(0),q(DistCmp(&dist_ptr)){
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


	void addSource(int s){
		assert(!sources.contains(s));
		sources.push(s);

		last_modification=-1;
		last_addition=-1;
		last_deletion=-1;
	}



	void setNodes(int n){

		check.capacity(n);

		INF=g.nodes+1;
		if(dist.size()<n){
			dist.growTo(n);

			for(int i =0;i<dist.size();i++)
				dist[i].growTo(n);

			prev.growTo(n);
			for(int i =0;i<dist.size();i++)
				prev[i].growTo(n);
		}
	}

	void update( ){
		static int iteration = 0;
		int local_it = ++iteration ;
		stats_full_updates++;

		if(last_modification>0 && g.modifications==last_modification){
			stats_skipped_updates++;
			return;
		}

		if(last_deletion==g.deletions){
			stats_num_skipable_deletions++;
		}
		INF=g.nodes+1;
		stats_full_updates++;
		double startdupdatetime = cpuTime();


		setNodes(g.nodes);

		q.clear();
		for(int i = 0;i<g.nodes;i++){
			for(int j = 0;j<g.nodes;j++){
				prev[i][j]=-1;
				dist[i][j]=INF;
			}
			dist[i][i]=0;
		}






		for(int s = 0;s<sources.size();s++){
			int source = sources[s];
			q.clear();
			//distcmp._dist=&dist[source];
			dist_ptr=&dist[source];
			dist[source][source]=0;
			q.insert(source);
			while(q.size()){
				int u = q.peakMin();
				if(dist[source][u]==INF)
					break;

				q.removeMin();
				for(int i = 0;i<g.adjacency[u].size();i++){
					if(!g.edgeEnabled( g.adjacency[u][i].id))
						continue;
					int v = g.adjacency[u][i].node;
					int alt = dist[source][u]+ 1;
					if(alt<dist[source][v]){
						dist[source][v]=alt;
						prev[source][v]=u;
						if(!q.inHeap(v))
							q.insert(v);
						else
							q.decrease(v);
					}
				}
			}
		}

		for(int i = 0;i<sources.size();i++){
			int s = sources[i];
			for(int u = 0;u<g.nodes;u++){
				if(dist[s][u]>=INF && reportPolarity<1){
					status.setReachable(s,u,false);
					status.setMininumDistance(s,u,false,INF);
				}else if(dist[s][u]<INF && reportPolarity>-1){
					status.setReachable(s,u,true);
					status.setMininumDistance(s,u,true,dist[s][u]);
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

	void getPath(int from, int to, vec<int> & path){
		update();
		assert(dist[from][to]<INF);
		int d = dist[from][to];
	/*	int intermediate = prev[from][to];
		if(intermediate>-1){
		getPath(from, intermediate, path);
		path.push(intermediate);
		getPath(intermediate,to,path);
		}*/
		path.clear();
		path.push(to);
		while(prev[from][to]!=-1){
			int p = prev[from][to];
			path.push(p);
			to=p;
		}

		for(int i = 0;i<path.size();i++){
			if(i>=path.size()-i-1)
				break;
			int t = path[i];
			int e = path[path.size()-i-1];
			path[i]=e;
			path[path.size()-i-1]=t;
		}
	}

	bool dbg_path(int from,int to){

		return true;
	}
	void drawFull(){
				/*printf("digraph{\n");
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

				printf("}\n");*/
			}
	bool dbg_uptodate(){

		return true;
	}

	bool connected_unsafe(int from,int t)const{
		return dist[from][t]<INF;
	}
	bool connected_unchecked(int from,int t)const{
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
		assert(t>=0 && t<prev.size());
		assert(prev[t]>=-1 && prev[t]<prev.size());
		return prev[t];
	}*/

};

#endif
