
#ifndef DISTANCE_H_
#define DISTANCE_H_

#include "mtl/Vec.h"
#include "mtl/Heap.h"
#include "DynamicGraph.h"
#include "core/Config.h"
#include "Reach.h"
namespace Minisat{

struct NullEdgeStatus{
	void setReachable(int u, bool reachable){

	}
	bool isReachable(int u) const{
		return false;
	}

	void setMininumDistance(int u, bool reachable, int distance){

	}
};
static NullEdgeStatus nullEdgeStatus;
/**
 * Detect connectivity within a number of steps in unweighted, directed graphs
 */
template<class Status,class EdgeStatus=DefaultEdgeStatus, bool undirected=false>
class Distance:public Reach{
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
	vec<int> dist;
//	vec<int> changed;


	vec<int> prev;

	//stats

	int stats_full_updates;
	int stats_fast_updates;
	int stats_fast_failed_updates;
	int stats_skip_deletes;
	int stats_skipped_updates;
	int stats_num_skipable_deletions;
	double mod_percentage;

	double stats_full_update_time;
	double stats_fast_update_time;

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


	Distance(int s,DynamicGraph<EdgeStatus> & graph, Status & _status, int _reportPolarity=0 ):g(graph), status(_status), last_modification(-1),last_addition(-1),last_deletion(-1),history_qhead(0),last_history_clear(0),source(s),INF(0),reportPolarity(_reportPolarity){

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

	void setNodes(int n){
		q.capacity(n);
		dist.growTo(n);
		prev.growTo(n);
		INF=g.nodes+1;
	}


	void update( ){
		static int iteration = 0;
		int local_it = ++iteration ;

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

		if(g.historyclears!=last_history_clear){
			last_history_clear=g.historyclears;
			history_qhead=0;
		}


		q.clear();
		for(int i = 0;i<g.nodes;i++){
			dist[i]=INF;
			prev[i]=-1;
		}
		auto & adjacency = undirected? g.adjacency_undirected:g.adjacency;
		dist[source]=0;
		q.push_(source);
		for (int i = 0;i<q.size();i++){
			int u = q[i];
			assert(dist[u]<INF);
			//if(reportPolarity>-1)
				status.setMininumDistance(u,true,dist[u]);
			int d = dist[u];
			for(int i = 0;i<adjacency[u].size();i++){
				if(!g.edgeEnabled(adjacency[u][i].id))
					continue;
				int edgeID = adjacency[u][i].id;
				int v = adjacency[u][i].node;
				int dv = dist[v];
				if(dist[v]>=INF){
					dist[v]=d+1;
					prev[v]=edgeID;
					q.push_(v);
				}else{
					assert(dist[v]<=d+1);
				}
			}
		}

		//if(reportPolarity<1){
			for(int u = 0;u<g.nodes;u++){
				if(dist[u]>=INF){
					status.setMininumDistance(u,dist[u]<INF,dist[u]);
				}
			}
		//}
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
		int p = previous(to);

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

					if(dist[i]<INF){
						printf("n%d [label=\"n%d %d \" fillcolor=blue style=filled]\n", i,i,dist[i]);
					}else{
						printf("n%d \n", i);
					}


				}
				auto & adjacency = undirected? g.adjacency_undirected:g.adjacency;
				for(int i = 0;i< adjacency.size();i++){
					for(int j =0;j<adjacency[i].size();j++){
					int id  =adjacency[i][j].id;
					int u =  adjacency[i][j].node;
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
		Dijkstra<EdgeStatus,undirected> d(source,g);
		d.update();
		//drawFull();

			for(int i = 0;i<g.nodes;i++){
				int distance = dist[i];
				int dbgdist = d.dist[i];
				assert(distance==dbgdist);
			}
#endif
		return true;
	}

	bool connected_unsafe(int t){
		return t<dist.size() && dist[t]<INF;
	}
	bool connected_unchecked(int t){
		assert(last_modification==g.modifications);
		return connected_unsafe(t);
	}
	bool connected(int t){
		if(last_modification!=g.modifications)
			update();

		assert(dbg_uptodate());

		return dist[t]<INF;
	}

	int distance(int t){
		if(connected(t)){
			return dist[t];
		}else{
			return INF;
		}
	}
	int distance_unsafe(int t){
		if(connected_unsafe(t))
			return dist[t];
		else
			return INF;
	}
	int incomingEdge(int t){

		assert(t>=0 && t<prev.size());
		assert(prev[t]>=-1 );
		return prev[t];
	}
	int previous(int t){

		if(incomingEdge(t)<0)
			return -1;
		if (undirected && g.all_edges[incomingEdge(t)].from==t){
			return g.all_edges[incomingEdge(t)].to;
		}
		assert(g.all_edges[incomingEdge(t)].to==t);
		return g.all_edges[incomingEdge(t)].from;
	}

};
};
#endif
