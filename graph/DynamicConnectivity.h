
#ifndef DYNAMIC_CONNECTIVITY_H_
#define DYNAMIC_CONNECTIVITY_H_

#include "mtl/Vec.h"
#include "mtl/Heap.h"
#include "DynamicGraph.h"
#include "EulerTree.h"
#include "TreapCustom.h"
#include "core/Config.h"
#include "Reach.h"
#include "ThorupDynamicConnectivity.h"
using namespace Minisat;



template<class EdgeStatus=DefaultEdgeStatus>
class DynamicConnectivity:public Reach{
public:

	DynamicGraph<EdgeStatus> & g;

	int last_modification;
	int last_addition;
	int last_deletion;
	int history_qhead;

	int last_history_clear;

	int source;
	int INF;

	ThorupDynamicConnectivity t;

	vec<int> check;
	const int reportPolarity;

	vec<char> seen;;


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

	vec<int> trail;
	vec<int> levels;
	vec<int> trail_lim;
public:


	DynamicConnectivity(int s,DynamicGraph<EdgeStatus> & graph, int _reportPolarity=0 ):g(graph), last_modification(-1),last_addition(-1),last_deletion(-1),history_qhead(0),last_history_clear(0),source(s),INF(0),reportPolarity(_reportPolarity){
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

		check.capacity(n);
		seen.growTo(n);
		prev.growTo(n);
		INF=g.nodes+1;

		for(int i = 0;i<n;i++){
			t.addNode();
		}
		for(int i = 0;i<g.edges;i++){
			t.addEdge(g.all_edges[i].id,g.all_edges[i].from,g.all_edges[i].to);
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
			double startdupdatetime = cpuTime();
			if(last_deletion==g.deletions){
				stats_num_skipable_deletions++;
			}


			setNodes(g.nodes);


	#ifndef NDEBUG
			dbg_sets.Reset();
			for(int i = 0;i<g.edges;i++){
				if(g.edgeEnabled(i)){
					int u = g.all_edges[i].from;
					int v = g.all_edges[i].to;
					dbg_sets.UnionElements(u,v);
				}
			}

	#endif

			if(g.historyclears!=last_history_clear){
				last_history_clear=g.historyclears;
				history_qhead=0;
				//start from scratch
			}else{
				//incremental/decremental update
				for(;history_qhead<g.history.size();history_qhead++){
					int edgeid = g.history[history_qhead].id;
					bool add = g.history[history_qhead].addition;
					int u =  g.history[history_qhead].u;
					int v =  g.history[history_qhead].v;
					t.setEdgeEnabled(edgeid,add);
				}
			}


			assert(dbg_sets.NumSets()== t.numComponents());


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
