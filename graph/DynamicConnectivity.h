
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
class DynamicConnectivity:public Reach, public AllPairs,public ConnectedComponents{
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

	const int reportPolarity;

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


	DynamicConnectivity(int s,DynamicGraph<EdgeStatus> & graph, int _reportPolarity=0 ):g(graph), last_modification(-1),last_addition(-1),last_deletion(-1),history_qhead(0),last_history_clear(0),source(s),INF(0),reportPolarity(_reportPolarity){

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
				t.clear();
				//start from scratch
				for(int i = 0;i<g.edges;i++){
					if(g.edgeEnabled(i)){
						t.setEdgeEnabled(i,true);
					}
				}

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





};

#endif
