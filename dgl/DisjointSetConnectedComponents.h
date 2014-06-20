/*
 * DisjointSetConnectedComponents.h
 *
 *  Created on: 2014-01-18
 *      Author: sam
 */

#ifndef DISJOINTSETSCONNECTEDCOMPONENTS_H_
#define DISJOINTSETSCONNECTEDCOMPONENTS_H_



#include <vector>
#include "alg/Heap.h"
#include "graph/DynamicGraph.h"
#include "core/Config.h"
#include "ConnectedComponents.h"
#include "alg/DisjointSets.h"
#include <limits>
namespace dgl{
template<class Status=ConnectedComponents::NullConnectedComponentsStatus>
class DisjointSetsConnectedComponents:public ConnectedComponents{
public:

	DynamicGraph & g;
	Status &  status;
	int last_modification;

	int last_addition;
	int last_deletion;
	int history_qhead;

	int last_history_clear;
	bool hasParents=false;
	int INF;
	DisjointSets sets;

	std::vector<int> q;
	std::vector<int> check;
	const int reportPolarity;
	struct ConnectCheck{
		int u;
		int v;
	};
	std::vector<ConnectCheck> connectChecks;
	//std::vector<char> old_seen;

	//stats

	int stats_full_updates=0;
	int stats_fast_updates=0;
	int stats_fast_failed_updates=0;
	int stats_skip_deletes=0;
	int stats_skipped_updates=0;
	int stats_num_skipable_deletions=0;
	double mod_percentage=0.2;

	double stats_full_update_time=0;
	double stats_fast_update_time=0;


public:
	DisjointSetsConnectedComponents(DynamicGraph & graph, Status & _status, int _reportPolarity=0 ):g(graph), status(_status), last_modification(-1),last_addition(-1),last_deletion(-1),history_qhead(0),last_history_clear(0),INF(0),reportPolarity(_reportPolarity){

	}

	DisjointSetsConnectedComponents(DynamicGraph & graph,  int _reportPolarity=0 ):g(graph), status(nullConnectedComponentsStatus), last_modification(-1),last_addition(-1),last_deletion(-1),history_qhead(0),last_history_clear(0),INF(0),reportPolarity(_reportPolarity){

	}

	void setNodes(int n){
		q.reserve(n);
		check.reserve(n);

		INF=std::numeric_limits<int>::max();
		sets.AddElements(n);

	}

	void addConnectedCheck(int u, int v){
		connectChecks.push_back({u,v});
	}

	void update( ){
		static int iteration = 0;
		int local_it = ++iteration ;

		if(last_modification>0 && g.modifications==last_modification){
			stats_skipped_updates++;
			return;
		}
		stats_full_updates++;

		if(last_deletion==g.deletions){
			stats_num_skipable_deletions++;
		}
		hasParents=false;
		sets.Reset();
		setNodes(g.nodes());

		for(int i = 0;i<g.edges();i++){
			if(g.edgeEnabled(i)){
				int u = g.all_edges[i].from;
				int v = g.all_edges[i].to;
				sets.UnionElements(u,v);
			}
		}

		status.setComponents(sets.NumSets());

		for(auto c:connectChecks){
			int u = c.u;
			int v = c.v;
			bool connected = sets.FindSet(u) == sets.FindSet(v);
			if(reportPolarity>=0 && connected){
				status.setConnected(u,v,true);
			}else if (reportPolarity<=0 && !connected){
				status.setConnected(u,v,false);
			}
		}


		last_modification=g.modifications;
		last_deletion = g.deletions;
		last_addition=g.additions;

		history_qhead=g.history.size();
		last_history_clear=g.historyclears;
	}

	bool connected(int from, int to){
		update();
		return sets.FindSet(from)==sets.FindSet(to);
	}


	 int numComponents(){
		 update();
		 return sets.NumSets();
	 }
	int getComponent(int node){
		update();
		return sets.FindSet(node);
	}
	int getElement(int set){
		update();
		return sets.GetElement(set);
	}

	bool dbg_uptodate(){
		return true;
	};


};
};

#endif /* DISJOINTSETCONNECTEDCOMPONENTS_H_ */
