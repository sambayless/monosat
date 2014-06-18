/*
 * Steiner.h
 *
 *  Created on: Jun 17, 2014
 *      Author: sam
 */





#ifndef STEINER_APPROX_H_
#define STEINER_APPROX_H_

#include <vector>
#include "alg/Heap.h"
#include "mtl/Sort.h"
#include "graph/DynamicGraph.h"
#include "core/Config.h"
#include "MinimumSpanningTree.h"
#include "SteinerTree.h"
#include "alg/DisjointSets.h"
#include <limits>
#include <algorithm>
#include "BFS.h"
#include "Kruskal.h"

namespace dgl{
template<class TerminalSet, class Status=SteinerTree::NullStatus>
class SteinerApprox:public SteinerTree{
public:

	DynamicGraph & g;
	TerminalSet & terminals;
	Status &  status;

	int last_modification;
	int min_weight;
	int last_addition;
	int last_deletion;
	int history_qhead;

	int last_history_clear;

	int INF;

	const int reportPolarity;

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

	SteinerApprox(DynamicGraph & graph,TerminalSet & terminals, Status & _status, int _reportPolarity=0 ):g(graph),terminals(terminals), status(_status), last_modification(-1),last_addition(-1),last_deletion(-1),history_qhead(0),last_history_clear(0),INF(0),reportPolarity(_reportPolarity){

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
		INF=std::numeric_limits<int>::max();
	}

	void update( ){
		static int iteration = 0;
		int local_it = ++iteration ;
#ifdef RECORD
		if(g.outfile){
			fprintf(g.outfile,"m\n");
			fflush(g.outfile);
		}
#endif
		if(last_modification>0 && g.modifications==last_modification){
			stats_skipped_updates++;
			return;
		}
		stats_full_updates++;

		if(last_deletion==g.deletions){
			stats_num_skipable_deletions++;
		}


		setNodes(g.nodes());

		min_weight=0;

		std::vector<Reach*> reaches;

		//construct the metric closure of GL on the set of terminal nodes.
		//i.e., find the shortest path between each terminal node; then form a graph with one edge for each such path...
		DynamicGraph induced;
		for(int i = 0;i<g.nodes();i++){
			induced.addNode();
		}
		//This can be made much more efficient...
		for(int i = 0;i<terminals.nodes();i++){
			if (terminals.nodeEnabled(i)){
				reaches.push_back(new BFSReachability<Reach::NullStatus,true>(i,g));
				reaches.back()->update();
				//now add in the corresponding edges to the subgraph
				for(int n = 0;n<terminals.nodes();n++){
					if (n!=i && terminals.nodeEnabled(n)  && reaches.back()->connected(n)){
						int dist = reaches.back()->distance(n);
						induced.addEdge(i,n,-1,dist);
					}
				}
			}
		}

		Kruskal<> mst(induced);
		min_weight = mst.weight();
		status.setMinimumSteinerTree(min_weight);
		last_modification=g.modifications;
		last_deletion = g.deletions;
		last_addition=g.additions;

		history_qhead=g.history.size();
		last_history_clear=g.historyclears;

		assert(dbg_uptodate());
	}



	int weight(){

		update();

		assert(dbg_uptodate());

		return min_weight;
	}


	bool dbg_uptodate(){
#ifndef NDEBUG

#endif
		return true;
	};




};
};
#endif
