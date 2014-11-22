/****************************************************************************************[Solver.h]
The MIT License (MIT)

Copyright (c) 2014, Sam Bayless

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/

#ifndef DFS_CYCLE_H_
#define DFS_CYCLE_H_

#include <vector>
#include "alg/Heap.h"
#include "DynamicGraph.h"
#include "core/Config.h"
#include "Reach.h"

namespace dgl{

class DFSCycle:public Cycle{
public:

	DynamicGraph & g;

	bool directed;
	int last_modification;
	int last_addition;
	int last_deletion;
	int history_qhead;

	int last_history_clear;


	int INF;


	std::vector<int> q;
	std::vector<int> check;
	std::vector<int> path;
	std::vector<int> cycle;

	const int reportPolarity;

	//std::vector<char> old_seen;
	std::vector<bool> seen;
	std::vector<bool> in_q;

//	std::vector<int> changed;

	bool undirected_cycle;
	bool directed_cycle;

	std::vector<int> prev;

public:


	DFSCycle(DynamicGraph & graph,bool _directed=true, int _reportPolarity=0 ):g(graph),directed(_directed), last_modification(-1),last_addition(-1),last_deletion(-1),history_qhead(0),last_history_clear(0),INF(0),reportPolarity(_reportPolarity){
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
		undirected_cycle=false;
		directed_cycle=false;
	}


	void setNodes(int n){
		q.reserve(n);
		check.reserve(n);
		seen.resize(n);
		prev.resize(n);
		INF=g.nodes()+1;
	}

	void update( ){
		static int iteration = 0;
		int local_it = ++iteration ;

		if(last_modification>0 && g.modifications==last_modification){
			stats_skipped_updates++;
			//reportPolarity(undirected_cycle,directed_cycle);
			return;
		}

		if(last_modification>0 &&  last_deletion==g.deletions){
			stats_num_skipable_deletions++;
			//reportPolarity(undirected_cycle,directed_cycle);
			return;
		}

		setNodes(g.nodes());

		stats_full_updates++;
		

		q.clear();
		for(int i = 0;i<g.nodes();i++){
			seen[i]=0;
			in_q[i]=0;
			prev[i]=-1;
		}

		cycle.clear();
		path.clear();
		for(int k = 0;k<g.nodes();k++){//to handle disconnected graph
			if(seen[k])
				continue;
			q.push_back(k);
			while(q.size()){//dfs
				int u = q.back();
				if(in_q[u]){
					q.pop_back();
					path.pop_back();
					in_q[u]=false;
					continue;
				}else{
					in_q[u]=true;

				}

				assert(seen[u]);


				for(int i = 0;i<g.nIncident(u);i++){
					if(!g.edgeEnabled( g.incident(u,i).id))
						continue;
					int v = g.incident(u,i).node;

					if(!seen[v]){
						seen[v]=1;
						prev[v]=u;
						q.push_back(v);
						path.push_back(g.incident(u,i).id);
					}else{
						if(!undirected_cycle){
							//path.copyTo(cycle);

							assert(path.size()==q.size()-1);
							for(int j = 1;j<q.size();j++){
								if(in_q[j]){
									cycle.push_back(path[j-1]);
								}
							}
							cycle.push_back(g.incident(u,i).id);

						}
						undirected_cycle=true;
						if(!directed){
							break;
						}else if(in_q[i]){
							directed_cycle=true;
							cycle.clear();
							assert(path.size()==q.size()-1);
							for(int j = 1;j<q.size();j++){
								if(in_q[j]){
									cycle.push_back(path[j-1]);
								}
							}
							cycle.push_back(g.incident(u,i).id);

							break;
						}
					}
				}
			}
		}
		last_modification=g.modifications;
		last_deletion = g.deletions;
		last_addition=g.additions;

		history_qhead=g.history.size();
		last_history_clear=g.historyclears;
		;
	}


	bool hasDirectedCycle(){
		update();
		return directed_cycle;
	}
	bool hasUndirectedCycle(){
		update();
		return undirected_cycle;
	}

	std::vector<int> & getUndirectedCycle(){
		update();
		return cycle;
	}
	std::vector<int> & getDirectedCycle(){
		update();
		return cycle;
	}
};
};
#endif
