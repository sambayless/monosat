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
#include "Cycle.h"
namespace dgl {
template <typename Weight,bool directed=true, bool undirected=false>
class DFSCycle: public Cycle {
public:

	DynamicGraph<Weight> & g;

	int last_modification=0;
	int last_addition=0;
	int last_deletion=0;
	int history_qhead=0;

	int last_history_clear=0;

	int INF;

	std::vector<int> q;
	std::vector<int> check;
	std::vector<int> path;
	std::vector<int> undirected_cycle;
	std::vector<int> directed_cycle;

	const int reportPolarity;

	//std::vector<char> old_seen;
	std::vector<int> processed;
	std::vector<bool> seen;
	std::vector<bool> ever_seen;


//	std::vector<int> changed;
	
	bool has_undirected_cycle=false;
	bool has_directed_cycle=false;

	std::vector<int> undirected_prev;

	void setNodes(int n) {
		q.reserve(n);
		check.reserve(n);
		seen.resize(n);
		processed.resize(n);
		ever_seen.resize(n);
		undirected_prev.resize(n);
		INF = g.nodes() + 1;
	}
	
	void computeCycles(){

		path.clear();
		q.clear();
		for (int i = 0; i < g.nodes(); i++) {
			seen[i] = false;
			ever_seen[i]=false;
			processed[i]=0;
			undirected_prev[i]=-1;
		}
		
		directed_cycle.clear();
		undirected_cycle.clear();
		if(directed){
			for (int k = 0; k < g.nodes(); k++) { //to handle disconnected graph
				if (ever_seen[k])
					continue;
				q.push_back(k);
				ever_seen[k]=true;
				seen[k]=true;
				while (q.size()) { //dfs
					int u = q.back();
					assert(processed[u]<=g.nIncident(u));
					if (processed[u]==g.nIncident(u)) {
						q.pop_back();
						if(u!=k)
							path.pop_back();
						seen[u]=false;

						continue;
					}


					assert(seen[u]);

					for (;u == q.back() && processed[u] < g.nIncident(u); processed[u]++) {
						int i = processed[u];

						int id = g.incident(u, i).id;
						if (!g.hasEdge(id) || !g.edgeEnabled(id)){
							continue;
						}
						int v = g.incident(u, i).node;
						if (!seen[v]) {
							seen[v] = true;//for directed cycles
							ever_seen[v]=true;
							q.push_back(v);
							path.push_back(id);
							continue;
						}else{
							assert(!has_directed_cycle);
							has_directed_cycle = true;
							directed_cycle.clear();
							directed_cycle.push_back(id);
							assert(path.size() == q.size() - 1);
							for (int j = 1; j < q.size(); j++) {
								if (seen[j]) {
									directed_cycle.push_back(path[j - 1]);
								}
							}
							if(undirected && !has_undirected_cycle){
								//a directed cycle is also an undirected cycle.
								has_undirected_cycle=true;
								undirected_cycle= directed_cycle;
							}
							return;
						}

					}
				}
			}
		}

		if(undirected && !has_undirected_cycle){
			path.clear();
			q.clear();
			for (int i = 0; i < g.nodes(); i++) {

				ever_seen[i]=false;
				processed[i]=0;
				undirected_prev[i]=-1;
			}

			for (int k = 0; k < g.nodes(); k++) { //to handle disconnected graph
				if (ever_seen[k])
					continue;
				q.push_back(k);
				ever_seen[k]=true;

				while (q.size()) { //dfs
					int u = q.back();
					assert(processed[u]<=g.nIncident(u,true));
					if (processed[u]==g.nIncident(u,true)) {
						q.pop_back();
						if(u!=k)
							path.pop_back();
						continue;
					}
					int fromEdge = -1;
					if(u!=k){
						fromEdge=path.back();
						assert((g.getEdge(fromEdge).to==u) || g.getEdge(fromEdge).from==u);
					}
					assert(ever_seen[u]);

					for (;u == q.back() && processed[u] < g.nIncident(u,true); processed[u]++) {
						int i = processed[u];

						int id = g.incident(u, i,true).id;
						if (!g.hasEdge(id) || !g.edgeEnabled(id) || id==fromEdge){
							continue;
						}
						int v = g.incident(u, i,true).node;
						if (!ever_seen[v]) {
							ever_seen[v]=true;
							q.push_back(v);
							path.push_back(id);
							continue;
						}else{
							//a directed cycle is also an undirected cycle.
							has_undirected_cycle=true;
							assert(path.size() == q.size() - 1);
							for (int j = 1; j < q.size(); j++) {
								if (seen[j]) {
									undirected_cycle.push_back(path[j - 1]);
								}
							}
							return;
						}
					}
				}
			}
		}
	}

public:

	DFSCycle(DynamicGraph<Weight> & graph, bool _directed = true, int _reportPolarity = 0) :
			g(graph),  last_modification(-1), last_addition(-1), last_deletion(-1), history_qhead(
					0), last_history_clear(0), INF(0), reportPolarity(_reportPolarity) {

	}


	void update() {
		static int iteration = 0;
		int local_it = ++iteration;

		if (last_modification > 0 && g.modifications == last_modification ) {
			stats_skipped_updates++;
			//reportPolarity(undirected_cycle,directed_cycle);
			return;
		}


		if(last_modification<=0 || g.changed()){
			setNodes(g.nodes());
		}
		has_undirected_cycle=false;
		has_directed_cycle=false;


		stats_full_updates++;

		computeCycles();

		last_modification = g.modifications;
		last_deletion = g.deletions;
		last_addition = g.additions;
		
		history_qhead = g.historySize();
		last_history_clear = g.historyclears;
		;
	}
	
	bool hasDirectedCycle() {
		update();
		return has_directed_cycle;
	}
	bool hasUndirectedCycle() {
		update();
		return has_undirected_cycle;
	}
	
	std::vector<int> & getUndirectedCycle() {
		update();
		return undirected_cycle;
	}
	std::vector<int> & getDirectedCycle() {
		update();
		return directed_cycle;
	}
};
}
;
#endif
