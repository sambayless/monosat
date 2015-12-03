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

#ifndef STEINER_APPROX_H_
#define STEINER_APPROX_H_

#include <vector>
#include "alg/Heap.h"
#include "mtl/Sort.h"
#include "DynamicGraph.h"
#include "core/Config.h"
#include "MinimumSpanningTree.h"
#include "SteinerTree.h"
#include "alg/DisjointSets.h"
#include <limits>
#include <algorithm>
#include "Dijkstra.h"
#include "Kruskal.h"
#include "Distance.h"

namespace dgl {
template<class TerminalSet, class Status, typename Weight = int>
class SteinerApprox: public SteinerTree<Weight> {
public:
	
	DynamicGraph<Weight> & g;

	TerminalSet & terminals;
	Status & status;

	int last_modification;
	Weight min_weight = 0;
	bool is_disconnected = false;
	int last_addition;
	int last_deletion;
	int history_qhead;

	int last_history_clear;

	Weight INF;

	const int reportPolarity;
	std::vector<bool> in_tree;
	std::vector<int> tree_edges;
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

	SteinerApprox(DynamicGraph<Weight> & graph,  TerminalSet & terminals, Status & _status,
			int _reportPolarity = 0) :
			g(graph),terminals(terminals), status(_status), last_modification(-1), last_addition(-1), last_deletion(
					-1), history_qhead(0), last_history_clear(0), INF(0), reportPolarity(_reportPolarity) {
		
		mod_percentage = 0.2;
		stats_full_updates = 0;
		stats_fast_updates = 0;
		stats_skip_deletes = 0;
		stats_skipped_updates = 0;
		stats_full_update_time = 0;
		stats_fast_update_time = 0;
		stats_num_skipable_deletions = 0;
		stats_fast_failed_updates = 0;
		min_weight = -1;
		
	}
	
	void setNodes(int n) {
		INF = std::numeric_limits<int>::max();
	}
	
	void update() {
		static int iteration = 0;
		int local_it = ++iteration;

		if (g.outfile) {
			fprintf(g.outfile, "m\n");
			fflush(g.outfile);
		}

		if (last_modification > 0 && g.modifications == last_modification) {
			stats_skipped_updates++;
			return;
		}
		stats_full_updates++;
		
		if (last_deletion == g.deletions) {
			stats_num_skipable_deletions++;
		}
		
		setNodes(g.nodes());
		
		in_tree.clear();
		in_tree.resize(g.nodes(), false);
		
		min_weight = 0;
		tree_edges.clear();
		is_disconnected = false;
		if (terminals.numEnabled() > 1) {
			std::vector<Distance<Weight>*> reaches;
			
			//construct the metric closure of GL on the set of terminal nodes.
			//i.e., find the shortest path between each terminal node; then form a graph with one edge for each such path...
			DynamicGraph<Weight> induced;
			
			for (int i = 0; i < g.nodes(); i++) {
				induced.addNode();
				
			}
			
			//This can be made much more efficient...
			for (int i = 0; i < terminals.nodes(); i++) {
				reaches.push_back(new Dijkstra<Weight>(i, g));
			};
			for (int i = 0; i < terminals.nodes(); i++) {
				if (terminals.nodeEnabled(i)) {
					reaches[i]->update();
					//now add in the corresponding edges to the subgraph
					for (int n = 0; n < terminals.nodes(); n++) {
						assert(reaches[i]->connected(n) == reaches[n]->connected(i));
						if (n != i && terminals.nodeEnabled(n)) {
							if (reaches[i]->connected(n)) {
								Weight & dist = reaches[i]->distance(n);
								//temporarily disabled - need to add this back!
								//induced.addEdge(i,n,-1,dist);
							} else {
								is_disconnected = true;
							}
						}
					}
				}
			}
			
			Kruskal<typename MinimumSpanningTree<Weight>::NullStatus, Weight> mst(induced);
			min_weight = 0;
			
			for (int & edgeID : mst.getSpanningTree()) {
				int u = induced.getEdge(edgeID).from;
				int v = induced.getEdge(edgeID).to;
				int p = v;
				std::vector<int> path;
				int first = -1;
				int last = -1;
				assert(reaches[v]->connected(u));
				assert(reaches[u]->connected(v));
				//Find the first and last edges on this path that are already in T
				while (p != u) {
					int pathEdge = reaches[u]->incomingEdge(p);
					if (in_tree[p]) {
						if (first < 0) {
							first = p;
						}
						last = p;
					} else if (first < 0) {
						in_tree[p] = true;
						min_weight += g.getWeight(pathEdge);
						tree_edges.push_back(pathEdge);
					}
					p = reaches[u]->previous(p);
				}
				//now that we have found the last edge on the path that is in the steiner tree, add the subpath u..last to the tree
				p = u;
				while (p != u) {
					assert(in_tree[p]);
					int pathEdge = reaches[u]->incomingEdge(p);
					tree_edges.push_back(pathEdge);
					p = reaches[u]->previous(p);
					assert(!in_tree[p]);
					in_tree[p] = true;
				}
			}
			assert(min_weight <= mst.forestWeight());
		}
		
		if (is_disconnected) {
			status.setMinimumSteinerTree(INF);
		} else {
			status.setMinimumSteinerTree(min_weight);
		}
		
		last_modification = g.modifications;
		last_deletion = g.deletions;
		last_addition = g.additions;
		
		history_qhead = g.historySize();
		last_history_clear = g.historyclears;
		//dbg_drawSteiner();
		assert(dbg_uptodate());
	}
	
	void dbg_drawSteiner() {
#ifndef NDEBUG
		
		printf("digraph{\n");
		for (int i = 0; i < g.nodes(); i++) {
			
			if (terminals.nodeEnabled(i)) {
				
				printf("n%d [fillcolor=blue,style=filled]\n", i);
			} else if (in_tree[i]) {
				printf("n%d [fillcolor=gray,style=filled]\n", i);
			} else {
				printf("n%d\n", i);
			}
			
		}
		
		for (int i = 0; i < g.nodes(); i++) {
			for (int j = 0; j < g.nIncident(i); j++) {
				int id = g.incident(i, j).id;
				int u = g.incident(i, j).node;
				
				const char * s = "black";
				
				if (g.edgeEnabled(id))
					s = "black";
				else
					s = "gray";
				
				if (std::count(tree_edges.begin(), tree_edges.end(), id)) {
					s = "blue";
				}
				
				//printf("n%d -> n%d [label=\"v%d w=%d\",color=\"%s\"]\n", i,u, id,g.getWeight(id), s);
				
			}
		}
		printf("}\n");
#endif
	}
	
	Weight &weight() {
		
		update();
		
		assert(dbg_uptodate());
		if (is_disconnected)
			return INF;
		return min_weight;
	}
	
	bool disconnected() {
		update();
		return is_disconnected;
	}
	
	void getSteinerTree(std::vector<int> & edges) {
		edges = tree_edges;
	}
	
	bool dbg_uptodate() {
#ifndef NDEBUG
		
#endif
		return true;
	}
	;
	
};
}
;
#endif
