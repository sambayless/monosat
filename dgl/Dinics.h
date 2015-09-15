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

#ifndef DINICS_H
#define DINICS_H

#include "DynamicGraph.h"
#include "MaxFlow.h"
#include <vector>
#include "core/Config.h"
#include "EdmondsKarpAdj.h"
#include <algorithm>
#include <climits>

namespace dgl {
template<typename Weight = int>
class Dinitz: public MaxFlow<Weight> {
	
public:
	
	std::vector<Weight> F;

	struct LocalEdge {
		int from;
		int id;
		bool backward = false;
		LocalEdge(int from = -1, int id = -1, bool backward = false) :
				from(from), id(id), backward(backward) {
			
		}
	};
	Weight curflow;
	int last_modification;
	int last_deletion;
	int last_addition;
	bool opt_dinics_recursive = false;
	int history_qhead;
	int last_history_clear;
	std::vector<LocalEdge> prev;
	std::vector<Weight> M;
	std::vector<int> dist;
	std::vector<int> pos; //position in the combined forward and backward adjacency list of each node in the DFS.
	std::vector<bool> changed;
	DynamicGraph<Weight>& g;


	int src=-1;
	int dst=-1;
	int source = -1;
	int sink = -1;
	Weight INF;
	std::vector<int> Q;
	long stats_augmenting_rounds = 0;
	long stats_rounds = 0;
#ifdef DEBUG_MAXFLOW
	std::vector<int> dbg_pos;
	EdmondsKarpAdj<Capacity,Weight> ek;
#endif
	
public:
	Dinitz(DynamicGraph<Weight> & _g,int source = -1, int sink = -1) :
			g(_g),  source(source), sink(sink), INF(0xF0F0F0)
#ifdef DEBUG_MAXFLOW
	,ek(_g,cap,source,sink)
#endif
	{
		curflow = 0;
		last_modification = -1;
		last_deletion = -1;
		last_addition = -1;
		
		history_qhead = -1;
		last_history_clear = -1;
		//setAllEdgeCapacities(1);
	}
	int getSource() const {
		return source;
	}
	int getSink() const {
		return sink;
	}
	void printStats() {
		printf("Dinics :\n");
		
		printf("Rounds: %ld, Augmenting Rounds: %ld\n", stats_rounds, stats_augmenting_rounds);
	}
	
	void setCapacity(int u, int w, Weight c) {
		//C.resize(g.edges());
		//C[ ]=c;
		
	}
	void setAllEdgeCapacities(Weight c) {
		
	}
	void dbg_print_graph(int from, int to) {
#ifndef NDEBUG
		return;
		static int it = 0;
		if (++it == 6) {
			int a = 1;
		}
		printf("Graph %d\n", it);
		printf("digraph{\n");
		for (int i = 0; i < g.nodes(); i++) {
			if (i == from) {
				printf("n%d [label=\"From\", style=filled, fillcolor=blue]\n", i);
			} else if (i == to) {
				printf("n%d [label=\"To\", style=filled, fillcolor=red]\n", i);
			} else
				printf("n%d\n", i);
		}
		
		for (int i = 0; i < g.edges(); i++) {
			if (g.edgeEnabled(i)) {
				auto & e = g.getEdge(i);
				const char * s = "black";
				if (dist[e.to] == dist[e.from] + 1) {
					s = "blue";
				}
				/*if(value(e.v)==l_True)
				 s="blue";
				 else if (value(e.v)==l_False)
				 s="red";*/
				std::cout << "n" << e.from << " -> n" << e.to << " [label=\"" << i << ": " << F[i] << "/" << g.getWeight(i)
						<< "\" color=\"" << s << "\"]\n";
			}
		}
		
		printf("}\n");
#endif
	}
	bool buildLevelGraph(int src, int dst) {
		dist.clear();
		dist.resize(g.nodes(), -1);
		dist[src] = 0;
		Q.push_back(src);
		//Build the level graph using a simple BFS
		for (int i = 0; i < Q.size(); i++) {
			int u = Q[i];
			for (int j = 0; j < g.nIncident(u); j++) {
				int edgeID = g.incident(u, j).id;
				if (!g.edgeEnabled(edgeID))
					continue;
				int v = g.incident(u, j).node;
				if (dist[v] < 0 && F[edgeID] < g.getWeight(edgeID)) {
					dist[v] = dist[u] + 1;
					Q.push_back(v);
				}
			}
			for (int j = 0; j < g.nIncoming(u); j++) {
				int edgeID = g.incoming(u, j).id;
				if (!g.edgeEnabled(edgeID))
					continue;
				int v = g.incoming(u, j).node;
				//this is a backward edge, so it has capacity exactly if the forward edge has flow
				if (dist[v] < 0 && F[edgeID]>0) {
					dist[v] = dist[u] + 1;
					Q.push_back(v);
				}
			}
		}
		Q.clear();
		return dist[dst] >= 0;
	}
	
	Weight findAugmentingPath(int u) {
		Weight m = 0;
		assert(Q.size() == 0);
		
		Q.push_back(src);
		M[src] = INT_MAX;    	//this isn't safe for other types than int, fix this...
		while (Q.size()) {
			int u = Q.back();
			if (u == dst)
				return M[u];
			bool found = false;
			for (; pos[u] < g.nIncident(u); pos[u]++) {
				//int edgeID = g.adjacency[u][pos[u]].id;
				int edgeID = g.incident(u, pos[u]).id;
				if (!g.edgeEnabled(edgeID))
					continue;
				int v = g.incident(u, pos[u]).node;
				if (dist[v] == dist[u] + 1 && F[edgeID] < g.getWeight(edgeID)) {
					//printf("%d\n",edgeID);
					found = true;
					Weight c = g.getWeight(edgeID) - F[edgeID];
					M[v] = std::min(M[u], c);
					prev[v] = LocalEdge(u, edgeID, false);
					if (v == dst) {
						//M[v] = min(M[u], g.getWeight(edgeID) - F[id]);
						m = M[dst];
						assert(Q.back() == u);
						Q.pop_back();
						while (Q.size()) {
							pos[Q.back()]--;
							Q.pop_back();
						}
						//pos[u]++;
						break;
					} else {
						Q.push_back(v);
						
						pos[u]++;
						break;
					}
				}
			}
			if (!found) {
				for (; pos[u] - g.nIncident(u) < g.nIncoming(u); pos[u]++) {
					//int edgeID = g.inverted_adjacency[u][pos[u]-g.nIncident(u)].id;
					int edgeID = g.incoming(u, pos[u] - g.nIncident(u)).id;
					if (!g.edgeEnabled(edgeID))
						continue;
					int v = g.incoming(u, pos[u] - g.nIncident(u)).node;
					//these are backwards edges, which have capacity exactly if the forward edge has non-zero flow
					if (dist[v] == dist[u] + 1 && F[edgeID]>0) {
						//printf("-%d\n",edgeID);
						found = true;
						M[v] = std::min(M[u], F[edgeID]);
						prev[v] = LocalEdge(u, edgeID, true);    	//this is a backward edge
						if (v == dst) {
							m = M[dst];
							assert(Q.back() == u);
							while (Q.size()) {
								pos[Q.back()]--;
								Q.pop_back();
							}
							break;
						} else {
							Q.push_back(v);
							
							pos[u]++;
							break;
						}
					}
				}
			}
			if (!found) {
				Q.pop_back();
			}
		}
		Q.clear();
		if (m > 0) {
			//we found an augmenting flow, so update all the edge flows correspondingly.
			int v = dst;
			while (v != src) {
				int u = prev[v].from;
				int id = prev[v].id;
				if (prev[v].backward) {
					F[id] = F[id] - m;
				} else
					F[id] = F[id] + m;
				v = u;
			}
		}
		return m;
	}
	Weight dbg_findAugmentingPath_recursive(int u, Weight f) {
#ifdef DEBUG_MAXFLOW
		if (u == dst)
		return f;

		for (;dbg_pos[u]<g.nIncident(u);dbg_pos[u]++) {
			int edgeID = g.incident(u,dbg_pos[u]).id;
			if(!g.edgeEnabled(edgeID))
			continue;
			int v = g.incident(u,dbg_pos[u]).node;
			if (dist[v] == dist[u] + 1 && F[edgeID] < g.getWeight(edgeID)) {
				int df = dbg_findAugmentingPath_recursive(v, min(f, g.getWeight(edgeID) - F[edgeID]));
				if (df > 0) {
					//F[edgeID] += df;
					return df;
				}
			}
		}
		for (;dbg_pos[u]-g.nIncident(u) <g.nIncoming(u);dbg_pos[u]++) {
			int edgeID = g.incoming(u,dbg_pos[u]-g.nIncident(u)).id;
			if(!g.edgeEnabled(edgeID))
			continue;
			int v = g.incoming(u,dbg_pos[u]-g.nIncident(u)).node;
			//these are backwards edges, which have capacity exactly if the forward edge has non-zero flow
			if (dist[v] == dist[u] + 1 && F[edgeID]>0) {
				int df = dbg_findAugmentingPath_recursive(v, std::min(f, F[edgeID]));
				if (df > 0) {
					//F[edgeID] -= df;
					return df;
				}
			}
		}
#endif
		return 0;
	}
	Weight findAugmentingPath_recursive(int u, Weight f) {
		if (u == dst)
			return f;
		
		for (; pos[u] < g.nIncident(u); pos[u]++) {
			//int edgeID = g.adjacency[u][pos[u]].id;
			int edgeID = g.incident(u, pos[u]).id;
			if (!g.edgeEnabled(edgeID))
				continue;
			int v = g.incident(u, pos[u]).node;
			if (dist[v] == dist[u] + 1 && F[edgeID] < g.getWeight(edgeID)) {
				//printf("%d\n",edgeID);
				Weight c = g.getWeight(edgeID) - F[edgeID];
				Weight df = findAugmentingPath_recursive(v, std::min(f, c));
				if (df > 0) {
					F[edgeID] += df;
					return df;
				}
			}
		}
		
		for (; pos[u] - g.nIncident(u) < g.nIncoming(u); pos[u]++) {
			//int edgeID = g.inverted_adjacency[u][pos[u]-g.nIncident(u)].id;
			int edgeID = g.incoming(u, pos[u] - g.nIncident(u)).id;
			if (!g.edgeEnabled(edgeID))
				continue;
			int v = g.incoming(u, pos[u] - g.nIncident(u)).node;
			//these are backwards edges, which have capacity exactly if the forward edge has non-zero flow
			if (dist[v] == dist[u] + 1 && F[edgeID]>0) {
				//printf("-%d\n",edgeID);
				Weight df = findAugmentingPath_recursive(v, std::min(f, F[edgeID]));
				if (df > 0) {
					F[edgeID] -= df;
					return df;
				}
			}
		}
		return 0;
	}
	long num_updates = 0;
	int numUpdates() const {
		return num_updates;
	}
	const Weight update() {
		return maxFlow(source, sink);
	}
	std::vector<int> changed_edges;
	std::vector<int> & getChangedEdges() {
		return changed_edges;
	}
	void clearChangedEdges() {
		for (int edgeID : changed_edges) {
			assert(changed[edgeID]);
			changed[edgeID] = false;
		}
		changed_edges.clear();
	}
	
private:
	
	void markChanged(int edgeID) {
		if (!changed[edgeID]) {
			changed[edgeID] = true;
			changed_edges.push_back(edgeID);
		}
	}
public:
	void setSource(int s) {
		if (source == s) {
			return;
		}
		source = s;
		last_modification = g.modifications - 1;
	}
	void setSink(int t) {
		if (sink == t) {
			return;
		}
		sink = t;
		last_modification = g.modifications - 1;
	}
	const Weight maxFlow(int s, int t) {
		Weight f = 0;
		if (g.outfile) {
			fprintf(g.outfile, "f %d %d\n", s, t);
			fflush(g.outfile);
		}

		if (last_modification > 0 && g.modifications == last_modification) {
			return curflow;
		}
		
		src = s;
		dst = t;
		F.clear();
		F.resize(g.edges());
		dist.clear();
		dist.resize(g.nodes());
		M.resize(g.nodes());
		prev.resize(g.nodes());
		changed.resize(g.nEdgeIDs());
		f = 0;
		dbg_print_graph(s, t);
		while (buildLevelGraph(s, t)) {
			dbg_print_graph(s, t);
			stats_rounds++;
			pos.clear();
			pos.resize(g.nodes());
#ifdef DEBUG_MAXFLOW
			dbg_pos.clear();dbg_pos.resize(g.nodes());
#endif
			if (opt_dinics_recursive) {
				while (Weight delta = findAugmentingPath_recursive(s, INF)) {
					stats_augmenting_rounds++;
					f += delta;
					dbg_print_graph(s, t);
				}
			} else {
				//int expect = dbg_findAugmentingPath_recursive(s,INT_MAX);
				while (Weight delta = findAugmentingPath(s)) {
					//assert(delta==expect);
					f += delta;
					stats_augmenting_rounds++;
					dbg_print_graph(s, t);
					//expect = dbg_findAugmentingPath_recursive(s,INT_MAX);
				}
			}
		}
		
#ifdef DEBUG_MAXFLOW
		Weight expected_flow =ek.maxFlow(s,t);
		assert(f==expected_flow);
#endif
		
		//dbg_print_graph(s,t);
		curflow = f;
		num_updates++;
		last_modification = g.modifications;
		last_deletion = g.deletions;
		last_addition = g.additions;
		
		history_qhead = g.historySize();
		last_history_clear = g.historyclears;
		return f;
	}
	
	std::vector<bool> seen;
	std::vector<bool> visited;
	const Weight minCut(std::vector<MaxFlowEdge> & cut) {
		return minCut(source, sink, cut);
	}
	const Weight minCut(int s, int t, std::vector<MaxFlowEdge> & cut) {
		const Weight f = maxFlow(s, t);
		//ok, now find the cut
		Q.clear();
		Q.push_back(s);
		seen.clear();
		seen.resize(g.nodes());
		seen[s] = true;
		cut.clear();
		/*		if(f==0)
		 return 0;*/
		//explore the residual graph
		for (int j = 0; j < Q.size(); j++) {
			int u = Q[j];
			
			for (int i = 0; i < g.nIncident(u); i++) {
				if (!g.edgeEnabled(g.incident(u, i).id))
					continue;
				int v = g.incident(u, i).node;
				int id = g.incident(u, i).id;
				if (g.getWeight(id) - F[id] == 0) {
					cut.push_back(MaxFlowEdge { u, v, id });    	//potential element of the cut
				} else if (!seen[v]) {
					Q.push_back(v);
					seen[v] = true;
				}
			}
			for (int i = 0; i < g.nIncoming(u); i++) {
				if (!g.edgeEnabled(g.incoming(u, i).id))
					continue;
				int v = g.incoming(u, i).node;
				int id = g.incoming(u, i).id;
				if (F[id] == 0) {
					
				} else if (!seen[v]) {
					Q.push_back(v);
					seen[v] = true;
				}
			}
		}
		//Now keep only the edges from a seen vertex to an unseen vertex
		int i, j = 0;
		for (i = 0; i < cut.size(); i++) {
			if (!seen[cut[i].v] && seen[cut[i].u]) {
				cut[j++] = cut[i];
			}
		}
		cut.resize(j);
#ifndef NDEBUG
		Weight dbg_sum = 0;
		for (int i = 0; i < cut.size(); i++) {
			int id = cut[i].id;
			assert(F[id] == g.getWeight(id));
			dbg_sum += F[id];
		}
		assert(dbg_sum == f);
#endif
		return f;
	}
	const Weight getEdgeCapacity(int id) {
		assert(g.edgeEnabled(id));
		return g.getWeight(id);
	}
	const Weight getEdgeFlow(int id) {
		assert(g.edgeEnabled(id));
		return F[id];    	// reserve(id);
	}
	const Weight getEdgeResidualCapacity(int id) {
		assert(g.edgeEnabled(id));
		return g.getWeight(id) - F[id];    	// reserve(id);
	}
};
}
;
#endif

