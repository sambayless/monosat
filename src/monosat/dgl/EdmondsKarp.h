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

#ifndef EDMONDS_KARP_H
#define EDMONDS_KARP_H

#include "MaxFlow.h"
#include <vector>
#include <algorithm>
#include "DynamicGraph.h"
namespace dgl {
template<typename Weight = int>
class EdmondsKarp: public MaxFlow<Weight> {
public:
	
	std::vector<std::vector<Weight> > F; //(Residual capacity from u to v is C[u,v] - F[u,v])
	std::vector<std::vector<Weight> > C;
	std::vector<int> P;
	std::vector<Weight> M;

	Weight curflow;

	int last_modification;
	int last_deletion;
	int last_addition;

	int history_qhead;
	int last_history_clear;
	int source = -1;
	int sink = -1;
	DynamicGraph<Weight>& g;
	Weight INF;

	std::vector<int> Q;
	std::vector<bool> changed;

	Weight BreadthFirstSearch(int s, int t) {
		for (int i = 0; i < g.nodes(); i++)
			P[i] = -1;
		P[s] = -2;
		Q.clear();
		Q.push_back(s);
		
		//while (Q.size() > 0){
		//int u =Q.back(); Q.pop_back();
		for (int j = 0; j < Q.size(); j++) {
			int u = Q[j];
			
			for (int i = 0; i < g.nIncident(u); i++) {
				if (!g.edgeEnabled(g.incident(u, i).id))
					continue;
				int v = g.incident(u, i).node;
				///(If there is available capacity, and v is not seen before in search)
				Weight c = C[u][v];
				Weight f = F[u][v];
				
				if (((C[u][v] - F[u][v]) > 0) && (P[v] == -1)) {
					P[v] = u;
					Weight b = C[u][v] - F[u][v];
					M[v] = std::min(M[u], b);
					if (v != t)
						Q.push_back(v);
					else
						return M[t];
				}
			}
			//Must also try reverse edges
			for (int i = 0; i < g.nIncoming(u); i++) {
				if (!g.edgeEnabled(g.incoming(u, i).id))
					continue;
				int v = g.incoming(u, i).node;
				///(If there is available capacity, and v is not seen before in search)
				Weight c = C[u][v];
				Weight f = F[u][v];
				
				if (((C[u][v] - F[u][v]) > 0) && (P[v] == -1)) {
					P[v] = u;
					Weight b = C[u][v] - F[u][v];
					M[v] = std::min(M[u], b);        	   //C[u][v] - F[u][v]);
					if (v != t)
						Q.push_back(v);
					else
						return M[t];
				}
			}
		}
		return 0;
		
	}
public:
	EdmondsKarp(DynamicGraph<Weight>& _g, int source = -1, int sink = -1) :
			g(_g), source(source), sink(sink), INF(0xF0F0F0) {
		curflow = -1;
		
		last_modification = -1;
		last_deletion = -1;
		last_addition = -1;
		
		history_qhead = -1;
		last_history_clear = -1;
		setAllEdgeCapacities(1);
	}
	int getSource() const {
		return source;
	}
	int getSink() const {
		return sink;
	}
	void setCapacity(int u, int w, Weight c) {
		if (C.size() < g.nodes()) {
			C.resize(g.nodes());
			for (int i = 0; i < g.nodes(); i++) {
				C[i].resize(g.nodes());
			}
		}
		C[u][w] = c;
	}
	void setAllEdgeCapacities(Weight c) {
		for (int i = 0; i < g.nodes(); i++) {
			for (int j = 0; j < g.nIncident(i); j++) {
				if (!g.edgeEnabled(g.incident(i, j).id))
					continue;
				setCapacity(i, g.incident(i, j).node, c);
			}
		}
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
		if (last_modification > 0 && g.modifications == last_modification) {
			
			return curflow;
		}
		Weight f = 0;
		C.resize(g.nodes());
		F.resize(g.nodes());
		P.resize(g.nodes());
		M.resize(g.nodes());
		
		for (int i = 0; i < g.nodes(); i++) {
			P[i] = -1;
			M[i] = 0;
			F[i].resize(g.nodes());
			C[i].resize(g.nodes());
			for (int j = 0; j < g.nodes(); j++) {
				F[i][j] = 0;
			}
		}
		P[s] = -2;
		M[s] = INF;
		while (true) {
			Weight m = BreadthFirstSearch(s, t);
			
			if (m == 0)
				break;
			
			f = f + m;
			
			int v = t;
			while (v != s) {
				int u = P[v];
				F[u][v] = F[u][v] + m;
				F[v][u] = F[v][u] - m;
				v = u;
			}
		}
		curflow = f;
		num_updates++;
		last_modification = g.modifications;
		last_deletion = g.deletions;
		last_addition = g.additions;
		
		history_qhead = g.historySize();
		last_history_clear = g.historyclears;
		return f;
	}
	const Weight maxFlow(int s, int t, int max_length) {
		Weight f = 0;
		C.resize(g.nodes());
		F.resize(g.nodes());
		P.resize(g.nodes());
		M.resize(g.nodes());
		changed.resize(g.nEdgeIDs());
		for (int i = 0; i < g.nodes(); i++) {
			P[i] = -1;
			M[i] = 0;
			F[i].resize(g.nodes());
			C[i].resize(g.nodes());
			for (int j = 0; j < g.nodes(); j++) {
				F[i][j] = 0;
			}
		}
		P[s] = -2;
		M[s] = INF;
		while (true) {
			int m = BreadthFirstSearch(s, t);
			
			if (m == 0)
				break;
			
			f = f + m;
			if (f > max_length) {
				return f;
			}
			int v = t;
			while (v != s) {
				int u = P[v];
				F[u][v] = F[u][v] + m;
				F[v][u] = F[v][u] - m;
				v = u;
			}
		}
		return f;
	}
	
	std::vector<bool> seen;
	std::vector<bool> visited;
	const Weight minCut(std::vector<MaxFlowEdge> & cut) {
		return minCut(source, sink, cut);
	}
	bool minCut(int s, int t, int max_length, std::vector<MaxFlowEdge> & cut) {
		Weight f = maxFlow(s, t, max_length);
		if (f > max_length) {
			return false;
		}
		//ok, now find the cut
		Q.clear();
		Q.push_back(s);
		seen.clear();
		seen.resize(g.nodes());
		seen[s] = true;
		cut.clear();
		//	visited.clear();
		//visited.resize(g.nodes());
		//	visited[s]=true;
		for (int j = 0; j < Q.size(); j++) {
			int u = Q[j];
			
			for (int i = 0; i < g.nIncident(u); i++) {
				if (!g.edgeEnabled(g.incident(u, i).id))
					continue;
				int v = g.incident(u, i).node;
				int id = g.incident(u, i).id;
				if (C[u][v] - F[u][v] == 0) {
					cut.push_back(MaxFlowEdge { u, v, id });
				} else if (!seen[v]) {
					Q.push_back(v);
					seen[v] = true;
				}
			}
		}
		//Now remove any edges that lead to vertices that we ended up visiting
		int i, j = 0;
		for (i = 0; i < cut.size(); i++) {
			if (!seen[cut[i].v]) {
				cut[j++] = cut[i];
			}
		}
		cut.resize(j);
		return true;
	}
	
	const Weight getEdgeFlow(int edgeid) {
		assert(g.edgeEnabled(edgeid));
		int u = g.all_edges[edgeid].from;
		int v = g.all_edges[edgeid].to;
		return F[u][v];
	}
	const Weight getEdgeCapacity(int id) {
		assert(g.edgeEnabled(id));
		int u = g.all_edges[id].from;
		int v = g.all_edges[id].to;
		return C[u][v];
	}
	
	const Weight getEdgeResidualCapacity(int id) {
		assert(g.edgeEnabled(id));
		int u = g.all_edges[id].from;
		int v = g.all_edges[id].to;
		return C[u][v] - F[u][v];
		
	}
	const Weight minCut(int s, int t, std::vector<MaxFlowEdge> & cut) {
		Weight f = maxFlow(s, t);
		//ok, now find the cut
		Q.clear();
		Q.push_back(s);
		seen.clear();
		seen.resize(g.nodes());
		seen[s] = true;
		
		//explore the residual graph
		for (int j = 0; j < Q.size(); j++) {
			int u = Q[j];
			
			for (int i = 0; i < g.nIncident(u); i++) {
				if (!g.edgeEnabled(g.incident(u, i).id))
					continue;
				int v = g.incident(u, i).node;
				int id = g.incident(u, i).id;
				if (C[u][v] - F[u][v] == 0) {
					cut.push_back(MaxFlowEdge { u, v, id });	//potential element of the cut
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
				if (F[v][u] == 0) {
					
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
			int u = cut[i].u;
			int v = cut[i].v;
			assert(F[u][v] == C[u][v]);
			dbg_sum += F[u][v];
		}
		assert(dbg_sum == f);
#endif
		return f;
	}
	
};

}
;
#endif

