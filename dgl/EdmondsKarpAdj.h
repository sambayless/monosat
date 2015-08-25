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

#ifndef EDMONDS_KARP_ADJ_H
#define EDMONDS_KARP_ADJ_H

#include "DynamicGraph.h"
#include "MaxFlow.h"
#include <vector>
#include <algorithm>
#include "core/Config.h"
#include "dgl/EdmondsKarp.h"

namespace dgl {
template<typename Weight = int>
class EdmondsKarpAdj: public MaxFlow<Weight> {
	
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
	int last_s = -1;
	int last_t = -1;

	int last_modification=-1;
	int last_deletion=0;
	int last_addition=0;

	int history_qhead=0;
	int last_history_clear=0;
	std::vector<LocalEdge> prev;
	std::vector<Weight> M;
	std::vector<bool> changed;
	DynamicGraph<Weight>& g;

	int source = -1;
	int sink = -1;
	Weight INF;
	/*
	 #ifdef DEBUG_MAXFLOW
	 EdmondsKarp<Weight> ek;
	 #endif
	 */

	std::vector<int> Q;

	Weight BreadthFirstSearch(int s, int t) {
		for (int i = 0; i < g.nodes(); i++)
			prev[i].from = -1;
		prev[s].from = -2;
		Q.clear();
		Q.push_back(s);
		
		for (int j = 0; j < Q.size(); j++) {
			int u = Q[j];
			
			for (int i = 0; i < g.nIncident(u); i++) {
				if (!g.edgeEnabled(g.incident(u, i).id))
					continue;
				int id = g.incident(u, i).id;
				int v = g.incident(u, i).node;
				///(If there is available capacity, and v is not seen before in search)
				
				Weight &f = F[id];
				Weight c = g.getWeight(id);
				
				//  int fr = F[id];
				if (((c - F[id]) > 0) && (prev[v].from == -1)) {
					prev[v] = LocalEdge(u, id, false);
					Weight b = c - F[id];
					M[v] = std::min(M[u], b);
					if (v != t)
						Q.push_back(v);
					else
						return M[t];
				}
			}
			
			for (int i = 0; i < g.nIncoming(u); i++) {
				int id = g.incoming(u, i).id;
				if (!g.edgeEnabled(id))
					continue;
				
				int v = g.incoming(u, i).node;
				
				Weight f = 0;
				Weight& c = F[id];
				
				if (((c - f) > 0) && (prev[v].from == -1)) {
					prev[v] = LocalEdge(u, id, true);
					Weight b = c - f;
					M[v] = std::min(M[u], b);
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
	EdmondsKarpAdj(DynamicGraph<Weight>& _g, int source = -1, int sink = -1) :
			g(_g),  source(source), sink(sink), INF(0xF0F0F0)
	/*
	 #ifdef DEBUG_MAXFLOW
	 ,ek(_g,source,sink)
	 #endif
	 */
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
	void setCapacity(int u, int w, Weight c) {
		//C.resize(g.edges());
		//C[ ]=c;
		
	}
	void setAllEdgeCapacities(Weight c) {
		
	}
	void dbg_print_graph(int from, int to) {
#ifndef NDEBUG
/*		return;
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
				std::cout << "n" << e.from << " -> n" << e.to << " [label=\"" << i << ": " << F[i] << "/" << g.getWeight(i)
						<< "\" color=\"" << s << "\"]\n";
				//printf("n%d -> n%d [label=\"%d: %d/%d\",color=\"%s\"]\n", e.from,e.to, i, F[i],capacity[i] , s);
			}
		}
		
		printf("}\n");*/
#endif
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
	const Weight maxFlow() {
		return this->maxFlow(source, sink);
	}
	const Weight maxFlow(int s, int t) {
		Weight f = 0;

		if (g.outfile) {
			fprintf(g.outfile, "f %d %d\n", s, t);
			fflush(g.outfile);
		}

		
		if (last_modification > 0 && g.modifications == last_modification && s==last_s && t==last_t) {

			return curflow;
		}
		last_s =s;
		last_t = t;

		changed.resize(g.nEdgeIDs());
		F.resize(g.edges());
		for(int edgeID = 0;edgeID<g.edges();edgeID++){
			//can do much better than this...
			if(g.hasEdge(edgeID) && F[edgeID]>0){
				if (F[edgeID]>0){
					markChanged(edgeID);
				}
			}
			F[edgeID]=0;
		}


		prev.resize(g.nodes());
		M.resize(g.nodes());
		
		for (int i = 0; i < g.nodes(); i++) {
			prev[i].from = -1;
			M[i] = 0;
		}
		prev[s].from = -2;
		M[s] = INF;
		while (true) {
			Weight m = BreadthFirstSearch(s, t);
			
			if (m == 0)
				break;
			
			f = f + m;
			
			int v = t;
			while (v != s) {
				int u = prev[v].from;
				int id = prev[v].id;
				if (prev[v].backward) {
					F[id] = F[id] - m;
				} else
					F[id] = F[id] + m;
				markChanged(id);
				v = u;
			}
			
		}

		//dbg_print_graph(s,t);
		curflow = f;
		num_updates++;
		last_modification = g.modifications;
		last_deletion = g.deletions;
		last_addition = g.additions;
		dbg_print_graph(s, t);
		history_qhead = g.historySize();
		last_history_clear = g.historyclears;
		assert(f>=0);
		return f;
	}
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
		/*	if(f==0)
		 return 0;*/
		dbg_print_graph(s, t);
		//explore the residual graph
		for (int j = 0; j < Q.size(); j++) {
			int u = Q[j];
			
			for (int i = 0; i < g.nIncident(u); i++) {
				if (!g.edgeEnabled(g.incident(u, i).id))
					continue;
				int v = g.incident(u, i).node;
				int id = g.incident(u, i).id;
				if ((g.getWeight(id) - F[id] == 0) && (F[id]>0)) {
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

		return g.getWeight(id);
	}
	const Weight getEdgeFlow(int id) {
		if(!g.edgeEnabled(id)){
			assert(F[id]==0);
		}
		return F[id];    	// reserve(id);
	}
	const Weight getEdgeResidualCapacity(int id) {
		if(!g.edgeEnabled(id)){
				assert(F[id]==0);
			}
		return g.getWeight(id) - F[id];    	// reserve(id);
	}
};
}
;
#endif

