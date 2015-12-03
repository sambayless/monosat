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

#ifndef RAMAL_REPS_H_
#define RAMAL_REPS_H_

#include <dgl/alg/Heap.h>
#include <dgl/Dijkstra.h>
#include <dgl/Distance.h>
#include <dgl/DynamicGraph.h>
#include <dgl/Reach.h>
//#include "core/Config.h"
//#include <algorithm>
#include <cassert>
#include <cstdio>
//#include <exception>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>


namespace dgl {
template<typename Weight = int, class Status = typename Distance<Weight>::NullStatus>
class RamalReps: public Distance<Weight>, public DynamicGraphAlgorithm {
public:
	static bool ever_warned_about_zero_weights;
	DynamicGraph<Weight> & g;
	std::vector<Weight> & weights;
	Status & status;
	int reportPolarity;
	bool reportDistance;
	int last_modification;
	int last_addition;
	int last_deletion;
	int history_qhead;

	int last_history_clear;

	int source;
	Weight INF;

	std::vector<Weight> old_dist;
	std::vector<int> changed;
	std::vector<bool> node_changed;
	std::vector<Weight> dist;
	std::vector<int> prev;
	struct DistCmp {
		std::vector<Weight> & _dist;
		bool operator()(int a, int b) const {
			return _dist[a] < _dist[b];
		}
		DistCmp(std::vector<Weight> & d) :
				_dist(d) {
		}
		;
	};
	Heap<DistCmp> q;

	std::vector<int> edgeInShortestPathGraph;
	std::vector<int> delta;
	std::vector<int> changeset;
	int alg_id;

	struct LocalDistanceStatus {
		RamalReps & outer;

		void setReachable(int u, bool reachable){

		}
		bool isReachable(int u) const {
			return false;
		}

		void setMininumDistance(int u, bool reachable, Weight& distance){
			if(reachable){
				if(outer.dist[u]!=distance){
					outer.dist[u]=distance;
					if(!outer.node_changed[u]){
						outer.node_changed[u]=true;
						outer.changed.push_back(u);
					}
				}
			}else{
				if(outer.dist[u]!=outer.INF){
					outer.dist[u]=outer.INF;
					if(!outer.node_changed[u]){
						outer.node_changed[u]=true;
						outer.changed.push_back(u);
					}
				}
			}
		}
		LocalDistanceStatus(RamalReps & _outer) :
			outer(_outer) {
		}
	} local_distance_status;
	Dijkstra<Weight,LocalDistanceStatus> dijkstras;
	bool has_zero_weights=false;
public:
	
	long stats_full_updates=0;
	long stats_fast_updates=0;
	long stats_fast_failed_updates=0;
	long stats_skip_deletes=0;
	long stats_skipped_updates=0;
	long stats_num_skipable_deletions=0;
	double mod_percentage=0;

	double stats_full_update_time=0;
	double stats_fast_update_time=0;
	RamalReps(int s, DynamicGraph<Weight> & graph,Status & status, int reportPolarity = 0,
			bool reportDistance = false) :
			g(graph), weights(g.getWeights()), status(status), reportPolarity(reportPolarity), reportDistance(reportDistance), last_modification(
					-1), last_addition(-1), last_deletion(-1), history_qhead(0), last_history_clear(0), source(s), INF(
					0), q(DistCmp(dist)),local_distance_status(*this),dijkstras(s,graph,local_distance_status,reportPolarity) {
		
		mod_percentage = 0.2;
		alg_id=g.addDynamicAlgorithm(this);
	}
	//Dijkstra(const Dijkstra& d):g(d.g), last_modification(-1),last_addition(-1),last_deletion(-1),history_qhead(0),last_history_clear(0),source(d.source),INF(0),q(DistCmp(dist)),stats_full_updates(0),stats_fast_updates(0),stats_skip_deletes(0),stats_skipped_updates(0),stats_full_update_time(0),stats_fast_update_time(0){marked=false;};
	
	void setSource(int s) {
		source = s;
		last_modification = -1;
		last_addition = -1;
		last_deletion = -1;
	}
	int getSource() {
		return source;
	}
	
	std::vector<int> & getChanged() {
		return changed;
	}
	void clearChanged() {
		changed.clear();
	}
	
	void drawFull() {
		
	}
	
	void dbg_delta() {
#ifdef DEBUG_RAMAL
		dbg_delta_lite();
		assert(delta.size() == g.nodes());
		
		for (int i = 0; i < g.nEdgeIDs(); i++) {
			if (!g.edgeEnabled(i)) {
				assert(!edgeInShortestPathGraph[i]);
			}
		}
		
		std::vector<int> dbg_delta;
		std::vector<Weight> dbg_dist;
		dbg_dist.resize(g.nodes(), INF);
		dbg_delta.resize(g.nodes());
		dbg_dist[getSource()] = 0;
		struct DistCmp {
			std::vector<Weight> & _dist;
			bool operator()(int a, int b) const {
				return _dist[a] < _dist[b];
			}
			DistCmp(std::vector<Weight> & d) :
					_dist(d) {
			}
			;
		};
		Heap<DistCmp> q(dbg_dist);
		
		q.insert(getSource());
		
		while (q.size()) {
			int u = q.removeMin();
			if (dbg_dist[u] == INF)
				break;
			dbg_delta[u] = 0;
			
			for (int i = 0; i < g.nIncoming(u); i++) {
				if (!g.edgeEnabled(g.incoming(u, i).id))
					continue;
				
				int edgeID = g.incoming(u, i).id;
				int v = g.getEdge(edgeID).from;
				Weight alt = dbg_dist[v] + weights[edgeID];
				assert(alt >= dbg_dist[u]);
				/*	if (alt==dbg_dist[u]){
				 dbg_delta[u]++;
				 }*/
			}
			
			for (int i = 0; i < g.nIncident(u); i++) {
				if (!g.edgeEnabled(g.incident(u, i).id))
					continue;
				
				int edgeID = g.incident(u, i).id;
				int v = g.getEdge(edgeID).to;
				Weight alt = dbg_dist[u] + weights[edgeID];
				if (alt < dbg_dist[v]) {
					
					dbg_dist[v] = alt;
					
					if (!q.inHeap(v))
						q.insert(v);
					else
						q.decrease(v);
				}/*else if (alt==dbg_dist[v]){
				 dbg_delta[v]++;
				 }*/
			}
		}
		
		for (int u = 0; u < g.nodes(); u++) {
			Weight & d = dist[u];
			
			assert(dbg_dist[u] == dist[u]);
			for (int i = 0; i < g.nIncoming(u); i++) {
				if (!g.edgeEnabled(g.incoming(u, i).id))
					continue;
				
				int edgeID = g.incoming(u, i).id;
				int v = g.getEdge(edgeID).from;
				
				Weight alt = dbg_dist[v] + weights[edgeID];
				assert(alt >= dbg_dist[u]);
				if (alt == dbg_dist[u]) {
					dbg_delta[u]++; //is this right?
					assert(edgeInShortestPathGraph[edgeID]);
				} else {
					assert(!edgeInShortestPathGraph[edgeID]);
				}
				
			}
			
		}
		for (int u = 0; u < g.nodes(); u++) {
			Weight& d = dist[u];
			Weight& d_expect = dbg_dist[u];
			assert(d == dbg_dist[u]);
			
		}
		for (int u = 0; u < g.nodes(); u++) {
			int d = delta[u];
			int d_expect = dbg_delta[u];
			assert(d == dbg_delta[u]);
			
		}
		dbg_delta_lite();
#endif
	}
	
	void GRRInc(int edgeID) {
		static int iter = 0;
		++iter;
		dbg_delta_lite();
		assert(g.edgeEnabled(edgeID));
		if (edgeInShortestPathGraph[edgeID])
			return;
		int ru = g.getEdge(edgeID).from;
		int rv = g.getEdge(edgeID).to;
		
		Weight & rdv = dist[rv];
		Weight & rdu = dist[ru];
		
		Weight& weight = weights[edgeID];
		if (dist[rv] < dist[ru] + weight)
			return;
		else if (dist[rv] == dist[ru] + weight) {
			assert(!edgeInShortestPathGraph[edgeID]);
			edgeInShortestPathGraph[edgeID] = true;
			delta[rv]++; //we have found an alternative shortest path to v
			return;
		}
		edgeInShortestPathGraph[edgeID] = true;
		delta[rv]++;
		dist[rv] = dist[ru] + weight;
		q.clear();
		q.insert(rv);
		
		while (q.size()) {
			int u = q.removeMin();
			
			if (!node_changed[u]) {
				node_changed[u] = true;
				changed.push_back(u);
			}
			delta[u] = 0;
			//for(auto & e:g.inverted_adjacency[u]){
			for (int i = 0; i < g.nIncoming(u); i++) {
				auto & e = g.incoming(u, i);
				int adjID = e.id;
				if (g.edgeEnabled(adjID)) {
					
					assert(g.getEdge(adjID).to == u);
					int v = g.getEdge(adjID).from;
					Weight & w = weights[adjID]; //assume a weight of one for now
					Weight & du = dist[u];
					Weight & dv = dist[v];
					if (dist[u] == dist[v] + w) {
						edgeInShortestPathGraph[adjID] = true;
						delta[u]++;
					} else if (dist[u] < (dist[v] + w)) {
						//This doesn't hold for us, because we are allowing multiple edges to be added at once.
						//assert(dist[u]<(dist[v]+w));
						
						edgeInShortestPathGraph[adjID] = false;
					} else {
						//don't do anything. This will get corrected in a future call to GRRInc.
						//assert(false);
						
					}
				} else {
					edgeInShortestPathGraph[adjID] = false;	//need to add this, because we may have disabled multiple edges at once.
				}
			}
			
			for (int i = 0; i < g.nIncident(u); i++) {
				auto & e = g.incident(u, i);
				int adjID = e.id;
				if (g.edgeEnabled(adjID)) {
					assert(g.getEdge(adjID).from == u);
					int s = g.getEdge(adjID).to;
					Weight & w = weights[adjID];							//assume a weight of one for now
					Weight & du = dist[u];
					Weight & ds = dist[s];
					if (dist[s] > dist[u] + w) {
						dist[s] = dist[u] + w;
						q.update(s);
					} else if (dist[s] == dist[u] + w && !edgeInShortestPathGraph[adjID]) {
						edgeInShortestPathGraph[adjID] = true;
						delta[s]++;
					}
				}
			}
		}
		dbg_delta_lite();
	}
	void dbg_delta_lite() {
#ifdef DEBUG_RAMAL
		for (int u = 0; u < g.nodes(); u++) {
			int del = delta[u];
			Weight d = dist[u];
			int num_in = 0;
			for (int i = 0; i < g.nIncoming(u); i++) {
				auto & e = g.incoming(u, i);
				int adjID = e.id;
				int from = g.getEdge(adjID).from;
				
				Weight dfrom = dist[from];
				if (edgeInShortestPathGraph[adjID])
					num_in++;
			}
			assert(del == num_in);
		}
#endif
		
	}
	
	void GRRDec(int edgeID) {
		dbg_delta_lite();
		assert(!g.edgeEnabled(edgeID));
		//First, check if this edge is actually in the shortest path graph
		if (!edgeInShortestPathGraph[edgeID])
			return;
		edgeInShortestPathGraph[edgeID] = false;						//remove this edge from the shortest path graph
				
		int ru = g.getEdge(edgeID).from;
		int rv = g.getEdge(edgeID).to;
		assert(delta[rv] > 0);
		delta[rv]--;
		if (delta[rv] > 0)
			return; //the shortest path hasn't changed in length, because there was an alternate route of the same length to this node.
			
		q.clear();
		changeset.clear();
		changeset.push_back(rv);
		
		//find all effected nodes whose shortest path lengths may now be increased (or that may have become unreachable)
		for (int i = 0; i < changeset.size(); i++) {
			int u = changeset[i];
			dist[u] = INF;
			for (int i = 0; i < g.nIncident(u); i++) {
				auto & e = g.incident(u, i);
				int adjID = e.id;
				if (g.edgeEnabled(adjID)) {
					if (edgeInShortestPathGraph[adjID]) {
						edgeInShortestPathGraph[adjID] = false;
						assert(g.getEdge(adjID).from == u);
						int s = g.getEdge(adjID).to;
						assert(delta[s] > 0);
						delta[s]--;
						if (delta[s] == 0) {
							changeset.push_back(s);
						}
					}
				}
			}
		}
		
		for (int i = 0; i < changeset.size(); i++) {
			int u = changeset[i];
			assert(dist[u] == INF);
			for (int i = 0; i < g.nIncoming(u); i++) {
				auto & e = g.incoming(u, i);
				int adjID = e.id;
				
				if (g.edgeEnabled(adjID)) {
					assert(g.getEdge(adjID).to == u);
					int v = g.getEdge(adjID).from;
					Weight & w = weights[adjID]; //assume a weight of one for now
					Weight alt = dist[v] + w;
					assert(!edgeInShortestPathGraph[adjID]);
					if (dist[u] > alt) {
						dist[u] = alt;
					}
				}
				
			}
			if (dist[u] != INF) {
				//q.insert(u);
				//dbg_Q_add(q,u);
				q.insert(u);
				
				if (!reportDistance && reportPolarity >= 0) {
					if (!node_changed[u]) {
						node_changed[u] = true;
						changed.push_back(u);
					}
				}
			} else if (reportPolarity <= 0) {
				//have to mark this change even if we are reporting distanec, as u has not been added to the queue.
				if (!node_changed[u]) {
					node_changed[u] = true;
					changed.push_back(u);
				}
			}
		}
		
		while (q.size() > 0) {
			int u = q.removeMin();
			if (reportDistance) {
				if (dist[u] != INF) {
					if (reportPolarity >= 0) {
						if (!node_changed[u]) {
							node_changed[u] = true;
							changed.push_back(u);
						}
					}
				} else if (reportPolarity <= 0) {
					if (!node_changed[u]) {
						node_changed[u] = true;
						changed.push_back(u);
					}
				}
			}
			for (int i = 0; i < g.nIncident(u); i++) {
				auto & e = g.incident(u, i);
				int adjID = e.id;
				if (g.edgeEnabled(adjID)) {
					assert(g.getEdge(adjID).from == u);
					int s = g.getEdge(adjID).to;
					Weight w = weights[adjID];				//assume a weight of one for now
					Weight alt = dist[u] + w;
					if (dist[s] > alt) {
						if (reportPolarity >= 0 && dist[s] >= 0) {
							//This check is needed (in addition to the above), because even if we are NOT reporting distances, it is possible for a node that was previously not reachable
							//to become reachable here. This is ONLY possible because we are batching multiple edge incs/decs at once (otherwise it would be impossible for removing an edge to decrease the distance to a node).
							if (!node_changed[s]) {
								node_changed[s] = true;
								changed.push_back(s);
							}
						}
						
						dist[s] = alt;
						q.update(s);
					} else if (dist[s] == alt && !edgeInShortestPathGraph[adjID]) {
						edgeInShortestPathGraph[adjID] = true;
						delta[s]++;							//added by sam... not sure if this is correct or not.
					}
				}
			}
			
			for (int i = 0; i < g.nIncoming(u); i++) {
				auto & e = g.incoming(u, i);
				int adjID = e.id;
				if (g.edgeEnabled(adjID)) {
					
					assert(g.getEdge(adjID).to == u);
					int v = g.getEdge(adjID).from;
					Weight & dv = dist[v];
					Weight & du = dist[u];
					bool edgeIn = edgeInShortestPathGraph[adjID];
					Weight & w = weights[adjID];							//assume a weight of one for now
					if (dist[u] == dist[v] + w && !edgeInShortestPathGraph[adjID]) {
						assert(!edgeInShortestPathGraph[adjID]);
						edgeInShortestPathGraph[adjID] = true;
						delta[u]++;
					} else if (dist[u] < dist[v] + w && edgeInShortestPathGraph[adjID]) {
						edgeInShortestPathGraph[adjID] = false;
						delta[u]--;
						assert(!edgeInShortestPathGraph[adjID]);
					} else if (dist[u] > dist[v] + w) {
						//assert(false);
					}
				}
			}
		}
		dbg_delta_lite();
	}
	
	long num_updates = 0;
	int numUpdates() const {
		return num_updates;
	}
	
	void update() {

		if (g.outfile) {
			fprintf(g.outfile, "r %d\n", getSource());
		}

		static int iteration = 0;
		int local_it = ++iteration;
		if (local_it == 7668) {
			int a = 1;
		}
		if (last_modification > 0 && g.modifications == last_modification)
			return;
		if (last_modification <= 0 || g.changed()) {
			INF = 1;							//g.nodes()+1;
			has_zero_weights=false;
			for (Weight & w : weights) {
				if (w <= 0) {
					//Note: in the future, we could implement the DFMN algorithm (Maintaining Shortest Paths in Digraphs with Arbitrary Arc Weights: An Experimental Study), which does support negative length weights, but is slower than RR.
					//throw std::invalid_argument("Ramalingham-Reps doesn't support zero-weight edges (select a different distance algorithm, such as dijkstra)");
					//for the moment: the _first_ time <0 weights are detected, simply fallback on dijkstra's, permanently.

					has_zero_weights=true;
				}
				INF += w;
			}
			dist.resize(g.nodes(), INF);
			dist[getSource()] = 0;
			delta.resize(g.nodes());
			node_changed.resize(g.nodes(), true);
			
			for (int i = 0; i < g.nodes(); i++) {
				if ((dist[i] >= INF && reportPolarity <= 0) || (dist[i] < INF && reportPolarity >= 0)) {
					node_changed[i] = true;
					changed.push_back(i);							//On the first round, report status of all nodes.
				}
			}
		}
		edgeInShortestPathGraph.resize(g.nEdgeIDs());
		if (last_history_clear != g.historyclears) {
			history_qhead = 0;
			last_history_clear = g.historyclears;
			
		}
		
		if(has_zero_weights){
			if(!ever_warned_about_zero_weights){
				ever_warned_about_zero_weights=true;
				fprintf(stderr,"Warning: Ramalingham-Reps doesn't support zero-weight edges; falling back on Dijkstra's (which is much slower)\n");
			}
			dijkstras.update();
		}else{
			for (int i = history_qhead; i < g.historySize(); i++) {
				int edgeid = g.getChange(i).id;
				if (g.getChange(i).addition && g.edgeEnabled(edgeid)) {
					GRRInc(edgeid);
				} else if (!g.getChange(i).addition && !g.edgeEnabled(edgeid)) {
					GRRDec(edgeid);
				}
			}
		}
			//for(int i = 0;i<g.nodes();i++){
			//	int u=i;
			for (int u : changed) {
				//int u = changed[i];
				node_changed[u] = false;
				//CANNOT clear the change flag here, because if we backtrack and then immediately re-propagate an edge before calling theoryPropagate, the change to this node may be missed.
				if (reportPolarity <= 0 && dist[u] >= INF) {
					status.setReachable(u, false);
					status.setMininumDistance(u, dist[u] < INF, dist[u]);
				} else if (reportPolarity >= 0 && dist[u] < INF) {
					status.setReachable(u, true);
					status.setMininumDistance(u, dist[u] < INF, dist[u]);
				}
			}
			changed.clear();


		assert(dbg_uptodate());
		num_updates++;
		last_modification = g.modifications;
		last_deletion = g.deletions;
		last_addition = g.additions;
		g.updateAlgorithmHistory(this,alg_id,history_qhead);
		history_qhead = g.historySize();
		last_history_clear = g.historyclears;
		
		;
	}
	void updateHistory(){
		update();
	}
	bool dbg_path(int to) {
#ifdef DEBUG_RAMAL
		/*	assert(connected(to));
		 if(to == source){
		 return true;
		 }
		 int p = previous(to);

		 if(p<0){
		 return false;
		 }
		 if(p==to){
		 return false;
		 }

		 return dbg_path(p);*/

#endif
		return true;
	}
	bool dbg_uptodate() {
#ifdef DEBUG_RAMAL
		/*if(last_modification<0)
		 return true;
		 dbg_delta();
		 Dijkstra<Weight> d(source,g,weights);

		 for(int i = 0;i<g.nodes();i++){
		 Weight dis = dist[i];
		 bool c = i<dist.size() && dist[i]<INF;
		 if(!c)
		 dis = this->unreachable();

		 Weight dbgdist = d.distance(i);

		 if(dis!=dbgdist){
		 assert(false);
		 throw std::logic_error();
		 }
		 }*/
//#endif
#endif
		return true;
	}
	
	bool connected_unsafe(int t) {
		dbg_uptodate();
		return t < dist.size() && dist[t] < INF;
	}
	bool connected_unchecked(int t) {
		assert(last_modification == g.modifications);
		return connected_unsafe(t);
	}
	bool connected(int t) {
		if (last_modification != g.modifications)
			update();

		assert(dbg_uptodate());
		
		return dist[t] < INF;
	}
	Weight & distance(int t) {
		if (last_modification != g.modifications)
			update();
		if (connected_unsafe(t))
			return dist[t];
		else
			return this->unreachable();
	}
	Weight &distance_unsafe(int t) {
		if (connected_unsafe(t))
			return dist[t];
		else
			return this->unreachable();
	}
	int incomingEdge(int t) {
		/*
		 assert(false);//not yet implemented...
		 assert(t>=0 && t<prev.size());
		 assert(prev[t]>=-1 );
		 return prev[t];*/
		//not supported
		 throw std::runtime_error("not implemented");
	}
	int previous(int t) {
		/*if(prev[t]<0)
		 return -1;

		 assert(g.all_edges[incomingEdge(t)].to==t);
		 return g.all_edges[incomingEdge(t)].from;*/
		 throw std::runtime_error("not implemented");
	}
};

template<typename Weight, class Status>
class UnweightedRamalReps: public Distance<int>, public DynamicGraphAlgorithm {
public:
	DynamicGraph<Weight> & g;
	Status & status;
	int reportPolarity;
	bool reportDistance;
	int last_modification;
	int last_addition;
	int last_deletion;
	int history_qhead;

	int last_history_clear;

	int source;
	int INF;

	int maxDistance = -1;

	std::vector<int> old_dist;
	std::vector<int> changed;
	std::vector<bool> node_changed;
	std::vector<int> dist;
	std::vector<int> prev;
	struct DistCmp {
		std::vector<int> & _dist;
		bool operator()(int a, int b) const {
			return _dist[a] < _dist[b];
		}
		DistCmp(std::vector<int> & d) :
				_dist(d) {
		}
		;
	};
	//Heap<DistCmp> q;
	std::vector<bool> in_queue;
	std::vector<bool> in_queue2;
	std::vector<int> q;
	std::vector<int> q2;
	std::vector<int> edgeInShortestPathGraph;
	std::vector<int> delta;
	std::vector<int> changeset;
	int alg_id;
public:
	
	long stats_full_updates=0;
	long stats_fast_updates=0;
	long stats_fast_failed_updates=0;
	long stats_skip_deletes=0;
	long stats_skipped_updates=0;
	long stats_num_skipable_deletions=0;
	double mod_percentage=0;

	double stats_full_update_time=0;
	double stats_fast_update_time=0;
	UnweightedRamalReps(int s, DynamicGraph<Weight> & graph, Status & status, int reportPolarity = 0,
			bool reportDistance = true) :
			g(graph), status(status), reportPolarity(reportPolarity), reportDistance(reportDistance), last_modification(
					-1), last_addition(-1), last_deletion(-1), history_qhead(0), last_history_clear(0), source(s), INF(
					-1) {
		maxDistance = -1;
		mod_percentage = 0.2;
		alg_id=g.addDynamicAlgorithm(this);
	}
	//Dijkstra(const Dijkstra& d):g(d.g), last_modification(-1),last_addition(-1),last_deletion(-1),history_qhead(0),last_history_clear(0),source(d.source),INF(0),q(DistCmp(dist)),stats_full_updates(0),stats_fast_updates(0),stats_skip_deletes(0),stats_skipped_updates(0),stats_full_update_time(0),stats_fast_update_time(0){marked=false;};
	void setMaxDistance(int &_maxDistance) {
		if (_maxDistance != maxDistance) {
			last_modification = -1;		//force the next update to recompute from scratch
			if (_maxDistance < 0) {
				maxDistance = INF;
			} else
				maxDistance = _maxDistance;
		}
	}
	
	void setSource(int s) {
		source = s;
		last_modification = -1;
		last_addition = -1;
		last_deletion = -1;
	}
	int getSource() {
		return source;
	}
	
	std::vector<int> & getChanged() {
		return changed;
	}
	void clearChanged() {
		changed.clear();
	}
	
	void drawFull() {
		
	}
	
	void dbg_delta() {
#ifdef DEBUG_RAMAL
		//g.drawFull();
		dbg_delta_lite();
		assert(delta.size() == g.nodes());
		
		for (int i = 0; i < g.nEdgeIDs(); i++) {
			if (!g.edgeEnabled(i)) {
				assert(!edgeInShortestPathGraph[i]);
				if (edgeInShortestPathGraph[i]) {
					throw std::runtime_error("");
				}
			}
		}
		
		std::vector<int> dbg_delta;
		std::vector<int> dbg_dist;
		dbg_dist.resize(g.nodes(), INF);
		dbg_delta.resize(g.nodes());
		dbg_dist[getSource()] = 0;
		
		struct DistCmp {
			std::vector<int> & _dist;
			bool operator()(int a, int b) const {
				return _dist[a] < _dist[b];
			}
			DistCmp(std::vector<int> & d) :
					_dist(d) {
			}
			;
		};
		Heap<DistCmp> q(dbg_dist);
		
		q.insert(getSource());
		
		while (q.size()) {
			int u = q.removeMin();
			if (dbg_dist[u] == INF)
				break;
			dbg_delta[u] = 0;
			
			for (int i = 0; i < g.nIncoming(u); i++) {
				if (!g.edgeEnabled(g.incoming(u, i).id))
					continue;
				
				int edgeID = g.incoming(u, i).id;
				int v = g.getEdge(edgeID).from;
				int alt = dbg_dist[v] + 1;
				if (maxDistance >= 0 && alt > maxDistance)
					alt = INF;
				assert(alt >= dbg_dist[u]);
				/*	if (alt==dbg_dist[u]){
				 dbg_delta[u]++;
				 }*/
			}
			
			for (int i = 0; i < g.nIncident(u); i++) {
				if (!g.edgeEnabled(g.incident(u, i).id))
					continue;
				
				int edgeID = g.incident(u, i).id;
				int v = g.getEdge(edgeID).to;
				int alt = dbg_dist[u] + 1;
				if (maxDistance >= 0 && alt > maxDistance)
					alt = INF;
				if (alt < dbg_dist[v]) {
					
					dbg_dist[v] = alt;
					
					if (!q.inHeap(v))
						q.insert(v);
					else
						q.decrease(v);
				}
			}
		}
		
		for (int u = 0; u < g.nodes(); u++) {
			int d = dist[u];
			
			int db = dbg_dist[u];
			assert(dbg_dist[u] == dist[u]);
			
			for (int i = 0; i < g.nIncoming(u); i++) {
				if (!g.edgeEnabled(g.incoming(u, i).id))
					continue;
				
				int edgeID = g.incoming(u, i).id;
				int v = g.getEdge(edgeID).from;
				int alt = dbg_dist[v] + 1;
				int du = dbg_dist[u];
				if (maxDistance >= 0 && alt > maxDistance)
					alt = INF;
				assert(alt >= dbg_dist[u]);
				if (alt == dbg_dist[u] && alt < INF) {
					dbg_delta[u]++;
					assert(edgeInShortestPathGraph[edgeID]);
				} else if (alt < INF) {
					assert(!edgeInShortestPathGraph[edgeID]);
				}
			}
		}
		
		for (int u = 0; u < g.nodes(); u++) {
			int d = dist[u];
			int d_expect = dbg_dist[u];
			assert(d == dbg_dist[u]);
		}
		for (int u = 0; u < g.nodes(); u++) {
			int du = dist[u];
			if (dist[u] < INF) {
				int d = delta[u];
				int d_expect = dbg_delta[u];
				assert(d == dbg_delta[u]);
			}
		}
		dbg_delta_lite();
#endif
	}
	
	void GRRInc(int edgeID) {
		static int iter = 0;
		++iter;
		dbg_delta_lite();
		assert(g.edgeEnabled(edgeID));
		if (edgeInShortestPathGraph[edgeID])
			return;
		int ru = g.getEdge(edgeID).from;
		int rv = g.getEdge(edgeID).to;
		
		int rdv = dist[rv];
		int rdu = dist[ru];
		
		int weight = 1;
		int altw = dist[ru] + weight;
		if (altw > maxDistance)
			altw = INF;
		if (dist[rv] < altw)
			return;
		else if (dist[rv] == altw && dist[rv] < INF) {
			assert(!edgeInShortestPathGraph[edgeID]);
			edgeInShortestPathGraph[edgeID] = true;
			delta[rv]++;		//we have found an alternative shortest path to v
			return;
		} else if (altw == INF) {
			return;		//don't do anything
		}
		assert(altw < INF);
		edgeInShortestPathGraph[edgeID] = true;
		delta[rv]++;
		dist[rv] = altw;
		
		q.clear();
		in_queue.clear();
		in_queue.resize(g.nodes());
		in_queue2.resize(g.nodes());
		q.push_back(rv);
		in_queue[rv] = true;
		
		for (int i = 0; i < q.size(); i++) {
			int u = q[i];
			dbg_Q_order(q);
			assert(dist[u] < INF);
			if (!node_changed[u]) {
				node_changed[u] = true;
				changed.push_back(u);
			}
			delta[u] = 0;
			//for(auto & e:g.inverted_adjacency[u]){
			for (int j = 0; j < g.nIncoming(u); j++) {
				auto & e = g.incoming(u, j);
				int adjID = e.id;
				if (g.edgeEnabled(adjID)) {
					
					assert(g.getEdge(adjID).to == u);
					int v = g.getEdge(adjID).from;
					int w = 1;		//assume a weight of one for now
					int du = dist[u];
					int dv = dist[v];
					int alt = dist[v] + w;
					if (alt > maxDistance)
						alt = INF;
					if (dist[u] == alt) {
						assert(alt < INF);
						edgeInShortestPathGraph[adjID] = true;
						delta[u]++;
					} else if (dist[u] < alt) {
						//This doesn't hold for us, because we are allowing multiple edges to be added at once.
						//assert(dist[u]<(dist[v]+w));
						
						edgeInShortestPathGraph[adjID] = false;
					} else {
						//don't do anything. This will get corrected in a future call to GRRInc.
						//assert(false);
						
					}
				} else {
					edgeInShortestPathGraph[adjID] = false;	//need to add this, because we may have disabled multiple edges at once.
				}
			}
			
			for (int j = 0; j < g.nIncident(u); j++) {
				auto & e = g.incident(u, j);
				int adjID = e.id;
				if (g.edgeEnabled(adjID)) {
					assert(g.getEdge(adjID).from == u);
					int s = g.getEdge(adjID).to;
					int w = 1;							//assume a weight of one for now
					int du = dist[u];
					int ds = dist[s];
					int alt = dist[u] + w;
					if (alt > maxDistance)
						alt = INF;
					if (dist[s] > alt) {
						dist[s] = alt;
						if (!in_queue[s]) {
							//q.update(s);
							dbg_Q_add(q, s);
							q.push_back(s);
							in_queue[s] = true;
						}
						dbg_not_seen_q(q, s, i);
					} else if (dist[s] == alt && !edgeInShortestPathGraph[adjID] && alt < INF) {
						assert(alt < INF);
						edgeInShortestPathGraph[adjID] = true;
						delta[s]++;
					}
				}
			}
		}
		dbg_delta_lite();
	}
	void dbg_not_seen_q(std::vector<int> & q, int u, int from) {
#ifdef DEBUG_RAMAL
		bool found = false;
		for (int i = from; i < q.size(); i++) {
			if (q[i] == u) {
				found = true;
				break;
			}
		}
		assert(found);
#endif
	}
	void dbg_Q_add(std::vector<int> & q, int u) {
#ifdef DEBUG_RAMAL
		//assert(!in_queue[u]);
		for (int v : q) {
			assert(u != v);
			assert(dist[v] <= dist[u]);
			//assert(in_queue[v]);
		}
		
#endif
	}
	void dbg_Q_order(std::vector<int> & _q) {
#ifdef DEBUG_RAMAL
		
		for (int i = 1; i < _q.size(); i++) {
			int v = _q[i];
			int u = _q[i - 1];
			
			if (&_q == &q) {
				assert(in_queue[v]);
				assert(in_queue[u]);
				if (!(in_queue2[u] || in_queue2[v])) {
					
					assert(dist[u] <= dist[v]);
				}
			} else {
				assert(in_queue2[u]);
				assert(in_queue2[v]);
				assert(dist[u] <= dist[v]);
			}
		}
		
#endif
	}
	
	void dbg_delta_lite() {
#ifdef DEBUG_RAMAL
		for (int u = 0; u < g.nodes(); u++) {
			int del = delta[u];
			int d = dist[u];
			int num_in = 0;
			for (int i = 0; i < g.nIncoming(u); i++) {
				auto & e = g.incoming(u, i);
				int adjID = e.id;
				int from = g.getEdge(adjID).from;
				
				int dfrom = dist[from];
				if (edgeInShortestPathGraph[adjID]) {
					
					assert(dist[u] < INF);
					num_in++;
				}
			}
			assert(del == num_in);
		}
#endif
		
	}
	
	void GRRDec(int edgeID) {
		dbg_delta_lite();
		assert(!g.edgeEnabled(edgeID));
		//First, check if this edge is actually in the shortest path graph
		if (!edgeInShortestPathGraph[edgeID])
			return;
		edgeInShortestPathGraph[edgeID] = false;						//remove this edge from the shortest path graph
				
		int ru = g.getEdge(edgeID).from;
		int rv = g.getEdge(edgeID).to;
		
		assert(delta[rv] > 0);
		delta[rv]--;
		if (delta[rv] > 0)
			return; //the shortest path hasn't changed in length, because there was an alternate route of the same length to this node.
			
		in_queue.clear();
		in_queue.resize(g.nodes());
		in_queue2.clear();
		in_queue2.resize(g.nodes());
		q.clear();
		q2.clear();
		
		changeset.clear();
		changeset.push_back(rv);
		
		//find all effected nodes whose shortest path lengths may now be increased (or that may have become unreachable)
		for (int i = 0; i < changeset.size(); i++) {
			int u = changeset[i];
			
			dist[u] = INF;
			for (int i = 0; i < g.nIncident(u); i++) {
				auto & e = g.incident(u, i);
				int adjID = e.id;
				if (g.edgeEnabled(adjID)) {
					if (edgeInShortestPathGraph[adjID]) {
						edgeInShortestPathGraph[adjID] = false;
						assert(g.getEdge(adjID).from == u);
						int s = g.getEdge(adjID).to;
						
						assert(delta[s] > 0);
						delta[s]--;
						if (delta[s] == 0) {
							changeset.push_back(s);
						}
					}
				}
			}
		}
		
		for (int i = 0; i < changeset.size(); i++) {
			int u = changeset[i];
			
			assert(dist[u] == INF);
			//for(auto & e:g.inverted_adjacency[u]){
			for (int i = 0; i < g.nIncoming(u); i++) {
				auto & e = g.incoming(u, i);
				int adjID = e.id;
				
				if (g.edgeEnabled(adjID)) {
					assert(g.getEdge(adjID).to == u);
					int v = g.getEdge(adjID).from;
					int w = 1; //assume a weight of one for now
					int alt = dist[v] + w;
					if (alt > maxDistance)
						alt = INF;
					assert(!edgeInShortestPathGraph[adjID]);
					if (dist[u] > alt) {
						dist[u] = alt;
					}
				}
				
			}
			
			if (dist[u] < INF) {
				//q.insert(u);
				//dbg_Q_add(q,u);
				q.push_back(u);
				in_queue[u] = true;
				
				if (reportPolarity >= 0) {
					if (!node_changed[u]) {
						node_changed[u] = true;
						changed.push_back(u);
					}
				}
			} else if (reportPolarity <= 0) {
				//call this even if we are reporting distance, because u hasn't been placed in the queue!
				if (!node_changed[u]) {
					node_changed[u] = true;
					changed.push_back(u);
				}
			}
			
		}
		std::sort(q.begin(), q.end(), DistCmp(dist));
		int i = 0, j = 0;
		while (i < q.size() || j < q2.size()) {
			int u;
			if (i == q.size()) {
				assert(j < q2.size());
				u = q2[j++];
			} else if (j == q2.size()) {
				assert(i < q.size());
				u = q[i++];
				if (in_queue2[u]) {
					continue;
				}
			} else if (dist[q[i]] < dist[q2[j]]) {
				u = q[i++];
				if (in_queue2[u]) {
					continue;
				}
			} else {
				assert(dist[q2[j]] <= dist[q[i]]);
				u = q2[j++];
			}
			assert(dist[u] < INF);
			if (reportDistance) {
				if (reportPolarity >= 0) {
					if (!node_changed[u]) {
						node_changed[u] = true;
						changed.push_back(u);
					}
				}
			}
			dbg_Q_order(q);
			dbg_Q_order(q2);
			
			for (int i = 0; i < g.nIncident(u); i++) {
				auto & e = g.incident(u, i);
				int adjID = e.id;
				if (g.edgeEnabled(adjID)) {
					assert(g.getEdge(adjID).from == u);
					int s = g.getEdge(adjID).to;
					int w = 1;				//assume a weight of one for now
					int alt = dist[u] + w;
					if (alt > maxDistance)
						alt = INF;
					if (dist[s] > alt) {
						assert(alt < INF);
						if (reportPolarity >= 0 && dist[s] >= 0) {
							//This check is needed (in addition to the above), because even if we are NOT reporting distances, it is possible for a node that was previously not reachable
							//to become reachable here. This is ONLY possible because we are batching multiple edge incs/decs at once (otherwise it would be impossible for removing an edge to decrease the distance to a node).
							if (!node_changed[s]) {
								node_changed[s] = true;
								changed.push_back(s);
							}
						}
						dist[s] = alt;
						if (!in_queue2[s]) {
							dbg_Q_add(q2, s);
							q2.push_back(s);
							in_queue2[s] = true;
						}
						
						dbg_Q_order(q2);
						//dbg_not_seen_q(q2,s,j);
					} else if (dist[s] == alt && !edgeInShortestPathGraph[adjID] && dist[s] < INF) {
						assert(dist[s] < INF);
						edgeInShortestPathGraph[adjID] = true;
						delta[s]++;								//added by sam... not sure if this is correct or not.
					}
				}
			}
			
			for (int i = 0; i < g.nIncoming(u); i++) {
				auto & e = g.incoming(u, i);
				int adjID = e.id;
				if (g.edgeEnabled(adjID)) {
					
					assert(g.getEdge(adjID).to == u);
					int v = g.getEdge(adjID).from;
					int dv = dist[v];
					int du = dist[u];
					bool edgeIn = edgeInShortestPathGraph[adjID];
					int w = 1;								//assume a weight of one for now
					int alt = dist[v] + w;
					if (alt > maxDistance)
						alt = INF;
					if (dist[u] == alt && !edgeInShortestPathGraph[adjID]) {
						assert(!edgeInShortestPathGraph[adjID]);
						assert(alt < INF);
						edgeInShortestPathGraph[adjID] = true;
						delta[u]++;
					} else if (dist[u] < alt && edgeInShortestPathGraph[adjID]) {
						edgeInShortestPathGraph[adjID] = false;
						delta[u]--;
						assert(!edgeInShortestPathGraph[adjID]);
					} else if (dist[u] > alt) {
						//assert(false);
					}
				}
			}
			
		}
		dbg_delta_lite();
	}
	
	long num_updates = 0;
	int numUpdates() const {
		return num_updates;
	}
	
	void update() {

		if (g.outfile) {
			fprintf(g.outfile, "r %d %d %d %d %d\n", getSource(),last_modification, g.modifications,g.changed(), g.historySize() );
		}

		if (last_modification > 0 && g.modifications == last_modification){
			return;
		}
		if (last_modification <= 0 || g.changed()) {//Note for the future: there is probably room to improve this further.
			stats_full_updates++;
			INF = g.nodes() + 1;
			dist.resize(g.nodes(), INF);
			dist[getSource()] = 0;
			delta.resize(g.nodes());
			edgeInShortestPathGraph.resize(g.nEdgeIDs());
			node_changed.resize(g.nodes());
			changed.clear();
			if (maxDistance < 0)
				maxDistance = INF;
			for (int i = 0; i < g.nodes(); i++) {
				if ((dist[i] >= INF && reportPolarity <= 0) || (dist[i] < INF && reportPolarity >= 0)) {
					node_changed[i] = true;
					changed.push_back(i);							//On the first round, report status of all nodes.
				}
			}
		}
		
		if (last_history_clear != g.historyclears) {
			history_qhead = 0;
			last_history_clear = g.historyclears;
			for (int edgeid = 0; edgeid < g.edges(); edgeid++) {
				if (g.edgeEnabled(edgeid)) {
					GRRInc(edgeid);
				} else {
					GRRDec(edgeid);
				}
			}
		}
		
		for (int i = history_qhead; i < g.historySize(); i++) {
			int edgeid = g.getChange(i).id;
			if (g.getChange(i).addition && g.edgeEnabled(edgeid)) {
				GRRInc(edgeid);
			} else if (!g.getChange(i).addition && !g.edgeEnabled(edgeid)) {
				GRRDec(edgeid);
			}
		}
		
		//for(int i = 0;i<g.nodes();i++){
		for (int u : changed) {
			//int u=i;
			//int u = changed[i];
			node_changed[u] = false;
			
			if (reportPolarity <= 0 && dist[u] >= INF) {
				status.setReachable(u, false);
				status.setMininumDistance(u, dist[u] < INF, dist[u]);
			} else if (reportPolarity >= 0 && dist[u] < INF) {
				status.setReachable(u, true);
				status.setMininumDistance(u, dist[u] < INF, dist[u]);
			}
		}
		changed.clear();
		//}
		dbg_delta();
		num_updates++;
		last_modification = g.modifications;
		last_deletion = g.deletions;
		last_addition = g.additions;
		
		history_qhead = g.historySize();
		g.updateAlgorithmHistory(this,alg_id,history_qhead);
		last_history_clear = g.historyclears;
		assert(dbg_uptodate());

	}
	void updateHistory(){
		update();
	}

	bool dbg_path(int to) {
#ifdef DEBUG_RAMAL
		assert(connected(to));
		if(to == source) {
			return true;
		}
		int p = previous(to);

		if(p<0) {
			return false;
		}
		if(p==to) {
			return false;
		}

		return dbg_path(p);

#endif
		return true;
	}
	bool dbg_uptodate() {
//#ifdef DEBUG_GRAPH
#ifdef DEBUG_RAMAL
		if(last_modification<0)
		return true;
		dbg_delta();
		UnweightedDijkstra<Weight> d(source,g);

		for(int i = 0;i<g.nodes();i++) {
			int dis = dist[i];
			bool c = i<dist.size() && dist[i]<INF;
			if(!c)
			dis = this->unreachable();
			if(maxDistance>=0 && dis>maxDistance) {
				dis=-1;
			}
			int dbgdist = d.distance(i);
			if(maxDistance>=0 && dbgdist>maxDistance) {
				dbgdist=-1;
			}
			if(dis!=dbgdist) {
				assert(false);
				throw std::runtime_error("");
			}
			if(d.connected(i) && d.distance(i)<maxDistance){

				int dd =d.dist[i];
				int mdis = dist[i];
				if(! (mdis< INF)){
					assert(false);
					throw std::runtime_error("");
				}
				if(dd!=mdis) {
					assert(false);
					throw std::runtime_error("");
				}
			}else{
				if(dist[i]<maxDistance){
					assert(false);
					throw std::runtime_error("");
				}
			}

		}
//#endif
#endif
		return true;
	}
	
	bool connected_unsafe(int t) {
		dbg_uptodate();
		return t < dist.size() && dist[t] < INF;
	}
	bool connected_unchecked(int t) {
		assert(last_modification == g.modifications);
		return connected_unsafe(t);
	}
	bool connected(int t) {
		if (last_modification < 0 ||  last_modification != g.modifications)
			update();


		assert(dbg_uptodate());
		
		return dist[t] < INF;
	}
	int& distance(int t) {
		if (last_modification < 0 ||  last_modification != g.modifications)
			update();

		if (connected_unsafe(t))
			return dist[t];
		else
			return this->unreachable();
	}
	int& distance_unsafe(int t) {
		if (connected_unsafe(t))
			return dist[t];
		else
			return this->unreachable();
	}
	int incomingEdge(int t) {
		/*
		 assert(false);//not yet implemented...
		 assert(t>=0 && t<prev.size());
		 assert(prev[t]>=-1 );
		 return prev[t];*/

		throw std::runtime_error("not implemented");
	}
	int previous(int t) {
		/*		if(prev[t]<0)
		 return -1;

		 assert(g.all_edges[incomingEdge(t)].to==t);
		 return g.all_edges[incomingEdge(t)].from;*/
		 throw std::runtime_error("not implemented");
	}
};
template<typename Weight, class Status>
bool RamalReps<Weight,Status>::ever_warned_about_zero_weights = 0;
}
;
#endif
