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
 NONinf()RINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
 OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 **************************************************************************************************/

#ifndef DIJKSTRA_H_
#define DIJKSTRA_H_

#include <vector>
#include "alg/Heap.h"
#include "DynamicGraph.h"
#include "Reach.h"
#include "Distance.h"
#include "core/Config.h"
#include <limits>
namespace dgl {

template<typename Weight = long, class Status = typename Distance<Weight>::NullStatus, bool undirected = false>
class Dijkstra: public Distance<Weight> {
	using Distance<Weight>::inf;
	using Distance<Weight>::unreachable;
public:
	DynamicGraph<Weight> & g;

	Status & status;
	int reportPolarity;

	int last_modification=-1;
	int last_addition=0;
	int last_deletion=0;
	int last_edge_inc =0;
	int last_edge_dec = 0;
	int history_qhead=0;

	int last_history_clear=0;

	int source;
	//Weight inf();
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

public:
	
	int stats_full_updates = 0;
	int stats_fast_updates = 0;
	int stats_fast_failed_updates = 0;
	int stats_skip_deletes = 0;
	int stats_skipped_updates = 0;
	int stats_num_skipable_deletions = 0;
	double mod_percentage;

	double stats_full_update_time = 0;
	double stats_fast_update_time = 0;
	Dijkstra(int s, DynamicGraph<Weight> & graph, Status & status, int reportPolarity = 0) :
			g(graph), status(status), reportPolarity(reportPolarity), source(s),  q(DistCmp(dist)) {
		
		mod_percentage = 0.2;
		
	}
	
	Dijkstra(int s, DynamicGraph<Weight>  & graph,  int reportPolarity = 0) :
			g(graph), status(Distance<Weight>::nullStatus), reportPolarity(reportPolarity),  source(s), q(DistCmp(dist)) {
		
		mod_percentage = 0.2;
		//inf()=std::numeric_limits<Weight>::max()/2;
		
	}
	//Dijkstra(const Dijkstra& d):g(d.g), last_modification(-1),last_addition(-1),last_deletion(-1),history_qhead(0),last_history_clear(0),source(d.source),inf()(0),q(DistCmp(dist)),stats_full_updates(0),stats_fast_updates(0),stats_skip_deletes(0),stats_skipped_updates(0),stats_full_update_time(0),stats_fast_update_time(0){marked=false;};

	void setSource(int s) {
		source = s;
		last_modification = -1;
		last_addition = -1;
		last_deletion = -1;
	}
	int getSource() {
		return source;
	}
	
	void drawFull() {
		
	}
	long num_updates = 0;
	int numUpdates() const {
		return num_updates;
	}
	void update() {
		static int iteration = 0;
		int local_it = ++iteration;
		
		if (last_modification > 0 && g.modifications == last_modification)
			return;
		
		if (last_addition == g.additions && last_edge_inc==g.edge_increases  && last_edge_dec==g.edge_decreases  && last_modification > 0) {
			//if none of the deletions were to edges that were the previous edge of some shortest path, then we don't need to do anything
			if (last_history_clear != g.historyclears) {
				history_qhead = 0;
				last_history_clear = g.historyclears;
			}
			bool need_recompute = false;
			//ok, now check if any of the added edges allow for a decrease in distance.
			for (int i = history_qhead; i < g.historySize(); i++) {
				assert(!g.getChange(i).addition);
				int edgeid = g.getChange(i).id;
				int u = g.getEdge(edgeid).from;
				int v = g.getEdge(edgeid).to;
				if (incomingEdge(u) == edgeid || incomingEdge(v) == edgeid) {
					history_qhead = i - 1;
					need_recompute = true;
					//this deletion matters, so we need to recompute.
					break;
				}
				/*if(previous(v)==u || (undirected && previous(u)==v )){
				 history_qhead = i-1;
				 need_recompute=true;
				 //this deletion matters, so we need to recompute.
				 break;
				 }*/
			}
			if (!need_recompute) {
				//none of these deletions touched any shortest paths, so we can ignore them.
				
				last_modification = g.modifications;
				last_deletion = g.deletions;
				last_addition = g.additions;
				last_edge_inc = g.edge_increases;
				last_edge_dec = g.edge_decreases;
				history_qhead = g.historySize();
				last_history_clear = g.historyclears;
				
				assert(dbg_uptodate());
				stats_skip_deletes++;
				return;
			}
		}
		
		stats_full_updates++;
		if (dist.size() != g.nodes()) {
			
/*			inf() = 1;
			for (int edgeID = 0;edgeID<g.edges();edgeID++) {
				if(g.hasEdge(edgeID))
					inf() += g.getWeight(edgeID);
			}*/
			dist.resize(g.nodes());
			prev.resize(g.nodes());
		}
		
		//old_dist.resize(g.nodes());
		q.clear();
		for (int i = 0; i < g.nodes(); i++) {
			//old_dist[i]=last_modification > 0 ? dist[i]:inf();//this won't work properly if we added nodes...
			dist[i] =  inf();
			prev[i] = -1;
		}
		
		dist[source] = 0;
		q.insert(source);
		while (q.size()) {
			int u = q.peakMin();
			if (dist[u] == inf())
				break;
			/*if(old_dist[u]>=inf()){
			 changed.push_back(u);
			 }*/
			q.removeMin();
			for (int i = 0; i < g.nIncident(u, undirected); i++) {
				if (!g.edgeEnabled(g.incident(u, i, undirected).id))
					continue;
				int edgeID = g.incident(u, i, undirected).id;
				int v = g.incident(u, i, undirected).node;
				Weight alt = dist[u] +  g.getWeight(edgeID);
				if (alt < dist[v]) {
					dist[v] = alt;
					prev[v] = edgeID;
					if (!q.inHeap(v))
						q.insert(v);
					else
						q.decrease(v);
				}
			}
		}
		
		assert(dbg_uptodate());
		for (int u = 0; u < g.nodes(); u++) {
			if (reportPolarity <= 0 && dist[u] >= inf()) {
				status.setReachable(u, false);
				status.setMininumDistance(u, dist[u] < inf(), dist[u]);
			} else if (reportPolarity >= 0 && dist[u] < inf()) {
				status.setReachable(u, true);
				status.setMininumDistance(u, dist[u] < inf(), dist[u]);
			}
		}
		num_updates++;
		last_modification = g.modifications;
		last_deletion = g.deletions;
		last_addition = g.additions;
		last_edge_inc = g.edge_increases;
		last_edge_dec = g.edge_decreases;
		history_qhead = g.historySize();
		last_history_clear = g.historyclears;
		
	}
	bool dbg_path(int to) {
#ifdef DEBUG_DIJKSTRA
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
#ifdef DEBUG_DIJKSTRA
		if(last_modification<=0)
		return true;
		/*	DynamicGraph gdbg;
		 for(int i = 0;i<g.nodes();i++){
		 gdbg.addNode();
		 }

		 for(int i = 0;i<g.nodes();i++){
		 for(int j = 0;j<g.nIncident(u);j++){
		 int u = g.incident(i,j);
		 gdbg.addEdge(i,u);
		 }
		 }*/

		/*	Dijkstra d(source,g);
		 d.update();
		 for(int i = 0;i<g.nodes();i++){
		 int distance = dist[i];
		 int dbgdist = d.dist[i];
		 assert(distance==dbgdist);
		 }*/
#endif
		return true;
	}
	
	bool connected_unsafe(int t) {
		return t < dist.size() && dist[t] < inf();
	}
	bool connected_unchecked(int t) {
		assert(last_modification == g.modifications);
		return connected_unsafe(t);
	}
	bool connected(int t) {
		if (last_modification != g.modifications)
			update();
		
		assert(dbg_uptodate());
		
		return dist[t] < inf();
	}
	Weight & distance(int t) {
		if (last_modification != g.modifications)
			update();
		if (connected_unsafe(t))
			return dist[t];
		return this->unreachable();
	}
	Weight &distance_unsafe(int t) {
		if (connected_unsafe(t))
			return dist[t];
		else
			return this->unreachable();
	}
	int incomingEdge(int t) {
		assert(t >= 0 && t < prev.size());
		assert(prev[t] >= -1);
		return prev[t];
	}
	int previous(int t) {
		if (prev[t] < 0)
			return -1;
		if (undirected && g.getEdge(incomingEdge(t)).from == t) {
			return g.getEdge(incomingEdge(t)).to;
		}
		assert(g.getEdge(incomingEdge(t)).to == t);
		return g.getEdge(incomingEdge(t)).from;
	}
};

template<typename Weight, class Status =typename Distance<int>::NullStatus, bool undirected = false>
class UnweightedDijkstra: public Distance<int> {
	using Distance<int>::inf;
	using Distance<int>::unreachable;
public:
	DynamicGraph<Weight>  & g;
	Status & status;
	int reportPolarity;

	int last_modification=-1;
	int last_addition=-1;
	int last_deletion=-1;
	int history_qhead=0;

	int last_history_clear=0;

	int source;

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
	Heap<DistCmp> q;

public:
	
	int stats_full_updates = 0;
	int stats_fast_updates = 0;
	int stats_fast_failed_updates = 0;
	int stats_skip_deletes = 0;
	int stats_skipped_updates = 0;
	int stats_num_skipable_deletions = 0;
	double mod_percentage;

	double stats_full_update_time = 0;
	double stats_fast_update_time = 0;
	UnweightedDijkstra(int s, DynamicGraph<Weight>  & graph, Status & status, int reportPolarity = 0) :
			g(graph), status(status), reportPolarity(reportPolarity),source(s),  q(DistCmp(dist)) {
		
		mod_percentage = 0.2;
		
	}
	
	UnweightedDijkstra(int s, DynamicGraph<Weight>  & graph, int reportPolarity = 0) :
			g(graph), status(Distance<int>::nullStatus), reportPolarity(reportPolarity), source(s),  q(DistCmp(dist)) {
		
		mod_percentage = 0.2;
		
	}
	//Dijkstra(const Dijkstra& d):g(d.g), last_modification(-1),last_addition(-1),last_deletion(-1),history_qhead(0),last_history_clear(0),source(d.source),inf()(0),q(DistCmp(dist)),stats_full_updates(0),stats_fast_updates(0),stats_skip_deletes(0),stats_skipped_updates(0),stats_full_update_time(0),stats_fast_update_time(0){marked=false;};

	void setSource(int s) {
		source = s;
		last_modification = -1;
		last_addition = -1;
		last_deletion = -1;
	}
	int getSource() {
		return source;
	}
	/*
	 void updateFast(){
	 stats_fast_updates++;
	 



	 for(int i = 0;i<g.nodes();i++)
	 changed.push_back(i);
	 assert(last_deletion==g.deletions);
	 last_modification=g.modifications;
	 last_addition=g.additions;

	 dist.resize(g.nodes());
	 prev.resize(g.nodes());
	 q.clear();
	 if(last_history_clear!=g.historyclears){
	 history_qhead=0;
	 last_history_clear=g.historyclears;
	 }
	 //ok, now check if any of the added edges allow for a decrease in distance.
	 for (int i = history_qhead;i<g.historySize();i++){
	 assert(g.getChange(i).addition); //NOTE: Currently, this is glitchy in some circumstances - specifically, ./modsat -rinc=1.05 -rnd-restart  -conflict-shortest-path  -no-conflict-min-cut   -rnd-init -rnd-seed=01231 -rnd-freq=0.01 /home/sam/data/gnf/unit_tests/unit_test_17_reduced.gnf can trigger this assertion!
	 int u=g.getChange(i).u;
	 int v=g.getChange(i).v;
	 int edgeID = g.getChange(i).id;
	 int alt = dist[u]+1 ;
	 if(alt< dist[v]){

	 if(dist[v]>=inf()){
	 //this was changed
	 changed.push_back(v);
	 }

	 dist[v]=alt;
	 prev[v]=edgeID;

	 if(!q.inHeap(v))
	 q.insert(v);
	 else
	 q.decrease(v);
	 }else if (undirected){
	 int u=g.getChange(i).v;
	 int v=g.getChange(i).u;
	 int alt = dist[u]+1 ;
	 if(alt< dist[v]){
	 if(dist[v]>=inf()){
	 //this was changed
	 changed.push_back(v);
	 }

	 dist[v]=alt;
	 prev[v]=edgeID;

	 if(!q.inHeap(v))
	 q.insert(v);
	 else
	 q.decrease(v);
	 }
	 }



	 *
	 *
	 //Is this altered code still correct? Well, not for dijkstras, but probably for connectivity
	 if(dist[v]>=inf()){
	 //this was changed
	 changed.push_back(v);

	 dist[v]=alt;
	 prev[v]=u;

	 if(!q.inHeap(v))
	 q.insert(v);
	 }
	 *


	 }
	 history_qhead=g.historySize();
	 auto & adjacency = undirected? g.adjacency_undirected:g.adjacency;
	 while(q.size()){
	 int u = q.removeMin();
	 if(dist[u]==inf())
	 break;

	 for(int i = 0;i<adjacency[u].size();i++){
	 if(!g.edgeEnabled( adjacency[u][i].id))
	 continue;
	 int v = adjacency[u][i].to;
	 int edgeID = adjacency[u][i].id;
	 int alt = dist[u]+ 1;
	 if(alt<dist[v]){
	 if(dist[v]>=inf()){
	 //this was changed
	 changed.push_back(v);
	 }

	 dist[v]=alt;
	 prev[v]=edgeID;
	 if(!q.inHeap(v))
	 q.insert(v);
	 else
	 q.decrease(v);
	 }

	 }
	 }

	 for(int u:changed){
	 if(reportPolarity<=0 && dist[u]>=inf()){
	 status.setReachable(u,false);
	 status.setMininumDistance(u,dist[u]<inf(),dist[u]);
	 }else if (reportPolarity>=0 && dist[u]<inf()){
	 status.setReachable(u,true);
	 status.setMininumDistance(u,dist[u]<inf(),dist[u]);
	 }
	 }
	 changed.clear();

	 
	 }*/
	/*	std::vector<int> & getChanged(){
	 return changed;
	 }
	 void clearChanged(){
	 changed.clear();
	 }*/

	void drawFull() {
		
	}
	long num_updates = 0;
	int numUpdates() const {
		return num_updates;
	}
	void update() {
		static int iteration = 0;
		int local_it = ++iteration;
		
		if (last_modification > 0 && g.modifications == last_modification)
			return;
		
		if (last_addition == g.additions && last_modification > 0) {
			//if none of the deletions were to edges that were the previous edge of some shortest path, then we don't need to do anything
			if (last_history_clear != g.historyclears) {
				history_qhead = 0;
				last_history_clear = g.historyclears;
			}
			bool need_recompute = false;
			//ok, now check if any of the added edges allow for a decrease in distance.
			for (int i = history_qhead; i <g.historySize();i++) {
				assert(!g.getChange(i).addition);
				int edgeid = g.getChange(i).id;
				int u = g.getEdge(edgeid).from;
				int v = g.getEdge(edgeid).to;
				if (incomingEdge(u) == edgeid || incomingEdge(v) == edgeid) {
					history_qhead = i - 1;
					need_recompute = true;
					//this deletion matters, so we need to recompute.
					break;
				}
				/*if(previous(v)==u || (undirected && previous(u)==v )){
				 history_qhead = i-1;
				 need_recompute=true;
				 //this deletion matters, so we need to recompute.
				 break;
				 }*/
			}
			if (!need_recompute) {
				//none of these deletions touched any shortest paths, so we can ignore them.
				
				last_modification = g.modifications;
				last_deletion = g.deletions;
				last_addition = g.additions;
				
				history_qhead = g.historySize();
				last_history_clear = g.historyclears;
				
				assert(dbg_uptodate());
				stats_skip_deletes++;
				return;
			}
		}
		
		/*if(last_deletion==g.deletions && last_modification>0  ){
		 //Don't need to do anything at all.
		 if(last_addition==g.additions){
		 last_modification = g.modifications;
		 stats_skipped_updates++;
		 assert(dbg_uptodate());
		 return;
		 }
		 //we can use a faster, simple dijkstra update method
		 updateFast();
		 assert(dbg_uptodate());
		 return;
		 }*/

		stats_full_updates++;

		dist.resize(g.nodes());
		prev.resize(g.nodes());
		//old_dist.resize(g.nodes());
		q.clear();
		for (int i = 0; i < g.nodes(); i++) {
			//old_dist[i]=last_modification > 0 ? dist[i]:inf();//this won't work properly if we added nodes...
			dist[i] = inf();
			prev[i] = -1;
		}
		
		dist[source] = 0;
		q.insert(source);
		while (q.size()) {
			int u = q.peakMin();
			if (dist[u] == inf())
				break;
			/*if(old_dist[u]>=inf()){
			 changed.push_back(u);
			 }*/
			q.removeMin();
			for (int i = 0; i < g.nIncident(u, undirected); i++) {
				if (!g.edgeEnabled(g.incident(u, i, undirected).id))
					continue;
				int edgeID = g.incident(u, i, undirected).id;
				int v = g.incident(u, i, undirected).node;
				int alt = dist[u] + 1;
				if (alt < dist[v]) {
					dist[v] = alt;
					prev[v] = edgeID;
					if (!q.inHeap(v))
						q.insert(v);
					else
						q.decrease(v);
				}
			}
		}
		
		/*	for(int u = 0;u<g.nodes();u++){
		 //while(q.size()){
		 //iterate through the unreached nodes and check which ones were previously reached

		 if(last_modification <=0  || (old_dist[u]<inf() && dist[u]>=inf())){
		 changed.push_back(u);
		 }
		 }*/
		//}
		assert(dbg_uptodate());
		for (int u = 0; u < g.nodes(); u++) {
			if (reportPolarity <= 0 && dist[u] >= inf()) {
				status.setReachable(u, false);
				status.setMininumDistance(u, dist[u] < inf(), dist[u]);
			} else if (reportPolarity >= 0 && dist[u] < inf()) {
				status.setReachable(u, true);
				status.setMininumDistance(u, dist[u] < inf(), dist[u]);
			}
		}
		num_updates++;
		last_modification = g.modifications;
		last_deletion = g.deletions;
		last_addition = g.additions;

		history_qhead = g.historySize();
		last_history_clear = g.historyclears;
		
	}
	bool dbg_path(int to) {
#ifdef DEBUG_DIJKSTRA
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
#ifdef DEBUG_DIJKSTRA
		if(last_modification<=0)
		return true;
		/*	DynamicGraph gdbg;
		 for(int i = 0;i<g.nodes();i++){
		 gdbg.addNode();
		 }

		 for(int i = 0;i<g.nodes();i++){
		 for(int j = 0;j<g.nIncident(u);j++){
		 int u = g.incident(i,j);
		 gdbg.addEdge(i,u);
		 }
		 }*/

		/*	Dijkstra d(source,g);
		 d.update();
		 for(int i = 0;i<g.nodes();i++){
		 int distance = dist[i];
		 int dbgdist = d.dist[i];
		 assert(distance==dbgdist);
		 }*/
#endif
		return true;
	}
	
	bool connected_unsafe(int t) {
		return t < dist.size() && dist[t] < inf();
	}
	bool connected_unchecked(int t) {
		assert(last_modification == g.modifications);
		return connected_unsafe(t);
	}
	bool connected(int t) {
		if (last_modification != g.modifications)
			update();
		
		assert(dbg_uptodate());
		
		return dist[t] < inf();
	}
	int & distance(int t) {
		if (last_modification != g.modifications)
			update();
		if (connected_unsafe(t))
			return dist[t];
		return this->unreachable();
	}
	int &distance_unsafe(int t) {
		if (connected_unsafe(t))
			return dist[t];
		else
			return this->unreachable(); //return inf();
	}
	int incomingEdge(int t) {
		assert(t >= 0 && t < prev.size());
		assert(prev[t] >= -1);
		return prev[t];
	}
	int previous(int t) {
		if (prev[t] < 0)
			return -1;
		if (undirected && g.getEdge(incomingEdge(t)).from == t) {
			return g.getEdge(incomingEdge(t)).to;
		}
		assert(g.getEdge(incomingEdge(t)).to == t);
		return g.getEdge(incomingEdge(t)).from;
	}
};
}
;
#endif /* DIJKSTRA_H_ */
