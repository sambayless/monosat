
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


#ifndef PKTOPOLOGICALSORT_H_
#define PKTOPOLOGICALSORT_H_

//Algorithm PK from "A Dynamic Topological Sort Algorithm for Directed Acyclic Graphs", 2006


#include <vector>
#include "alg/Heap.h"
#include "DynamicGraph.h"
#include "core/Config.h"
#include "Reach.h"
#include <algorithm>
#include "TarjansSCC.h"

namespace dgl {

class PKToplogicalSort: public Cycle {
public:

	DynamicGraph & g;
	TarjansSCC<> scc;

	std::vector<bool> is_strict_scc;
	std::vector<int> strict_sccs;
	//int num_invalidated_sccs=0;
	int num_strict_sccs=0;
	bool directed;
	int last_modification;
	int last_addition;
	int last_deletion;
	int history_qhead;

	int last_history_clear;

	int INF;

	std::vector<bool> in_cycle;
	std::vector<int> cycle;
	std::vector<int> prev;

	const int reportPolarity;

	//If true, force the (analyzed) graph to remain a DAG.
	//This is only valid if it is known that any discovered cycles will be eliminated.
	//In practice, what will happen is that as soon as a cycle is created on edge addition, that edge addition will be cancelled (but the cycle will be stored in "cycle" for analysis)
	bool force_dag=false;

	std::vector<int> ignore;
	std::vector<int> ord;
	std::vector<bool> visited;
	std::vector<bool> tmp_mark;
	std::vector<int> L;
	std::vector<int> R;
	std::vector<int> l_xy_B;

	std::vector<int> l_xy_F;
	//bool might_not_have_cycle=false;
	bool has_cycle = false;
	bool cycleComputed=false;
	bool has_topo=false;
	int nextOrd =0;
public:

	void forceDAG(){
		//If true, force the (analyzed) graph to remain a DAG.
		//This is only valid if it is known that any discovered cycles will be eliminated.
		//In practice, what will happen is that as soon as a cycle is created on edge addition, that edge addition will be cancelled (but the cycle will be stored in "cycle" for analysis)
		forceDAG=true;
	}

	PKToplogicalSort(DynamicGraph & graph, bool _directed = true, int _reportPolarity = 0) :
			g(graph), directed(_directed), last_modification(-1), last_addition(-1), last_deletion(-1), history_qhead(
					0), last_history_clear(0), INF(0), reportPolarity(_reportPolarity) {
		marked = false;
		mod_percentage = 0.2;
		stats_full_updates = 0;
		stats_fast_updates = 0;
		stats_skip_deletes = 0;
		stats_skipped_updates = 0;
		stats_full_update_time = 0;
		stats_fast_update_time = 0;
		stats_num_skipable_deletions = 0;
		stats_fast_failed_updates = 0;

	}

	void setNodes(int n) {
		in_cycle.resize(n);
		visited.resize(n);
		ord.resize(n);//how should this be initialised?
		is_strict_scc.resize(n);
		tmp_mark.resize(n);
		INF = g.nodes() + 1;
	}
private:

	bool inSCC(int node){
		assert(has_cycle);
		if(in_cycle[node])
			return true;
		int sccID = scc.getComponentUnsafe(node);
		return is_strict_scc[sccID];
	}

	void invalidateSCC(int node){
		assert(has_cycle);
		if(in_cycle[node]){
			num_strict_sccs=0;
			for(int n:cycle){
				assert(in_cycle[n]);
				in_cycle[n]=false;
			}
			cycle.clear();
			cycleComputed=false;
		}else{

			int sccID =scc.getComponentUnsafe(node);
			if (is_strict_scc[sccID]){
				is_strict_scc[sccID]=false;
				num_strict_sccs--;
				if(num_strict_sccs==0){
					cycleComputed=false;
				}
			}
			assert(num_strict_sccs>=0);
		}
	}

	void updateSCCs(){
		if(has_cycle && !cycleComputed){
			for(int n:cycle){
				assert(in_cycle[n]);
				in_cycle[n]=false;
			}
			cycle.clear();
			scc.update();
			for(int id:strict_sccs){
				is_strict_scc[id]=false;
			}
			num_strict_sccs=0;
			strict_sccs.clear();
			for(int id: scc.getStrictSCCs()){
				is_strict_scc[id]=true;
				num_strict_sccs++;
				strict_sccs.push_back(id);
			}
			has_cycle=num_strict_sccs>0;
			if(has_cycle){
				cycleComputed=true;
			}
		}
		if(has_cycle){
			assert(num_strict_sccs>0);
			assert(cycleComputed);
		}
	}

	void removeEdge(int edgeID){
		if(!has_cycle){
			//don't need to do anything at all.
			return;
		}else{
			//int from = g.getEdge(edgeID).from;
			int to = g.getEdge(edgeID).to;
			if(!cycleComputed){
				updateSCCs();
			}
			if(inSCC(to)){

				//then we have just broken an existing cycle, and need to check if the graph is still cyclic.
				invalidateSCC(to);
				//what to do in this case? for now, first: check if tarjan's scc still has any sccs that are not broken.
				if(num_strict_sccs>0){
					//then we still have at least one scc, somewhere. mark the current cycle invalid in scc.
				}else{
					updateSCCs();
				}
			}else{
				//there is still at least one cycle in this graph.
				return;
			}
		}
	}

	//PK Algorithm:

	int lower_bound;
	int upper_bound;
	void addEdge(int edgeID){
		//once the PK algorithm has a cycle, it is in an invalid state.
		if(has_cycle && !forceDAG)
			return;
		int from = g.getEdge(edgeID).from;
		int to = g.getEdge(edgeID).to;
		 lower_bound = ord[to];
		 upper_bound = ord[from];
		if(lower_bound<upper_bound){
			dfs_forward(to);
			if(!has_cycle){
				dfs_backward(from);
				reorder();
			}else{
				l_xy_F.clear();
				l_xy_B.clear();
				if(force_dag){
					has_cycle=false;
				}
			}
		}

	}

	void dfs_forward(int n){
		assert(!has_cycle);
		visited[n]=true;
		l_xy_F.push_back(n);
		for (int i = 0;!has_cycle && i<g.nIncident(n);i++){
			int edgeID = g.incident(n,i);
			if(g.edgeEnabled(edgeID)){
				int w = g.getEdge(edgeID).to;
				if(ord[w]==upper_bound){
					//cycle detected.
					has_cycle=true;
					//assert(cycleComputed==false);//cycle will be computed next time it is relevant.
					cycleComputed=true;
					num_strict_sccs=1;
					cycle.clear();
					scc.getSCC(w,g,cycle);

					return;
				}
				if(!visited[w] && ord[w]<upper_bound){
					dfs_forward(w);
				}
			}
		}
	}


	void dfs_backward(int n){
		assert(!has_cycle);
		visited[n]=true;
		l_xy_B.push_back(n);
		for (int i = 0;i<g.nIncoming(n);i++){
			int edgeID = g.incoming(n,i);
			if(g.edgeEnabled(edgeID)){
				int w = g.getEdge(edgeID).to;

				if(!visited[w] && lower_bound<ord[w]){
					dfs_backward(w);
				}
			}
		}
	}

	void merge(std::vector<int> & left, std::vector<int> & right,std::vector<int> & store ){
		store.clear();
		int lpos=0;
		int rpos=0;
		while(lpos< left.size() && rpos< right.size()){
			if(left[lpos]<=right[rpos]){
				store.push_back(left[lpos++]);
			}else{
				store.push_back(right[rpos++]);
			}
		}
		while(lpos< left.size()){
			store.push_back(left[lpos++]);
		}
		while(rpos< right.size()){
			store.push_back(right[rpos++]);
		}
	}

	void reorder(){
		assert(!has_cycle);
		//sort l_xy arrays with respect to their current ordinal
		std::sort(l_xy_B.begin(),l_xy_B.end(),ord);
		std::sort(l_xy_F.begin(),l_xy_F.end(),ord);


		L.clear();
		R.clear();
		for(int i =0;i<l_xy_B.size();i++){
			int w = l_xy_B[i];
			l_xy_B[i]=ord[w];
			visited[w]=false;
			L.push_back(w);
		}
		for(int i =0;i<l_xy_F.size();i++){
			int w = l_xy_F[i];
			l_xy_F[i]=ord[w];
			visited[w]=false;
			L.push_back(w);
		}

		merge(l_xy_B, l_xy_F,R);
		for(int i = 0;i<L.size();i++){
			ord[L[i]]=R[i];
		}
		l_xy_F.clear();
		l_xy_B.clear();
	}
public:
	void update() {
		static int iteration = 0;
		int local_it = ++iteration;

		if (last_modification > 0 && g.modifications == last_modification) {
			stats_skipped_updates++;
			//reportPolarity(undirected_cycle,directed_cycle);
			return;
		}

		if(last_modification<=0 || g.changed()){
			cycleComputed=false;
			has_cycle=false;
			has_topo=false;
			if(!topologicalSort()){
				updateSCCs();
			}
		}

		for (int i = history_qhead; i < g.history.size(); i++) {
			int edgeid = g.history[i].id;
			if (g.selfLoop(edgeid))
				continue; //skip self loops
			if (g.history[i].addition && g.edgeEnabled(edgeid) ) {
				addEdge(edgeID);
			} else if (!g.history[i].addition && !g.edgeEnabled(edgeid)) {
				removeEdge(edgeID);
			}
		}

		last_modification = g.modifications;
		last_deletion = g.deletions;
		last_addition = g.additions;

		history_qhead = g.history.size();
		last_history_clear = g.historyclears;

	}
private:
	bool sortVisit(int node){
		if(tmp_mark[node]){
			has_cycle=true;
			//assert(cycleComputed==false);//cycle will be computed next time it is relevant.
			cycleComputed=true;
			num_strict_sccs=1;
			cycle.clear();
			scc.getSCC(node,g,cycle);

			return false;
		}
		if(!visited[node]){
			tmp_mark[node];
			for (int j = 0;j<g.nIncident(node);j++){
				int edgeID = g.incident(node,j).id;
				if(g.edgeEnabled(edgeID)){
					if(!sortVisit(g.incident(node,j).node)){
						assert(has_cycle);
						return false;
					}
				}
			}
			visited[node]=true;
			tmp_mark[node]=false;
			ord[node]=nextOrd++;
		}
		return true;
	}
	//topologically sorts the graph, or returns false if it has a directed cycle.
	bool topologicalSort(){
		if(has_topo)
			return true;
		has_topo=true;
		nextOrd=0;
		//from wikipedia's pseudocode
		L.clear();
		for(int n = 0;n<g.nodes();n++){
			if(!visited[n]){
				if(!sortVisit(n)){
					assert(has_cycle);
					return has_topo=false;
					return false;
				}
			}
		}

		return has_topo;
	}
public:
	bool hasDirectedCycle() {
		update();
		return has_cycle;
	}

	//get _any_ directed cycle from this graph (must be cyclic)
	void getDirectedCycle(std::vector<int> & store) {
		update();
		store.clear();

		if(cycle.size()){
			store= cycle;
			return;
		}

		if(!hasDirectedCycle()){
			assert(cycle.size()==0);
			return ;
		}


		for(int id:strict_sccs){
			int element = scc.getElement(id);
			store.clear();
			scc.getSCC(element,g,store);
			if(store.size()>1){
				return;
			}
			store.clear();
			/*if(cycle.size()>1){
				for(int n:cycle){
					in_cycle[n]=true;
				}
				return cycle;
			}*/
		}
	}

	bool hasUndirectedCycle(){
		assert(false);
		return false;//not implemented
	}

	std::vector<int> & getUndirectedCycle(){
		assert(false);
		return ignore;
	}


};
}
;




#endif /* PKTOPOLOGICALSORT_H_ */
