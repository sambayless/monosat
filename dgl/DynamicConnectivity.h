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

#ifndef DYNAMIC_CONNECTIVITY_H_
#define DYNAMIC_CONNECTIVITY_H_

#include <vector>
#include "alg/Heap.h"
#include "DynamicGraph.h"
#include "alg/EulerTree.h"
#include "alg/TreapCustom.h"
#include "core/Config.h"
#include "Reach.h"
#include "mtl/Sort.h"
#include "ThorupDynamicConnectivity.h"

namespace dgl {
template<typename Weight, class Status>
class DynamicConnectivity: public Reach, public AllPairs, public ConnectedComponents {
public:
	
	DynamicGraph<Weight> & g;
	Status & status;
	int last_modification;
	int last_addition;
	int last_deletion;
	int history_qhead;

	int last_history_clear;
	bool hasPrev;
	std::vector<int> sources;
	//the transitive closure of the graph (only for components connected to one of the sources)
	struct ClosureData {
		//	int incomingEdge;
		char reachable :1;
		char changed :1;
		ClosureData() :
				reachable(false), changed(false) {
		}
	};
	std::vector<std::vector<ClosureData>> transitive_closure;
	struct Change {
		int sourceNum;
		int node;
	};
	std::vector<Change> changes;
	std::vector<int> prev;
	std::vector<int> path;
	std::vector<int> u_component;
	std::vector<int> v_component;
	std::vector<int> component;
	int INF;

	ThorupDynamicConnectivity t;
	int default_source=0;
	int default_source_index=0;
	const int reportPolarity;

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

	DynamicConnectivity(DynamicGraph<Weight> & graph, Status & status, int _reportPolarity = 0) :
			g(graph), status(status), last_modification(-1), last_addition(-1), last_deletion(-1), history_qhead(0), last_history_clear(
					0), INF(0), reportPolarity(_reportPolarity) {
		default_source = -1;
		hasPrev = false;
		iteration = 0;
	}
	
	void addSource(int s) {
		int sourceNum = sources.size();
		sources.push_back(s);
		changes.push_back( { sourceNum, s });
		
		transitive_closure.push_back( { });
		transitive_closure[sourceNum].resize(g.nodes());
		transitive_closure[sourceNum][s].reachable = true;
		transitive_closure[sourceNum][s].changed = true;
		last_modification = -1;
		last_addition = -1;
		last_deletion = -1;
	}
	
	void setSource(int s) {
		if (!std::count(sources.begin(), sources.end(), s)) {
			addSource(s);
		}
		default_source_index = sources.size() - 1;
		default_source = s;
		
	}
	int getSource() {
		return default_source;
	}
	
	void setNodes(int n) {
		
		INF = g.nodes() + 1;
		
		while (t.nNodes() < g.nodes()) {
			t.addNode();
			prev.push_back(-1);
		}
		for (int i = 0; i < transitive_closure.size(); i++) {
			transitive_closure[i].resize(g.nodes());
		}
		
	}
	
	void updateEdge(int u, int v, int edgeid, bool add) {
		if (edgeid == 2) {
			int a = 1;
		}
		if (u == 16 || v == 16) {
			int a = 1;
		}
		assert(transitive_closure[0][sources[0]].reachable);
		assert(
				(g.all_edges[edgeid].to == u && g.all_edges[edgeid].from == v)
						|| (g.all_edges[edgeid].to == v && g.all_edges[edgeid].from == u));
		if (add) {
			bool already_connected = false;
			if (!t.connected(u, v)) {
				already_connected = false;
				u_component.clear();
				v_component.clear();
				for (int i = 0; i < sources.size(); i++) {
					
					//then it may be the case that new nodes are connected to one of the sources we are tracking
					if (transitive_closure[i][u].reachable && !transitive_closure[i][v].reachable) {
						//ok, explore the new connected component and enclose it
						
						if (!v_component.size())
							t.getConnectedComponentEdges(v, v_component, true);
						assert(!transitive_closure[i][v].reachable);
						
						transitive_closure[i][v].reachable = true;
						if (!transitive_closure[i][v].changed) {
							transitive_closure[i][v].changed = true;
							changes.push_back( { i, v });
						}
						
						for (int edgeid : v_component) {
							
							int u = g.all_edges[edgeid].from;
							int v = g.all_edges[edgeid].to;
							if (u == 16 || v == 16) {
								int a = 1;
							}
							assert(!(transitive_closure[i][u].reachable && transitive_closure[i][v].reachable));
							assert((transitive_closure[i][u].reachable || transitive_closure[i][v].reachable));
							int n = transitive_closure[i][u].reachable ? v : u;
							assert(!transitive_closure[i][n].reachable);
							transitive_closure[i][n].reachable = true;
							
							if (!transitive_closure[i][n].changed) {
								transitive_closure[i][n].changed = true;
								changes.push_back( { i, n });
							}
						}
					} else if (transitive_closure[i][v].reachable && !transitive_closure[i][u].reachable) {
						//ok, explore the new connected component and enclose it
						if (!u_component.size())
							t.getConnectedComponentEdges(u, u_component, true);
						assert(!transitive_closure[i][u].reachable);
						
						transitive_closure[i][u].reachable = true;
						if (!transitive_closure[i][u].changed) {
							transitive_closure[i][u].changed = true;
							changes.push_back( { i, u });
						}
						
						for (int edgeid : u_component) {
							int u = g.all_edges[edgeid].from;
							int v = g.all_edges[edgeid].to;
							if (u == 16 || v == 16) {
								int a = 1;
							}
							assert(!(transitive_closure[i][u].reachable && transitive_closure[i][v].reachable));
							assert((transitive_closure[i][u].reachable || transitive_closure[i][v].reachable));
							int n = transitive_closure[i][u].reachable ? v : u;
							assert(!transitive_closure[i][n].reachable);
							transitive_closure[i][n].reachable = true;
							
							if (!transitive_closure[i][n].changed) {
								transitive_closure[i][n].changed = true;
								changes.push_back( { i, n });
							}
						}
					}
				}
			} else {
				already_connected = true;
			}
			//avoiding an extra connected check here by using the unchecked variatn...
			t.setEdgeEnabledUnchecked(u, v, edgeid, already_connected);
		} else {
			if (t.setEdgeEnabled(u, v, edgeid, add)) {
				assert(!t.connected(u, v));
				u_component.clear();
				v_component.clear();
				for (int i = 0; i < sources.size(); i++) {
					int source = sources[i];
					//then it may be the case that new nodes are connected to one of the sources we are tracking
					if (transitive_closure[i][u].reachable) {
						assert(transitive_closure[i][v].reachable);
						//ok, explore the new connected component and enclose it
						if (t.connected(source, u)) {
							assert(!t.connected(source, v));
							if (!v_component.size())
								t.getConnectedComponent(v, v_component);
							for (int n : v_component) {
								if (n == 16) {
									int a = 1;
								}
								transitive_closure[i][n].reachable = false;
								
								if (!transitive_closure[i][n].changed) {
									transitive_closure[i][n].changed = true;
									changes.push_back( { i, n });
								}
							}
						} else {
							assert(t.connected(source, v));
							assert(!t.connected(source, u));
							if (!u_component.size())
								t.getConnectedComponent(u, u_component);
							for (int n : u_component) {
								if (n == 16) {
									int a = 1;
								}
								transitive_closure[i][n].reachable = false;
								
								if (!transitive_closure[i][n].changed) {
									transitive_closure[i][n].changed = true;
									changes.push_back( { i, n });
								}
							}
						}
					}
				}
			}
		}
		assert(transitive_closure[0][sources[0]].reachable);
	}
	
#ifndef NDEBUG
	DisjointSets dbg_sets;
#endif
	int iteration;
	long num_updates = 0;
	int numUpdates() const {
		return num_updates;
	}
	void update() {
		
		int local_it = ++iteration;
		
		if (last_modification > 0 && g.modifications == last_modification) {
			stats_skipped_updates++;
			return;
		}
		stats_full_updates++;
		
		setNodes(g.nodes());
		hasPrev = false;
		
#ifndef NDEBUG
		dbg_sets.Reset();
		
		dbg_sets.AddElements(g.nodes());
		
		for (int i = 0; i < g.all_edges.size(); i++) {
			if (g.edgeEnabled(i) && g.all_edges[i].id >= 0) {
				int u = g.all_edges[i].from;
				int v = g.all_edges[i].to;
				dbg_sets.UnionElements(u, v);
			}
		}
		
#endif
		dbg_transitive_closure();
		
		assert(transitive_closure[0][sources[0]].reachable);
		
		if (last_modification <= 0 || g.historyclears != last_history_clear) {
			last_history_clear = g.historyclears;
			history_qhead = 0;
			
			//initialize the transitive closure.
			for (int s = 0; s < sources.size(); s++) {
				int source = sources[s];
				for (int n = 0; n < g.nodes(); n++) {
					bool connected = t.connected(source, n);
					if (connected && reportPolarity >= 0) {
						transitive_closure[s][n].reachable = true;
						if (!transitive_closure[s][n].changed) {
							transitive_closure[s][n].changed = true;
							changes.push_back( { s, n });
						}
					} else if (!connected && reportPolarity <= 0) {
						transitive_closure[s][n].reachable = false;
						if (!transitive_closure[s][n].changed) {
							transitive_closure[s][n].changed = true;
							changes.push_back( { s, n });
						}
					}
				}
			}
			
			//start from scratch
			for (int i = 0; i < g.all_edges.size(); i++) {
				if (g.all_edges[i].id >= 0) {
					bool add = g.edgeEnabled(i);
					int u = g.all_edges[i].from;
					int v = g.all_edges[i].to;
					updateEdge(u, v, i, add);
				}
			}
			
		} else {
			//incremental/decremental update
			for (; history_qhead < g.historySize(); history_qhead++) {
				int edgeid = g.getChange(history_qhead).id;
				bool add = g.getChange(history_qhead).addition;
				int u = g.all_edges[edgeid].from;
				int v = g.all_edges[edgeid].to;
				updateEdge(u, v, edgeid, add);
			}
			
		}
#ifndef NDEBUG
		for (int i = 0; i < g.edges(); i++) {
			if (g.all_edges[i].id >= 0) {
				assert(t.edges[i].edgeID == g.all_edges[i].id);
				assert(t.edges[i].from == g.all_edges[i].from);
				assert(t.edges[i].to == g.all_edges[i].to);
				assert(t.edgeEnabled(i) == g.edgeEnabled(i));
			}
		}
		for (int i = 0; i < g.nodes(); i++) {
			for (int j = 0; j < g.nodes(); j++) {
				bool dbg_connected = dbg_sets.FindSet(i) == dbg_sets.FindSet(j);
				assert(dbg_connected == t.connected(i, j));
			}
		}
#endif
		for (Change c : changes) {
			int source = sources[c.sourceNum];
			assert(transitive_closure[c.sourceNum][c.node].changed);
			transitive_closure[c.sourceNum][c.node].changed = false;
			if (source == default_source) {
				bool reachable = transitive_closure[c.sourceNum][c.node].reachable;
				if (reachable && reportPolarity >= 0) {
					assert(t.connected(source, c.node));
					status.setReachable(c.node, true);
				} else if (!reachable && reportPolarity <= 0) {
					assert(!t.connected(source, c.node));
					status.setReachable(c.node, false);
				}
			}
		}
		changes.clear();
#ifndef NDEBUG
		{
			for (std::vector<ClosureData> & v : transitive_closure)
				for (ClosureData & c : v)
					assert(!c.changed);
		}
#endif
		dbg_transitive_closure();
		/*

		 if (default_source>=0){
		 prev.clear();
		 prev.resize(g.nodes(),-1);
		 //ok, traverse the nodes connected to this component
		 component.clear();
		 static int iter=0;
		 //this is NOT the right way to do this.
		 //need to only see check from t!
		 t.getConnectedComponent(default_source,component);

		 if(reportPolarity>=0){
		 for(int v:component){
		 status.setReachable(v,true);
		 }
		 }
		 if(reportPolarity<=0){
		 sort(component);
		 int nextpos = 0;
		 for(int v = 0;v<g.nodes();v++){
		 if(nextpos<component.size() && v==component[nextpos]){
		 nextpos++;
		 assert(dbg_sets.FindSet(v)==dbg_sets.FindSet(default_source));
		 //status.setReachable(v,true);

		 }else{
		 assert(dbg_sets.FindSet(v)!=dbg_sets.FindSet(default_source));
		 status.setReachable(v,false);

		 }
		 }
		 }
		 }
		 */

		assert(dbg_sets.NumSets() == t.numComponents());
		
		num_updates++;
		last_modification = g.modifications;
		last_deletion = g.deletions;
		last_addition = g.additions;
		
		history_qhead = g.historySize();
		last_history_clear = g.historyclears;
		
		//stats_full_update_time+=cpuTime()-startdupdatetime;;
	}
	
	void dbg_transitive_closure() {
#ifndef NDEBUG
		for (int i = 0; i < sources.size(); i++) {
			int source = sources[i];
			assert(transitive_closure[i][source].reachable);
			for (int n = 0; n < transitive_closure[i].size(); n++) {
				bool reachable = t.connected(source, n);
				assert(((bool )transitive_closure[i][n].reachable) == reachable);
				
			}
		}
		
#endif
	}
	void dbg_path_edges(int from, int to, std::vector<int> & path) {
#ifndef NDEBUG
		assert(path.size());
		int n = from;
		
		for (int i = 0; i < path.size(); i++) {
			int edgeid = path[i];
			
			int v = g.all_edges[edgeid].from;
			int u = g.all_edges[edgeid].to;
			assert(v == n || u == n);
			assert(g.hasEdgeUndirected(u, v));
			if (v == n) {
				n = u;
			} else {
				n = v;
			}
		}
		assert(n == to);
#endif
	}
	void dbg_path(int from, int to, std::vector<int> & path) {
#ifndef NDEBUG
		assert(path.size());
		assert(path[0] == from);
		assert(path.back() == to);
		for (int i = 1; i < path.size(); i++) {
			int v = path[i];
			int u = path[i - 1];
			assert(g.hasEdgeUndirected(u, v));
			
		}
#endif
	}
	
	bool connected(int from, int to) {
		update();
		return t.connected(from, to);
	}
	
	int numComponents() {
		update();
		return t.numComponents();
	}
	int getComponent(int node) {
		update();
		return t.findRoot(node);
	}
	
	bool connected_unsafe(int from, int to) {
		return t.connected(from, to);
	}
	bool connected_unchecked(int from, int to) {
		return t.connected(from, to);
	}
	
	int distance(int from, int to) {
		update();
		return t.connected(from, to) ? 0 : INF;
	}
	int distance_unsafe(int from, int to) {
		return t.connected(from, to) ? 0 : INF;
	}
	void getPath(int source, int to, std::vector<int> & path_store) {
		return t.getPath(source, to, path_store);
	}
	
	bool connected_unsafe(int t) {
		return connected_unsafe(getSource(), t);
	}
	bool connected_unchecked(int t) {
		return connected_unchecked(getSource(), t);
	}
	bool connected(int to) {
		update();
		assert(default_source_index >= 0);
		assert(to >= 0);
		dbg_transitive_closure();
		assert(transitive_closure[default_source_index].size() > to);
		assert((bool )(transitive_closure[default_source_index][to].reachable) == t.connected(getSource(), to));
		return transitive_closure[default_source_index][to].reachable;
		//return connected(getSource(),t);
	}
	int distance(int t) {
		return distance(getSource(), t);
	}
	int distance_unsafe(int t) {
		return distance_unsafe(getSource(), t);
	}
	/*		 int previous( int to){
	 update();
	 if(prev[to]<0){

	 t.getPath(getSource(),to,path);
	 for(int i:path){
	 printf("%d->",i);
	 }
	 printf("\n");

	 dbg_path(getSource(),to,path);

	 prev.clear();
	 prev.resize(g.nodes(),-1);
	 assert(path[0]==getSource());
	 for(int i = 1;i<path.size();i++){
	 prev[path[i]]=path[i-1];
	 }

	 }
	 assert(prev[to]!=to);
	 assert(prev[to]>=0);
	 return prev[to];
	 }*/

	int incomingEdge(int to) {
		update();
		assert(to >= 0);
		if (!hasPrev) {
			hasPrev = true;
			prev.clear();
			prev.resize(g.nodes(), -1);
		}
		if (prev[to] < 0) {
			assert(t.connected(getSource(), to));
			t.getPathEdges(getSource(), to, path);
			/*		for(int i:path){
			 printf("%d->",i);
			 }
			 printf("\n");
			 */
			dbg_path_edges(getSource(), to, path);
			
			int n = getSource();
			
			for (int edge : path) {
				int v = g.all_edges[edge].from;
				assert(v == n || g.all_edges[edge].to == n);
				if (v == n) {
					v = g.all_edges[edge].to;
				}
				assert(v != getSource());
				prev[v] = edge;
				n = v;
			}
			assert(n == to);
		}
		dbg_transitive_closure();
		return prev[to];
		//return prev[to];
	}
	int previous(int t) {
		int edgeID = incomingEdge(t);
		if (edgeID < 0)
			return -1;
		assert(transitive_closure[default_source_index][t].reachable);
		
		if (g.all_edges[edgeID].from == t) {
			return g.all_edges[edgeID].to;
		}
		assert(g.all_edges[edgeID].to == t);
		return g.all_edges[edgeID].from;
	}
	
	void getPath(int t, std::vector<int> & path_store) {
		getPath(getSource(), t, path_store);
	}
};
}
;
#endif
