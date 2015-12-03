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

#include "ConnectedComponentsDetector.h"
#include "GraphTheory.h"
#include "dgl/BFS.h"
#include "dgl/Reach.h"
#include <limits>
using namespace Monosat;
template<typename Weight>
ConnectedComponentsDetector<Weight>::ConnectedComponentsDetector(int _detectorID, GraphTheorySolver<Weight> * _outer,
		DynamicGraph<Weight>  &_g, DynamicGraph<Weight>  &_antig, double seed) :
		Detector(_detectorID), outer(_outer), g_under(_g), g_over(_antig), rnd_seed(seed), underapprox_component_detector(
				NULL), overapprox_component_detector(NULL), positiveReachStatus(NULL), negativeReachStatus(NULL) {
	
	positiveReachStatus = new ConnectedComponentsDetector<Weight>::ConnectedComponentsStatus(*this, true);
	negativeReachStatus = new ConnectedComponentsDetector<Weight>::ConnectedComponentsStatus(*this, false);
	//Note: these are _intentionalyl_ swapped
	overapprox_component_detector = new DisjointSetsConnectedComponents<Weight,
			ConnectedComponentsDetector<Weight>::ConnectedComponentsStatus>(_g, *(negativeReachStatus), 1);
	underapprox_component_detector = new DisjointSetsConnectedComponents<Weight,
			ConnectedComponentsDetector<Weight>::ConnectedComponentsStatus>(_antig, *(positiveReachStatus), 1);
	
	components_low_marker = outer->newReasonMarker(getID());
	components_high_marker = outer->newReasonMarker(getID());
	
	connected_marker = outer->newReasonMarker(getID());
	not_connected_marker = outer->newReasonMarker(getID());
	
}
template<typename Weight>
Lit ConnectedComponentsDetector<Weight>::getConnectLit(int u, int v) {
	if (reachLits.size() > u && reachLits[u].size() > v && reachLits[u][v] != lit_Undef) {
		return reachLits[u][v];
	}
	if (reachLits.size() < g_under.nodes()) {
		reachLits.growTo(g_under.nodes());
		for (int i = 0; i < reachLits.size(); i++) {
			reachLits[i].growTo(g_under.nodes(), lit_Undef);
		}
	}
	Var connect_var = outer->newVar(getID(), true);
	reachLits[u][v] = mkLit(connect_var);
	reachLits[v][u] = reachLits[u][v];
	reach_map.insert(connect_var, { u, v });
	return reachLits[u][v];
}
template<typename Weight>
void ConnectedComponentsDetector<Weight>::addConnectedLit(Var outer_reach_var, int node1, int node2) {
	g_under.invalidate();
	g_over.invalidate();
	if (reachLits.size() < g_under.nodes()) {
		reachLits.growTo(g_under.nodes());
		for (int i = 0; i < reachLits.size(); i++) {
			reachLits[i].growTo(g_under.nodes(), lit_Undef);
		}
	}
	
	if (reachLits[node1][node2] == lit_Undef) {
		assert(reachLits[node2][node1] == lit_Undef);
		Var connect_var = outer->newVar(outer_reach_var, getID());
		reachLits[node1][node2] = mkLit(connect_var);
		reachLits[node2][node1] = reachLits[node1][node2];
		reach_map.insert(connect_var, { node1, node2 });
	} else {
		Lit r = reachLits[node1][node2];
		//force equality between the new lit and the old reach lit, in the SAT solver
		outer->makeEqualInSolver(outer->toSolver(r), mkLit(outer_reach_var));
	}
	assert(reachLits[node1][node2] == reachLits[node2][node1]);
}
template<typename Weight>
void ConnectedComponentsDetector<Weight>::addConnectedComponentsLit(Var outer_weight_var, int min_components) {
	
	g_under.invalidate();
	g_over.invalidate();
	//while( dist_lits[to].size()<=within_steps)
	//	dist_lits[to].push({lit_Undef,-1});
	
	Var weight_var = outer->newVar(outer_weight_var, getID());
	
	/*	while(outer->S->nVars()<=weight_var)
	 outer->S->newVar();*/

	Lit reachLit = mkLit(weight_var, false);
	bool found = false;
	for (int i = 0; i < connected_components_lits.size(); i++) {
		if (connected_components_lits[i].min_components == min_components) {
			found = true;
			Lit r = connected_components_lits[i].l;
			//force equality between the new lit and the old reach lit, in the SAT solver
			//outer->S->addClause(~r, reachLit);
			//outer->S->addClause(r, ~reachLit);
			outer->makeEqual(r, reachLit);
		}
	}
	if (!found) {
		connected_components_lits.push();
		connected_components_lits.last().l = reachLit;
		connected_components_lits.last().min_components = min_components;
		
		//weight_lit_map.insert(min_weight,weight_lits.size()-1);
	}
	
}

/*void ConnectedComponentsDetector<Weight>::ConnectedComponentsStatus::setConnected(int edgeid, bool in_tree){
 if(edgeid<detector.tree_edge_lits.size()){
 Lit l = detector.tree_edge_lits[edgeid].l;
 if(l!=lit_Undef)
 detector.changed_edges.push({in_tree? l:~l,edgeid});
 }
 }*/
template<typename Weight>
void ConnectedComponentsDetector<Weight>::ConnectedComponentsStatus::setConnected(int u, int v, bool connected) {
	if (u > v) {
		std::swap(u, v);
	}
	if (detector.reachLits[u][v] != lit_Undef) {
		Lit l = detector.reachLits[u][v];
		lbool assign = detector.outer->value(l);
		if (assign == l_True && connected) {
			//do nothing
		} else if (assign == l_False && !connected) {
			//do nothing
		} else
			detector.changed.push( { connected ? l : ~l, u, v });
	}
}
template<typename Weight>
void ConnectedComponentsDetector<Weight>::ConnectedComponentsStatus::setComponents(int components) {
	
	for (int i = 0; i < detector.connected_components_lits.size(); i++) {
		int min_components = detector.connected_components_lits[i].min_components;
		
		Lit l = detector.connected_components_lits[i].l;
		if (l != lit_Undef) {
			
			assert(l != lit_Undef);
			if (components <= min_components && !polarity) {
				lbool assign = detector.outer->value(l);
				if (assign != l_False) {
					detector.changed_weights.push( { ~l, min_components });
				}
			} else if (components > min_components && polarity) {
				lbool assign = detector.outer->value(l);
				if (assign != l_True) {
					detector.changed_weights.push( { l, min_components });
				}
			}
		}
	}
	
}

template<typename Weight>
void ConnectedComponentsDetector<Weight>::buildMinComponentsTooLowReason(int min_components, vec<Lit> & conflict) {
	
	double starttime = rtime(2);
	
	seen.clear();
	seen.growTo(g_under.nodes());
	
	visit.clear();
	//ok,construct spanning forest from among the connected elements using dfs (we don't need to use kruskal's, as the graph is unweighted), and learn that at least one edge in that each must be disabled.
	for (int k = 0; k < overapprox_component_detector->numComponents(); k++) {
		//learn that at least one edge in the tree must be disabled (else, number of connected components cannot be increased)
		int root = overapprox_component_detector->getComponent(k);
		seen[root] = true;
		visit.push(root);
		
		while (visit.size()) {
			int u = visit.last();
			visit.pop();
			
			for (int i = 0; i < g_under.nIncident(u); i++) {
				int v = g_under.incident(u, i).node;
				int edgeid = g_under.incident(u, i).id;
				if (g_under.edgeEnabled(edgeid) && !seen[v]) {
					seen[v] = true;
					visit.push(v);
					//this edge is in the spanning tree
					Var v = outer->edge_list[edgeid].v;
					assert(outer->value(v)==l_True);
					conflict.push(mkLit(v, true));
				}
			}
		}
	}
	outer->num_learnt_paths++;
	outer->learnt_path_clause_length += (conflict.size() - 1);
	double elapsed = rtime(2) - starttime;
	outer->pathtime += elapsed;
}
template<typename Weight>
void ConnectedComponentsDetector<Weight>::buildMinComponentsTooHighReason(int min_components, vec<Lit> & conflict) {
	static int it = 0;
	++it;
	
	//drawFull( non_reach_detectors[detector]->getSource(),u);
	//assert(outer->dbg_distance( source,u));
	
	double starttime = rtime(2);
	int INF = std::numeric_limits<int>::max();
	underapprox_component_detector->update();
	int numComponents = underapprox_component_detector->numComponents();
	assert(numComponents >= min_components);
	
	//IF the ConnectedComponents is disconnected, the reason is a separating cut between any two disconnected components.
	//we can find one of these by identifying any two roots
	
	int min_conflict = INF;
	
	//walk back down from the each root to find a separating cut of disabled edge.
	//return the smallest such cut.
	
	edge_in_clause.clear();
	edge_in_clause.growTo(g_under.nEdgeIDs());
	assert(conflict.size() == 1);
	for (int i = 0; i < underapprox_component_detector->numComponents(); i++) {
		int root = underapprox_component_detector->getComponent(i);
		
		//now explore that component.
		
		visit.clear();
		//ok, now traverse the connected component in the tree below this root.
		seen.clear();
		seen.growTo(g_under.nodes());
		
		visit.push(root);
		seen[root] = true;
		while (visit.size()) {
			int u = visit.last();
			visit.pop();
			
			for (int i = 0; i < g_over.nIncident(u, true); i++) {
				int edgeid = g_over.incident(u, i).id;
				if (g_over.edgeEnabled(edgeid)) {
					int v = g_over.incident(u, i).node;
					
					if (!seen[v]) {
						//u is a child of node in the minimum spanning tree
						seen[v] = true;
						visit.push(v);
					}
				} else if (!g_over.edgeEnabled(edgeid)) {
					int v = g_over.incident(u, i).node;
					if (!edge_in_clause[edgeid]) {
						edge_in_clause[edgeid] = true;
						Var e = outer->edge_list[edgeid].v;
						
						assert(outer->value(e)==l_False);
						conflict.push(mkLit(e, false));
					}
				}
			}
		}
		
	}
	
	if (opt_components_learn_connect) {
		vec<Lit> c;
		c.push(conflict[0]);
		//for each component, find the lowest valued node; learn that at least one of these lowest nodes must be connected to each other.
		for (int i = 0; i < underapprox_component_detector->numComponents(); i++) {
			int root1 = underapprox_component_detector->getComponent(i);
			for (int j = i + 1; j < underapprox_component_detector->numComponents(); j++) {
				int root2 = underapprox_component_detector->getComponent(j);
				
				Lit connect = getConnectLit(root1, root1);
				c.push(connect);
				
			}
			
		}
		outer->addClauseSafely(c);
	}
	
}
template<typename Weight>
void ConnectedComponentsDetector<Weight>::buildNodesConnectedReason(int source, int node, vec<Lit> & conflict) {
	if (source > node) {
		std::swap(source, node);
	}
	
	UnweightedBFS<Weight,Distance<int>::NullStatus, true> d(source, g_under);
	double starttime = rtime(2);
	d.update();
	
	assert(outer->dbg_reachable(source, node));
	vec<Lit> & reach_lits = reachLits[source];
	
	assert(d.connected_unchecked(node));
	if (opt_learn_reaches == 0 || opt_learn_reaches == 2) {
		int u = node;
		int p;
		while ((p = d.previous(u)) != -1) {
			Edge & edg = outer->edge_list[d.incomingEdge(u)]; //outer->edges[p][u];
			Var e = edg.v;
			lbool val = outer->value(e);
			assert(outer->value(e)==l_True);
			conflict.push(mkLit(e, true));
			u = p;
			
		}
	} else {
		//Instead of a complete path, we can learn reach variables, if they exist
		int u = node;
		int p;
		while ((p = d.previous(u)) != -1) {
			Edge & edg = outer->edge_list[d.incomingEdge(u)]; //outer->edges[p][u];
			Var e = edg.v;
			lbool val = outer->value(e);
			assert(outer->value(e)==l_True);
			conflict.push(mkLit(e, true));
			u = p;
			if (u < reach_lits.size() && reach_lits[u] != lit_Undef && outer->value(reach_lits[u]) == l_True
					&& outer->level(var(reach_lits[u])) < outer->decisionLevel()) {
				//A potential (fixed) problem with the above: reach lit can be false, but have been assigned after r in the trail, messing up clause learning if this is a reason clause...
				//This is avoided by ensuring that L is lower level than the conflict.
				Lit l = reach_lits[u];
				assert(outer->value(l)==l_True);
				conflict.push(~l);
				break;
			}
		}
		
	}
	
	outer->num_learnt_paths++;
	outer->learnt_path_clause_length += (conflict.size() - 1);
	double elapsed = rtime(2) - starttime;
	outer->pathtime += elapsed;
}
template<typename Weight>
void ConnectedComponentsDetector<Weight>::buildNodesNotConnectedReason(int source, int node, vec<Lit> & conflict) {
	if (source > node) {
		std::swap(source, node);
	}
	assert(outer->dbg_notreachable(source, node));
	//assert(!negative_reach_detector->connected_unchecked(node));
	vec<Lit> & reach_lits = reachLits[source];
	double starttime = rtime(2);
	outer->cutGraph.clearHistory();
	outer->stats_mc_calls++;
	/*	if(opt_conflict_min_cut){
	 if(mincutalg!= MinCutAlg::ALG_EDKARP_ADJ){
	 //ok, set the weights for each edge in the cut graph.
	 //Set edges to infinite weight if they are undef or true, and weight 1 otherwise.
	 for(int u = 0;u<outer->cutGraph.nodes();u++){
	 for(int j = 0;j<outer->cutGraph.nIncident(u);j++){
	 int v = outer->cutGraph.incident(u,j).node;
	 int edgeid =  outer->cutGraph.incident(u,j).id;
	 Var var = outer->getEdgeVar(edgeid);
	 if(S->value(var)==l_False){
	 mc.setCapacity(u,v,1);
	 }else{
	 outer->mc->setCapacity(u,v,0xF0F0F0);
	 //}
	 }
	 }

	 //find any edges assigned to false, and set their capacity to 1
	 for(int i =0;i<outer->trail.size();i++){
	 if(outer->trail[i].isEdge && !outer->trail[i].assign){
	 outer->mc->setCapacity(outer->trail[i].from, outer->trail[i].to,1);
	 }
	 }
	 }
	 outer->cut.clear();

	 int f =outer->mc->minCut(source,node,outer->cut);
	 assert(f<0xF0F0F0); assert(f==outer->cut.size());//because edges are only ever infinity or 1
	 for(int i = 0;i<outer->cut.size();i++){
	 MaxFlowEdge e = outer->cut[i];

	 Lit l = mkLit(outer->getEdgeVar(e.id),false);
	 assert(outer->value(l)==l_False);
	 conflict.push(l);
	 }
	 }else*/{
		//We could learn an arbitrary (non-infinite) cut here, or just the whole set of false edges
		//or perhaps we can learn the actual 1-uip cut?
		
		vec<int>& to_visit = outer->to_visit;
		vec<char>& seen = outer->seen;
		
		to_visit.clear();
		to_visit.push(node);
		seen.clear();
		seen.growTo(outer->nNodes());
		seen[node] = true;
		
		do {
			
			assert(to_visit.size());
			int u = to_visit.last();
			assert(u != source);
			to_visit.pop();
			assert(seen[u]);
			assert(!outer->dbg_reachable(source, u, false));
			//assert(!negative_reach_detector->connected_unsafe(u));
			//Ok, then add all its incoming disabled edges to the cut, and visit any unseen, non-disabled incoming.edges()
			for (int i = 0; i < outer->inv_adj[u].size(); i++) {
				int v = outer->inv_adj[u][i].v;
				int from = outer->inv_adj[u][i].from;
				int edge_num = outer->getEdgeID(v);								// v-outer->min_edge_var;
				if (from == u) {
					assert(outer->edge_list[edge_num].to == u);
					assert(outer->edge_list[edge_num].from == u);
					continue;				//Self loops are allowed, but just make sure nothing got flipped around...
				}
				assert(from != u);
				assert(outer->inv_adj[u][i].to == u);
				//Note: the variable has to not only be assigned false, but assigned false earlier in the trail than the reach variable...
				
				if (outer->value(v) == l_False) {
					//note: we know we haven't seen this edge variable before, because we know we haven't visited this node before
					//if we are already planning on visiting the from node, then we don't need to include it in the conflict (is this correct?)
					//if(!seen[from])
					conflict.push(mkLit(v, false));
				} else {
					assert(from != source);
					//even if it is undef? probably...
					if (!seen[from]) {
						seen[from] = true;
						if ((opt_learn_reaches == 2 || opt_learn_reaches == 3) && from < reach_lits.size()
								&& reach_lits[from] != lit_Undef && outer->value(reach_lits[from]) == l_False
								&& outer->level(var(reach_lits[from])) < outer->decisionLevel()) {
							//The problem with the above: reach lit can be false, but have been assigned after r in the trail, messing up clause learning if this is a reason clause...
							Lit r = reach_lits[from];
							assert(var(r) < outer->nVars());
							assert(outer->value(r)==l_False);
							conflict.push(r);
						} else
							to_visit.push(from);
					}
				}
			}
		} while (to_visit.size());
		
	}
}
template<typename Weight>
void ConnectedComponentsDetector<Weight>::buildReason(Lit p, vec<Lit> & reason, CRef marker) {
	
	if (marker == components_low_marker) {
		reason.push(p);
		
		Var v = var(p);
		int weight = -1;
		//could swap this out for a map if there are lots of lits..
		for (int i = 0; i < connected_components_lits.size(); i++) {
			if (var(connected_components_lits[i].l) == v) {
				weight = connected_components_lits[i].min_components;
				break;
			}
		}
		assert(weight >= 0);
		buildMinComponentsTooHighReason(weight, reason);
		
		//double elapsed = rtime(2)-startpathtime;
		//	pathtime+=elapsed;
	} else if (marker == components_high_marker) {
		reason.push(p);
		//the reason is a cut separating p from s;
		//We want to find a min-cut in the full graph separating, where activated edges (ie, those still in antig) are weighted infinity, and all others are weighted 1.
		
		//This is a cut that describes a minimal set of edges which are disabled in the current graph, at least one of which would need to be activated in order for s to reach p
		//assign the mincut edge weights if they aren't already assigned.
		Var v = var(p);
		
		int weight = -1;
		//could swap this out for a map if there are lots of lits..
		for (int i = 0; i < connected_components_lits.size(); i++) {
			if (var(connected_components_lits[i].l) == v) {
				weight = connected_components_lits[i].min_components;
				break;
			}
		}
		assert(weight >= 0);
		buildMinComponentsTooLowReason(weight, reason);
		
	} else if (marker == connected_marker) {
		reason.push(p);
		Var va = var(p);
		auto pair = reach_map[va];
		int u = pair.u;
		int v = pair.v;
		
		buildNodesConnectedReason(u, v, reason);
	} else if (marker == not_connected_marker) {
		reason.push(p);
		Var va = var(p);
		auto pair = reach_map[va];
		int u = pair.u;
		int v = pair.v;
		buildNodesNotConnectedReason(u, v, reason);
	} else {
		assert(false);
	}
}
template<typename Weight>
bool ConnectedComponentsDetector<Weight>::propagate(vec<Lit> & conflict) {
	
	changed.clear();
	changed_weights.clear();
	if (!opt_detect_pure_theory_lits || unassigned_positives > 0) {
		double startdreachtime = rtime(2);
		underapprox_component_detector->update();
		double reachUpdateElapsed = rtime(2) - startdreachtime;
		outer->reachupdatetime += reachUpdateElapsed;
	} else {
		outer->stats_pure_skipped++;
	}
	
	if (!opt_detect_pure_theory_lits || unassigned_negatives > 0) {
		double startunreachtime = rtime(2);
		overapprox_component_detector->update();
		double unreachUpdateElapsed = rtime(2) - startunreachtime;
		outer->unreachupdatetime += unreachUpdateElapsed;
	} else {
		outer->stats_pure_skipped++;
	}
	
	for (int j = 0; j < changed_weights.size(); j++) {
		Lit l = changed_weights[j].l;
		int components = changed_weights[j].min_components;
		bool reach = !sign(l);
		if (outer->value(l) == l_True) {
			//do nothing
		} else if (outer->value(l) == l_Undef) {
			//trail.push(Assignment(false,reach,detectorID,0,var(l)));
			if (reach)
				outer->enqueue(l, components_low_marker);
			else
				outer->enqueue(l, components_high_marker);
			
		} else if (outer->value(l) == l_False) {
			conflict.push(l);
			
			if (reach) {
				
				buildMinComponentsTooHighReason(components, conflict);
				
			} else {
				buildMinComponentsTooLowReason(components, conflict);
			}
			
			return false;
		} else {
			int a = 1;
		}
		
	}
	for (int j = 0; j < changed.size(); j++) {
		Lit l = changed[j].l;
		int u = changed[j].u;
		int v = changed[j].v;
		int components = changed_weights[j].min_components;
		bool reach = !sign(l);
		if (outer->value(l) == l_True) {
			//do nothing
		} else if (outer->value(l) == l_Undef) {
			if (reach)
				outer->enqueue(l, components_low_marker);
			else
				outer->enqueue(l, components_high_marker);
		} else if (outer->value(l) == l_False) {
			conflict.push(l);
			
			if (reach) {
				buildNodesConnectedReason(u, v, conflict);
			} else {
				buildNodesNotConnectedReason(u, v, conflict);
			}
			return false;
		}
		
	}
	return true;
}

template<typename Weight>
void ConnectedComponentsDetector<Weight>::printSolution(std::ostream & write_to) {
	if (opt_verb > 0) {
		int numComponents = underapprox_component_detector->numComponents();
		printf("Number of connected components (graph %d) is: %d\n", outer->getGraphID(), numComponents);
	}
}

template<typename Weight>
bool ConnectedComponentsDetector<Weight>::checkSatisfied() {
	int numConnected = underapprox_component_detector->numComponents();
	int numConnectedOver = overapprox_component_detector->numComponents();
	for (int k = 0; k < connected_components_lits.size(); k++) {
		Lit l = connected_components_lits[k].l;
		int moreThanThisManyComponents = connected_components_lits[k].min_components;
		
		if (l != lit_Undef) {
			
			if (outer->value(l) == l_True) {
				if (overapprox_component_detector->numComponents() <= moreThanThisManyComponents) {
					return false;
				}
			} else if (outer->value(l) == l_False) {
				if (underapprox_component_detector->numComponents() > moreThanThisManyComponents) {
					return false;
				}
			} else {
				if (overapprox_component_detector->numComponents() > moreThanThisManyComponents) {
					return false;
				}
				if (underapprox_component_detector->numComponents() <= moreThanThisManyComponents) {
					return false;
				}
			}
		}
	}
	
	return true;
}
template<typename Weight>
Lit ConnectedComponentsDetector<Weight>::decide() {
	/*ConnectedComponentsDetector *r =this;
	 MinimumSpanningTree<ConnectedComponentsDetector<Weight>::ConnectedComponentsStatus> * over = (MinimumSpanningTree<ConnectedComponentsDetector<Weight>::ConnectedComponentsStatus>*) r->negative_reach_detector;

	 MinimumSpanningTree<ConnectedComponentsDetector<Weight>::ConnectedComponentsStatus> * under = (MinimumSpanningTree<ConnectedComponentsDetector<Weight>::ConnectedComponentsStatus>*) r->positive_reach_detector;

	 //we can probably also do something similar, but with cuts, for nodes that are decided to be unreachable.

	 //ok, for each node that is assigned reachable, but that is not actually reachable in the under approx, decide an edge on a feasible path

	 //this can be obviously more efficient
	 //for(int j = 0;j<nNodes();j++){
	 if(opt_decide_graph_neg){

	 }

	 for(int k = 0;k<dist_lits.size();k++){
	 for(int n = 0;n<dist_lits[k].size();n++){
	 Lit l = dist_lits[k][n].l;
	 int min_weight = dist_lits[k][n].min_weight;
	 if(l==lit_Undef)
	 continue;
	 int j = r->getNode(var(l));
	 if(outer->value(l)==l_True){
	 if(opt_decide_graph_pos){
	 //if(S->level(var(l))>0)
	 //	continue;

	 assert(over->distance(j)<=min_dist);//else we would already be in conflict before this decision was attempted!
	 if(under->distance(j)>min_dist){
	 //then lets try to connect this
	 //static vec<bool> print_path;

	 assert(over->connected(j));//Else, we would already be in conflict


	 int p =j;
	 int last=j;
	 if(!opt_use_random_path_for_decisions)
	 {
	 //ok, read back the path from the over to find a candidate edge we can decide
	 //find the earliest unconnected node on this path
	 over->update();
	 p = j;
	 last = j;
	 int dist = 0;
	 while(under->distance(p)>=min_dist-dist){

	 last=p;
	 assert(p!=r->source);
	 int prev = over->previous(p);
	 assert(over->distance(p)<=min_dist-dist);
	 dist+=1;//should really be weighted
	 p = prev;

	 }
	 }else{
	 //This won't work (without modification) because we need to constrain these paths to ones of maximum real distance < min_dist.
	 //Randomly re-weight the graph sometimes
	 if(drand(rnd_seed)<opt_decide_graph_re_rnd){

	 for(int i=0;i<outer->g.nodes();i++){
	 double w = drand(rnd_seed);
	 w-=0.5;
	 w*=w;
	 //printf("%f (%f),",w,rnd_seed);
	 rnd_path->setWeight(i,w);
	 }
	 }
	 rnd_path->update();
	 //derive a random path in the graph
	 p = j;
	 last = j;
	 assert( rnd_path->connected(p));
	 while(!under->connected(p)){

	 last=p;
	 assert(p!=source);
	 int prev = rnd_path->previous(p);
	 p = prev;
	 assert(p>=0);
	 }

	 }



	 //ok, now pick some edge p->last that will connect p to last;
	 assert(!under->connected(last));
	 assert(under->connected(p));

	 assert(over->connected(last));
	 assert(over->connected(p));
	 assert(p>-1);
	 if(p>-1){
	 Var v = outer->edges[p][last].v;
	 if(outer->value(v)==l_Undef){
	 return mkLit(v,false);
	 }else{
	 assert(outer->value(v)!=l_True);
	 }
	 }
	 for(int k = 0;k<outer->antig.adjacency[p].size();k++){
	 int to = outer->antig.adjacency[p][k].node;
	 if (to==last){
	 Var v =outer->edge_list[ outer->antig.adjacency[p][k].id].v;
	 if(outer->value(v)==l_Undef){
	 return mkLit(v,false);
	 }else{
	 assert(outer->value(v)!=l_True);
	 }
	 }
	 }

	 }
	 }
	 }else if(outer->value(l)==l_False){
	 if(opt_decide_graph_neg){

	 //assert(over->distance(j)<=min_dist);//else we would already be in conflict before this decision was attempted!


	 if(over->distance(j)<=min_dist && under->distance(j)>min_dist){
	 //then lets try to disconnect this node from source by walking back along the path in the over approx, and disabling the first unassigned edge we see.
	 //(there must be at least one such edge, else the variable would be connected in the under approximation as well - in which case it would already have been propagated.
	 int p = j;
	 int last = j;
	 int dist = 0;
	 // tmp_nodes.clear();
	 while(under->distance(p)>min_dist-dist){
	 //int d = under->distance(p);
	 //int d_over = over->distance(p);
	 last=p;
	 assert(p!=source);
	 int prev = over->previous(p);
	 Var v = outer->edges[prev][p].v;
	 if(outer->value(v)==l_Undef){
	 //if(opt_use_random_path_for_decisions)
	 //	tmp_nodes.push(v);
	 //else
	 return mkLit(v,true);
	 }else{
	 assert(outer->value(v)!=l_False);
	 }
	 assert(over->distance(p)<=min_dist-dist);
	 dist+=1;//should really be weighted
	 p = prev;

	 }
	 assert(opt_use_random_path_for_decisions);
	 assert(tmp_nodes.size()>0);
	 int i = irand(rnd_seed,tmp_nodes.size());
	 Var v= tmp_nodes[i];
	 return mkLit(v,true);

	 }
	 }
	 }
	 }
	 }*/
	return lit_Undef;
}
;
template class Monosat::ConnectedComponentsDetector<int> ;
template class Monosat::ConnectedComponentsDetector<long> ;
template class Monosat::ConnectedComponentsDetector<double> ;
#include <gmpxx.h>
template class Monosat::ConnectedComponentsDetector<mpq_class> ;

