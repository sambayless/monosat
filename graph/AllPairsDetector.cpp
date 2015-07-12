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

#include "AllPairsDetector.h"
#include "GraphTheory.h"
#include "dgl/FloydWarshall.h"
#include "dgl/DijkstraAllPairs.h"
#include "dgl/ThorupDynamicConnectivity.h"
using namespace Monosat;
template<typename Weight>
AllPairsDetector<Weight>::AllPairsDetector(int _detectorID, GraphTheorySolver<Weight> * _outer, DynamicGraph<Weight> &_g,
		DynamicGraph<Weight> &_antig, double seed) :
		Detector(_detectorID), outer(_outer), g_under(_g), g_over(_antig), rnd_seed(seed), underapprox_reach_detector(NULL), overapprox_reach_detector(
				NULL), underapprox_path_detector(NULL), positiveReachStatus(NULL), negativeReachStatus(NULL) {
	
	positiveReachStatus = new AllPairsDetector<Weight>::ReachStatus(*this, true);
	negativeReachStatus = new AllPairsDetector<Weight>::ReachStatus(*this, false);
	if (allpairsalg == AllPairsAlg::ALG_FLOYDWARSHALL) {
		underapprox_reach_detector = new FloydWarshall<Weight,AllPairsDetector<Weight>::ReachStatus>(_g,
				*(positiveReachStatus), 1);
		overapprox_reach_detector = new FloydWarshall<Weight,AllPairsDetector<Weight>::ReachStatus>(_antig,
				*(negativeReachStatus), -1);
		underapprox_path_detector = underapprox_reach_detector;
	}/*else if (allpairsalg==ALG_THORUP_ALLPAIRS){
	 positive_reach_detector = new DynamicConnectivity<AllPairsDetector<Weight>::ReachStatus>(_g,*(positiveReachStatus),1);
	 negative_reach_detector = new DynamicConnectivity<AllPairsDetector<Weight>::ReachStatus>(_antig,*(negativeReachStatus),-1);
	 positive_path_detector = positive_reach_detector;
	 }*/else {
		underapprox_reach_detector = new DijkstraAllPairs<Weight,AllPairsDetector<Weight>::ReachStatus>(_g,
				*(positiveReachStatus), 1);
		overapprox_reach_detector = new DijkstraAllPairs<Weight,AllPairsDetector<Weight>::ReachStatus>(_antig,
				*(negativeReachStatus), -1);
		underapprox_path_detector = underapprox_reach_detector;
	}
	/*	if(opt_conflict_shortest_path)
	 reach_detectors.last()->positive_dist_detector = new Dijkstra<PositiveEdgeStatus>(from,g);*/
#ifdef DEBUG_ALLPAIRS
	{	
		dbg_positive_reach_detector = new DijkstraAllPairs<AllPairsDetector<Weight>::IgnoreStatus>(_g,ignoreStatus,1);
		dbg_negative_reach_detector = new DijkstraAllPairs<AllPairsDetector<Weight>::IgnoreStatus>(_antig,ignoreStatus,-1);
	}
#endif
	first_reach_var = var_Undef;
	underprop_marker = outer->newReasonMarker(getID());
	overprop_marker = outer->newReasonMarker(getID());
	
}

template<typename Weight>
void AllPairsDetector<Weight>::addLit(int from, int to, Var outer_reach_var, int within_steps) {
	g_under.invalidate();
	g_over.invalidate();
	
	Var reach_var = outer->newVar(outer_reach_var, getID());
	
	if (first_reach_var == var_Undef) {
		first_reach_var = reach_var;
	} else {
		assert(reach_var > first_reach_var);
	}
	if (within_steps < 0)
		within_steps = outer->nNodes();
	if (within_steps > outer->nNodes())
		within_steps = outer->nNodes();
	
	dist_lits.growTo(outer->nNodes());
	
	dist_lits[from].growTo(outer->nNodes());
	installed_sources.growTo(outer->nNodes());
	
	//while( dist_lits[to].size()<=within_steps)
	//	dist_lits[to].push({lit_Undef,-1});
	
	/*while(outer->S->nVars()<=reach_var)
	 outer->S->newVar();*/

	if (!installed_sources[from]) {
		installed_sources[from] = true;
		sources.push(from);
		underapprox_reach_detector->addSource(from);
		overapprox_reach_detector->addSource(from);
#ifdef DEBUG_ALLPAIRS
		dbg_positive_reach_detector->addSource(from);
		dbg_negative_reach_detector->addSource(from);
#endif
	}
	
	Lit reachLit = mkLit(reach_var, false);
	bool found = false;
	for (int i = 0; i < dist_lits[from][to].size(); i++) {
		if (dist_lits[from][to][i].min_distance == within_steps) {
			found = true;
			Lit r = dist_lits[from][to][i].l;
			//force equality between the new lit and the old reach lit, in the SAT solver
			outer->makeEqual(r, reachLit);
		}
	}
	if (!found) {
		dist_lits[from][to].push();
		dist_lits[from][to].last().l = reachLit;
		dist_lits[from][to].last().min_distance = within_steps;
		dist_lits[from][to].last().source = from;
		while (reach_lit_map.size() <= reach_var - first_reach_var) {
			reach_lit_map.push( { lit_Undef, -1, -1, -1 });
		}
		
		reach_lit_map[reach_var-first_reach_var]= {reachLit,within_steps,from,to};
	}

}
template<typename Weight>
void AllPairsDetector<Weight>::ReachStatus::setReachable(int from, int u, bool reachable) {
	/*	if(polarity==reachable && u<detector.reach_lits.size()){
	 Lit l = detector.reach_lits[u];
	 if(l!=lit_Undef){
	 lbool assign = detector.outer->value(l);
	 if(assign!= (reachable? l_True:l_False )){
	 detector.changed.push({reachable? l:~l,u});
	 }
	 }
	 }*/
}
template<typename Weight>
void AllPairsDetector<Weight>::ReachStatus::setMininumDistance(int from, int u, bool reachable, int distance) {
	assert(reachable == (distance < detector.g_under.nodes()));
	if (distance <= detector.g_under.nodes()) {
		setReachable(from, u, reachable);
	}
	
	if (u < detector.dist_lits.size()) {
		assert(distance >= 0);
		
		for (int i = 0; i < detector.dist_lits[from][u].size(); i++) {
			int d = detector.dist_lits[from][u][i].min_distance;
			
			Lit l = detector.dist_lits[from][u][i].l;
			if (l != lit_Undef) {
				
				assert(l != lit_Undef);
				if (d < distance && !polarity) {
					lbool assign = detector.outer->value(l);
					if (assign != l_False) {
						detector.changed.push( { ~l, u, from });
					}
				} else if (d >= distance && polarity) {
					lbool assign = detector.outer->value(l);
					if (assign != l_True) {
						detector.changed.push( { l, u, from });
					}
				}
			}
		}
	}
	
}

template<typename Weight>
void AllPairsDetector<Weight>::buildReachReason(int source, int to, vec<Lit> & conflict) {
	//drawFull();
	AllPairs & d = *underapprox_path_detector;
	static int iter = 0;
	if (++iter == 3) {
		int a = 1;
	}
	
	double starttime = rtime(2);
	d.update();
	
	tmp_path.clear();
	
	assert(d.connected_unchecked(source, to));
#ifdef DEBUG_ALLPAIRS
	if(!dbg_positive_reach_detector->connected(source,to)) {
		assert(false);
		exit(4);
	}
#endif
	//d.update();
	d.getPath(source, to, tmp_path);
	//if(opt_learn_reaches ==0 || opt_learn_reaches==2)
	{
		int u = to;
		
		//while(( p = d.previous(u)) != -1){
		for (int i = tmp_path.size() - 2; i >= 0; i--) {
			int edge_id = tmp_path[i];
			int p = outer->edge_list[edge_id].from;
			
			Var e = outer->getEdgeVar(edge_id);
			lbool val = outer->value(e);
			assert(outer->value(e)==l_True);
			
			conflict.push(mkLit(e, true));
			u = p;
			
		}
	}
	outer->num_learnt_paths++;
	outer->learnt_path_clause_length += (conflict.size() - 1);
	double elapsed = rtime(2) - starttime;
	outer->pathtime += elapsed;
	
}
template<typename Weight>
void AllPairsDetector<Weight>::buildNonReachReason(int source, int node, vec<Lit> & conflict) {
	static int it = 0;
	++it;
	int u = node;
	//drawFull( non_reach_detectors[detector]->getSource(),u);
	//assert(outer->dbg_distance( source,u));
	double starttime = rtime(2);
	outer->cutGraph.clearHistory();
	outer->stats_mc_calls++;
	{
		
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
			//assert(negative_reach_detector->distance_unsafe(u)>d);
			//Ok, then add all its incoming disabled edges to the cut, and visit any unseen, non-disabled incoming.edges()
			for (int i = 0; i < outer->inv_adj[u].size(); i++) {
				int v = outer->inv_adj[u][i].v;
				int from = outer->inv_adj[u][i].from;
				assert(from != u);
				assert(outer->inv_adj[u][i].to == u);
				//Note: the variable has to not only be assigned false, but assigned false earlier in the trail than the reach variable...
				int edge_num = outer->getEdgeID(v);				    	// v-outer->min_edge_var;
						
				if (outer->value(v) == l_False) {
					//note: we know we haven't seen this edge variable before, because we know we haven't visited this node before
					//if we are already planning on visiting the from node, then we don't need to include it in the conflict (is this correct?)
					//if(!seen[from])
					conflict.push(mkLit(v, false));
				} else if (from != source) {
					//for distance analysis, we _can_ end up reaching source.
					
					//even if it is undef? probably...
					if (!seen[from]) {
						seen[from] = true;
						to_visit.push(from);
					}
				}
			}
		} while (to_visit.size());
		
	}
	
	outer->num_learnt_cuts++;
	outer->learnt_cut_clause_length += (conflict.size() - 1);
	
	double elapsed = rtime(2) - starttime;
	outer->mctime += elapsed;
	
}

template<typename Weight>
void AllPairsDetector<Weight>::buildReason(Lit p, vec<Lit> & reason, CRef marker) {
	
	if (marker == underprop_marker) {
		reason.push(p);
		//	double startpathtime = rtime(2);
		
		/*Dijkstra & detector = *reach_detectors[d]->positive_dist_detector;
		 //the reason is a path from s to p(provided by d)
		 //p is the var for a reachability detector in dijkstra, and corresponds to a node
		 detector.update();
		 Var v = var(p);
		 int u =reach_detectors[d]->getNode(v); //reach_detectors[d]->reach_lit_map[v];
		 assert(detector.connected(u));
		 int w;
		 while(( w= detector.previous(u)) > -1){
		 Lit l = mkLit( edges[w][u].v,true );
		 assert(S->value(l)==l_False);
		 reason.push(l);
		 u=w;
		 }*/
		Var v = var(p);
		int u = getNode(v);
		int source = getSource(v);
		buildReachReason(source, u, reason);
		
		//double elapsed = rtime(2)-startpathtime;
		//	pathtime+=elapsed;
	} else if (marker == overprop_marker) {
		reason.push(p);
		
		//the reason is a cut separating p from s;
		//We want to find a min-cut in the full graph separating, where activated edges (ie, those still in antig) are weighted infinity, and all others are weighted 1.
		
		//This is a cut that describes a minimal set of edges which are disabled in the current graph, at least one of which would need to be activated in order for s to reach p
		//assign the mincut edge weights if they aren't already assigned.
		
		Var v = var(p);
		int t = getNode(v); // v- var(reach_lits[d][0]);
		int source = getSource(v);
		buildNonReachReason(source, t, reason);
		
	} else {
		assert(false);
	}
}
template<typename Weight>
bool AllPairsDetector<Weight>::propagate(vec<Lit> & conflict) {
	
	double startdreachtime = rtime(2);
	getChanged().clear();
	underapprox_reach_detector->update();
	double reachUpdateElapsed = rtime(2) - startdreachtime;
	outer->reachupdatetime += reachUpdateElapsed;
	
	double startunreachtime = rtime(2);
	overapprox_reach_detector->update();
	double unreachUpdateElapsed = rtime(2) - startunreachtime;
	outer->unreachupdatetime += unreachUpdateElapsed;
	
	for (int j = 0; j < getChanged().size(); j++) {
		Lit l = getChanged()[j].l;
		int u = getChanged()[j].u;
		int source = getChanged()[j].source;
		bool reach = !sign(l);
		
		if (outer->value(l) == l_True) {
			//do nothing
		} else if (outer->value(l) == l_Undef) {
#ifdef DEBUG_GRAPH
			assert(outer->dbg_propgation(l));
#endif
#ifdef DEBUG_SOLVER
			if(S->dbg_solver)
			S->dbg_check_propagation(l);
#endif
			//trail.push(Assignment(false,reach,detectorID,0,var(l)));
			if (reach)
				outer->enqueue(l, underprop_marker);
			else
				outer->enqueue(l, overprop_marker);
			
		} else if (outer->value(l) == l_False) {
			conflict.push(l);
			
			if (reach) {
				
				//conflict
				//The reason is a path in g from to s in d
				buildReachReason(source, u, conflict);
				//add it to s
				//return it as a conflict
				
			} else {
				//The reason is a cut separating s from t
				buildNonReachReason(source, u, conflict);
				
			}
#ifdef DEBUG_GRAPH
			for(int i = 0;i<conflict.size();i++)
			assert(outer->value(conflict[i])==l_False);
#endif
#ifdef DEBUG_SOLVER
			if(S->dbg_solver)
			S->dbg_check(conflict);
#endif
			
			return false;
		} else {
			int a = 1;
		}
		
	}
	
#ifdef DEBUG_ALLPAIRS
	for(int m = 0;m<sources.size();m++) {
		int s = sources[m];
		for(int i = 0;i<dist_lits[s].size();i++) {
			for(int j = 0;j<dist_lits[s][i].size();j++) {
				Lit l = dist_lits[s][i][j].l;
				int dist = dist_lits[s][i][j].min_distance;
				if(l!=lit_Undef) {
					int u = getNode(var(l));
					int pos_d = underapprox_reach_detector->distance_unsafe(s,u);
					int neg_d = overapprox_reach_detector->distance_unsafe(s,u);
					int dbg_pos_d=dbg_positive_reach_detector->distance(s,u);
					int dbg_neg_d=dbg_negative_reach_detector->distance(s,u);
					if(dbg_pos_d!=pos_d) {
						assert(false);
						exit(2);
					}
					if(dbg_neg_d!=neg_d) {
						assert(false);
						exit(7);
					}

				}
			}
		}
	}
#endif
	
#ifdef DEBUG_GRAPH
	for(int m = 0;m<sources.size();m++) {
		int s = sources[m];
		for(int i = 0;i<dist_lits[s].size();i++) {
			for(int j = 0;j<dist_lits[s][i].size();j++) {
				Lit l = dist_lits[s][i][j].l;
				int dist = dist_lits[s][i][j].min_distance;
				if(l!=lit_Undef) {
					int u = getNode(var(l));
					if(underapprox_reach_detector->distance_unsafe(s,u)<=dist) {
						assert(outer->value(l)==l_True);
					} else if (overapprox_reach_detector->distance_unsafe(s,u)>dist) {
						int d =overapprox_reach_detector->distance_unsafe(s,u);
						assert(outer->value(l)==l_False);
					}
				}
			}
		}
	}
#endif
	return true;
}
template<typename Weight>
bool AllPairsDetector<Weight>::checkSatisfied() {
	for (int m = 0; m < sources.size(); m++) {
		int s = sources[m];
		for (int j = 0; j < dist_lits[s].size(); j++) {
			for (int k = 0; k < dist_lits[s][j].size(); k++) {
				Lit l = dist_lits[s][j][k].l;
				int dist = dist_lits[s][j][k].min_distance;
				DistLit dl = dist_lits[s][j][k];
				assert(s == dist_lits[s][j][k].source);
				if (l != lit_Undef) {
					int node = getNode(var(l));
					
					if (outer->value(l) == l_True) {
						if (underapprox_reach_detector->distance(s, node) > dist) {
							return false;
						}
					} else if (outer->value(l) == l_False) {
						if (overapprox_reach_detector->distance(s, node) <= dist) {
							return false;
						}
					} else {
						if (underapprox_reach_detector->distance(s, node) <= dist) {
							return false;
						}
						if (overapprox_reach_detector->distance(s, node) > dist) {
							return false;
						}
					}
				}
			}
		}
	}
	return true;
}
template<typename Weight>
Lit AllPairsDetector<Weight>::decide() {
	
	//AllPairs * over = (FloydWarshall<AllPairsDetector<Weight>::ReachStatus>*) negative_reach_detector;
	
	//FloydWarshall<AllPairsDetector<Weight>::ReachStatus> * under = (FloydWarshall<AllPairsDetector<Weight>::ReachStatus>*)positive_reach_detector;
	
	//we can probably also do something similar, but with cuts, for nodes that are decided to be unreachable.
	
	//ok, for each node that is assigned reachable, but that is not actually reachable in the under approx, decide an edge on a feasible path
	
	//this can be obviously more efficient
	//for(int j = 0;j<nNodes();j++){
	for (int m = 0; m < sources.size(); m++) {
		int s = sources[m];
		for (int k = 0; k < dist_lits[s].size(); k++) {
			for (int n = 0; n < dist_lits[s][k].size(); n++) {
				Lit l = dist_lits[s][k][n].l;
				int min_dist = dist_lits[s][k][n].min_distance;
				if (l == lit_Undef)
					continue;
				int j = getNode(var(l));
				if (outer->value(l) == l_True) {
					//if(S->level(var(l))>0)
					//	continue;
					
					assert(overapprox_reach_detector->distance(s, j) <= min_dist);//else we would already be in conflict before this decision was attempted!
					if (underapprox_reach_detector->distance(s, j) > min_dist) {
						//then lets try to connect this
						static vec<bool> print_path;
						
						assert(overapprox_reach_detector->connected(s, j));		//Else, we would already be in conflict
						int p = j;
						int last = j;
						//if(!opt_use_random_path_for_decisions)
						{
							//ok, read back the path from the over to find a candidate edge we can decide
							//find the earliest unconnected node on this path
							overapprox_reach_detector->update();
							assert(overapprox_reach_detector->distance(s, j) <= min_dist);
#ifdef DEBUG_ALLPAIRS
							if(!overapprox_reach_detector->connected(s,j)) {
								assert(false);
								exit(4);
							}

							if(!(overapprox_reach_detector->distance(s,j)<=min_dist)) {
								assert(false);
								exit(4);
							}
#endif
							tmp_path.clear();
							overapprox_reach_detector->getPath(s, j, tmp_path);
							
							p = j;
							last = j;
							
							assert(tmp_path.size());
							int d = 0;
							while (underapprox_reach_detector->distance(s, p) > (min_dist - d)) {
								
								tmp_path.pop_back();
								if (tmp_path.size() == 0)
									break;
								last = p;
								assert(p != s);
								assert(tmp_path.size());
								int prev = tmp_path.back(); //negative_reach_detector->previous(s,p);
								d += 1; //for weighted graphs, this needs to be the weight of the edge.
								p = prev;
#ifdef DEBUG_ALLPAIRS
								if(!overapprox_reach_detector->connected(s,p)) {
									assert(false);

								}
								if(!(overapprox_reach_detector->distance(s,p)<=min_dist-d)) {
									assert(false);

								}
#endif
							}
							
						}
						
						for (int k = 0; k < g_over.nIncident(p); k++) {
							int to = g_over.incident(p, k).node;
							if (to == last) {
								Var v = outer->edge_list[g_over.incident(p, k).id].v;
								if (outer->value(v) == l_Undef) {
									return mkLit(v, false);
								} else {
									assert(outer->value(v)!=l_True);
								}
							}
						}
						
					}
				}
			}
		}
	}
	return lit_Undef;
}
;
template class Monosat::AllPairsDetector<int> ;
template class Monosat::AllPairsDetector<long> ;
template class Monosat::AllPairsDetector<double> ;
#include <gmpxx.h>
template class Monosat::AllPairsDetector<mpq_class> ;

