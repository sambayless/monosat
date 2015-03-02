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

#include "MaxflowDetector.h"
#include "GraphTheory.h"
#include "dgl/EdmondsKarpAdj.h"
#include "dgl/KohliTorr.h"
#include "dgl/EdmondsKarpDynamic.h"
#include "dgl/Dinics.h"
#include "dgl/DinicsLinkCut.h"
using namespace Monosat;

template<>
void MaxflowDetector<int>::buildDinitzLinkCut() {
	underapprox_detector = new DinitzLinkCut<std::vector<int>>(g_under, capacities, source, target);
	overapprox_detector = new DinitzLinkCut<std::vector<int>>(g_over, capacities, source, target);
}

template<typename Weight>
void MaxflowDetector<Weight>::buildDinitzLinkCut() {
	underapprox_detector = nullptr;
	overapprox_detector = nullptr;
	assert(false);
	fprintf(stderr,
			"Warning: Dinitz Link/Cut Tree implementation only supports unweighted or integer weight graphs, aborting!\n");
	exit(1);
}

template<typename Weight>
MaxflowDetector<Weight>::MaxflowDetector(int _detectorID, GraphTheorySolver<Weight> * _outer,
		std::vector<Weight> & capacities, DynamicGraph &_g, DynamicGraph &_antig, int from, int _target, double seed) :
		LevelDetector(_detectorID), outer(_outer), capacities(capacities), g_under(_g), g_over(_antig), source(
				from), target(_target), rnd_seed(seed), cutStatus(*this) {
	
	if (mincutalg == MinCutAlg::ALG_EDKARP_DYN) {
		underapprox_detector = new EdmondsKarpDynamic<std::vector<Weight>, Weight>(_g, capacities, source, target);
		overapprox_detector = new EdmondsKarpDynamic<std::vector<Weight>, Weight>(_antig, capacities, source, target);
		underapprox_conflict_detector = underapprox_detector; //new EdmondsKarpAdj<std::vector<Weight>,Weight>(_g,capacities);
		overapprox_conflict_detector = overapprox_detector; // new EdmondsKarpAdj<std::vector<Weight>,Weight>(_antig,capacities);
		if (opt_conflict_min_cut_maxflow || opt_adaptive_conflict_mincut)
			learn_cut = new EdmondsKarpDynamic<CutStatus, long>(learn_graph, cutStatus, source, target);
		
		/*if(opt_conflict_min_cut_maxflow)
		 learn_cut = new EdmondsKarpAdj<std::vector<int>,int>(learn_graph,learn_caps,source,target);*/
	} else if (mincutalg == MinCutAlg::ALG_EDKARP_ADJ) {
		underapprox_detector = new EdmondsKarpAdj<std::vector<Weight>, Weight>(_g, capacities, source, target);
		overapprox_detector = new EdmondsKarpAdj<std::vector<Weight>, Weight>(_antig, capacities, source, target);
		underapprox_conflict_detector = underapprox_detector;
		overapprox_conflict_detector = overapprox_detector;
		if (opt_conflict_min_cut_maxflow || opt_adaptive_conflict_mincut)
			learn_cut = new EdmondsKarpAdj<CutStatus, long>(learn_graph, cutStatus, source, target);
	} else if (mincutalg == MinCutAlg::ALG_DINITZ) {
		underapprox_detector = new Dinitz<std::vector<Weight>, Weight>(_g, capacities, source, target);
		overapprox_detector = new Dinitz<std::vector<Weight>, Weight>(_antig, capacities, source, target);
		underapprox_conflict_detector = underapprox_detector; // new EdmondsKarpAdj<std::vector<Weight>,Weight>(_g,capacities);
		overapprox_conflict_detector = overapprox_detector; //new EdmondsKarpAdj<std::vector<Weight>,Weight>(_antig,capacities);
		if (opt_conflict_min_cut_maxflow || opt_adaptive_conflict_mincut)
			learn_cut = new Dinitz<CutStatus, long>(learn_graph, cutStatus, source, target);
	} else if (mincutalg == MinCutAlg::ALG_DINITZ_LINKCUT) {
		//link-cut tree currently only supports ints (enforcing this using tempalte specialization...).
		buildDinitzLinkCut();
		underapprox_conflict_detector = new EdmondsKarpAdj<std::vector<Weight>, Weight>(_g, capacities, source, target);
		overapprox_conflict_detector = new EdmondsKarpAdj<std::vector<Weight>, Weight>(_antig, capacities, source,
				target);
		if (opt_conflict_min_cut_maxflow || opt_adaptive_conflict_mincut)
			learn_cut = new EdmondsKarpAdj<CutStatus, long>(learn_graph, cutStatus, source, target);
	} else if (mincutalg == MinCutAlg::ALG_KOHLI_TORR) {
		underapprox_detector = new KohliTorr<std::vector<Weight>, Weight>(_g, capacities, source, target,
				opt_kt_preserve_order);
		overapprox_detector = new KohliTorr<std::vector<Weight>, Weight>(_antig, capacities, source, target,
				opt_kt_preserve_order);
		if (opt_use_kt_for_conflicts) {
			underapprox_conflict_detector = underapprox_detector; //new EdmondsKarpDynamic<std::vector<Weight>,Weight>(_g,capacities);
			overapprox_conflict_detector = overapprox_detector; //new EdmondsKarpDynamic<std::vector<Weight>,Weight>(_antig,capacities);
		} else {
			//for reasons I don't yet understand, kohli-torr seems to produce maxflows that work very poorly as theory-decisions for some problems.
			//Possibly, this is because Kohli-torr will sometimes produces flows that use multiple paths even when a maxflow with a single path is possible.
			underapprox_conflict_detector = new EdmondsKarpDynamic<std::vector<Weight>, Weight>(_g, capacities, source,
					target);
			overapprox_conflict_detector = new EdmondsKarpDynamic<std::vector<Weight>, Weight>(_antig, capacities,
					source, target);
		}
		if (opt_conflict_min_cut_maxflow || opt_adaptive_conflict_mincut)
			learn_cut = new KohliTorr<CutStatus, long>(learn_graph, cutStatus, source, target, opt_kt_preserve_order);
		
	} else {
		underapprox_detector = new EdmondsKarpAdj<std::vector<Weight>, Weight>(_g, capacities, source, target);
		overapprox_detector = new EdmondsKarpAdj<std::vector<Weight>, Weight>(_antig, capacities, source, target);
		underapprox_conflict_detector = underapprox_detector;
		overapprox_conflict_detector = overapprox_detector;
		if (opt_conflict_min_cut_maxflow || opt_adaptive_conflict_mincut)
			learn_cut = new EdmondsKarpAdj<CutStatus, long>(learn_graph, cutStatus, source, target);
	}
	
#ifdef RECORD
	{
		char t[30];
		sprintf(t, "LEARN_GRAPH%d", getID());
		learn_graph.outfile = fopen(t, "w");
	}
	
#endif
	
	first_reach_var = var_Undef;
	underprop_marker = outer->newReasonMarker(getID());
	overprop_marker = outer->newReasonMarker(getID());
	
}

template<typename Weight>
void MaxflowDetector<Weight>::addFlowLit(Weight maxflow, Var outer_reach_var, bool inclusive) {
	g_under.invalidate();
	g_over.invalidate();
	Var reach_var = outer->newVar(outer_reach_var, getID());
	if (first_reach_var == var_Undef) {
		first_reach_var = reach_var;
	} else {
		assert(reach_var > first_reach_var);
	}
	
	assert(maxflow >= 0);
	
	//while( dist_lits[to].size()<=within_steps)
	//	dist_lits[to].push({lit_Undef,-1});
	
	/*	while(outer->S->nVars()<=reach_var)
	 outer->S->newVar();*/

	Lit reachLit = mkLit(reach_var, false);
	bool found = false;
	
	flow_lits.push();
	flow_lits.last().l = reachLit;
	flow_lits.last().max_flow = maxflow;
	flow_lits.last().inclusive=inclusive;
	while (reach_lit_map.size() <= reach_var - first_reach_var) {
		reach_lit_map.push(-1);
	}
	
	reach_lit_map[reach_var - first_reach_var] = flow_lits.size() - 1;
	
}

template<typename Weight>
void MaxflowDetector<Weight>::buildMaxFlowTooHighReason(Weight flow, vec<Lit> & conflict) {
	//drawFull();
	double starttime = rtime(2);
	tmp_cut.clear();
	Weight actual_flow = underapprox_conflict_detector->maxFlow();
	
	//just collect the set of edges which have non-zero flow, and return them
	//(Or, I could return a cut, probably)
	for (int i = 0; i < g_under.edges(); i++) {
		if (g_under.edgeEnabled(i)) {
			if (underapprox_conflict_detector->getEdgeFlow(i) > 0) {
				Var v = outer->edge_list[i].v;
				assert(outer->value(v)==l_True);
				conflict.push(mkLit(v, true));
			}
		}
	}
	
#ifdef RECORD
	std::sort(conflict.begin(), conflict.end());
	if (g_under.outfile) {
		fprintf(g_under.outfile, "toohigh ");
		for (int i = 0; i < conflict.size(); i++) {
			fprintf(g_under.outfile, "%d,", dimacs(conflict[i]));
		}
		fprintf(g_under.outfile, "\n");
		fflush(g_under.outfile);
	}
	if (g_over.outfile) {
		fprintf(g_over.outfile, "toohigh ");
		for (int i = 0; i < conflict.size(); i++) {
			fprintf(g_over.outfile, "%d,", dimacs(conflict[i]));
		}
		fprintf(g_over.outfile, "\n");
		fflush(g_over.outfile);
	}
#endif
	stats_under_conflicts++;
	outer->num_learnt_paths++;
	outer->learnt_path_clause_length += (conflict.size() - 1);
	double elapsed = rtime(2) - starttime;
	outer->pathtime += elapsed;
	
	stats_under_conflict_time += elapsed;
}
/*

 template<typename Weight>
 int MaxflowDetector<Weight>::dbg_minconflict(){
 #ifndef NDEBUG
 Weight foundflow = negative_conflict_detector->maxFlow();

 int INF=0x0FF0F0;
 //for each edge in the original graph, need to add a forward and backward edge here.
 if(learn_graph.nodes()<g.nodes()){
 while(learn_graph.nodes()<g.nodes())
 learn_graph.addNode();

 for (auto & e:g.all_edges){
 learn_graph.addEdge(e.from,e.to);
 }
 back_edges.growTo(g.edges());
 for (auto & e:g.all_edges){
 back_edges[e.id] = learn_graph.addEdge(e.to,e.from);
 }
 learn_caps.resize(learn_graph.edges());
 }

 //now, set learn_graph to the residual graph
 for (auto & e:g.all_edges){
 int from = e.from;
 int to = e.to;
 int v = outer->getEdgeVar(e.id);
 lbool val =outer->value(v);
 if(val!=l_False){
 Weight flow = negative_conflict_detector->getEdgeFlow(e.id);
 Weight capacity =  negative_conflict_detector->getEdgeCapacity(e.id);
 if(capacity==0){
 learn_graph.disableEdge(e.id);
 learn_graph.disableEdge(back_edges[e.id]);
 }else{
 if(flow>0){
 //then there is capacity in the backward edge in the residual graph
 int back_edge = back_edges[e.id];
 learn_graph.enableEdge(back_edges[e.id]);
 learn_caps[back_edge] = INF;
 }else{
 learn_graph.disableEdge(back_edges[e.id]);
 }
 if (flow<capacity){
 //then there is capacity in the forward edge in the residual graph
 learn_caps[e.id] = INF;
 learn_graph.enableEdge(e.id);
 }else{
 learn_graph.disableEdge(e.id);
 }
 }
 }else{

 learn_graph.enableEdge(e.id);
 learn_graph.disableEdge(back_edges[e.id]);
 learn_caps[e.id] = 1;
 }
 }
 learn_graph.invalidate();
 cut.clear();

 auto check_cut = new EdmondsKarpAdj<std::vector<int>,int>(learn_graph,learn_caps,source,target);
 int f =check_cut->minCut(cut);
 return cut.size();
 #endif
 return 0;
 }
 */
void bassert(bool condition) {
#ifndef NDEBUG
	assert(condition);
	if (!condition) {
		exit(3);
	}
#endif
}

template<typename Weight>
void MaxflowDetector<Weight>::buildMaxFlowTooLowReason(Weight maxflow, vec<Lit> & conflict, bool force_maxflow) {
	static int it = 0;
	++it;
	if (it == 136) {
		int a = 1;
	}
	//printf("%d\n",it);
	double starttime = rtime(2);
	if (force_maxflow || opt_conflict_min_cut_maxflow) {
		Weight foundflow = overapprox_conflict_detector->maxFlow();
		
		/*	{

		 EdmondsKarpAdj<std::vector<Weight>,Weight> ek(antig,capacities,source,target);
		 std::vector<MaxFlowEdge> tmpcut;
		 Weight tf = ek.maxFlow();

		 if(tf != foundflow){
		 exit(3);
		 }


		 }*/

		collectChangedEdges();
		collectDisabledEdges();
#ifndef NDEBUG
		for (auto & e : g_over.all_edges) {
			int from = e.from;
			int to = e.to;
			int v = outer->getEdgeVar(e.id);
			lbool val = outer->value(v);
			if (val != l_False) {
				Weight flow = overapprox_conflict_detector->getEdgeFlow(e.id);
				Weight capacity = overapprox_conflict_detector->getEdgeCapacity(e.id);
				if (capacity == 0) {
					bassert(!learn_graph.edgeEnabled(e.id * 2));
					bassert(!learn_graph.edgeEnabled(e.id * 2 + 1));
					bassert(!learn_graph.edgeEnabled(back_edges[e.id * 2]));
					bassert(!learn_graph.edgeEnabled(back_edges[e.id * 2 + 1]));
					
				} else {
					if (flow > 0) {
						//then there is capacity in the backward edge in the residual graph
						int back_edge = back_edges[e.id];
						/*	learn_graph.enableEdge(back_edges[e.id]);
						 learn_caps[back_edge] = INF;*/
						bassert(!learn_graph.edgeEnabled(back_edges[e.id * 2]));
						bassert(learn_graph.edgeEnabled(back_edges[e.id * 2 + 1]));
					} else {
						//learn_graph.disableEdge(back_edges[e.id]);
						bassert(!learn_graph.edgeEnabled(back_edges[e.id * 2]));
						bassert(!learn_graph.edgeEnabled(back_edges[e.id * 2 + 1]));
					}
					if (flow < capacity) {
						//then there is capacity in the forward edge in the residual graph
						/*	learn_caps[e.id] = INF;
						 learn_graph.enableEdge(e.id);*/
						bassert(!learn_graph.edgeEnabled(e.id * 2));
						bassert(learn_graph.edgeEnabled(e.id * 2 + 1));
					} else {
						//learn_graph.disableEdge(e.id);
						bassert(!learn_graph.edgeEnabled(e.id * 2));
						bassert(!learn_graph.edgeEnabled(e.id * 2 + 1));
					}
				}
			} else {
				
				/*			learn_graph.enableEdge(e.id);
				 learn_graph.disableEdge(back_edges[e.id]);
				 learn_caps[e.id] = 1;*/
				bassert(learn_graph.edgeEnabled(e.id * 2));
				bassert(!learn_graph.edgeEnabled(e.id * 2 + 1));
				bassert(!learn_graph.edgeEnabled(back_edges[e.id * 2]));
				bassert(!learn_graph.edgeEnabled(back_edges[e.id * 2 + 1]));
			}
		}
#endif
		
		long f = learn_cut->minCut(cut);
		learn_graph.clearChanged();
		learn_graph.clearHistory();
		/*				{

		 EdmondsKarpAdj<CutStatus,long> ek(learn_graph, cutStatus,source,target);
		 std::vector<MaxFlowEdge> tmpcut;
		 long tf = ek.minCut(tmpcut);
		 printf("cut size:%d, %d, expected: %d, %d \n",cut.size(),f, tmpcut.size(), tf);
		 if(f != tf || cut.size()!= tmpcut.size()){
		 exit(3);
		 }


		 }*/
		if (f < cutStatus.inf) {
			assert(f < 0xF0F0F0);
			assert(f == cut.size());			//because edges are only ever infinity or 1
			for (int i = 0; i < cut.size(); i++) {
				MaxFlowEdge e = cut[i];
				int edgeID = e.id / 2;
				Lit l = mkLit(outer->getEdgeVar(edgeID), false);
				bassert(outer->value(l) == l_False);
				conflict.push(l);
			}
		} else {
			int a = 1;
			
			//there is no way to increase the max flow.
		}
		
		outer->num_learnt_cuts++;
		outer->learnt_cut_clause_length += (conflict.size() - 1);
		stats_over_conflicts++;
		double elapsed = rtime(2) - starttime;
		stats_over_conflict_time += elapsed;
		
		return;
	}
	
	//drawFull( non_reach_detectors[detector]->getSource(),u);
	//assert(outer->dbg_distance( source,u));
	
	//The reason why we can't reach this assignment is a cut through the disabled edges in the residual graph from the overapprox.
	//we could search for a min cut, but instead we will just step back in the RESIDUAL graph, from u to s, collecting disabled edges.
	seen.clear();
	seen.growTo(outer->nNodes());
	visit.clear();
	Weight foundflow = overapprox_conflict_detector->maxFlow();
	/*			std::vector<MaxFlowEdge> ignore;
	 negative_conflict_detector->minCut(ignore);*/
	visit.push(target);
	for (int k = 0; k < visit.size(); k++) {
		int u = visit[k];
		for (int i = 0; i < g_under.nIncoming(u); i++) {
			int p = g_under.incoming(u, i).node;
			if (!seen[p]) {
				
				int edgeid = g_under.incoming(u, i).id;
				int v = outer->getEdgeVar(edgeid);
				
				//assert( g.incoming(u,i).to==u);
				if (outer->value(v) != l_False) {
					//this is an enabled edge in the overapprox
					
					Weight residual_capacity = overapprox_conflict_detector->getEdgeResidualCapacity(edgeid);
					if (residual_capacity > 0) {
						seen[p] = true;
						visit.push(p);
					}
				} else {
					//this is a disabled edge, and we can add it to the cut.
					//we're going to assume the edge has non-zero capacity here, otherwise we could exclude it (but it shouldn't really even be in this graph in that case, anyways).
					conflict.push(mkLit(v, false));
				}
			}
			//pred
		}
		//we are walking back in the RESIDUAL graph, which can mean backwards traversal of edges with flow!
		for (int i = 0; i < g_under.nIncident(u); i++) {
			int p = g_under.incident(u, i).node;
			if (!seen[p]) {
				int edgeid = g_under.incident(u, i).id;
				int v = outer->getEdgeVar(edgeid);
				
				//assert( g.incoming(u,i).to==u);
				if (outer->value(v) != l_False) {
					//this is the residual capacity of the backwards edge in the residual graph - which is equal to the forwards flow on this edge!
					Weight residual_capacity = overapprox_conflict_detector->getEdgeFlow(edgeid);
					if (residual_capacity) {
						seen[p] = true;
						visit.push(p);
					}
				} else {
					//this is a disabled edge, and we can add it to the cut.
					//we're going to assume the edge has non-zero capacity here, otherwise we could exclude it (but it shouldn't really even be in this graph in that case, anyways).
					conflict.push(mkLit(v, false));
				}
			}
			//pred
		}
	}
	
	if (!force_maxflow && opt_adaptive_conflict_mincut > 0 && (conflict.size() - 1 > opt_adaptive_conflict_mincut)) { //-1 to ignore the predicate's literal stored at position 0
		conflict.shrink(conflict.size() - 1);
		assert(conflict.size() == 1);
		buildMaxFlowTooLowReason(maxflow, conflict, true);
		return;
	}
	
	/*if(conflict.size()<dbg_minconflict()){
	 exit(4);
	 }*/
#ifdef RECORD
	if (overapprox_detector != overapprox_conflict_detector) {
		Weight foundflow2 = overapprox_detector->maxFlow();
		assert(foundflow == foundflow2);
	}
	std::sort(conflict.begin(), conflict.end());
	if (g_under.outfile) {
		fprintf(g_under.outfile, "toolow ");
		for (int i = 0; i < conflict.size(); i++) {
			fprintf(g_under.outfile, "%d,", dimacs(conflict[i]));
		}
		fprintf(g_under.outfile, "\n");
		fflush(g_under.outfile);
	}
	if (g_over.outfile) {
		fprintf(g_over.outfile, "toolow ");
		for (int i = 0; i < conflict.size(); i++) {
			fprintf(g_over.outfile, "%d,", dimacs(conflict[i]));
		}
		fprintf(g_over.outfile, "\n");
		fflush(g_over.outfile);
	}
#endif
	outer->num_learnt_cuts++;
	outer->learnt_cut_clause_length += (conflict.size() - 1);
	stats_over_conflicts++;
	double elapsed = rtime(2) - starttime;
	stats_over_conflict_time += elapsed;
	
}

template<typename Weight>
void MaxflowDetector<Weight>::buildReason(Lit p, vec<Lit> & reason, CRef marker) {
	
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
		Weight flow = flow_lits[reach_lit_map[v - first_reach_var]].max_flow;
		//int u =getNode(v);
		buildMaxFlowTooHighReason(flow, reason);
		
		//double elapsed = rtime(2)-startpathtime;
		//	pathtime+=elapsed;
	} else if (marker == overprop_marker) {
		reason.push(p);
		
		//the reason is a cut separating p from s;
		//We want to find a min-cut in the full graph separating, where activated edges (ie, those still in antig) are weighted infinity, and all others are weighted 1.
		
		//This is a cut that describes a minimal set of edges which are disabled in the current graph, at least one of which would need to be activated in order for s to reach p
		//assign the mincut edge weights if they aren't already assigned.
		
		Var v = var(p);
		Weight flow = flow_lits[reach_lit_map[v - first_reach_var]].max_flow;
		//int t = getNode(v); // v- var(reach_lits[d][0]);
		buildMaxFlowTooLowReason(flow, reason);
		
	} else {
		assert(false);
	}
}
template<typename Weight>
bool MaxflowDetector<Weight>::propagate(vec<Lit> & conflict) {
	if (flow_lits.size() == 0) {
		return true;
	}
	static int iter1 = 0;
	if (++iter1 == 135) {
		int a = 1;
	}
#ifdef RECORD
	if (g_under.outfile) {
		fprintf(g_under.outfile, "iter %d\n", iter1);
		fflush(g_under.outfile);
	}
	if (g_over.outfile) {
		fprintf(g_over.outfile, "iter %d\n", iter1);
		fflush(g_over.outfile);
	}
#endif
	
	Weight over_maxflow = -1;
	Weight under_maxflow = -1;
	
	if (underapprox_detector && (!opt_detect_pure_theory_lits || unassigned_positives > 0)) {
		double startdreachtime = rtime(2);
		stats_under_updates++;
		under_maxflow = underapprox_detector->maxFlow();
		assert(under_maxflow == underapprox_conflict_detector->maxFlow());
		double reachUpdateElapsed = rtime(2) - startdreachtime;
		stats_under_update_time += reachUpdateElapsed;
	} else
		stats_skipped_under_updates++;
	
	if (overapprox_detector && (!opt_detect_pure_theory_lits || unassigned_negatives > 0)) {
		double startunreachtime = rtime(2);
		stats_over_updates++;
		over_maxflow = overapprox_detector->maxFlow();
		assert(over_maxflow == overapprox_conflict_detector->maxFlow());
		double unreachUpdateElapsed = rtime(2) - startunreachtime;
		stats_over_update_time += unreachUpdateElapsed;
	} else
		stats_skipped_over_updates++;
	
	for (int j = 0; j < flow_lits.size(); j++) {
		
		DistLit f = flow_lits[j];
		Lit l = f.l;
		bool inclusive = f.inclusive;
		Weight& maxflow = f.max_flow;
		//int u = getNode(var(l));
		
		if ((inclusive && under_maxflow >= maxflow) || (!inclusive && under_maxflow > maxflow)) {
			if (outer->value(l) == l_True) {
				//do nothing
			} else if (outer->value(l) == l_Undef) {
				//trail.push(Assignment(false,true,detectorID,0,var(l)));
				outer->enqueue(l, underprop_marker);
				
			} else if (outer->value(l) == l_False) {
				conflict.push(l);
				buildMaxFlowTooHighReason(maxflow, conflict);
				return false;
			}
			
		} else if ((inclusive && over_maxflow < maxflow) || (!inclusive && over_maxflow <= maxflow)) {
			if (outer->value(l) == l_False) {
				//do nothing
			} else if (outer->value(l) == l_Undef) {
				//trail.push(Assignment(false,false,detectorID,0,var(l)));
				outer->enqueue(~l, underprop_marker);
				
			} else if (outer->value(l) == l_True) {
				conflict.push(~l);
				buildMaxFlowTooLowReason(maxflow, conflict);
				return false;
			}
			
		}
		
	}
	
	return true;
}
template<typename Weight>
bool MaxflowDetector<Weight>::checkSatisfied() {
	EdmondsKarpAdj<std::vector<Weight>, Weight> positiveCheck(g_under, capacities, source, target);
	EdmondsKarpAdj<std::vector<Weight>, Weight> negativeCheck(g_over, capacities, source, target);
	for (int j = 0; j < flow_lits.size(); j++) {
		
		Lit l = flow_lits[j].l;
		Weight dist = flow_lits[j].max_flow;
		
		if (l != lit_Undef) {
			//int node =getNode(var(l));
			
			if (outer->value(l) == l_True) {
				if (positiveCheck.maxFlow() < dist) {
					return false;
				}
			} else if (outer->value(l) == l_False) {
				if (negativeCheck.maxFlow() >= dist) {
					return false;
				}
			} else {
				if (positiveCheck.maxFlow() >= dist) {
					return false;
				}
				if (negativeCheck.maxFlow() < dist) {
					return false;
				}
			}
		}
		
	}
	return true;
}
template<typename Weight>
void MaxflowDetector<Weight>::printSolution(std::ostream & write_to) {
	Weight f = underapprox_conflict_detector->maxFlow();
	write_to << "Graph " << outer->getGraphID() << " maxflow " << source << " to " << target << " is " << f << "\n";
	
	int maxw = log10(outer->nNodes()) + 1;
	int width = sqrt(outer->nNodes());
	if (opt_width > 0) {
		width = opt_width;
	}
	int height = width;
	if (opt_height > 0) {
		height = opt_height;
	}
	
	int lasty = 0;
	int extra = outer->nNodes() % width ? (width - outer->nNodes() % width) : 0;
	for (int n = 0; n < outer->nNodes(); n++) {
		int x = n % width;
		
		int y = (n) / width;
		if (y > lasty)
			write_to << "\n";
		Weight total_flow = 0;
		for (int e = 0; e < g_under.edges(); e++) {
			if (g_under.getEdge(e).to == n && g_under.edgeEnabled(e)) {
				total_flow += underapprox_conflict_detector->getEdgeFlow(e);
				
			}
		}
		
		//printf("%*d ",maxw,total_flow);
		write_to << total_flow << " ";
		lasty = y;
	}
	write_to << "\n";
	
	for (int n = 0; n < outer->nNodes(); n++) {
		
		int total_flow = 0;
		for (int e = 0; e < g_under.edges(); e++) {
			if (g_under.getEdge(e).to == n && g_under.edgeEnabled(e)) {
				Weight flow = underapprox_conflict_detector->getEdgeFlow(e);
				if (flow > 0) {
					write_to << "Graph " << outer->getGraphID() << " maxflow " << source << " to " << target
							<< " assigns edge " << g_under.getEdge(e).from << " -> " << g_under.getEdge(e).to << " flow " << flow
							<< "\n";
				}
			}
		}
		
	}
	write_to << "\n";
}

template<typename Weight>
void MaxflowDetector<Weight>::collectDisabledEdges() {
	if (opt_conflict_min_cut_maxflow) {
		if (learn_graph.nodes() < g_under.nodes()) {
			while (learn_graph.nodes() < g_under.nodes())
				learn_graph.addNode();
			back_edges.growTo(g_under.edges() * 2, -1);
			
			//add every edge twice
			for (auto & e : g_under.all_edges) {
				if (back_edges[e.id * 2] == -1) {
					learn_graph.addEdge(e.from, e.to);
					learn_graph.addEdge(e.from, e.to);
				}
			}
			
			for (auto & e : g_under.all_edges) {
				if (back_edges[e.id * 2] == -1) {
					back_edges[e.id * 2] = learn_graph.addEdge(e.to, e.from);
					back_edges[e.id * 2 + 1] = learn_graph.addEdge(e.to, e.from);
					
				}
			}
			for (auto & e : learn_graph.all_edges) {
				learn_graph.disableEdge(e.id);
			}
			learn_graph.invalidate();
			
#ifdef RECORD
			if (learn_graph.outfile) {
				for (int edgeID = 0; edgeID < learn_graph.edges(); edgeID++) {
					fprintf(learn_graph.outfile, "edge_weight %d %ld\n", edgeID, cutStatus[edgeID]);
				}
				fflush(learn_graph.outfile);
			}
#endif
		}
		
		if (learngraph_history_clears != g_over.historyclears || g_over.changed()) {
			//refresh
			for (int edgeid = 0; edgeid < g_over.edges(); edgeid++) {
				if (edgeid == 34) {
					int a = 1;
				}
				if (g_over.hasEdge(edgeid) && !g_over.edgeEnabled(edgeid)) {
					learn_graph.enableEdge(edgeid * 2);
					learn_graph.disableEdge(back_edges[edgeid * 2]);
					learn_graph.disableEdge(edgeid * 2 + 1);
					learn_graph.disableEdge(back_edges[edgeid * 2 + 1]);
				} else {
					Weight flow = overapprox_conflict_detector->getEdgeFlow(edgeid);
					Weight capacity = overapprox_conflict_detector->getEdgeCapacity(edgeid);
					if (capacity == 0) {
						//learn_graph.disableEdge(e.id);
						//learn_graph.disableEdge(back_edges[e.id]);
						
						learn_graph.disableEdge(edgeid * 2);
						learn_graph.disableEdge(back_edges[edgeid * 2]);
						learn_graph.disableEdge(edgeid * 2 + 1);
						learn_graph.disableEdge(back_edges[edgeid * 2 + 1]);
					} else {
						if (flow > 0) {
							//then there is capacity in the backward edge in the residual graph
							/*					int back_edge = back_edges[e.id];
							 learn_graph.enableEdge(back_edges[e.id]);
							 learn_caps[back_edge] = INF;*/

							learn_graph.disableEdge(back_edges[edgeid * 2]);
							learn_graph.enableEdge(back_edges[edgeid * 2 + 1]);
						} else {
							learn_graph.disableEdge(back_edges[edgeid * 2]);
							learn_graph.disableEdge(back_edges[edgeid * 2 + 1]);
						}
						if (flow < capacity) {
							//then there is capacity in the forward edge in the residual graph
							/*	learn_caps[e.id] = INF;
							 learn_graph.enableEdge(e.id);
							 */
							learn_graph.disableEdge(edgeid * 2);
							//learn_graph.disableEdge(back_edges[edgeid*2]);
							learn_graph.enableEdge(edgeid * 2 + 1);
							//learn_graph.disableEdge(back_edges[edgeid*2+1]);
						} else {
							assert(flow == capacity);
							//learn_graph.disableEdge(e.id);
							
							learn_graph.disableEdge(edgeid * 2);
							//learn_graph.disableEdge(back_edges[edgeid*2]);
							learn_graph.disableEdge(edgeid * 2 + 1);
							//learn_graph.disableEdge(back_edges[edgeid*2+1]);
						}
					}
				}
			}
			learngraph_history_clears = g_over.historyclears;
			learngraph_history_qhead = g_over.history.size();
		} else {
			for (int i = learngraph_history_qhead; i < g_over.history.size(); i++) {
				int edgeid = g_over.history[i].id;
				if (edgeid == 34) {
					int a = 1;
				}
				/*	if(antig.selfLoop(edgeid))
				 continue;//skip self loops*/
				if (g_over.history[i].addition && g_over.edgeEnabled(edgeid)) {
					Weight flow = overapprox_conflict_detector->getEdgeFlow(edgeid);
					Weight capacity = overapprox_conflict_detector->getEdgeCapacity(edgeid);
					if (capacity == 0) {
						//learn_graph.disableEdge(e.id);
						//learn_graph.disableEdge(back_edges[e.id]);
						
						learn_graph.disableEdge(edgeid * 2);
						learn_graph.disableEdge(back_edges[edgeid * 2]);
						learn_graph.disableEdge(edgeid * 2 + 1);
						learn_graph.disableEdge(back_edges[edgeid * 2 + 1]);
					} else {
						if (flow > 0) {
							//then there is capacity in the backward edge in the residual graph
							/*					int back_edge = back_edges[e.id];
							 learn_graph.enableEdge(back_edges[e.id]);
							 learn_caps[back_edge] = INF;*/

							learn_graph.disableEdge(back_edges[edgeid * 2]);
							learn_graph.enableEdge(back_edges[edgeid * 2 + 1]);
						} else {
							learn_graph.disableEdge(back_edges[edgeid * 2]);
							learn_graph.disableEdge(back_edges[edgeid * 2 + 1]);
						}
						if (flow < capacity) {
							//then there is capacity in the forward edge in the residual graph
							/*	learn_caps[e.id] = INF;
							 learn_graph.enableEdge(e.id);
							 */
							learn_graph.disableEdge(edgeid * 2);
							//learn_graph.disableEdge(back_edges[edgeid*2]);
							learn_graph.enableEdge(edgeid * 2 + 1);
							//learn_graph.disableEdge(back_edges[edgeid*2+1]);
						} else {
							assert(flow == capacity);
							//learn_graph.disableEdge(e.id);
							
							learn_graph.disableEdge(edgeid * 2);
							//learn_graph.disableEdge(back_edges[edgeid*2]);
							learn_graph.disableEdge(edgeid * 2 + 1);
							//learn_graph.disableEdge(back_edges[edgeid*2+1]);
						}
					}
				} else if (!g_over.history[i].addition && !g_over.edgeEnabled(edgeid)) {
					//any edge that is disabled has a capacity of 1 in the learn graph
					learn_graph.enableEdge(edgeid * 2);
					learn_graph.disableEdge(back_edges[edgeid * 2]);
					learn_graph.disableEdge(edgeid * 2 + 1);
					learn_graph.disableEdge(back_edges[edgeid * 2 + 1]);
				}
			}
			learngraph_history_qhead = g_over.history.size();
		}
	}
}

template<typename Weight>
void MaxflowDetector<Weight>::collectChangedEdges() {
	
	if (is_potential_decision.size() < g_under.edges()) {
		is_potential_decision.growTo(g_under.edges(), false);
	}
	dbg_decisions();
	overapprox_conflict_detector->update();
	
	if (opt_conflict_min_cut_maxflow) {
		if (learn_graph.nodes() < g_under.nodes()) {
			while (learn_graph.nodes() < g_under.nodes())
				learn_graph.addNode();
			back_edges.growTo(g_under.edges() * 2, -1);
			
			//add every edge twice
			for (auto & e : g_under.all_edges) {
				if (back_edges[e.id * 2] == -1) {
					learn_graph.addEdge(e.from, e.to);
					learn_graph.addEdge(e.from, e.to);
				}
			}
			
			for (auto & e : g_under.all_edges) {
				if (back_edges[e.id * 2] == -1) {
					back_edges[e.id * 2] = learn_graph.addEdge(e.to, e.from);
					back_edges[e.id * 2 + 1] = learn_graph.addEdge(e.to, e.from);
				}
			}
			for (auto & e : learn_graph.all_edges) {
				learn_graph.disableEdge(e.id);
			}
			learn_graph.invalidate();
			
#ifdef RECORD
			if (learn_graph.outfile) {
				for (int edgeID = 0; edgeID < learn_graph.edges(); edgeID++) {
					fprintf(learn_graph.outfile, "edge_weight %d %ld\n", edgeID, cutStatus[edgeID]);
				}
				fflush(learn_graph.outfile);
			}
#endif
		}
	}
	
	std::vector<int> & changed_edges = overapprox_conflict_detector->getChangedEdges();
	if (opt_rnd_order_graph_decisions) {
		/*			static vec<int> tmp_changed;
		 tmp_changed.clear();
		 for(int edge:changed_edges):
		 tmp_changed.p*/
		randomShuffle(rnd_seed, changed_edges);
	}
	//for(int j = 0;j<changed_edges.size();j++){
	//			int edgeid = changed_edges[j];
	for (int j = changed_edges.size() - 1; j >= 0; j--) {
		int edgeid = changed_edges[j];
		
		if (!is_potential_decision[edgeid]) {
			Lit l = mkLit(outer->getEdgeVar(edgeid), false);
			if (overapprox_conflict_detector->getEdgeFlow(edgeid) > 0) {
				is_potential_decision[edgeid] = true;
				if (opt_maxflow_decisions_q == 0) {
					potential_decisions_q.insertBack(edgeid);
				} else if (opt_maxflow_decisions_q == 1) {
					potential_decisions_q.insert(edgeid);		//insertBack
				} else if (opt_maxflow_decisions_q == 2) {
					potential_decisions_q.insert(edgeid);
				} else if (opt_maxflow_decisions_q == 3) {
					potential_decisions_q.insertBack(edgeid);
				} else if (opt_maxflow_decisions_q == 4) {
					potential_decisions_q.insertBack(edgeid);
				}
			}
		}
		
		if (opt_conflict_min_cut_maxflow) {
			
			Lit l = mkLit(outer->getEdgeVar(edgeid), false);
			//maintain the residual graph for conflict analysis
			lbool val = outer->value(l);
			if (val != l_False) {
				Weight flow = overapprox_conflict_detector->getEdgeFlow(edgeid);
				Weight capacity = overapprox_conflict_detector->getEdgeCapacity(edgeid);
				if (capacity == 0) {
					//learn_graph.disableEdge(e.id);
					//learn_graph.disableEdge(back_edges[e.id]);
					
					learn_graph.disableEdge(edgeid * 2);
					learn_graph.disableEdge(back_edges[edgeid * 2]);
					learn_graph.disableEdge(edgeid * 2 + 1);
					learn_graph.disableEdge(back_edges[edgeid * 2 + 1]);
				} else {
					if (flow > 0) {
						//then there is capacity in the backward edge in the residual graph
						/*					int back_edge = back_edges[e.id];
						 learn_graph.enableEdge(back_edges[e.id]);
						 learn_caps[back_edge] = INF;*/

						learn_graph.disableEdge(back_edges[edgeid * 2]);
						learn_graph.enableEdge(back_edges[edgeid * 2 + 1]);
					} else {
						learn_graph.disableEdge(back_edges[edgeid * 2]);
						learn_graph.disableEdge(back_edges[edgeid * 2 + 1]);
					}
					if (flow < capacity) {
						//then there is capacity in the forward edge in the residual graph
						/*	learn_caps[e.id] = INF;
						 learn_graph.enableEdge(e.id);
						 */
						learn_graph.disableEdge(edgeid * 2);
						//learn_graph.disableEdge(back_edges[edgeid*2]);
						learn_graph.enableEdge(edgeid * 2 + 1);
						//learn_graph.disableEdge(back_edges[edgeid*2+1]);
					} else {
						assert(flow == capacity);
						//learn_graph.disableEdge(e.id);
						
						learn_graph.disableEdge(edgeid * 2);
						//learn_graph.disableEdge(back_edges[edgeid*2]);
						learn_graph.disableEdge(edgeid * 2 + 1);
						//learn_graph.disableEdge(back_edges[edgeid*2+1]);
					}
				}
			} else {
				learn_graph.enableEdge(edgeid * 2);
				learn_graph.disableEdge(back_edges[edgeid * 2]);
				learn_graph.disableEdge(edgeid * 2 + 1);
				learn_graph.disableEdge(back_edges[edgeid * 2 + 1]);
			}
		}
	}
	
	overapprox_conflict_detector->clearChangedEdges();
	dbg_decisions();
}

template<typename Weight>
void MaxflowDetector<Weight>::dbg_decisions() {
#ifndef NDEBUG
	static int iter = 0;
	++iter;
	if (iter == 911) {
		int a = 1;
	}
	//printf("Potential Decisions %d: ",iter);
	for (int edgeID = 0; edgeID < g_under.edges(); edgeID++) {
		
		if (is_potential_decision[edgeID]) {
			/*if(potential_decisions.contains(edgeID) || potential_decisions_q.contains(edgeID)){
			 printf("%d ",edgeID);
			 }*/
			Lit l = mkLit(outer->getEdgeVar(edgeID));
			//assert(!decisions.contains(l));
			//assert(!decisions.contains(~l));
			bassert(potential_decisions.count(edgeID) <= 1);
			bassert(potential_decisions_q.count(edgeID) <= 1);
			bassert(
					potential_decisions.contains(edgeID) || potential_decisions_q.contains(edgeID)
							|| decisions.contains(l) || decisions.contains(~l));
		}
	}
	for (Lit l : decisions) {
		Var v = var(l);
		int edgeID = outer->getEdgeID(v);
		assert(is_potential_decision[edgeID]);
		bassert(!potential_decisions.contains(edgeID));
		bassert(!potential_decisions_q.contains(edgeID));
		//printf("%d ",edgeID);
	}
	/*printf("\n");
	 printf("Decisions %d: ",iter);

	 printf("\n");*/
#endif
}

template<typename Weight>
Lit MaxflowDetector<Weight>::decide(int level) {
	static int it = 0;
	if (++it == 59) {
		int a = 1;
	}
	double startdecidetime = rtime(2);
	auto * over = overapprox_conflict_detector;
	auto * under = underapprox_conflict_detector;
	
	if (opt_lazy_maxflow_decisions) {
		collectChangedEdges();
		
#ifndef NDEBUG
		{
			bool found = false;
			for (int edgeID = 0; edgeID < g_under.edges(); edgeID++) {
				if (edgeID == 62) {
					int a = 1;
				}
				if (g_under.edgeEnabled(edgeID)) {
					bool hasFlow = over->getEdgeFlow(edgeID) > 0;
					Lit l = mkLit(outer->getEdgeVar(edgeID));
					if (hasFlow) {
						assert(is_potential_decision[edgeID] || decisions.contains(l));
					}
				}
			}
		}
#endif
		
		Lit decision = lit_Undef;
		/*if(opt_maxflow_decisions_q==0){
		 //should probably look into ordering heuristics for which edge to decide first, here...
		 while(potential_decisions.size()>0){
		 int edgeID = potential_decisions.last();
		 assert(is_potential_decision[edgeID]);


		 potential_decisions.pop();
		 Lit l = mkLit(outer->getEdgeVar(edgeID),false);
		 if(outer->value(l)==l_Undef && over->getEdgeFlow(edgeID)>0){
		 decideEdge(edgeID,level,true);
		 decision = l;
		 break;
		 }else if (outer->value(l)==l_True){
		 if(over->getEdgeFlow(edgeID)>0)//this check is optional
		 decideEdge(edgeID,level,true);
		 else{
		 is_potential_decision[edgeID]=false;
		 }
		 }else{
		 assert(over->getEdgeFlow(edgeID)==0);
		 is_potential_decision[edgeID]=false;
		 //decideEdge(edgeID,level,false);
		 }
		 }
		 }else{*/
		//should probably look into ordering heuristics for which edge to decide first, here...
		while (potential_decisions_q.size() > 0) {
			int edgeID;
			
			if (opt_maxflow_decisions_q == 0) {
				edgeID = potential_decisions_q.peekBack();
				potential_decisions_q.popBack();
			} else {
				edgeID = potential_decisions_q.peek();
				potential_decisions_q.pop();
			}
			assert(is_potential_decision[edgeID]);
			Lit l = mkLit(outer->getEdgeVar(edgeID), false);
			if (outer->value(l) == l_Undef && over->getEdgeFlow(edgeID) > 0) {
				decideEdge(edgeID, level, true);
				decision = l;
				break;
			} else if (outer->value(l) == l_True) {
				if (over->getEdgeFlow(edgeID) > 0) {			//this check is optional
					decideEdge(edgeID, level, true);
				} else {
					is_potential_decision[edgeID] = false;
				}
			} else {
				assert(over->getEdgeFlow(edgeID) == 0);
				is_potential_decision[edgeID] = false;
			}
		}
		//}
#ifndef NDEBUG
		{
			bool found = false;
			for (int edgeID = 0; edgeID < g_under.edges(); edgeID++) {
				if (g_under.edgeEnabled(edgeID)) {
					bool hasFlow = over->getEdgeFlow(edgeID) > 0;
					Lit l = mkLit(outer->getEdgeVar(edgeID));
					if (hasFlow) {
						assert(is_potential_decision[edgeID] || decisions.contains(l));
					}
				}
			}
		}
#endif
		double post_time = rtime(2);
		stats_decide_time += post_time - startdecidetime;
		stats_redecide_time += post_time - startdecidetime;
		return decision;
	} else if (opt_old_lazy_maxflow_decisions) {
		if (last_decision_status != over->numUpdates()) {
			last_decision_status = over->numUpdates();
			q.clear();
			last_decision_q_pos = 0;
			seen.clear();
			seen.growTo(g_over.nodes(), false);
			
			if (opt_conflict_from_source)
				q.push_back(source);
			else
				q.push_back(target);
		}
		
		stats_decision_calculations++;
		
		double post_calc_time = rtime(2);
		
		if (opt_conflict_dfs) {
			//do a dfs to find that edge. Could start from either the source or the target.
			if (opt_conflict_from_source) {
				
				while (q.size()) {
					int u = q.back();
					q.pop_back();
					int qs = q.size();
					for (int i = 0; i < g_over.nIncident(u); i++) {
						if (!g_over.edgeEnabled(g_over.incident(u, i).id))
							continue;
						int v = g_over.incident(u, i).node;
						int edgeID = g_over.incident(u, i).id;
						if (over->getEdgeFlow(edgeID) > 0) {
							Var var = outer->getEdgeVar(edgeID);
							if (outer->value(var) == l_Undef) {
								q.resize(qs);
								q.push_back(u);
								return (mkLit(var, false));
							}
							if (!seen[v]) {
								seen[v] = true;
								q.push_back(v);
							}
							
						}
					}
				}
			} else {
				
				while (q.size()) {
					int u = q.back();
					int qs = q.size();
					q.pop_back();
					for (int i = 0; i < g_over.nIncoming(u); i++) {
						if (!g_over.edgeEnabled(g_over.incoming(u, i).id))
							continue;
						int v = g_over.incoming(u, i).node;
						int edgeID = g_over.incoming(u, i).id;
						if (over->getEdgeFlow(edgeID) > 0) {
							Var var = outer->getEdgeVar(edgeID);
							if (outer->value(var) == l_Undef) {
								q.resize(qs);
								q.push_back(u);
								return (mkLit(var, false));
							}
							if (!seen[v]) {
								seen[v] = true;
								q.push_back(v);
							}
						}
					}
				}
			}
		} else {
			if (opt_conflict_from_source) {
				
				//do a bfs to find that edge. Could start from either the source or the target.
				for (; last_decision_q_pos < q.size(); last_decision_q_pos++) {
					int u = q[last_decision_q_pos];
					for (int i = 0; i < g_over.nIncident(u); i++) {
						if (!g_over.edgeEnabled(g_over.incident(u, i).id))
							continue;
						int edgeID = g_over.incident(u, i).id;
						int v = g_over.incident(u, i).node;
						if (over->getEdgeFlow(edgeID) > 0) {
							Var var = outer->getEdgeVar(edgeID);
							assert(outer->value(var)!=l_False);
							if (outer->value(var) == l_Undef) {
								return (mkLit(var, false));
							}
							if (!seen[v]) {
								seen[v] = true;
								q.push_back(v);
							}
						}
					}
				}
			} else {
				
				//do a bfs to find that edge. Could start from either the source or the target.
				for (; last_decision_q_pos < q.size(); last_decision_q_pos++) {
					int u = q[last_decision_q_pos];
					for (int i = 0; i < g_over.nIncoming(u); i++) {
						if (!g_over.edgeEnabled(g_over.incoming(u, i).id))
							continue;
						int edgeID = g_over.incoming(u, i).id;
						int v = g_over.incoming(u, i).node;
						if (over->getEdgeFlow(edgeID) > 0) {
							Var var = outer->getEdgeVar(edgeID);
							assert(outer->value(var)!=l_False);
							if (outer->value(var) == l_Undef) {
								return (mkLit(var, false));
							}
							if (!seen[v]) {
								seen[v] = true;
								q.push_back(v);
							}
						}
					}
				}
			}
		}
		stats_flow_recalc_time += rtime(2) - post_calc_time;
		
	} else {
		
		//Weight under_flow = under->maxFlow(source,target) ;
		if (last_decision_status != over->numUpdates())
			to_decide.clear();
		
		if (to_decide.size()) {
			while (to_decide.size()) {
				Lit l = to_decide.last();
				to_decide.pop();
				if (outer->value(l) == l_Undef) {
					double post_time = rtime(2);
					stats_decide_time += post_time - startdecidetime;
					stats_redecide_time += post_time - startdecidetime;
					return l;
				}
			}
		}
		//Weight under_flow = positive_detector->maxFlow(source,target);//intentionally not using the conflict detector; it isn't required.
		Weight over_flow = overapprox_detector->maxFlow();
		
		for (int k = 0; k < flow_lits.size(); k++) {
			Lit l = flow_lits[k].l;
			if (l == lit_Undef)
				continue;
			
			Weight & required_flow = flow_lits[k].max_flow;
			
			if (outer->value(l) == l_True && opt_decide_graph_pos) {
				assert(over_flow >= required_flow);
				
#ifndef NDEBUG
				static vec<bool> dbg_expect;
				int dbg_count = 0;
				dbg_expect.clear();
				dbg_expect.growTo(g_under.edges());
				for (int edgeID = 0; edgeID < g_under.edges(); edgeID++) {
					lbool val = outer->value(outer->getEdgeVar(edgeID));
					if (val == l_Undef) {
						if (over->getEdgeFlow(edgeID) > 0) {
							dbg_expect[edgeID] = true;
							dbg_count++;
						}
					} else if (val == l_False) {
						//assert(over->getEdgeFlow(edgeID)==0);
					}
				}
#endif
				//in principle we can do this check, and it can avoid un-needed decisions - but if we are detecting pure theory literals,
				//then we might not have computed the underflow yet, and it is probably too expensive to compute here in that case.
				
				//if(under_flow <required_flow)
				{
					stats_decision_calculations++;
					//then decide an unassigned edge of the currently selected flow
					to_decide.clear();
					over->maxFlow();
					double st = rtime(2);
					over->getEdgeFlow(0);
					double post_calc_time = rtime(2);
					stats_flow_calc_time += post_calc_time - st;
					last_decision_status = over->numUpdates();
					
					seen.clear();
					seen.growTo(g_over.nodes(), false);
					q.clear();
					
					if (opt_conflict_dfs) {
						//do a dfs to find that edge. Could start from either the source or the target.
						if (opt_conflict_from_source) {
							
							q.push_back(source);
							
							while (q.size()) {
								int u = q.back();
								q.pop_back();
								for (int i = 0; i < g_over.nIncident(u); i++) {
									if (!g_over.edgeEnabled(g_over.incident(u, i).id))
										continue;
									int v = g_over.incident(u, i).node;
									int edgeID = g_over.incident(u, i).id;
									if (over->getEdgeFlow(edgeID) > 0) {
										Var var = outer->getEdgeVar(edgeID);
										if (outer->value(var) == l_Undef) {
											to_decide.push(mkLit(var, false));
										}
										if (!seen[v]) {
											seen[v] = true;
											q.push_back(v);
										}
									}
								}
							}
						} else {
							q.push_back(target);
							while (q.size()) {
								int u = q.back();
								q.pop_back();
								for (int i = 0; i < g_over.nIncoming(u); i++) {
									if (!g_over.edgeEnabled(g_over.incoming(u, i).id))
										continue;
									int v = g_over.incoming(u, i).node;
									int edgeID = g_over.incoming(u, i).id;
									if (over->getEdgeFlow(edgeID) > 0) {
										Var var = outer->getEdgeVar(edgeID);
										if (outer->value(var) == l_Undef) {
											to_decide.push(mkLit(var, false));
										}
										if (!seen[v]) {
											seen[v] = true;
											q.push_back(v);
										}
									}
								}
							}
						}
					} else {
						if (opt_conflict_from_source) {
							q.push_back(source);
							//do a bfs to find that edge. Could start from either the source or the target.
							for (int i = 0; i < q.size(); i++) {
								int u = q[i];
								for (int i = 0; i < g_over.nIncident(u); i++) {
									if (!g_over.edgeEnabled(g_over.incident(u, i).id))
										continue;
									int edgeID = g_over.incident(u, i).id;
									int v = g_over.incident(u, i).node;
									if (over->getEdgeFlow(edgeID) > 0) {
										Var var = outer->getEdgeVar(edgeID);
										assert(outer->value(var)!=l_False);
										if (outer->value(var) == l_Undef) {
											to_decide.push(mkLit(var, false));
										}
										if (!seen[v]) {
											seen[v] = true;
											q.push_back(v);
										}
									}
								}
							}
						} else {
							q.push_back(target);
							//do a bfs to find that edge. Could start from either the source or the target.
							for (int i = 0; i < q.size(); i++) {
								int u = q[i];
								for (int i = 0; i < g_over.nIncoming(u); i++) {
									if (!g_over.edgeEnabled(g_over.incoming(u, i).id))
										continue;
									int edgeID = g_over.incoming(u, i).id;
									int v = g_over.incoming(u, i).node;
									if (over->getEdgeFlow(edgeID) > 0) {
										Var var = outer->getEdgeVar(edgeID);
										assert(outer->value(var)!=l_False);
										if (outer->value(var) == l_Undef) {
											to_decide.push(mkLit(var, false));
										}
										if (!seen[v]) {
											seen[v] = true;
											q.push_back(v);
										}
									}
								}
							}
						}
					}
					stats_flow_recalc_time += rtime(2) - post_calc_time;
					//assert(to_decide.size()==dbg_count);
				}
			} else if (outer->value(l) == l_False && opt_decide_graph_neg) {
				
			}
			if (to_decide.size() && last_decision_status == over->numUpdates()) {
				while (to_decide.size()) {
					Lit l = to_decide.last();
					to_decide.pop();
					if (outer->value(l) == l_Undef) {
						stats_decide_time += rtime(2) - startdecidetime;
						return l;
					}
				}
			}
		}
	}
	stats_decide_time += rtime(2) - startdecidetime;
	return lit_Undef;
}
;

template class MaxflowDetector<int> ;
template class MaxflowDetector<long> ;
template class MaxflowDetector<double> ;
#include <gmpxx.h>
template class MaxflowDetector<mpq_class> ;
