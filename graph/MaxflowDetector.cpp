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

#include <core/Config.h>
#include <dgl/Dinics.h>
#include <dgl/DinicsLinkCut.h>
#include <dgl/EdmondsKarpAdj.h>
#include <dgl/EdmondsKarpDynamic.h>
#include <dgl/Reach.h>
#include <dgl/BFS.h>
#include <graph/GraphTheory.h>
#include <graph/MaxflowDetector.h>
#include <mtl/Rnd.h>
#include <mtl/Vec.h>
#include <utils/Options.h>
//#include "dgl/KohliTorr.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
using namespace Monosat;

template<typename Weight>
MaxflowDetector<Weight>::MaxflowDetector(int _detectorID, GraphTheorySolver<Weight> * _outer,
		 DynamicGraph<Weight>  &_g, DynamicGraph<Weight>  &_antig, int from, int _target, double seed) :
		Detector(_detectorID), outer(_outer),  g_under(_g), g_over(_antig), source(
				from), target(_target), rnd_seed(seed),order_heap(EdgeOrderLt(activity)) {
	var_decay =opt_var_decay;

	MinCutAlg alg = mincutalg;
	if(outer->hasBitVectorEdges()){
		if (alg!= MinCutAlg::ALG_EDKARP_ADJ && alg != MinCutAlg::ALG_KOHLI_TORR){
			printf("Note: falling back on kohli-torr for maxflow, because edge weights are bitvectors\n");
			alg=MinCutAlg::ALG_KOHLI_TORR;
		}
	}

	if (alg == MinCutAlg::ALG_EDKARP_DYN) {
		underapprox_detector = new EdmondsKarpDynamic<Weight>(_g, source, target);
		overapprox_detector = new EdmondsKarpDynamic<Weight>(_antig, source, target);
		underapprox_conflict_detector = underapprox_detector; //new EdmondsKarpAdj<std::vector<Weight>,Weight>(_g,capacities);
		overapprox_conflict_detector = overapprox_detector; // new EdmondsKarpAdj<std::vector<Weight>,Weight>(_antig,capacities);
		if (opt_conflict_min_cut_maxflow || opt_adaptive_conflict_mincut)
			learn_cut = new EdmondsKarpDynamic<Weight>(learn_graph, source, target);

		/*if(opt_conflict_min_cut_maxflow)
		 learn_cut = new EdmondsKarpAdj<std::vector<int>,int>(learn_graph,learn_caps,source,target);*/
	} else if (alg == MinCutAlg::ALG_EDKARP_ADJ) {
		underapprox_detector = new EdmondsKarpAdj<Weight>(_g, source, target);
		overapprox_detector = new EdmondsKarpAdj<Weight>(_antig, source, target);
		underapprox_conflict_detector = underapprox_detector;
		overapprox_conflict_detector = overapprox_detector;
		if (opt_conflict_min_cut_maxflow || opt_adaptive_conflict_mincut)
			learn_cut = new EdmondsKarpAdj<Weight>(learn_graph, source, target);
	} else if (alg == MinCutAlg::ALG_DINITZ) {
		underapprox_detector = new Dinitz<Weight>(_g, source, target);
		overapprox_detector = new Dinitz<Weight>(_antig, source, target);
		underapprox_conflict_detector = underapprox_detector; // new EdmondsKarpAdj<std::vector<Weight>,Weight>(_g,capacities);
		overapprox_conflict_detector = overapprox_detector; //new EdmondsKarpAdj<std::vector<Weight>,Weight>(_antig,capacities);
		if (opt_conflict_min_cut_maxflow || opt_adaptive_conflict_mincut)
			learn_cut = new Dinitz<Weight>(learn_graph, source, target);
	} else if (alg == MinCutAlg::ALG_DINITZ_LINKCUT) {
		//link-cut tree currently only supports ints (enforcing this using tempalte specialization...).
		underapprox_detector = new DinitzLinkCut<Weight>(g_under, source, target);
		overapprox_detector = new DinitzLinkCut<Weight>(g_over, source, target);
		underapprox_conflict_detector = underapprox_detector;
			overapprox_conflict_detector = overapprox_detector;

		if (opt_conflict_min_cut_maxflow || opt_adaptive_conflict_mincut)
			learn_cut = new EdmondsKarpAdj<Weight>(learn_graph, source, target);
	} else if (alg == MinCutAlg::ALG_KOHLI_TORR) {
		underapprox_detector = new KohliTorr<Weight>(_g, source, target,
				opt_kt_preserve_order);
		overapprox_detector = new KohliTorr<Weight>(_antig, source, target,
				opt_kt_preserve_order);
		if (opt_use_kt_for_conflicts) {
			underapprox_conflict_detector = underapprox_detector; //new EdmondsKarpDynamic<std::vector<Weight>,Weight>(_g,capacities);
			overapprox_conflict_detector = overapprox_detector; //new EdmondsKarpDynamic<std::vector<Weight>,Weight>(_antig,capacities);
		} else {
			//for reasons I don't yet understand, kohli-torr seems to produce maxflows that work very poorly as theory-decisions for some problems.
			//Possibly, this is because Kohli-torr will sometimes produces flows that use multiple paths even when a maxflow with a single path is possible.
			underapprox_conflict_detector = new EdmondsKarpDynamic<Weight>(_g, source,
					target);
			overapprox_conflict_detector = new EdmondsKarpDynamic<Weight>(_antig,
					source, target);
		}
		if (opt_conflict_min_cut_maxflow || opt_adaptive_conflict_mincut)
			learn_cut = new KohliTorr<Weight>(learn_graph, source, target, opt_kt_preserve_order);
		
	} else {
		underapprox_detector = new EdmondsKarpAdj<Weight>(_g, source, target);
		overapprox_detector = new EdmondsKarpAdj<Weight>(_antig, source, target);
		underapprox_conflict_detector = underapprox_detector;
		overapprox_conflict_detector = overapprox_detector;
		if (opt_conflict_min_cut_maxflow || opt_adaptive_conflict_mincut)
			learn_cut = new EdmondsKarpAdj<Weight>(learn_graph, source, target);
	}
	if (opt_adaptive_history_clear > 0) {
		learn_graph.adaptive_history_clear = true;
		learn_graph.historyClearInterval = opt_adaptive_history_clear;
	} else {

		learn_graph.historyClearInterval = opt_history_clear;
	}
	learn_graph.dynamic_history_clears=opt_dynamic_history_clear;
	learn_graph.disable_history_clears=disable_history_clears;


	if(learn_cut){
		alg_id = g_over.addDynamicAlgorithm(this);
	}

	if (opt_learn_acyclic_flows){
		acyclic_flow=new AcyclicFlow<Weight>(g_under);
	}

	if(opt_record){
		std::string t = (const char*)opt_record_file;
		t+="/LOG_LEARN_GRAPH" +std::to_string(getID());
		learn_graph.outfile = fopen(t.c_str(), "w");
	}

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
void MaxflowDetector<Weight>::addFlowBVLessThan(const BitVector<Weight>  &bv, Var outer_reach_var, bool inclusive) {
	g_under.invalidate();
	g_over.invalidate();
	Var reach_var = outer->newVar(outer_reach_var, getID());
	if (first_reach_var == var_Undef) {
		first_reach_var = reach_var;
	} else {
		assert(reach_var > first_reach_var);
	}

	Lit reachLit = mkLit(reach_var, false);
	bool found = false;

	flow_lits.push();
	flow_lits.last().l = reachLit;
	flow_lits.last().bv = bv;
	flow_lits.last().max_flow = -1;
	flow_lits.last().inclusive=inclusive;
	while (reach_lit_map.size() <= reach_var - first_reach_var) {
		reach_lit_map.push(-1);
	}

	reach_lit_map[reach_var - first_reach_var] = flow_lits.size() - 1;

}


template<typename Weight>
Lit MaxflowDetector<Weight>::findFirstReasonTooHigh(Weight flow) {

	Weight actual_flow = underapprox_conflict_detector->maxFlow();

	//just collect the set of edges which have non-zero flow, and return them
	//(Or, I could return a cut, probably)
	for (int i = 0; i < g_under.edges(); i++) {
		if (g_under.edgeEnabled(i)) {
			if (underapprox_conflict_detector->getEdgeFlow(i) > 0) {
				Var v = outer->edge_list[i].v;
				Lit l = mkLit(v,true);
				assert(outer->value(v)==l_True);
				if(!g_under.isConstant(i) && outer->decidable(l)){
					return l;
				}
				if(outer->hasBitVector(i)){
					throw std::runtime_error("Error in maxflow");
					//outer->buildBVReason(outer->getEdgeBV(i).getID(),Comparison::geq,underapprox_conflict_detector->getEdgeFlow(i),conflict);
					//outer->buildBVReason(bv.getID(),inclusive ? Comparison::gt:Comparison::geq,bv.getUnder(),conflict);
				}
			}
		}
	}
	return lit_Undef;
}
template<typename Weight>
Lit MaxflowDetector<Weight>::findFirstReasonTooLow(Weight flow) {
	Weight foundflow = overapprox_conflict_detector->maxFlow();
	assert(!seen.contains(true));

	visit.clear();
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
					}else{
						//if the edge _was_ enabled, and all of its capacity was used, then the reason that it didn't have more capacity must be included.
						if(outer->hasBitVector(edgeid)){
							throw std::runtime_error("Error in maxflow");
						}
					}
				} else {
					//this is a disabled edge, and we can add it to the cut.
					//we're going to assume the edge has non-zero capacity here, otherwise we could exclude it (but it shouldn't really even be in this graph in that case, anyways).
					if( !g_over.isConstant(edgeid)){
						Lit l = mkLit(v);
						if(outer->decidable(l)){
							while(visit.size()){
								seen[visit.last()]=false;
								visit.pop();
							}
							visit.clear();

							return l;
						}
					}
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
					if (residual_capacity>0) {
						seen[p] = true;
						visit.push(p);
					}else if (g_over.getWeight(edgeid)==0){
						//if the edge _was_ enabled, and it had no capacity capacity was used, then the reason that it didn't have more capacity must be included.
						//does this really hold for backwards edges?
						if(outer->hasBitVector(edgeid)){
							throw std::runtime_error("Error in maxflow");
							//outer->buildBVReason(bv.getID(),inclusive ? Comparison::gt:Comparison::geq,bv.getUnder(),conflict);
						}
					}
				} else {
					//this is a disabled edge, and we can add it to the cut.
					//we're going to assume the edge has non-zero capacity here, otherwise we could exclude it (but it shouldn't really even be in this graph in that case, anyways).
					if( !g_over.isConstant(edgeid)){
						Lit l = mkLit(v);
						if(outer->decidable(l)){
							while(visit.size()){
								seen[visit.last()]=false;
								visit.pop();
							}
							visit.clear();

							return l;
						}
					}
				}
			}
			//pred
		}
	}
	while(visit.size()){
		seen[visit.last()]=false;
		visit.pop();
	}
	return lit_Undef;
}


template<typename Weight>
void MaxflowDetector<Weight>::buildMaxFlowTooHighReason(Weight flow, vec<Lit> & conflict) {
	//drawFull();
	double starttime = rtime(2);
	tmp_cut.clear();
	Weight actual_flow = underapprox_conflict_detector->maxFlow();
	if(opt_learn_acyclic_flows){

		refined_flow.resize(g_under.edges());
		for(int i = 0;i<g_under.edges();i++){
			if(g_under.hasEdge(i) && g_under.edgeEnabled(i)){
				refined_flow[i]=underapprox_conflict_detector->getEdgeFlow(i);
			}else{
				refined_flow[i]=0;
			}
		}

		acyclic_flow->getAcyclicFlow(source,target,refined_flow);
		for (int i = 0; i < g_under.edges(); i++) {
			if (g_under.edgeEnabled(i)) {
				if (refined_flow[i] > 0) {
					Var v = outer->edge_list[i].v;
					assert(outer->value(v)==l_True);
					if(!g_under.isConstant(i)){
						conflict.push(mkLit(v, true));
					}
					if(outer->hasBitVector(i) && ! outer->getEdgeBV(i).isConst()){
						outer->buildBVReason(outer->getEdgeBV(i).getID(),Comparison::geq,refined_flow[i],conflict);
					}
				}

			}
		}
	}else{
		//just collect the set of edges which have non-zero flow, and return them
		for (int i = 0; i < g_under.edges(); i++) {
			if (g_under.edgeEnabled(i)) {
				if (underapprox_conflict_detector->getEdgeFlow(i) > 0) {
					Var v = outer->edge_list[i].v;
					assert(outer->value(v)==l_True);
					if(!g_under.isConstant(i)){
						conflict.push(mkLit(v, true));
					}
					if(outer->hasBitVector(i)){
						outer->buildBVReason(outer->getEdgeBV(i).getID(),Comparison::geq,underapprox_conflict_detector->getEdgeFlow(i),conflict);
						//outer->buildBVReason(bv.getID(),inclusive ? Comparison::gt:Comparison::geq,bv.getUnder(),conflict);
					}
				}
	
			}
		}
	}
	bumpConflictEdges(conflict);


	if (g_under.outfile) {
		std::sort(conflict.begin(), conflict.end());
		fprintf(g_under.outfile, "toohigh ");
		for (int i = 0; i < conflict.size(); i++) {
			fprintf(g_under.outfile, "%d,", dimacs(conflict[i]));
		}
		fprintf(g_under.outfile, "\n");
		fflush(g_under.outfile);
	}
	if (g_over.outfile) {
		std::sort(conflict.begin(), conflict.end());
		fprintf(g_over.outfile, "toohigh ");
		for (int i = 0; i < conflict.size(); i++) {
			fprintf(g_over.outfile, "%d,", dimacs(conflict[i]));
		}
		fprintf(g_over.outfile, "\n");
		fflush(g_over.outfile);
	}

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
		throw std::runtime_error("Assertion error");
	}
#endif
}

template<typename Weight>
void MaxflowDetector<Weight>::buildMaxFlowTooLowReason(Weight maxflow, vec<Lit> & conflict, bool force_maxflow) {
	static int it = 0;
	++it;
	if (it == 3) {
		int a = 1;
	}
	if(opt_verb>1){
		printf("Maxflow conflict %d, graph %d\n", it, outer->getTheoryIndex());
	}
	if(g_over.edges()==0)
		return;
	//printf("%d\n",it);
	double starttime = rtime(2);
	if (force_maxflow || opt_conflict_min_cut_maxflow) {
		Weight foundflow = overapprox_conflict_detector->maxFlow();
		collectChangedEdges();
		collectDisabledEdges();
		//g_over.drawFull(true);
		//learn_graph.drawFull(true);
#ifndef NDEBUG
		for (auto & e : g_over.getEdges()) {
			int from = e.from;
			int to = e.to;
			int v = outer->getEdgeVar(e.id);
			lbool val = outer->value(v);
			int edgeid = e.id;
			if (g_over.hasEdge(edgeid) && g_over.edgeEnabled(edgeid) ){

				Weight flow = overapprox_conflict_detector->getEdgeFlow(edgeid);
				Weight capacity = overapprox_conflict_detector->getEdgeCapacity(edgeid);
				assert(capacity<=g_over.getWeight(edgeid));
				Weight level0capacity =  outer->hasBitVector(edgeid)?	 outer->getEdgeBV(edgeid).getOver(true):capacity;

				assert(learn_graph.getWeight(back_edges[edgeid])==0x0FF0F0);//capacity of this edge in the backward direction is the same as the forward flow.
				assert(learn_graph.edgeEnabled(back_edges[edgeid]) == (flow>0));
				if (flow<capacity){
					//then this edge is not fully utilized, so it cannot be a limiting factor in the maxflow. Set its weight to infinity.
					assert(learn_graph.getWeight(edgeid)==0x0FF0F0);
				}else{
					assert(learn_graph.edgeEnabled(edgeid) == ((level0capacity-capacity)>0));
					assert(learn_graph.getWeight(edgeid)==(((level0capacity-capacity) > 0) ? 1:0));
				}
			}else{
				//flow is 0.
				assert(!learn_graph.edgeEnabled(back_edges[edgeid]));
				assert(learn_graph.getWeight(back_edges[edgeid])==0);//capacity of this edge in the backward direction is the same as the forward flow, which is 0.
				Weight w = (outer->hasBitVector(edgeid)?	 outer->getEdgeBV(edgeid).getOver(true):outer->edge_weights[edgeid]);
				assert(learn_graph.getWeight(edgeid)==(w>0?1:0));
			}

		}
#endif
		//g_over.drawFull(true);
		//learn_graph.drawFull(true);

		Weight f = learn_cut->minCut(cut);
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

		if (f < 0x0FF0F0) {
			assert(f < 0x0FF0F0);

			for (int i = 0; i < cut.size(); i++) {
				MaxFlowEdge e = cut[i];
				int edgeID = e.id;//because we've doubled the set of edges in the learn graph relative to g_over/g_under.
				Lit l = mkLit(outer->getEdgeVar(edgeID), false);
				if(outer->value(l)==l_False){//it is possible for the edge to be enabled, but to be set to capacity 0.
					bassert(outer->value(l) == l_False);
					conflict.push(l);

				}else if(outer->hasBitVector(edgeID) && ! outer->getEdgeBV(edgeID).isConst()){
					Weight residual = overapprox_conflict_detector->getEdgeResidualCapacity(edgeID);
					assert(overapprox_conflict_detector->getEdgeResidualCapacity(edgeID)==0);//the edge has no residual capacity, so its capacity must be increased.
					outer->buildBVReason(outer->getEdgeBV(edgeID).getID(),Comparison::leq,g_over.getWeight(edgeID),conflict);

				}
			}
		} else {
			int a = 1;

			//there is no way to increase the max flow.
		}

	/*	printf("conflict in theory %d ", outer->getTheoryIndex());
		if (f < 0x0FF0F0) {
			assert(f < 0x0FF0F0);

			for (int i = 0; i < cut.size(); i++) {
				MaxFlowEdge e = cut[i];
				int edgeID = e.id;//because we've doubled the set of edges in the learn graph relative to g_over/g_under.
				Lit l = mkLit(outer->getEdgeVar(edgeID), false);
				if(outer->value(l)==l_False){//it is possible for the edge to be enabled, but to be set to capacity 0.

					std::cout << "edge " << edgeID << " enabled ";
				}else if(outer->hasBitVector(edgeID)){
					std::cout << "edge " << edgeID << " bvid " << outer->getEdgeBV(edgeID).getID() << " > " << g_over.getWeight(edgeID)<< " ";
				}
			}
		}
		printf("\n");*/

		bumpConflictEdges(conflict);
		outer->num_learnt_cuts++;
		outer->learnt_cut_clause_length += (conflict.size() - 1);
		stats_over_conflicts++;
		double elapsed = rtime(2) - starttime;
		stats_over_conflict_time += elapsed;

		return;
	}
	
	//drawFull( non_reach_detectors[detector]->getSource(),u);
	//assert(outer->dbg_distance( source,u));
	//g_over.drawFull(true);
	//The reason why we can't reach this assignment is a cut through the disabled edges in the residual graph from the overapprox.
	//we could search for a min cut, but instead we will just step back in the RESIDUAL graph, from u to s, collecting disabled edges.
	assert(!seen.contains(true));

	visit.clear();
	Weight foundflow = overapprox_conflict_detector->maxFlow();
	/*			std::vector<MaxFlowEdge> ignore;
	 negative_conflict_detector->minCut(ignore);*/
	visit.push(target);
	seen[target]=true;
	for (int k = 0; k < visit.size(); k++) {
		int u = visit[k];
		for (int i = 0; i < g_over.nIncoming(u); i++) {
			int p = g_over.incoming(u, i).node;
			if(p==u)
				continue;//skip self loop edges, as they cannot be required for the maxflow
			if (!seen[p]) {
				
				int edgeid = g_over.incoming(u, i).id;


				int v = outer->getEdgeVar(edgeid);
				
				//assert( g.incoming(u,i).to==u);
				if (outer->value(v) != l_False) {
					//this is an enabled edge in the overapprox
					assert(overapprox_conflict_detector->getEdgeCapacity(edgeid)== g_over.getWeight(edgeid));
					
					Weight residual_capacity = overapprox_conflict_detector->getEdgeResidualCapacity(edgeid);
					if (residual_capacity > 0) {
						seen[p] = true;
						visit.push(p);
					}else{
						//if the edge _was_ enabled, and all of its capacity was used, then the reason that it didn't have more capacity must be included.
						if(outer->hasBitVector(edgeid) && ! outer->getEdgeBV(edgeid).isConst()){
							assert(g_over.getWeight(edgeid)==outer->getEdgeBV(edgeid).getOver());
							outer->buildBVReason(outer->getEdgeBV(edgeid).getID(),Comparison::leq,overapprox_conflict_detector->getEdgeFlow(edgeid),conflict);
							//outer->buildBVReason(bv.getID(),inclusive ? Comparison::gt:Comparison::geq,bv.getUnder(),conflict);
						}
					}
				} else {
					//this is a disabled edge, and we can add it to the cut.
					//we're going to assume the edge has non-zero capacity here, otherwise we could exclude it (but it shouldn't really even be in this graph in that case, anyways).
					if( !g_over.isConstant(edgeid)){
						//if(g_over.getWeight(edgeid)>0){
							conflict.push(mkLit(v, false));
				/*		}else if (!outer->constantWeight(edgeid)){
							if(outer->hasBitVector(edgeid)){
								outer->buildBVReason(outer->getEdgeBV(edgeid).getID(),Comparison::leq,0,conflict);
							}
						}*/
					}
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
					if (residual_capacity>0) {
						seen[p] = true;
						visit.push(p);
					}else if (g_over.getWeight(edgeid)==0){
						//if the edge _was_ enabled, and it had no capacity capacity was used, then the reason that it didn't have more capacity must be included.
						//does this really hold for backwards edges?
						if(outer->hasBitVector(edgeid)){
							outer->buildBVReason(outer->getEdgeBV(edgeid).getID(),Comparison::leq,0,conflict);
							//outer->buildBVReason(bv.getID(),inclusive ? Comparison::gt:Comparison::geq,bv.getUnder(),conflict);
						}
					}
				} else {
					//this is a disabled edge, and we can add it to the cut.
					//we're going to assume the edge has non-zero capacity here, otherwise we could exclude it (but it shouldn't really even be in this graph in that case, anyways).
					if( !g_over.isConstant(edgeid)){
						conflict.push(mkLit(v, false));
					}
				}
			}
			//pred
		}
	}
	
	while(visit.size()){
		seen[visit.last()]=false;
		visit.pop();
	}

	if (!force_maxflow && opt_adaptive_conflict_mincut > 0 && (conflict.size() - 1 > opt_adaptive_conflict_mincut)) { //-1 to ignore the predicate's literal stored at position 0
		double elapsed = rtime(2) - starttime;
		stats_over_conflict_time += elapsed;
		conflict.shrink(conflict.size() - 1);
		assert(conflict.size() == 1);
		buildMaxFlowTooLowReason(maxflow, conflict, true);
		return;
	}
	bumpConflictEdges(conflict);
	/*if(conflict.size()<dbg_minconflict()){
	 exit(4);
	 }*/
#ifndef NDEBUG
	if (overapprox_detector != overapprox_conflict_detector) {
		Weight foundflow2 = overapprox_detector->maxFlow();
		assert(foundflow == foundflow2);
	}
#endif



	if (g_under.outfile) {
		std::sort(conflict.begin(), conflict.end());
		fprintf(g_under.outfile, "toolow ");
		for (int i = 0; i < conflict.size(); i++) {
			fprintf(g_under.outfile, "%d,", dimacs(conflict[i]));
		}
		fprintf(g_under.outfile, "\n");
		fflush(g_under.outfile);
	}
	if (g_over.outfile) {
		std::sort(conflict.begin(), conflict.end());
		fprintf(g_over.outfile, "toolow ");
		for (int i = 0; i < conflict.size(); i++) {
			fprintf(g_over.outfile, "%d,", dimacs(conflict[i]));
		}
		fprintf(g_over.outfile, "\n");
		fflush(g_over.outfile);
	}

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
		DistLit & f = flow_lits[reach_lit_map[v - first_reach_var]];
		if (f.max_flow>=0){
			Weight flow = f.max_flow;
			//int u =getNode(v);
			buildMaxFlowTooHighReason(flow, reason);
		}else{
			auto & bv= f.bv;
			outer->buildBVReason(bv.getID(),f.inclusive ? Comparison::gt:Comparison::geq,bv.getUnder(),reason);
			buildMaxFlowTooHighReason(bv.getOver(), reason);
		}
		//double elapsed = rtime(2)-startpathtime;
		//	pathtime+=elapsed;
	} else if (marker == overprop_marker) {
		reason.push(p);
		
		//the reason is a cut separating p from s;
		//We want to find a min-cut in the full graph separating, where activated edges (ie, those still in antig) are weighted infinity, and all others are weighted 1.
		
		//This is a cut that describes a minimal set of edges which are disabled in the current graph, at least one of which would need to be activated in order for s to reach p
		//assign the mincut edge weights if they aren't already assigned.
		
		Var v = var(p);
		DistLit & f = flow_lits[reach_lit_map[v - first_reach_var]];
		if (f.max_flow>=0){
			Weight flow = f.max_flow;
			//int t = getNode(v); // v- var(reach_lits[d][0]);
			buildMaxFlowTooLowReason(flow, reason);
		}else{
			auto & bv = f.bv;
			outer->buildBVReason(bv.getID(),f.inclusive ? Comparison::gt:Comparison::geq,bv.getUnder(),reason);
			//int t = getNode(v); // v- var(reach_lits[d][0]);
			buildMaxFlowTooLowReason(bv.getUnder(), reason);
		}
		
	} else {
		assert(false);
	}
}
template<typename Weight>
bool MaxflowDetector<Weight>::propagate(vec<Lit> & conflict, bool backtrackOnly, Lit & conflictLit) {
	if (flow_lits.size() == 0) {
		return true;
	}
	static int iter1 = 0;
	if (++iter1 == 79) {
		int a = 1;
	}

	if (g_under.outfile) {
		fprintf(g_under.outfile, "iter %d\n", iter1);
		fflush(g_under.outfile);
	}
	if (g_over.outfile) {
		fprintf(g_over.outfile, "iter %d\n", iter1);
		fflush(g_over.outfile);
	}

	
	double start_prop_time = rtime(2);

	Weight over_maxflow = -1;
	Weight under_maxflow = -1;
	bool computed_under=false;
	bool computed_over=false;
	
	if (underapprox_detector && (!opt_detect_pure_theory_lits || unassigned_positives > 0)) {
		double startdreachtime = rtime(2);
		stats_under_updates++;
		computed_under=true;
		under_maxflow = underapprox_detector->maxFlow();
		assert(under_maxflow == underapprox_conflict_detector->maxFlow());
		double reachUpdateElapsed = rtime(2) - startdreachtime;
		stats_under_update_time += reachUpdateElapsed;
	} else
		stats_skipped_under_updates++;
	
	if (overapprox_detector && (!opt_detect_pure_theory_lits || unassigned_negatives > 0)) {
		double startunreachtime = rtime(2);
		stats_over_updates++;
		computed_over=true;
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
		if(f.max_flow>=0){
			Weight& maxflow = f.max_flow;
			//int u = getNode(var(l));
		
			if (computed_under && ((inclusive && under_maxflow >= maxflow) || (!inclusive && under_maxflow > maxflow))) {
				if (outer->value(l) == l_True) {
					//do nothing
				} else if (outer->value(l) == l_Undef) {
					//trail.push(Assignment(false,true,detectorID,0,var(l)));
					outer->enqueue(l, underprop_marker);

				} else if (outer->value(l) == l_False) {
					stats_total_prop_time += rtime(2)-start_prop_time;
					if(backtrackOnly){
						return false;
						/*conflictLit = findFirstReasonTooHigh(maxflow);
						if(conflictLit!=lit_Undef){
							return false;
						}*/
					}
					conflict.push(l);
					buildMaxFlowTooHighReason(maxflow, conflict);

					return false;
				}
				
			} else if (computed_over && ((inclusive && over_maxflow < maxflow) || (!inclusive && over_maxflow <= maxflow))) {
				if (outer->value(l) == l_False) {
					//do nothing
				} else if (outer->value(l) == l_Undef) {
					//trail.push(Assignment(false,false,detectorID,0,var(l)));
					outer->enqueue(~l, underprop_marker);

				} else if (outer->value(l) == l_True) {
					stats_total_prop_time += rtime(2)-start_prop_time;
					if(backtrackOnly){
						return false;
					/*	conflictLit= findFirstReasonTooLow(maxflow);
						if(conflictLit!=lit_Undef){
							return false;
						}*/
					}

					conflict.push(~l);
					buildMaxFlowTooLowReason(maxflow, conflict);
					return false;
				}
				
			}
		}else{
			Lit l = f.l;

			bool inclusive = f.inclusive;
			BitVector<Weight> & bv = f.bv;
			Weight & max_flow_under = bv.getUnder();
			Weight & max_flow_over = bv.getOver();

			if  (computed_under && ((inclusive && under_maxflow >= max_flow_over) || (!inclusive && under_maxflow > max_flow_over))) {
				if (outer->value(l) == l_True) {
					//do nothing
				} else if (outer->value(l) == l_Undef) {
					//trail.push(Assignment(false,true,detectorID,0,var(l)));
					outer->enqueue(l, underprop_marker);
					//should also enqueue that the flow is >= under->flow, and <= over->flow...
				} else if (outer->value(l) == l_False) {
					stats_total_prop_time += rtime(2)-start_prop_time;
					if(backtrackOnly)
						return false;

					conflict.push(l);
					outer->buildBVReason(bv.getID(),inclusive ? Comparison::leq:Comparison::lt,under_maxflow,conflict);
					buildMaxFlowTooHighReason(bv.getOver(), conflict);
					return false;
				}

			} else if  (computed_over && ((inclusive && over_maxflow < max_flow_under) || (!inclusive && over_maxflow <= max_flow_under))) {
				if (outer->value(l) == l_False) {
					//do nothing
				} else if (outer->value(l) == l_Undef) {
					//trail.push(Assignment(false,false,detectorID,0,var(l)));
					outer->enqueue(~l, underprop_marker);
					//should also enqueue that the flow is >= under->flow, and <= over->flow...
				} else if (outer->value(l) == l_True) {
					stats_total_prop_time += rtime(2)-start_prop_time;
					if(backtrackOnly)
						return false;

					conflict.push(~l);
					outer->buildBVReason(bv.getID(),inclusive ? Comparison::gt:Comparison::geq,over_maxflow,conflict);
					buildMaxFlowTooLowReason(bv.getUnder(), conflict);
					return false;
				}

			}
		}
	}
	stats_total_prop_time += rtime(2)-start_prop_time;
	return true;
}
template<typename Weight>
bool MaxflowDetector<Weight>::checkSatisfied() {
	g_under.drawFull(true);
	EdmondsKarpAdj<Weight> underCheck(g_under, source, target);
	EdmondsKarpAdj<Weight> overCheck(g_over, source, target);
	for (int j = 0; j < flow_lits.size(); j++) {
		Lit l = flow_lits[j].l;
		if(flow_lits[j].max_flow>=0){
			Weight dist = flow_lits[j].max_flow;
			
			if (l != lit_Undef) {
				//int node =getNode(var(l));

				if (outer->value(l) == l_True) {
					if (underCheck.maxFlow() < dist) {
						return false;
					}
				} else if (outer->value(l) == l_False) {
					if (overCheck.maxFlow() >= dist) {
						return false;
					}
				} else {
					if (underCheck.maxFlow() >= dist) {
						return false;
					}
					if (overCheck.maxFlow() < dist) {
						return false;
					}
				}
			}
		}else{
			Weight dist_under = flow_lits[j].bv.getUnder();
			Weight dist_over = flow_lits[j].bv.getOver();

			if (l != lit_Undef) {
				//int node =getNode(var(l));

				if (outer->value(l) == l_True) {
					if (underCheck.maxFlow() < dist_under) {
						return false;
					}
				} else if (outer->value(l) == l_False) {
					if (overCheck.maxFlow() >= dist_over) {
						return false;
					}
				} else {
					if (underCheck.maxFlow() >= dist_under) {
						return false;
					}
					if (overCheck.maxFlow() < dist_over) {
						return false;
					}
				}
			}
		}
		
	}
	return true;
}
template<typename Weight>
void MaxflowDetector<Weight>::printSolution(std::ostream & write_to) {
	Weight f = underapprox_conflict_detector->maxFlow();
	BFSReachability<Weight> r = BFSReachability<Weight>(source,g_under);
	r.update();
	//Filter the proposed flow to exclude superfluous cycles (which, though legal parts of a flow, don't contribute to the overall flow
	//and are usually not a useful thing (and also, can always be safely removed without affecting the logical correctness of the solution))
	//note: it _IS_ possible
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
				if(opt_maxflow_allow_cycles || r.connected(n)){
					if (flow > 0) {
						write_to << "Graph " << outer->getGraphID() << " maxflow " << source << " to " << target
								<< " assigns edge " << g_under.getEdge(e).from << " -> " << g_under.getEdge(e).to << " flow " << flow
								<< "\n";
					}
				}
			}
		}
		
	}
	write_to << "\n";
}

template<typename Weight>
void MaxflowDetector<Weight>::collectDisabledEdges() {
	if (opt_conflict_min_cut_maxflow) {
		buildLearnGraph();
		
		if (learngraph_history_clears != g_over.historyclears || g_over.changed()) {
			//refresh
			overapprox_conflict_detector->update();
			for (int edgeid = 0; edgeid < g_over.edges(); edgeid++) {
				if (g_over.hasEdge(edgeid) && g_over.edgeEnabled(edgeid) ){

					Weight flow = overapprox_conflict_detector->getEdgeFlow(edgeid);
					Weight capacity = overapprox_conflict_detector->getEdgeCapacity(edgeid);
					assert(capacity<=g_over.getWeight(edgeid));
					Weight level0capacity =  outer->hasBitVector(edgeid)?	 outer->getEdgeBV(edgeid).getOver(true):capacity;

					learn_graph.setEdgeWeight(back_edges[edgeid],0x0FF0F0);//capacity of this edge in the backward direction is the same as the forward flow.
					if(flow>0){
						learn_graph.enableEdge(back_edges[edgeid]);
					}else{
						learn_graph.disableEdge(back_edges[edgeid]);
					}

					if (flow<capacity){
						//then this edge is not fully utilized, so it cannot be a limiting factor in the maxflow. Set its weight to infinity.
						learn_graph.setEdgeWeight(edgeid,0x0FF0F0);
						learn_graph.enableEdge(edgeid);
					}else{
						learn_graph.setEdgeWeight(edgeid,(level0capacity-capacity > 0) ? 1:0);
						if(learn_graph.getWeight(edgeid)>0){
							learn_graph.enableEdge(edgeid);
						}else{
							learn_graph.disableEdge(edgeid);
						}
					}
				}else{
					//flow is 0.
					learn_graph.setEdgeWeight(back_edges[edgeid],0);//capacity of this edge in the backward direction is the same as the forward flow, which is 0.
					Weight w =outer->hasBitVector(edgeid)?	 outer->getEdgeBV(edgeid).getOver(true):outer->edge_weights[edgeid];
					learn_graph.setEdgeWeight(edgeid,  (w>0) ? 1:0 );
					learn_graph.disableEdge(back_edges[edgeid]);
					if(learn_graph.getWeight(edgeid)>0){
						learn_graph.enableEdge(edgeid);
					}else{
						learn_graph.disableEdge(edgeid);
					}
				}
			}
			learngraph_history_clears = g_over.historyclears;
			learngraph_history_qhead = g_over.historySize();
		} else {
			for (int i = learngraph_history_qhead; i < g_over.historySize(); i++) {
				int edgeid = g_over.getChange(i).id;
				if (g_over.hasEdge(edgeid) && g_over.edgeEnabled(edgeid) ){

					Weight flow = overapprox_conflict_detector->getEdgeFlow(edgeid);
					Weight capacity = overapprox_conflict_detector->getEdgeCapacity(edgeid);
					assert(capacity<=g_over.getWeight(edgeid));
					Weight level0capacity =  outer->hasBitVector(edgeid)?	 outer->getEdgeBV(edgeid).getOver(true):capacity;

					learn_graph.setEdgeWeight(back_edges[edgeid],0x0FF0F0);//capacity of this edge in the backward direction is the same as the forward flow.
					if(flow>0){
						learn_graph.enableEdge(back_edges[edgeid]);
					}else{
						learn_graph.disableEdge(back_edges[edgeid]);
					}
					if (flow<capacity){
						//then this edge is not fully utilized, so it cannot be a limiting factor in the maxflow. Set its weight to infinity.
						learn_graph.setEdgeWeight(edgeid,0x0FF0F0);
						learn_graph.enableEdge(edgeid);
					}else{
						learn_graph.setEdgeWeight(edgeid,(level0capacity-capacity > 0) ? 1:0);
						if(learn_graph.getWeight(edgeid)>0){
							learn_graph.enableEdge(edgeid);
						}else{
							learn_graph.disableEdge(edgeid);
						}
					}
				}else{
					//flow is 0.
					learn_graph.setEdgeWeight(back_edges[edgeid],0);//capacity of this edge in the backward direction is the same as the forward flow, which is 0.
					Weight w =outer->hasBitVector(edgeid)?	 outer->getEdgeBV(edgeid).getOver(true):outer->edge_weights[edgeid];
					learn_graph.setEdgeWeight(edgeid,  (w>0) ? 1:0 );

					learn_graph.disableEdge(back_edges[edgeid]);
					if(learn_graph.getWeight(edgeid)>0){
						learn_graph.enableEdge(edgeid);
					}else{
						learn_graph.disableEdge(edgeid);
					}
				}
			}
			learngraph_history_qhead = g_over.historySize();
		}
	}
	g_over.updateAlgorithmHistory(this,alg_id,learngraph_history_qhead);
}

template<typename Weight>
void MaxflowDetector<Weight>::collectChangedEdges() {
	if(!(opt_conflict_min_cut_maxflow || opt_decide_theories))
		return;

	dbg_decisions();
	overapprox_conflict_detector->update();

	if (opt_conflict_min_cut_maxflow) {

		buildLearnGraph();
	}
	
	std::vector<int> & changed_edges = overapprox_conflict_detector->getChangedEdges();


	if (opt_decide_theories &&  opt_rnd_order_graph_decisions) {
		/*			static vec<int> tmp_changed;
		 tmp_changed.clear();
		 for(int edge:changed_edges):
		 tmp_changed.p*/
		randomShuffle(rnd_seed, changed_edges);
	}
	//for(int j = 0;j<changed_edges.size();j++){
	//			int edgeid = changed_edges[j];

	//need to deal with changes to bv edge weights, also!
	for (int j = changed_edges.size() - 1; j >= 0; j--) {
		static int iter = 0;
		++iter;
		int edgeid = changed_edges[j];
		if(opt_theory_internal_vsids){
			insertEdgeOrder(edgeid);
		}
		if(opt_decide_theories){
		//if (!is_potential_decision[edgeid]) {
			Lit l = mkLit(outer->getEdgeVar(edgeid), false);
			if((outer->decidable(l) || outer->level(var(l))>0) || (outer->edgeWeightDecidable(edgeid, DetectorComparison::geq,  overapprox_conflict_detector->getEdgeFlow(edgeid) )) ){
				if (overapprox_conflict_detector->getEdgeFlow(edgeid) > 0) {

					is_potential_decision[edgeid] = true;
					if(!in_decision_q[edgeid]){
						in_decision_q[edgeid]=true;
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
					/*	printf("g%d %d mf collect change q (size %d): ", outer->id,iter,potential_decisions_q.size());
						for(int i =0;i<potential_decisions_q.size();i++){
							printf("%d, ",(int) potential_decisions_q[i]);
							break;
						}
						printf("%d, ",(int) potential_decisions_q[potential_decisions_q.size()-1]);
						printf("\n");*/
					}
				}
			}
		}
		//}
		if (opt_conflict_min_cut_maxflow) {
			
			Lit l = mkLit(outer->getEdgeVar(edgeid), false);
			//maintain the residual graph for conflict analysis
			lbool val = outer->value(l);
			if (g_over.hasEdge(edgeid) && g_over.edgeEnabled(edgeid) ){

				Weight flow = overapprox_conflict_detector->getEdgeFlow(edgeid);
				Weight capacity = overapprox_conflict_detector->getEdgeCapacity(edgeid);
				assert(capacity<=g_over.getWeight(edgeid));
				Weight level0capacity =  outer->hasBitVector(edgeid)?	 outer->getEdgeBV(edgeid).getOver(true):capacity;

				learn_graph.setEdgeWeight(back_edges[edgeid],0x0FF0F0);//capacity of this edge in the backward direction is the same as the forward flow.
				if(flow>0){
					learn_graph.enableEdge(back_edges[edgeid]);
				}else{
					learn_graph.disableEdge(back_edges[edgeid]);
				}

				if (flow<capacity){
					//then this edge is not fully utilized, so it cannot be a limiting factor in the maxflow. Set its weight to infinity.
					learn_graph.setEdgeWeight(edgeid,0x0FF0F0);
					learn_graph.enableEdge(edgeid);
				}else{
					learn_graph.setEdgeWeight(edgeid,(level0capacity-capacity > 0) ? 1:0);
					if(learn_graph.getWeight(edgeid)>0){
						learn_graph.enableEdge(edgeid);
					}else{
						learn_graph.disableEdge(edgeid);
					}
				}
			}else{
				//flow is 0.
				learn_graph.setEdgeWeight(back_edges[edgeid],0);//capacity of this edge in the backward direction is the same as the forward flow, which is 0.
				Weight w =outer->hasBitVector(edgeid)?	 outer->getEdgeBV(edgeid).getOver(true):outer->edge_weights[edgeid];
				learn_graph.setEdgeWeight(edgeid,  (w>0) ? 1:0 );

				learn_graph.disableEdge(back_edges[edgeid]);
				if(learn_graph.getWeight(edgeid)>0){
					learn_graph.enableEdge(edgeid);
				}else{
					learn_graph.disableEdge(edgeid);
				}
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
	/*for (int edgeID = 0; edgeID < g_under.edges(); edgeID++) {
		
		if (is_potential_decision[edgeID]) {
			if(potential_decisions.contains(edgeID) || potential_decisions_q.contains(edgeID)){
			 printf("%d ",edgeID);
			 }
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
	}*/
	/*printf("\n");
	 printf("Decisions %d: ",iter);

	 printf("\n");*/
#endif
}

/*template<typename Weight>
Lit MaxflowDetector<Weight>::decideByPath(int level){
	//given an edgeID, decide any undecided edge contributing to its flow.
	if(current_decision_edge<0)
		return lit_Undef;
	if (overapprox_conflict_detector->getEdgeFlow(current_decision_edge)<=0){
		current_decision_edge=-1;
	}
	bool forward = (opt_maxflow_decision_paths==1);
	//explore forward and backwards from this edge to find other edges with flow (in the over approx)
	int lastDecisionEdge = -1;


	if(lastDecisionEdge==-1){
		int node = g_over.getEdge(current_decision_edge).from;
		while(node!=(forward? target:source)){
			for (int i = 0;i<g_over.nDirectedEdges(node,!forward);i++){
				int edgeID = g_over.directedEdges(node,i,!forward).id;
				if(g_over.edgeEnabled(edgeID)){
					if(overapprox_conflict_detector->getEdgeFlow(edgeID)>0){
						//if(!g_under.edgeEnabled(edgeID)){
						if(outer->value(outer->getEdgeVar(edgeID))==l_Undef){
							lastDecisionEdge=edgeID;//in the forward direction, decide the first encountered node, not the last.
							node = forward? target:source;
							break;
						}
						node = g_over.directedEdges(node,i,!forward).node;
						break;
					}
				}
			}
		}
	}
	current_decision_edge=lastDecisionEdge;
	if(lastDecisionEdge>-1){
		Lit l = mkLit(outer->getEdgeVar(lastDecisionEdge), false);
		decideEdge(lastDecisionEdge, level, true);
		return l;
	}
	current_decision_edge=-1;
	return lit_Undef;
}*/

template<typename Weight>
bool MaxflowDetector<Weight>::decideEdgeWeight(int edgeID, Weight & store, DetectorComparison & op){
	store =  overapprox_conflict_detector->getEdgeFlow(edgeID);
	op=DetectorComparison::geq;
	return store>0;
}

template<typename Weight>
void MaxflowDetector<Weight>::undecideEdgeWeight(int edgeid){
	if(opt_theory_internal_vsids){
		insertEdgeOrder(edgeid);
	}
	if(!opt_decide_theories)
		return;
	if((outer->hasBitVector(edgeid) || is_potential_decision[edgeid]) && !in_decision_q[edgeid]){

		//problem: this check is applied before any backtracking might occur (if lazy backtracking is not applied).
		if (overapprox_conflict_detector->getEdgeFlow(edgeid) > 0) {			//this check is optional

			if(!in_decision_q[edgeid]){
				in_decision_q[edgeid]=true;
				if (opt_maxflow_decisions_q == 0)
					potential_decisions_q.insertBack(edgeid);
				else if (opt_maxflow_decisions_q == 1) {
					potential_decisions_q.insertBack(edgeid);//insert in LIFO order, no FIFO, because we are unwinding the decisions
				} else if (opt_maxflow_decisions_q == 2) {
					potential_decisions_q.insert(edgeid);			//insert in FIFO order instead
				} else if (opt_maxflow_decisions_q == 3) {
					potential_decisions_q.insertBack(edgeid);//insert in LIFO order, no FIFO, because we are unwinding the decisions
				} else if (opt_maxflow_decisions_q == 4) {
					potential_decisions_q.insert(edgeid);//insert in LIFO order, no FIFO, because we are unwinding the decisions
				}
			}
		/*	printf("g%d mf backtrack q (size %d): ", outer->id,potential_decisions_q.size());
			for(int i =0;i<potential_decisions_q.size();i++){
				printf("%d, ",(int) potential_decisions_q[i]);
				break;
			}
			if(potential_decisions_q.size()){
				printf("%d, ",(int) potential_decisions_q[potential_decisions_q.size()-1]);
			}
			printf("\n");*/
		} else {
			is_potential_decision[edgeid] = false;			//discard this edge from the set of potential decisions
			//printf("undecide remove from q %d\n", edgeid);
		}
	}
}

template<typename Weight>
void MaxflowDetector<Weight>::undecide(Lit l) {
	static int iter = 0;

	if(outer->isEdgeVar(var(l))){
		++iter;
		int edgeid = outer->getEdgeID(var(l));
		undecideEdgeWeight(edgeid);


	}



	/*printf("g%d %d undecide %d q: ", outer->id,iter,dimacs(l));
	for(int i =0;i<potential_decisions_q.size();i++){
		printf("%d, ",(int) potential_decisions_q[i]);
	}
	printf("\n");*/
}
template<typename Weight>
void MaxflowDetector<Weight>::suggestDecision(Lit l){
	if(outer->isEdgeVar(var(l))){
		int edgeID = outer->getEdgeID(var(l));
		priority_decisions.push(edgeID);
	}
}

template<typename Weight>
Lit MaxflowDetector<Weight>::decide() {
	static int it = 0;
	if (++it == 22) {
		int a = 1;
	}
	double startdecidetime = rtime(2);
	auto * over = overapprox_conflict_detector;
	auto * under = underapprox_conflict_detector;
	

	//typically, only one of the priority decisions will be selected (and, in doing so, satisfy the learn clause that filled the priority queue)
	while (priority_decisions.size() > 0) {
		int edgeID =priority_decisions.last();
		priority_decisions.pop();

		Lit l = mkLit(outer->getEdgeVar(edgeID), false);
		if ((outer->decidable(l) || outer->edgeWeightDecidable(edgeID, DetectorComparison::geq,  overapprox_conflict_detector->getEdgeFlow(edgeID)) ) && over->getEdgeFlow(edgeID)>0) {
			n_stats_priority_decisions++;
			double post_time = rtime(2);
			stats_decide_time += post_time - startdecidetime;
			stats_redecide_time += post_time - startdecidetime;
			if(opt_theory_priority_clear)
				priority_decisions.clear();//at most one decision per conflict clause.
			return l;
		}
	}

	while(order_heap.size()){
		int edgeID = order_heap.removeMin();
		Lit l = mkLit(outer->getEdgeVar(edgeID), false);
		if ((outer->decidable(l) || outer->edgeWeightDecidable(edgeID, DetectorComparison::geq,  overapprox_conflict_detector->getEdgeFlow(edgeID)) ) && over->getEdgeFlow(edgeID) > 0) {
			n_stats_vsids_decisions++;
			double post_time = rtime(2);
			stats_decide_time += post_time - startdecidetime;
			stats_redecide_time += post_time - startdecidetime;
			return l;
		}
	}

	if (opt_lazy_maxflow_decisions) {

		collectChangedEdges();
		//DEBUG:check that the decision_q contains the right elements
#ifndef NDEBUG
		if(outer->hasBitVectorEdges()){
		for(int edgeID = 0;edgeID<g_over.edges();edgeID++){
			if(g_over.hasEdge(edgeID) && g_over.edgeEnabled(edgeID)){
				Var v = outer->getEdgeVar(edgeID);
				lbool val = outer->value(v);
				Weight flow =  overapprox_conflict_detector->getEdgeFlow(edgeID);
				int bvID = outer->getEdgeBV(edgeID).getID();
				Weight over_weight = outer->getEdgeBV(edgeID).getOver();
				Weight under_weight = outer->getEdgeBV(edgeID).getUnder();
				if(opt_decide_graph_bv){
					if((val==l_Undef && flow>0) || (flow>0 && flow>under_weight)){
						if(!in_decision_q[edgeID]){
							throw std::runtime_error("Error in maxflow");
						}
					}
				}else{
					if((val==l_Undef && flow>0)){
						if(!in_decision_q[edgeID]){
							throw std::runtime_error("Error in maxflow");
						}
					}
				}
			}
		}
		}
#endif


		Lit decision = lit_Undef;

		while (potential_decisions_q.size() > 0) {
			int edgeID;
			//bool path_decision=false;
			if (opt_maxflow_decisions_q == 0) {
				edgeID = potential_decisions_q.peekBack();
				//path_decision = potential_decisions_q.peekBack().path_decision;
				potential_decisions_q.popBack();
			} else {
				edgeID = potential_decisions_q.peek();
				//path_decision = potential_decisions_q.peek().path_decision;
				potential_decisions_q.pop();
			}
			assert(in_decision_q[edgeID]);
			in_decision_q[edgeID]=false;
			assert(outer->hasBitVector(edgeID) || is_potential_decision[edgeID]);
			Lit l = mkLit(outer->getEdgeVar(edgeID), false);
			if ((outer->decidable(l) || outer->edgeWeightDecidable(edgeID, DetectorComparison::geq,  overapprox_conflict_detector->getEdgeFlow(edgeID)))  && over->getEdgeFlow(edgeID) > 0) {
				//decideEdge(edgeID, true);
				if(edgeID>446){
					int a =1;
				}else{
					int a =1;
				}
				decision = l;
				//printf("decide edge %d\n", edgeID);
				break;
			} else if (outer->value(l) != l_False) {
				if (over->getEdgeFlow(edgeID) > 0) {			//this check is optional
					//decideEdge(edgeID, true);
				//	printf("skip decision keep edge %d\n", edgeID);
				} else {
					is_potential_decision[edgeID] = false;
				//	printf("skip decision remove edge %d\n", edgeID);
				}
			} else {
				assert(over->getEdgeFlow(edgeID) == 0);
				is_potential_decision[edgeID] = false;
				//printf("skip decision remove edge %d\n", edgeID);
			}
		}
		//}

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

template class Monosat::MaxflowDetector<int> ;
template class Monosat::MaxflowDetector<long> ;
template class Monosat::MaxflowDetector<double> ;
#include <gmpxx.h>
template class Monosat::MaxflowDetector<mpq_class> ;
