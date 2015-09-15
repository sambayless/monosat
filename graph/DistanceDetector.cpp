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
#include <dgl/RamalReps.h>
#include "dgl/EdmondsKarp.h"
#include "dgl/EdmondsKarpAdj.h"
#include "dgl/KohliTorr.h"
#include "dgl/EdmondsKarpDynamic.h"
#include "dgl/Dinics.h"
#include "dgl/DinicsLinkCut.h"
#include <graph/DistanceDetector.h>
#include <graph/GraphTheory.h>
#include <mtl/Rnd.h>
#include <mtl/Vec.h>
#include <utils/Options.h>
//#include "dgl/UnweightedDistance.h"
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <iomanip>
using namespace Monosat;

template<typename Weight>
DistanceDetector<Weight>::DistanceDetector(int _detectorID, GraphTheorySolver<Weight> * outer,
		DynamicGraph<Weight>  &_g, DynamicGraph<Weight>  &_antig, int from, double seed) :
		Detector(_detectorID), outer(outer), g_under(_g), g_over(_antig), source(from), rnd_seed(seed) {
	max_unweighted_distance = 0;
	rnd_path = NULL;
	
	constraintsBuilt = -1;
	first_reach_var = var_Undef;
	stats_pure_skipped = 0;
	if (distalg == DistAlg::ALG_SAT) {
		positiveReachStatus = nullptr;
		negativeReachStatus = nullptr;
		positiveDistanceStatus=nullptr;
		negativeDistanceStatus=nullptr;
		underapprox_unweighted_distance_detector = nullptr;
		overapprox_unweighted_distance_detector = nullptr;
		underapprox_path_detector = nullptr;
		unweighted_underprop_marker = CRef_Undef;
		unweighted_overprop_marker = CRef_Undef;
		
		weighted_underprop_marker = CRef_Undef;
		weighted_overprop_marker = CRef_Undef;
		//we are just going to directly enforce these constraints in the original SAT solver, so _nothing_ will end up happening in this detector (except for creating the clauses needed to enforce these constraints).
		
		return;
	}
	
	if (opt_use_random_path_for_decisions) {
		rnd_weight.clear();
		rnd_path = new WeightedDijkstra<Weight,double>(from, _antig, rnd_weight);
		for (int i = 0; i < outer->edge_list.size(); i++) {
			double w = drand(rnd_seed);
			
			rnd_weight.push_back(w);
		}
		
	}
	/* if(opt_use_optimal_path_for_decisions){
	 opt_path = new WeightedDijkstra< OptimalWeightEdgeStatus >(from,_antig,opt_weight);
	 }*/
	positiveReachStatus = new DistanceDetector<Weight>::ReachStatus(*this, true);
	negativeReachStatus = new DistanceDetector<Weight>::ReachStatus(*this, false);
	
	//select the unweighted distance detectors
	if (distalg == DistAlg::ALG_DISTANCE) {
		if (outer->all_edges_unit) {
			if (!opt_encode_dist_underapprox_as_sat)
				underapprox_unweighted_distance_detector = new UnweightedBFS<Weight,typename DistanceDetector<Weight>::ReachStatus>(from,
					_g, *(positiveReachStatus), 0);
			overapprox_unweighted_distance_detector = new UnweightedBFS<Weight,typename DistanceDetector<Weight>::ReachStatus>(from,
					_antig, *(negativeReachStatus), 0);
			if(underapprox_unweighted_distance_detector)
				underapprox_path_detector = underapprox_unweighted_distance_detector;
			else{
				underapprox_path_detector = new UnweightedBFS<Weight,typename DistanceDetector<Weight>::ReachStatus>(from,
									_g, *(positiveReachStatus), 0);
			}

		} else {
			if (!opt_encode_dist_underapprox_as_sat){
				underapprox_unweighted_distance_detector = new UnweightedDijkstra<Weight,typename DistanceDetector<Weight>::ReachStatus>(
						from, _g, *positiveReachStatus, 0);
			}
			overapprox_unweighted_distance_detector = new UnweightedDijkstra<Weight,typename DistanceDetector<Weight>::ReachStatus>(
					from, _antig, *negativeReachStatus, 0);
			underapprox_path_detector = new UnweightedBFS<Weight,Distance<int>::NullStatus>(from, _g, Distance<int>::nullStatus, 0);
		}
		
		/*	if(opt_conflict_shortest_path)
		 reach_detectors.last()->positive_dist_detector = new Dijkstra<PositiveEdgeStatus>(from,g);*/
	} else if (distalg == DistAlg::ALG_RAMAL_REPS) {
		if (!opt_encode_dist_underapprox_as_sat){
			/*underapprox_unweighted_distance_detector = new UnweightedRamalReps<Weight,
					typename DistanceDetector<Weight>::ReachStatus>(from, _g, *(positiveReachStatus), 0);*/
			underapprox_unweighted_distance_detector = new UnweightedDijkstra<Weight,typename DistanceDetector<Weight>::ReachStatus>(
							from, _g, *positiveReachStatus, 0);
		}
		overapprox_unweighted_distance_detector =
				new UnweightedRamalReps<Weight,typename DistanceDetector<Weight>::ReachStatus>(from, _antig,
						*(negativeReachStatus), 0);
		//fix this? Can RamalReps report paths?
		underapprox_path_detector = new UnweightedBFS<Weight,Distance<int>::NullStatus>(from, _g, Distance<int>::nullStatus, 0);
	} else {
		if (!opt_encode_dist_underapprox_as_sat){
			underapprox_unweighted_distance_detector =
					new UnweightedDijkstra<Weight,typename DistanceDetector<Weight>::ReachStatus>(from, _g, *positiveReachStatus,
							0);
		}
		overapprox_unweighted_distance_detector =
				new UnweightedDijkstra<Weight,typename DistanceDetector<Weight>::ReachStatus>(from, _antig,
						*negativeReachStatus, 0);

		if(underapprox_unweighted_distance_detector)
			underapprox_path_detector = underapprox_unweighted_distance_detector;
		else{
			underapprox_path_detector =	new UnweightedDijkstra<Weight,typename DistanceDetector<Weight>::ReachStatus>(from, _g, *positiveReachStatus,0);
		}
		//reach_detectors.last()->positive_dist_detector = new Dijkstra(from,g);
	}
	
	//select the _weighted_ distance detectors
	positiveDistanceStatus = new DistanceDetector<Weight>::DistanceStatus(*this, true);
	negativeDistanceStatus = new DistanceDetector<Weight>::DistanceStatus(*this, false);
	
	if (outer->hasBitVectorEdges()){
		printf("Note: falling back on Dijkstra for shortest paths, because edge weights are bitvectors\n");
		//ramel reps doesn't support bvs yet
		underapprox_weighted_distance_detector =
		new Dijkstra<Weight, typename DistanceDetector<Weight>::DistanceStatus>(from, _g,
				*positiveDistanceStatus, 0);
		overapprox_weighted_distance_detector = new Dijkstra<Weight, typename DistanceDetector<Weight>::DistanceStatus>(
				from, _antig,  *negativeDistanceStatus, 0);
		underapprox_weighted_path_detector = underapprox_weighted_distance_detector;
	}else if (  distalg == DistAlg::ALG_RAMAL_REPS) {

		underapprox_weighted_distance_detector =
				new RamalReps<Weight, typename DistanceDetector<Weight>::DistanceStatus>(from, _g,
						*(positiveDistanceStatus), 0);
		overapprox_weighted_distance_detector =
				new RamalReps<Weight, typename DistanceDetector<Weight>::DistanceStatus>(from, _antig,
						*(negativeDistanceStatus), 0);
		underapprox_weighted_path_detector = new Dijkstra<Weight>(from, _g);
	} else {
		underapprox_weighted_distance_detector =
				new Dijkstra<Weight, typename DistanceDetector<Weight>::DistanceStatus>(from, _g,
						*positiveDistanceStatus, 0);
		overapprox_weighted_distance_detector = new Dijkstra<Weight, typename DistanceDetector<Weight>::DistanceStatus>(
				from, _antig,  *negativeDistanceStatus, 0);
		underapprox_weighted_path_detector = underapprox_weighted_distance_detector;
	}
	
	if (opt_conflict_min_cut) {
		if (mincutalg == MinCutAlg::ALG_EDKARP_DYN) {
			conflict_flow = new EdmondsKarpDynamic<long>(outer->cutGraph,  source, 0);
		} else if (mincutalg == MinCutAlg::ALG_EDKARP_ADJ) {
			conflict_flow = new EdmondsKarpAdj<long>(outer->cutGraph,  source, 0);
		} else if (mincutalg == MinCutAlg::ALG_DINITZ) {
			conflict_flow = new Dinitz<long>(outer->cutGraph,  source, 0);
		} else if (mincutalg == MinCutAlg::ALG_DINITZ_LINKCUT) {
			//link-cut tree currently only supports ints
			conflict_flow = new Dinitz<long>(outer->cutGraph,  source, 0);
			
		} else if (mincutalg == MinCutAlg::ALG_KOHLI_TORR) {
			if (opt_use_kt_for_conflicts) {
				conflict_flow = new KohliTorr<long>(outer->cutGraph, source, 0,
						opt_kt_preserve_order);
			} else
				conflict_flow = new EdmondsKarpDynamic<long>(outer->cutGraph,  source, 0);
		} else {
			conflict_flow = new EdmondsKarpAdj<long>(outer->cutGraph,  source, 0);
		}
	}
	
	unweighted_underprop_marker = outer->newReasonMarker(getID());
	unweighted_overprop_marker = outer->newReasonMarker(getID());
	
	weighted_underprop_marker = outer->newReasonMarker(getID());
	weighted_overprop_marker = outer->newReasonMarker(getID());
}

template<typename Weight>
void DistanceDetector<Weight>::buildUnweightedSATConstraints(bool onlyUnderApprox, int within_steps) {
	if (within_steps < 0)
		within_steps = g_under.nodes();
	if (within_steps > g_under.nodes())
		within_steps = g_under.nodes();
	if (within_steps > g_under.edges())
		within_steps = g_under.edges();
	if (constraintsBuilt >= within_steps)
		return;
	
	if (onlyUnderApprox && opt_encode_dist_underapprox_as_sat==2){
		//this is only correct for directed graphs.
		//To do the same for undirected graps, also non-deterministically set the direction of each edge.
		Lit True = mkLit(outer->newVar());
		outer->addClause(True);
		BitVector<Weight> one = outer->comparator->newBitvector(-1,outer->getEdgeWeightBitWidth() ,1);
		for(int  to = 0;to<unweighted_dist_lits.size();to++){
			if(unweighted_dist_lits[to].size() && to!=source){

				//For each target node, create an Erez-Nadel style nondeterministic graph
				//then non-deterministically pick an acyclic path
				//and assert that it's length is less than the target length
				vec<Lit> induced_graph_edges;
				induced_graph_edges.growTo(g_under.edges(),lit_Undef);
				for (int edgeID = 0;edgeID<g_under.edges();edgeID++){
					if(g_over.hasEdge(edgeID) ){//&& g_over.edgeEnabled(edgeID)
						int from = g_under.getEdge(edgeID).from;
						int to = g_under.getEdge(edgeID).to;
						Lit edge_enabled = mkLit(outer->getEdgeVar(edgeID),false);
						Lit induced_edge_enabled = mkLit(outer->newVar(), false);
						induced_graph_edges[edgeID]=induced_edge_enabled;
						//If the edge is disabled, then our induced graph edge must also be disabled (but _not_ vice versa).
						outer->addClause(edge_enabled,~induced_edge_enabled);
					}
				}

				vec<BitVector<Weight>> node_distances;
				//introduce a bitvector distance for each node
				for(int n = 0;n<g_under.nodes();n++){
					if (n==source){
						node_distances.push(outer->newBV(0));//source has constant 0 distance
					}else{
						BitVector<Weight> bv = outer->newBV();
						node_distances.push(bv);
					}
				}

				//now, non-deterministically select an acyclic path from source out of this graph, asserting that it ends at the target node
				Lit reaches_to=lit_Undef;
				//every vertex (other than source or to) has at most one in-coming edge enabled, and if it has such an edge, then it has exactly one outgoing edge enabled (else no outgoing edges)
				//source has no outgoing edge, and exactly one outgoing edge
				for(int n = 0;n<g_under.nodes();n++){
					{
						//Lit any_incoming_enabled= mkLit( outer->newVar(),false); //~True;
						//vec<Lit> any_incoming_clause;
						//any_incoming_clause.push(~any_incoming_enabled);
						Lit any_incoming_enabled= ~True;
						BitVector<Weight> dist = node_distances[n];
						if(n!=source){
							BitVector<Weight> dist =opt_sat_distance_encoding_unconstrained_default?  outer->newBV() : outer->newBV(0);//should this be unconstrained, or zero?
							for(int i = 0;i<g_over.nIncoming(n);i++){
								int edgeID = g_over.incoming(n,i).id;
								int from = g_over.incoming(n,i).node;
								if(from==n)
									continue;//no need to consider self loops
								Lit l= induced_graph_edges[edgeID];
								//If any prior incoming edge is enabled, this one must be disabled
								outer->addClause(~any_incoming_enabled,~l);

								//if any prior incoming edge was enabled, or l is enabled, then next_incoming is enabled
								Lit next_incoming = mkLit(outer->newVar(),false);

								outer->addClause(~any_incoming_enabled,next_incoming);
								outer->addClause(~l,next_incoming);
								//if both incoming and l are false, then next incoming is false
								outer->addClause(l,any_incoming_enabled,~next_incoming);
								any_incoming_enabled=next_incoming;

								BitVector<Weight> sum = outer->newBV();

								//this is an unweighted constraint, so just add 1
								outer->comparator->newAdditionBV(sum.getID(), node_distances[from].getID(),one.getID());
								BitVector<Weight> new_dist = outer->newBV();
								//If this incoming edge is enabled, then the distance is the cost of this edge
								outer->comparator->newConditionalBV(l,new_dist.getID(),sum.getID(),dist.getID());

								dist = new_dist;
							}
							outer->comparator->makeEquivalent(dist.getID(),node_distances[to].getID());
						}else{
							//all incoming edges must be disabled
							for(int i = 0;i<g_over.nIncoming(n);i++){
								int edgeID = g_over.incoming(n,i).id;
								outer->addClause(~induced_graph_edges[edgeID]);
							}
						}
						//any incoming enabled is true if any incoming edges are enabled
						//also, at most one incoming edge can possibly be enabled.
						if(n==to){
							reaches_to=any_incoming_enabled;
						}
						Lit any_outgoing_enabled=~True;
						//exactly one outgoing edge can be enabled, but _only_ if any incoming edges are enabled.
						if(n!=to){
							for(int i = 0;i<g_over.nIncident(n);i++){
								int edgeID = g_over.incident(n,i).id;
								Lit l= induced_graph_edges[edgeID];

								//if no incoming edge was enabled, then all outgoing edges must be disabled
								outer->addClause(any_incoming_enabled,~l);

								//If any prior incoming edge is enabled, this one must be disabled
								outer->addClause(~any_outgoing_enabled,~l);

								//if any prior incoming edge was enabled, or l is enabled, then next_incoming is enabled
								Lit next_outgoing = mkLit(outer->newVar(),false);

								outer->addClause(~any_outgoing_enabled,next_outgoing);
								outer->addClause(~l,next_outgoing);
								//if both incoming nor l are false, then next incoming is false
								outer->addClause(l,any_outgoing_enabled,~next_outgoing);
								any_outgoing_enabled=next_outgoing;
							}
						}else{
							//all outgoing edges must be disabled
							for(int i = 0;i<g_over.nIncident(n);i++){
								int edgeID = g_over.incident(n,i).id;
								outer->addClause(~induced_graph_edges[edgeID]);
							}
						}
						if(n==source){
							outer->addClause(any_outgoing_enabled);
						}else if (n==to){
							outer->addClause(any_incoming_enabled);
						}else{
							//if and only if any incoming edge is enabled, then exactly one outgoing edge must be enabled
							outer->addClause(~any_incoming_enabled,any_outgoing_enabled);
							outer->addClause(any_incoming_enabled,~any_outgoing_enabled);
						}
					}
				}
				//for each distance, assert that it's lit is true iff the distance on the path is <= 'distance'
				//following Erez and Nadel 2015's non-graph aware encoding, use bitvectors
				//Var lit = unweighted_dist_lits[to].l;
				//int distance = unweighted_dist_lits[to].min_unweighted_distance;

				BitVector<Weight> dist = node_distances[to];
				for(int i = 0;i<unweighted_dist_lits[to].size();i++){
					Lit l = unweighted_dist_lits[to][i].l;
					int min_unweighted_distance = unweighted_dist_lits[to][i].min_unweighted_distance;
					BitVector<Weight> comparison =  outer->comparator->newBitvector(-1,outer->getEdgeWeightBitWidth() ,min_unweighted_distance);
					Lit conditional = outer->comparator->newComparison(Comparison::leq, dist.getID() ,comparison.getID());

					Lit s_l = outer->toSolver(l);
					Lit s_conditional =  outer->comparator->toSolver(conditional);
					Lit s_reaches = outer->toSolver(reaches_to);

					//l is true if the distance is <= comparison, AND their exists a path to node
					outer->addClauseToSolver(~s_l,s_conditional);
					outer->addClauseToSolver(~s_l,s_reaches);
					outer->addClauseToSolver(s_l,~s_conditional,~s_reaches);
				}
			}
		}
	}else{
		assert(outer->decisionLevel() == 0);
		vec<Lit> c;

		if (constraintsBuilt <= 0) {
			constraintsBuilt = 0;
			unweighted_dist_lits.push();
			Lit True = mkLit(outer->newVar());
			outer->addClause(True);
			assert(outer->value(True)==l_True);
			Lit False = ~True;
			for (int i = 0; i < g_under.nodes(); i++) {
				unweighted_sat_lits[0].push(False);
			}
			unweighted_sat_lits[0][source] = True;
		}

		vec<Lit> reaches;

		//bellman-ford:
		for (int i = constraintsBuilt; i < within_steps; i++) {
			unweighted_sat_lits.last().copyTo(reaches);
			unweighted_sat_lits.push();
			reaches.copyTo(unweighted_sat_lits.last());
			assert(outer->value( reaches[source])==l_True);
			//For each edge:
			for (int j = 0; j < g_under.nodes(); j++) {
				Lit r_cur = reaches[j];

				for (Edge & e : outer->inv_adj[j]) {
					if (outer->value(unweighted_sat_lits.last()[e.to]) == l_True) {
						//do nothing
					} else if (outer->value(reaches[e.from]) == l_False) {
						//do nothing
					} else {
						Lit l = mkLit(e.v, false);
						Lit r = mkLit(outer->newVar(), false);

						c.clear();
						c.push(~r);
						c.push(reaches[e.to]);
						c.push(l); //r -> (e.l or reaches[e.to])
						outer->addClause(c);
						c.clear();
						c.push(~r);
						c.push(reaches[e.to]);
						c.push(reaches[e.from]); //r -> (reaches[e.from]) or reaches[e.to])
						outer->addClause(c);
						c.clear();
						c.push(r);
						c.push(~reaches[e.to]); //~r -> ~reaches[e.to]
						outer->addClause(c);
						c.clear();
						c.push(r);
						c.push(~reaches[e.from]);
						c.push(~l); //~r -> (~reaches[e.from] or ~e.l)
						outer->addClause(c);
						r_cur = r;

					}
				}
				unweighted_sat_lits.last()[j] = r_cur; //reaches[e.to] == (var & reaches[e.from])| reaches[e.to];
			}

		}
		
		constraintsBuilt = within_steps;
	}
}

template<typename Weight>
void DistanceDetector<Weight>::addWeightedShortestPathLit(int from, int to, Var outer_reach_var,
		Weight within_distance, bool strictComparison) {
	g_under.invalidate();
	g_over.invalidate();
	Var reach_var = outer->newVar(outer_reach_var, getID());
	assert(from == source);
	weighted_dist_lits.push(WeightedDistLit { mkLit(reach_var), to, within_distance,strictComparison });
	//sort(weighted_dist_lits);
}

template<typename Weight>
void DistanceDetector<Weight>::addWeightedShortestPathBVLit(int from, int to, Var outer_reach_var,
		const BitVector<Weight>  &bv, bool strictComparison) {
	g_under.invalidate();
	g_over.invalidate();
	Var reach_var = outer->newVar(outer_reach_var, getID());
	assert(from == source);
	weighted_dist_bv_lits.push(WeightedDistBVLit { mkLit(reach_var), to, bv, strictComparison });
	//sort(weighted_dist_lits);
}


template<typename Weight>
void DistanceDetector<Weight>::addUnweightedShortestPathLit(int from, int to, Var outer_reach_var, int within_steps) {
	g_under.invalidate();
	g_over.invalidate();
	Var reach_var = outer->newVar(outer_reach_var, getID());
	
	if (first_reach_var == var_Undef) {
		first_reach_var = reach_var;
	} else {
		assert(reach_var > first_reach_var);
	}
	assert(from == source);
	assert(within_steps >= 0);
	if (within_steps > g_under.nodes())
		within_steps = g_under.nodes();
	if (within_steps < 0)
		within_steps = g_under.nodes();
	
	if (within_steps >= max_unweighted_distance) {
		max_unweighted_distance = within_steps;
		if (opt_compute_max_distance) {
			underapprox_unweighted_distance_detector->setMaxDistance(max_unweighted_distance);
			overapprox_unweighted_distance_detector->setMaxDistance(max_unweighted_distance);
			//positive_path_detector->setMaxDistance(max_distance);
		}
	}
	
	Lit reachLit = mkLit(reach_var, false);

	
	while (unweighted_dist_lits.size() <= to)
		unweighted_dist_lits.push();
	
	bool found = false;
	for (int i = 0; i < unweighted_dist_lits[to].size(); i++) {
		if (unweighted_dist_lits[to][i].min_unweighted_distance == within_steps) {
			found = true;
			Lit r = unweighted_dist_lits[to][i].l;
			//force equality between the new lit and the old reach lit, in the SAT solver
			outer->makeEqual(r, reachLit);
			/*outer->S->addClause(~r, reachLit);
			 outer->S->addClause(r, ~reachLit);*/
		}
	}
	
	if (!found) {
		unweighted_dist_lits[to].push();
		unweighted_dist_lits[to].last().l = reachLit;
		unweighted_dist_lits[to].last().min_unweighted_distance = within_steps;

		while (reach_lit_map.size() <= reach_var - first_reach_var) {
			reach_lit_map.push({-1,-1});
		}
		reach_lit_map[reach_var - first_reach_var] = {to,within_steps};
	}
	
}
template<typename Weight>
void DistanceDetector<Weight>::ReachStatus::setReachable(int u, bool reachable) {
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
void DistanceDetector<Weight>::ReachStatus::setMininumDistance(int u, bool reachable, int distance) {
	assert(reachable == (distance < detector.g_under.nodes()));
	if (distance <= detector.g_under.nodes()) {
		setReachable(u, reachable);
	}
	
	if (u < detector.unweighted_dist_lits.size() && detector.unweighted_dist_lits[u].size()) {
		assert(distance >= 0);
		
		//for(int i = 0;i<detector.unweighted_dist_lits[u].size();i++){
		//int& min_distance =  detector.unweighted_dist_lits[u][i].min_unweighted_distance;
		
		//Lit l = detector.unweighted_dist_lits[u][i].l;
		if (!detector.is_changed[u]) {
			
			if (!polarity) {				// && (!reachable || (distance>min_distance))){
				//lbool assign = detector.outer->value(l);
				//if( assign!= l_False ){
				detector.is_changed[u] = true;
				detector.changed.push( { u });						//{var(l),u,min_distance});
				//}
			} else if (polarity) {						// && reachable && (distance<=min_distance)){
				//lbool assign = detector.outer->value(l);
				//if( assign!= l_True ){
				detector.is_changed[u] = true;
				detector.changed.push( { u });						//{var(l),u,min_distance});
				//}
			}
		}
		//}
	}
	
}

template<typename Weight>
void DistanceDetector<Weight>::DistanceStatus::setReachable(int u, bool reachable) {
	
}
template<typename Weight>
void DistanceDetector<Weight>::DistanceStatus::setMininumDistance(int u, bool reachable, Weight & distance) {
	
}

template<typename Weight>
void DistanceDetector<Weight>::buildUnweightedDistanceLEQReason(int node, vec<Lit> & conflict) {
	//drawFull();
	Reach & d = *underapprox_path_detector;
	stats_unweighted_leq_reasons++;
	stats_under_conflicts++;
	double starttime = rtime(2);
	d.update();
	assert(outer->dbg_reachable(d.getSource(), node));
	/*
	 if(!outer->dbg_reachable(d.getSource(),node)){
	 outer->drawFull();

	 d.drawFull();

	 assert(false);
	 }*/
	assert(underapprox_unweighted_distance_detector->connected(node));
	/*if(!d.connected_unchecked(node)){
	 exit(3);
	 }*/
	assert(d.connected_unchecked(node));
	//if(opt_learn_reaches ==0 || opt_learn_reaches==2)
	{
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
	}/*else{
	 //Instead of a complete path, we can learn reach variables, if they exist
	 int u = node;
	 int p;
	 while(( p = d.previous(u)) != -1){
	 Edge edg = outer->edges[p][u];
	 Var e =outer->edges[p][u].v;
	 lbool val = outer->value(e);
	 assert(outer->value(e)==l_True);
	 conflict.push(mkLit(e, true));
	 u = p;
	 if( u< reach_lits.size() && reach_lits[u]!=lit_Undef && outer->value(reach_lits[u])==l_True && outer->S->level(var(reach_lits[u]))< outer->S->decisionLevel()){
	 //A potential (fixed) problem with the above: reach lit can be false, but have been assigned after r in the trail, messing up clause learning if this is a reason clause...
	 //This is avoided by ensuring that L is lower level than the conflict.
	 Lit l =reach_lits[u];
	 assert(outer->value(l)==l_True);
	 conflict.push(~l);
	 break;
	 }
	 }

	 }*/
	outer->num_learnt_paths++;
	outer->learnt_path_clause_length += (conflict.size() - 1);
	
	stats_under_clause_length += (conflict.size() - 1);
	double elapsed = rtime(2) - starttime;
	stats_under_conflict_time += elapsed;
	
}
template<typename Weight>
void DistanceDetector<Weight>::buildUnweightedDistanceGTReason(int node, int within_steps, vec<Lit> & conflict) {
	static int it = 0;
	stats_unweighted_gt_reasons++;
	stats_over_conflicts++;
	double starttime = rtime(2);
	++it;
	int u = node;
	bool reaches = overapprox_unweighted_distance_detector->connected(node);
	
	if (!reaches && within_steps>=g_over.nodes() && opt_conflict_min_cut && conflict_flow) {

		g_over.drawFull();
		cut.clear();
		long f;
		
		assert(conflict_flow->getSource() == source);
		conflict_flow->setSink(node);
		f = conflict_flow->minCut(cut);
		assert(f == cut.size()); //because edges are only ever infinity or 1
		if(f != cut.size()){
			throw  std::runtime_error("Bad learnt cut");
		}
		for (int i = 0; i < cut.size(); i++) {
			MaxFlowEdge e = cut[i];
			int cut_id = e.id;
			assert(cut_id % 2 == 0);
			Lit l = mkLit(outer->getEdgeVar(cut_id / 2), false);
			assert(outer->value(l)==l_False);
			conflict.push(l);
		}
		outer->learnt_cut_clause_length += (conflict.size() - 1);
		
		stats_over_clause_length += (conflict.size() - 1);
		double elapsed = rtime(2) - starttime;
		stats_over_conflict_time += elapsed;
		return;
	}
	
	//drawFull( non_reach_detectors[detector]->getSource(),u);
	//assert(outer->dbg_distance( source,u));
	
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
				assert(outer->inv_adj[u][i].to == u);
				//Note: the variable has to not only be assigned false, but assigned false earlier in the trail than the reach variable...
				int edge_num = outer->getEdgeID(v);				// v-outer->min_edge_var;
						
				if (from == u) {
					assert(outer->edge_list[edge_num].to == u);
					assert(outer->edge_list[edge_num].from == u);
					continue;				//Self loops are allowed, but just make sure nothing got flipped around...
				}
				assert(from != u);
				
				if (has_unweighted_shortest_paths_overapprox && reaches) {
					//This is an optional optimization: if we know that even with all possible edges enabled, the shortest path to from + 1 is >= than the current distance to this node, enabling this edge cannot decrease the shortest path,
					//and so we don't need to consider this edge
					int current_dist = overapprox_unweighted_distance_detector->distance(u);
					assert(current_dist > 0);
					int over_approx_dist = unweighted_over_approx_shortest_paths[from] + 1;
					if (over_approx_dist >= current_dist) {
						//following this edge cannot shorten this path, so skip it.
						stats_gt_unweighted_edges_skipped++;
						continue;
					}
				}
				
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
	
	stats_over_clause_length += (conflict.size() - 1);
	double elapsed = rtime(2) - starttime;
	stats_over_conflict_time += elapsed;
	
}

template<typename Weight>
void DistanceDetector<Weight>::printSolution(std::ostream& write_to) {
	
	vec<bool> to_show;
	to_show.growTo(g_under.nodes());
	for (auto & w : weighted_dist_lits) {
		to_show[w.u] = true;
	}
	for (auto & w : weighted_dist_bv_lits) {
		to_show[w.u] = true;
	}
	underapprox_weighted_path_detector->update();
	for (int to = 0; to < g_under.nodes(); to++) {
		if (!to_show[to])
			continue;
		
		Distance<Weight> & d = *underapprox_weighted_path_detector;
		d.update();
		Weight & actual_dist = d.distance(to);
		write_to << "Shortest Weighted Path " << source << "->" << to << " is " << actual_dist << ": ";
		//std::cout<< "Weighted Distance Constraint " <<dimacs(w.l) << " (" << source <<"->" << w.u << ") <=" << w.min_distance ;
		
		if (d.connected(to)) {
			vec<int> path;
			int u = to;
			path.push(u);
			int p;
			while ((p = d.previous(u)) != -1) {
				Edge & edg = outer->edge_list[d.incomingEdge(u)]; //outer->edges[p][u];
				path.push(p);
				u = p;
			}
			
			for (int i = path.size() - 1; i >= 0; i--) {
				write_to << path[i] << ",";
			}
			write_to << '\n';
		} else {
			write_to << ": FALSE\n";
		}
	}
	
	int width = sqrt(outer->nNodes());
	int maxw = log10(outer->nNodes()) + 1; //highestbit(bits);
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
		
		int y = (n + extra) / width;
		if (y > lasty)
			write_to << "\n";
		
		int d = underapprox_unweighted_distance_detector->distance(n);
		//printf("%*d ",maxw,d);
		write_to << std::setw(maxw) << d;
		
		lasty = y;
	}
	write_to << "\n";
	
}

template<typename Weight>
void DistanceDetector<Weight>::buildDistanceLEQReason(int to, Weight & min_distance, vec<Lit> & conflict, bool strictComparison) {
	stats_distance_leq_reasons++;
	Distance<Weight> & d = *underapprox_weighted_path_detector;
	underapprox_weighted_distance_detector->update();
	Weight & actual_dist = underapprox_weighted_distance_detector->distance(to);
	double starttime = rtime(2);
	d.update();
	stats_under_conflicts++;
	assert(outer->dbg_reachable(d.getSource(), to));
	assert(underapprox_unweighted_distance_detector->connected(to));
	assert(
			underapprox_weighted_distance_detector->distance(to) <= min_distance
					&& underapprox_weighted_distance_detector->distance(to)
							!= underapprox_weighted_distance_detector->unreachable());
	assert(d.connected_unchecked(to));
	if (!d.connected(to) || d.distance(to) > min_distance) {
		fprintf(stderr, "Error in shortest path detector, aborting\n");
		exit(4);
	}

	//the reason that the distance is less than or equal to min_distance is because the shortest path is less than this weight
	{
		int u = to;
		int p;
		while ((p = d.previous(u)) != -1) {
			Edge & edg = outer->edge_list[d.incomingEdge(u)]; //outer->edges[p][u];
			int edgeID = edg.edgeID;
			Var e = edg.v;
			lbool val = outer->value(e);
			assert(outer->value(e)==l_True);
			if(!g_under.isConstant(edgeID))
				conflict.push(mkLit(e, true));
			if(outer->hasBitVector(edgeID)){
	/*			Lit leq;
				Weight w = g_under.getWeight(edgeID);
				if(strictComparison){
					leq= outer->getEdgeWeightGT(edgeID,w);
				}else{
					leq= ~outer->getEdgeWeightLEQ(edgeID,w);
				}
				lbool val = outer->dbg_value(leq);
				assert(val!=l_True);
				conflict.push(leq);*/
				outer->buildBVReason(outer->getEdgeBV(edgeID).getID(),Comparison::leq,g_under.getWeight(edgeID),conflict);
			}

			u = p;
		}
	}
	outer->num_learnt_paths++;
	outer->learnt_path_clause_length += (conflict.size() - 1);
	
	stats_under_clause_length += (conflict.size() - 1);
	double elapsed = rtime(2) - starttime;
	stats_under_conflict_time += elapsed;
	
}

template<typename Weight>
void DistanceDetector<Weight>::buildDistanceGTReason(int to, Weight & min_distance, vec<Lit> & conflict, bool strictComparison) {
	static int it = 0;
	stats_distance_gt_reasons++;
	stats_over_conflicts++;
	++it;
	double starttime = rtime(2);
	int u = to;
	bool reaches = overapprox_weighted_distance_detector->connected(to);
	if (!reaches && opt_conflict_min_cut && conflict_flow) {
		g_over.drawFull();
		cut.clear();
		long f;
		
		assert(conflict_flow->getSource() == source);
		conflict_flow->setSink(to);
		f = conflict_flow->minCut(cut);
		
		assert(f == cut.size()); //because edges are only ever infinity or 1

		for (int i = 0; i < cut.size(); i++) {
			MaxFlowEdge e = cut[i];
			int cut_id = e.id;
			assert(cut_id % 2 == 0);
			int edgeID = cut_id/2;
			if(!g_over.isConstant(edgeID)){
				Lit l = mkLit(outer->getEdgeVar(edgeID), false);
				assert(outer->value(l)==l_False);
				conflict.push(l);
			}
		}
		outer->learnt_cut_clause_length += (conflict.size() - 1);
		
		stats_over_clause_length += (conflict.size() - 1);
		double elapsed = rtime(2) - starttime;
		outer->mctime += elapsed;
		return;
	}
	Weight & actual_dist = overapprox_weighted_distance_detector->distance(to);
	/*for(auto & e:antig.all_edges){
	 if(antig.edgeEnabled(e.id)){
	 printf("nxg.add_edge(%d,%d,weight=",e.from,e.to);
	 std::cout<<outer->getWeight(e.id)<<")\n";
	 }
	 }*/
	bool connected = overapprox_weighted_distance_detector->connected(to);
#ifndef NDEBUG
	Dijkstra<Weight> d(source, g_over);
	Weight & expected = d.distance(to);
	assert(expected == actual_dist);
	assert(!overapprox_weighted_distance_detector->connected(to) || (strictComparison? actual_dist > min_distance: actual_dist >= min_distance));
#endif
	
	{
		
		vec<int>& to_visit = outer->to_visit;
		vec<char>& seen = outer->seen;
		
		to_visit.clear();
		to_visit.push(to);
		seen.clear();
		seen.growTo(outer->nNodes());
		seen[to] = true;
		
		do {
			
			assert(to_visit.size());
			int u = to_visit.last();
			assert(u != source);
			to_visit.pop();
			assert(seen[u]);
			//add all of this node's incoming disabled edges to the cut, and visit any unseen, non-disabled incoming.edges()
			//for (int i = 0; i < outer->inv_adj[u].size(); i++) {
			for(int i = 0;i<g_over.nIncoming(u);i++){
				//int v = outer->inv_adj[u][i].v;
				//int from = outer->inv_adj[u][i].from;
				//assert(outer->inv_adj[u][i].to == u);
				int from =  g_over.incoming(u,i).node;
				int edge_num =  g_over.incoming(u,i).id;

				//Note: the variable has to not only be assigned false, but assigned false earlier in the trail than the reach variable...
				Var edge_enabled = outer->getEdgeVar(edge_num);
				if (from == u) {
					assert(outer->edge_list[edge_num].to == u);
					assert(outer->edge_list[edge_num].from == u);
					continue; //Self loops are allowed, but just make sure nothing got flipped around...
				}
				assert(from != u);
				
				if (has_weighted_shortest_paths_overapprox && reaches) {
					//This is an optional optimization: if we know that even with all possible edges enabled, the shortest path to from + 1 is >= than the current distance to this node, enabling this edge cannot decrease the shortest path,
					//and so we don't need to consider this edge
					Weight current_dist = overapprox_weighted_distance_detector->distance(u);
					Weight over_approx_dist = over_approx_shortest_paths[from] + g_over.getWeight(edge_num);
					if (over_approx_dist >= current_dist) {
						stats_gt_weighted_edges_skipped++;
						//following this edge cannot shorten this path, so skip it.
						continue;
					}
				}
				
				if (outer->value(edge_enabled) == l_False) {
					//note: we know we haven't seen this edge variable before, because we know we haven't visited this node before
					//if we are already planning on visiting the from node, then we don't need to include it in the conflict (is this correct?)
					//if(!seen[from])
					if(!g_over.isConstant(edge_num)){
						conflict.push(mkLit(edge_enabled,false));
					}
				} else{

					if(reaches){
						//if the edge _is_ enabled, and the node _is_ reachable, and the edge weight is symbolic
						//then part of the reason the shortest path is too long is that this edge is not less than its current weight.
						if(!outer->constantWeight(edge_num)){
							Weight w = g_over.getWeight(edge_num);

							Lit gt;
							if(strictComparison){
								gt= ~outer->getEdgeWeightGEQ(edge_num,w);
								//gt= outer->getEdgeWeightLT(edge_num,w);
							}else{
								gt= ~outer->getEdgeWeightGEQ(edge_num,w);
								//gt= outer->getEdgeWeightLEQ(edge_num,w);
								//gt= ~outer->getEdgeWeightGT(edge_num,w);
							}
							lbool val = outer->dbg_value(gt);
							assert(val!=l_True);
							conflict.push(gt);
						}
					}
					//for distance analysis, we _can_ end up reaching source.
					if (from != source) {
						//even if it is undef? probably...
						if (!seen[from]) {
							seen[from] = true;
							to_visit.push(from);
						}
					}
				}
			}
		} while (to_visit.size());
		
	}
	
	outer->learnt_cut_clause_length += (conflict.size() - 1);
	
	stats_over_clause_length += (conflict.size() - 1);
	double elapsed = rtime(2) - starttime;
	outer->mctime += elapsed;
	
}
template<typename Weight>
void DistanceDetector<Weight>::buildReason(Lit p, vec<Lit> & reason, CRef marker) {
	
	if (marker == unweighted_underprop_marker) {
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
		buildUnweightedDistanceLEQReason(u, reason);
		
		//double elapsed = rtime(2)-startpathtime;
		//	pathtime+=elapsed;
	} else if (marker == unweighted_overprop_marker) {
		reason.push(p);
		
		//the reason is a cut separating p from s;
		//We want to find a min-cut in the full graph separating, where activated edges (ie, those still in antig) are weighted infinity, and all others are weighted 1.
		
		//This is a cut that describes a minimal set of edges which are disabled in the current graph, at least one of which would need to be activated in order for s to reach p
		//assign the mincut edge weights if they aren't already assigned.
		
		Var v = var(p);
		int t = getNode(v); // v- var(reach_lits[d][0]);
		int within = getMaximumDistance(v);
		buildUnweightedDistanceGTReason(t,within, reason);
		
	} else {
		exit(3);
		assert(false);
	}
}

template<typename Weight>
void DistanceDetector<Weight>::updateShortestPaths(bool unweighted) {
	if (opt_shortest_path_prune_dist && outer->decisionLevel() == 0) {
		//only update these distances at level 0, to ensure they are a valid over approximate of the shortest possible path to each node
		if (unweighted) {
			has_unweighted_shortest_paths_overapprox = true;
			unweighted_over_approx_shortest_paths.growTo(g_under.nodes());
			for (int i = 0; i < g_over.nodes(); i++) {
				unweighted_over_approx_shortest_paths[i] = overapprox_unweighted_distance_detector->distance(i);
			}
		} else {
			has_weighted_shortest_paths_overapprox = true;
			over_approx_shortest_paths.growTo(g_under.nodes());
			for (int i = 0; i < g_over.nodes(); i++) {
				over_approx_shortest_paths[i] = overapprox_weighted_distance_detector->distance(i);
			}
		}
	}
}
template<typename Weight>
void DistanceDetector<Weight>::preprocess() {
	is_changed.growTo(g_under.nodes());

	if(!overapprox_unweighted_distance_detector){
		buildUnweightedSATConstraints(false);
	}else if (!underapprox_unweighted_distance_detector) {
		//can optimize
		buildUnweightedSATConstraints(true);
/*
		if (unweighted_sat_lits.size() >= within_steps) {
			//then we are using the sat encoding and can directly look up its corresponding distance lit

			Lit r = unweighted_sat_lits[within_steps][to];
			//force equality between the new lit and the old reach lit, in the SAT solver
			outer->makeEqual(r, reachLit);
		}*/

	}


}




template<typename Weight>
bool DistanceDetector<Weight>::propagate(vec<Lit> & conflict) {
	if (!underapprox_unweighted_distance_detector)
		return true;
	
	static int iter = 0;
	if (++iter == 267) { //18303
		int a = 1;
	}

	//printf("iter %d\n",iter);
	bool skipped_positive = false;
	//getChanged().clear();
	if (!opt_detect_pure_theory_lits || unassigned_positives > 0) {
		double startdreachtime = rtime(2);
		stats_under_updates++;
		underapprox_unweighted_distance_detector->update();
		double reachUpdateElapsed = rtime(2) - startdreachtime;
		stats_under_update_time += reachUpdateElapsed;
	} else {
		skipped_positive = true;
		stats_skipped_under_updates++;
		//stats_pure_skipped++;
	}
	
	//  outer->reachupdatetime+=reachUpdateElapsed;
	bool skipped_negative = false;
	
	if (!opt_detect_pure_theory_lits || unassigned_negatives > 0) {
		double startunreachtime = rtime(2);
		stats_over_updates++;
		overapprox_unweighted_distance_detector->update();
		double unreachUpdateElapsed = rtime(2) - startunreachtime;
		stats_over_update_time += unreachUpdateElapsed;
		updateShortestPaths(true); //only needed for the shortest path theory
	} else {
		skipped_negative = true;
		stats_skipped_over_updates++;
	}
	//outer->unreachupdatetime+=unreachUpdateElapsed;
	
	if (opt_rnd_shuffle) {
		randomShuffle(rnd_seed, changed);
	}
	
	while (changed.size()) {
		int sz = changed.size();
		int u = changed.last().u;
		if (u == 3) {
			int a = 1;
		}
		assert(is_changed[u]);
		for (int i = 0; i < unweighted_dist_lits[u].size(); i++) {
			int& min_distance = unweighted_dist_lits[u][i].min_unweighted_distance;
			
			Var v = var(unweighted_dist_lits[u][i].l);
			Lit l;
			
			if (underapprox_unweighted_distance_detector && !skipped_positive
					&& underapprox_unweighted_distance_detector->connected(u)
					&& underapprox_unweighted_distance_detector->distance(u) <= min_distance) {
				l = mkLit(v, false);
			} else if (overapprox_unweighted_distance_detector && !skipped_negative
					&& (!overapprox_unweighted_distance_detector->connected(u)
							|| overapprox_unweighted_distance_detector->distance(u) > min_distance)) {
				l = mkLit(v, true);
			} else {
				continue;
			}
			
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
				
				if (reach)
					outer->enqueue(l, unweighted_underprop_marker);
				else
					outer->enqueue(l, unweighted_overprop_marker);
				
			} else if (outer->value(l) == l_False) {
				conflict.push(l);
				
				if (reach) {
					
					//conflict
					//The reason is a path in g from to s in d
					buildUnweightedDistanceLEQReason(u, conflict);
					//add it to s
					//return it as a conflict
					
				} else {
					//The reason is a cut separating s from t
					buildUnweightedDistanceGTReason(u,getMaximumDistance(v), conflict);
					
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
		assert(sz == changed.size());
		assert(changed.last().u == u);
		is_changed[u] = false;
		changed.pop();
	}
	
	if (opt_rnd_shuffle && weighted_dist_lits.size()) {
		randomShuffle(rnd_seed, weighted_dist_lits);
	}
	if (opt_rnd_shuffle && weighted_dist_bv_lits.size()) {
		randomShuffle(rnd_seed, weighted_dist_bv_lits);
	}
	if (weighted_dist_lits.size() || weighted_dist_bv_lits.size()) {
		updateShortestPaths(false);						//only needed for the shortest path theory
	}
	//now, check for weighted distance lits
	for (auto & dist_lit : weighted_dist_lits) {
		bool strictComparison = dist_lit.strictComparison;
		Lit l = dist_lit.l;
		int to = dist_lit.u;
		Weight & min_dist = dist_lit.min_distance;
		Weight & over_dist = underapprox_weighted_distance_detector->distance(to);
		Weight & under_dist = overapprox_weighted_distance_detector->distance(to);
		if (underapprox_weighted_distance_detector->connected(to)
				&& (!strictComparison ? (underapprox_weighted_distance_detector->distance(to) <= min_dist) : (underapprox_weighted_distance_detector->distance(to) < min_dist))) {
			if (outer->value(l) == l_True) {
				//do nothing
			} else if (outer->value(l) == l_Undef) {
				outer->enqueue(l, weighted_underprop_marker);
			} else if (outer->value(l) == l_False) {
				//conflict
				conflict.push(l);

				buildDistanceLEQReason(to, min_dist, conflict,strictComparison);

				return false;
			}
		}
		if (!overapprox_weighted_distance_detector->connected(to)
				|| (!strictComparison ? (overapprox_weighted_distance_detector->distance(to) > min_dist) :  (overapprox_weighted_distance_detector->distance(to) >= min_dist))) {
			if (outer->value(~l) == l_True) {
				//do nothing
			} else if (outer->value(~l) == l_Undef) {
				outer->enqueue(~l, weighted_overprop_marker);
			} else if (outer->value(~l) == l_False) {
				//conflict
				conflict.push(~l);

				buildDistanceGTReason(to, min_dist, conflict,strictComparison);

				return false;
			}
		}
		
	}
	//g_under.drawFull(true);
	//g_over.drawFull(true);
	for (auto & dist_lit : weighted_dist_bv_lits) {
		Lit l = dist_lit.l;
		int to = dist_lit.u;
		bool strictComparison = dist_lit.strictComparison;
		BitVector<Weight> & bv = dist_lit.bv;
		Weight & min_dist_under = bv.getUnder();
		Weight & min_dist_over = bv.getOver();
		Weight & under_dist = underapprox_weighted_distance_detector->distance(to);
		Weight & over_dist = overapprox_weighted_distance_detector->distance(to);


		if (underapprox_weighted_distance_detector->connected(to)
				&& (strictComparison? (under_dist < min_dist_under) :  (under_dist <= min_dist_under))) {
			if (outer->value(l) == l_True) {
				//do nothing
			} else if (outer->value(l) == l_Undef) {
				outer->enqueue(l, weighted_underprop_marker);
			} else if (outer->value(l) == l_False) {
				//conflict

				conflict.push(l);
				//conflict.push(~outer->getBV_LEQ(bv.getID(),min_dist_under));
				if(strictComparison){
					//std::cout<<"distance conflict " << under_dist << " < " << bv.getID() << "\n";
					/*Lit d = outer->getBV_LEQ(bv.getID(),under_dist);
					conflict.push(d);
					lbool val = outer->value(d);
					lbool dbgval = outer->dbg_value(d);
					assert(dbgval!=l_True);*/
					outer->buildBVReason(bv.getID(),Comparison::gt,under_dist,conflict);
					buildDistanceLEQReason(to, min_dist_under, conflict,true);
				}else{
					//std::cout<<"distance conflict " << under_dist << " <= " << bv.getID() << "\n";
/*					conflict.push(outer->getBV_LT(bv.getID(),under_dist));
					lbool val = outer->value(outer->getBV_LT(bv.getID(),under_dist));
					lbool dbgval = outer->dbg_value(outer->getBV_LT(bv.getID(),under_dist));
					assert(dbgval!=l_True);*/
					outer->buildBVReason(bv.getID(),Comparison::geq,under_dist,conflict);
					buildDistanceLEQReason(to, min_dist_under, conflict,false);
				}

				return false;
			}
		}
		if (!overapprox_weighted_distance_detector->connected(to)
				|| (strictComparison? (over_dist >= min_dist_over) :  (over_dist > min_dist_over))){//over_dist > min_dist_over) {
			if (outer->value(~l) == l_True) {
				//do nothing
			} else if (outer->value(~l) == l_Undef) {
				outer->enqueue(~l, weighted_overprop_marker);
			} else if (outer->value(~l) == l_False) {
				//conflict
				conflict.push(~l);
				if(strictComparison){
					//std::cout<<"distance conflict " << over_dist << " > " << bv.getID() << "\n";
					//conflict.push(outer->getBV_LT(bv.getID(),min_dist_over));
					if(overapprox_weighted_distance_detector->connected(to)){
			/*			conflict.push(outer->getBV_GT(bv.getID(),over_dist));
						lbool val = outer->value(outer->getBV_GEQ(bv.getID(),over_dist));
						lbool dbgval = outer->dbg_value(outer->getBV_GEQ(bv.getID(),over_dist));
						assert(dbgval!=l_True);*/
						outer->buildBVReason(bv.getID(),Comparison::leq,over_dist,conflict);
					}
					buildDistanceGTReason(to, min_dist_over, conflict,false);
				}else{
					//std::cout<<"distance conflict " << over_dist << " >= " << bv.getID() << "\n";
					//conflict.push(outer->getBV_LT(bv.getID(),min_dist_over));
					if(overapprox_weighted_distance_detector->connected(to)){
		/*				conflict.push(outer->getBV_GEQ(bv.getID(),over_dist));
						lbool val = outer->value(outer->getBV_GEQ(bv.getID(),over_dist));
						lbool dbgval = outer->dbg_value(outer->getBV_GEQ(bv.getID(),over_dist));
						assert(dbgval!=l_True);*/
						outer->buildBVReason(bv.getID(),Comparison::lt,over_dist,conflict);
					}
					buildDistanceGTReason(to, min_dist_over, conflict,true);
				}
				return false;
			}
		}

	}
	return true;
}
template<typename Weight>
bool DistanceDetector<Weight>::checkSatisfied() {
	UnweightedDijkstra<Weight> under(source, g_under);
	UnweightedDijkstra<Weight> over(source, g_over);
	under.update();
	over.update();
	for (int j = 0; j < unweighted_dist_lits.size(); j++) {
		for (int k = 0; k < unweighted_dist_lits[j].size(); k++) {
			Lit l = unweighted_dist_lits[j][k].l;
			int dist = unweighted_dist_lits[j][k].min_unweighted_distance;
			
			if (l != lit_Undef) {
				int node = getNode(var(l));
				
				if (outer->value(l) == l_True) {
					if (!under.connected(node) || under.distance(node) > dist) {
						return false;
					}
				} else if (outer->value(l) == l_False) {
					if (under.connected(node) && over.distance(node) <= dist) {
						return false;
					}
				} else {
					if (under.connected(node) && under.distance(node) <= dist) {
						return false;
					}
					if (!over.connected(node) || over.distance(node) > dist) {
						return false;
					}
				}
			}
		}
	}

	{
		Dijkstra<Weight> under(source, g_under);
		Dijkstra<Weight> over(source, g_over);
		g_under.drawFull(true);
		g_over.drawFull(true);
		under.update();
		over.update();
		//now, check for weighted distance lits
		for (auto & dist_lit : weighted_dist_lits) {
			Lit l = dist_lit.l;
			int to = dist_lit.u;
			Weight & min_dist = dist_lit.min_distance;
			Weight & over_dist = under.distance(to);
			Weight & under_dist = over.distance(to);
			if (under.connected(to) &&  under.distance(to) <= min_dist) {
				if (outer->value(l) == l_True) {
					//do nothing
				} else if (outer->value(l) == l_Undef) {
					outer->enqueue(l, weighted_underprop_marker);
				} else if (outer->value(l) == l_False) {
					return false;
				}
			}
			if (!over.connected(to) || over.distance(to) > min_dist) {
				if (outer->value(~l) == l_True) {
					//do nothing
				} else if (outer->value(~l) == l_Undef) {
					outer->enqueue(~l, weighted_overprop_marker);
				} else if (outer->value(~l) == l_False) {
					
					return false;
				}
			}
		}
		
	}
	
	//	}
	return true;
}

/*
 int DistanceDetector<Weight>::OptimalWeightEdgeStatus::operator [] (int edge) const {
 Var v = detector.outer->edge_list[edge].v;
 lbool val = detector.outer->value(v);
 if(val==l_False){
 assert(false);
 return detector.outer->edge_list.size()*2;
 }else if (val==l_True){
 return 0;
 }else{
 assert(val==l_Undef);
 return 1;
 }
 }

 int DistanceDetector<Weight>::OptimalWeightEdgeStatus::size()const{
 return detector.outer->edge_list.size();
 }
 */

template<typename Weight>
Lit DistanceDetector<Weight>::decide() {
	if (!opt_decide_graph_distance || !overapprox_unweighted_distance_detector)
		return lit_Undef;
	DistanceDetector *r = this;
	auto * over = r->overapprox_unweighted_distance_detector;
	
	auto * under = r->underapprox_unweighted_distance_detector;
	
	//we can probably also do something similar, but with cuts, for nodes that are decided to be unreachable.
	
	//ok, for each node that is assigned reachable, but that is not actually reachable in the under approx, decide an edge on a feasible path
	
	//this can be obviously more efficient
	//for(int j = 0;j<nNodes();j++){
	/*	if(opt_decide_graph_neg){

	 }*/

	for (int k = 0; k < unweighted_dist_lits.size(); k++) {
		for (int n = 0; n < unweighted_dist_lits[k].size(); n++) {
			Lit l = unweighted_dist_lits[k][n].l;
			int min_dist = unweighted_dist_lits[k][n].min_unweighted_distance;
			if (l == lit_Undef)
				continue;
			int j = r->getNode(var(l));
			if (outer->value(l) == l_True) {
				if (opt_decide_graph_pos) {
					//if(S->level(var(l))>0)
					//	continue;
					
					assert(over->distance(j) <= min_dist);//else we would already be in conflict before this decision was attempted!
					if (under->distance(j) > min_dist) {
						//then lets try to connect this
						//static vec<bool> print_path;
						
						assert(over->connected(j));					//Else, we would already be in conflict
								
						int p = j;
						int last = j;
						int last_edge = -1;
						if (!opt_use_random_path_for_decisions) {
							//ok, read back the path from the over to find a candidate edge we can decide
							//find the earliest unconnected node on this path
							over->update();
							p = j;
							last = j;
							int dist = 0;
							while (under->distance(p) >= min_dist - dist) {
								
								last = p;
								assert(p != r->source);
								int prev = over->previous(p);
								last_edge = over->incomingEdge(p);
								assert(over->distance(p) <= min_dist - dist);
								dist += 1;						//should really be weighted
								p = prev;
								
							}
						} else {
							//This won't work (without modification) because we need to constrain these paths to ones of maximum real distance < min_dist.
							//Randomly re-weight the graph sometimes
							if (drand(rnd_seed) < opt_decide_graph_re_rnd) {
								
								for (int i = 0; i < outer->edge_list.size(); i++) {
									double w = drand(rnd_seed);
									/* w-=0.5;
									 w*=w;*/
									//printf("%f (%f),",w,rnd_seed);
									//rnd_path->setWeight(i,w);
									rnd_weight[i] = w;
								}
							}
							rnd_path->update();
							//derive a random path in the graph
							p = j;
							last = j;
							assert(rnd_path->connected(p));
							while (!under->connected(p)) {
								
								last = p;
								
								assert(p != source);
								int prev = rnd_path->previous(p);
								last_edge = rnd_path->incomingEdge(p);
								p = prev;
								assert(p >= 0);
							}
							
						}
						
						//ok, now pick some edge p->last that will connect p to last;
						/*	assert(!under->connected(last));
						 assert(under->connected(p));

						 assert(over->connected(last));
						 assert(over->connected(p));*/
						assert(p > -1);
						if (p > -1) {
							assert(outer->edge_list[last_edge].from == p);
							assert(outer->edge_list[last_edge].to == last);
							Var v = outer->getEdgeVar(last_edge); ///outer->edges[p][last].v;
							if (outer->value(v) == l_Undef) {
								return mkLit(v, false);
							} else {
								assert(outer->value(v)!=l_True);
							}
						}
						/*		for(int k = 0;k<outer->antig.adjacency[p].size();k++){
						 int to = outer->antig.adjacency[p][k].node;
						 if (to==last){
						 Var v =outer->edge_list[ outer->antig.adjacency[p][k].id].v;
						 if(outer->value(v)==l_Undef){
						 return mkLit(v,false);
						 }else{
						 assert(outer->value(v)!=l_True);
						 }
						 }
						 }*/

					}
				}
			} else if (outer->value(l) == l_False) {
				if (opt_decide_graph_neg) {
					
					//assert(over->distance(j)<=min_dist);//else we would already be in conflict before this decision was attempted!
					
					if (over->distance(j) <= min_dist && under->distance(j) > min_dist) {
						//then lets try to disconnect this node from source by walking back along the path in the over approx, and disabling the first unassigned edge we see.
						//(there must be at least one such edge, else the variable would be connected in the under approximation as well - in which case it would already have been propagated.
						int p = j;
						int last = j;
						int dist = 0;
						// tmp_nodes.clear();
						while (under->distance(p) > min_dist - dist) {
							//int d = under->distance(p);
							//int d_over = over->distance(p);
							last = p;
							assert(p != source);
							int prev = over->previous(p);
							int edgeid = over->incomingEdge(p);
							assert(outer->edge_list[edgeid].from == prev);
							assert(outer->edge_list[edgeid].to == p);
							Var v = outer->getEdgeVar(edgeid); //outer->edges[prev][p].v;
							if (outer->value(v) == l_Undef) {
								//if(opt_use_random_path_for_decisions)
								//	tmp_nodes.push(v);
								//else
								return mkLit(v, true);
							} else {
								assert(outer->value(v)!=l_False);
							}
							assert(over->distance(p) <= min_dist - dist);
							dist += 1;								//should really be weighted
							p = prev;
							
						}
						/*		assert(opt_use_random_path_for_decisions);
						 assert(tmp_nodes.size()>0);
						 int i = irand(rnd_seed,tmp_nodes.size());
						 Var v= tmp_nodes[i];
						 return mkLit(v,true);*/

					}
				}
			}
		}
	}
	return lit_Undef;
}
;

template class Monosat::DistanceDetector<int> ;
template class Monosat::DistanceDetector<long> ;
template class Monosat::DistanceDetector<double> ;
#include <gmpxx.h>
template class Monosat::DistanceDetector<mpq_class> ;
