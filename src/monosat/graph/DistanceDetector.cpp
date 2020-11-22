/**************************************************************************************************
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

#include "monosat/core/Config.h"
#include "monosat/dgl/RamalReps.h"
#include "monosat/dgl/RamalRepsBatched.h"
#include "monosat/dgl/RamalRepsBatchedUnified.h"

#include "monosat/dgl/EdmondsKarp.h"
#include "monosat/dgl/EdmondsKarpAdj.h"
#include "monosat/dgl/KohliTorr.h"
#include "monosat/dgl/EdmondsKarpDynamic.h"
#include "monosat/dgl/Dinics.h"
#include "monosat/dgl/DinicsLinkCut.h"
#include "monosat/graph/DistanceDetector.h"
#include "monosat/graph/GraphTheory.h"
#include "monosat/mtl/Rnd.h"
#include "monosat/mtl/Vec.h"
#include "monosat/utils/Options.h"
#include <monosat/graph/GraphHeuristic.h>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <iomanip>

using namespace Monosat;

template<typename Weight, typename Graph>
DistanceDetector<Weight, Graph>::DistanceDetector(int _detectorID, GraphTheorySolver<Weight>* outer,
                                                  Graph& g_under, Graph& g_over, Graph& cutGraph, int from, double seed)
        :
        Detector(_detectorID), outer(outer), g_under(g_under), g_over(g_over), cutGraph(cutGraph), source(from),
        rnd_seed(seed){
    max_unweighted_distance = -1;
    rnd_path = NULL;

    constraintsBuilt = -1;
    first_reach_var = var_Undef;
    stats_pure_skipped = 0;
    if(distalg == DistAlg::ALG_SAT){
        positiveReachStatus = nullptr;
        negativeReachStatus = nullptr;

        underapprox_unweighted_distance_detector = nullptr;
        overapprox_unweighted_distance_detector = nullptr;
        underapprox_path_detector = nullptr;

        //we are just going to directly enforce these constraints in the original SAT solver, so _nothing_ will end up happening in this detector (except for creating the clauses needed to enforce these constraints).

        return;
    }

    if(opt_use_random_path_for_decisions){
        rnd_weight.clear();
        rnd_path = new WeightedDijkstra<Weight, Graph, double>(from, g_over, rnd_weight);
        for(int i = 0; i < outer->nEdges(); i++){
            double w = drand(rnd_seed);
            rnd_weight.push_back(w);
        }

    }
    positiveReachStatus = new DistanceDetector<Weight, Graph>::ReachStatus(*this, true);
    negativeReachStatus = new DistanceDetector<Weight, Graph>::ReachStatus(*this, false);

    //select the unweighted distance detectors
    if(distalg == DistAlg::ALG_DISTANCE){
        if(outer->all_edges_unit){
            if(!opt_encode_dist_underapprox_as_sat)
                underapprox_unweighted_distance_detector = new UnweightedBFS<Weight, Graph, typename DistanceDetector<Weight, Graph>::ReachStatus>(
                        from,
                        g_under, *(positiveReachStatus), 0);
            overapprox_unweighted_distance_detector = new UnweightedBFS<Weight, Graph, typename DistanceDetector<Weight, Graph>::ReachStatus>(
                    from,
                    g_over, *(negativeReachStatus), 0);
            if(underapprox_unweighted_distance_detector)
                underapprox_path_detector = underapprox_unweighted_distance_detector;
            else{
                underapprox_path_detector = new UnweightedBFS<Weight, Graph, typename DistanceDetector<Weight, Graph>::ReachStatus>(
                        from,
                        g_under, *(positiveReachStatus), 0);
            }

        }else{
            if(!opt_encode_dist_underapprox_as_sat){
                underapprox_unweighted_distance_detector = new UnweightedDijkstra<Weight, Graph, typename DistanceDetector<Weight, Graph>::ReachStatus>(
                        from, g_under, *positiveReachStatus, 0);
            }
            overapprox_unweighted_distance_detector = new UnweightedDijkstra<Weight, Graph, typename DistanceDetector<Weight, Graph>::ReachStatus>(
                    from, g_over, *negativeReachStatus, 0);
            underapprox_path_detector = new UnweightedBFS<Weight, Graph, Distance<int>::NullStatus>(from, g_under,
                                                                                                    Distance<int>::nullStatus,
                                                                                                    0);
        }

    }else if(distalg == DistAlg::ALG_RAMAL_REPS){
        if(!opt_encode_dist_underapprox_as_sat){
            underapprox_unweighted_distance_detector = new UnweightedRamalReps<Weight, Graph,
                    typename DistanceDetector<Weight, Graph>::ReachStatus>(from, g_under, *(positiveReachStatus), 0);
        }
        overapprox_unweighted_distance_detector =
                new UnweightedRamalReps<Weight, Graph, typename DistanceDetector<Weight, Graph>::ReachStatus>(from,
                                                                                                              g_over,
                                                                                                              *(negativeReachStatus),
                                                                                                              0);

        if(underapprox_unweighted_distance_detector){
            underapprox_path_detector = underapprox_unweighted_distance_detector;
        }else{
            underapprox_path_detector = underapprox_unweighted_distance_detector = new UnweightedRamalReps<Weight, Graph,
                    typename DistanceDetector<Weight, Graph>::ReachStatus>(from, g_under, *(positiveReachStatus), 0);
        }
    }else if(distalg == DistAlg::ALG_RAMAL_REPS_BATCHED){
        if(!opt_encode_dist_underapprox_as_sat){
            underapprox_unweighted_distance_detector = new UnweightedRamalRepsBatched<Weight, Graph,
                    typename DistanceDetector<Weight, Graph>::ReachStatus>(from, g_under, *(positiveReachStatus), 0);
        }
        overapprox_unweighted_distance_detector =
                new UnweightedRamalRepsBatched<Weight, Graph, typename DistanceDetector<Weight, Graph>::ReachStatus>(
                        from, g_over,
                        *(negativeReachStatus), 0);

        if(underapprox_unweighted_distance_detector){
            underapprox_path_detector = underapprox_unweighted_distance_detector;
        }else{
            underapprox_path_detector = underapprox_unweighted_distance_detector = new UnweightedRamalRepsBatched<Weight, Graph,
                    typename DistanceDetector<Weight, Graph>::ReachStatus>(from, g_under, *(positiveReachStatus), 0);
        }
    }else if(distalg == DistAlg::ALG_RAMAL_REPS_BATCHED2){
        if(!opt_encode_dist_underapprox_as_sat){
            underapprox_unweighted_distance_detector = new UnweightedRamalRepsBatchedUnified<Weight, Graph,
                    typename DistanceDetector<Weight, Graph>::ReachStatus>(from, g_under, *(positiveReachStatus), 0);
        }
        overapprox_unweighted_distance_detector =
                new UnweightedRamalRepsBatchedUnified<Weight, Graph, typename DistanceDetector<Weight, Graph>::ReachStatus>(
                        from, g_over,
                        *(negativeReachStatus), 0);

        if(underapprox_unweighted_distance_detector){
            underapprox_path_detector = underapprox_unweighted_distance_detector;
        }else{
            underapprox_path_detector = underapprox_unweighted_distance_detector = new UnweightedRamalRepsBatchedUnified<Weight, Graph,
                    typename DistanceDetector<Weight, Graph>::ReachStatus>(from, g_under, *(positiveReachStatus), 0);
        }
    }else{
        if(!opt_encode_dist_underapprox_as_sat){
            underapprox_unweighted_distance_detector =
                    new UnweightedDijkstra<Weight, Graph, typename DistanceDetector<Weight, Graph>::ReachStatus>(from,
                                                                                                                 g_under,
                                                                                                                 *positiveReachStatus,
                                                                                                                 0);
        }
        overapprox_unweighted_distance_detector =
                new UnweightedDijkstra<Weight, Graph, typename DistanceDetector<Weight, Graph>::ReachStatus>(from,
                                                                                                             g_over,
                                                                                                             *negativeReachStatus,
                                                                                                             0);

        if(underapprox_unweighted_distance_detector)
            underapprox_path_detector = underapprox_unweighted_distance_detector;
        else{
            underapprox_path_detector = new UnweightedDijkstra<Weight, Graph, typename DistanceDetector<Weight, Graph>::ReachStatus>(
                    from, g_under, *positiveReachStatus, 0);
        }
    }


    if(opt_conflict_min_cut){
        if(mincutalg == MinCutAlg::ALG_EDKARP_DYN){
            conflict_flow = new EdmondsKarpDynamic<Weight>(cutGraph, source, 0);
        }else if(mincutalg == MinCutAlg::ALG_EDKARP_ADJ){
            conflict_flow = new EdmondsKarpAdj<Weight>(cutGraph, source, 0);
        }else if(mincutalg == MinCutAlg::ALG_DINITZ){
            conflict_flow = new Dinitz<Weight>(cutGraph, source, 0);
        }else if(mincutalg == MinCutAlg::ALG_DINITZ_LINKCUT){
            //link-cut tree currently only supports ints
            conflict_flow = new Dinitz<Weight>(cutGraph, source, 0);

        }else if(mincutalg == MinCutAlg::ALG_KOHLI_TORR){
            if(opt_use_kt_for_conflicts){
                conflict_flow = new KohliTorr<Weight>(cutGraph, source, 0,
                                                      opt_kt_preserve_order);
            }else
                conflict_flow = new EdmondsKarpDynamic<Weight>(cutGraph, source, 0);
        }else{
            conflict_flow = new EdmondsKarpAdj<Weight>(cutGraph, source, 0);
        }
    }

    unweighted_underprop_marker = outer->newReasonMarker(getID());
    unweighted_overprop_marker = outer->newReasonMarker(getID());

    default_heuristic = new GraphHeuristic<Weight>(outer, this);


}

template<typename Weight, typename Graph>
void DistanceDetector<Weight, Graph>::buildUnweightedSATConstraints(bool onlyUnderApprox, int within_steps){
    if(within_steps < 0)
        within_steps = g_under.nodes();
    if(within_steps > g_under.nodes())
        within_steps = g_under.nodes();
    if(within_steps > g_under.edges())
        within_steps = g_under.edges();
    if(constraintsBuilt >= within_steps)
        return;

    if(onlyUnderApprox && opt_encode_dist_underapprox_as_sat == 2){
        //this is only correct for directed graphs.
        //To do the same for undirected graps, also non-deterministically set the direction of each edge.
        Lit True = mkLit(outer->newVar());
        outer->addClause(True);
        BitVector<Weight> one = outer->bvTheory->newBitvector(-1, outer->getEdgeWeightBitWidth(), 1);
        for(int to = 0; to < unweighted_dist_lits.size(); to++){
            if(unweighted_dist_lits[to].size() && to != source){

                //For each target node, create an Erez-Nadel style nondeterministic graph
                //then non-deterministically pick an acyclic path
                //and assert that it's length is less than the target length
                vec<Lit> induced_graph_edges;
                induced_graph_edges.growTo(g_under.edges(), lit_Undef);
                for(int edgeID = 0; edgeID < g_under.edges(); edgeID++){
                    if(g_over.hasEdge(edgeID)){//&& g_over.edgeEnabled(edgeID)
                        int from = g_under.getEdge(edgeID).from;
                        int to = g_under.getEdge(edgeID).to;
                        Lit edge_enabled = mkLit(outer->getEdgeVar(edgeID), false);
                        Lit induced_edge_enabled = mkLit(outer->newVar(), false);
                        induced_graph_edges[edgeID] = induced_edge_enabled;
                        //If the edge is disabled, then our induced graph edge must also be disabled (but _not_ vice versa).
                        outer->addClause(edge_enabled, ~induced_edge_enabled);
                    }
                }

                vec<BitVector<Weight>> node_distances;
                //introduce a bitvector distance for each node
                for(int n = 0; n < g_under.nodes(); n++){
                    if(n == source){
                        node_distances.push(outer->newBV(0));//source has constant 0 distance
                    }else{
                        BitVector<Weight> bv = outer->newBV();
                        node_distances.push(bv);
                    }
                }

                //now, non-deterministically select an acyclic path from source out of this graph, asserting that it ends at the target node
                Lit reaches_to = lit_Undef;
                //every vertex (other than source or to) has at most one in-coming edge enabled, and if it has such an edge, then it has exactly one outgoing edge enabled (else no outgoing edges)
                //source has no outgoing edge, and exactly one outgoing edge
                for(int n = 0; n < g_under.nodes(); n++){
                    {
                        Lit any_incoming_enabled = ~True;
                        BitVector<Weight> dist = node_distances[n];
                        if(n != source){
                            BitVector<Weight> dist = opt_sat_distance_encoding_unconstrained_default ? outer->newBV()
                                                                                                     : outer->newBV(
                                            0);//should this be unconstrained, or zero?
                            for(int i = 0; i < g_over.nIncoming(n); i++){
                                int edgeID = g_over.incoming(n, i).id;
                                int from = g_over.incoming(n, i).node;
                                if(from == n)
                                    continue;//no need to consider self loops
                                Lit l = induced_graph_edges[edgeID];
                                //If any prior incoming edge is enabled, this one must be disabled
                                outer->addClause(~any_incoming_enabled, ~l);

                                //if any prior incoming edge was enabled, or l is enabled, then next_incoming is enabled
                                Lit next_incoming = mkLit(outer->newVar(), false);

                                outer->addClause(~any_incoming_enabled, next_incoming);
                                outer->addClause(~l, next_incoming);
                                //if both incoming and l are false, then next incoming is false
                                outer->addClause(l, any_incoming_enabled, ~next_incoming);
                                any_incoming_enabled = next_incoming;

                                BitVector<Weight> sum = outer->newBV();

                                //this is an unweighted constraint, so just add 1
                                outer->bvTheory->newAdditionBV(sum.getID(), node_distances[from].getID(), one.getID());
                                BitVector<Weight> new_dist = outer->newBV();
                                //If this incoming edge is enabled, then the distance is the cost of this edge
                                outer->bvTheory->newConditionalBV(l, new_dist.getID(), sum.getID(), dist.getID());

                                dist = new_dist;
                            }
                            outer->bvTheory->makeEquivalent(dist.getID(), node_distances[to].getID());
                        }else{
                            //all incoming edges must be disabled
                            for(int i = 0; i < g_over.nIncoming(n); i++){
                                int edgeID = g_over.incoming(n, i).id;
                                outer->addClause(~induced_graph_edges[edgeID]);
                            }
                        }
                        //any incoming enabled is true if any incoming edges are enabled
                        //also, at most one incoming edge can possibly be enabled.
                        if(n == to){
                            reaches_to = any_incoming_enabled;
                        }
                        Lit any_outgoing_enabled = ~True;
                        //exactly one outgoing edge can be enabled, but _only_ if any incoming edges are enabled.
                        if(n != to){
                            for(int i = 0; i < g_over.nIncident(n); i++){
                                int edgeID = g_over.incident(n, i).id;
                                Lit l = induced_graph_edges[edgeID];

                                //if no incoming edge was enabled, then all outgoing edges must be disabled
                                outer->addClause(any_incoming_enabled, ~l);

                                //If any prior incoming edge is enabled, this one must be disabled
                                outer->addClause(~any_outgoing_enabled, ~l);

                                //if any prior incoming edge was enabled, or l is enabled, then next_incoming is enabled
                                Lit next_outgoing = mkLit(outer->newVar(), false);

                                outer->addClause(~any_outgoing_enabled, next_outgoing);
                                outer->addClause(~l, next_outgoing);
                                //if both incoming nor l are false, then next incoming is false
                                outer->addClause(l, any_outgoing_enabled, ~next_outgoing);
                                any_outgoing_enabled = next_outgoing;
                            }
                        }else{
                            //all outgoing edges must be disabled
                            for(int i = 0; i < g_over.nIncident(n); i++){
                                int edgeID = g_over.incident(n, i).id;
                                outer->addClause(~induced_graph_edges[edgeID]);
                            }
                        }
                        if(n == source){
                            outer->addClause(any_outgoing_enabled);
                        }else if(n == to){
                            outer->addClause(any_incoming_enabled);
                        }else{
                            //if and only if any incoming edge is enabled, then exactly one outgoing edge must be enabled
                            outer->addClause(~any_incoming_enabled, any_outgoing_enabled);
                            outer->addClause(any_incoming_enabled, ~any_outgoing_enabled);
                        }
                    }
                }
                //for each distance, assert that it's lit is true iff the distance on the path is <= 'distance'
                //following Erez and Nadel 2015's non-graph aware encoding, use bitvectors
                //Var lit = unweighted_dist_lits[to].l;
                //int distance = unweighted_dist_lits[to].min_unweighted_distance;

                BitVector<Weight> dist = node_distances[to];
                for(int i = 0; i < unweighted_dist_lits[to].size(); i++){
                    Lit l = unweighted_dist_lits[to][i].l;
                    int min_unweighted_distance = unweighted_dist_lits[to][i].min_unweighted_distance;
                    BitVector<Weight> comparison = outer->bvTheory->newBitvector(-1, outer->getEdgeWeightBitWidth(),
                                                                                 min_unweighted_distance);
                    Lit conditional = outer->bvTheory->newComparison(Comparison::leq, dist.getID(), comparison.getID());

                    Lit s_l = outer->toSolver(l);
                    Lit s_conditional = outer->bvTheory->toSolver(conditional);
                    Lit s_reaches = outer->toSolver(reaches_to);

                    //l is true if the distance is <= comparison, AND their exists a path to node
                    outer->addClauseToSolver(~s_l, s_conditional);
                    outer->addClauseToSolver(~s_l, s_reaches);
                    outer->addClauseToSolver(s_l, ~s_conditional, ~s_reaches);
                }
            }
        }
    }else{
        assert(outer->decisionLevel() == 0);
        vec<Lit> c;

        if(constraintsBuilt <= 0){
            constraintsBuilt = 0;
            unweighted_dist_lits.push();
            Lit True = mkLit(outer->newVar());
            outer->addClause(True);
            assert(outer->value(True) == l_True);
            Lit False = ~True;
            for(int i = 0; i < g_under.nodes(); i++){
                unweighted_sat_lits[0].push(False);
            }
            unweighted_sat_lits[0][source] = True;
        }

        vec<Lit> reaches;

        //bellman-ford:
        for(int i = constraintsBuilt; i < within_steps; i++){
            unweighted_sat_lits.last().copyTo(reaches);
            unweighted_sat_lits.push();
            reaches.copyTo(unweighted_sat_lits.last());
            assert(outer->value(reaches[source]) == l_True);
            //For each edge:
            for(int j = 0; j < g_under.nodes(); j++){
                Lit r_cur = reaches[j];
                int to = j;

                for(int n = 0; n < g_over.nIncoming(j); n++){
                    int from = g_over.incoming(j, n).node;
                    int v = outer->getEdgeVar(g_over.incoming(j, n).id);
                    if(outer->value(unweighted_sat_lits.last()[to]) == l_True){
                        //do nothing
                    }else if(outer->value(reaches[from]) == l_False){
                        //do nothing
                    }else{
                        Lit l = mkLit(v, false);
                        Lit r = mkLit(outer->newVar(), false);

                        c.clear();
                        c.push(~r);
                        c.push(reaches[to]);
                        c.push(l); //r -> (e.l or reaches[e.to])
                        outer->addClause(c);
                        c.clear();
                        c.push(~r);
                        c.push(reaches[to]);
                        c.push(reaches[from]); //r -> (reaches[e.from]) or reaches[e.to])
                        outer->addClause(c);
                        c.clear();
                        c.push(r);
                        c.push(~reaches[to]); //~r -> ~reaches[e.to]
                        outer->addClause(c);
                        c.clear();
                        c.push(r);
                        c.push(~reaches[from]);
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

template<typename Weight, typename Graph>
void
DistanceDetector<Weight, Graph>::addUnweightedShortestPathLit(int from, int to, Var outer_reach_var, int within_steps){
    g_under.invalidate();
    g_over.invalidate();
    Var reach_var = outer->newVar(outer_reach_var, getID());

    if(first_reach_var == var_Undef){
        first_reach_var = reach_var;
    }else{
        assert(reach_var > first_reach_var);
    }
    assert(from == source);
    assert(within_steps >= 0);
    if(within_steps > g_under.nodes())
        within_steps = g_under.nodes();
    if(within_steps < 0)
        within_steps = g_under.nodes();

    if(within_steps > max_unweighted_distance){
        max_unweighted_distance = within_steps;
        if(opt_compute_max_distance){
            underapprox_unweighted_distance_detector->setMaxDistance(max_unweighted_distance);
            overapprox_unweighted_distance_detector->setMaxDistance(max_unweighted_distance);
        }
    }

    Lit reachLit = mkLit(reach_var, false);


    while(unweighted_dist_lits.size() <= to)
        unweighted_dist_lits.push();

    bool found = false;
    for(int i = 0; i < unweighted_dist_lits[to].size(); i++){
        if(unweighted_dist_lits[to][i].min_unweighted_distance == within_steps){
            found = true;
            Lit r = unweighted_dist_lits[to][i].l;
            //force equality between the new lit and the old reach lit, in the SAT solver
            outer->makeEqual(reachLit, r, true);
        }
    }

    if(!found){
        unweighted_dist_lits[to].push();
        unweighted_dist_lits[to].last().l = reachLit;
        unweighted_dist_lits[to].last().min_unweighted_distance = within_steps;

        while(reach_lit_map.size() <= reach_var - first_reach_var){
            reach_lit_map.push({-1, -1});
        }
        reach_lit_map[reach_var - first_reach_var] = {to, within_steps};
    }

}

template<typename Weight, typename Graph>
void DistanceDetector<Weight, Graph>::ReachStatus::setReachable(int u, bool reachable){

}

template<typename Weight, typename Graph>
void DistanceDetector<Weight, Graph>::ReachStatus::setMininumDistance(int u, bool reachable, int distance){
    assert(reachable == (distance < detector.g_under.nodes()));
    if(distance <= detector.g_under.nodes()){
        setReachable(u, reachable);
    }

    if(u < detector.unweighted_dist_lits.size() && detector.unweighted_dist_lits[u].size()){
        assert(distance >= 0);

        if(!polarity && !detector.is_changed_over[u]){
            detector.is_changed_over[u] = true;
            detector.changed.push({u, polarity});
        }else if(polarity &&
                 !detector.is_changed_under[u]){
            detector.is_changed_under[u] = true;
            detector.changed.push({u, polarity});
        }
    }

}

template<typename Weight, typename Graph>
void DistanceDetector<Weight, Graph>::DistanceStatus::setReachable(int u, bool reachable){

}

template<typename Weight, typename Graph>
void DistanceDetector<Weight, Graph>::DistanceStatus::setMininumDistance(int u, bool reachable, Weight& distance){

}

template<typename Weight, typename Graph>
void DistanceDetector<Weight, Graph>::buildUnweightedDistanceLEQReason(int node, vec<Lit>& conflict){
    Reach& d = *underapprox_path_detector;
    stats_unweighted_leq_reasons++;
    stats_under_conflicts++;
    double starttime = rtime(2);
    d.update();
    assert(outer->dbg_reachable(d.getSource(), node));

    assert(underapprox_unweighted_distance_detector->connected(node));

    assert(d.connected_unchecked(node));
    {
        int u = node;
        int p;
        while((p = d.previous(u)) != -1){
            Var e = outer->getEdgeVar(d.incomingEdge(u));
            lbool val = outer->value(e);
            assert(outer->value(e) == l_True);
            conflict.push(mkLit(e, true));
            u = p;

        }
    }
    outer->num_learnt_paths++;
    outer->learnt_path_clause_length += (conflict.size() - 1);

    stats_under_clause_length += (conflict.size() - 1);
    double elapsed = rtime(2) - starttime;
    stats_under_conflict_time += elapsed;

}

template<typename Weight, typename Graph>
void DistanceDetector<Weight, Graph>::buildUnweightedDistanceGTReason(int node, int within_steps, vec<Lit>& conflict){
    stats_unweighted_gt_reasons++;
    stats_over_conflicts++;
    double starttime = rtime(2);

    int u = node;
    bool reaches = overapprox_unweighted_distance_detector->connected(node);


    if(!reaches && within_steps >= g_over.nodes() && opt_conflict_min_cut && conflict_flow){

        g_over.drawFull(false);
        cut.clear();
        Weight f;

        assert(conflict_flow->getSource() == source);
        conflict_flow->setSink(node);
        f = conflict_flow->minCut(cut);
        assert(f == cut.size()); //because edges are only ever infinity or 1
        if(f != cut.size()){
            throw std::runtime_error("Bad learnt cut");
        }
        for(int i = 0; i < cut.size(); i++){
            MaxFlowEdge e = cut[i];
            int cut_id = e.id;
            assert(cut_id % 2 == 0);
            Lit l = mkLit(outer->getEdgeVar(cut_id / 2), false);
            assert(outer->value(l) == l_False);
            conflict.push(l);
        }
        outer->learnt_cut_clause_length += (conflict.size() - 1);

        stats_over_clause_length += (conflict.size() - 1);
        double elapsed = rtime(2) - starttime;
        stats_over_conflict_time += elapsed;
        return;
    }
    cutGraph.clearHistory();
    outer->stats_mc_calls++;
    {


        to_visit.clear();
        to_visit.push(node);
        seen.clear();
        seen.growTo(outer->nNodes());
        seen[node] = true;

        do{

            assert(to_visit.size());
            int u = to_visit.last();
            assert(u != source);
            to_visit.pop();
            assert(seen[u]);
            // then add all its incoming disabled edges to the cut, and visit any unseen, non-disabled incoming.edges()
            for(int i = 0; i < g_over.nIncoming(u); i++){
                int v = outer->getEdgeVar(g_over.incoming(u, i).id);
                int from = g_over.incoming(u, i).node;

                //Note: the variable has to not only be assigned false, but assigned false earlier in the trail than the reach variable...
                int edge_num = outer->getEdgeID(v);                // v-outer->min_edge_var;

                if(from == u){
                    assert(g_over.getEdge(edge_num).to == u);
                    assert(g_over.getEdge(edge_num).from == u);
                    continue;                //Self loops are allowed, but just make sure nothing got flipped around...
                }
                assert(from != u);

                if(has_unweighted_shortest_paths_overapprox && reaches){
                    //This is an optional optimization: if we know that even with all possible edges enabled, the shortest path to from + 1 is >= than the current distance to this node, enabling this edge cannot decrease the shortest path,
                    //and so we don't need to consider this edge
                    int current_dist = overapprox_unweighted_distance_detector->distance(u);
                    assert(current_dist > 0);
                    int over_approx_dist = unweighted_over_approx_shortest_paths[from] + 1;
                    if(over_approx_dist >= current_dist){
                        //following this edge cannot shorten this path, so skip it.
                        stats_gt_unweighted_edges_skipped++;
                        continue;
                    }
                }

                if(outer->value(v) == l_False){
                    //note: we know we haven't seen this edge variable before, because we know we haven't visited this node before
                    //if we are already planning on visiting the from node, then we don't need to include it in the conflict (is this correct?)
                    conflict.push(mkLit(v, false));
                }else if(from != source){
                    //for distance analysis, we _can_ end up reaching source.

                    //even if it is undef? probably...
                    if(!seen[from]){
                        seen[from] = true;
                        to_visit.push(from);
                    }
                }
            }
        }while(to_visit.size());

    }

    outer->num_learnt_cuts++;
    outer->learnt_cut_clause_length += (conflict.size() - 1);

    stats_over_clause_length += (conflict.size() - 1);
    double elapsed = rtime(2) - starttime;
    stats_over_conflict_time += elapsed;

}

template<typename Weight, typename Graph>
void DistanceDetector<Weight, Graph>::printSolution(std::ostream& write_to){


}


template<typename Weight, typename Graph>
void DistanceDetector<Weight, Graph>::buildReason(Lit p, vec<Lit>& reason, CRef marker){

    if(marker == unweighted_underprop_marker){
        reason.push(p);
        Var v = var(p);
        int u = getNode(v);
        buildUnweightedDistanceLEQReason(u, reason);
    }else if(marker == unweighted_overprop_marker){
        reason.push(p);

        //the reason is a cut separating p from s;
        //We want to find a min-cut in the full graph separating, where activated edges (ie, those still in antig) are weighted infinity, and all others are weighted 1.

        //This is a cut that describes a minimal set of edges which are disabled in the current graph, at least one of which would need to be activated in order for s to reach p
        //assign the mincut edge weights if they aren't already assigned.

        Var v = var(p);

        int t = getNode(v); // v- var(reach_lits[d][0]);
        int within = getMaximumDistance(v);
        buildUnweightedDistanceGTReason(t, within, reason);

    }else{
        throw std::runtime_error("Internal error in distance detector");
        assert(false);
    }
    outer->toSolver(reason);
}

template<typename Weight, typename Graph>
void DistanceDetector<Weight, Graph>::updateShortestPaths(){
    if(opt_shortest_path_prune_dist && outer->decisionLevel() == 0){
        //only update these distances at level 0, to ensure they are a valid over approximate of the shortest possible path to each node

        has_unweighted_shortest_paths_overapprox = true;
        unweighted_over_approx_shortest_paths.growTo(g_under.nodes());
        for(int i = 0; i < g_over.nodes(); i++){
            unweighted_over_approx_shortest_paths[i] = overapprox_unweighted_distance_detector->distance(i);
        }

    }
}

template<typename Weight, typename Graph>
void DistanceDetector<Weight, Graph>::preprocess(){
    is_changed_under.growTo(g_under.nodes());
    is_changed_over.growTo(g_over.nodes());

    if(!overapprox_unweighted_distance_detector){
        buildUnweightedSATConstraints(false);
    }else if(!underapprox_unweighted_distance_detector){
        //can optimize
        buildUnweightedSATConstraints(true);
    }


}


template<typename Weight, typename Graph>
bool DistanceDetector<Weight, Graph>::propagate(vec<Lit>& conflict){
    if(!underapprox_unweighted_distance_detector)
        return true;

    bool skipped_positive = false;
    if(!opt_detect_pure_theory_lits || unassigned_positives > 0){
        double startdreachtime = rtime(2);
        stats_under_updates++;
        underapprox_unweighted_distance_detector->update();
        double reachUpdateElapsed = rtime(2) - startdreachtime;
        stats_under_update_time += reachUpdateElapsed;
    }else{
        skipped_positive = true;
        stats_skipped_under_updates++;
    }

    bool skipped_negative = false;

    if(!opt_detect_pure_theory_lits || unassigned_negatives > 0){
        double startunreachtime = rtime(2);
        stats_over_updates++;
        overapprox_unweighted_distance_detector->update();
        double unreachUpdateElapsed = rtime(2) - startunreachtime;
        stats_over_update_time += unreachUpdateElapsed;
        updateShortestPaths(); //only needed for the shortest path theory
    }else{
        skipped_negative = true;
        stats_skipped_over_updates++;
    }

    if(opt_rnd_shuffle){
        randomShuffle(rnd_seed, changed);
    }

    while(changed.size()){
        int sz = changed.size();
        int u = changed.last().u;
        bool polarity = changed.last().polarity;
        assert(polarity ? is_changed_under[u] : is_changed_over[u]);
        for(int i = 0; i < unweighted_dist_lits[u].size(); i++){
            int& min_distance = unweighted_dist_lits[u][i].min_unweighted_distance;

            Var v = var(unweighted_dist_lits[u][i].l);
            Lit l;

            if(underapprox_unweighted_distance_detector && polarity
               && underapprox_unweighted_distance_detector->connected(u)
               && underapprox_unweighted_distance_detector->distance(u) <= min_distance){
                l = mkLit(v, false);
            }else if(overapprox_unweighted_distance_detector && !polarity
                     && (!overapprox_unweighted_distance_detector->connected(u)
                         || overapprox_unweighted_distance_detector->distance(u) > min_distance)){
                l = mkLit(v, true);
            }else{
                continue;
            }

            bool reach = !sign(l);
            if(outer->value(l) == l_True){
                //do nothing
            }else if(outer->value(l) == l_Undef){
#ifdef DEBUG_GRAPH
                assert(outer->dbg_propgation(l));
#endif
#ifdef DEBUG_SOLVER
                if(S->dbg_solver)
                S->dbg_check_propagation(l);
#endif

                if(reach)
                    outer->enqueue(l, unweighted_underprop_marker);
                else
                    outer->enqueue(l, unweighted_overprop_marker);

            }else if(outer->value(l) == l_False){
                conflict.push(l);

                if(reach){

                    //conflict
                    //The reason is a path in g from to s in d
                    buildUnweightedDistanceLEQReason(u, conflict);
                    //add it to s
                    //return it as a conflict

                }else{
                    //The reason is a cut separating s from t
                    buildUnweightedDistanceGTReason(u, getMaximumDistance(v), conflict);

                }
#ifdef DEBUG_GRAPH
                for(int i = 0;i<conflict.size();i++)
                assert(outer->value(conflict[i])==l_False);
#endif
#ifdef DEBUG_SOLVER
                if(S->dbg_solver)
                S->dbg_check(conflict);
#endif
                outer->toSolver(conflict);
                return false;
            }
        }
        if(sz == changed.size()){
            //it is possible, in rare cases, for literals to have been added to the chagned list while processing the changed list.
            //in that case, don't pop anything off the changed list, and instead accept that this literal may be re-processed
            //(should only be possible to have this happen at most twice per propagation call, so shouldn't matter too much)
            assert(sz == changed.size());
            assert(changed.last().u == u);
            if(polarity){
                is_changed_under[u] = false;
            }else{
                is_changed_over[u] = false;
            }
            changed.pop();
        }
    }


    return true;
}

//Return the path (in terms of nodes)
template<typename Weight, typename Graph>
bool DistanceDetector<Weight, Graph>::getModel_Path(int node, std::vector<int>& store_path){
    store_path.clear();
    Reach& d = *underapprox_path_detector;
    d.update();
    if(!d.connected(node))
        return false;
    int u = node;
    int p;
    while((p = d.previous(u)) != -1){
        Var e = outer->getEdgeVar(d.incomingEdge(u));
        lbool val = outer->value(e);
        assert(outer->value(e) == l_True);
        store_path.push_back(u);
        u = p;
    }
    assert(u == this->source);
    store_path.push_back(u);
    std::reverse(store_path.begin(), store_path.end());
    return true;
}

template<typename Weight, typename Graph>
bool DistanceDetector<Weight, Graph>::getModel_PathByEdgeLit(int node, std::vector<Lit>& store_path){
    store_path.clear();
    Reach& d = *underapprox_path_detector;
    d.update();
    if(!d.connected(node))
        return false;
    int u = node;
    int p;
    while((p = d.previous(u)) != -1){
        Var e = outer->getEdgeVar(d.incomingEdge(u));
        assert(outer->value(e) == l_True);
        lbool val = outer->value(e);
        assert(outer->value(e) == l_True);
        store_path.push_back(mkLit(outer->toSolver(e), false));
        u = p;
    }
    assert(u == this->source);
    std::reverse(store_path.begin(), store_path.end());
    return true;
}

template<typename Weight, typename Graph>
bool DistanceDetector<Weight, Graph>::checkSatisfied(){
    UnweightedDijkstra<Weight, Graph> under(source, g_under);
    UnweightedDijkstra<Weight, Graph> over(source, g_over);
    under.update();
    over.update();
    for(int j = 0; j < unweighted_dist_lits.size(); j++){
        for(int k = 0; k < unweighted_dist_lits[j].size(); k++){
            Lit l = unweighted_dist_lits[j][k].l;
            int dist = unweighted_dist_lits[j][k].min_unweighted_distance;
            lbool val = outer->value(l);
            if(l != lit_Undef){
                int node = getNode(var(l));

                if(outer->value(l) == l_True){
                    if(!under.connected(node) || under.distance(node) > dist){
                        return false;
                    }
                }else if(outer->value(l) == l_False){
                    if(under.connected(node) && over.distance(node) <= dist){
                        return false;
                    }
                }else{
                    if(under.connected(node) && under.distance(node) <= dist){
                        return false;
                    }
                    if(!over.connected(node) || over.distance(node) > dist){
                        return false;
                    }
                }
            }
        }
    }

    return true;
}

template<typename Weight, typename Graph>
Lit DistanceDetector<Weight, Graph>::decide(CRef& decision_reason){
    if(!opt_decide_graph_distance || !overapprox_unweighted_distance_detector)
        return lit_Undef;
    DistanceDetector* r = this;
    auto* over = r->overapprox_unweighted_distance_detector;

    auto* under = r->underapprox_unweighted_distance_detector;

    //we can probably also do something similar, but with cuts, for nodes that are decided to be unreachable.

    //ok, for each node that is assigned reachable, but that is not actually reachable in the under approx, decide an edge on a feasible path

    //this can be obviously more efficient


    for(int k = 0; k < unweighted_dist_lits.size(); k++){
        for(int n = 0; n < unweighted_dist_lits[k].size(); n++){
            Lit l = unweighted_dist_lits[k][n].l;
            int min_dist = unweighted_dist_lits[k][n].min_unweighted_distance;
            if(l == lit_Undef)
                continue;
            int j = r->getNode(var(l));
            if(outer->value(l) == l_True){
                if(opt_decide_graph_pos){

                    assert(over->distance(j) <=
                           min_dist);//else we would already be in conflict before this decision was attempted!
                    if(under->distance(j) > min_dist){
                        //then lets try to connect this
                        //static vec<bool> print_path;

                        assert(over->connected(j));                    //Else, we would already be in conflict

                        int p = j;
                        int last = j;
                        int last_edge = -1;
                        if(!opt_use_random_path_for_decisions){
                            //ok, read back the path from the over to find a candidate edge we can decide
                            //find the earliest unconnected node on this path
                            over->update();
                            p = j;
                            last = j;
                            int dist = 0;
                            while(under->distance(p) >= min_dist - dist){

                                last = p;
                                assert(p != r->source);
                                int prev = over->previous(p);
                                last_edge = over->incomingEdge(p);
                                assert(over->distance(p) <= min_dist - dist);
                                dist += 1;                        //should really be weighted
                                p = prev;

                            }
                        }else{
                            //This won't work (without modification) because we need to constrain these paths to ones of maximum real distance < min_dist.
                            //Randomly re-weight the graph sometimes
                            if(drand(rnd_seed) < opt_decide_graph_re_rnd){

                                for(int i = 0; i < outer->nEdges(); i++){
                                    double w = drand(rnd_seed);
                                    rnd_weight[i] = w;
                                }
                            }
                            rnd_path->update();
                            //derive a random path in the graph
                            p = j;
                            last = j;
                            assert(rnd_path->connected(p));
                            while(!under->connected(p)){

                                last = p;

                                assert(p != source);
                                int prev = rnd_path->previous(p);
                                last_edge = rnd_path->incomingEdge(p);
                                p = prev;
                                assert(p >= 0);
                            }

                        }

                        //ok, now pick some edge p->last that will connect p to last;

                        assert(p > -1);
                        if(p > -1){
                            assert(g_over.getEdge(last_edge).from == p);
                            assert(g_over.getEdge(last_edge).to == last);
                            Var v = outer->getEdgeVar(last_edge); ///outer->edges[p][last].v;
                            if(outer->value(v) == l_Undef){
                                return mkLit(v, false);
                            }else{
                                assert(outer->value(v) != l_True);
                            }
                        }

                    }
                }
            }else if(outer->value(l) == l_False){
                if(opt_decide_graph_neg){

                    //assert(over->distance(j)<=min_dist);//else we would already be in conflict before this decision was attempted!

                    if(over->distance(j) <= min_dist && under->distance(j) > min_dist){
                        //then lets try to disconnect this node from source by walking back along the path in the over approx, and disabling the first unassigned edge we see.
                        //(there must be at least one such edge, else the variable would be connected in the under approximation as well - in which case it would already have been propagated.
                        int p = j;
                        int last = j;
                        int dist = 0;
                        while(under->distance(p) > min_dist - dist){
                            last = p;
                            assert(p != source);
                            int prev = over->previous(p);
                            int edgeid = over->incomingEdge(p);
                            assert(g_over.getEdge(edgeid).from == prev);
                            assert(g_over.getEdge(edgeid).to == p);
                            Var v = outer->getEdgeVar(edgeid); //outer->edges[prev][p].v;
                            if(outer->value(v) == l_Undef){
                                return mkLit(v, true);
                            }else{
                                assert(outer->value(v) != l_False);
                            }
                            assert(over->distance(p) <= min_dist - dist);
                            dist += 1;                                //should really be weighted
                            p = prev;

                        }
                    }
                }
            }
        }
    }
    return lit_Undef;
};

#include <gmpxx.h>

template
class Monosat::DistanceDetector<int>;

template
class Monosat::DistanceDetector<int64_t>;

template
class Monosat::DistanceDetector<double>;

template
class Monosat::DistanceDetector<mpq_class>;

template
class Monosat::DistanceDetector<int, DynamicBackGraph<int>>;

template
class Monosat::DistanceDetector<int64_t, DynamicBackGraph<int64_t>>;

template
class Monosat::DistanceDetector<double, DynamicBackGraph<double>>;

template
class Monosat::DistanceDetector<mpq_class, DynamicBackGraph<mpq_class>>;
