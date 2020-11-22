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
#include "monosat/dgl/EdmondsKarp.h"
#include "monosat/dgl/EdmondsKarpAdj.h"
#include "monosat/dgl/KohliTorr.h"
#include "monosat/dgl/EdmondsKarpDynamic.h"
#include "monosat/dgl/Dinics.h"
#include "monosat/dgl/DynamicGraph.h"
#include "monosat/dgl/Graph.h"
#include "monosat/dgl/DinicsLinkCut.h"
#include "monosat/graph/WeightedDistanceDetector.h"
#include "monosat/graph/GraphTheory.h"
#include "monosat/mtl/Rnd.h"
#include "monosat/mtl/Vec.h"
#include "monosat/utils/Options.h"
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <iomanip>

using namespace Monosat;

template<typename Weight, typename Graph>
WeightedDistanceDetector<Weight, Graph>::WeightedDistanceDetector(int _detectorID, GraphTheorySolver<Weight>* outer,
                                                                  Graph& _g, Graph& _antig, Graph& cutGraph, int from,
                                                                  double seed) :
        Detector(_detectorID), outer(outer), g_under(_g), g_over(_antig), cutGraph(cutGraph), source(from),
        rnd_seed(seed){
    max_unweighted_distance = 0;
    rnd_path = NULL;
    bvTheory = outer->bvTheory;
    constraintsBuilt = -1;
    first_reach_var = var_Undef;
    stats_pure_skipped = 0;


    if(opt_use_random_path_for_decisions){
        rnd_weight.clear();
        rnd_path = new WeightedDijkstra<Weight, Graph, double>(from, _antig, rnd_weight);
        for(int i = 0; i < outer->nEdges(); i++){
            double w = drand(rnd_seed);

            rnd_weight.push_back(w);
        }

    }


    //select the _weighted_ distance detectors
    positiveDistanceStatus = new WeightedDistanceDetector<Weight, Graph>::DistanceStatus(*this, true);
    negativeDistanceStatus = new WeightedDistanceDetector<Weight, Graph>::DistanceStatus(*this, false);

    if(outer->hasBitVectorEdges()){
        printf("Note: falling back on Dijkstra for shortest paths, because edge weights are bitvectors\n");
        //ramel reps doesn't support bvs yet
        underapprox_weighted_distance_detector =
                new Dijkstra<Weight, Graph, typename WeightedDistanceDetector<Weight, Graph>::DistanceStatus>(from, _g,
                                                                                                              *positiveDistanceStatus,
                                                                                                              0);
        overapprox_weighted_distance_detector = new Dijkstra<Weight, Graph, typename WeightedDistanceDetector<Weight, Graph>::DistanceStatus>(
                from, _antig, *negativeDistanceStatus, 0);
        underapprox_weighted_path_detector = underapprox_weighted_distance_detector;
    }else if(distalg == DistAlg::ALG_RAMAL_REPS){

        underapprox_weighted_distance_detector =
                new RamalReps<Weight, Graph, typename WeightedDistanceDetector<Weight, Graph>::DistanceStatus>(from, _g,
                                                                                                               *(positiveDistanceStatus),
                                                                                                               -2);
        overapprox_weighted_distance_detector =
                new RamalReps<Weight, Graph, typename WeightedDistanceDetector<Weight, Graph>::DistanceStatus>(from,
                                                                                                               _antig,
                                                                                                               *(negativeDistanceStatus),
                                                                                                               -2);
        underapprox_weighted_path_detector = underapprox_weighted_distance_detector; //new Dijkstra<Weight>(from, _g);
    }else{
        underapprox_weighted_distance_detector =
                new Dijkstra<Weight, Graph, typename WeightedDistanceDetector<Weight, Graph>::DistanceStatus>(from, _g,
                                                                                                              *positiveDistanceStatus,
                                                                                                              0);
        overapprox_weighted_distance_detector = new Dijkstra<Weight, Graph, typename WeightedDistanceDetector<Weight, Graph>::DistanceStatus>(
                from, _antig, *negativeDistanceStatus, 0);
        underapprox_weighted_path_detector = underapprox_weighted_distance_detector;
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


    weighted_underprop_marker = outer->newReasonMarker(getID());
    weighted_overprop_marker = outer->newReasonMarker(getID());

    weighted_underprop_bv_marker = outer->newReasonMarker(getID());
    weighted_overprop_bv_marker = outer->newReasonMarker(getID());
}


template<typename Weight, typename Graph>
void WeightedDistanceDetector<Weight, Graph>::addWeightedShortestPathLit(int from, int to, Var outer_reach_var,
                                                                         Weight within_distance, bool strictComparison){
    g_under.invalidate();
    g_over.invalidate();
    Var reach_var = outer->newVar(outer_reach_var, getID());
    assert(from == source);
    weighted_dist_lits.push(WeightedDistLit{mkLit(reach_var), to, within_distance, strictComparison});

    if(first_reach_var == var_Undef){
        first_reach_var = reach_var;
    }else{
        assert(reach_var > first_reach_var);
    }
    while(reach_lit_map.size() <= reach_var - first_reach_var){
        reach_lit_map.push({-1, -1, None});
    }
    reach_lit_map[reach_var - first_reach_var] = {to, weighted_dist_lits.size() - 1, WeightedConstLit};
}

template<typename Weight, typename Graph>
void WeightedDistanceDetector<Weight, Graph>::addWeightedShortestPathBVLit(int from, int to, Var outer_reach_var,
                                                                           const BitVector<Weight>& bv,
                                                                           bool strictComparison){
    g_under.invalidate();
    g_over.invalidate();

    Var reach_var = outer->newVar(outer_reach_var, getID());
    assert(from == source);

    weighted_dist_bv_lits.push(WeightedDistBVLit{mkLit(reach_var), to, bv, strictComparison, nullptr});
    //sort(weighted_dist_lits);
    if(opt_graph_bv_prop){
        DistanceOp* distOp = new DistanceOp(*bvTheory, this, bv.getID(), to, strictComparison, mkLit(reach_var));
        weighted_dist_bv_lits.last().op = distOp;
    }

    if(first_reach_var == var_Undef){
        first_reach_var = reach_var;
    }else{
        assert(reach_var > first_reach_var);
    }
    while(reach_lit_map.size() <= reach_var - first_reach_var){
        reach_lit_map.push({-1, -1, None});
    }
    reach_lit_map[reach_var - first_reach_var] = {to, weighted_dist_bv_lits.size() - 1, WeightedBVLit};
}


template<typename Weight, typename Graph>
void WeightedDistanceDetector<Weight, Graph>::DistanceStatus::setReachable(int u, bool reachable){

}

template<typename Weight, typename Graph>
void
WeightedDistanceDetector<Weight, Graph>::DistanceStatus::setMininumDistance(int u, bool reachable, Weight& distance){

}


template<typename Weight, typename Graph>
void WeightedDistanceDetector<Weight, Graph>::printSolution(std::ostream& write_to){

}

template<typename Weight, typename Graph>
void WeightedDistanceDetector<Weight, Graph>::DistanceOp::analyzeReason(bool compareOver, Comparison op, Weight to,
                                                                        vec<Lit>& conflict){
//watch out - might need to backtrack the graph theory appropriately, here...
    GraphTheorySolver<Weight>::GraphTheoryOp::analyzeReason(compareOver, op, to, conflict);
    if(!compareOver){
        Lit l = this->comparisonLit;
        assert(outer->outer->value(l) == l_True);

        conflict.push(outer->outer->toSolver(~l));
        assert(op == Comparison::gt || op == Comparison::geq);
        assert(this->strictCompare == (op == Comparison::gt));
        outer->analyzeDistanceGTReason(this->to, to, conflict, strictCompare);
    }else{
        Lit l = this->comparisonLit;
        assert(outer->outer->value(l) == l_False);
        conflict.push(outer->outer->toSolver(l));
        assert(op == Comparison::lt || op == Comparison::leq);
        assert(this->strictCompare == (op == Comparison::lt));
        outer->analyzeDistanceLEQReason(this->to, to, conflict, strictCompare);

    }
    GraphTheorySolver<Weight>::GraphTheoryOp::completeAnalysis();
}

template<typename Weight, typename Graph>
void WeightedDistanceDetector<Weight, Graph>::analyzeDistanceLEQReason(int to, Weight& min_distance, vec<Lit>& conflict,
                                                                       bool strictComparison){
    Distance<Weight>& d = *underapprox_weighted_path_detector;
    underapprox_weighted_distance_detector->update();
    Weight& actual_dist = underapprox_weighted_distance_detector->distance(to);

    d.update();
    stats_under_conflicts++;
    assert(outer->dbg_reachable(d.getSource(), to));
    assert(underapprox_weighted_distance_detector->connected(to));
    assert(
            underapprox_weighted_distance_detector->distance(to) <= min_distance
            && underapprox_weighted_distance_detector->distance(to)
               != underapprox_weighted_distance_detector->unreachable());
    assert(d.connected_unchecked(to));
    if(!d.connected(to) || d.distance(to) > min_distance){
        throw std::runtime_error("Internal error in shortest path detector");
    }

    //the reason that the distance is less than or equal to min_distance is because the shortest path is less than this weight
    {
        int u = to;
        int p;
        while((p = d.previous(u)) != -1){

            int edgeID = d.incomingEdge(u);
            Var e = outer->getEdgeVar(edgeID);
            lbool val = outer->value(e);
            assert(outer->value(e) == l_True);
            if(!g_under.isConstant(edgeID))
                conflict.push(outer->toSolver(mkLit(e, true)));
            if(outer->hasBitVector(edgeID)){


                outer->bvTheory->addAnalysis(Comparison::leq, outer->getEdgeBV(edgeID).getID(),
                                             g_under.getWeight(edgeID));
                //outer->buildBVReason(outer->getEdgeBV(edgeID).getID(),Comparison::leq,g_under.getWeight(edgeID),conflict);
            }

            u = p;
        }
    }
}

template<typename Weight, typename Graph>
void WeightedDistanceDetector<Weight, Graph>::buildDistanceLEQReason(int to, Weight& min_distance, vec<Lit>& conflict,
                                                                     bool strictComparison){
    stats_distance_leq_reasons++;
    double starttime = rtime(2);

    analyzeDistanceLEQReason(to, min_distance, conflict, strictComparison);
    if(outer->bvTheory){
        outer->bvTheory->analyze(conflict);
    }
    outer->num_learnt_paths++;
    outer->learnt_path_clause_length += (conflict.size() - 1);

    stats_under_clause_length += (conflict.size() - 1);
    double elapsed = rtime(2) - starttime;
    stats_under_conflict_time += elapsed;

}

template<typename Weight, typename Graph>
void WeightedDistanceDetector<Weight, Graph>::analyzeDistanceGTReason(int to, Weight& min_distance, vec<Lit>& conflict,
                                                                      bool strictComparison){
    bool reaches = overapprox_weighted_distance_detector->connected(to);
    if(!reaches && opt_conflict_min_cut && conflict_flow){

        cut.clear();
        Weight f;

        assert(conflict_flow->getSource() == source);
        conflict_flow->setSink(to);
        f = conflict_flow->minCut(cut);

        assert(f == cut.size()); //because edges are only ever infinity or 1

        for(int i = 0; i < cut.size(); i++){
            MaxFlowEdge e = cut[i];
            int cut_id = e.id;
            assert(cut_id % 2 == 0);
            int edgeID = cut_id / 2;
            if(!g_over.isConstant(edgeID)){
                Lit l = mkLit(outer->getEdgeVar(edgeID), false);
                assert(outer->value(l) == l_False);
                conflict.push(outer->toSolver(l));
            }
        }

        return;
    }
    Weight& actual_dist = overapprox_weighted_distance_detector->distance(to);

    bool connected = overapprox_weighted_distance_detector->connected(to);
#ifdef DEBUG_GRAPH
    Dijkstra<Weight> d(source, g_over);
    Weight & expected = d.distance(to);
    assert(expected == actual_dist);
    assert(!overapprox_weighted_distance_detector->connected(to) || (strictComparison? actual_dist > min_distance: actual_dist >= min_distance));
#endif

    {


        to_visit.clear();
        to_visit.push(to);
        seen.clear();
        seen.growTo(outer->nNodes());
        seen[to] = true;

        do{

            assert(to_visit.size());
            int u = to_visit.last();
            assert(u != source);
            to_visit.pop();
            assert(seen[u]);
            //add all of this node's incoming disabled edges to the cut, and visit any unseen, non-disabled incoming.edges()

            for(int i = 0; i < g_over.nIncoming(u); i++){
                int from = g_over.incoming(u, i).node;
                int edge_num = g_over.incoming(u, i).id;

                //Note: the variable has to not only be assigned false, but assigned false earlier in the trail than the reach variable...
                Var edge_enabled = outer->getEdgeVar(edge_num);
                if(from == u){

                    continue; //Self loops are allowed, but just make sure nothing got flipped around...
                }
                assert(from != u);

                if(has_weighted_shortest_paths_overapprox && reaches){
                    //This is an optional optimization: if we know that even with all possible edges enabled, the shortest path to from + 1 is >= than the current distance to this node, enabling this edge cannot decrease the shortest path,
                    //and so we don't need to consider this edge
                    Weight current_dist = overapprox_weighted_distance_detector->distance(u);
                    Weight over_approx_dist = over_approx_shortest_paths[from] + g_over.getWeight(edge_num);
                    if(over_approx_dist >= current_dist){
                        stats_gt_weighted_edges_skipped++;
                        //following this edge cannot shorten this path, so skip it.
                        continue;
                    }
                }

                if(outer->value(edge_enabled) == l_False){
                    //note: we know we haven't seen this edge variable before, because we know we haven't visited this node before
                    //if we are already planning on visiting the from node, then we don't need to include it in the conflict (is this correct?)
                    if(!g_over.isConstant(edge_num)){
                        conflict.push(outer->toSolver(mkLit(edge_enabled)));
                    }
                }else{

                    if(reaches){
                        //if the edge _is_ enabled, and the node _is_ reachable, and the edge weight is symbolic
                        //then part of the reason the shortest path is too int64_t is that this edge is not less than its current weight.
                        if(!outer->constantWeight(edge_num)){
                            Weight w = g_over.getWeight(edge_num);


                            if(strictComparison){
                                outer->bvTheory->addAnalysis(Comparison::geq, outer->getEdgeBV(edge_num).getID(), w);
                            }else{
                                outer->bvTheory->addAnalysis(Comparison::geq, outer->getEdgeBV(edge_num).getID(), w);
                            }

                        }
                    }
                    //for distance analysis, we _can_ end up reaching source.
                    if(from != source){
                        //even if it is undef? probably...
                        if(!seen[from]){
                            seen[from] = true;
                            to_visit.push(from);
                        }
                    }
                }
            }
        }while(to_visit.size());
    }
}

template<typename Weight, typename Graph>
void WeightedDistanceDetector<Weight, Graph>::buildDistanceGTReason(int to, Weight& min_distance, vec<Lit>& conflict,
                                                                    bool strictComparison){
    stats_distance_gt_reasons++;
    stats_over_conflicts++;
    double starttime = rtime(2);
    analyzeDistanceGTReason(to, min_distance, conflict, strictComparison);
    if(outer->bvTheory){
        outer->bvTheory->analyze(conflict);
    }
    outer->learnt_cut_clause_length += (conflict.size() - 1);

    stats_over_clause_length += (conflict.size() - 1);
    double elapsed = rtime(2) - starttime;
    outer->mctime += elapsed;

}

template<typename Weight, typename Graph>
void WeightedDistanceDetector<Weight, Graph>::buildReason(Lit p, vec<Lit>& reason, CRef marker){

    if(marker == weighted_underprop_marker){
        reason.push(outer->toSolver(p));
        assert(getLitType(p) == WeightedConstLit);
        Var v = var(p);
        WeightedDistLit& dist_lit = getDistLit(v);
        bool strictComparison = dist_lit.strictComparison;
        Lit l = dist_lit.l;
        lbool val = outer->value(l);

        int to = dist_lit.u;
        Weight& min_dist = dist_lit.min_distance;
        Weight& over_dist = underapprox_weighted_distance_detector->distance(to);
        Weight& under_dist = overapprox_weighted_distance_detector->distance(to);

        buildDistanceLEQReason(to, min_dist, reason, strictComparison);

    }else if(marker == weighted_overprop_marker){
        reason.push(outer->toSolver(p));
        assert(sign(p));
        assert(getLitType(p) == WeightedConstLit);
        Var v = var(p);
        WeightedDistLit& dist_lit = getDistLit(v);
        bool strictComparison = dist_lit.strictComparison;
        Lit l = dist_lit.l;
        lbool val = outer->value(l);

        int to = dist_lit.u;
        Weight& min_dist = dist_lit.min_distance;
        Weight& over_dist = underapprox_weighted_distance_detector->distance(to);
        Weight& under_dist = overapprox_weighted_distance_detector->distance(to);

        buildDistanceGTReason(to, min_dist, reason, strictComparison);

    }else if(marker == weighted_underprop_bv_marker){
        reason.push(outer->toSolver(p));
        assert(getLitType(p) == WeightedBVLit);

        Var v = var(p);
        WeightedDistBVLit& dist_lit = getBVDistLit(v);
        Lit l = dist_lit.l;
        int to = dist_lit.u;
        bool strictComparison = dist_lit.strictComparison;
        BitVector<Weight>& bv = dist_lit.bv;
        Weight& min_dist_under = bv.getUnder();
        Weight& min_dist_over = bv.getOver();
        Weight& under_dist = underapprox_weighted_distance_detector->distance(to);
        Weight& over_dist = overapprox_weighted_distance_detector->distance(to);

        if(strictComparison){

            outer->bvTheory->addAnalysis(Comparison::gt, bv.getID(), under_dist);

            buildDistanceLEQReason(to, min_dist_under, reason, true);
        }else{
            outer->bvTheory->addAnalysis(Comparison::geq, bv.getID(), under_dist);
            buildDistanceLEQReason(to, min_dist_under, reason, false);
        }


    }else if(marker == weighted_overprop_bv_marker){
        reason.push(outer->toSolver(p));
        assert(sign(p));
        assert(getLitType(p) == WeightedBVLit);

        Var v = var(p);
        WeightedDistBVLit& dist_lit = getBVDistLit(v);
        Lit l = dist_lit.l;
        int to = dist_lit.u;
        bool strictComparison = dist_lit.strictComparison;
        BitVector<Weight>& bv = dist_lit.bv;
        Weight& min_dist_under = bv.getUnder();
        Weight& min_dist_over = bv.getOver();
        Weight& under_dist = underapprox_weighted_distance_detector->distance(to);
        Weight& over_dist = overapprox_weighted_distance_detector->distance(to);


        if(strictComparison){
            if(overapprox_weighted_distance_detector->connected(to)){
                outer->bvTheory->addAnalysis(Comparison::leq, bv.getID(), over_dist);
            }
            buildDistanceGTReason(to, min_dist_over, reason, false);
        }else{
            if(overapprox_weighted_distance_detector->connected(to)){
                outer->bvTheory->addAnalysis(Comparison::lt, bv.getID(), over_dist);
            }
            buildDistanceGTReason(to, min_dist_over, reason, true);
        }


    }else{
        throw std::runtime_error("Internal error in distance detector");
        assert(false);
    }

}

template<typename Weight, typename Graph>
void WeightedDistanceDetector<Weight, Graph>::updateShortestPaths(){
    if(opt_shortest_path_prune_dist && outer->decisionLevel() == 0){
        //only update these distances at level 0, to ensure they are a valid over approximate of the shortest possible path to each node

        has_weighted_shortest_paths_overapprox = true;
        over_approx_shortest_paths.growTo(g_under.nodes());
        for(int i = 0; i < g_over.nodes(); i++){
            over_approx_shortest_paths[i] = overapprox_weighted_distance_detector->distance(i);
        }

    }
}

template<typename Weight, typename Graph>
void WeightedDistanceDetector<Weight, Graph>::preprocess(){
    is_changed_under.growTo(g_under.nodes());
    is_changed_over.growTo(g_over.nodes());


}


template<typename Weight, typename Graph>
bool WeightedDistanceDetector<Weight, Graph>::propagate(vec<Lit>& conflict){
    bool skipped_positive = false;
    if(!opt_detect_pure_theory_lits || unassigned_positives > 0){
        double startdreachtime = rtime(2);
        stats_under_updates++;
        underapprox_weighted_distance_detector->update();
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
        overapprox_weighted_distance_detector->update();
        double unreachUpdateElapsed = rtime(2) - startunreachtime;
        stats_over_update_time += unreachUpdateElapsed;
        updateShortestPaths(); //only needed for the shortest path theory
    }else{
        skipped_negative = true;
        stats_skipped_over_updates++;
    }
    if(opt_rnd_shuffle && weighted_dist_lits.size()){
        randomShuffle(rnd_seed, weighted_dist_lits);
    }
    if(opt_rnd_shuffle && weighted_dist_bv_lits.size()){
        randomShuffle(rnd_seed, weighted_dist_bv_lits);
    }

    //now, check for weighted distance lits
    for(auto& dist_lit : weighted_dist_lits){
        bool strictComparison = dist_lit.strictComparison;
        Lit l = dist_lit.l;
        lbool val = outer->value(l);

        int to = dist_lit.u;
        Weight& min_dist = dist_lit.min_distance;
        Weight& over_dist = underapprox_weighted_distance_detector->distance(to);
        Weight& under_dist = overapprox_weighted_distance_detector->distance(to);
        if(underapprox_weighted_distance_detector->connected(to)
           && (!strictComparison ? (underapprox_weighted_distance_detector->distance(to) <= min_dist) : (
                underapprox_weighted_distance_detector->distance(to) < min_dist))){
            if(outer->value(l) == l_True){
                //do nothing
            }else if(outer->value(l) == l_Undef){
                outer->enqueue(l, weighted_underprop_marker);
            }else if(outer->value(l) == l_False){
                //conflict
                conflict.push(outer->toSolver(l));

                buildDistanceLEQReason(to, min_dist, conflict, strictComparison);

                return false;
            }
        }
        if(!overapprox_weighted_distance_detector->connected(to)
           || (!strictComparison ? (overapprox_weighted_distance_detector->distance(to) > min_dist) : (
                overapprox_weighted_distance_detector->distance(to) >= min_dist))){
            if(outer->value(~l) == l_True){
                //do nothing
            }else if(outer->value(~l) == l_Undef){
                outer->enqueue(~l, weighted_overprop_marker);
            }else if(outer->value(~l) == l_False){
                //conflict
                conflict.push(outer->toSolver(~l));

                buildDistanceGTReason(to, min_dist, conflict, strictComparison);

                return false;
            }
        }
    }
    for(auto& dist_lit : weighted_dist_bv_lits){
        Lit l = dist_lit.l;
        int to = dist_lit.u;
        bool strictComparison = dist_lit.strictComparison;
        BitVector<Weight>& bv = dist_lit.bv;
        Weight& min_dist_under = bv.getUnder();
        Weight& min_dist_over = bv.getOver();
        Weight& under_dist = underapprox_weighted_distance_detector->distance(to);
        Weight& over_dist = overapprox_weighted_distance_detector->distance(to);


        if(underapprox_weighted_distance_detector->connected(to)
           && (strictComparison ? (under_dist < min_dist_under) : (under_dist <= min_dist_under))){
            if(outer->value(l) == l_True){
                //do nothing
            }else if(outer->value(l) == l_Undef){
                outer->enqueue(l, weighted_underprop_bv_marker);
            }else if(outer->value(l) == l_False){
                //conflict

                conflict.push(outer->toSolver(l));
                if(strictComparison){
                    outer->bvTheory->addAnalysis(Comparison::gt, bv.getID(), under_dist);
                    buildDistanceLEQReason(to, min_dist_under, conflict, true);
                }else{

                    outer->bvTheory->addAnalysis(Comparison::geq, bv.getID(), under_dist);
                    buildDistanceLEQReason(to, min_dist_under, conflict, false);
                }

                return false;
            }
        }
        if(!overapprox_weighted_distance_detector->connected(to)
           || (strictComparison ? (over_dist >= min_dist_over) : (over_dist >
                                                                  min_dist_over))){
            if(outer->value(~l) == l_True){
                //do nothing
            }else if(outer->value(~l) == l_Undef){
                outer->enqueue(~l, weighted_overprop_bv_marker);
            }else if(outer->value(~l) == l_False){
                //conflict
                conflict.push(outer->toSolver(~l));
                if(strictComparison){
                    if(overapprox_weighted_distance_detector->connected(to)){
                        outer->bvTheory->addAnalysis(Comparison::leq, bv.getID(), over_dist);

                    }
                    buildDistanceGTReason(to, min_dist_over, conflict, false);
                }else{
                    if(overapprox_weighted_distance_detector->connected(to)){
                        outer->bvTheory->addAnalysis(Comparison::lt, bv.getID(), over_dist);
                    }
                    buildDistanceGTReason(to, min_dist_over, conflict, true);
                }

                return false;
            }
        }


        assert(min_dist_under <= min_dist_over);
        assert(under_dist == -1 || over_dist <= under_dist);
        if(opt_graph_bv_prop){

            DistanceOp* distOp = dist_lit.op;
            if(distOp && !bv.isConst()){
                if(outer->value(l) == l_True){

                    //propagate over/under constraints on the output lit back to the bitvector theory
                    if(!overapprox_weighted_distance_detector->connected(to)){
                        //this should have triggered a conflict above
                        throw std::runtime_error("Internal error in distance detector");
                    }else{
                        if(!outer->assignBV(bv.getID(), strictComparison ? Comparison::gt : Comparison::geq, over_dist,
                                            distOp)){
                            throw std::runtime_error("Bad bv assignment in distance detector");
                        }

                    }
                }else if(outer->value(l) == l_False){
                    if(!underapprox_weighted_distance_detector->connected(to)){
                        //this is allowed, but what value should the bv take in this case?
                    }else{
                        if(!strictComparison && under_dist == 0){
                            //then we can't state that the bv is less than zero (since the distance on the edge must be positive), so
                            //actually this should have been caught as a conflict above
                            throw std::runtime_error("Bad bv assignment in distance detector");
                        }

                        if(!outer->assignBV(bv.getID(), strictComparison ? Comparison::leq : Comparison::lt, under_dist,
                                            distOp)){
                            throw std::runtime_error("Bad bv assignment in distance detector");
                        }
                    }
                }

            }
        }
    }
    return true;
}

//Return the path (in terms of nodes)
template<typename Weight, typename Graph>
bool WeightedDistanceDetector<Weight, Graph>::getModel_Path(int node, std::vector<int>& store_path){
    store_path.clear();
    Reach& d = *underapprox_weighted_path_detector;
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
bool WeightedDistanceDetector<Weight, Graph>::getModel_PathByEdgeLit(int node, std::vector<Lit>& store_path){
    store_path.clear();
    Reach& d = *underapprox_weighted_path_detector;
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
bool WeightedDistanceDetector<Weight, Graph>::checkSatisfied(){

    Dijkstra<Weight> under(source, g_under);
    Dijkstra<Weight> over(source, g_over);
    under.update();
    over.update();
    //now, check for weighted distance lits
    for(auto& dist_lit : weighted_dist_lits){
        Lit l = dist_lit.l;
        lbool val = outer->value(l);
        int to = dist_lit.u;
        Weight& min_dist = dist_lit.min_distance;
        Weight& over_dist = under.distance(to);
        Weight& under_dist = over.distance(to);
        bool strictComparison = dist_lit.strictComparison;
        if(!strictComparison){
            if(under.connected(to) && under.distance(to) <= min_dist){
                if(outer->value(l) == l_True){
                    //do nothing
                }else if(outer->value(l) == l_Undef){
                    outer->enqueue(l, weighted_underprop_marker);
                }else if(outer->value(l) == l_False){
                    return false;
                }
            }
            if(!over.connected(to) || over.distance(to) > min_dist){
                if(outer->value(~l) == l_True){
                    //do nothing
                }else if(outer->value(~l) == l_Undef){
                    outer->enqueue(~l, weighted_overprop_marker);
                }else if(outer->value(~l) == l_False){

                    return false;
                }
            }
        }else{
            if(under.connected(to) && under.distance(to) < min_dist){
                if(outer->value(l) == l_True){
                    //do nothing
                }else if(outer->value(l) == l_Undef){
                    outer->enqueue(l, weighted_underprop_marker);
                }else if(outer->value(l) == l_False){
                    return false;
                }
            }
            if(!over.connected(to) || over.distance(to) >= min_dist){
                if(outer->value(~l) == l_True){
                    //do nothing
                }else if(outer->value(~l) == l_Undef){
                    outer->enqueue(~l, weighted_overprop_marker);
                }else if(outer->value(~l) == l_False){

                    return false;
                }
            }
        }
    }



    //	}
    return true;
}


template<typename Weight, typename Graph>
Lit WeightedDistanceDetector<Weight, Graph>::decide(CRef& decision_reason){
    if(!opt_decide_graph_distance || !overapprox_weighted_distance_detector)
        return lit_Undef;
    WeightedDistanceDetector* r = this;
    auto* over = r->overapprox_weighted_distance_detector;

    auto* under = r->underapprox_weighted_distance_detector;

    //we can probably also do something similar, but with cuts, for nodes that are decided to be unreachable.

    //ok, for each node that is assigned reachable, but that is not actually reachable in the under approx, decide an edge on a feasible path

    //this can be obviously more efficient
    for(int type = 0; type < 2; type++){

        for(int k = 0; k < (type ? weighted_dist_bv_lits.size() : weighted_dist_lits.size()); k++){

            Lit l;
            if(!type)
                l = weighted_dist_lits[k].l;
            else
                l = weighted_dist_bv_lits[k].l;

            if(l == lit_Undef)
                continue;

            Weight min_dist;
            if(!type)
                min_dist = weighted_dist_lits[k].min_distance;
            else{
                min_dist = weighted_dist_bv_lits[k].bv.getUnder();
            }


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

                    if(over->distance(j) <= min_dist && under->distance(j) > min_dist){
                        //then lets try to disconnect this node from source by walking back aint64_t the path in the over approx, and disabling the first unassigned edge we see.
                        //(there must be at least one such edge, else the variable would be connected in the under approximation as well - in which case it would already have been propagated.
                        int p = j;
                        int last = j;
                        int dist = 0;
                        while(under->distance(p) > min_dist - dist){
                            last = p;
                            assert(p != source);
                            int prev = over->previous(p);
                            int edgeid = over->incomingEdge(p);

                            Var v = outer->getEdgeVar(edgeid);
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

template
class Monosat::WeightedDistanceDetector<int>;

template
class Monosat::WeightedDistanceDetector<int64_t>;

template
class Monosat::WeightedDistanceDetector<double>;

template
class Monosat::WeightedDistanceDetector<mpq_class>;
