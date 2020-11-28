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

#include "GraphTheory.h"
#include "monosat/dgl/SteinerApprox.h"
#include "monosat/dgl/FloydWarshall.h"

using namespace Monosat;

template<typename Weight, typename Graph>
SteinerDetector<Weight, Graph>::SteinerDetector(int detectorID, GraphTheorySolver<Weight>* outer,
                                                Graph& g, Graph& antig, double seed) :
        Detector(detectorID), outer(outer), g_under(g), g_over(antig), rnd_seed(seed){
    checked_unique = false;
    all_unique = true;
    positiveStatus = new SteinerDetector<Weight, Graph>::SteinerStatus(*this, true);
    negativeStatus = new SteinerDetector<Weight, Graph>::SteinerStatus(*this, false);

    //NOTE: the terminal sets are intentionally swapped, in order to preserve monotonicity
    underapprox_detector = new SteinerApprox<Weight, Graph, DynamicNodes, SteinerDetector<Weight, Graph>::SteinerStatus>(
            g,
            overTerminalSet, *positiveStatus,
            1); //new SpiraPan<SteinerDetector<Weight,Graph>::MSTStatus>(_g,*(positiveReachStatus),1);
    overapprox_detector = new SteinerApprox<Weight, Graph, DynamicNodes, SteinerDetector<Weight, Graph>::SteinerStatus>(
            antig,
            underTerminalSet, *negativeStatus, -1);

    underprop_marker = outer->newReasonMarker(getID());
    overprop_marker = outer->newReasonMarker(getID());

    underprop_edge_marker = outer->newReasonMarker(getID());
    overprop_edge_marker = outer->newReasonMarker(getID());
    first_reach_var = var_Undef;
}

template<typename Weight, typename Graph>
void SteinerDetector<Weight, Graph>::addTerminalNode(int node, Var outer_Var){
    Var var = outer->newVar(outer_Var, getID());
    underTerminalSet.addNode(node);
    underTerminalSet.setNodeEnabled(node, false);
    overTerminalSet.addNode(node);
    overTerminalSet.setNodeEnabled(node, true);
    terminal_map.growTo(g_under.nodes(), var_Undef);
    terminal_map[node] = var;
    terminal_var_map.growTo(var + 1, -1);
    terminal_var_map[var] = node;
}

template<typename Weight, typename Graph>
void SteinerDetector<Weight, Graph>::addWeightLit(Weight& min_weight, Var outer_weight_var){
    g_under.invalidate();
    g_over.invalidate();

    Var weight_var = outer->newVar(outer_weight_var, getID());

    Lit reachLit = mkLit(weight_var, false);
    bool found = false;
    for(int i = 0; i < weight_lits.size(); i++){
        if(weight_lits[i].min_weight == min_weight){
            found = true;
            Lit r = weight_lits[i].l;
            //force equality between the new lit and the old reach lit, in the SAT solver
            outer->makeEqual(r, reachLit);

        }
    }
    if(!found){
        weight_lits.push();
        weight_lits.last().l = reachLit;
        weight_lits.last().min_weight = min_weight;

    }

}

template<typename Weight, typename Graph>
void SteinerDetector<Weight, Graph>::SteinerStatus::setMinimumSteinerTree(Weight& weight){

    for(int i = 0; i < detector.weight_lits.size(); i++){
        Weight& min_weight = detector.weight_lits[i].min_weight;
        Lit l = detector.weight_lits[i].l;
        if(l != lit_Undef){
            assert(l != lit_Undef);
            if(min_weight < weight && !polarity){
                lbool assign = detector.outer->value(l);
                if(assign != l_False){
                    detector.changed_weights.push({~l, min_weight});
                }
            }else if(min_weight >= weight && polarity){
                lbool assign = detector.outer->value(l);
                if(assign != l_True){
                    detector.changed_weights.push({l, min_weight});
                }
            }
        }
    }

}

template<typename Weight, typename Graph>
void SteinerDetector<Weight, Graph>::buildMinWeightTooSmallReason(Weight& weight, vec<Lit>& conflict){
    //if the weight is too small, then either an edge has to be enabled, or a terminal node that is currently disabled has to be enabled.
    for(int i = 0; i < underTerminalSet.nodes(); i++){
        if(!underTerminalSet.nodeEnabled(i) && terminal_map[i] != var_Undef){
            conflict.push(mkLit(terminal_map[i], false));
        }
    }
    std::vector<int> edges;
    underapprox_detector->getSteinerTree(edges);
    for(int edgeID : edges){
        assert(g_under.edgeEnabled(edgeID));
        conflict.push(mkLit(outer->getEdgeVar(edgeID), true));
    }

}

template<typename Weight, typename Graph>
void SteinerDetector<Weight, Graph>::buildMinWeightTooLargeReason(Weight& weight, vec<Lit>& conflict){
    assert(overapprox_detector->disconnected() || overapprox_detector->weight() > weight);
    //if the weight is too large, then either an edge has to be enabled, or a terminal node that is currently enabled has to be disabled.
    for(int i = 0; i < overTerminalSet.nodes(); i++){
        if(overTerminalSet.nodeEnabled(i) && terminal_map[i] != var_Undef){
            conflict.push(mkLit(terminal_map[i], true));
        }
    }

    if(overapprox_detector->disconnected()){
        //Then find a separating cut between any two terminals.

        double starttime = rtime(2);
        int INF = std::numeric_limits<int>::max();

        //IF the mst is disconnected, then we define it's weight to be infinite. In this case, the reason is a separating cut between any two disconnected components.
        //we can find one of these by identifying any two roots

        int min_conflict = INF;

        //walk back down from the each root to find a separating cut of disabled edge.
        //return the smallest such cut.

        DisjointSets sets;
        sets.AddElements(g_over.nodes());

        for(int i = 0; i < g_over.edges(); i++){
            if(g_over.edgeEnabled(i)){
                int u = g_over.getEdge(i).from;
                int v = g_over.getEdge(i).to;
                sets.UnionElements(u, v);
            }
        }
        assert(sets.NumSets() > 1);
        vec<bool> visited;
        visited.growTo(g_under.nodes());

        for(int i = 0; i < underTerminalSet.nodes(); i++){
            if(!underTerminalSet.nodeEnabled(i) || visited[i])
                continue;
            visited[i] = true;

            int root = i;
            int set = sets.FindSet(root);
            tmp_conflict.clear();
            visit.clear();
            //ok, now traverse the connected component in the tree below this root.
            visit.push(root);
            seen.clear();

            seen.growTo(g_under.nodes());
            seen[root] = true;
            while(visit.size()){
                int u = visit.last();
                visit.pop();
                assert(u > -1);
                for(int i = 0; i < g_over.nIncident(u, true); i++){
                    int edgeid = g_over.incident(u, i, true).id;
                    if(g_over.edgeEnabled(edgeid)){
                        int v = g_over.incident(u, i, true).node;
                        assert(sets.FindSet(v) == set);
                        if(!seen[v]){
                            //u is a child of node in the minimum spanning tree
                            seen[v] = true;
                            visit.push(v);
                        }
                    }else if(!g_over.edgeEnabled(edgeid)){
                        int v = g_over.incident(u, i, true).node;
                        int s = sets.FindSet(v);
                        if(s !=
                           set){    //&& hasTerminal[s] this isn't correct, because there can be a separating component between two terminals that has no terminal itself, and yet an edge has to be added to it to make the terminals connected.
                            //this node is on a separating cut.
                            Var e = outer->getEdgeVar(edgeid);
                            assert(outer->value(e) == l_False);
                            tmp_conflict.push(mkLit(e, false));
                        }
                    }
                }
            }

            if(tmp_conflict.size() < min_conflict){
                min_conflict = tmp_conflict.size();
                conflict.shrink(conflict.size() - 1);//keep only the first conflict element, which is the constraint lit
                assert(conflict.size() == 1);
                for(int i = 0; i < tmp_conflict.size(); i++){
                    conflict.push(tmp_conflict[i]);
                }
            }
            if(!opt_mst_min_cut){
                break;                //stop looking as soon as we find a cut, instead of trying to find the smallest one.
            }
        }

        return;

    }

    //all pairs shortest paths
    FloydWarshall<Weight> fw(g_over);

    for(int i = 0; i < g_over.edges(); i++){
        if(g_over.isEdge(i) && !g_over.edgeEnabled(i)){
            int u = g_over.getEdge(i).from;
            int v = g_over.getEdge(i).to;
            if(g_over.getWeight(i) < fw.distance(u, v)){
                //then adding this edge can potentially decrease the steiner tree
                //(there is probably an additional test we can do to include fewer edges...)
                conflict.push(mkLit(outer->getEdgeVar(i), false));
            }
        }
    }
}

template<typename Weight, typename Graph>
void SteinerDetector<Weight, Graph>::buildReason(Lit p, vec<Lit>& reason, CRef marker){

    if(marker == underprop_marker){
        reason.push(p);

        Var v = var(p);
        Weight weight = -1;
        //could swap this out for a map if there are lots of lits..
        for(int i = 0; i < weight_lits.size(); i++){
            if(var(weight_lits[i].l) == v){
                weight = weight_lits[i].min_weight;
                break;
            }
        }
        assert(weight >= 0);
        buildMinWeightTooSmallReason(weight, reason);

    }else if(marker == overprop_marker){
        reason.push(p);

        //the reason is a cut separating p from s;
        //We want to find a min-cut in the full graph separating, where activated edges (ie, those still in antig) are weighted infinity, and all others are weighted 1.

        //This is a cut that describes a minimal set of edges which are disabled in the current graph, at least one of which would need to be activated in order for s to reach p
        //assign the mincut edge weights if they aren't already assigned.

        Var v = var(p);

        Weight weight = -1;
        //could swap this out for a map if there are lots of lits..
        for(int i = 0; i < weight_lits.size(); i++){
            if(var(weight_lits[i].l) == v){
                weight = weight_lits[i].min_weight;
                break;
            }
        }
        assert(weight >= 0);

        buildMinWeightTooLargeReason(weight, reason);

    }else{
        assert(false);
    }
    outer->toSolver(reason);
}

template<typename Weight, typename Graph>
bool SteinerDetector<Weight, Graph>::propagate(vec<Lit>& conflict){
    changed_weights.clear();
    double startdreachtime = rtime(2);
    stats_under_updates++;
    underapprox_detector->update();
    double reachUpdateElapsed = rtime(2) - startdreachtime;
    stats_under_update_time += reachUpdateElapsed;

    double startunreachtime = rtime(2);
    stats_over_updates++;
    overapprox_detector->update();
    double unreachUpdateElapsed = rtime(2) - startunreachtime;
    stats_over_update_time += unreachUpdateElapsed;

    for(int j = 0; j < changed_weights.size(); j++){
        Lit l = changed_weights[j].l;
        Weight weight = changed_weights[j].weight;

        bool reach = !sign(l);
        if(outer->value(l) == l_True){
            //do nothing
        }else if(outer->value(l) == l_Undef){
            if(reach)
                outer->enqueue(l, underprop_marker);
            else
                outer->enqueue(l, overprop_marker);

        }else if(outer->value(l) == l_False){
            conflict.push(l);

            if(reach){

                //conflict
                //The reason is a path in g from to s in d
                buildMinWeightTooSmallReason(weight, conflict);
                //add it to s
                //return it as a conflict

            }else{

                buildMinWeightTooLargeReason(weight, conflict);

            }
            outer->toSolver(conflict);
            return false;
        }

    }

    return true;
}

template<typename Weight, typename Graph>
bool SteinerDetector<Weight, Graph>::checkSatisfied(){

    assert(underTerminalSet.numEnabled() == overTerminalSet.numEnabled());

    SteinerApprox<Weight, Graph, DynamicNodes, typename SteinerTree<Weight>::NullStatus> positive_checker(g_under,
                                                                                                          overTerminalSet,
                                                                                                          SteinerTree<Weight>::nullStatus,
                                                                                                          0);
    SteinerApprox<Weight, Graph, DynamicNodes, typename SteinerTree<Weight>::NullStatus> negative_checker(g_over,
                                                                                                          underTerminalSet,
                                                                                                          SteinerTree<Weight>::nullStatus,
                                                                                                          0);
    positive_checker.update();
    negative_checker.update();
    Weight& w = positive_checker.weight();
    for(int k = 0; k < weight_lits.size(); k++){
        Lit l = weight_lits[k].l;
        Weight dist = weight_lits[k].min_weight;

        if(l != lit_Undef){

            if(outer->value(l) == l_True){
                if(positive_checker.weight() > dist){
                    return false;
                }
            }else if(outer->value(l) == l_False){
                if(negative_checker.weight() <= dist){
                    return false;
                }
            }else{
                if(positive_checker.weight() <= dist){
                    return false;
                }
                if(negative_checker.weight() > dist){
                    return false;
                }
            }
        }
    }

    return true;
}

template<typename Weight, typename Graph>
Lit SteinerDetector<Weight, Graph>::decide(CRef& decision_reason){

    return lit_Undef;
};

template
class Monosat::SteinerDetector<int>;

template
class Monosat::SteinerDetector<int64_t>;

template
class Monosat::SteinerDetector<double>;

template
class Monosat::SteinerDetector<mpq_class>;
