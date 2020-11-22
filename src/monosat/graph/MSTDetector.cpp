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
#include "monosat/dgl/Prim.h"
#include "monosat/dgl/SpiraPan.h"

#include <set>

using namespace Monosat;

template<typename Weight, typename Graph>
MSTDetector<Weight, Graph>::MSTDetector(int detectorID, GraphTheorySolver<Weight>* outer, Graph& g,
                                        Graph& antig, double seed) :
        Detector(detectorID), outer(outer), g_under(g), g_over(antig), rnd_seed(seed){
    checked_unique = false;
    all_unique = true;
    positiveReachStatus = new MSTDetector<Weight, Graph>::MSTStatus(*this, true);
    negativeReachStatus = new MSTDetector<Weight, Graph>::MSTStatus(*this, false);

    if(mstalg == MinSpanAlg::ALG_KRUSKAL){
        underapprox_detector = new Kruskal<MSTDetector<Weight, Graph>::MSTStatus, Weight>(g,
                                                                                          *(positiveReachStatus), 1);
        overapprox_detector = new Kruskal<MSTDetector<Weight, Graph>::MSTStatus, Weight>(antig,
                                                                                         *(negativeReachStatus), -1);
        underapprox_conflict_detector = underapprox_detector;
        overapprox_conflict_detector = overapprox_detector;
    }else if(mstalg == MinSpanAlg::ALG_PRIM){
        underapprox_detector = new Prim<MSTDetector<Weight, Graph>::MSTStatus, Weight>(g, *(positiveReachStatus),
                                                                                       1);
        overapprox_detector = new Prim<MSTDetector<Weight, Graph>::MSTStatus, Weight>(antig,
                                                                                      *(negativeReachStatus), -1);
        underapprox_conflict_detector = new Kruskal<typename MinimumSpanningTree<Weight>::NullStatus, Weight>(g,
                                                                                                              MinimumSpanningTree<Weight>::nullStatus,
                                                                                                              1);
        overapprox_conflict_detector = new Kruskal<typename MinimumSpanningTree<Weight>::NullStatus, Weight>(antig,
                                                                                                             MinimumSpanningTree<Weight>::nullStatus,
                                                                                                             -1);

    }else if(mstalg == MinSpanAlg::ALG_SPIRA_PAN){

        underapprox_detector = new SpiraPan<MSTDetector<Weight, Graph>::MSTStatus, Weight>(g,
                                                                                           *(positiveReachStatus),
                                                                                           1); //new SpiraPan<MSTDetector<Weight,Graph>::MSTStatus>(_g,*(positiveReachStatus),1);
        overapprox_detector = new SpiraPan<MSTDetector<Weight, Graph>::MSTStatus, Weight>(antig,
                                                                                          *(negativeReachStatus), -1);
        underapprox_conflict_detector = new Kruskal<typename MinimumSpanningTree<Weight>::NullStatus, Weight>(g,
                                                                                                              MinimumSpanningTree<Weight>::nullStatus,
                                                                                                              1);
        overapprox_conflict_detector = new Kruskal<typename MinimumSpanningTree<Weight>::NullStatus, Weight>(antig,
                                                                                                             MinimumSpanningTree<Weight>::nullStatus,
                                                                                                             -1);
    }

    underprop_marker = outer->newReasonMarker(getID());
    overprop_marker = outer->newReasonMarker(getID());

    underprop_edge_marker = outer->newReasonMarker(getID());
    overprop_edge_marker = outer->newReasonMarker(getID());
    first_reach_var = var_Undef;
}

template<typename Weight, typename Graph>
void MSTDetector<Weight, Graph>::addWeightLit(Var outer_weight_var, Weight& min_weight, bool inclusive){
    g_under.invalidate();
    g_over.invalidate();

    Var weight_var = outer->newVar(outer_weight_var, getID());
    if(lowest_weight_lit < 0 || min_weight < lowest_weight_lit){
        lowest_weight_lit = min_weight;
    }
    if(min_weight > highest_weight_lit){
        highest_weight_lit = min_weight;
    }
    Lit reachLit = mkLit(weight_var, false);
    bool found = false;
    for(int i = 0; i < weight_lits.size(); i++){
        if(weight_lits[i].min_weight == min_weight && weight_lits[i].inclusive == inclusive){
            found = true;
            Lit r = weight_lits[i].l;
            //force equality between the new lit and the old reach lit, in the SAT solver
            outer->makeEqual(reachLit, r, true);
            /*outer->S->addClause(~r, reachLit);
             outer->S->addClause(r, ~reachLit);*/
        }
    }
    if(!found){
        weight_lits.push();
        weight_lits.last().l = reachLit;
        weight_lits.last().inclusive = inclusive;
        weight_lits.last().min_weight = min_weight;
    }
}

template<typename Weight, typename Graph>
void MSTDetector<Weight, Graph>::addTreeEdgeLit(int edge_id, Var outer_reach_var){
    g_under.invalidate();
    g_over.invalidate();
    Var reach_var = outer->newVar(outer_reach_var, getID());
    /*while(outer->S->nVars()<=reach_var)
     outer->S->newVar();*/

    if(!checked_unique){
        checked_unique = true;

        std::set<Weight> seen;
        for(int i = 0; i < g_over.getWeights().size(); i++){
            Weight w = g_over.getWeight(i);
            if(seen.count(w) > 0){
                all_unique = false;
                std::cout
                        << "Minimum spanning tree edge constraints can only be enforced in graphs in which all edges have unique weights.\nTree edge constraints will be left unconstrained because multiple edges have weight "
                        << w << "\n";
                return;
            }else{
                seen.insert(w);
            }
        }

    }
    if(!all_unique)
        return;
    tree_edge_lits.growTo(outer->nEdges());

    //while( dist_lits[to].size()<=within_steps)
    //	dist_lits[to].push({lit_Undef,-1});

    Lit reachLit = mkLit(reach_var, false);
    bool found = false;

    if(tree_edge_lits[edge_id].l != lit_Undef){
        found = true;
        Lit r = tree_edge_lits[edge_id].l;
        //force equality between the new lit and the old reach lit, in the SAT solver
        outer->makeEqual(reachLit, r, true);
        /*		outer->S->addClause(~r, reachLit);
         outer->S->addClause(r, ~reachLit);*/
    }else{
        if(first_reach_var == var_Undef){
            first_reach_var = reach_var;
        }
        assert(first_reach_var <= reach_var);
        tree_edge_lits[edge_id].l = reachLit;
        tree_edge_lits[edge_id].edgeID = edge_id;

        Var edgeVar = outer->getEdgeVar(edge_id);
        Lit edgeEnabled = mkLit(edgeVar, false);
        outer->addClause(edgeEnabled,
                         reachLit);//If the edge is not enabled, then we artificially enforce that the edge counts as being in the tree.
        //this is required to enforce monotonicity of the in the tree constraint.
        while(tree_edge_lits_map.size() <= reach_var - first_reach_var){
            tree_edge_lits_map.push(-1);
        }

        tree_edge_lits_map[reach_var - first_reach_var] = edge_id;
    }

}

template<typename Weight, typename Graph>
void MSTDetector<Weight, Graph>::MSTStatus::inMinimumSpanningTree(int edgeid, bool in_tree){
    if(edgeid < detector.tree_edge_lits.size()){
        Lit l = detector.tree_edge_lits[edgeid].l;
        //Note: for the tree edge detector, polarity is effectively reversed.
        if(l != lit_Undef){
            if(!polarity && in_tree){
                if(!detector.is_edge_changed[edgeid]){
                    detector.is_edge_changed[edgeid] = true;
                    detector.changed_edges.push({var(l), edgeid});
                }
            }else if(polarity && !in_tree){
                if(!detector.is_edge_changed[edgeid]){
                    detector.is_edge_changed[edgeid] = true;
                    //Var edgevar = detector.outer->getEdgeVar(edgeid);
                    //assert(detector.outer->value(edgevar)!=l_False);		//else the edge counts as in the tree
                    detector.changed_edges.push({var(l), edgeid});
                }
            }
        }
    }
}

template<typename Weight, typename Graph>
void MSTDetector<Weight, Graph>::MSTStatus::setMinimumSpanningTree(Weight& weight, bool connected){

    /*for(int i = 0;i<detector.weight_lits.size();i++){
     Weight & min_weight =  detector.weight_lits[i].min_weight;
     Lit l = detector.weight_lits[i].l;
     if(l!=lit_Undef){
     assert(l!=lit_Undef);
     if((!connected ||   weight>min_weight) && !polarity){
     lbool assign = detector.outer->value(l);
     if( assign!= l_False ){
     detector.changed_weights.push({~l,min_weight});
     }
     }else if((connected && weight<=min_weight) && polarity){
     lbool assign = detector.outer->value(l);
     if( assign!= l_True ){
     detector.changed_weights.push({l,min_weight});
     }
     }
     }
     }*/

}

template<typename Weight, typename Graph>
void MSTDetector<Weight, Graph>::buildMinWeightTooSmallReason(Weight& weight, vec<Lit>& conflict){

    MinimumSpanningTree<Weight>& d = *underapprox_conflict_detector;

    double starttime = rtime(2);
    d.update();
    Weight& minweight = d.weight();
    assert(d.weight() <= weight);
    //learn that at least one edge in the tree must be disabled (else, the minimum weight cannot be higher than the current weight)
    std::vector<int>& mst = d.getSpanningTree();
    for(int i = 0; i < mst.size(); i++){
        int edgeid = mst[i];
        Var v = outer->getEdgeVar(edgeid);
        assert(outer->value(v) == l_True);
        conflict.push(mkLit(v, true));
    }

    outer->num_learnt_paths++;
    outer->learnt_path_clause_length += (conflict.size() - 1);
    double elapsed = rtime(2) - starttime;
    outer->pathtime += elapsed;

}

template<typename Weight, typename Graph>
bool MSTDetector<Weight, Graph>::walkback(Weight& min_weight, int from, int to){
    int u = from;
    while(u != to && u != -1){
        int p = overapprox_conflict_detector->getParent(u);
        if(p != -1){
            int edgeid = overapprox_conflict_detector->getParentEdge(u);
            Weight weight = g_over.getWeight(edgeid);
            if(weight > min_weight){
                return true;
            }
        }
        u = p;
    }
    return false;
}

template<typename Weight, typename Graph>
void MSTDetector<Weight, Graph>::TarjanOLCA(int node, vec<Lit>& conflict){
    ancestors[node] = node;
    for(int i = 0; i < g_over.nIncident(node, true); i++){
        int edgeid = g_over.incident(node, i, true).id;
        if(overapprox_conflict_detector->edgeInTree(edgeid)){
            assert(g_over.edgeEnabled(edgeid));
            int v = g_over.incident(node, i, true).node;
            if(overapprox_conflict_detector->getParent(v) == node){
                //u is a child of node in the minimum spanning tree
                TarjanOLCA(v, conflict);
                sets.UnionElements(node, v);
                int set = sets.FindSet(node);
                ancestors[set] = node;

            }
        }
    }
    black[node] = true;
    //now visit all _disabled_ edges of node
    for(int i = 0; i < g_over.nIncident(node, true); i++){
        int edgeid = g_over.incident(node, i, true).id;
        if(!g_over.edgeEnabled(edgeid)){
            //this is a disabled edge
            int v = g_under.incident(node, i, true).node;
            if(black[v]){
                int set = sets.FindSet(v);
                int lowest_common_ancestor = ancestors[set];

                //ok, now walk back from u and v in the minimum spanning tree until either the lowest common ancestor is seen, or an edge larger than the weight of this disabled edge is found
                Weight weight = g_over.getWeight(edgeid);
                bool any_larger_weights = walkback(weight, node, lowest_common_ancestor)
                                          || walkback(weight, v, lowest_common_ancestor);
                //if any larger edge was found in either path from u or v to their common ancestor in the minimum spanning tree, then enabling this edge
                //would have replaced that larger edge in the minimum spanning tree, resulting in a smaller minimum spanning tree.
                if(any_larger_weights){
                    Var e = outer->getEdgeVar(edgeid);
                    assert(outer->value(e) == l_False);
                    conflict.push(mkLit(e, false));
                }
            }
        }
    }
}

template<typename Weight, typename Graph>
void MSTDetector<Weight, Graph>::buildMinWeightTooLargeReason(Weight& weight, vec<Lit>& conflict){

    double starttime = rtime(2);

    overapprox_conflict_detector->update();

    Weight& mstweight = overapprox_conflict_detector->weight();

    if(overapprox_conflict_detector->numComponents() > 1){
        //IF the mst is disconnected, then we define it's weight to be infinite. In this case, the reason is a separating cut between any two disconnected components.
        //we can find one of these by identifying any two roots

        int min_conflict = -1;

        //walk back down from the each root to find a separating cut of disabled edge.
        //return the smallest such cut.

        assert(conflict.size() == 1);
        for(int i = 0; i < overapprox_conflict_detector->numComponents(); i++){
            int root = overapprox_conflict_detector->getRoot(i);
            tmp_conflict.clear();
            int set = overapprox_conflict_detector->getComponent(root); //sets.FindSet(root);
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
                        assert(overapprox_conflict_detector->getComponent(v) == set);
                        if(!seen[v]){
                            //u is a child of node in the minimum spanning tree
                            seen[v] = true;
                            visit.push(v);
                        }
                    }else if(!g_over.edgeEnabled(edgeid)){
                        int v = g_over.incident(u, i, true).node;
                        if(overapprox_conflict_detector->getComponent(v) != set){
                            //Then this edge is on the cut between this component and the other components
                            Var e = outer->getEdgeVar(edgeid);
                            assert(outer->value(e) == l_False);
                            tmp_conflict.push(mkLit(e, false));
                        }
                    }
                }
            }

            if(min_conflict == -1 || tmp_conflict.size() < min_conflict){
                min_conflict = tmp_conflict.size();
                conflict.shrink(
                        conflict.size() - 1); //keep only the first conflict element, which is the constraint lit
                assert(conflict.size() == 1);
                for(int i = 0; i < tmp_conflict.size(); i++){
                    conflict.push(tmp_conflict[i]);
                }
            }
            if(!opt_mst_min_cut){
                break; //stop looking as soon as we find a cut, instead of trying to find the smallest one.
            }
        }

        return;
    }

    ancestors.clear();
    ancestors.growTo(g_under.nodes(), -1);
    black.clear();
    black.growTo(g_under.nodes());
    sets.Reset();
    sets.AddElements(g_under.nodes());

    //Noah's algorithm for finding a minimal necessary set of edges to add to decrease the weight of the tree
    //all we are doing here is visiting each disabled edge, then examining the cycle that would be created in the minimum spanning tree if we were to
    //add that edge to the tree. If any edge in that cycle is larger than the disabled edge, then enabling that edge would result in a smaller minimum spanning tree.
    //so we have to add that edge to the conflict as part of the reason that the minimum spanning tree is larger than it should be.

    //conceptually, thats really simple. The problem is that visiting that cycle requires computing the lowest common ancestor of the two nodes on either side of the disabled edge,
    //and to make that efficient we are going to use Tarjan's offline lca algorithm. This results in the convoluted code below.

    int root = 0;
    while(overapprox_conflict_detector->getParent(root) != -1){
        root = overapprox_conflict_detector->getParent(root);
    }
    TarjanOLCA(root, conflict);
#ifdef DEBUG_GRAPH
    for (int i = 0; i < black.size(); i++)
        assert(black[i]);
#endif

    outer->num_learnt_cuts++;
    outer->learnt_cut_clause_length += (conflict.size() - 1);

    double elapsed = rtime(2) - starttime;
    outer->mctime += elapsed;

}

template<typename Weight, typename Graph>
Weight MSTDetector<Weight, Graph>::walkback_edge(Weight& min_weight, int edge_id, int from, int to, bool& found){
    int u = from;
    Weight maxweight = 0;
    while(u != to && u != -1){
        int p = overapprox_conflict_detector->getParent(u);
        if(p != -1){
            int edgeid = overapprox_conflict_detector->getParentEdge(u);
            Weight weight = g_over.getWeight(edgeid);
            if(edgeid == edge_id){
                assert(!found);    //can't enounter the edge twice while traversing the cycle
                found = true;
            }
            if(weight > maxweight){
                maxweight = weight;
            }
        }
        u = p;
    }
    return maxweight;
}

template<typename Weight, typename Graph>
void MSTDetector<Weight, Graph>::TarjanOLCA_edge(int node, int check_edgeid, int lowest_endpoint, vec<Lit>& conflict){
    ancestors[node] = node;
    for(int i = 0; i < g_over.nIncident(node, true); i++){
        int edgeid = g_over.incident(node, i, true).id;
        if(!g_over.edgeEnabled(edgeid))
            continue;
        if(overapprox_conflict_detector->edgeInTree(edgeid)){
            int v = g_over.incident(node, i, true).node;
            if(overapprox_conflict_detector->getParent(v) == node){
                //u is a child of node in the minimum spanning tree
                TarjanOLCA_edge(v, check_edgeid, lowest_endpoint, conflict);
                sets.UnionElements(node, v);
                int set = sets.FindSet(node);
                assert(sets.FindSet(v) == set);

                ancestors[set] = node;

            }
        }
    }
    black[node] = true;
    //now visit all _disabled_ edges of node
    for(int i = 0; i < g_under.nIncident(node, true); i++){
        int edgeid = g_under.incident(node, i, true).id;
        if(!g_over.edgeEnabled(edgeid)){

            int v = g_under.incident(node, i, true).node;
            if(black[v]){
                int set = sets.FindSet(v);
                int lowest_common_ancestor = ancestors[set];

                //ok, now walk back from u and v in the minimum spanning tree until either the lowest common ancestor is seen, or an edge larger than the weight of this disabled edge is found
                Weight weight = g_under.getWeight(edgeid);
                bool found = false;

                Weight largest_weight = walkback_edge(weight, check_edgeid, node, lowest_common_ancestor, found);

                Weight largest_weight_other = walkback_edge(weight, check_edgeid, v, lowest_common_ancestor, found);
                if(found && (largest_weight > weight || largest_weight_other > weight)){
                    //if any edge larger than the disabled edge's weight was found in either path from u or v to their common ancestor in the minimum spanning tree, then enabling this edge
                    //would lead  to replacing that larger edge in the minimum spanning tree, resulting in a smaller minimum spanning tree, possibly resulting in the edge we really care about being removed from the tree.

                    Var e = outer->getEdgeVar(edgeid);
                    assert(outer->value(e) == l_False);
                    conflict.push(mkLit(e, false));

                }
            }
        }
    }
}

template<typename Weight, typename Graph>
void MSTDetector<Weight, Graph>::buildEdgeNotInTreeReason(int edgeid, vec<Lit>& conflict){
    Var edgevar = outer->getEdgeVar(edgeid);
    assert(outer->value(edgevar) != l_False);                            //else the edge counts as in the tree
    //what about if the mst is disconnected?
    //the reason that an edge is NOT in the minimum spanning tree is the paths to lca from either edge. So long as all those edges are in the tree, each of which is <= weight to this edge,
    //this edge cannot be in the mst.
    seen.clear();
    seen.growTo(g_under.nodes());
    int u = g_under.getEdge(edgeid).from;
    int v = g_under.getEdge(edgeid).to;
    underapprox_conflict_detector->update();
    //assert(!g.edgeEnabled(edgeid));
    assert(!underapprox_conflict_detector->edgeInTree(edgeid));
    int r = u;
    while(r != -1){
        seen[r] = true;
        r = underapprox_conflict_detector->getParent(r);

    }
    r = v;
    while(!seen[r] && r != -1){
        r = underapprox_conflict_detector->getParent(r);

    }
    if(r == -1){
        //then the mst is disconnected, but we can pretend it is connected with an edge of infinite weight directly to the root node of the other component.
        assert(underapprox_conflict_detector->numComponents() > 1);
        assert(underapprox_conflict_detector->getComponent(u) != underapprox_conflict_detector->getComponent(v));
    }
    int lca = r;
    r = u;
    while(r != lca){
        assert(seen[r]);
        seen[r] = false;
        int p = underapprox_conflict_detector->getParent(r);
        if(p != -1){
            int edge = underapprox_conflict_detector->getParentEdge(r);
            assert(g_under.edgeEnabled(edge));
            Var v = outer->getEdgeVar(edge);
            assert(outer->value(v) == l_True);
            conflict.push(mkLit(v, true));
        }
        r = p;
    }

    r = v;
    while(r != lca){
        assert(!seen[r]);
        int p = underapprox_conflict_detector->getParent(r);
        if(p != -1){
            int edge = underapprox_conflict_detector->getParentEdge(r);
            assert(g_under.edgeEnabled(edge));
            Var v = outer->getEdgeVar(edge);
            assert(outer->value(v) == l_True);
            conflict.push(mkLit(v, true));
        }
        r = p;
    }

    r = lca;
    while(r != -1){
        assert(seen[r] || r == lca);
        seen[r] = false;
        r = underapprox_conflict_detector->getParent(r);
    }

}

template<typename Weight, typename Graph>
void MSTDetector<Weight, Graph>::buildEdgeInTreeReason(int edgeid, vec<Lit>& conflict){
    //if the edge is disabled, then the reason for the edge being in the tree is that we have defined disabled edges to be in the mst.
    if(!g_over.edgeEnabled(edgeid)){
        Var v = outer->getEdgeVar(edgeid);
        assert(outer->value(v) == l_False);
        conflict.push(mkLit(v, false));
        return;
    }
    Var vt = outer->getEdgeVar(edgeid);
    assert(vt > 0);

    ancestors.clear();
    ancestors.growTo(g_under.nodes(), -1);

    black.clear();
    black.growTo(g_under.nodes());
    double starttime = rtime(2);
    int INF = std::numeric_limits<int>::max();
    overapprox_conflict_detector->update();
    assert(overapprox_conflict_detector->edgeInTree(edgeid));
    sets.Reset();
    sets.AddElements(g_under.nodes());

    //Noah's algorithm for finding a minimal necessary set of edges to add to decrease the weight of the tree
    //all we are doing here is visiting each disabled edge, then examining the cycle that would be created in the minimum spanning tree if we were to
    //add that edge to the tree. If any edge in that cycle is larger than the disabled edge, then enabling that edge would result in a smaller minimum spanning tree.
    //so we have to add that edge to the conflict as part of the reason that the minimum spanning tree is larger than it should be.

    //conceptually, thats really simple. The problem is that visiting that cycle requires computing the lowest common ancestor of the two nodes on either side of the disabled edge,
    //and to make that efficient we are going to use Tarjan's offline lca algorithm. This results in the convoluted code below.
    int u = g_under.getEdge(edgeid).from;
    int v = g_under.getEdge(edgeid).to;
    int lower_endpoint = u;
    if(overapprox_conflict_detector->getParent(v) == u){
        lower_endpoint = v;
    }else{
        assert(overapprox_conflict_detector->getParent(u) == v);
    }

    int root = v;
    while(overapprox_conflict_detector->getParent(root) != -1){
        root = overapprox_conflict_detector->getParent(root);
    }
    TarjanOLCA_edge(root, edgeid, lower_endpoint,
                    conflict);//run tarjan's off-line lowest common ancestor query from node 0, arbitrarily.


    outer->num_learnt_cuts++;
    outer->learnt_cut_clause_length += (conflict.size() - 1);

    double elapsed = rtime(2) - starttime;
    outer->mctime += elapsed;
}

template<typename Weight, typename Graph>
void MSTDetector<Weight, Graph>::buildReason(Lit p, vec<Lit>& reason, CRef marker){

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

        //double elapsed = rtime(2)-startpathtime;
        //	pathtime+=elapsed;
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

    }else if(marker == underprop_edge_marker){
        reason.push(p);

        Var v = var(p);
        assert(v >= first_reach_var);
        int edgeid = tree_edge_lits_map[v - first_reach_var];

        buildEdgeInTreeReason(edgeid, reason);

        //double elapsed = rtime(2)-startpathtime;
        //	pathtime+=elapsed;
    }else if(marker == overprop_edge_marker){
        reason.push(p);

        //the reason is a cut separating p from s;
        //We want to find a min-cut in the full graph separating, where activated edges (ie, those still in antig) are weighted infinity, and all others are weighted 1.

        //This is a cut that describes a minimal set of edges which are disabled in the current graph, at least one of which would need to be activated in order for s to reach p
        //assign the mincut edge weights if they aren't already assigned.

        Var v = var(p);
        assert(v >= first_reach_var);
        int edgeid = tree_edge_lits_map[v - first_reach_var];
        buildEdgeNotInTreeReason(edgeid, reason);

    }else{
        assert(false);
    }
    outer->toSolver(reason);
}


template<typename Weight, typename Graph>
void MSTDetector<Weight, Graph>::preprocess(){
    is_edge_changed.growTo(g_under.edges());
}

template<typename Weight, typename Graph>
bool MSTDetector<Weight, Graph>::propagate(vec<Lit>& conflict){
    if(outer->has_any_bitvector_edges){
        throw std::runtime_error("MST constraints don't yet support bitvector weight edges");
    }
//NOTE! Cannot use pure theory lits here, because the edge literals and the mst weight literals are computed
//in opposed polarity by each detector! So unless we separate those unassigned counts out, we can't skip either...
//Fix this later if needed...
//if(positive_reach_detector && (!opt_detect_pure_theory_lits || unassigned_positives>0)){
    double startdreachtime = rtime(2);
    stats_under_updates++;
    underapprox_detector->update();
    double reachUpdateElapsed = rtime(2) - startdreachtime;
    stats_under_update_time += reachUpdateElapsed;
//}else
//	stats_skipped_under_updates++;

//if(negative_reach_detector && (!opt_detect_pure_theory_lits || unassigned_negatives>0)){
    double startunreachtime = rtime(2);
    stats_over_updates++;
    overapprox_detector->update();
    double unreachUpdateElapsed = rtime(2) - startunreachtime;
    stats_over_update_time += unreachUpdateElapsed;
//}else
//	stats_skipped_over_updates++;

    Weight& under_weight = underapprox_detector->weight();
    Weight& over_weight = overapprox_detector->weight();
    if(weight_lits.size()
       && ((under_weight >= lowest_weight_lit || underapprox_detector->numComponents() > 1)
           || (over_weight <= highest_weight_lit && overapprox_detector->numComponents() <= 1))){
        //Probably should do a binary search here.
        for(int j = 0; j < weight_lits.size(); j++){
            Lit l = weight_lits[j].l;
            bool inclusive = weight_lits[j].inclusive;
            //printf("mst: %d\n",dimacs(l));
            Weight min_weight = weight_lits[j].min_weight;

            if(inclusive && under_weight <= min_weight && underapprox_detector->numComponents() <= 1){
                //dont flip sign
            }else if(!inclusive && under_weight < min_weight && underapprox_detector->numComponents() <= 1){
                //dont flip sign
            }else if(inclusive && (over_weight > min_weight || overapprox_detector->numComponents() > 1)){
                l = ~l;            //flip sign
            }else if(!inclusive && (over_weight >= min_weight || overapprox_detector->numComponents() > 1)){
                l = ~l;            //flip sign
            }else{
                //try the next edge weight
                continue;
            }

            bool reach = !sign(l);
            if(outer->value(l) == l_True){
                //do nothing
            }else if(outer->value(l) == l_Undef){
                //trail.push(Assignment(false,reach,detectorID,0,var(l)));
                if(reach)
                    outer->enqueue(l, underprop_marker);
                else
                    outer->enqueue(l, overprop_marker);

            }else if(outer->value(l) == l_False){
                conflict.push(l);

                if(reach){

                    //conflict
                    //The reason is a path in g from to s in d
                    buildMinWeightTooSmallReason(min_weight, conflict);
                    //add it to s
                    //return it as a conflict

                }else{

                    buildMinWeightTooLargeReason(min_weight, conflict);

                }
                outer->toSolver(conflict);
                return false;
            }

        }
    }

    while(changed_edges.size()){
        int sz = changed_edges.size();
        Var v = changed_edges.last().v;
        int edgeID = changed_edges.last().edgeID;
        assert(is_edge_changed[edgeID]);
        Lit l;
        Var edgevar = outer->getEdgeVar(edgeID);
        lbool edge_val = outer->value(edgevar);
        if((!g_over.edgeEnabled(edgeID) || overapprox_detector->edgeInTree(edgeID))){
            l = mkLit(v, false);
        }else if(g_under.edgeEnabled(edgeID) && !underapprox_detector->edgeInTree(edgeID)){

            assert(outer->value(edgevar) != l_False);            //else the edge counts as in the tree
            l = mkLit(v, true);
        }else{
            assert(changed_edges.size() == sz);
            assert(changed_edges.last().edgeID == edgeID);
            is_edge_changed[edgeID] = false;
            changed_edges.pop();
            //this can happen if the changed node's reachability status was reported before a backtrack in the solver.
            continue;
        }

        bool reach = !sign(l);
        if(outer->value(l) == l_True){
            //do nothing
        }else if(outer->value(l) == l_Undef){
            //trail.push(Assignment(false,reach,detectorID,0,var(l)));
            if(reach){
                outer->enqueue(l, underprop_edge_marker);
            }else
                outer->enqueue(l, overprop_edge_marker);

        }else if(outer->value(l) == l_False){
            conflict.push(l);

            if(reach){

                //conflict
                //The reason is a path in g from to s in d
                buildEdgeInTreeReason(edgeID, conflict);

                //add it to s
                //return it as a conflict

            }else{

                buildEdgeNotInTreeReason(edgeID, conflict);

            }
            outer->toSolver(conflict);
            return false;
        }
        assert(changed_edges.size() == sz);
        assert(changed_edges.last().edgeID == edgeID);
        is_edge_changed[edgeID] = false;
        changed_edges.pop();
    }
    return true;
}

template<typename Weight, typename Graph>
bool MSTDetector<Weight, Graph>::checkSatisfied(){
    Kruskal<typename MinimumSpanningTree<Weight>::NullStatus, Weight> positive_checker(g_under,
                                                                                       MinimumSpanningTree<Weight>::nullStatus,
                                                                                       0);
    Kruskal<typename MinimumSpanningTree<Weight>::NullStatus, Weight> negative_checker(g_over,
                                                                                       MinimumSpanningTree<Weight>::nullStatus,
                                                                                       0);
    positive_checker.update();
    negative_checker.update();
    for(int k = 0; k < weight_lits.size(); k++){
        Lit l = weight_lits[k].l;
        Weight& dist = weight_lits[k].min_weight;
        bool connected = positive_checker.numComponents() <= 1;
        if(l != lit_Undef){

            if(outer->value(l) == l_True){
                if(!connected || positive_checker.weight() > dist){
                    return false;
                }
            }else if(outer->value(l) == l_False){
                if(connected && negative_checker.weight() <= dist){
                    return false;
                }
            }else{
                if(connected && positive_checker.weight() <= dist){
                    return false;
                }
                if(!connected || negative_checker.weight() > dist){
                    return false;
                }
            }
        }
    }

    for(int k = 0; k < tree_edge_lits.size(); k++){
        Lit l = tree_edge_lits[k].l;
        int edgeid = tree_edge_lits[k].edgeID;

        if(l != lit_Undef){
            Var v = outer->getEdgeVar(edgeid);
            bool edgedisabled = true;
            bool edgeenabled = false;
            if(v >= 0){
                edgedisabled = outer->value(v) == l_False;
                edgeenabled = outer->value(v) == l_True;
            }
            if(outer->value(l) == l_True){
                assert(edgedisabled || positive_checker.edgeInTree(edgeid));
                if(!(edgedisabled || positive_checker.edgeInTree(edgeid))){
                    return false;
                }
            }else if(outer->value(l) == l_False){
                if(edgedisabled)
                    return false;

                if(negative_checker.edgeInTree(edgeid)){
                    return false;
                }
            }else{
                if(edgedisabled)
                    return false;
                if(positive_checker.edgeInTree(edgeid)){
                    return false;
                }
                if(!negative_checker.edgeInTree(edgeid))
                    return false;
            }
        }
    }
    Weight sum_weight = 0;
    for(int edgeid = 0; edgeid < g_under.edges(); edgeid++){
        if(positive_checker.edgeInTree(edgeid)){
            if(!negative_checker.edgeInTree(edgeid)){
                return false;
            }
            sum_weight += g_under.getWeight(edgeid);

        }else if(negative_checker.edgeInTree(edgeid)){
            return false;
        }
    }

    if(sum_weight != positive_checker.forestWeight()){
        return false;
    }
    return true;
}

template<typename Weight, typename Graph>
void MSTDetector<Weight, Graph>::printSolution(std::ostream& write_to){

    if(underapprox_detector->numComponents() > 1){
        printf("Min Spanning Tree is disconnected (%d components)\n", underapprox_detector->numComponents());
        Weight min_weight = underapprox_detector->forestWeight();
        write_to << "Min Spanning Forest Weight: " << min_weight << "\n";

    }else{
        Weight min_weight = underapprox_detector->weight();
        write_to << "Min Spanning Tree Weight: " << min_weight << "\n";
    }

    int maxw = log10(outer->nNodes()) + 1;
    int width = sqrt(outer->nNodes());
    if(opt_width > 0){
        width = opt_width;
    }
    int height = width;
    if(opt_height > 0){
        height = opt_height;
    }

    int lasty = 0;
    int extra = outer->nNodes() % width ? (width - outer->nNodes() % width) : 0;
    for(int n = 0; n < outer->nNodes(); n++){
        int x = n % width;

        int y = (n) / width;
        if(y > lasty)
            write_to << "\n";
        bool in_tree = false;
        for(int e = 0; e < g_under.edges(); e++){
            if(g_under.getEdge(e).to == n && g_under.edgeEnabled(e) && this->underapprox_detector->edgeInTree(e)){
                in_tree = true;
                break;
            }
        }
        if(in_tree){
            write_to << "-";
        }else{
            write_to << " ";
        }
    }
    write_to << "\n";

}

template<typename Weight, typename Graph>
Lit MSTDetector<Weight, Graph>::decide(CRef& decision_reason){
    /*MSTDetector *r =this;
     MinimumSpanningTree<MSTDetector<Weight,Graph>::MSTStatus> * over = (MinimumSpanningTree<MSTDetector<Weight,Graph>::MSTStatus>*) r->negative_reach_detector;

     MinimumSpanningTree<MSTDetector<Weight,Graph>::MSTStatus> * under = (MinimumSpanningTree<MSTDetector<Weight,Graph>::MSTStatus>*) r->positive_reach_detector;

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
};

template
class Monosat::MSTDetector<int>;

template
class Monosat::MSTDetector<int64_t>;

template
class Monosat::MSTDetector<double>;

template
class Monosat::MSTDetector<mpq_class>;
