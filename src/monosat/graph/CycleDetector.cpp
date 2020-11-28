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
#include "monosat/dgl/PKTopologicalSort.h"

using namespace Monosat;

template<typename Weight, typename Graph>
CycleDetector<Weight, Graph>::CycleDetector(int _detectorID, GraphTheorySolver<Weight>* _outer, Graph& g_under,
                                            Graph& g_over, bool detect_directed_cycles, double seed) :
        Detector(_detectorID), outer(_outer), g_under(g_under), g_over(g_over), rnd_seed(seed),
        underapprox_directed_cycle_detector(NULL), overapprox_directed_cycle_detector(
        NULL){

    undirected_acyclic_lit = lit_Undef;
    directed_acyclic_lit = lit_Undef;


    if(cyclealg == CycleAlg::ALG_DFS_CYCLE){
        underapprox_directed_cycle_detector = new DFSCycle<Weight, Graph, true, true>(g_under, detect_directed_cycles,
                                                                                      1);
        overapprox_directed_cycle_detector = new DFSCycle<Weight, Graph, true, true>(g_over, detect_directed_cycles, 1);

        overapprox_undirected_cycle_detector = overapprox_directed_cycle_detector;
        underapprox_undirected_cycle_detector = underapprox_directed_cycle_detector;

    }else if(cyclealg == CycleAlg::ALG_PK_CYCLE){
        underapprox_directed_cycle_detector = new PKToplogicalSort<Weight, Graph>(g_under, 1);
        overapprox_directed_cycle_detector = new PKToplogicalSort<Weight, Graph>(g_over, 1);

        underapprox_undirected_cycle_detector = new DFSCycle<Weight, Graph, false, true>(g_under, false, 1);
        overapprox_undirected_cycle_detector = new DFSCycle<Weight, Graph, false, true>(g_over, false, 1);
    }
    directed_cycle_marker = outer->newReasonMarker(getID());
    no_directed_cycle_marker = outer->newReasonMarker(getID());

    undirected_cycle_marker = outer->newReasonMarker(getID());
    no_undirected_cycle_marker = outer->newReasonMarker(getID());
}

template<typename Weight, typename Graph>
Lit CycleDetector<Weight, Graph>::addAcyclicLit(bool directed, Var outer_reach_var){
    g_under.invalidate();
    g_over.invalidate();
    if(!directed){
        if(undirected_acyclic_lit == lit_Undef){

            Var v = outer->newVar(outer_reach_var, getID());
            Lit l = mkLit(v, false);
            undirected_acyclic_lit = l;
            return l;
        }else{
            if(outer_reach_var == var_Undef){
                return undirected_acyclic_lit;
            }else{
                Var v = outer->newVar(outer_reach_var, getID());
                Lit l = mkLit(v, false);
                outer->makeEqual(l, undirected_acyclic_lit, true);
                return l;
            }
        }
    }else{
        if(directed_acyclic_lit == lit_Undef){
            Var v = outer->newVar(outer_reach_var, getID());
            Lit l = mkLit(v, false);
            directed_acyclic_lit = l;
            return l;
        }else{
            if(outer_reach_var == var_Undef){
                return directed_acyclic_lit;
            }else{
                Var v = outer->newVar(outer_reach_var, getID());
                Lit l = mkLit(v, false);
                outer->makeEqual(l, directed_acyclic_lit, true);
                return l;
            }
        }
    }

}

template<typename Weight, typename Graph>
void CycleDetector<Weight, Graph>::buildNoUndirectedCycleReason(vec<Lit>& conflict){
    //its clear that we can do better than this, but its also not clear how to do so efficiently...
    //for now, learn the trivial clause...

    assert(!overapprox_undirected_cycle_detector->hasUndirectedCycle());
    for(int edgeID = 0; edgeID < g_over.edges(); edgeID++){
        if(g_over.hasEdge(edgeID) && !g_over.edgeEnabled(edgeID)){
            Var v = outer->getEdgeVar(edgeID);
            if(outer->value(v) == l_False){
                conflict.push(mkLit(v, false));
            }
        }
    }
}

template<typename Weight, typename Graph>
void CycleDetector<Weight, Graph>::buildNoDirectedCycleReason(vec<Lit>& conflict){
    //its clear that we can do better than this, but its also not clear how to do so efficiently...
    //for now, learn the trivial clause...
    //One thing you could do would be to first do an over-approx cycle detection at level 0, and exclude from here any edges that can't possibly be part of any scc.
    //Is that the best one can do?
    assert(!overapprox_directed_cycle_detector->hasDirectedCycle());
    for(int edgeID = 0; edgeID < g_over.edges(); edgeID++){
        if(g_over.hasEdge(edgeID) && !g_over.edgeEnabled(edgeID)){
            Var v = outer->getEdgeVar(edgeID);
            if(outer->value(v) == l_False){
                conflict.push(mkLit(v, false));
            }
        }
    }

}

template<typename Weight, typename Graph>
void CycleDetector<Weight, Graph>::buildUndirectedCycleReason(vec<Lit>& conflict){
    assert(underapprox_undirected_cycle_detector->hasUndirectedCycle());

    std::vector<int>& cycle = underapprox_undirected_cycle_detector->getUndirectedCycle();
    for(int i = 0; i < cycle.size(); i++){
        int edgeID = cycle[i];
        Lit l = mkLit(outer->getEdgeVar(edgeID), false);
        assert(outer->value(l) == l_True);
        conflict.push(~l);
    }

}

template<typename Weight, typename Graph>
void CycleDetector<Weight, Graph>::buildDirectedCycleReason(vec<Lit>& conflict){
    assert(underapprox_directed_cycle_detector->hasDirectedCycle());

    std::vector<int>& cycle = underapprox_directed_cycle_detector->getDirectedCycle();

    for(int i = 0; i < cycle.size(); i++){
        int edgeID = cycle[i];
        Lit l = mkLit(outer->getEdgeVar(edgeID), false);
        assert(outer->value(l) == l_True);
        conflict.push(~l);
    }

}

template<typename Weight, typename Graph>
void CycleDetector<Weight, Graph>::buildReason(Lit p, vec<Lit>& reason, CRef marker){

    if(marker == directed_cycle_marker){
        reason.push(p);

        buildDirectedCycleReason(reason);

    }else if(marker == no_directed_cycle_marker){
        reason.push(p);

        buildNoDirectedCycleReason(reason);

    }else if(marker == undirected_cycle_marker){
        reason.push(p);

        buildUndirectedCycleReason(reason);

    }else if(marker == no_undirected_cycle_marker){
        reason.push(p);

        buildNoUndirectedCycleReason(reason);

    }else{
        assert(false);
    }
    outer->toSolver(reason);
}

template<typename Weight, typename Graph>
bool CycleDetector<Weight, Graph>::propagate(vec<Lit>& conflict){
    if(directed_acyclic_lit != lit_Undef){

        if(outer->value(directed_acyclic_lit) != l_False && underapprox_directed_cycle_detector->hasDirectedCycle()){

            Lit l = ~directed_acyclic_lit;

            if(outer->value(l) == l_True){
                //do nothing
            }else if(outer->value(l) == l_Undef){
                //trail.push(Assignment(false,false,detectorID,0,var(l)));
                outer->enqueue(l, no_directed_cycle_marker);
            }else if(outer->value(l) == l_False){
                conflict.push(l);
                buildDirectedCycleReason(conflict);
                outer->toSolver(conflict);
                return false;
            }
        }else if(outer->value(directed_acyclic_lit) != l_True &&
                 !overapprox_directed_cycle_detector->hasDirectedCycle()){
            Lit l = directed_acyclic_lit;

            if(outer->value(l) == l_True){
                //do nothing
            }else if(outer->value(l) == l_Undef){
                outer->enqueue(l, directed_cycle_marker);
            }else if(outer->value(l) == l_False){
                conflict.push(l);
                buildNoDirectedCycleReason(conflict);
                outer->toSolver(conflict);
                return false;
            }
        }

    }

    if(undirected_acyclic_lit != lit_Undef){

        if(outer->value(undirected_acyclic_lit) != l_False &&
           underapprox_undirected_cycle_detector->hasUndirectedCycle()){

            Lit l = ~undirected_acyclic_lit;

            if(outer->value(l) == l_True){
                //do nothing
            }else if(outer->value(l) == l_Undef){
                outer->enqueue(l, no_undirected_cycle_marker);
            }else if(outer->value(l) == l_False){
                conflict.push(l);
                buildUndirectedCycleReason(conflict);
                outer->toSolver(conflict);
                return false;
            }
        }else if(outer->value(undirected_acyclic_lit) != l_True &&
                 !overapprox_undirected_cycle_detector->hasUndirectedCycle()){
            Lit l = undirected_acyclic_lit;

            if(outer->value(l) == l_True){
                //do nothing
            }else if(outer->value(l) == l_Undef){
                outer->enqueue(l, undirected_cycle_marker);
            }else if(outer->value(l) == l_False){
                conflict.push(l);
                buildNoUndirectedCycleReason(conflict);
                outer->toSolver(conflict);
                return false;
            }
        }
    }
    return true;
}

template<typename Weight, typename Graph>
bool CycleDetector<Weight, Graph>::checkSatisfied(){
    DFSCycle<Weight, Graph, true, false> checkDirected(g_under, true, 1);
    DFSCycle<Weight, Graph, false, true> checkUndirected(g_under, false, 1);
    if(directed_acyclic_lit != lit_Undef){
        if(outer->value(directed_acyclic_lit) == l_True && checkDirected.hasDirectedCycle()){
            return false;
        }else if(outer->value(directed_acyclic_lit) == l_False && !checkDirected.hasDirectedCycle()){
            return false;
        }

    }
    if(undirected_acyclic_lit != lit_Undef){
        if(outer->value(undirected_acyclic_lit) == l_True && checkUndirected.hasUndirectedCycle()){
            return false;
        }else if(outer->value(undirected_acyclic_lit) == l_False && !checkUndirected.hasUndirectedCycle()){
            return false;
        }

    }
    return true;
}

template<typename Weight, typename Graph>
Lit CycleDetector<Weight, Graph>::decide(CRef& decision_reason){

    return lit_Undef;
};

template
class Monosat::CycleDetector<int>;

template
class Monosat::CycleDetector<int64_t>;

template
class Monosat::CycleDetector<double>;

template
class Monosat::CycleDetector<mpq_class>;
