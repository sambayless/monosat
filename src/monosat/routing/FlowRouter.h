/****************************************************************************************[Solver.h]
 The MIT License (MIT)

 Copyright (c) 2017, Sam Bayless

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

#ifndef MONOSAT_FLOWROUTER_H
#define MONOSAT_FLOWROUTER_H

#include "monosat/graph/GraphTheory.h"
#include "monosat/graph/ReachDetector.h"
#include "monosat/mtl/Vec.h"
#include "monosat/mtl/IntMap.h"
#include "monosat/mtl/Heap.h"
#include "monosat/mtl/Alg.h"
#include "monosat/mtl/Rnd.h"
#include "monosat/utils/Options.h"
#include "monosat/core/SolverTypes.h"
#include "monosat/core/Theory.h"
#include "monosat/core/TheorySolver.h"
#include "monosat/core/Config.h"
using namespace dgl;
namespace Monosat {
template<typename Weight>
class FlowRouter : public Theory {
    IntSet<> edgevars;
    int theory_index=-1;
    struct Net{
        //vec<int> members;
        vec<ReachDetector<Weight> * >detectors;
        vec<Lit> reach_lits;
        //vec<int> dest_edges;
        vec<Lit> dest_edgelits;
        //int disconnectedEdge=-1;
        Lit disconnectedEdgeLit=lit_Undef;
    };
    Solver * S;
    GraphTheorySolver<Weight> *g_theory=nullptr;
    MaxflowDetector<Weight> * maxflow_detector=nullptr;
    vec<Net> nets;

    int sourceNode;
    int destNode;
    Lit flowLit;
    int routerID=-1;
    vec<Lit> enabled_routing_lits;
    vec<Lit> disabled_routing_lits;
    vec<Lit> inner_conflict;
public:
    FlowRouter(Solver * S,GraphTheorySolver<Weight> * g,int sourceNode,int destNode,Lit maxflowLit);

    int getRouterID(){
        return routerID;
    }

    void addNet(Lit disabledEdge,vec<Lit> & dest_edges, vec<Lit> & reach_lits);

    bool propagateTheory(vec<Lit> & conflict) override{
        return propagateTheory(conflict,false);
    }
    bool propagateTheory(vec<Lit> & conflict, bool solve);
    bool solveTheory(vec<Lit> & conflict) override{
        return propagateTheory(conflict,true);
    }

    inline int getTheoryIndex()const {
        return theory_index;
    }
    inline void setTheoryIndex(int id) {
        theory_index = id;
    }
    inline void newDecisionLevel() {

    }
    inline void backtrackUntil(int untilLevel){

    }
    inline int decisionLevel() {
        return S->decisionLevel();
    }
    inline void undecideTheory(Lit l){

    }
    void enqueueTheory(Lit l) {

    }
};


template<typename Weight>
FlowRouter<Weight>::FlowRouter(Solver * S,GraphTheorySolver<Weight> * g,int sourceNode,int destNode,Lit maxflowLit):S(S),g_theory(g),sourceNode(sourceNode),destNode(destNode),flowLit(maxflowLit){

    S->addTheory(this);
    routerID = getTheoryIndex();
    Theory * t = S->theories[S->getTheoryID(maxflowLit)];
   //this is really a mess...
   //and should be done more safely in the future...
   //This code is recovering the reachability detector associated with each net member
   GraphTheorySolver<Weight> * g2 = (GraphTheorySolver<Weight> *) t;
   Detector * d =  g2->detectors[g2->getDetector(var(S->getTheoryLit(maxflowLit)))];
   assert(d);
   assert(strcmp(d->getName(),"Max-flow Detector")==0);
   MaxflowDetector<Weight> * rd = (MaxflowDetector<Weight> *) d;
   this->maxflow_detector = rd;
}

template<typename Weight>
void FlowRouter<Weight>::addNet(Lit disabledEdge,vec<Lit> & dest_edges, vec<Lit> & reach_lits){
    /*
     *   struct Net{
        //vec<int> members;
        ReachDetector<Weight> * det;
        vec<Lit> member_lits;
        //vec<int> dest_edges;
        vec<Lit> dest_edgelits;
        //int disconnectedEdge=-1;
        Lit disconnectedEdgeLit=lit_Undef;
    };
     */
    nets.push();
    for(Lit l:reach_lits){
        nets.last().reach_lits.push(l);
        assert(S->hasTheory(l));
        Theory * t = S->theories[S->getTheoryID(l)];
        //this is really a mess...
        //and should be done more safely in the future...
        //This code is recovering the reachability detector associated with each net member
        GraphTheorySolver<Weight> * g = (GraphTheorySolver<Weight> *) t;
        Detector * d =  g->detectors[g->getDetector(var(S->getTheoryLit(l)))];
        assert(d);
        assert(strcmp(d->getName(),"Reachability Detector")==0);
        ReachDetector<Weight> * rd = (ReachDetector<Weight> *) d;
        nets.last().detectors.push(rd);
    }
    nets.last().disconnectedEdgeLit = disabledEdge;
    dest_edges.copyTo(nets.last().dest_edgelits);
    for(Lit l:dest_edges){
    	S->setDecisionVar(var(l),false);
    }
    S->setDecisionVar(var(disabledEdge),false);
}

template<typename Weight>
bool FlowRouter<Weight>::propagateTheory(vec<Lit> &conflict, bool solve) {
    //for each net to be routed, pick one unrouted endpoint (if any).
    //connect it to destination in g.
    if(S->decisionLevel()==0)
        return true;//don't do anything at level 0

    if(!maxflow_detector->propagate(conflict)){
        return false;
    }
    bool has_level_decision=false;
    int lev = S->decisionLevel();
    enabled_routing_lits.clear();
    disabled_routing_lits.clear();
    inner_conflict.clear();

    DynamicGraph<Weight> & g = g_theory->g_over;
    //vec<int> routing_edges;

    //nets that are completely routed are connected directly to the dest node; if all are connected then do nothing
    for(int j = 0;j<nets.size();j++){
        Net & net = nets[j];

        //check if any members of this detector are unrouted
        Lit unrouted = lit_Undef;
        int unrouted_n = -1;
        for(int i = 0;i<net.dest_edgelits.size();i++){
            ReachDetector<Weight> * r = net.detectors[i];
            Lit l = net.reach_lits[i];
            if(S->value(l)==l_True) {
                Lit edge_lit = net.dest_edgelits[i];
                //find an endpoint that is not yet connected in the under approx graph
                if (!r->isConnected(S->getTheoryLit(l),false) && unrouted==lit_Undef){
                	unrouted = edge_lit;
					unrouted_n = i;
                	if(S->value(edge_lit)==l_Undef){
						enabled_routing_lits.push(edge_lit);
                        if(!has_level_decision){
                            has_level_decision=true;
                            S->newDecisionLevel();
                        }
						S->enqueue((edge_lit),CRef_Undef);
                	}
                }else{
                	if(S->value(edge_lit)==l_Undef){
                    disabled_routing_lits.push(~edge_lit);
                    //g_theory->enqueue(~S->getTheoryLit(edge_lit),CRef_Undef);
                        if(!has_level_decision){
                            has_level_decision=true;
                            S->newDecisionLevel();
                        }
                        S->enqueue(~(edge_lit),CRef_Undef);
                	}
                }
            }
        }
        if(unrouted==lit_Undef){
            //connect source to dest directly
            //g.enableEdge(net.disconnectedEdge);
        	if(S->value(net.disconnectedEdgeLit)==l_Undef){
                if(!has_level_decision){
                    has_level_decision=true;
                    S->newDecisionLevel();
                }
        		S->enqueue((net.disconnectedEdgeLit),CRef_Undef);
        		enabled_routing_lits.push(net.disconnectedEdgeLit);
        	}
            //routing_edges.push(net.disconnectedEdgeLit);
        }else{
            //connect unrouted to dest
            //int edge_ID = net.dest_edges[unrouted_n];
            //g.enableEdge(edgeID);
            //Lit l = net.dest_edgelits[unrouted_n];
            //routing_edges.push(edge_ID);
        	if(S->value(net.disconnectedEdgeLit)==l_Undef){
                if(!has_level_decision){
                    has_level_decision=true;
                    S->newDecisionLevel();
                }
        		S->enqueue(~(net.disconnectedEdgeLit),CRef_Undef);
        		disabled_routing_lits.push(~net.disconnectedEdgeLit);
        	}
        }
    }
    inner_conflict.clear();
    conflict.clear();
    assert(lev>=1);
    if(!maxflow_detector->propagate(inner_conflict)){
        //if there is a conflict, remove any of the edge lits that were assigned here
        for(int i = 0;i<inner_conflict.size();i++){
            Lit l = inner_conflict[i];
            assert(l!=lit_Undef);
            if(edgevars.has(var(l))){
                //drop this literal
            }else{
                conflict.push(l);
            }
        }
        S->cancelUntil(lev);
        for(Lit l : enabled_routing_lits){
            //g_theory->backtrackAssign(S->getTheoryLit(l));
            //S->unsafeUnassign(l);
            assert(S->value(l)==l_Undef);
        }

        for(Lit l : disabled_routing_lits){
            //g_theory->backtrackAssign(S->getTheoryLit(l));
            //S->unsafeUnassign(l);
            assert(S->value(l)==l_Undef);
        }

        return false;
    }
    if(!solve) {
        S->cancelUntil(lev);

        for (Lit l : enabled_routing_lits) {
            //S->unsafeUnassign(l);
            //g_theory->backtrackAssign(S->getTheoryLit(l));
            assert(S->value(l) == l_Undef);
        }

        for (Lit l : disabled_routing_lits) {
            //S->unsafeUnassign(l);
            //g_theory->backtrackAssign(S->getTheoryLit(l));
            assert(S->value(l) == l_Undef);
        }
    }
    return true;

}

};

#endif //MONOSAT_FLOWROUTER_H
