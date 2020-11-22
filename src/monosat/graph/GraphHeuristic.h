/**************************************************************************************************
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

#ifndef MONOSAT_GRAPHHEURISTIC_H
#define MONOSAT_GRAPHHEURISTIC_H

#include "monosat/graph/Detector.h"
#include "monosat/graph/GraphTheory.h"
#include "monosat/core/SolverTypes.h"
#include "monosat/core/Heuristic.h"

namespace Monosat {

template<typename Weight>
class GraphHeuristic : public Heuristic {
    Detector* detector;
    GraphTheorySolver<Weight>* outer;
    Solver* S;
    CRef local_decision_reason = CRef_Undef;

public:
    GraphHeuristic(GraphTheorySolver<Weight>* outer, Detector* detector) : outer(outer), S(outer->getSolver()),
                                                                           detector(detector){
        outer->addHeuristic(this);
        local_decision_reason = outer->newReasonMarker(this, true);//normal_detector->detectorID
    }

    int getTheoryIndex() const override{
        return outer->getTheoryIndex();
    }

    virtual Lit decide(CRef& decision_reason){
        return detector->decide(decision_reason);
    }

    Lit decideTheory(CRef& decision_reason) override{
        //first, give the main graph theory a chance to make a decision
        outer->dbg_full_sync();
        if(opt_lazy_backtrack && outer->supportsLazyBacktracking() && opt_lazy_backtrack_decisions &&
           outer->detectors.size()){//the detectors.size() check is a hack, to prevent empty graphs from forcing the decisions that they didn't originally contribute to.

            //when redeciding a literal, should check to see whether it would still be recomended as a decision by its detector...
            if(outer->lazy_trail_head != var_Undef){
                assert(outer->value(outer->lazy_trail_head) != l_Undef);
                assert(S->value(outer->toSolver(outer->lazy_trail_head)) == l_Undef);
                Lit d = mkLit(outer->lazy_trail_head, outer->value(outer->lazy_trail_head) == l_False);
                Lit solverLit = outer->toSolver(d);

                outer->stats_lazy_decisions++;
                outer->stats_decisions++;
                decision_reason = outer->graph_decision_reason;
                return solverLit;
            }
        }


        decision_reason = CRef_Undef;

        Lit l = lit_Undef;
        decision_reason = local_decision_reason;

        if(outer->satisfied_detectors[detector->getID()])
            return lit_Undef;
        l = decide(decision_reason);


        if(l != lit_Undef){

            if(opt_decide_graph_bv && !sign(l) && outer->isEdgeVar(var(l)) &&
               outer->hasBitVector(outer->getEdgeID(var(l))) && detector->supportsEdgeDecisions()){
                int edgeID = outer->getEdgeID(var(l));
                EdgeDecider<Weight>* d = dynamic_cast<EdgeDecider<Weight>*>(detector); //(EdgeDecider<Weight>*)r;
                Weight edgeWeight = -1;

                DetectorComparison op;
                if(d->decideEdgeWeight(edgeID, edgeWeight, op)){
                    BVTheorySolver<Weight>* bvTheory = outer->bvTheory;
                    assert(edgeWeight >= 0);
                    Lit bv_decision = lit_Undef;

                    Comparison bvOp;
                    if(op == DetectorComparison::leq){
                        bvOp = Comparison::leq;
                        bv_decision = bvTheory->decideBV(bvOp, outer->getEdgeBV(edgeID).getID(), edgeWeight);


                    }else if(op == DetectorComparison::lt){
                        bvOp = Comparison::lt;
                        bv_decision = bvTheory->decideBV(bvOp, outer->getEdgeBV(edgeID).getID(), edgeWeight);
                    }else if(op == DetectorComparison::geq){
                        bvOp = Comparison::geq;
                        bv_decision = bvTheory->decideBV(bvOp, outer->getEdgeBV(edgeID).getID(), edgeWeight);
                    }else if(op == DetectorComparison::gt){
                        bvOp = Comparison::gt;
                        bv_decision = bvTheory->decideBV(bvOp, outer->getEdgeBV(edgeID).getID(), edgeWeight);
                    }else if(op == DetectorComparison::eq){
                        bv_decision = bvTheory->decideBV(Comparison::leq, outer->getEdgeBV(edgeID).getID(), edgeWeight);
                        if(bv_decision == lit_Undef)
                            bv_decision = bvTheory->decideBV(Comparison::geq, outer->getEdgeBV(edgeID).getID(),
                                                             edgeWeight);
                    }else{
                        throw std::runtime_error("error in decision heuristic: not supported");
                    }
                    if(bv_decision != lit_Undef){
                        assert(S->value(bv_decision) == l_Undef);
                        outer->stats_decisions++;
                        detector->undecide(l);
                        if(S->value(bv_decision) != l_Undef){
                            throw std::runtime_error("error in decision heuristic");
                        }
                        return bv_decision;
                    }
                }
            }
            assert(l == lit_Undef || S->value(outer->toSolver(l)) == l_Undef);
            outer->stats_decisions++;
            detector->stats_decisions++;
            assert(l == lit_Undef || S->value(outer->toSolver(l)) == l_Undef);
            return outer->toSolver(l);
        }
        return l;
    }

};
}
#endif //MONOSAT_GRAPHHEURISTIC_H
