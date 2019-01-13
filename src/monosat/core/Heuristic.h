/**************************************************************************************************
 The MIT License (MIT)

 Copyright (c) 2016, Sam Bayless

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

#ifndef MONOSAT_HEURISTIC_H
#define MONOSAT_HEURISTIC_H

#include "monosat/core/SolverTypes.h"

namespace Monosat {
/**
 * Abstract interface for decision heuristics for theory solvers
 */
class Heuristic {
    int priority = 0;
    double activity = 0;
    int heuristic_index = -1;
    int heuristic_order = 0;
    Heuristic* parent = nullptr;
public:


    virtual ~Heuristic(){

    }

    int getPriority() const{
        return priority;
    }

    void setPriority(int p){
        priority = p;
    }

    int getHeuristicOrder() const{
        return heuristic_order;
    }

    void setHeuristicOrder(int p){
        heuristic_order = p;
    }

    double& getActivity(){
        return activity;
    }

    void setActivity(double p){
        activity = p;
    }

    //the theory index of the theory this heuristic belongs to, or -1 if there is no theory.
    virtual int getTheoryIndex() const{
        return -1;
    }

    int getHeuristicIndex() const{
        return heuristic_index;
    }

    void setHeuristicIndex(int id){
        assert(id > 0);
        heuristic_index = id;
    }

    virtual Lit decideTheory(CRef& decision_reason){
        decision_reason = CRef_Undef;
        return lit_Undef;
    }

    void setParentHeuristic(Heuristic* parent){
        this->parent = parent;
    }

    Heuristic* getParentHeuristic() const{
        return this->parent;
    }
};
}

#endif //MONOSAT_HEURISTIC_H
