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

#ifndef FSM_DETECTOR_H_
#define FSM_DETECTOR_H_

#include "monosat/core/SolverTypes.h"

#include "monosat/mtl/Vec.h"
#include "monosat/core/Config.h"
#include <cstdio>
#include <iostream>

namespace Monosat {

class FSMDetector {
public:

    int detectorID;

    int unassigned_positives = 0;
    int unassigned_negatives = 0;
    int64_t last_negated_check = 0;

    //Stats
    double stats_under_update_time = 0;
    double stats_over_update_time = 0;
    int64_t stats_under_updates = 0;
    int64_t stats_over_updates = 0;
    int64_t stats_under_conflicts = 0;
    int64_t stats_over_conflicts = 0;
    double stats_under_conflict_time = 0;
    double stats_over_conflict_time = 0;
    int64_t stats_skipped_under_updates = 0;
    int64_t stats_skipped_over_updates = 0;
    int64_t stats_decisions = 0;
    double stats_decide_time = 0;
    int64_t stats_under_clause_length = 0;
    int64_t stats_over_clause_length = 0;

    int getID(){
        return detectorID;
    }

    virtual const char* getName(){
        return "<unknown>";
    }

    virtual void printStats(){
        if(opt_verb > 0){
            printf("Detector %d (%s):\n", getID(), getName());
            //printf("Updates: %d (under), %d over\n", stats_under_updates, stats_over_updates);
            printf("\tUnder-approx updates: %" PRId64 " (%" PRId64 " skipped) (%f s total, %f s avg)\n",
                   stats_under_updates,
                   stats_skipped_under_updates, (double) stats_under_update_time,
                   (double) stats_under_update_time / (double) (stats_under_updates + 1));
            printf("\tOver-approx updates: %" PRId64 " (%" PRId64 " skipped)  (%f s total, %f s avg)\n",
                   stats_over_updates,
                   stats_skipped_over_updates, (double) stats_over_update_time,
                   (double) stats_over_update_time / (double) (stats_over_updates + 1));
            printf("\tTheory Decisions: %" PRId64 " (%f s total, %f s avg)\n", stats_decisions,
                   (double) stats_decide_time,
                   (double) stats_decide_time / (double) (stats_decisions + 1));
            printf(
                    "\tConflicts (under,over): %" PRId64 " (clause literals: %" PRId64 "), %" PRId64 ", (clause literals: %" PRId64 "), (under time %f s, over time %f s)\n",
                    stats_under_conflicts, stats_under_clause_length, stats_over_conflicts, stats_over_clause_length,
                    stats_under_conflict_time, stats_over_conflict_time);

        }
    }

    virtual void printSolution(std::ostream& write_to = std::cout){
    }

    virtual bool propagate(vec<Lit>& conflict) = 0;

    virtual void buildReason(Lit p, vec<Lit>& reason, CRef marker) = 0;

    virtual bool checkSatisfied() = 0;

    virtual void preprocess(){

    }

    virtual void backtrack(int level){
        //do nothing
    }

    virtual Lit decide(int level){
        return lit_Undef;
    }

    virtual void setOccurs(Lit l, bool occurs){

        if(!occurs){
            if(sign(l))
                unassigned_negatives--;
            else
                unassigned_positives--;
        }else{
            if(sign(l))
                unassigned_negatives++;
            else
                unassigned_positives++;
        }
        assert(unassigned_positives >= 0);
        assert(unassigned_negatives >= 0);
    }

    virtual void assign(Lit l){
        if(!sign(l))
            unassigned_negatives++;
        else
            unassigned_positives++;
        assert(unassigned_positives >= 0);
        assert(unassigned_negatives >= 0);

    }

    virtual void unassign(Lit l){
        if(!sign(l))
            unassigned_negatives--;
        else
            unassigned_positives--;
        assert(unassigned_positives >= 0);
        assert(unassigned_negatives >= 0);
    }

    virtual void setSatisfied(Lit l, bool isSatisfied){

    }

    virtual bool checkNegatedPolarity(){
        if(opt_detect_satisfied_predicates <= 0)
            return false;

        if(++last_negated_check >= opt_detect_satisfied_predicates){
            last_negated_check = 0;
            return true;
        }else{
            return false;
        }
    }

    FSMDetector(int detectorID) :
            detectorID(detectorID), unassigned_positives(0), unassigned_negatives(0){
    };

    virtual ~FSMDetector(){

    }


    virtual void addVar(Var v){
    }
};

class FSMLevelDetector : public FSMDetector {

    vec<int> level_trail;

public:
    FSMLevelDetector(int detectorID) :
            FSMDetector(detectorID){
    };

    virtual ~FSMLevelDetector(){
    }

    virtual void backtrack(int level){
        while(level_trail.size() && level < level_trail.last()){
            level_trail.pop();
            localBacktrack();
        }
    }

    virtual void localBacktrack(){
        //do nothing
    }

    //local decision level within the detector
    int decisionLevel(){
        return level_trail.size();
    }

    void newDecisionLevel(int outer_level){

        level_trail.push(outer_level);
    }

};

};

#endif /* FSM_DETECTOR_H_ */
