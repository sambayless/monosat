/**************************************************************************************************
 The MIT License (MIT)

 Copyright (c) 2015, Sam Bayless

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

#ifndef AMOTHEORY_H_
#define AMOTHEORY_H_

#include "monosat/mtl/Vec.h"
#include "monosat/core/SolverTypes.h"
#include "monosat/core/Theory.h"

namespace Monosat {


//At-Most-One theory. This is a special case of PB constraints, for handling at-most-one constraints.
//Each instance of this theory supports a _single_ at-most-one constraint; to implement multiple such constraints, instantiate multiple copies of the theory.
class AMOTheory : public Theory {
    Solver* S;
    int theory_index = -1;

public:

    CRef assign_false_reason;
    //CRef assign_true_reason;

    vec<Lit> amo;//list of variables, at most one of which should be true.
    VMap<Lit> positive;
    vec<Lit> tmp_clause;

    double propagationtime = 0;
    int64_t stats_propagations = 0;
    int64_t stats_lit_propagations = 0;
    int64_t stats_propagations_skipped = 0;
    int64_t stats_shrink_removed = 0;
    int64_t stats_reasons = 0;
    int64_t stats_conflicts = 0;

    Lit true_lit = lit_Undef;
    Lit conflict_lit = lit_Undef;
    bool needs_propagation = false;
    bool clausified = false;
public:
    const char* getTheoryType() override{
        return "AMO";
    }

    AMOTheory(Solver* S) :
            S(S){
        S->addTheory(this);
        assign_false_reason = S->newReasonMarker(this);
        //assign_true_reason=S->newReasonMarker(this);

    }

    ~AMOTheory(){
    };

    static bool clausify_amo(Solver* S, const vec<Lit>& lits){
        vec<Lit> set;

        assert(S->decisionLevel() == 0);
        Lit constant_true_lit = lit_Undef;

        for(int i = 0; i < lits.size(); i++){
            Lit l = lits[i];
            if(S->value(l) == l_False && S->level(var(l)) == 0){
                //don't add this to the set
            }else if(S->value(l) == l_True && S->level(var(l)) == 0){
                if(constant_true_lit == lit_Undef){
                    constant_true_lit = l;
                }else{
                    //this is a conflict
                    S->addClause(~constant_true_lit, ~l);//this is a conflict
                    return false;
                }
            }else{
                set.push(l);
            }
        }

        if(constant_true_lit == lit_Undef){
            for(int i = 0; i < set.size(); i++){
                for(int j = i + 1; j < set.size(); j++){
                    S->addClause(~set[i], ~set[j]);
                }
            }
        }else{
            //all remaining elements of the set must be false, because constant_true_lit is true.
            for(int i = 0; i < set.size(); i++){
                S->addClause(~constant_true_lit,
                             ~set[i]);//technically don't need to include ~constant_true_lit here, but it might make things cleaner to reason about elsewhere.
                // (It will be eliminated by the solver anyhow, so this doesn't cost anything.)
            }
        }
        return true;
    }


    //Add a variable (not literal!) to the set of which at most one may be true.
    void addVar(Var solverVar){
        addLit(mkLit(solverVar));

    }

    //Add a literal to the set of which at most one may be true.
    void addLit(Lit solverLit){
        S->newTheoryVar(var(solverLit), getTheoryIndex(),
                        var(solverLit));//using same variable indices in the theory as out of the theory
        amo.push(solverLit);
        positive.insert(var(solverLit), solverLit, lit_Undef);
    }


    inline int getTheoryIndex() const override{
        return theory_index;
    }

    inline void setTheoryIndex(int id) override{
        theory_index = id;
    }

    inline void newDecisionLevel() override{

    }

    inline void backtrackUntil(int untilLevel) override{

    }

    inline int decisionLevel(){
        return S->decisionLevel();
    }

    inline void undecideTheory(Lit l) override{
        if(l == true_lit){
            needs_propagation = false;
            true_lit = lit_Undef;
            assert(conflict_lit == lit_Undef);
        }
        if(l == conflict_lit){
            conflict_lit = lit_Undef;
        }
    }

    void enqueueTheory(Lit l) override{
        if(clausified)
            return;
        if(conflict_lit == lit_Undef){
            Lit expected = positive[var(l)];
            if(l == expected){
                if(true_lit == lit_Undef){
                    true_lit = l;
                    assert(!needs_propagation);
                    if(opt_amo_eager_prop){
                        //enqueue all of the remaining lits in the solver, now.
                        stats_propagations++;
                        for(Lit l:amo){
                            if(l != true_lit){
                                S->enqueue(~l, assign_false_reason);
                            }
                        }
                    }else{
                        needs_propagation = true;
                    }
                }else if(l == true_lit){
                    //we already knew this lit was assigned to true, do nothing.
                }else{
                    //there is a conflict - both conflict_lit and true_lit are assigned true, which is not allowed.
                    conflict_lit = l;
                }
            }else{
                //it is always safe to assign a var to false.
            }
        }
    };

    bool propagateTheory(vec<Lit>& conflict) override{
        if(clausified){
            S->setTheorySatisfied(this);
            return true;
        }
        S->theoryPropagated(this);
        if(decisionLevel() == 0){
            //remove constants from the set
            int n_nonconstants = 0;
            bool has_true_lit = false;
            int i, j = 0;
            for(i = 0; i < amo.size(); i++){
                Lit l = amo[i];
                if(S->value(l) == l_False && S->level(var(l)) == 0){
                    //drop this literal from the set
                }else if(S->value(l) == l_True && S->level(var(l)) == 0){
                    amo[j++] = l;
                    has_true_lit = true;
                }else{
                    n_nonconstants++;
                    amo[j++] = l;
                }
            }
            amo.shrink(i - j);
            if(has_true_lit || amo.size() == 0 || amo.size() <= opt_clausify_amo){
                clausified = true;
                vec<Lit> amoLits;
                for(Lit l:amo){
                    amoLits.push(l);
                }
                if(opt_verb > 1){
                    printf("Clausifying amo theory %d with %d lits\n", this->getTheoryIndex(), amo.size());
                }
                S->setTheorySatisfied(this);
                return clausify_amo(S, amoLits);
            }
        }

        if(conflict_lit != lit_Undef){
            conflict.clear();
            assert(true_lit != lit_Undef);
            assert(true_lit != conflict_lit);
            conflict.push(~conflict_lit);
            conflict.push(~true_lit);
            needs_propagation = false;
            stats_conflicts++;
            return false;
        }else if(true_lit != lit_Undef && needs_propagation){
            stats_propagations++;
            needs_propagation = false;
            assert(!opt_amo_eager_prop);
            //enqueue all of the remaining lits in the solver, now.
            for(Lit l:amo){
                if(l != true_lit){
                    stats_lit_propagations++;
                    S->enqueue(~l, assign_false_reason);
                }
            }
        }
        return true;
    }

    void printStats(int detailLevel) override{
        if(!clausified){
            printf("AMO Theory %d stats:\n", this->getTheoryIndex());

            printf("Propagations: %" PRId64 " (%f s, avg: %f s, %" PRId64 " skipped,  %" PRId64 " lits)\n",
                   stats_propagations, propagationtime,
                   (propagationtime) / ((double) stats_propagations + 1), stats_propagations_skipped,
                   stats_lit_propagations);

            printf("Conflicts: %" PRId64 "\n", stats_conflicts);
            printf("Reasons: %" PRId64 "\n", stats_reasons);

            fflush(stdout);
        }
    }

    inline bool solveTheory(vec<Lit>& conflict) override{
        return propagateTheory(conflict);
    }

    inline void buildReason(Lit p, vec<Lit>& reason, CRef reason_marker) override{
        stats_reasons++;
        assert(reason_marker == assign_false_reason);
        if(p != true_lit){
            assert(sign(p));
            assert(S->value(p) == l_True);
            assert(S->value(true_lit) == l_True);
            reason.push(p);
            reason.push(~true_lit);//either true_lit (currently assigned true) must be false, or p must be false
        }else{
            assert(false);
        }
    }

    bool check_solved() override{
        int n_true = 0;
        for(Lit l:amo){
            if(S->value(l) == l_True){
                n_true += 1;
                if(n_true > 1){
                    return false;
                }
            }
        }

        return true;
    }

private:


};

};

#endif /* AMOTheory_H_ */
