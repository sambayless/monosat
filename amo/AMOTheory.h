/****************************************************************************************[Solver.h]
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

#include "mtl/Vec.h"
#include "core/SolverTypes.h"
#include "core/Theory.h"

namespace Monosat {


//At-Most-One theory. This is a special case of PB constraints, for handling at-most-one constraints.
//Each instance of this theory supports a _single_ at-most-one constraint; to implement multiple such constraints, instantiate multiple copies of the theory.
class AMOTheory: public Theory {
	Solver * S;
	int theory_index=-1;

public:

	CRef assign_false_reason;
	//CRef assign_true_reason;

	vec<Var> amo;//list of variables, at most one of which should be true.

	vec<Lit> tmp_clause;

	double propagationtime=0;
	long stats_propagations=0;
	long stats_propagations_skipped=0;
	long stats_shrink_removed = 0;
	long stats_reasons = 0;
	long stats_conflicts = 0;

	Var true_var=var_Undef;
	Var conflict_var=var_Undef;
	bool needs_propagation=false;

public:
	
	AMOTheory(Solver * S) :
			S(S) {
		S->addTheory(this);
		assign_false_reason=S->newReasonMarker(this);
		//assign_true_reason=S->newReasonMarker(this);

	}
	~AMOTheory() {
	}
	;

	//Add a variable (not literal!) to the set of which at most one may be true.
	void addVar(Var solverVar){
		S->newTheoryVar(solverVar, getTheoryIndex(),solverVar);//using same variable indices in the theory as out of the theory
		amo.push(solverVar);
	}


	inline int getTheoryIndex() {
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
		if(var(l)==true_var){
			needs_propagation=false;
			true_var=var_Undef;
			assert(conflict_var==var_Undef);
		}
		if(var(l)==conflict_var){
			conflict_var=var_Undef;
		}
	}
	void enqueueTheory(Lit l) {
		if(conflict_var==var_Undef){
			if (!sign(l)){
				if (true_var==var_Undef){
					true_var=var(l);
					assert(!needs_propagation);
					if(opt_amo_eager_prop){
						//enqueue all of the remaining lits in the solver, now.
						for(Var v:amo){
							if(v!=true_var){
								S->enqueue(mkLit(v,true),assign_false_reason);
							}
						}
					}else{
						needs_propagation=true;
					}
				}else if (var(l)==true_var){
					//we already knew this lit was assigned to true, do nothing.
				}else{
					//there is a conflict - both conflict_var and true_var are assigned true, which is not allowed.
					conflict_var=var(l);
				}
			}else{
				//it is always safe to assign a var to false.
			}
		}
	}
	;
	bool propagateTheory(vec<Lit> & conflict) {
		if(conflict_var!=var_Undef){
			conflict.clear();
			assert(true_var!=var_Undef);
			assert(true_var!=conflict_var);
			conflict.push(mkLit(conflict_var,true));
			conflict.push(mkLit(true_var,true));
			needs_propagation=false;
			return false;
		}else if (true_var!=var_Undef && needs_propagation){
			needs_propagation=false;
			assert(!opt_amo_eager_prop);
			//enqueue all of the remaining lits in the solver, now.
			for(Var v:amo){
				if(v!=true_var){
					S->enqueue(mkLit(v,true),assign_false_reason);
				}
			}
		}
		return true;
	}
	inline bool solveTheory(vec<Lit> & conflict){
		return propagateTheory(conflict);
	}
	inline void buildReason(Lit p, vec<Lit> & reason, CRef reason_marker){
		assert(reason_marker==assign_false_reason);
		if(var(p)!=true_var){
			assert(sign(p));
			assert(S->value(p)==l_True);
			assert(S->value(true_var)==l_True);
			reason.push(mkLit( var(p),true));
			reason.push(mkLit(true_var,true));//either true_var (currently assigned true) must be false, or var(p) must be false
		}else{
			assert(false);
		}
	}
	bool check_solved() {
		int n_true=0;
		for (Var v:amo){
			if(S->value(v)==l_True){
				n_true+=1;
				if(n_true>1){
					return false;
				}
			}
		}

		return true;
	}
private:

	
};

}
;

#endif /* AMOTheory_H_ */
