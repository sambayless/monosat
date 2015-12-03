/****************************************************************************************[Solver.h]
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

#ifndef THEORYSOLVER_H_
#define THEORYSOLVER_H_
#include "SolverTypes.h"
namespace Monosat {

/**
 * Abstract interface to SMT theory solvers, as accessed by components of the theory solver.
 */
class TheorySolver {
public:
	virtual ~TheorySolver() {
		
	}
	virtual void addTheory(Theory * t)=0;
	virtual lbool value(Lit l)const=0;
	virtual lbool value(Var v)const=0;
	virtual bool isConstant(Var v)const=0;
	inline bool isConstant(Lit l) const{
		return isConstant(var(l));
	}
	virtual void needsPropagation(int theoryID)=0;
	virtual void instantiateLazyDecision(Lit l,int atLevel, CRef reason=CRef_Undef){

	}
	/*virtual void prependToTrail(Lit solverLit, int atLevel){

	}*/
	virtual Lit True()=0;
	//Set whether a variable can be a chosen as a decision in the SAT solver or not
	virtual void setDecisionVar(Var solverVar, bool decidable)=0;
	virtual Var newTheoryVar(Var solverVar, int theoryID, Var theoryVar)=0;
	virtual Var newVar(bool polarity = true, bool dvar = true)=0;
	virtual int nVars()const=0;
	virtual bool enqueue(Lit l, CRef reason=CRef_Undef)=0;
	virtual CRef newReasonMarker(Theory * forTheory)=0;
	virtual CRef reason(Var v)const=0;

	virtual void addClauseSafely(vec<Lit> & clause)=0;
	virtual bool addConflictClause(vec<Lit> & ps, CRef & confl_out, bool permanent)=0;
	virtual int level(Var v)const=0;
	virtual int decisionLevel()const=0;
	virtual Lit theoryDecisionLit(int theoryID){
		return mkLit(newVar(true,false));
	}
};
}
;

#endif /* THEORYSOLVER_H_ */
