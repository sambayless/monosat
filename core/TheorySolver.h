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
	virtual lbool value(Lit l)=0;
	virtual lbool value(Var v)=0;
	virtual bool isConstant(Var v)=0;
	inline bool isConstant(Lit l) {
		return isConstant(var(l));
	}
	virtual Var newVar(Var v, int detector_id)=0;
	virtual bool enqueue(Lit l, CRef reason)=0;
	virtual CRef newReasonMarker(int detectorID)=0;
	virtual void needsPropagation(int detector_id)=0;
};
}
;

#endif /* THEORYSOLVER_H_ */
