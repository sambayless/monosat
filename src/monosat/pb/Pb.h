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


#ifndef MONOSAT_PB_H
#define MONOSAT_PB_H

#include "monosat/core/SolverTypes.h"

namespace Monosat {
namespace PB {
enum Ineq {
    LT = -2,
    LEQ = -1,
    EQ = 0,
    GEQ = 1,
    GT = 2
};

//abstract interface to Psuedoboolean constraint solver
class PBConstraintSolver {
public:
    virtual ~PBConstraintSolver(){

    }

    virtual bool addConstr(const vec<Lit>& ps, const vec<int>& Cs, int rhs, Ineq ineq) = 0;

    virtual Lit
    addConditionalConstr(const vec<Lit>& ps, const vec<int>& Cs, int rhs, Ineq ineq, Lit cond = lit_Undef) = 0;

    virtual void convert(){

    }
};
}
}
#endif //MONOSAT_PB_H
