/**************************************************************************[PbSolver_convertAdd.cc]
Copyright (c) 2005-2010, Niklas Een, Niklas Sorensson

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

#include "PbSolver.h"
#include "Hardware.h"

namespace Monosat {
namespace PB {

// Write 'd' in binary, then substitute 0 with '_0_', 1 with 'f'. This is the resulting 'out' vector.
//
static inline void bitAdder(Int d, Formula f, vec<Formula>& out){
    out.clear();
    for(; d != 0; d >>= 1)
        out.push(((d & 1) != 0) ? f : _0_);
}

// Produce a conjunction of formulas that forces the constraint 'xs <= ys'. The formulas will
// be pushed onto 'out'.
//
static void lte(vec<Formula>& xs, vec<Formula>& ys, vec<Formula>& out){
    for(int i = 0; i < xs.size(); i++){
        Formula c = _0_;
        for(int j = i + 1; j < max(xs.size(), ys.size()); j++){
            Formula x = j < xs.size() ? xs[j] : _0_;
            Formula y = j < ys.size() ? ys[j] : _0_;
            c = c | (x ^ y);
        }
        c = c | ~xs[i] | (i < ys.size() ? ys[i] : _0_);
        assert(c != _0_);
        if(c != _1_)
            out.push(c);
    }

}

void linearAddition(const Linear& l, vec<Formula>& out){
    vec<Formula> sum;
    vec<Formula> inp;
    vec<Int> cs;

    for(int i = 0; i < l.size; i++){
        inp.push(id(var(var(l[i])), sign(l[i])));
        cs.push(l(i));
    }

    Int maxlim = (l.hi != Int_MAX) ? l.hi : (l.lo - 1);
    int bits = 0;
    for(Int i = maxlim; i != 0; i >>= 1)
        bits++;

    int nodes = FEnv::nodes.size();

    addPb(inp, cs, sum, bits);
    if(opt_verbosity >= 1){
        char* tmp = toString(maxlim);
        reportf("Adder-cost: %d   maxlim: %s   bits: %d/%d\n", FEnv::nodes.size() - nodes, tmp, sum.size(), bits);
        xfree(tmp);
    }


    if(l.lo != Int_MIN){
        //reportf("lower limit\n");
        bitAdder(l.lo, _1_, inp);
        lte(inp, sum, out);
    }
    if(l.hi != Int_MAX){
        //reportf("upper limit\n");
        bitAdder(l.hi, _1_, inp);
        lte(sum, inp, out);
    }
}


//-------------------------------------------------------------------------------------------------
// (old)
}
}
#if 0
namespace Monosat{
namespace PB{
bool PbSolver::convertPbs(void)
{
    vec<Formula> fs;
    vec<Lit>     ts;

    if (!rewriteAlmostClauses())
        return false;

    // construct formulas for all constraints
    for (int i = 0; i < constrs.size(); i++){
        if (constrs[i] == NULL) continue;

        /**/reportf("---[%4d]---> ", constrs.size() - 1 - i);
        linearAddition2(*constrs[i], fs);
    }

    clausify(sat_solver,fs,ts);

    for (int i = 0; i < ts.size(); i++)
        sat_solver.addUnit(ts[i]);

    constrs.clear();
    mem.clear();

    return sat_solver.okay();
}


static void linearAddition(const Linear& l, vec<Formula>& out)
{
    vec<Formula> sum;
    vec<Formula> x, y;

    //reportf("linearAddition\n");
    for (int i = 0; i < l.size; i++){
        sum.moveTo(x);
        bitAdder(l(i), id(var(var(l[i])),sign(l[i])), y);
        //reportf("step %d, x = %d, y = %d\n", i, x.size(), y.size());
        rippleAdder(x,y,sum);
    }

    //reportf("size: %d\n", sum.size());

    if (l.lo != Int_MIN){
        //reportf("lower limit\n");
        bitAdder(l.lo,_1_,x);
        lte(x,sum,out);
    }
    if (l.hi != Int_MAX){
        //reportf("upper limit\n");
        bitAdder(l.lo,_1_,x);
        lte(sum,x,out);
    }
}
}
}
#endif
