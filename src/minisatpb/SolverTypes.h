/***********************************************************************************[SolverTypes.h]
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

/**************************************************************************************************

Contains the solver specific types: Var, Lit, Clause

**************************************************************************************************/


#ifndef SolverTypes_h
#define SolverTypes_h

#include "Int.h"


//=================================================================================================
// Variables, literals:


// NOTE! Variables are just integers. No abstraction here. They should be chosen from 0..N,
// so that they can be used as array indices.

typedef int Var;
#define var_Undef (-1)


class Lit {
    int     x;
public:
    Lit(void)   /* unspecifed value allowed for efficiency */      { }
    explicit Lit(Var var, bool sign = false) : x((var+var) + sign) { }
    friend Lit operator ~ (Lit p);

    friend bool sign (Lit p);
    friend int  var  (Lit p);
    friend int  index(Lit p);
    friend Lit  toLit(int i);

    friend bool operator == (Lit p, Lit q);
    friend bool operator <  (Lit p, Lit q);
};
inline Lit operator ~ (Lit p) { Lit q; q.x = p.x ^ 1; return q; }
inline bool sign (Lit p) { return p.x & 1; }
inline int  var  (Lit p) { return p.x >> 1; }
inline int  index(Lit p) { return p.x; }        // A "toInt" method that guarantees small, positive integers suitable for array indexing.
inline Lit  toLit(int i) { Lit p; p.x = i; return p; }
inline bool operator == (Lit p, Lit q) { return index(p) == index(q); }
inline bool operator <  (Lit p, Lit q) { return index(p)  < index(q); }  // '<' guarantees that p, ~p are adjacent in the ordering.

const Lit lit_Undef(var_Undef, false);  // }- Useful special constants.
const Lit lit_Error(var_Undef, true );  // }

macro uint64 abstLit(Lit p) { return ((uint64)1) << (index(p) & 63); }

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#define L_LIT    "%sx%d"
#define L_lit(p) sign(p)?"~":"", var(p)

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

struct BasicSolverStats {
    int64   starts, decisions, propagations, inspects, conflicts;
    BasicSolverStats(void) : starts(0), decisions(0), propagations(0), inspects(0), conflicts(0) { }
};

//=================================================================================================
#endif
