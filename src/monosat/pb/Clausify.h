/*****************************************************************************[Hardware_sorters.cc]
 Copyright (c) 2016, Sam Bayless
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
#ifndef MONOSAT_CLAUSIFY_CONTEXT_H_H
#define MONOSAT_CLAUSIFY_CONTEXT_H_H

#include "monosat/core/SolverTypes.h"
#include "monosat/pb/ADTs/FEnv.h"

namespace Monosat {
namespace PB {

struct ClausifyContext {
    CMap<int> occ;
    CMap<Var> vmap;
    CMap<Lit, true> vmapp;

    ClausifyContext() : occ(0), vmap(var_Undef), vmapp(lit_Undef){

    }
    /*CMap<int>      Clausifier::occ(0);
CMap<Var>      Clausifier::vmap(var_Undef);
CMap<Lit, true> Clausifier::vmapp(lit_Undef);*/
};
}
}

#endif //MONOSAT_CLAUSIFY_CONTEXT_H_H
