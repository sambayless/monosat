/****************************************************************************************[Solver.h]
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

#ifndef OPTIMIZE_H_
#define OPTIMIZE_H_

#include "core/Solver.h"

#include "simp/SimpSolver.h"
#include "bv/BVTheorySolver.h"
#include "core/SolverTypes.h"

#include "mtl/Vec.h"


namespace Monosat{

long optimize_linear(Monosat::SimpSolver * S, Monosat::BVTheorySolver<long> * bvTheory,const vec<Lit> & assume,int bvID, bool & hit_cutoff, long & n_solves);

long optimize_binary(Monosat::SimpSolver * S, Monosat::BVTheorySolver<long> * bvTheory,const vec<Lit> & assume,int bvID, bool & hit_cutoff, long & n_solves);

lbool optimize_and_solve(Monosat::SimpSolver & S,const vec<Lit> & assume,const vec<int> & bvs, bool & found_optimal);
};
#endif /* OPTIMIZE_H_ */
