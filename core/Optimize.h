/*
 * Opt.h
 *
 *  Created on: Nov 12, 2015
 *      Author: sam
 */

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

lbool optimize_and_solve(Monosat::SimpSolver & S,const vec<Lit> & assume,const vec<int> & bvs);
};
#endif /* OPTIMIZE_H_ */
