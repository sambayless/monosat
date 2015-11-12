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

long optimize_linear(Monosat::SimpSolver * S, Monosat::BVTheorySolver<long> * bvTheory,const vec<Lit> & assume,int bvID, int time_cutoff, bool & hit_cutoff, long & n_solves);

long optimize_binary(Monosat::SimpSolver * S, Monosat::BVTheorySolver<long> * bvTheory,const vec<Lit> & assume,int bvID, int time_cutoff, bool & hit_cutoff, long & n_solves);
};
#endif /* OPTIMIZE_H_ */
