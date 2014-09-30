/*
 * Aiger.cpp
 *
 *  Created on: 2012-10-18
 *      Author: sam
 */

//Utility functions for using AIG with Minisat

extern "C" {
#include "aiger/aiger.h"
}
#include "core/Solver.h"

Monosat::Lit unroll(Monosat::Solver& target,aiger * mgr, int iterations, Monosat::vec<Monosat::Var> & in_latches, Monosat::vec<Monosat::Var> & out_latches);

Monosat::Lit unroll(Monosat::Solver& target,aiger * mgr, Monosat::vec<Monosat::Var> & in_latches, Monosat::vec<Monosat::Var> & out_latches);
bool isReset(Monosat::vec<Monosat::Lit> & assignment);
void prepare(Monosat::Solver& target,aiger * mgr, Monosat::vec<Monosat::Var> & latches);
int zero(Monosat::Solver& target,aiger * mgr, Monosat::vec<Monosat::Var> & latches );
void add_passthrough(aiger*mgr, bool add_passthrough_latch);
