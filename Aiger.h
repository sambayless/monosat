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

Minisat::Lit unroll(Minisat::Solver& target,aiger * mgr, int iterations, Minisat::vec<Minisat::Var> & in_latches, Minisat::vec<Minisat::Var> & out_latches);

Minisat::Lit unroll(Minisat::Solver& target,aiger * mgr, Minisat::vec<Minisat::Var> & in_latches, Minisat::vec<Minisat::Var> & out_latches);
bool isReset(Minisat::vec<Minisat::Lit> & assignment);
void prepare(Minisat::Solver& target,aiger * mgr, Minisat::vec<Minisat::Var> & latches);
int zero(Minisat::Solver& target,aiger * mgr, Minisat::vec<Minisat::Var> & latches );
void add_passthrough(aiger*mgr, bool add_passthrough_latch);
