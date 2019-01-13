/**************************************************************************************************
 The MIT License (MIT)

 Copyright (c) 2018, Sam Bayless

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

#ifndef MONOSAT_MONOSATINTERNAL_H
#define MONOSAT_MONOSATINTERNAL_H


#include "monosat/utils/ParseUtils.h"
#include "monosat/utils/Options.h"
#include "monosat/core/Solver.h"
#include "monosat/simp/SimpSolver.h"
#include "monosat/graph/GraphTheory.h"
#include "monosat/fsm/FSMTheory.h"
#include "monosat/pb/PbTheory.h"
#include "monosat/amo/AMOTheory.h"
#include "Monosat.h"
#include "monosat/core/Dimacs.h"
#include "monosat/bv/BVParser.h"
#include "monosat/graph/GraphParser.h"
#include "monosat/pb/PbParser.h"
#include "monosat/amo/AMOParser.h"
#include "monosat/core/Optimize.h"
#include "monosat/pb/PbSolver.h"
#include "monosat/routing/FlowRouter.h"
#include "monosat/api/Circuit.h"
#include <string>
#include <cstdio>
#include <ctime>

//Helper data structures for the Monosat API, intended for internal use only.

#ifndef __APPLE__
typedef timer_t posix_timer;
#else
//apple doesn't define timer_t
typedef void* posix_timer;
#endif

struct MonosatData {
    Monosat::Circuit<Monosat::SimpSolver> circuit;
    Monosat::BVTheorySolver<int64_t>* bv_theory = nullptr;
    Monosat::FSMTheorySolver* fsm_theory = nullptr;
    PB::PbSolver* pbsolver = nullptr;
    int time_limit = -1;
    bool has_timer = false;
    posix_timer solver_timer;
    vec<Monosat::GraphTheorySolver<int64_t>*> graphs;
    bool last_solution_optimal = true;
    bool has_conflict_clause_from_last_solution = false;
    vec<Objective> optimization_objectives;
    Dimacs<StreamBuffer, SimpSolver>* parser = nullptr;
    FILE* outfile = nullptr;
    std::string args = "";

    MonosatData(SimpSolver* solver) : circuit(*solver){

    }

    ~MonosatData(){
        for(auto* p:parser->getParsers()){
            delete (p);
        }
        delete (parser);
    }
};

#endif //MONOSAT_MONOSATINTERNAL_H
