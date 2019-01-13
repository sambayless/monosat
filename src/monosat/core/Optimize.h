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

#ifndef OPTIMIZE_H_
#define OPTIMIZE_H_

#include "monosat/core/Solver.h"

#include "monosat/simp/SimpSolver.h"
#include "monosat/bv/BVTheorySolver.h"
#include "monosat/core/SolverTypes.h"

#include "monosat/mtl/Vec.h"


namespace Monosat {

//An individual goal during optimization. If there are multiple goals, they will be optimized in lexicographic order.
struct Objective {
    enum class Type {
        BV, PB
    };
    bool maximize = false;//else, minimize
    Type type;
    int bvID = -1;//used if type is BV
    //Used if type is PB
    vec<Lit> pb_lits;
    vec<int> pb_weights; //if empty, weights are treated as '1'
    bool isPB() const{
        return type == Type::PB;
    }

    bool isBV() const{
        return type == Type::BV;
    }

    Objective(const Objective& from) : maximize(from.maximize), type(from.type), bvID(from.bvID){
        from.pb_lits.copyTo(pb_lits);
        from.pb_weights.copyTo(pb_weights);
    }

    Objective(Objective&& from) : maximize(from.maximize), type(from.type), bvID(from.bvID){
        from.pb_lits.copyTo(pb_lits);
        from.pb_weights.copyTo(pb_weights);
        from.pb_lits.clear();
        from.pb_weights.clear();
        from.bvID = -1;
    }

    Objective& operator=(Objective&& from){
        if(this != &from){
            maximize = from.maximize;
            type = from.type;
            bvID = from.bvID;
            from.pb_lits.copyTo(pb_lits);
            from.pb_weights.copyTo(pb_weights);
            from.pb_lits.clear();
            from.pb_weights.clear();
            from.bvID = -1;
        }
        return *this;
    }

    Objective() : maximize(false), type(Type::BV), bvID(-1){
    }

    Objective(int bvID, bool maximize) : maximize(maximize), type(Type::BV), bvID(bvID){
    }

    Objective(const vec<Lit>& lits, bool maximize) : maximize(maximize), type(Type::PB), bvID(-1){
        lits.copyTo(pb_lits);
        pb_weights.growTo(pb_lits.size(), 1);
    }

    Objective(const vec<Lit>& lits, const vec<int>& weights, bool maximize) : maximize(maximize), type(Type::PB),
                                                                              bvID(-1){
        lits.copyTo(pb_lits);
        weights.copyTo(pb_weights);
        if(pb_weights.size() > pb_lits.size())
            pb_weights.shrink(pb_weights.size() - pb_lits.size());
        pb_weights.growTo(pb_lits.size(), 1);
    }

};

int64_t
optimize_linear(Monosat::SimpSolver* S, Monosat::BVTheorySolver<int64_t>* bvTheory, const vec<Lit>& assume, int bvID,
                bool& hit_cutoff, int64_t& n_solves);

int64_t
optimize_binary(Monosat::SimpSolver* S, Monosat::BVTheorySolver<int64_t>* bvTheory, const vec<Lit>& assume, int bvID,
                bool& hit_cutoff, int64_t& n_solves);

lbool optimize_and_solve(Monosat::SimpSolver& S, const vec<Lit>& assume, const vec<Objective>& objectives, bool do_simp,
                         bool& found_optimal);


//Reduce the given assumptions to a (locally) minimal unsat core, if they are mutually unsat.
//Returns l_False if this method succeeds in reducing the assumptions to a provably minimal unsat core (the resulting unsat core will
// be stored in the supplied assumptions vector).
//Returns l_True if the assumptions are satisfiable.
//Returns l_Undef if solve time constraints prevent the assumptions from being reduced to a provably locally minimal unsat core
lbool minimizeCore(SimpSolver& S, vec<Lit>& assumptions, bool do_simp = false);
};
#endif /* OPTIMIZE_H_ */
