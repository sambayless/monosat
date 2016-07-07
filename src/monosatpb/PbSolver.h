/**************************************************************************************[PbSolver.h]
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

#ifndef PbSolver_h
#define PbSolver_h

#include "monosat/mtl/Vec.h"
#include "monosat/simp/SimpSolver.h"

#include "monosatpb/ADTs/Map.h"
#include "monosatpb/ADTs/StackAlloc.h"
#include "monosatpb/ADTs/Int.h"
#include "monosatpb/Config_pb.h"

namespace Monosat {
namespace PB{
using Monosat::Var;
using Monosat::Lit;
using Monosat::SimpSolver;
using Monosat::lbool;
using Monosat::mkLit;
using Monosat::lit_Undef;
/*using Monosat::l_Undef;
using Monosat::l_False;
using Monosat::l_True;
using Monosat::var_Undef;*/
//=================================================================================================
// Linear -- a class for storing pseudo-boolean constraints:


class Linear {
    int orig_size;  // Allocated terms in constraint.
public:
    int size;       // Terms in constraint.
    Int lo, hi;     // Sum should be in interval [lo,hi] (inclusive).
private:
    char data[0];    // (must be last element of the struct)
public:
    // NOTE: Cannot be used by normal 'new' operator!
    Linear(const vec<Lit> &ps, const vec<Int> &Cs, Int low, Int high) {
        orig_size = size = ps.size(), lo = low, hi = high;
        char *p = data;
        for (int i = 0; i < ps.size(); i++) *(Lit *) p = ps[i], p += sizeof(Lit);
        for (int i = 0; i < Cs.size(); i++) new((Int *) p) Int(Cs[i]), p += sizeof(Int);
    }

    ~Linear() {
        for (int i = 0; i < size; i++)
            (*this)(i).~Int();
    }

    Lit operator[](int i) const { return *(Lit *) (data + sizeof(Lit) * i); }

    Int operator()(int i) const { return *(Int *) (data + sizeof(Lit) * orig_size + sizeof(Int) * i); }

    Lit &operator[](int i) { return *(Lit *) (data + sizeof(Lit) * i); }

    Int &operator()(int i) { return *(Int *) (data + sizeof(Lit) * orig_size + sizeof(Int) * i); }
};


//=================================================================================================
// PbSolver -- Pseudo-boolean solver (linear boolean constraints):


class PbSolver {
protected:
    SimpSolver & sat_solver;     // Underlying SAT solver.
    vec<Lit> trail;          // Chronological assignment stack.

    StackAlloc<char *> mem;            // Used to allocate the 'Linear' constraints stored in 'constrs' (other 'Linear's, such as the goal function, are allocated with 'xmalloc()')

public:
    vec<Linear *> constrs;        // Vector with all constraints.
    Linear *goal;           // Non-normalized goal function (used in optimization). NULL means no goal function specified. NOTE! We are always minimizing.
protected:
    vec<int> n_occurs;       // Lit -> int: Number of occurrences.
    vec<vec<int> > occur;          // Lit -> vec<int>: Occur lists. Left empty until 'setupOccurs()' is called.

    int propQ_head;     // Head of propagation queue (index into 'trail').
    Monosat::vec<Lit> tmp_clause;

    // Main internal methods:
    //
    bool propagate(Linear &c);

    void propagate();

    bool addUnit(Lit p) {
        if (value(p) == l_Undef) trail.push(p);
        return sat_solver.addClause(p);
    }

    bool addClause(const vec<Lit> &ps) {
        tmp_clause.clear();
        for (int i = 0; i < ps.size(); i++) tmp_clause.push(ps[i]);
        return sat_solver.addClause_(tmp_clause);
    }

    bool normalizePb(vec<Lit> &ps, vec<Int> &Cs, Int &C);

    void storePb(const vec<Lit> &ps, const vec<Int> &Cs, Int lo, Int hi);

    void setupOccurs();   // Called on demand from 'propagate()'.
    void findIntervals();

    bool rewriteAlmostClauses();

    bool convertPbs(bool first_call);   // Called from 'solve()' to convert PB constraints to clauses.

public:
    PbSolver(SimpSolver & sat_solver, bool use_preprocessing = false)
            : sat_solver(sat_solver), goal(NULL), propQ_head(0)
            //, stats(sat_solver.stats_ref())
            , declared_n_vars(-1), declared_n_constrs(-1), best_goalvalue(Int_MAX) {
        // Turn off preprocessing if wanted.
        if (!use_preprocessing)
            sat_solver.eliminate(true);
    }

    // Helpers (semi-internal):
    //
    lbool value(Var x) const { return sat_solver.value(x); }

    lbool value(Lit p) const { return sat_solver.value(p); }

    int nVars() const { return sat_solver.nVars(); }

    int nConstrs() const { return constrs.size(); }

    // Public variables:
    //BasicSolverStats& stats;
    void printStats();

    int declared_n_vars;            // Number of variables declared in file header (-1 = not specified).
    int declared_n_constrs;         // Number of constraints declared in file header (-1 = not specified).
    int pb_n_vars;                  // Actual number of variables (before clausification).
    int pb_n_constrs;               // Actual number of constraints (before clausification).

    Map<cchar *, int> name2index;
    vec<cchar *> index2name;
    vec<bool> best_model;     // Best model found (size is 'pb_n_vars').
    Int best_goalvalue; // Value of goal function for that model (or 'Int_MAX' if no models were found).

    // Problem specification:
    //
    int getVar(cchar *name);

    void allocConstrs(int n_vars, int n_constrs);

    void addGoal(const vec<Lit> &ps, const vec<Int> &Cs);

    bool addConstr(const vec<Lit> &ps, const vec<Int> &Cs, Int rhs, int ineq);

    // Solve:
    //
    bool okay(void) { return sat_solver.okay(); }

    enum solve_Command {
        sc_Minimize, sc_FirstSolution, sc_AllSolutions
    };

    void solve(solve_Command cmd = sc_Minimize);    // Returns best/first solution found or Int_MAX if UNSAT.
};
}
}
//=================================================================================================
#endif
