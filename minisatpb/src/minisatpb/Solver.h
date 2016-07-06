/****************************************************************************************[Solver.h]
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

#ifndef Solver_h
#define Solver_h

#include "MiniSat.h"
#include "SatELite.h"
#include "Global.h"
#include "Main.h"
//=================================================================================================


class Solver {
    MiniSat ::Solver*   minisat;
    SatELite::Solver*   satelite;
public:
    bool&             ok_ref     () { return (minisat != NULL) ? minisat->ok      : satelite->ok     ; }
    vec<int>&         assigns_ref() { return (minisat != NULL) ? minisat->assigns : satelite->assigns; }
    vec<Lit>&         trail_ref  () { return (minisat != NULL) ? minisat->trail   : satelite->trail  ; }
    BasicSolverStats& stats_ref  () { return (minisat != NULL) ? (BasicSolverStats&)minisat->stats : (BasicSolverStats&)satelite->stats  ; }

    void setVerbosity(int level) {
        if (minisat != NULL)
            minisat->verbosity = level;
        else
            satelite->verbosity = level; }

    Var         newVar         (bool dvar = true)   { return (minisat != NULL) ? minisat->newVar(dvar) : satelite->newVar(dvar); }
    bool        addClause      (const vec<Lit>& ps) { return (minisat != NULL) ? minisat->addClause(ps) : (satelite->addClause(ps), satelite->okay()); }
    bool        addUnit        (Lit p)              { return (minisat != NULL) ? minisat->addUnit(p) : (satelite->addUnit(p), satelite->okay()); }
    void        freeze         (Var x)              { if (minisat == NULL) satelite->freeze(x); }
    void        suggestPolarity(Var x, lbool value) { if (minisat != NULL) minisat->polarity_sug[x] = toInt(value); else satelite->polarity_sug[x] = toInt(value); }
    bool        solve     (const vec<Lit>& assumps) { return (minisat != NULL) ? minisat->solve(assumps) : satelite->solve(assumps); }
    bool        solve          ()                   { vec<Lit> tmp; return solve(tmp); }
    vec<lbool>& model          ()                   { return (minisat != NULL) ? minisat->model : satelite->model; }
    bool        varElimed      (Var x)              { return (minisat != NULL) ? false : satelite->var_elimed[x]; }
    bool        okay           ()                   { return (minisat != NULL) ? minisat->okay() : satelite->okay(); }
    int         nVars          ()                   { return (minisat != NULL) ? minisat->nVars() : satelite->nVars(); }
    void        exportCnf      (cchar* filename)    {
        if (minisat != NULL){
            minisat->simplifyDB();
            minisat->exportClauses(filename);
        }else{

            SatELite::opt_pre_sat = true;
            SatELite::output_file = filename;
            if (opt_verbosity >= 1) reportf("=================================[SATELITE+]==================================\n");
            satelite->simplifyDB(true);
        } }
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    Solver(bool use_minisat) : minisat(use_minisat ? new MiniSat::Solver : NULL), satelite(use_minisat ? NULL : new SatELite::Solver) {}
};


//=================================================================================================
#endif
