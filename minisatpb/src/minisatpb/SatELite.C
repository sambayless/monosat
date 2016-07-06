/**************************************************************************************[SatELite.C]
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

/**************************************************************************************************

A simple Chaff-like SAT-solver with support for incremental SAT.

**************************************************************************************************/

#include "SatELite.h"
#include "Sort.h"
#include <cmath>
#include "Global.h"
#include "Main.h"
namespace SatELite {

//#################################################################################################
// "Main.C" BEGIN
//#################################################################################################
bool opt_confl_1sub   = true;
bool opt_confl_ksub   = false;
bool opt_var_elim     = true;
bool opt_0sub         = true;
bool opt_1sub         = true;
bool opt_2sub         = false;
bool opt_repeated_sub = false;
bool opt_def_elim     = true;
bool opt_unit_def     = true;
bool opt_hyper1_res   = true;
bool opt_pure_literal = false;
bool opt_asym_branch  = false;
bool opt_keep_all     = false;
bool opt_no_random    = false;
bool opt_pre_sat      = false;
bool opt_ext_sat      = false;
bool opt_niver        = false;
cchar* input_file    = NULL;
cchar* output_file   = NULL;
cchar* varmap_file   = NULL;
cchar* elimed_file   = NULL;
cchar* model_file    = NULL;
extern int verbosity;
//#################################################################################################
// "Main.C" END
//#################################################################################################


//#################################################################################################
// "BcnfWrite.iC" BEGIN
//#################################################################################################
//
// NOTE! Must include "SolverTypes.h" before including this file.
//

#define BCNF_CHUNK_LIMIT 1048576


//=================================================================================================
// BCNF Writer:


class BcnfWriter {
    FILE*   out;
    int     n_vars;
    int     n_clauses;
    int     chunk_sz;
    int*    chunk;

public:
    int     stated_n_vars;      // }- The "p cnf" line.
    int     stated_n_clauses;   // }

    BcnfWriter(cchar* output_file);
   ~BcnfWriter(void);

    void addClause(vec<Lit>& c);
    int  nVars   (void) { return n_vars; }
    int  nClauses(void) { return n_clauses; }
};


BcnfWriter::BcnfWriter(cchar* output_file)
    : n_vars(0), n_clauses(0), chunk_sz(1), stated_n_vars(-1), stated_n_clauses(-1)
{
    out = fopen(output_file, "w+b");
    if (out == NULL) fprintf(stderr, "ERROR! Could not open file for writing: %s\n", output_file), exit(2);   // <<= Exception handling
    fputc(0, out); fputc(0, out); fputc(0, out); fputc(0, out);     // Room for header: "BCNF"
    fputc(0, out); fputc(0, out); fputc(0, out); fputc(0, out);     // Room for byte-order: 1,2,3,4
    fputc(0, out); fputc(0, out); fputc(0, out); fputc(0, out);     // Room for #variables
    fputc(0, out); fputc(0, out); fputc(0, out); fputc(0, out);     // Room for #clauses

    chunk = xmalloc<int>(BCNF_CHUNK_LIMIT);
}


BcnfWriter::~BcnfWriter(void)
{
    chunk[0] = chunk_sz;
    chunk[chunk_sz++] = -1;
    fwrite(chunk, 4, chunk_sz, out);
    xfree(chunk);

    fflush(out);
    rewind(out);
    int byte_order = 0x01020304;
    fputc('B', out); fputc('C', out); fputc('N', out); fputc('F', out);
    fwrite(&byte_order, 1, 4, out);
    fwrite(&n_vars    , 1, 4, out);
    fwrite(&n_clauses , 1, 4, out);
    fclose(out);
}


void BcnfWriter::addClause(vec<Lit>& c)
{
    n_clauses++;

    if (chunk_sz + 3 + c.size() >= BCNF_CHUNK_LIMIT){  // leave room for final terminator size ("-1") and the chunk size itself  (not really part of the chunk, but just in case).
        chunk[0] = chunk_sz;
        chunk[chunk_sz++] = -1;
        fwrite(chunk, 4, chunk_sz, out);
        chunk_sz = 1; }

    assert(chunk_sz + 3 + c.size() < BCNF_CHUNK_LIMIT);
    chunk[chunk_sz++] = c.size();
    for (int i = 0; i < c.size(); i++)
        chunk[chunk_sz++] = index(c[i]),
        n_vars = max(n_vars, var(c[i]) + 1);
}
//#################################################################################################
// "BcnfWrite.iC" END
//#################################################################################################


//#################################################################################################
// "TmpFiles.C" BEGIN
//#################################################################################################
static vec<FILE*>  tmp_fps;
static vec<cchar*> tmp_files;


// 'out_name' should NOT be freed by caller.
FILE* createTmpFile(cchar* prefix, cchar* mode, char*& out_name)
{
    char*   name = xmalloc<char>(strlen(prefix) + 6 + 1);
    strcpy(name, prefix);
    strcat(name, "XXXXXX");

    int fd = mkstemp(name);
    if (fd == -1){
        fprintf(stderr, "ERROR! Could not create temporary file with prefix: %s\n", prefix);
        exit(1); }
    FILE* fp = fdopen(fd, mode); assert(fp != NULL);

    tmp_fps  .push(fp);
    tmp_files.push(name);
    if (&out_name != NULL) out_name = name;
    return fp;
}


// If 'exact' is set, 'prefix' is the full name returned by 'createTmpFile()'.
void deleteTmpFile(cchar* prefix, bool exact)
{
    uint    len = strlen(prefix);
    for (int i = 0; i < tmp_files.size();){
        if (!exact && (strlen(tmp_files[i]) == len + 6 && strncmp(tmp_files[i], prefix, len) == 0)
        ||   exact && (strcmp(tmp_files[i], prefix) == 0))
        {
            fclose(tmp_fps[i]);
            ::remove(tmp_files[i]);
            xfree(tmp_files[i]);
            tmp_fps[i] = tmp_fps.last();
            tmp_fps.pop();
            tmp_files[i] = tmp_files.last();
            tmp_files.pop();
        }else
            i++;
    }
}


// May be used to delete all temporaries. Will also be called automatically upon (normal) program exit.
void deleteTmpFiles(void)
{
    for (int i = 0; i < tmp_fps.size(); i++){
        fclose(tmp_fps[i]);
        ::remove(tmp_files[i]);
        xfree(tmp_files[i]);
    }
}


FINALIZER(TmpFiles) {
    deleteTmpFiles();
}
//#################################################################################################
// "TmpFiles.C" END
//#################################################################################################


//#################################################################################################
// "Solver_clause.iC" BEGIN
//#################################################################################################


//=================================================================================================
// Allocation:


#if 0
// HACKISH OPTIMIZATION OF MEMORY ALLOCATION:

#include "VecAlloc.h"
struct Size24 { char dummy[24]; };
struct Size28 { char dummy[28]; };
VecAlloc<Size24> mem24;
VecAlloc<Size28> mem28;

template <class T> macro void* ymalloc(int size);
template <> macro void* ymalloc<char>(int size)
{
    if      (size == 24) return (void*)mem24.alloc();
    else if (size == 28) return (void*)mem28.alloc();
    else                 return (void*)xmalloc<char>(size);
}


macro void yfree(void* ptr)
{
    Clause c((Clause_t*)ptr);
    assert(!c.dynamic());
    int size = sizeof(Clause_t) + sizeof(uint)*(c.size() + (int)c.learnt());
    if      (size == 24) mem24.free((Size24*)ptr);
    else if (size == 28) mem28.free((Size28*)ptr);
    else                 xfree(ptr);
}

#else
#define ymalloc xmalloc
#define yfree   xfree
#endif

//=================================================================================================
// Solver methods operating on clauses:


// Will allocate space in 'constrs' or 'learnts' and return the position to put new clause at.
int Solver::allocClauseId(bool learnt)
{
    int     id;
    if (learnt){
        if (learnts_free.size() > 0)
            id = learnts_free.last(),
            learnts_free.pop();
        else
            id = learnts.size(),
            learnts.push();
    }else{
        if (constrs_free.size() > 0)
            id = constrs_free.last(),
            constrs_free.pop();
        else
            id = constrs.size(),
            constrs.push();
    }
    return id;
}


void Solver::freeClauseId(int id, bool learnt)
{
    if (learnt){
        learnts[id] = Clause_NULL;
        learnts_free.push(id);
    }else{
        constrs[id] = Clause_NULL;
        constrs_free.push(id);
    }
}


Clause Solver::allocClause(const vec<Lit>& ps, bool learnt, Clause overwrite)
{
    assert(sizeof(Lit)   == sizeof(uint));
    assert(sizeof(float) == sizeof(uint));
    int       id  = overwrite.null() ? allocClauseId(learnt) : overwrite.id();
    void*     mem = overwrite.null() ? ymalloc<char>(sizeof(Clause_t) + sizeof(uint)*(ps.size() + (int)learnt)) : (void*)overwrite.ptr();
    Clause_t* c   = new (mem) Clause_t;

    c->id_         = id;
    c->abst_       = 0;
    c->size_learnt = (int)learnt | (ps.size() << 1);
    for (int i = 0; i < ps.size(); i++){
        c->data[i] = ps[i];
        c->abst_  |= abstLit(ps[i]);
    }
    if (learnt) Clause(c).activity() = 0.0;

    return Clause(c);
}


// PUBLIC: May return NULL if clause is already satisfied by the top-level assignment.
//
Clause Solver::addClause(const vec<Lit>& ps_, bool learnt, Clause overwrite)
{
    if (!ok) return Clause_NULL;

    vec<Lit>    qs;

    // Contains no eliminated variables?
    for (int i = 0; i < ps_.size(); i++)
        assert(!var_elimed[var(ps_[i])]);

    if (!learnt){
        assert(decisionLevel() == 0);
        ps_.copyTo(qs);             // Make a copy of the input vector.

        // Remove false literals:
        for (int i = 0; i < qs.size();){
            if (value(qs[i]) != l_Undef){
                if (value(qs[i]) == l_True){
                    if (!overwrite.null()) deallocClause(overwrite);
                    return Clause_NULL; }  // Clause always true -- don't add anything.
                else
                    qs[i] = qs.last(),
                    qs.pop();
            }else
                i++;
        }

        // Remove duplicates:
        sortUnique(qs);
        for (int i = 0; i < qs.size()-1; i++){
            if (qs[i] == ~qs[i+1]){
                if (!overwrite.null()) deallocClause(overwrite);
                return Clause_NULL;        // Clause always true -- don't add anything.
            }
        }
    }
    const vec<Lit>& ps = learnt ? ps_ : qs;     // 'ps' is now the (possibly) reduced vector of literals.
    //*HACK!*/if (ps.size() <= 3) learnt = false;
    if (opt_keep_all) learnt = false;

    if (ps.size() == 0){
        if (!overwrite.null()) deallocClause(overwrite);
        ok = false;
        return Clause_NULL;
    }else if (ps.size() == 1){
        if (!overwrite.null()) deallocClause(overwrite);
        if (decisionLevel() > 0) units.push(ps[0]);
        if (!enqueue(ps[0])){
            assert(decisionLevel() == 0);
            ok = false; }
        return Clause_NULL; }
    else{
        // Allocate clause:
        Clause c = allocClause(ps, learnt, overwrite);

        // Subsumption:
        if (occur_mode != occ_Off){
            if (fwd_subsump && isSubsumed(c)){
                deallocClause(c);
                return Clause_NULL; }
            if (!learnt) subsume0(c);
        }

        // Occur lists:
        if (updateOccur(c)){
            for (int i = 0; i < ps.size(); i++)
                occur[index(ps[i])].push(c),
                touch(var(ps[i]));
            if (overwrite.null())
                cl_added.add(c);
            else
                cl_touched.add(c);
        }

        // For learnt clauses only:
        if (learnt){
            // Put the second watch on the literal with highest decision level:
            int     max_i = 1;
            int     max   = level[var(ps[1])];
            for (int i = 2; i < ps.size(); i++)
                if (level[var(ps[i])] > max)
                    max   = level[var(ps[i])],
                    max_i = i;
            c[1]     = ps[max_i];
            c[max_i] = ps[1];

            // Bumping:
            claBumpActivity(c); // (newly learnt clauses should be considered active)
        }
        //*TEST*/for (int i = 0; i < c.size(); i++) varBumpActivity(c[i]);

        // Attach clause:
        if (watches_setup)
            watch(c, ~c[0]),
            watch(c, ~c[1]);

        if (learnt) learnts[c.id()] = c;
        else        constrs[c.id()] = c;

        if (!learnt){
            n_literals += c.size();
            for (int i = 0; i < c.size(); i++)
                n_occurs[index(c[i])]++;
        }

        return c;
    }
}


void Solver::deallocClause(Clause c, bool quick)    // (quick means only free memory; used in destructor of 'Solver')
{
    if (quick)
        yfree(c.ptr());
    else{
        // Free resources:
        freeClauseId(c.id(), c.learnt());
        yfree(c.ptr());
    }
}


void Solver::unlinkClause(Clause c, Var elim)
{
    if (elim != var_Undef){
        assert(!c.learnt());
        io_tmp.clear();
        io_tmp.push(toLit(c.size()));
        for (int i = 0; i < c.size(); i++)
            io_tmp.push(c[i]);
        fwrite((Lit*)io_tmp, 4, io_tmp.size(), elim_out);
    }

    if (updateOccur(c)){
        for (int i = 0; i < c.size(); i++){
            maybeRemove(occur[index(c[i])], c);
          #ifndef TOUCH_LESS
            touch(c[i]);
          #endif
        }
    }

    if (watches_setup)
        unwatch(c, ~c[0]),
        unwatch(c, ~c[1]);
    if (c.learnt()) learnts[c.id()] = Clause_NULL;
    else            constrs[c.id()] = Clause_NULL;

    if (!c.learnt()){
        n_literals -= c.size();
        for (int i = 0; i < c.size(); i++)
            n_occurs[index(c[i])]--;
    }

    // Remove from iterator vectors/sets:
    for (int i = 0; i < iter_vecs.size(); i++){
        vec<Clause>& cs = *iter_vecs[i];
        for (int j = 0; j < cs.size(); j++)
            if (cs[j] == c)
                cs[j] = NULL;
    }
    for (int i = 0; i < iter_sets.size(); i++){
        CSet& cs = *iter_sets[i];
        cs.exclude(c);
    }

    // Remove clause from clause touched set:
    if (updateOccur(c))
        cl_touched.exclude(c),
        cl_added  .exclude(c);

}


// 'p' is the literal that became TRUE. Returns FALSE if conflict found. 'keep_watch' should be set
//  to FALSE before call. It will be re-set to TRUE if the watch should be kept in the current
//  watcher list.
//
bool Solver::propagateClause(Clause c, Lit p, bool& keep_watch)
{
    assert(watches_setup);

    // Make sure the false literal is c[1]:
    Lit     false_lit = ~p;
    if (c(0) == false_lit)
        c(0) = c(1), c(1) = false_lit;
    assert(c(1) == false_lit);

    // If 0th watch is true, then clause is already satisfied.
    if (value(c(0)) == l_True){
        keep_watch = true;
        return true; }

    // Look for new watch:
    for (int i = 2; i < c.size(); i++){
        if (value(c(i)) != l_False){
            c(1) = c(i), c(i) = false_lit;
            watches[index(~c(1))].push(c);
            return true; } }

    // Clause is unit under assignment:
    keep_watch = true;
    return enqueue(c(0), c);
}


// Can assume 'out_reason' to be empty.
// Calculate reason for 'p'. If 'p == lit_Undef', calculate reason for conflict.
//
void Solver::calcReason(Clause c, Lit p, vec<Lit>& out_reason)
{
    assert(p == lit_Undef || p == c[0]);
    for (int i = ((p == lit_Undef) ? 0 : 1); i < c.size(); i++)
        assert(value(c[i]) == l_False),
        out_reason.push(~c[i]);
    if (c.learnt()) claBumpActivity(c);
}


// Will remove clause and add a new shorter, potentially unit and thus adding facts to propagation queue.
//
void Solver::strengthenClause(Clause c, Lit p)
{
    assert(c.size() > 1);

    vec<Lit>    new_clause;
    for (int i = 0; i < c.size(); i++)
        if (c[i] != p && (value(c[i]) != l_False || level[var(c[i])] > 0))
            new_clause.push(c[i]);

    unlinkClause(c);    // (will touch all variables)
    addClause(new_clause, c.learnt(), c);
}


//=================================================================================================
// Optimized occur table builder:


void Solver::setOccurMode(OccurMode new_occur_mode)
{
    assert(new_occur_mode != occ_All);      // (not implemented)
    occur_mode = new_occur_mode;

    if (occur_mode == occ_Off)
        for (int i = 0; i < occur.size(); i++)
            occur[i].clear(true);
    else{
        assert(occur_mode == occ_Permanent);
        // Allocate vectors of right capacities:
        for (int i = 0; i < nVars()*2; i++){
            vec<Clause> tmp(xmalloc<Clause>(n_occurs[i]), n_occurs[i]); tmp.clear();
            tmp.moveTo(occur[i]); }
        // Fill vectors:
        for (int i = 0; i < constrs.size(); i++){
            if (constrs[i].null()) continue;
            for (int j = 0; j < constrs[i].size(); j++)
                assert(occur[index(constrs[i][j])].size() < n_occurs[index(constrs[i][j])]),
                occur[index(constrs[i][j])].push(constrs[i]);
        }
    }
}


void Solver::setupWatches(void)
{
    if (watches_setup) return;
    assert(learnts.size() == 0);

    for (int i = 0; i < constrs.size(); i++)
        if (!constrs[i].null())
            watch(constrs[i], ~constrs[i][0]),
            watch(constrs[i], ~constrs[i][1]);
    watches_setup = true;
}
//#################################################################################################
// "Solver_clause.iC" END
//#################################################################################################



//=================================================================================================
// Minor methods:


// Creates a new SAT variable in the solver. If 'decision_var' is cleared, variable will not be
// used as a decision variable (NOTE! This has effects on the meaning of a SATISFIABLE result).
//
Var Solver::newVar(bool dvar)
{
    int     index;
    index = nVars();
    watches     .push();        // (list for positive literal)
    watches     .push();        // (list for negative literal)
    occur       .push();
    occur       .push();
    n_occurs    .push(0);
    n_occurs    .push(0);
    reason      .push(Clause_NULL);
    assigns     .push(toInt(l_Undef));
    level       .push(-1);
    activity    .push(0);
    polarity_sug.push(toInt(l_Undef));
    order       .newVar(dvar);
    seen_tmp    .push(0);       // (one for each polarity)
    seen_tmp    .push(0);
    touched     .push(1);
    touched_list.push(index);
    touched_tmp .push(0);
    var_elimed  .push(0);
    frozen      .push(0);
    return index;
}


// Returns FALSE if immediate conflict.
bool Solver::assume(Lit p) {
    assert(propQ.size() == 0);
    if (verbosity >= 2) reportf(L_IND"assume("L_LIT")\n", L_ind, L_lit(p));
    trail_lim.push(trail.size());
    return enqueue(p); }


// Revert one variable binding on the trail.
//
inline void Solver::undoOne(void)
{
    if (verbosity >= 2){ Lit p = trail.last(); reportf(L_IND"unbind("L_LIT")\n", L_ind, L_lit(p)); }
    Lit     p  = trail.last(); trail.pop();
    Var     x  = var(p);
    assigns[x] = toInt(l_Undef);
    reason [x] = Clause_NULL;
    level  [x] = -1;
    order.undo(x);
}


// Reverts to the state before last 'assume()'.
//
void Solver::cancel(void)
{
    assert(propQ.size() == 0);
    if (verbosity >= 2){ if (trail.size() != trail_lim.last()){ Lit p = trail[trail_lim.last()]; reportf(L_IND"cancel("L_LIT")\n", L_ind, L_lit(p)); } }
    for (int c = trail.size() - trail_lim.last(); c != 0; c--)
        undoOne();
    trail_lim.pop();
}


// Revert to the state at given level.
//
void Solver::cancelUntil(int level) {
    while (decisionLevel() > level) cancel(); }


// Record a clause and drive backtracking. 'clause[0]' must contain the asserting literal.
//
void Solver::record(const vec<Lit>& clause)
{
    assert(clause.size() != 0);
    Clause c = addClause(clause, true); assert(ok);
    check(enqueue(clause[0], c));
}


//=================================================================================================
// Major methods:


Solver::~Solver(void)
{
    assert(iter_vecs.size() == 0); assert(iter_sets.size() == 0);
    for (int i = 0; i < constrs.size(); i++) if (!constrs[i].null()) deallocClause(constrs[i], true);
    for (int i = 0; i < learnts.size(); i++) if (!learnts[i].null()) deallocClause(learnts[i], true);
    deleteTmpFiles();
}


/*_________________________________________________________________________________________________
|                                                                                                  
|  analyze : (confl : Clause) (out_learnt : vec<Lit>&) (out_btlevel : int&)  ->  [void]            
|                                                                                                  
|  Description:                                                                                    
|    Analyze conflict and produce a reason clause.                                                 
|                                                                                                  
|    Pre-conditions:                                                                               
|      * 'out_learnt' is assumed to be cleared.                                                    
|      * Current decision level must be greater than root level.                                   
|                                                                                                  
|    Post-conditions:                                                                              
|      * 'out_learnt[0]' is the asserting literal at level 'out_btlevel'.                          
|                                                                                                  
|  Effect:                                                                                         
|    Will undo part of the trail, upto but not beyond the assumption of the current decision level.
|________________________________________________________________________________________________@*/
void Solver::analyze(Clause confl, vec<Lit>& out_learnt, int& out_btlevel)
{
    vec<char>&  seen = seen_tmp;
    int         pathC    = 0;
    Lit         p = lit_Undef;
    vec<Lit>    p_reason;

    // Generate conflict clause:
    //
#if 1
// Niklas S�rensson's version
    // Generate conflict clause:
    //
    out_learnt.push();      // (leave room for the asserting literal)
    out_btlevel = 0;
    int index = trail.size()-1;
    do{
        assert(confl != Clause_NULL);          // (otherwise should be UIP)
        if (confl.learnt()) claBumpActivity(confl);

        for (int j = p == lit_Undef ? 0 : 1; j < confl.size(); j++){
            Lit q = confl(j);
            if (!seen[var(q)] && level[var(q)] > 0){
                seen[var(q)] = 1;
                varBumpActivity(q);
                if (level[var(q)] == decisionLevel())
                    pathC++;
                else{
                    out_learnt.push(q),
                    out_btlevel = max(out_btlevel, level[var(q)]);
                }
            }
        }

        // Select next clause to look at:
        while (!seen[var(trail[index--])]);
        p     = trail[index+1];
        confl = reason[var(p)];
        seen[var(p)] = 0;
        pathC--;

    }while (pathC > 0);
    out_learnt[0] = ~p;

#else
    out_learnt.push();      // (leave room for the asserting literal)
    out_btlevel = 0;
    do{
        assert(!confl.null());          // (otherwise should be UIP)

        p_reason.clear();
        calcReason(confl, p, p_reason);

        for (int j = 0; j < p_reason.size(); j++){
            Lit q = p_reason[j];
            if (!seen[var(q)] && level[var(q)] > 0){
              #ifdef BUMP_MORE
                varBumpActivity(q);
              #endif
                seen[var(q)] = 1;
                if (level[var(q)] == decisionLevel())
                    pathC++;
                else{
                    out_learnt.push(~q),
                    out_btlevel = max(out_btlevel, level[var(q)]);
                }
            }
        }

        // Select next clause to look at:
        do{
            p = trail.last();
            confl = reason[var(p)];
            undoOne();
        }while (!seen[var(p)]);
        pathC--;
        seen[var(p)] = 0;
    }while (pathC > 0);
    out_learnt[0] = ~p;
#endif

  #ifndef BUMP_MORE
    // Bump variables:
    for (int i = 0; i < out_learnt.size(); i++)
        varBumpActivity(out_learnt[i]);
  #endif

    // Remove literals:
//#define DEBUG_CSK
  #ifdef DEBUG_CSK
    vec<char>   seen_copy; seen.copyTo(seen_copy);
    vec<Lit>    out_learnt_copy; out_learnt.copyTo(out_learnt_copy);
    vec<Lit>    candidate;
  #endif

    vec<Var>    to_clear(out_learnt.size()); for (int i = 0; i < out_learnt.size(); i++) to_clear[i] = var(out_learnt[i]);
    if (opt_confl_ksub){
        vec<int>    levels;
        for (int i = 1; i < out_learnt.size(); i++)
            levels.push(level[var(out_learnt[i])]);
        sortUnique(levels);

        for (int k = 0; k < levels.size(); k++){
            assert(levels[k] > 0);
            int from = trail_lim[levels[k]-1];
            int to   = (levels[k] >= trail_lim.size()) ? trail.size() : trail_lim[levels[k]];
            for (int i = from; i < to; i++){
                Var     x = var(trail[i]);
                Clause  r = reason[x];
                if (!r.null()){
                    for (int j = 1; j < r.size(); j++)
                        if (seen[var(r[j])] == 0 && level[var(r[j])] != 0)
                            goto NoConsequence;
                    seen[x] = 2;
                    to_clear.push(x);
                  NoConsequence:;
                }
            }
        }
        for (int i = 0; i < out_learnt.size(); i++)
            if (seen[var(out_learnt[i])] == 2)
                out_learnt[i] = lit_Undef;

      #ifdef DEBUG_CSK
        for (int i = 0; i < out_learnt.size(); i++)
            if (out_learnt[i] != lit_Undef)
                candidate.push(out_learnt[i]);
        seen_copy.copyTo(seen);
        out_learnt_copy.copyTo(out_learnt);

        for (int i = 0; i < out_learnt.size(); i++){
            Clause c = reason[var(out_learnt[i])];
            if (!c.null()){
                assert(c[0] == ~out_learnt[i]);
                for (int j = 1; j < c.size(); j++){
                    if (seen[var(c[j])] != 1)
                        goto Done2;
                }
                out_learnt[i] = lit_Undef;
              Done2:;
            }
        }
      #endif

    }else if (opt_confl_1sub){
        for (int i = 0; i < out_learnt.size(); i++){
            Clause c = reason[var(out_learnt[i])];
            if (!c.null()){
                assert(c[0] == ~out_learnt[i]);
                for (int j = 1; j < c.size(); j++){
                    if (seen[var(c[j])] != 1)
                        goto Done;
                }
                out_learnt[i] = lit_Undef;
              Done:;
            }
        }
    }
    if (opt_confl_1sub || opt_confl_ksub){
        int new_sz = 0;
        for (int i = 0; i < out_learnt.size(); i++)
            if (out_learnt[i] != lit_Undef)
                out_learnt[new_sz++] = out_learnt[i];
        out_learnt.shrink(out_learnt.size() - new_sz);
    }
  #ifdef DEBUG_CSK
    if (candidate.size() > out_learnt.size()){
        reportf("csk: "), dump(candidate);
        reportf("cs1: "), dump(out_learnt);
        exit(0);
    }
  #endif


    // Clear 'seen':
    for (int j = 0; j < to_clear.size(); j++) seen[to_clear[j]] = 0;    // ('seen[]' is now cleared)

    if (verbosity >= 2){
        reportf(L_IND"Learnt {", L_ind);
        for (int i = 0; i < out_learnt.size(); i++) reportf(" "L_LIT, L_lit(out_learnt[i]));
        reportf(" } at level %d\n", out_btlevel); }
}


/*_________________________________________________________________________________________________
|                                                                                                  
|  enqueue : (p : Lit) (from : Clause)  ->  [bool]                                                 
|                                                                                                  
|  Description:                                                                                    
|    Puts a new fact on the propagation queue as well as immediately updating the variable's value.
|    Should a conflict arise, FALSE is returned.                                                   
|                                                                                                  
|  Input:                                                                                          
|    p    - The fact to enqueue                                                                    
|    from - [Optional] Fact propagated from this (currently) unit clause. Stored in 'reason[]'.    
|           Default value is NULL (no reason).                                                     
|                                                                                                  
|  Output:                                                                                         
|    TRUE if fact was enqueued without conflict, FALSE otherwise.                                  
|________________________________________________________________________________________________@*/
bool Solver::enqueue(Lit p, Clause from)
{
    if (value(p) != l_Undef){
      #ifdef RELEASE
        return value(p) != l_False;
      #else
        if (value(p) == l_False){
            // Conflicting enqueued assignment
            assert(decisionLevel() > 0);
            return false;
        }else if (value(p) == l_True){
            // Existing consistent assignment -- don't enqueue
            return true;
        }else{
            assert(value(p) == l_Error);
            // Do nothing -- clause should be removed.
            return true;
        }
      #endif
    }else{
        // New fact -- store it.
      #ifndef RELEASE
        if (verbosity >= 2) reportf(L_IND"bind("L_LIT")\n", L_ind, L_lit(p));
      #endif
        assigns[var(p)] = toInt(lbool(!sign(p)));
        level  [var(p)] = decisionLevel();
        reason [var(p)] = from;
        trail.push(p);
        propQ.insert(p);
        return true;
    }
}


/*_________________________________________________________________________________________________
|                                                                                                  
|  propagateToplevel : [void]  ->  [void]                                                          
|                                                                                                  
|  Description:                                                                                    
|    Destructively update clause database with enqueued top-level facts.                           
|________________________________________________________________________________________________@*/
void Solver::propagateToplevel(void)
{
    assert(decisionLevel() == 0);

    for (int i = 0; i < units.size(); i++){
        if (!enqueue(units[i])){
            propQ.clear();
            return; } }
    units.clear();

    while (propQ.size() > 0){
        Lit            p = propQ.dequeue();
        vec<Clause>    cs;

        // Remove satisfied clauses:
        occur[index(p)].moveTo(cs);
        for (int i = 0; i < cs.size(); i++)
            removeClause(cs[i]);

        // Remove false literals from clauses:
        occur[index(~p)].moveTo(cs);
        registerIteration(cs);
        for (int i = 0; i < cs.size(); i++){
            if (!cs[i].null())
                strengthenClause(cs[i], ~p);    // (may enqueue new facts to propagate)
        }
        unregisterIteration(cs);
    }
    propQ.clear();
}


/*_________________________________________________________________________________________________
|                                                                                                  
|  propagate : [void]  ->  [Clause]                                                                
|                                                                                                  
|  Description:                                                                                    
|    Propagates all enqueued facts. If a conflict arises, the conflicting clause is returned,      
|    otherwise NULL.                                                                               
|                                                                                                  
|    Post-conditions:                                                                              
|      * the propagation queue is empty, even if there was a conflict.                             
|________________________________________________________________________________________________@*/
#if 0
// Standard
Clause Solver::propagate(void)
{
    if (decisionLevel() == 0 && occur_mode != occ_Off){
        propagateToplevel();
        return Clause_NULL;
    }

    Clause confl;
    while (propQ.size() > 0){
        stats.propagations++;
        Lit           p  = propQ.dequeue();        // 'p' is enqueued fact to propagate.
        vec<Clause>&  ws = watches[index(p)];
        bool          keep_watch;
        int           i, j;
        for (i = j = 0; confl.null() && i < ws.size(); i++){
            stats.inspects++;
            keep_watch = false;
            if (!propagateClause(ws[i], p, keep_watch))
                confl = ws[i],
                propQ.clear();
            if (keep_watch)
                ws[j++] = ws[i];
        }

        // Copy the remaining watches:
        while (i < ws.size())
            ws[j++] = ws[i++];

        watches[index(p)].shrink(i - j);
    }

    return confl;
}
#endif

#if 0
// Optimized
Clause Solver::propagate(void)
{
    if (decisionLevel() == 0 && occur_mode != occ_Off){
        propagateToplevel();
        return Clause_NULL;
    }

    Clause confl;
    while (propQ.size() > 0){
        stats.propagations++;
        Lit           p  = propQ.dequeue();        // 'p' is enqueued fact to propagate.
        vec<Clause>&  ws = watches[index(p)];
        int           i, j;
        for (i = j = 0; i < ws.size(); i++){
            stats.inspects++;

            // Make sure the false literal is c[1]:
            Clause& c = ws[i];
            Lit     false_lit = ~p;
            if (c(0) == false_lit)
                c(0) = c(1), c(1) = false_lit;
            assert(c(1) == false_lit);

            // If 0th watch is true, then clause is already satisfied.
            if (value(c(0)) == l_True){
                ws[j++] = ws[i];
                goto Continue; }

            // Look for new watch:
            for (int k = 2; k < c.size(); k++){
                if (value(c(k)) != l_False){
                    c(1) = c(k), c(k) = false_lit;
                    watches[index(~c(1))].push(c);
                    goto Continue; } }

            // Clause is unit under assignment:
            ws[j++] = ws[i];
            if (!enqueue(c(0), c)){
                confl = ws[i],
                propQ.clear();
                i++; break;
            }

          Continue:;
        }

        // Copy the remaining watches:
        while (i < ws.size())
            ws[j++] = ws[i++];

        watches[index(p)].shrink(i - j);
    }

    return confl;
}
#endif

#if 1
// Borrowed from Niklas S�rensson -- uses "unsafe" type casts to achieve maximum performance
Clause Solver::propagate(void)
{
    if (decisionLevel() == 0 && occur_mode != occ_Off){
        propagateToplevel();
        return Clause_NULL;
    }
    assert(watches_setup);

    Clause confl = Clause_NULL;
    while (propQ.size() > 0){
        stats.propagations++;
        Lit           p  = propQ.dequeue();        // 'p' is enqueued fact to propagate.
        vec<Clause>&  ws = watches[index(p)];
        Clause_t      **i, **j, **end = (Clause_t**)(Clause*)ws + ws.size();
        for (i = j = (Clause_t**)(Clause*)ws; confl == NULL && i < end; ){
            stats.inspects++;
            Clause_t& c = **i;

            // Make sure the false literal is data[1]:
            Lit false_lit = ~p;
            if (c(0) == false_lit)
                c(0) = c(1), c(1) = false_lit;

            assert(c(1) == false_lit);

            // If 0th watch is true, then clause is already satisfied.
            if (value(c(0)) == l_True)
                *j++ = *i;
            else{
                // Look for new watch:
                for (int k = 2; k < c.size(); k++)
                    if (value(c(k)) != l_False){
                        c(1) = c(k); c(k) = false_lit;
                        watches[index(~c(1))].push(&c);
                        goto next; }

                // Clause is unit under assignment:
                *j++ = *i;
                if (!enqueue(c(0), &c)){
                    confl = *(Clause*)i;
                    propQ.clear();
                }
            }
        next:
            i++;
        }

        // Copy the remaining watches:
        while (i < end)
            *j++ = *i++;

        ws.shrink(i - j);
    }

    return confl;
}
#endif


/*_________________________________________________________________________________________________
|                                                                                                  
|  reduceDB : ()  ->  [void]                                                                       
|                                                                                                  
|  Description:                                                                                    
|    Remove half of the learnt clauses, minus the clauses locked by the current assignment. Locked 
|    clauses are clauses that are reason to a some assignment.                                     
|________________________________________________________________________________________________@*/
bool satisfied(Solver& S, Clause c)
{
    for (int i = 0; i < c.size(); i++){
        if ((S.value(c[i]) == l_True && S.level[var(c[i])] == 0) || S.value(c[i]) == l_Error)    // (l_Error means variable is eliminated)
            return true; }
    return false;
}

struct reduceDB_lt { bool operator () (Clause x, Clause y) { return x.activity() < y.activity(); } };
void Solver::reduceDB(void)
{
/**/if (occur_mode == occ_All) return;  // <<= Temporary fix

    // <<= BAD IMPLEMENTATION! FIX!
    stats.reduceDBs++;
    int     i;
    double  extra_lim = cla_inc / learnts.size();    // Remove any clause below this activity
    vec<Clause> ls;
    for (int i = 0; i < learnts.size(); i++)
        if (!learnts[i].null())
             ls.push(learnts[i]);

    sort(ls,
    reduceDB_lt());
    for (i = 0; i < ls.size() / 3; i++){
        if (!locked(ls[i]))
            removeClause(ls[i]);
    }
    for (; i < ls.size(); i++){
//      if (!locked(ls[i]) && ls[i].activity() < extra_lim)
        if (!locked(ls[i]) && (ls[i].activity() < extra_lim || satisfied(*this, ls[i])))
            removeClause(ls[i]);
    }
}


void Solver::compressDB(void)
{
    // UNFINISHED
    vec<Pair<int,Var> > n_occurs(2*nVars());

    for (int i = 0; i < n_occurs.size(); i++)
        n_occurs[i].fst = 0,
        n_occurs[i].snd = i;

    for (int i = 0; i < learnts.size(); i++){
        Clause  c = learnts[i];
        if (c == Clause_NULL) goto Skip;
        for (int j = 0; j < c.size(); j++)
            if (value(c[j]) == l_Error) goto Skip;

     for (int j = 0; j < c.size(); j++)
            n_occurs[index(c[j])].fst--;
        Skip:;
    }

    sort(n_occurs);

    /**/for (int i = 0; i < n_occurs.size(); i++) if (n_occurs[i].fst != 0) reportf("%d ", n_occurs[i].fst); reportf("\n"); exit(0);
}


/*_________________________________________________________________________________________________
|                                                                                                  
|  simplifyDB : [void]  ->  [bool]                                                                 
|                                                                                                  
|  Description:                                                                                    
|    Simplify all constraints according to the current top-level assigment (redundant constraints  
|    may be removed altogether).                                                                   
|________________________________________________________________________________________________@*/
void Solver::simplifyDB(bool subsume)
{
    if (!ok) return;    // GUARD (public method)
    assert(decisionLevel() == 0);

    // temporary placement -- put at end of solve?  <<= flytta till 'propagateToplevel()'
    // end

    if (!propagate().null()){   // (cannot use 'propagateToplevel()' here since it behaves different for 'occur_mode == occ_Off')
        assert(ok == false);
        return; }

    //**/if (occur_mode != occ_Off) clauseReduction();

    if (nAssigns() == last_simplify)
        return;
    last_simplify = nAssigns();

    // Subsumption simplification:
    //*HACK!*/if (!subsume && opt_var_elim){ opt_var_elim = false; if (opt_repeated_sub) reportf("                                                                                (var.elim. off)\r"); }
    if (occur_mode != occ_Off){
        if (!opt_repeated_sub){
            if (subsume) simplifyBySubsumption();
        }else
            simplifyBySubsumption();
    }

    // Removed satisfied clauses from the learnt clause database:
    if (stats.inspects - last_inspects > nLearnts() * 32){
        last_inspects = stats.inspects;
        for (int i = 0; i < learnts.size(); i++) if (!learnts[i].null() && satisfied(*this, learnts[i])) removeClause(learnts[i]);
    }

    //**/if (occur_mode != occ_Off) clauseReduction();
}


/*_________________________________________________________________________________________________
|                                                                                                  
|  search : (nof_conflicts : int) (nof_learnts : int) (params : const SearchParams&)  ->  [lbool]  
|                                                                                                  
|  Description:                                                                                    
|    Search for a model the specified number of conflicts, keeping the number of learnt clauses    
|    below the provided limit. NOTE! Use negative value for 'nof_conflicts' or 'nof_learnts' to    
|    indicate infinity.                                                                            
|                                                                                                  
|  Output:                                                                                         
|    'l_True' if a partial assigment that is consistent with respect to the clauseset is found. If 
|    all variables are decision variables, this means that the clause set is satisfiable. 'l_False'
|    if the clause set is unsatisfiable. 'l_Undef' if the bound on number of conflicts is reached. 
|________________________________________________________________________________________________@*/
lbool Solver::search(int nof_conflicts, int nof_learnts, const SearchParams& params)
{
    if (!ok) return l_False;    // GUARD (public method)
    assert(root_level == decisionLevel());

    stats.starts++;
    int     conflictC = 0;
    var_decay = 1 / params.var_decay;
    cla_decay = 1 / params.clause_decay;
    model.clear();

    for (;;){
        Clause confl = propagate();
        if (!confl.null()){
            // CONFLICT

            if (verbosity >= 2) reportf(L_IND"**CONFLICT**\n", L_ind);
            stats.conflicts++; conflictC++;
            vec<Lit>    learnt_clause;
            int         backtrack_level;
            if (decisionLevel() == root_level)
                return l_False;
            analyze(confl, learnt_clause, backtrack_level);
            cancelUntil(max(backtrack_level, root_level));
            record(learnt_clause);
            varDecayActivity(); claDecayActivity();

        }else{
            // NO CONFLICT

            if (nof_conflicts >= 0 && conflictC >= nof_conflicts){
                // Reached bound on number of conflicts:
                progress_estimate = progressEstimate();
                propQ.clear();
                cancelUntil(root_level);
                return l_Undef; }

            if (decisionLevel() == 0 && params.simplify){
                // Simplify the set of problem clauses:
                simplifyDB();
                if (!ok) return l_False;
            }

            if (nof_learnts >= 0 && nLearnts()-nAssigns() >= nof_learnts)
                // Reduce the set of learnt clauses:
                reduceDB();

            // New variable decision:
            stats.decisions++;
            Var next = order.select(params.random_var_freq);

            if (next == var_Undef){
                // Model found:
                model.growTo(nVars());
                if (occur_mode == occ_Off)
                    for (int i = 0; i < nVars(); i++) model[i] = value(i);
                else{
                  #ifdef WATCH_OPTIMIZATION
                    watches.clear(true); watches_setup = false;
                  #endif
                    extendModel();
                }
                cancelUntil(root_level);
                return l_True;
            }

            if (polarity_sug[next] == toInt(l_False))
                check(assume(~Lit(next)));
            else if (polarity_sug[next] == toInt(l_True))
                check(assume(Lit(next)));
            else
                check(assume(~Lit(next)));  // Arbitrarly default to negative polarity...
        }
    }
}


/*_________________________________________________________________________________________________
|                                                                                                  
|  extendModel : [void]  ->  [void]                                                                
|                                                                                                  
|  Description:                                                                                    
|    Extend the partial model of the current SAT environment to a full model, reading back         
|    eliminated clauses from the temporary file.                                                   
|________________________________________________________________________________________________@*/
void Solver::extendModel(void)
{
    if (verbosity >= 1){
        reportf("==============================================================================\n");
        reportf("(current CPU-time: %g s)\n", cpuTime()); }
    //    reportf("Extending model [cpu-time: %g s]\n", cpuTime()); }

    //**/int64 mem0 = memUsed();
  #if 1
    // READING FILE FORWARDS (simple)
    Solver  S(occ_Off);
    for (int i = 0; i < nVars(); i++){
        S.newVar();
        if (!var_elimed[i] && value(i) != l_Undef)
            S.addUnit(((value(i) == l_True) ? Lit(i) : ~Lit(i)));
    }

    fflush(elim_out);
    rewind(elim_out);
    for(;;){
        Lit     p;
        int     n = fread(&p, 4, 1, elim_out);
        if (n == 0)
            break;
        int     size = index(p);

        io_tmp.clear();
        io_tmp.growTo(size);
        fread((Lit*)io_tmp, 4, size, elim_out);

        S.addClause(io_tmp);
        /**/if (!S.ok) reportf("PANIC! False clause read back: "), dump(S, io_tmp);
        assert(S.ok);
    }
    fflush(elim_out);

    check(S.solve());
    S.model.moveTo(model);

  #else
    // READING FILE BACKWARDS
    Solver  S(occ_Off);
    S.setupWatches();
    for (int i = 0; i < nVars(); i++){
        S.newVar();
        if (!var_elimed[i] && value(i) != l_Undef)
            S.addUnit(((value(i) == l_True) ? Lit(i) : ~Lit(i)));
    }
    check(S.propagate() == Clause_NULL);

    fflush(elim_out);
    const int       chunk_size = 1000;
    vec<long>       offsets;
    vec<vec<Lit> >  tmps(chunk_size);
    Lit     p;
    int     n, size, c;

    rewind(elim_out);
    offsets.push(ftell(elim_out));
    c = 0;
    for(;;){
        n = fread(&p, 4, 1, elim_out); if (n == 0) break; size = index(p);              // (read size of clause or abort if no more clauses)
        io_tmp.clear(); io_tmp.growTo(size); fread((Lit*)io_tmp, 4, size, elim_out);    // (read clause)
        c++;
        if (c == chunk_size)
            offsets.push(ftell(elim_out)),
            c = 0;
    }

    assert(S.constrs.size() == 0);
    assert(S.propQ.size() == 0);
    rewind(elim_out);
    for (int i = offsets.size() - 1; i >= 0; i--){
        fseek(elim_out, SEEK_SET, offsets[i]);
        for(c = 0; c < chunk_size; c++){
            n = fread(&p, 4, 1, elim_out); if (n == 0) break; size = index(p);              // (read size of clause or abort if no more clauses)
            tmps[c].clear(); tmps[c].growTo(size); fread((Lit*)tmps[c], 4, size, elim_out); // (read clause)
        }
        for (; c > 0; c--){
            assert(S.watches_setup);
            S.addClause(tmps[c-1]); /*DEBUG*/if (!S.ok) reportf("PANIC! False clause read back: "), dump(S, io_tmp); assert(S.ok);
            check(S.propagate() == Clause_NULL);
        }
    }
    fseek(elim_out, SEEK_END, 0);
    fflush(elim_out);

    check(S.solve());
    S.model.moveTo(model);
  #endif

    //**/reportf("MEM USED for extending model: %g MB\n", (memUsed() - mem0) / (1024*1024.0));
}


// Return search-space coverage. Not extremely reliable.
//
double Solver::progressEstimate(void)
{
    int n_vars = 0; for (int i = 0; i < nVars(); i++) n_vars += !var_elimed[i];
    double  progress = 0;
    double  F = 1.0 / n_vars;
    for (int i = 0; i < nVars(); i++)
        if (value(i) != l_Undef && !var_elimed[i])
            progress += pow(F, level[i]);
    return progress / n_vars;
}


// Divide all variable activities by 1e100.
//
void Solver::varRescaleActivity(void)
{
    for (int i = 0; i < nVars(); i++)
        activity[i] *= 1e-100;
    var_inc *= 1e-100;
}


// Divide all constraint activities by 1e100.
//
void Solver::claRescaleActivity(void)
{
    for (int i = 0; i < learnts.size(); i++)
        if (!learnts[i].null())
            learnts[i].activity() *= 1e-20;
    cla_inc *= 1e-20;
}


/*_________________________________________________________________________________________________
|                                                                                                  
|  solve : (assumps : const vec<Lit>&)  ->  [bool]                                                 
|                                                                                                  
|  Description:                                                                                    
|    Top-level solve. If using assumptions (non-empty 'assumps' vector), you must call             
|    'simplifyDB()' first to see that no top-level conflict is present (which would put the solver 
|    in an undefined state).                                                                       
|________________________________________________________________________________________________@*/
bool Solver::solve(const vec<Lit>& assumps)
{
    if (verbosity >= 1){
        reportf("=================================[SATELITE+]==================================\n");
        reportf("|           |     ORIGINAL     |              LEARNT              |          |\n");
        reportf("| Conflicts | Clauses Literals |   Limit Clauses Literals  Lit/Cl | Progress |\n");
        reportf("==============================================================================\n");
        // Hack:
        double nof_learnts = 0.3;
        int learnt_literals = 0; for (int i = 0; i < learnts.size(); i++) if (!learnts[i].null()) learnt_literals += learnts[i].size();
        reportf("| %9d | %7d %8d | %7d %7d %8d %7.1f | %6.3f %% |\n", (int)stats.conflicts, nClauses(), nLiterals(), (int)(nof_learnts*nClauses()), nLearnts(), learnt_literals, learnt_literals / (double)nLearnts(), progress_estimate*100);
    }

    if (occur_mode == occ_Off || opt_asym_branch) setupWatches();
    simplifyDB(true);
    if (!ok) return false;
    setupWatches();

    SearchParams    params(0.95, 0.999, opt_no_random ? 0 : 0.02);
    double  nof_conflicts = 100;
    double  nof_learnts   = 0.4; //0.3;
    lbool   status        = l_Undef;

    for (int i = 0; i < assumps.size(); i++){
        assert(!var_elimed[var(assumps[i])]);
        if (!propagate().null() || !assume(assumps[i]) || !propagate().null()){
            propQ.clear();
            cancelUntil(0);
            return false; }
    }
    root_level = decisionLevel();

    while (status == l_Undef){
        //if (verbosity >= 1) reportf("RESTART -- conflicts=%d   clauses=%d   learnts=%d/%d   progress=%.4f %%\n", (int)nof_conflicts, constrs.size(), learnts.size(), (int)(nof_learnts*nClauses()), progress_estimate*100);
        if (verbosity >= 1){
            int learnt_literals = 0; for (int i = 0; i < learnts.size(); i++) if (!learnts[i].null()) learnt_literals += learnts[i].size();
            reportf("| %9d | %7d %8d | %7d %7d %8d %7.1f | %6.3f %% |\n", (int)stats.conflicts, nClauses(), nLiterals(), (int)(nof_learnts*nClauses()), nLearnts(), learnt_literals, learnt_literals / (double)nLearnts(), progress_estimate*100);
        }
        status = search((int)nof_conflicts, (int)(nof_learnts*nClauses()), params);
        nof_conflicts *= 1.5;
        nof_learnts   *= 1.1;
        //**/if (learnts.size() > 1000) compressDB();
    }

    if (verbosity >= 1) reportf("==============================================================================\n");

    cancelUntil(0);
    return (status == l_True);
}


//#################################################################################################
// "Solver_subsume.C" BEGIN
//#################################################################################################

#define Report(format, args...) ((verbosity >= 1) ? reportf(format , ## args) : (void)0)
//#define Report(format, args...)
#define Report2(format, args...) ((verbosity >= 2) ? reportf(format , ## args) : (void)0)


macro bool has(Clause c, Lit p) {
    for (int i = 0; i < c.size(); i++)
        if (c[i] == p) return true;
    return false; }


#if 1
// Assumes 'seen' is cleared (will leave it cleared)
static
bool subset(Clause A, Clause B, vec<char>& seen)
{
    for (int i = 0; i < B.size(); i++) seen[index(B[i])] = 1;
    for (int i = 0; i < A.size(); i++){
        if (!seen[index(A[i])]){
            for (int i = 0; i < B.size(); i++) seen[index(B[i])] = 0;
            return false;
        }
    }
    for (int i = 0; i < B.size(); i++) seen[index(B[i])] = 0;
    return true;
}
#else
static
bool subset(Clause A, Clause B, vec<char>& seen)
{
    for (int i = 0; i < A.size(); i++)
        if (!has(B, A[i]))
            return false;
    return true;
}
#endif


macro bool subset(uint64 A, uint64 B) { return (A & ~B) == 0; }


// Assumes 'seen' is cleared (will leave it cleared)
static
bool selfSubset(Clause A, Clause B, vec<char>& seen)
{
    for (int i = 0; i < B.size(); i++) seen[index(B[i])] = 1;

    bool    flip = false;
    for (int i = 0; i < A.size(); i++){
        if (!seen[index(A[i])]){
            if (flip == true || !seen[index(~A[i])]){
                for (int i = 0; i < B.size(); i++) seen[index(B[i])] = 0;
                return false;
            }
            flip = true;
        }
    }
    for (int i = 0; i < B.size(); i++) seen[index(B[i])] = 0;
    return flip;
}


macro bool selfSubset(uint64 A, uint64 B)
{
    uint64 B_tmp = B | ((B & 0xAAAAAAAAAAAAAAAALL) >> 1) | ((B & 0x5555555555555555LL) << 1);
    if ((A & ~B_tmp) == 0){
        uint64 C = A & ~B;
        return (C & (C-1)) == 0;
    }else
        return false;
}


void Solver::findSubsumed(Clause ps, vec<Clause>& out_subsumed)
{
    uint64  abst;
    if (ps.dynamic()){
        abst = 0;
        for (int i = 0; i < ps.size(); i++)
            abst |= abstLit(ps[i]);
    }else
        abst = ps.abst();

    int min_i = 0;
    for (int i = 1; i < ps.size(); i++){
        if (occur[index(ps[i])].size() < occur[index(ps[min_i])].size())
            min_i = i;
    }

    vec<Clause>& cs = occur[index(ps[min_i])];
    for (int i = 0; i < cs.size(); i++){
        if (cs[i] != ps && ps.size() <= cs[i].size() && subset(abst, cs[i].abst()) && subset(ps, cs[i], seen_tmp))
            out_subsumed.push(cs[i]);
    }
}


bool Solver::isSubsumed(Clause ps)
{
    uint64  abst;
    if (ps.dynamic()){
        abst = 0;
        for (int i = 0; i < ps.size(); i++)
            abst |= abstLit(ps[i]);
    }else
        abst = ps.abst();

    for (int j = 0; j < ps.size(); j++){
        vec<Clause>& cs = occur[index(ps[j])];

        for (int i = 0; i < cs.size(); i++){
            if (cs[i] != ps && cs[i].size() <= ps.size() && subset(cs[i].abst(), abst)){
                if (subset(cs[i], ps, seen_tmp))
                    return true;
            }
        }
    }
    return false;
}


// Will put NULL in 'cs' if clause removed.
void Solver::subsume0(Clause ps, int& counter)
{
    if (!opt_0sub) return;

    vec<Clause>    subs;
    findSubsumed(ps, subs);
    for (int i = 0; i < subs.size(); i++){
        if (&counter != NULL) counter++;
        removeClause(subs[i]);
    }
}


#if 1
// With queue
void Solver::subsume1(Clause ps, int& counter)
{
    if (!opt_1sub) return;

    vec<Clause>    Q;
    vec<Clause>    subs;
    Clause_t       qs;
    int            q;

    registerIteration(Q);
    registerIteration(subs);

    Q.push(ps);
    q = 0;
    while (q < Q.size()){
        if (Q[q].null()) { q++; continue; }
        qs.clear();
        for (int i = 0; i < Q[q].size(); i++)
            qs.push(Q[q][i]);

        for (int i = 0; i < qs.size(); i++){
            qs[i] = ~qs[i];
            findSubsumed(qs, subs);
            for (int j = 0; j < subs.size(); j++){
                /*DEBUG*/
              #ifndef NDEBUG
                if (&counter != NULL && counter == -1){
                    dump(subs[j]);
                    qs[i] = ~qs[i];
                    dump(qs);
                    printf(L_LIT"\n", L_lit(qs[i]));
                    exit(0);
                }
              #endif
                /*END*/
                if (subs[j].null()) continue;
                if (&counter != NULL) counter++;
                strengthenClause(subs[j], qs[i]);
                Q.push(subs[j]);
            }

            qs[i] = ~qs[i];
            subs.clear();
        }
        q++;
    }

    unregisterIteration(Q);
    unregisterIteration(subs);
}
#else
// Without queue
void Solver::subsume1(Clause ps, int& counter)
{
    if (!opt_1sub) return;

    vec<Clause>    subs;
    Clause_t       qs;

    registerIteration(subs);

    assert(!ps.null());
    for (int i = 0; i < ps.size(); i++)
        qs.push(ps[i]);

    for (int i = 0; i < qs.size(); i++){
        qs[i] = ~qs[i];
        findSubsumed(qs, subs);
        for (int j = 0; j < subs.size(); j++){
            if (subs[j].null()) continue;
            if (&counter != NULL) counter++;
            strengthenClause(subs[j], qs[i]);
        }

        qs[i] = ~qs[i];
        subs.clear();
    }
    unregisterIteration(subs);
}
#endif


void Solver::simplifyBySubsumption(bool with_var_elim)
{
    propagateToplevel(); if (!ok){ Report("(contradiction during subsumption)\n"); return; }

    int     orig_n_clauses  = nClauses();
    int     orig_n_literals = nLiterals();
    do{
        // SUBSUMPTION:
        //
      #ifndef SAT_LIVE
        Report("  -- subsuming                       \r");
      #endif
        int     clauses_subsumed = 0, literals_removed = 0;
        if (!opt_0sub && !opt_1sub){
            cl_added  .clear();
            cl_touched.clear();
            goto NoSubsumption; }

        if (cl_added.size() > nClauses() / 2){
            // Optimized variant when virtually whole database is involved:
            cl_added  .clear();
            cl_touched.clear();

            for (int i = 0; i < constrs.size(); i++) if (!constrs[i].null()) subsume1(constrs[i], literals_removed);
            propagateToplevel(); if (!ok){ Report("(contradiction during subsumption)\n"); return; }

            CSet s1;
            registerIteration(s1);
            while (cl_touched.size() > 0){
                for (int i = 0; i < cl_touched.size(); i++)
                    if (!cl_touched[i].null())
                        s1.add(cl_touched[i]);
                cl_touched.clear();

                for (int i = 0; i < s1.size(); i++) if (!s1[i].null()) subsume1(s1[i], literals_removed);
                s1.clear();

                propagateToplevel();
                if (!ok){ Report("(contradiction during subsumption)\n"); unregisterIteration(s1); return; }
            }
            unregisterIteration(s1);

            for (int i = 0; i < constrs.size(); i++) if (!constrs[i].null()) subsume0(constrs[i], clauses_subsumed);

        }else{
            //  Set used in 1-subs:
            //      (1) clauses containing a negated literal of an added clause.
            //      (2) all added or strengthened ("touched") clauses.
            //
            //  Set used in 0-subs:
            //      (1) clauses containing a (non-negated) literal of an added clause, including the added clause itself.
            //      (2) all strenghtened clauses -- REMOVED!! We turned on eager backward subsumption which supersedes this.

//Report("  PREPARING\n");
            CSet        s0, s1;     // 's0' is used for 0-subsumption, 's1' for 1-subsumption
            vec<char>   ol_seen(nVars()*2, 0);
            for (int i = 0; i < cl_added.size(); i++){
                Clause  c = cl_added[i]; if (c.null()) continue;
                s1.add(c);
                for (int j = 0; j < c.size(); j++){
                    if (ol_seen[index(c[j])]) continue;
                    ol_seen[index(c[j])] = 1;

                    vec<Clause>& n_occs = occur[index(~c[j])];
                    for (int k = 0; k < n_occs.size(); k++)     // <<= B�ttra p�. Beh�ver bara kolla 'n_occs[k]' mot 'c'
                        if (n_occs[k] != c && n_occs[k].size() <= c.size() && selfSubset(n_occs[k].abst(), c.abst()) && selfSubset(n_occs[k], c, seen_tmp))
                            s1.add(n_occs[k]);

                    vec<Clause>& p_occs = occur[index(c[j])];
                    for (int k = 0; k < p_occs.size(); k++)     // <<= B�ttra p�. Beh�ver bara kolla 'p_occs[k]' mot 'c'
                        if (subset(p_occs[k].abst(), c.abst()))
                            s0.add(p_occs[k]);
                }
            }
            cl_added.clear();

            registerIteration(s0);
            registerIteration(s1);

//Report("  FIXED-POINT\n");
            // Fixed-point for 1-subsumption:
            while (s1.size() > 0 || cl_touched.size() > 0){
                for (int i = 0; i < cl_touched.size(); i++)
                    if (!cl_touched[i].null())
                        s1.add(cl_touched[i]),
                        s0.add(cl_touched[i]);

                cl_touched.clear();
                assert(propQ.size() == 0);
//Report("s1.size()=%d  cl_touched.size()=%d\n", s1.size(), cl_touched.size());
                for (int i = 0; i < s1.size(); i++) if (!s1[i].null()){ subsume1(s1[i], literals_removed); }
                s1.clear();

                propagateToplevel();
                if (!ok){
                    Report("(contradiction during subsumption)\n");
                    unregisterIteration(s0);
                    unregisterIteration(s1);
                    return; }
                assert(cl_added.size() == 0);
            }
            unregisterIteration(s1);

            // Iteration pass for 0-subsumption:
            for (int i = 0; i < s0.size(); i++) if (!s0[i].null()) subsume0(s0[i], clauses_subsumed);
            s0.clear();
            unregisterIteration(s0);
        }

        if (literals_removed > 0 || clauses_subsumed > 0)
            Report2("  #literals-removed: %d    #clauses-subsumed: %d\n", literals_removed, clauses_subsumed);


        // VARIABLE ELIMINATION:
        //
      NoSubsumption:
        if (!with_var_elim || !opt_var_elim) break;

//Report("VARIABLE ELIMINIATION\n");
        vec<Var>    init_order;
        orderVarsForElim(init_order);   // (will untouch all variables)

        for (bool first = true;; first = false){
            int vars_elimed = 0;
            int clauses_before = nClauses();
            vec<Var>    order;

            if (opt_pure_literal){
                for (int i = 0; i < touched_list.size(); i++){
                    Var x = touched_list[i];
                    if (n_occurs[index(Lit(x))] == 0 && value(x) == l_Undef && !frozen[x] && !var_elimed[x] && n_occurs[index(~Lit(x))] > 0){
                        enqueue(~Lit(x));
                    }else if (n_occurs[index(~Lit(x))] == 0 && value(x) == l_Undef && !frozen[x] && !var_elimed[x] && n_occurs[index(Lit(x))] > 0){
                        enqueue(Lit(x));
                    }
                }
                propagateToplevel(); assert(ok);
            }

            if (first)
                init_order.copyTo(order);
            else{
                for (int i = 0; i < touched_list.size(); i++)
                    if (!var_elimed[touched_list[i]])
                        order.push(touched_list[i]),
                    touched[touched_list[i]] = 0;
                touched_list.clear();
            }

            assert(propQ.size() == 0);
            for (int i = 0; i < order.size(); i++){
              #ifndef SAT_LIVE
                if (i % 1000 == 999 || i == order.size()-1) Report("  -- var.elim.:  %d/%d          \r", i+1, order.size());
              #endif
                if (maybeEliminate(order[i])){
                    vars_elimed++;
                    propagateToplevel(); if (!ok){ Report("(contradiction during subsumption)\n"); return; }
                }
            }
            assert(propQ.size() == 0);

            if (vars_elimed == 0)
                break;

            Report2("  #clauses-removed: %-8d #var-elim: %d\n", clauses_before - nClauses(), vars_elimed);
        }
        //assert(touched_list.size() == 0);     // <<= No longer true, due to asymmetric branching. Ignore it for the moment...
//  }while (cl_added.size() > 0);
    }while (cl_added.size() > 100);

    if (orig_n_clauses != nClauses() || orig_n_literals != nLiterals())
        //Report("#clauses: %d -> %d    #literals: %d -> %d    (#learnts: %d)\n", orig_n_clauses, nClauses(), orig_n_literals, nLiterals(), nLearnts());
        //Report("#clauses:%8d (%6d removed)    #literals:%8d (%6d removed)    (#learnts:%8d)\n", nClauses(), orig_n_clauses-nClauses(), nLiterals(), orig_n_literals-nLiterals(), nLearnts());
        Report("| %9d | %7d %8d | %7s %7d %8s %7s | %6s   | %d/%d\n", (int)stats.conflicts, nClauses(), nLiterals(), "--", nLearnts(), "--", "--", "--", nClauses() - orig_n_clauses, nLiterals() - orig_n_literals);


    if (opt_pre_sat){
        // Compact variables:
        vec<Var>    vmap(nVars(), -1);
        int         c = 0;
        for (int i = 0; i < nVars(); i++)
            if (!var_elimed[i] && (occur[index(Lit(i))].size() != 0 || occur[index(~Lit(i))].size() != 0))
                vmap[i] = c++;

        Report("==============================================================================\n");
        Report("Result  :   #vars: %d   #clauses: %d   #literals: %d\n", c, nClauses(), nLiterals());
        Report("CPU time:   %g s\n", cpuTime());
        Report("==============================================================================\n");

        // Write CNF or BCNF file:
        cchar*  filename = (output_file == NULL) ? "pre-satelited.cnf" : output_file;
        int     len = strlen(filename);
        if (len >= 5 && strcmp(filename+len-5, ".bcnf") == 0){
            BcnfWriter w(filename);
            vec<Lit>   lits;
            for (int i = 0; i < constrs.size(); i++){
                Clause c = constrs[i]; if (c.null()) continue;
                lits.clear(); lits.growTo(c.size());
                for (int j = 0; j < c.size(); j++)
                    lits[j] = Lit(vmap[var(c[j])], sign(c[j]));
                w.addClause(lits);
            }
        }else{
            FILE*   out = (strcmp(filename, "-") == 0) ? stdout : fopen(filename, "wb");
            if (out == NULL) fprintf(stderr, "ERROR! Could not open output file: %s\n", filename), exit(1);

            fprintf(out, "p cnf %d %d\n", c, nClauses());
            fprintf(out, "c\n");
            fprintf(out, "c   #vars: %d   #clauses: %d   #literals: %d\n"  , c, nClauses(), nLiterals());
            fprintf(out, "c\n");
            for (int i = 0; i < constrs.size(); i++){
                Clause c = constrs[i]; if (c.null()) continue;
                for (int j = 0; j < c.size(); j++)
                    assert(vmap[var(c[j])]+1 != -1),
                    fprintf(out, "%s%d ", sign(c[j])?"-":"", vmap[var(c[j])]+1);
                fprintf(out, "0\n");
            }
            fclose(out);
        }

        // Write varible map:
        if (varmap_file != NULL){
            FILE*   out = fopen(varmap_file, "wb");
            if (out == NULL) fprintf(stderr, "ERROR! Could not open output file: %s\n", varmap_file), exit(1);

            fprintf(out, "%d 0\n", nVars());

            c = 0;
            for (int i = 0; i < nVars(); i++){
                if (!var_elimed[i] && value(i) != l_Undef){
                    c++;
                    fprintf(out, "%s%s%d", (c > 1)?" ":"", (value(i)==l_True)?"":"-", i+1);
                }
            }
            fprintf(out, " 0\n");

            c = 0;
            for (int i = 0; i < nVars(); i++){
                if (vmap[i] != -1){
                    assert(vmap[i] == c);
                    c++;
                    fprintf(out, "%s%d", (c > 1)?" ":"", i+1);
                }
            }
            fprintf(out, " 0\n");
            fclose(out);
        }
        fflush(elim_out);
        exit(0);
    }

//Report("DONE!\n");
}





//=================================================================================================
// Eliminate variables:


// Side-effect: Will untouch all variables.
void Solver::orderVarsForElim(vec<Var>& order)
{
    order.clear();
    vec<Pair<int, Var> > cost_var;
    for (int i = 0; i < touched_list.size(); i++){
        Var x = touched_list[i];
        touched[x] = 0;
        cost_var.push(Pair_new( occur[index(Lit(x))].size() * occur[index(~Lit(x))].size() , x ));
    }
    touched_list.clear();
    sort(cost_var);

    for (int x = 0; x < cost_var.size(); x++)
        if (cost_var[x].fst != 0)
            order.push(cost_var[x].snd);
}


#if 0
// Returns FALSE if clause is always satisfied ('out_clause' should not be used).
static
bool merge(Clause ps, Clause qs, Lit without_p, Lit without_q, vec<Lit>& out_clause)
{
    int     i = 0, j = 0;
    while (i < ps.size() && j < qs.size()){
        if      (ps[i] == without_p) i++;
        else if (qs[j] == without_q) j++;
        else if (ps[i] == ~qs[j])    return false;
        else if (ps[i] < qs[j])      out_clause.push(ps[i++]);
        else if (ps[i] > qs[j])      out_clause.push(qs[j++]);
        else                         out_clause.push(ps[i++]), j++;
    }
    while (i < ps.size()){
        if (ps[i] == without_p) i++;
        else                    out_clause.push(ps[i++]); }
    while (j < qs.size()){
        if (qs[j] == without_q) j++;
        else                    out_clause.push(qs[j++]); }
    return true;
}
#endif

// Returns FALSE if clause is always satisfied ('out_clause' should not be used). 'seen' is assumed to be cleared.
static
bool merge(Clause ps, Clause qs, Lit without_p, Lit without_q, vec<char>& seen, vec<Lit>& out_clause)
{
    for (int i = 0; i < ps.size(); i++){
        if (ps[i] != without_p){
            seen[index(ps[i])] = 1;
            out_clause.push(ps[i]);
        }
    }

    for (int i = 0; i < qs.size(); i++){
        if (qs[i] != without_q){
            if (seen[index(~qs[i])]){
                for (int i = 0; i < ps.size(); i++) seen[index(ps[i])] = 0;
                return false; }
            if (!seen[index(qs[i])])
                out_clause.push(qs[i]);
        }
    }

    for (int i = 0; i < ps.size(); i++) seen[index(ps[i])] = 0;
    return true;
}


Lit Solver::findUnitDef(Var x, vec<Clause>& poss, vec<Clause>& negs)
{
    vec<Lit>    imps;   // All p:s s.t. "~x -> p"
    for (int i = 0; i < poss.size(); i++){
        if (poss[i].size() == 2){
            if (var(poss[i][0]) == x)
                imps.push(poss[i][1]);
            else
                assert(var(poss[i][1]) == x),
                imps.push(poss[i][0]);
        }
    }

    // Quadratic algorithm; maybe should improve?
    for (int i = 0; i < negs.size(); i++){
        if (negs[i].size() == 2){
            Lit p = (var(negs[i][0]) == x ) ? ~negs[i][1] : (assert(var(negs[i][1]) == x), ~negs[i][0]);
            for (int j = 0; j < imps.size(); j++)
                if (imps[j] == p)
                    return ~imps[j];
                //**/else if (imps[j] == ~p)
                //**/    enqueue(~p);
        }
    }

    return lit_Undef;
}


// If returns TRUE 'out_def' is the definition of 'x' (always a disjunction).
// An empty definition means '~x' is globally true.
//
bool Solver::findDef(Lit x, vec<Clause>& poss, vec<Clause>& negs, Clause out_def)
{
    assert(out_def.size() == 0);
    Clause_t    imps;       // (negated implied literals, actually)
    uint64      abst = 0;

    for (int i = 0; i < negs.size(); i++){
        if (negs[i].size() == 2){
            Lit imp;
            if (negs[i][0] == ~x)
                imp = ~negs[i][1];
            else
                assert(negs[i][1] == ~x),
                imp = ~negs[i][0];

            imps.push(imp);
            abst |= abstLit(imp);
        }
    }

    if (opt_hyper1_res && isSubsumed(imps))
        return true;
    if (!opt_def_elim)
        return false;

    imps.push(x);
    abst |= abstLit(x);

    for (int i = 0; i < poss.size(); i++){
        if (subset(poss[i].abst(), abst) && subset(poss[i], imps, seen_tmp)){
            Clause c = poss[i];

            if (out_def.size() > 0 && c.size()-1 < out_def.size())
                //**/printf("----------\nImps  : "), dump(imps),printf("Before: "), dump(out_def), printf("Now   : "), dump(c, false), printf("  -- minus "L_LIT"\n", L_lit(x));
                out_def.clear();    // (found a better definition)
            if (out_def.size() == 0){
                for (int j = 0; j < c.size(); j++)
                    if (c[j] != x)
                        out_def.push(~c[j]);
                assert(out_def.size() == c.size()-1);
            }
        }
    }
    return (out_def.size() > 0);
}

/*
 - 1 2 3                 -3 -4
  1 -2 3  resolved with  -3 -5
  4  5 3                 -3 6 7
                         -3 -6 -7

3 -> ~4       == { ~3, ~4 }   (ger ~4, ~5 + ev mer som superset; negera detta och l�gg till 3:an)
3 -> ~5       == { ~3, ~5 }

~4 & ~5 -> 3  == { 4, 5, 3 }
*/

int Solver::substitute(Lit x, Clause def, vec<Clause>& poss, vec<Clause>& negs, vec<Clause>& new_clauses = *(vec<Clause>*)NULL)
{
    vec<Lit>    tmp;
    int         counter = 0;


    // Positives:
    //**/if (hej) printf("    --from POS--\n");
    for (int i = 0; i < def.size(); i++){
        for (int j = 0; j < poss.size(); j++){
            if (poss[j].null()) continue;

            Clause c = poss[j];
            for (int k = 0; k < c.size(); k++){
                if (c[k] == ~def[i])
                    goto Skip;
            }
            tmp.clear();
            for (int k = 0; k < c.size(); k++){
                if (c[k] == x)
                    tmp.push(def[i]);
                else
                    tmp.push(c[k]);
            }
            //**/if (hej) printf("    "), dump(*this, tmp);
            if (&new_clauses != NULL){
                Clause tmp_c;
                tmp_c = addClause(tmp);
                if (!tmp_c.null())
                    new_clauses.push(tmp_c);
            }else{
                sortUnique(tmp);
                for (int i = 0; i < tmp.size()-1; i++)
                    if (tmp[i] == ~tmp[i+1])
                        goto Skip;
                counter++;
            }
          Skip:;
        }
    }
    //**/if (hej) printf("    --from NEG--\n");

    // Negatives:
    for (int j = 0; j < negs.size(); j++){
        if (negs[j].null()) continue;

        Clause c = negs[j];
        // If any literal from 'def' occurs in 'negs[j]', it is satisfied.
        for (int i = 0; i < def.size(); i++){
            if (has(c, def[i]))
                goto Skip2;
        }

        tmp.clear();
        for (int i = 0; i < c.size(); i++){
            if (c[i] == ~x){
                for (int k = 0; k < def.size(); k++)
                    tmp.push(~def[k]);
            }else
                tmp.push(c[i]);
        }
        //**/if (hej) printf("    "), dump(*this, tmp);
        if (&new_clauses != NULL){
            Clause tmp_c;
            tmp_c = addClause(tmp);
            if (!tmp_c.null())
                new_clauses.push(tmp_c);
        }else{
            sortUnique(tmp);
            for (int i = 0; i < tmp.size()-1; i++)
                if (tmp[i] == ~tmp[i+1])
                    goto Skip2;
            counter++;
        }
      Skip2:;
    }

    //**/if (counter != 0 && def.size() >= 2) exit(0);
    return counter;
}


void Solver::asymmetricBranching(Lit p)
{
    if (!opt_asym_branch) return;

    assert(decisionLevel() == 0);
    propagateToplevel(); if (!ok){ Report("(contradiction during asymmetric branching)\n"); return; }

    vec<Lit>    learnt;
    //vec<char>&  seen = seen_tmp;
    vec<Clause> cs;
    occur[index(p)].copyTo(cs);
    registerIteration(cs);

    for (int i = 0; i < cs.size(); i++){
        Clause  c = cs[i]; if (c.null()) continue;

#if 1
        Clause  confl;
        bool    analyze    = false;

        //**/printf("Inspecting: "); dump(*this, c);
        for (int j = 0; j < c.size(); j++){
            if (c[j] == p) continue;
            Lit q = ~c[j];
            //**/printf("Assuming: "L_LIT"\n", L_lit(q));
            if (!assume(q)){
                //**/printf("Assumption failed!\n");
                learnt.clear();
                if (level[var(q)] > 0)
                    learnt.push(~q);
                Clause rs = reason[var(q)];
                for (int k = 1; k < rs.size(); k++){
                    assert(value(rs[k]) == l_False);
                    seen_tmp[index(~rs[k])] = 1; }
                analyze = true;
                break;
            }

            confl = propagate();
            if (!confl.null()){
                //**/printf("Lead to conflict! "); dump(*this, confl);
                learnt.clear();
                for (int k = 0; k < confl.size(); k++){
                    assert(value(confl[k]) == l_False);
                    seen_tmp[index(~confl[k])] = 1; }
                analyze = true;
                break;
            }
        }

        if (analyze){
            // Analyze conflict:
            for (int j = trail.size()-1; j >= 0; j--){
                Lit q = trail[j];
                if (seen_tmp[index(q)]){
                    Clause rs = reason[var(q)];
                    if (rs.null()){
                        if (level[var(q)] > 0)
                            learnt.push(~q);
                    }else{
                        for (int k = 1; k < rs.size(); k++)
                            seen_tmp[index(~rs[k])] = 1;
                    }
                    seen_tmp[index(q)] = 0;
                }
            }
        }

        cancelUntil(0);

        /*Test 'seen[]'*/
        //for (int i = 0; i < seen.size(); i++) assert(seen[i] == 0);
        /*End*/

        if (analyze){
            // Update clause:
            //**/printf("%d ", c.size() - learnt.size()); fflush(stdout);
            //**/{static int counter = 0; static const char chr[] = { '|', '/', '-', '\\' }; printf("%c\r", chr[counter++ & 3]); fflush(stdout); }
            assert(propQ.size() == 0);
            /*Test Subset*/
          #ifndef NDEBUG
            for (int i = 0; i < learnt.size(); i++){
                for (int j = 0; j < c.size(); j++)
                    if (c[j] == learnt[i])
                        goto Found;
                printf("\n");
                printf("asymmetricBranching("L_LIT" @ %d)\n", L_lit(p), level[var(p)]);
                printf("learnt: "); dump(*this, learnt);
                printf("c     : "); dump(*this, c     );
                assert(false);  // no subset!
              Found:;
            }
          #endif
            /*End*/
            unlinkClause(c);    // (will touch all variables)
            addClause(learnt, false, c);
            propagateToplevel(); if (!ok){ Report("(contradiction during asymmetric branching after learning unit fact)\n"); return; }
        }

#else
        assert(propQ.size() == 0);
        bool conflict = false;
        for (int j = 0; j < c.size(); j++){
            if (c[j] == p) continue;
            if (!assume(~c[j]))      { conflict = true; break; }
            if (!propagate().null()) { conflict = true; break; }
        }
        assert(propQ.size() == 0);
        cancelUntil(0);

        if (conflict){
            /**/putchar('+'); fflush(stdout);
            //**/printf("----------------------------------------\n");
            //**/printf("asymmetricBranching("L_LIT" @ %d)\n", L_lit(p), level[var(p)]);
            //**/dump(*this, c);
            //**/for (int j = 0; j < c.size(); j++) if (c[j] == p) goto Found; assert(false); Found:;
            assert(!cs[i].null());
            strengthenClause(c, p);
            propagateToplevel(); if (!ok){ Report("(contradiction during asymmetric branching after learning unit fact)\n"); return; }
        }
#endif
    }

    unregisterIteration(cs);
    assert(propQ.size() == 0);
}


// Returns TRUE if variable was eliminated.
bool Solver::maybeEliminate(Var x)
{
    assert(propQ.size() == 0);
    assert(!var_elimed[x]);
    if (frozen[x]) return false;
    if (value(x) != l_Undef) return false;
    if (occur[index(Lit(x))].size() == 0 && occur[index(~Lit(x))].size() == 0) return false;

    vec<Clause>&   poss = occur[index( Lit(x))];
    vec<Clause>&   negs = occur[index(~Lit(x))];
    vec<Clause>    new_clauses;

    int before_clauses  = -1;
    int before_literals = -1;

    bool    elimed     = false;
    bool    tried_elim = false;

    #define MigrateToPsNs vec<Clause> ps; poss.moveTo(ps); vec<Clause> ns; negs.moveTo(ns); for (int i = 0; i < ps.size(); i++) unlinkClause(ps[i], x); for (int i = 0; i < ns.size(); i++) unlinkClause(ns[i], x);
    #define DeallocPsNs   for (int i = 0; i < ps.size(); i++) deallocClause(ps[i]); for (int i = 0; i < ns.size(); i++) deallocClause(ns[i]);

    // Find 'x <-> p':
    if (opt_unit_def){
        Lit p = findUnitDef(x, poss, negs);
        /**/if (p == lit_Undef) propagateToplevel(); if (!ok) return true;
        if (p != lit_Undef){
            //**/printf("DEF: x%d = "L_LIT"\n", x, L_lit(p));
            /*BEG
            printf("POS:\n"); for (int i = 0; i < poss.size(); i++) printf("  "), dump(poss[i]);
            printf("NEG:\n"); for (int i = 0; i < negs.size(); i++) printf("  "), dump(negs[i]);
            printf("DEF: x%d = "L_LIT"\n", x, L_lit(p));
            hej = true;
            END*/
            Clause_t    def; def.push(p);
            MigrateToPsNs
            substitute(Lit(x), def, ps, ns, new_clauses);
            /*BEG
            hej = false;
            printf("NEW:\n");
            for (int i = 0; i < new_clauses.size(); i++)
                printf("  "), dump(new_clauses[i]);
            END*/
            propagateToplevel(); if (!ok) return true;
            DeallocPsNs
            goto Eliminated;
        }
    }

    // ...
#if 0
    if (poss.size() < 10){ asymmetricBranching( Lit(x)); if (!ok) return true; }
    if (negs.size() < 10){ asymmetricBranching(~Lit(x)); if (!ok) return true; }
    if (value(x) != l_Undef) return false;
#endif

    // Heuristic:
    if (poss.size() >= 8 && negs.size() >= 8)      // <<== CUT OFF
//  if (poss.size() >= 7 && negs.size() >= 7)      // <<== CUT OFF
//  if (poss.size() >= 6 && negs.size() >= 6)      // <<== CUT OFF
        return false;

    // Count clauses/literals before elimination:
    before_clauses  = poss.size() + negs.size();
    before_literals = 0;
    for (int i = 0; i < poss.size(); i++) before_literals += poss[i].size();
    for (int i = 0; i < negs.size(); i++) before_literals += negs[i].size();

    if (poss.size() >= 3 && negs.size() >= 3 && before_literals > 300)  // <<== CUT OFF
        return false;

    // Check for definitions:
    if (opt_def_elim || opt_hyper1_res){
        //if (poss.size() > 1 || negs.size() > 1){
        if (poss.size() > 2 || negs.size() > 2){
            Clause_t def;
            int      result_size;
            if (findDef(Lit(x), poss, negs, def)){
                if (def.size() == 0){ enqueue(~Lit(x)); return true; }  // Hyper-1-resolution
                result_size = substitute(Lit(x), def, poss, negs);
                if (result_size <= poss.size() + negs.size()){  // <<= elimination threshold (maybe subst. should return literal count as well)
                    MigrateToPsNs
                    substitute(Lit(x), def, ps, ns, new_clauses);
                    propagateToplevel(); if (!ok) return true;
                    DeallocPsNs
                    goto Eliminated;
                }else
                    tried_elim = true,
                    def.clear();
            }
            if (!elimed && findDef(~Lit(x), negs, poss, def)){
                if (def.size() == 0){ enqueue(Lit(x)); return true; }  // Hyper-1-resolution
                result_size = substitute(~Lit(x), def, negs, poss);
                if (result_size <= poss.size() + negs.size()){  // <<= elimination threshold
                    MigrateToPsNs
                    substitute(~Lit(x), def, ns, ps, new_clauses);
                    propagateToplevel(); if (!ok) return true;
                    DeallocPsNs
                    goto Eliminated;
                }else
                    tried_elim = true;
            }
        }
    }

    if (!tried_elim){
        // Count clauses/literals after elimination:
        int after_clauses  = 0;
        int after_literals = 0;
        Clause_t  dummy;
        for (int i = 0; i < poss.size(); i++){
            for (int j = 0; j < negs.size(); j++){
                // Merge clauses. If 'y' and '~y' exist, clause will not be created.
                dummy.clear();
                bool ok = merge(poss[i], negs[j], Lit(x), ~Lit(x), seen_tmp, dummy.asVec());
                if (ok){
                    after_clauses++;
                    /**/if (after_clauses > before_clauses) goto Abort;
                    after_literals += dummy.size(); }
            }
        }
      Abort:;

        // Maybe eliminate:
        if ((!opt_niver && after_clauses  <= before_clauses)
        ||  ( opt_niver && after_literals <= before_literals)
        ){
            MigrateToPsNs
            for (int i = 0; i < ps.size(); i++){
                for (int j = 0; j < ns.size(); j++){
                    dummy.clear();
                    bool ok = merge(ps[i], ns[j], Lit(x), ~Lit(x), seen_tmp, dummy.asVec());
                    if (ok){
                        Clause c = addClause(dummy.asVec());
                        if (!c.null()){
                            new_clauses.push(c); }
                        propagateToplevel(); if (!ok) return true;
                    }
                }
            }
            DeallocPsNs
            goto Eliminated;
        }

        /*****TEST*****/
#if 1
        // Try to remove 'x' from clauses:
        bool    ran = false;
        if (poss.size() < 10){ ran = true; asymmetricBranching( Lit(x)); if (!ok) return true; }
        if (negs.size() < 10){ ran = true; asymmetricBranching(~Lit(x)); if (!ok) return true; }
        if (value(x) != l_Undef) return false;
        if (!ran) return false;

        {
            // Count clauses/literals after elimination:
            int after_clauses  = 0;
            int after_literals = 0;
            Clause_t  dummy;
            for (int i = 0; i < poss.size(); i++){
                for (int j = 0; j < negs.size(); j++){
                    // Merge clauses. If 'y' and '~y' exist, clause will not be created.
                    dummy.clear();
                    bool ok = merge(poss[i], negs[j], Lit(x), ~Lit(x), seen_tmp, dummy.asVec());
                    if (ok){
                        after_clauses++;
                        after_literals += dummy.size(); }
                }
            }

            // Maybe eliminate:
            if ((!opt_niver && after_clauses  <= before_clauses)
            ||  ( opt_niver && after_literals <= before_literals)
            ){
                MigrateToPsNs
                for (int i = 0; i < ps.size(); i++){
                    for (int j = 0; j < ns.size(); j++){
                        dummy.clear();
                        bool ok = merge(ps[i], ns[j], Lit(x), ~Lit(x), seen_tmp, dummy.asVec());
                        if (ok){
                            Clause c = addClause(dummy.asVec());
                            if (!c.null()){
                                new_clauses.push(c); }
                            propagateToplevel(); if (!ok) return true;
                        }
                    }
                }
                DeallocPsNs
                goto Eliminated;
            }
        }
#endif
        /*****END TEST*****/
    }

    return false;

  Eliminated:
    assert(occur[index(Lit(x))].size() + occur[index(~Lit(x))].size() == 0);
    var_elimed[x] = 1;
    assigns   [x] = toInt(l_Error);
    return true;
}


// Inefficient test implementation! (pre-condition: satisfied clauses have been removed)
//
void Solver::clauseReduction(void)
{
    /**/reportf("clauseReduction() -- BEGIN\n");
    assert(decisionLevel() == 0);
    propagateToplevel(); if (!ok){ Report("(contradiction during clause reduction)\n"); return; }

    vec<Lit>    learnt;
    vec<char>&  seen = seen_tmp;

    for (int i = 0; i < constrs.size(); i++){
        Clause  c = constrs[i]; if (c.null()) continue;
        Clause  confl;
        bool    analyze = false;

        //**/dump(*this, c);
        unwatch(c, ~c[0]);
        unwatch(c, ~c[1]);
        for (int j = 0; j < c.size(); j++){
            //**/printf(L_LIT" = %c\n", L_lit(c[j]), name(value(c[j])));

            if (!assume(~c[j])){
                learnt.clear();
                if (level[var(c[j])] > 0)
                    learnt.push(c[j]);
                Clause rs = reason[var(c[j])];
                for (int k = 1; k < rs.size(); k++){
                    assert(value(rs[k]) == l_False);
                    seen[index(~rs[k])] = 1; }
                analyze = true;
                break;
            }

            confl = propagate();
            if (!confl.null()){
                learnt.clear();
                for (int k = 0; k < confl.size(); k++){
                    assert(value(confl[k]) == l_False);
                    seen[index(~confl[k])] = 1; }
                analyze = true;
                break;
            }
        }

        if (analyze){
            // Analyze conflict:
            for (int j = trail.size()-1; j >= 0; j--){
                Lit p = trail[j];
                if (seen[index(p)]){
                    Clause rs = reason[var(p)];
                    if (rs.null()){
                        if (level[var(p)] > 0)
                            learnt.push(~p);
                    }else{
                        for (int k = 1; k < rs.size(); k++)
                            seen[index(~rs[k])] = 1;
                    }
                    seen[index(p)] = 0;
                }
            }

            // (kolla att seen[] �r nollad korrekt h�r)

            /**/if (learnt.size() < c.size()){
                putchar('*'); fflush(stdout);
                //**/printf("Original: "); dump(c);
                //**/printf("Conflict:");
                //**/for (int j = 0; j < learnt.size(); j++) printf(" "L_LIT, L_lit(learnt[j]));
                //**/printf("\n");
                //**/if (learnt.size() == 0) exit(0);
            /**/}
        }

        cancelUntil(0);
        watch(c, ~c[0]);
        watch(c, ~c[1]);

        if (analyze && learnt.size() < c.size()){
            // Update clause:
            //**/printf("{"); for (int i = 0; i < c.size(); i++) printf(" %d", occur[index(c[i])].size()); printf(" }\n");

            /**/assert(propQ.size() == 0);
            /*Test Subset*/
            for (int i = 0; i < learnt.size(); i++){
                for (int j = 0; j < c.size(); j++)
                    if (c[j] == learnt[i])
                        goto Found;
                assert(false);  // no subset!
              Found:;
            }
            /*End*/
            unlinkClause(c);    // (will touch all variables)
            addClause(learnt, false, c);
            propagateToplevel(); if (!ok){ Report("(contradiction during clause reduction after learning unit fact)\n"); return; }
            /**/assert(propQ.size() == 0);
        }
        //**/else{ printf("                    {"); for (int i = 0; i < c.size(); i++) printf(" %d", occur[index(c[i])].size()); printf(" }\n"); }
    }
    /**/reportf("clauseReduction() -- END\n");
}
//#################################################################################################
// "Solver_subsume.C" END
//#################################################################################################

}// end namespace SatELite
