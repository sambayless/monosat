/***************************************************************************************[Solver.cc]
 The MIT License (MIT)

 Copyright (c) 2014-2017, Sam Bayless
 Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
 Copyright (c) 2007-2010, Niklas Sorensson

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

#include <cmath>
#include "monosat/mtl/Alg.h"
#include <algorithm>
#include "monosat/mtl/Sort.h"
#include "monosat/graph/GraphTheory.h"
#include <ctype.h>

using namespace Monosat;
#ifndef NDEBUG
#define DEBUG_SOLUTION
#endif
//=================================================================================================
// Options:
// Collected in Config.h
//=================================================================================================
// Constructor/Destructor:
bool Solver::shown_warning = false;

Solver::Solver() :

// Parameters (user settable):
//
        verbosity(opt_verb), var_decay(opt_var_decay), clause_decay(opt_clause_decay), theory_decay(opt_var_decay),
        random_var_freq(
                opt_random_var_freq), random_seed(opt_random_seed), luby_restart(opt_luby_restart), ccmin_mode(
        opt_ccmin_mode), phase_saving(opt_phase_saving), rnd_pol(false), rnd_init_act(opt_rnd_init_act), garbage_frac(
        opt_garbage_frac), restart_first(opt_restart_first), restart_inc(opt_restart_inc)

        // Parameters (the rest):
        //
        , learntsize_factor((double) 1 / (double) 3), learntsize_inc(1.1)

        // Parameters (experimental):
        //
        , learntsize_adjust_start_confl(100), learntsize_adjust_inc(1.5)

        // Statistics: (formerly in 'SolverStats')
        //
        , solves(0), starts(0), decisions(0), rnd_decisions(0), propagations(0), conflicts(0), stats_pure_lits(
        0), stats_pure_theory_lits(0), pure_literal_detections(0), stats_removed_clauses(0), dec_vars(0),
        clauses_literals(
                0), learnts_literals(0), max_literals(0), tot_literals(0), stats_pure_lit_time(0), ok(
        true), cla_inc(1), var_inc(1), theory_inc(1), watches(WatcherDeleted(ca)), qhead(0), simpDB_assigns(-1),
        simpDB_props(
                0), order_heap(VarOrderLt(activity, priority)), theory_order_heap(HeuristicOrderLt(), HeuristicToInt()),
        progress_estimate(0), remove_satisfied(true) //lazy_heap( LazyLevelLt(this)),

        // Resource constraints:
        //
        , conflict_budget(-1), propagation_budget(-1){
    all_decision_heuristics.push(nullptr);//prevent any decision heuristic from getting a heuristic id of 0
    if(opt_vsids_solver_as_theory){
        decisionTheory = new SolverDecisionTheory(*this);
        this->addHeuristic(decisionTheory);
    }
}

Solver::~Solver(){
    for(Theory* t : theories){
        delete (t);
    }
    delete pbsolver;
}

const std::string Solver::empty_name = "";
//=================================================================================================
// Minor methods:

// Creates a new SAT variable in the solver. If 'decision' is cleared, variable will not be
// used as a decision variable (NOTE! This has effects on the meaning of a SATISFIABLE result).
//
Var Solver::newVar(bool sign, bool dvar){

    int v;
    if(free_vars.size() > 0){
        v = free_vars.last();
        free_vars.pop();
    }else{
        v = nVars();
        assigns.push(l_Undef);
        vardata.push();
        priority.push();
        theory_vars.push();
        activity.push();
        seen.push(0);
        polarity.push();
        decision.push();
        trail.capacity(v + 1);
    }
    watches.init(mkLit(v, false));
    watches.init(mkLit(v, true));
    assigns[v] = l_Undef;
    vardata[v] = mkVarData(CRef_Undef, 0);
    int p = 0;
    if(max_decision_var > 0 && v > max_decision_var)
        p = 1;
    if(v < min_decision_var)
        p = 1;
    priority[v] = p;
    theory_vars[v].theoryMap.clear(true);
    theory_vars[v].theories.clear(true);
    activity[v] = (rnd_init_act ? drand(random_seed) * 0.00001 : 0);
    seen[v] = 0;
    polarity[v] = opt_init_rnd_phase ? irand(random_seed, 1) : sign;
    //decision.push();//set below

    if(max_decision_var > 0 && v > max_decision_var)
        dvar = false;
    if(v < min_decision_var)
        dvar = false;
    setDecisionVar(v, dvar);

    return v;
}


// Note: at the moment, only unassigned variable will be released (this is to avoid duplicate
// releases of the same variable).
void Solver::releaseVar(Lit l){
    assert(!hasTheory(
            var(l)));//theory vars should never be released (unless the theory is already detached from the variable?)
    if(value(l) == l_Undef && !hasTheory(l)){//dont ever release theory atoms
        released_vars.push(var(l));
    }
    addClause(l);
}

void Solver::checkClause(vec<Lit>& clause){
#ifdef DEBUG_SOLUTION
    //start with an easy check - does this clause trivially violate the known solution
    for(LSet& solution:known_solutions){
        bool found = false;
        for(Lit l:clause){
            if(!solution.contains(~l)){
                found = true;
                break;
            }
        }
        if(!found){
            throw std::runtime_error("Learned clause violates known solution");
        }
    }
    //slow check
    // Solver dbg_Solver;
#endif
}

bool Solver::addClause_(vec<Lit>& ps, bool is_derived_clause){


    assert(decisionLevel() == 0);
    if(!ok)
        return false;
    for(Lit l:ps){
        if(hasTheory(l)){
            resetInitialPropagation();    //Ensure that super solver call propagate on this solver at least once.
        }
    }
    // Check if clause is satisfied and remove false/duplicate literals:
    sort(ps);
    Lit p;
    int i, j;
    for(i = j = 0, p = lit_Undef; i < ps.size(); i++)
        if(value(ps[i]) == l_True || ps[i] == ~p)
            return true;
        else if(value(ps[i]) != l_False && ps[i] != p)
            ps[j++] = p = ps[i];
    ps.shrink(i - j);
    checkClause(ps);
    if(ps.size() == 0)
        return ok = false;
    else if(ps.size() == 1){
        uncheckedEnqueue(ps[0]);
        ok = (propagate(false) ==
              CRef_Undef); //do NOT propagate theory solvers here, or else adding unit clauses can become very expensive in some circumstances (such as when constructing the initial CNF for example)

        return ok;
    }else{
        CRef cr = ca.alloc(ps, false);
        //mark this clause as being either a learned clause, or a clause derived from a theory, or a clause derived from simplification/preprocessing.
        ca[cr].setDerived(is_derived_clause);
        clauses.push(cr);
        attachClause(cr);
    }

    return true;
}

CRef Solver::attachReasonClause(Lit r, vec<Lit>& ps){
    assert(value(r) == l_True);

    if(opt_write_learnt_clauses){
        fprintf(opt_write_learnt_clauses, "learnt ");
        for(Lit l:ps){
            fprintf(opt_write_learnt_clauses, " %d", dimacs(unmap(l)));
        }
        fprintf(opt_write_learnt_clauses, " 0\n");
        fflush(opt_write_learnt_clauses);
    }

    //sort(ps);
    Lit p;
    int i, j;
    for(i = j = 0, p = lit_Undef; i < ps.size(); i++)
        if(((value(ps[i]) == l_True && level(var(ps[i])) == 0)) || ps[i] == ~p)
            return CRef_Undef;
        else if((value(ps[i]) != l_False || level(var(ps[i])) != 0) && ps[i] != p)
            ps[j++] = p = ps[i];
    ps.shrink(i - j);

    CRef confl_out = CRef_Undef;
    if(ps.size() == 0){
        ok = false;
        //cancelUntil(0);
        return CRef_Undef;
    }else if(ps.size() == 1){
        //cancelUntil(0);
        assert(var(ps[0]) < nVars());

        enqueueLazy(ps[0], 0);

        return CRef_Undef;
    }else{

        int nfalse = 0;
        int max_lev = 0;
        //bool satisfied = false;
        int notFalsePos1 = -1;
        int notFalsePos2 = -1;
        for(int j = 0; j < ps.size(); j++){
            assert(var(ps[j]) < nVars());
            //assert(value(ps[j])==l_False);
            if(value(ps[j]) == l_False){
                nfalse++;
            }else{

                //if (value(ps[j]) == l_True)
                //	satisfied = true;
                if(notFalsePos1 < 0)
                    notFalsePos1 = j;
                else if(notFalsePos2 < 0){
                    notFalsePos2 = j;
                }
            }
            if(value(ps[j]) != l_Undef){
                int l = level(var(ps[j]));
                if(l > max_lev){
                    max_lev = l;
                }
            }
        }

        if(notFalsePos1 < 0){
            //all literals in this clause are false - this is the usual case for reason clauses.
            //so just treat the first literal as the non-false one
            notFalsePos1 = 0;
        }

        assert(notFalsePos1 >= 0);
        if(notFalsePos1 >= 0 && notFalsePos2 >= 0){
            assert(notFalsePos1 != notFalsePos2);
            if(notFalsePos1 == 1){
                std::swap(ps[0], ps[notFalsePos2]);
            }else{
                std::swap(ps[0], ps[notFalsePos1]);
                std::swap(ps[1], ps[notFalsePos2]);
            }
        }else{
            std::swap(ps[0], ps[notFalsePos1]);
        }
        assert(value(ps[0]) != l_False);
        if(notFalsePos2 >= 0){
            assert(value(ps[1]) != l_False);
        }


        CRef cr = ca.alloc(ps);
        ca[cr].setDerived(true);
        clauses.push(cr);
        attachClause(cr);
        enqueueLazy(ps[0], max_lev, cr);

        return cr;
    }
}

void Solver::attachClause(CRef cr){
    const Clause& c = ca[cr];
    assert(c.size() > 1);
#ifdef DEBUG_CORE
    for(int p = 0; p <= 1; p++)
        if(value(c[p]) == l_False && level(var(c[p])) == 0){
            //then c[0] must not have been propagated yet - it must be after qhead. Otherwise this clause will never be enforced
            if(decisionLevel() > 0){
                assert(false); //because we have already propagated all the level 0 elements
            }else{
                bool found = false;
                for(int i = qhead; i < trail.size(); i++){
                    if(var(trail[i]) == var(c[p])){
                        found = true;
                        break;
                    }
                }
                assert(found);
            }
        }
#endif
    watches[~c[0]].push(Watcher(cr, c[1]));
    watches[~c[1]].push(Watcher(cr, c[0]));
    if(c.learnt())
        learnts_literals += c.size();
    else
        clauses_literals += c.size();
}

void Solver::detachClause(CRef cr, bool strict){
    const Clause& c = ca[cr];
    assert(c.size() > 1);

    if(strict){
        remove(watches[~c[0]], Watcher(cr, c[1]));
        remove(watches[~c[1]], Watcher(cr, c[0]));
    }else{
        // Lazy detaching: (NOTE! Must clean all watcher lists before garbage collecting this clause)
        watches.smudge(~c[0]);
        watches.smudge(~c[1]);
    }

    if(c.learnt())
        learnts_literals -= c.size();
    else
        clauses_literals -= c.size();
}

void Solver::removeClause(CRef cr){
    CRef remove_clause = cr;
    Clause& c = ca[cr];
    detachClause(cr);
    // Don't leave pointers to free'd memory!
    if(locked(c))
        vardata[var(c[0])].reason = CRef_Undef;
    c.mark(1);
    ca.free(cr);
}

bool Solver::satisfied(const Clause& c) const{
    for(int i = 0; i < c.size(); i++)
        if(value(c[i]) == l_True)
            return true;
    return false;
}

// Revert to the state at given level (keeping all assignment at 'level' but not beyond).
//
void Solver::cancelUntil(int lev){
    if(decisionLevel() > lev){

        for(int i = 0; i < theories.size(); i++){
            if(opt_lazy_backtrack && theories[i]->supportsLazyBacktracking()){
                //if we _are_ backtracking lazily, then the assumption is that the theory solver will, after backtracking, mostly re-assign the same literals.
                //so instead, we will backtrack the theory lazily, in the future, if it encounters an apparent conflict (and this backtracking may alter or eliminate that conflict.)
            }else{
                theories[i]->backtrackUntil(lev);
            }
        }

        //printf("s: cancel %d\n",lev);
        for(int c = trail.size() - 1; c >= trail_lim[lev]; c--){
            Var x = var(trail[c]);
            int xlev = level(x);

            if(xlev <= lev){
                to_reenqueue.push(trail[c]);
            }else{
                //if(hasTheory(x)){
                for(int n = 0; n < getNTheories(x); n++){
                    int theoryID = getTheoryID(x, n);
                    Lit l = getTheoryLit(trail[c], n);

                    if(c <= satisfied_theory_trail_pos[theoryID]){
                        satisfied_theory_trail_pos[theoryID] = -1;
                        post_satisfied_theory_trail_pos[theoryID] = -1;
                        //printf("theory %d no longer sat at lev %d\n",theoryID, decisionLevel());
                    }

                    if(satisfied_theory_trail_pos[theoryID] < 0){
                        theories[theoryID]->undecideTheory(l);
                    }else if(c > satisfied_theory_trail_pos[theoryID] &&
                             c <= post_satisfied_theory_trail_pos[theoryID]){
                        theories[theoryID]->undecideTheory(l);
                        post_satisfied_theory_trail_pos[theoryID] = c - 1;
                    }
                    assert(satisfied_theory_trail_pos[theoryID] < c);
                    assert(post_satisfied_theory_trail_pos[theoryID] < c);
                }
                //}
                assigns[x] = l_Undef;
                if(phase_saving > 1 || ((phase_saving == 1) && c > trail_lim.last()))
                    polarity[x] = sign(trail[c]);
                insertVarOrder(x);
            }
        }
        qhead = trail_lim[lev];
        if(local_qhead > qhead){
            local_qhead = qhead;
        }
        if(S && super_qhead > S->qhead){
            super_qhead = S->qhead;
        }

        trail.shrink(trail.size() - trail_lim[lev]);
        trail_lim.shrink(trail_lim.size() - lev);

        if(qhead > trail.size()){
            throw std::runtime_error("Internal error in backtracking");
        }

        while(decision_heuristic_trail.size() &&
              first_heuristic_decision_level[decision_heuristic_trail.last()->getHeuristicIndex()] > decisionLevel()){
            first_heuristic_decision_level[decision_heuristic_trail.last()->getHeuristicIndex()] = -1;
            decision_heuristic_trail.pop();
        }

        for(int i = 0; i < theories.size(); i++){
            int theoryID = theories[i]->getTheoryIndex();
            int c = trail.size() - 1;
            if(c <= satisfied_theory_trail_pos[theoryID]){
                satisfied_theory_trail_pos[theoryID] = -1;
                post_satisfied_theory_trail_pos[theoryID] = -1;
                //printf("theory %d no longer sat at lev %d\n",theoryID, decisionLevel());
            }else if(c > satisfied_theory_trail_pos[theoryID] && c <= post_satisfied_theory_trail_pos[theoryID]){
                theories[theoryID]->undecideTheory(getTheoryLit(trail[c], theories[theoryID]));
                post_satisfied_theory_trail_pos[theoryID] = c - 1;
            }
        }

        /*while(theory_sat_queue.size() && theory_sat_queue.last().trail_size>trail.size()){
			int theoryID = theory_sat_queue.last().theoryID;
			assert(theorySatisfied(theories[theoryID]));
			satisfied_theory_trail_pos[theoryID]=-1;
			theory_sat_queue.pop();
		}*/

        //remove any lits from the lazy heap that are now unassigned.
/*		while(lazy_heap.size() && value(toLit( lazy_heap.peekMin())) == l_Undef ){
			lazy_heap.pop();
		}*/

        if(decisionLevel() < track_min_level){
            track_min_level = decisionLevel();
        }
        //you cannot necessarily remove a theory from the theory queue on backtracking; it should only be removed following a successful propagate.
        /*for (int i : theory_queue) {
			in_theory_queue[i] = false;
		}
		theory_queue.clear();*/

        int lowest_re_enqueue = -1;
        for(int i = 0; i < theories.size(); i++){
            int theoryID = theories[i]->getTheoryIndex();
            if(qhead < theory_init_prop_trail_pos[theoryID]){
                theory_init_prop_trail_pos[theoryID] = -1;
                theory_reprop_trail_pos[theoryID] = -1;
            }else if(theory_init_prop_trail_pos[theoryID] >= 0 && qhead < theory_reprop_trail_pos[theoryID]){
                theory_reprop_trail_pos[theoryID] = qhead;
                assert(theory_init_prop_trail_pos[theoryID] >= 0);
                int init_pos = theory_init_prop_trail_pos[theoryID];
                if(trail_lim.size() > 0 && init_pos < trail_lim[decisionLevel()]){
                    init_pos = trail_lim[decisionLevel()];
                }
                Lit p = trail[init_pos];
                assert(p != lit_Undef);
                int back_lev = level(var(p)) - 1;
                if(back_lev < 0)
                    back_lev = 0;
                if(lowest_re_enqueue < 0 || trail_lim[back_lev] < lowest_re_enqueue){
                    lowest_re_enqueue = trail_lim[back_lev];
                }
                theories[theoryID]->backtrackUntil(back_lev);//or back_lev -1?
            }

        }
        if(lowest_re_enqueue > -1){
            for(int q = lowest_re_enqueue; q < trail.size(); q++){ //should this be q<qhead, or q<=qhead?
                Lit p = trail[q];
                for(int n = 0; n < getNTheories(var(p)); n++){
                    int theoryID = getTheoryID(p, n);
                    Lit theoryLit = getTheoryLit(p, n);
                    theories[theoryID]->backtrackUntil(decisionLevel());
                    if(theory_reprop_trail_pos[theoryID] == -1 && q >= theory_init_prop_trail_pos[theoryID] &&
                       !theorySatisfied(theories[theoryID])){
                        needsPropagation(theoryID);
                        //theories[theoryID]->backtrackUntil(level(var(p)));
                        theories[theoryID]->enqueueTheory(theoryLit);
                    }
                }
            }
        }

        //re-enqueue any lazy lits as appropriate
        while(to_reenqueue.size()){
            Lit p = to_reenqueue.last();
            to_reenqueue.pop();
            assert(level(var(p)) <= lev);
            trail.push(p);
            //is this really needed?
            for(int n = 0; n < getNTheories(var(p)); n++){
                int theoryID = getTheoryID(p, n);
                Lit l = getTheoryLit(p, n);
                if(!theorySatisfied(theories[theoryID])){
                    needsPropagation(theoryID);
                    theories[theoryID]->enqueueTheory(l);
                }
            }
        }
        //should qhead be adjusted, here? Or do we want to repropagate these literals? (currently, re-propagating these literals).

        //If opt_theory_order_vsids is enabled, need to add back into the theory decision heap theories that may now be able to make decisions again.
        while(theory_decision_trail.size() && theory_decision_trail.last().second >= decisionLevel()){
            Heuristic* h = theory_decision_trail.last().first;
            if(!theory_order_heap.inHeap(h)){
                theory_order_heap.insert(h);
            }
            theory_decision_trail.pop();
        }


    }

    while(partially_propagated_levels.size() && partially_propagated_levels.last() > decisionLevel()){
        partially_propagated_levels.pop();
    }
    if(partially_propagated_levels.size() && partially_propagated_levels.last() == decisionLevel()){
        //mark all theories as needing propagation. in principle, we could try to keep track of exactly which theories might need propagation, but that could be expensive (memory-wise)
        partially_propagated_levels.pop();
        for(int i = 0; i < theories.size(); i++){
            if(theories[i] && !theories[i]->unskipable())
                needsPropagation(i);
        }

    }
    assert(partially_propagated_levels.size() == 0 || partially_propagated_levels.last() < decisionLevel());
}

void Solver::backtrackUntil(int lev){
    if(S){
        if(S->trail.size() < super_qhead)
            super_qhead = S->trail.size();
    }
    cancelUntil(0);
}

//=================================================================================================
// Major methods:

Lit Solver::pickBranchLit(){
    Var next = var_Undef;
    decisions++;
    // Random decision:
    if(drand(random_seed) < random_var_freq && !order_heap.empty()){
        next = order_heap[irand(random_seed, order_heap.size())];
        if(value(next) == l_Undef && decision[next])
            rnd_decisions++;
    }

    // Activity based decision:
    while(next == var_Undef || value(next) != l_Undef || !decision[next])
        if(order_heap.empty()){
            next = var_Undef;
            break;
        }else
            next = order_heap.removeMin();

    return next == var_Undef ? lit_Undef : mkLit(next, rnd_pol ? drand(random_seed) < 0.5 : polarity[next]);
}

void Solver::instantiateLazyDecision(Lit p, int atLevel, CRef reason){
    assert(value(p) == l_Undef);
    assigns[var(p)] = lbool(!sign(p));
    vardata[var(p)] = mkVarData(reason, atLevel);
    assert(atLevel <= decisionLevel());
    assert(atLevel > 0);
    int trail_pos = trail_lim[atLevel - 1];
    Lit curDec = trail[trail_pos];
    if(curDec == p)//this lit was already the decision at this level
        return;
    assert(curDec == theoryDecision);
    if(curDec != theoryDecision){
        throw std::runtime_error("Critical error: bad decision");
    }
    trail[trail_pos] = p;
    for(int n = 0; n < getNTheories(var(p)); n++){
        int theoryID = getTheoryID(p, n);
        Lit l = getTheoryLit(p, n);
        if(!theorySatisfied(theories[theoryID])){
            needsPropagation(theoryID);
            theories[theoryID]->enqueueTheory(l);
        }
    }
}

void Solver::analyzeHeuristicDecisions(CRef confl, IntSet<int>& conflicting_heuristics, int max_involved_heuristics,
                                       int minimimum_involved_decision_priority){
    if(max_involved_heuristics > -1 && conflicting_heuristics.size() >= max_involved_heuristics)
        return;

    int pathC = 0;
    CRef original_confl = confl;
    Lit p = lit_Undef;
    assert(confl != CRef_Undef);

#ifdef DEBUG_CORE
    assert(!seen.contains(1));
#endif

    int index = trail.size() - 1;
    bool found_any = false;
    do{
        if(confl != CRef_Undef){
            assert(!isTheoryCause(confl));
            Clause& c = ca[confl];

            for(int j = (p == lit_Undef) ? 0 : 1; j < c.size(); j++){
                Lit q = c[j];
                if(!seen[var(q)] && level(var(q)) > 0){
                    seen[var(q)] = 1;
                    pathC++;
                }
            }
        }
        confl = CRef_Undef;


        // Select next clause to look at:
        while(index >= 0 && (!seen[var(trail[index--])]));
        assert(index >= -1);
        p = trail[index + 1];

        assert(var(p) != var(theoryDecision));

        confl = reasonOrDecision(var(p));
        if(confl == CRef_Undef){
            //do nothing
        }else if(isDecisionReason(confl)){
            Heuristic* h = getHeuristic(confl);
            if(!conflicting_heuristics.has(h->getHeuristicIndex()) &&
               h->getPriority() >= minimimum_involved_decision_priority){
                conflicting_heuristics.insert(h->getHeuristicIndex());
                if(max_involved_heuristics > -1 && conflicting_heuristics.size() >= max_involved_heuristics){
                    for(int i = 0; i <= index + 1; i++){
                        seen[var(trail[i])] = false;
                    }
#ifdef DEBUG_CORE
                    assert(!seen.contains(1));
#endif
                    return;
                }
            }
            confl = CRef_Undef;
        }else if(isTheoryCause(confl)){
            Theory* t = getTheory(confl);
            assert(t);
            //this is not quite right...
            if(t->getHeuristicIndex() >= 0){
                if(!conflicting_heuristics.has(t->getHeuristicIndex()) &&
                   t->getPriority() >= minimimum_involved_decision_priority){
                    conflicting_heuristics.insert(t->getHeuristicIndex());
                    if(max_involved_heuristics > -1 && conflicting_heuristics.size() >= max_involved_heuristics){
                        for(int i = 0; i <= index + 1; i++){
                            seen[var(trail[i])] = false;
                        }
#ifdef DEBUG_CORE
                        assert(!seen.contains(1));
#endif
                        return;
                    }
                }
            }
            confl = CRef_Undef;
        }

        seen[var(p)] = 0;
        pathC--;
    }while(pathC > 0);
    assert(max_involved_heuristics < 0 || conflicting_heuristics.size() < max_involved_heuristics);
    if(p != lit_Undef){
        confl = reasonOrDecision(var(p));
        if(confl == CRef_Undef){
            //do nothing

        }else if(isDecisionReason(confl)){
            Heuristic* h = getHeuristic(confl);
            if(!conflicting_heuristics.has(h->getHeuristicIndex()) &&
               h->getPriority() >= minimimum_involved_decision_priority){
                conflicting_heuristics.insert(h->getHeuristicIndex());
            }

        }else{
            if(isTheoryCause(confl)){
                Theory* t = getTheory(confl);
                assert(t);
                //this is not quite right...
                if(t->getHeuristicIndex() >= 0){
                    if(!conflicting_heuristics.has(t->getHeuristicIndex()) &&
                       t->getPriority() >= minimimum_involved_decision_priority){
                        conflicting_heuristics.insert(t->getHeuristicIndex());
                    }
                }

            }
        }
    }

#ifdef DEBUG_CORE
    assert(!seen.contains(1));
#endif
}

/*_________________________________________________________________________________________________
 |
 |  analyze : (confl : Clause*) (out_learnt : vec<Lit>&) (out_btlevel : int&)  ->  [void]
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
 |      * If out_learnt.size() > 1 then 'out_learnt[1]' has the greatest decision level of the
 |        rest of literals. There may be others from the same level though.
 |
 |________________________________________________________________________________________________@*/
void Solver::analyze(CRef confl, vec<Lit>& out_learnt, int& out_btlevel){
    int pathC = 0;
    CRef original_confl = confl;
    Lit p = lit_Undef;
    assert(confl != CRef_Undef);
    int maxlev = 0;
    {
        Clause& check = ca[confl];
        for(Lit p:check){
            int lev = level(var(p));
            if(lev > maxlev){
                maxlev = lev;
            }
        }
    }

    cancelUntil(maxlev);//use of lazily enqueued literals can trigger conflicts at earlier decision levels
    // Generate conflict clause:
    //
#ifdef DEBUG_CORE
    assert(!seen.contains(1));
#endif
    bool possibly_missed_1uip = false;
    to_analyze.clear();
    out_learnt.push(lit_Undef);      // (leave room for the asserting literal)
    int index = trail.size() - 1;
    int stop = trail_lim.last();
    do{
        if(confl != CRef_Undef){
            assert(!isTheoryCause(confl));
            Clause& c = ca[confl];

            if(c.learnt())
                claBumpActivity(c);

            for(int j = (p == lit_Undef) ? 0 : 1; j < c.size(); j++){
                Lit q = c[j];
                if(!seen[var(q)] && level(var(q)) > 0){
                    assert(value(q) == l_False);
                    assert(var(q) != var(theoryDecision));
                    varBumpActivity(var(q));
                    seen[var(q)] = 1;

                    if(level(var(q)) >= decisionLevel())
                        pathC++;
                    else
                        out_learnt.push(q);
                }
            }
        }else{
            out_learnt.push(~p);
        }
        bool searching = true;
        while(searching){
            searching = false;

            // Select next clause to look at:
            while(index >= stop && (!seen[var(trail[index--])] || (level(var(trail[index + 1])) < decisionLevel())));
            assert(index >= -1);
            p = trail[index + 1];

            assert(var(p) != var(theoryDecision));
            confl = reason(var(p));
            int was_at_level = level(var(p));
            if(isTheoryCause(confl)){
                //lazily construct the reason for this theory propagation now that we need it
                confl = constructReason(p);
                //for some theories, we may discover while constructing the cause that p is at a lower level than we thought.
            }
            if(level(var(p)) < decisionLevel()){
                assert(value(p) == l_True);
                if(was_at_level == decisionLevel()){
                    //the level of this variable changed while deriving a reason for it.
                    if(level(var(p)) > 0){
                        out_learnt.push(~p);
                    }
                    pathC--;
                    possibly_missed_1uip = true;
                }
                seen[var(p)] = 0;
                searching = pathC > 0;
            }
        }

        seen[var(p)] = 0;
        pathC--;

    }while(pathC > 0);
    out_learnt[0] = ~p;
    assert(var(p) != var(theoryDecision));

    if(possibly_missed_1uip){
        //because of literals that were enqueued lazily, at the wrong level, and only discovered during lazy reason construction,
        //in rare circumstances the 1uip may be missed. If that happens, re-start clause learning now that all the relevant reasons have
        //been constructed and all relevant levels are corrected.
        for(Lit l:out_learnt){
            seen[var(l)] = 0;
        }
        assert(!seen.contains(1));
        out_learnt.clear();
        analyze(original_confl, out_learnt, out_btlevel);
        return;
    }

#ifdef DEBUG_CORE
    for(Lit p : out_learnt)
        assert(value(p) == l_False);

#endif
    // Simplify conflict clause:
    //
    int i, j;
    out_learnt.copyTo(analyze_toclear);
    if(ccmin_mode == 2){
        uint32_t abstract_level = 0;
        for(i = 1; i < out_learnt.size(); i++)
            abstract_level |= abstractLevel(
                    var(out_learnt[i])); // (maintain an abstraction of levels involved in conflict)

        for(i = j = 1; i < out_learnt.size(); i++)
            if(!ca.isClause(reason(var(out_learnt[i]))) || !litRedundant(out_learnt[i], abstract_level))
                out_learnt[j++] = out_learnt[i];

    }else if(ccmin_mode == 1){
        for(i = j = 1; i < out_learnt.size(); i++){
            Var x = var(out_learnt[i]);

            if(reason(x) == CRef_Undef)
                out_learnt[j++] = out_learnt[i];
            else{
                Clause& c = ca[reason(var(out_learnt[i]))];
                for(int k = 1; k < c.size(); k++)
                    if(!seen[var(c[k])] && level(var(c[k])) > 0){
                        out_learnt[j++] = out_learnt[i];
                        break;
                    }
            }
        }
    }else
        i = j = out_learnt.size();

    max_literals += out_learnt.size();
    out_learnt.shrink(i - j);
    tot_literals += out_learnt.size();

    // Find correct backtrack level:
    //
    if(out_learnt.size() == 1)
        out_btlevel = 0;
    else{
        int max_i = 1;
        // Find the first literal assigned at the next-highest level:
        for(int i = 2; i < out_learnt.size(); i++)
            if(level(var(out_learnt[i])) > level(var(out_learnt[max_i])))
                max_i = i;
        // Swap-in this literal at index 1:
        Lit p = out_learnt[max_i];
        out_learnt[max_i] = out_learnt[1];
        out_learnt[1] = p;
        out_btlevel = level(var(p));
    }
#ifdef DEBUG_CORE
    for(Lit p : out_learnt){
        assert(var(p) != var(theoryDecision));
        assert(value(p) == l_False);
    }
#endif

    for(int j = 0; j < analyze_toclear.size(); j++)
        seen[var(analyze_toclear[j])] = 0;    // ('seen[]' is now cleared)
#ifdef DEBUG_CORE
    assert(!seen.contains(1));
#endif
}

// Check if 'p' can be removed. 'abstract_levels' is used to abort early if the algorithm is
// visiting literals at levels that cannot be removed later.
bool Solver::litRedundant(Lit p, uint32_t abstract_levels){
    analyze_stack.clear();
    analyze_stack.push(p);
    int top = analyze_toclear.size();
    while(analyze_stack.size() > 0){
        assert(reason(var(analyze_stack.last())) != CRef_Undef);
        if(!ca.isClause(reason(var(analyze_stack.last())))){
            //Don't pursue this if it is a subsolver reason
            for(int j = top; j < analyze_toclear.size(); j++)
                seen[var(analyze_toclear[j])] = 0;
            analyze_toclear.shrink(analyze_toclear.size() - top);
            analyze_stack.clear();
            return false;
        }

        Clause& c = ca[reason(var(analyze_stack.last()))];
        analyze_stack.pop();

        for(int i = 1; i < c.size(); i++){
            Lit p = c[i];
            if(!seen[var(p)] && level(var(p)) > 0){
                if(reason(var(p)) != CRef_Undef && (abstractLevel(var(p)) & abstract_levels) != 0){
                    seen[var(p)] = 1;
                    analyze_stack.push(p);
                    analyze_toclear.push(p);
                }else{
                    for(int j = top; j < analyze_toclear.size(); j++)
                        seen[var(analyze_toclear[j])] = 0;
                    analyze_toclear.shrink(analyze_toclear.size() - top);
                    return false;
                }
            }
        }
    }

    return true;
}

/*_________________________________________________________________________________________________
 |
 |  analyzeFinal : (p : Lit)  ->  [void]
 |
 |  Description:
 |    Specialized analysis procedure to express the final conflict in terms of assumptions.
 |    Calculates the (possibly empty) set of assumptions that led to the assignment of 'p', and
 |    stores the result in 'out_conflict'.
 |________________________________________________________________________________________________@*/
void Solver::analyzeFinal(Lit p, vec<Lit>& out_conflict){
    out_conflict.clear();
    out_conflict.push(p);

    if(decisionLevel() == 0)
        return;

    seen[var(p)] = 1;

    for(int i = trail.size() - 1; i >= trail_lim[0]; i--){
        Var x = var(trail[i]);
        if(seen[x]){
            if(reason(x) == CRef_Undef){
                assert(level(x) > 0);
                out_conflict.push(~trail[i]);
            }else{
                if(isTheoryCause(reason(x))){
                    constructReason(trail[i]);
                }
                //this can happen, if the reason was a theory reason of size 1 or 0
                if(reason(x) == CRef_Undef){
                    //note that this is NOT an assumption; it is really a theory implication that should be at level 0
                }else{

                    Clause& c = ca[reason(x)];
                    assert(var(c[0]) == x);
                    for(int j = 1; j < c.size(); j++)
                        if(level(var(c[j])) > 0)
                            seen[var(c[j])] = 1;
                }
            }
            seen[x] = 0;
        }
    }

    seen[var(p)] = 0;
}

void Solver::enqueueLazy(Lit p, int lev, CRef from){
    assert(value(p) != l_False);
    if(value(p) == l_True && lev < level(var(p))){
        //then the lit was already implied, but needs to be (lazily) moved to an earlier level.
        vardata[var(p)] = mkVarData(from, lev);
    }else if(value(p) == l_Undef && lev < decisionLevel()){
        assert(value(p) == l_Undef);
        assigns[var(p)] = lbool(!sign(p));
        vardata[var(p)] = mkVarData(from, lev);
        trail.push_(p);
        //lazy_heap.insert(toInt(p));

        for(int n = 0; n < getNTheories(var(p)); n++){
            int theoryID = getTheoryID(p, n);
            Lit l = getTheoryLit(p, n);
            if(!theorySatisfied(theories[theoryID])){
                needsPropagation(theoryID);
                theories[theoryID]->enqueueTheory(l);
            }
        }
    }else if(value(p) == l_Undef){
        uncheckedEnqueue(p, from);
    }else{
        //do nothing
    }
}

void Solver::uncheckedEnqueue(Lit p, CRef from){
    assert(value(p) == l_Undef);
    assigns[var(p)] = lbool(!sign(p));
    vardata[var(p)] = mkVarData(from, decisionLevel());
    trail.push_(p);
    for(int n = 0; n < getNTheories(var(p)); n++){
        int theoryID = getTheoryID(p, n);
        Lit l = getTheoryLit(p, n);
        if(!theorySatisfied(theories[theoryID])){
            needsPropagation(theoryID);
            theories[theoryID]->enqueueTheory(l);
        }
    }
}

void Solver::analyzeFinal(CRef confl, Lit skip_lit, vec<Lit>& out_conflict){
    out_conflict.clear();
    if(decisionLevel() == 0)
        return;

    Clause& c = ca[confl];

    for(int i = 0; i < c.size(); i++){
        Var x = var(c[i]);
        if(x == var(skip_lit))
            continue;

        assert(x >= 0);
        assert(x < nVars());
        if(level(x) > 0)
            seen[x] = 1;
        assert(value(x) != l_Undef);
    }

    int start = trail.size() - 1;
    int i;
    for(i = start; i >= trail_lim[0]; i--){
        Var x = var(trail[i]);
        assert(value(x) != l_Undef);
        assert(x >= 0);
        assert(x < nVars());

        if(seen[x]){
            assert(x != var(skip_lit));
            CRef r = reason(x);
            if(r == CRef_Undef){
                //Var v = var(trail[i]);
                //int lev = level(v);
                out_conflict.push(~trail[i]);
            }else{
                if(isTheoryCause(r)){
                    r = constructReason(trail[i]);
                }
                if(r == CRef_Undef){
                    //note that this is NOT an assumption; it is really a theory implication that should be at level 0
                }else{
                    Clause& c = ca[r];
                    for(int j = 0; j < c.size(); j++){

                        if(level(var(c[j])) > 0){
                            seen[var(c[j])] = 1;
                        }
                    }
                }
            }
            seen[x] = 0;
        }
    }
}

void Solver::buildReason(Lit p, vec<Lit>& reason){
    Lit local_l = fromSuper(p);    // mkLit(var(p)-super_offset, sign(p));
    analyzeFinal(local_l, reason);
    toSuper(reason, reason);
    interpolant.push();
    reason.copyTo(interpolant.last());    //Add this clause into the interpolant vector
}


void Solver::enqueueTheory(Lit l){

}

//Propagate assignments from the super solver's interface variables to this solver (and, if this solver makes further assignments to the interface, pass those back to the super solver)
bool Solver::propagateTheory(vec<Lit>& conflict_out){
    //backtrack as needed
    if(S->trail.size() && super_qhead > 0){
        Lit last_super = S->trail[super_qhead - 1];
        cancelUntil(S->level(var(last_super)));
    }else{
        cancelUntil(0);
    }
    assert(decisionLevel() <= S->decisionLevel());

    CRef confl = CRef_Undef;

    while(confl == CRef_Undef && super_qhead < S->qhead){
        Lit out_l = S->trail[super_qhead++];

        if(var(out_l) < min_super || var(out_l) > max_super)
            continue;    //this lit is not on the interface

        int lev = S->level(var(out_l));
        assert(decisionLevel() <= lev);
        while(decisionLevel() < lev){

            //We are going to start a new decision level; make sure the current one is fully propagated first
            initial_level = decisionLevel();
            track_min_level = initial_level;

            confl = propagate();
            if(confl != CRef_Undef || track_min_level < initial_level){
                cancelUntil(track_min_level);
                S->cancelUntil(track_min_level);
                goto conflict;
            }

            newDecisionLevel();
        }

        Lit local_l = fromSuper(out_l);

        if(!enqueue(local_l, CRef_Undef)){
            //this is a conflict
            conflict_out.clear();
            analyzeFinal(~local_l, conflict_out);
            toSuper(conflict_out, conflict_out);
            interpolant.push();
            conflict_out.copyTo(interpolant.last());    //Add this clause into the interpolant vector
            return false;
        }

    }
    initial_level = decisionLevel();
    track_min_level = initial_level;
    confl = propagate();
    if(confl != CRef_Undef || track_min_level < initial_level){
        cancelUntil(track_min_level);
        S->cancelUntil(track_min_level);
    }
    conflict:

    if(confl != CRef_Undef){
        //then we have a conflict which we need to instantiate in S
        conflict_out.clear();
        analyzeFinal(confl, lit_Undef, conflict_out);
        toSuper(conflict_out, conflict_out);
        interpolant.push();
        conflict_out.copyTo(interpolant.last());    //Add this clause into the interpolant vector
        return false;
    }else{
        //find lits to prop
        while(local_qhead < qhead){
            Lit local_l = trail[local_qhead++];
            if(var(local_l) < min_local || var(local_l) > max_local)
                continue;    //this lit is not on the interface
            Lit out_l = toSuper(local_l);
            lbool v = S->value(out_l);

            if(S->value(out_l) != l_True || S->level(var(out_l)) > level(var(local_l))){

                if(S->decisionLevel() > level(var(local_l))){
                    cancelUntil(level(var(local_l)));
                    S->cancelUntil(level(var(local_l)));

                }
                if(!S->enqueue(out_l, cause_marker)){    //can probably increment super_qhead here...
                    //this is a conflict
                    conflict_out.clear();
                    analyzeFinal(local_l, conflict_out);
                    toSuper(conflict_out, conflict_out);
                    interpolant.push();
                    conflict_out.copyTo(interpolant.last());            //Add this clause into the interpolant vector
                    return false;
                }
            }
        }

        while(decisionLevel() < S->decisionLevel())
            newDecisionLevel();

        assert(decisionLevel() == S->decisionLevel());
        return true;
    }
}

void Solver::enqueueAnyUnqueued(){
    //in some solving modes, a theory can sometimes delay enqueing some atoms. Check for any such atoms here.
    //so far, this should only happen if the theory was marked as 'satisfied'

    int startingPos = -1;
    for(Theory* t:theories){
        assert(t);
        if(theorySatisfied(t)){
            int start = post_satisfied_theory_trail_pos[t->getTheoryIndex()];
            assert(post_satisfied_theory_trail_pos[t->getTheoryIndex()] >=
                   satisfied_theory_trail_pos[t->getTheoryIndex()]);
            if(start >= 0 && (start < startingPos || startingPos < 0)){
                startingPos = start;
            }
        }
    }
    if(startingPos < 0){
        for(Theory* t:theories){
            t->enqueueAnyUnqueued();
        }
        return;
    }
    int lev_0_pos;
    if(trail_lim.size()){
        lev_0_pos = trail_lim[0];
    }else{
        lev_0_pos = trail.size();
    }
    for(int i = 0; i < trail.size(); i++){
        Lit p = trail[i];

        for(int n = 0; n < getNTheories(var(p)); n++){
            int theoryID = getTheoryID(p, n);
            Lit l = getTheoryLit(p, n);


            int start = post_satisfied_theory_trail_pos[theoryID];
            if(start >= 0 && start <= i){
                Theory* theory = theories[theoryID];

                theory->enqueueTheory(l);
                assert(post_satisfied_theory_trail_pos[theoryID] <= i);
                post_satisfied_theory_trail_pos[theoryID] = i;
                assert(post_satisfied_theory_trail_pos[theoryID] >= satisfied_theory_trail_pos[theoryID]);
            }
        }
    }
    for(Theory* t:theories){
        assert(t);
        assert(post_satisfied_theory_trail_pos[t->getTheoryIndex()] >= satisfied_theory_trail_pos[t->getTheoryIndex()]);
        if(theorySatisfied(t)){
            int start = satisfied_theory_trail_pos[t->getTheoryIndex()];
            if(start >= 0 && start < lev_0_pos - 1){
                satisfied_theory_trail_pos[t->getTheoryIndex()] = lev_0_pos - 1;
            }
        }
    }
    for(Theory* t:theories){
        t->enqueueAnyUnqueued();
    }
}


/*_________________________________________________________________________________________________
 |
 |  propagate : [void]  ->  [Clause*]
 |
 |  Description:
 |    Propagates all enqueued facts. If a conflict arises, the conflicting clause is returned,
 |    otherwise CRef_Undef.
 |
 |    Post-conditions:
 |      * the propagation queue is empty, even if there was a conflict.
 |________________________________________________________________________________________________@*/
CRef Solver::propagate(bool propagate_theories){
    if(qhead == trail.size() && (!initialPropagate || decisionLevel() > 0) && !unskippable_theory_q.size() &&
       (!propagate_theories ||
        !theory_queue.size())){//it is possible that the theory solvers need propagation, even if the sat solver has an empty queue.
        return CRef_Undef;
    }
    conflicting_heuristic = nullptr;
    CRef confl = CRef_Undef;
    int num_props = 0;
    int initial_qhead = qhead;
    watches.cleanAll();
    if(decisionLevel() == 0 && !propagate_theories){
        initialPropagate = true;//we will need to propagate this assignment to the theories at some point in the future.
    }
    do{

        while(qhead < trail.size()){
            Lit p = trail[qhead++];     // 'p' is enqueued fact to propagate.

            vec<Watcher>& ws = watches[p];
            Watcher* i, * j, * end;
            num_props++;
            for(i = j = (Watcher*) ws, end = i + ws.size(); i != end;){
                // Try to avoid inspecting the clause:
                Lit blocker = i->blocker;
                if(value(blocker) == l_True){
                    *j++ = *i++;
                    continue;
                }

                // Make sure the false literal is data[1]:
                CRef cr = i->cref;
                Clause& c = ca[cr];
                Lit false_lit = ~p;
                if(c[0] == false_lit)
                    c[0] = c[1], c[1] = false_lit;
                assert(c[1] == false_lit);
                i++;

                // If 0th watch is true, then clause is already satisfied.
                Lit first = c[0];
                Watcher w = Watcher(cr, first);
                if(first != blocker && value(first) == l_True){
                    *j++ = w;
                    continue;
                }

                // Look for new watch:
                for(int k = 2; k < c.size(); k++)
                    if(value(c[k]) != l_False){
                        c[1] = c[k];
                        c[k] = false_lit;
                        watches[~c[1]].push(w);
                        goto NextClause;
                    }

                // Did not find watch -- clause is unit under assignment:
                *j++ = w;
                if(value(first) == l_False){
                    confl = cr;
                    qhead = trail.size();
                    // Copy the remaining watches:
                    while(i < end)
                        *j++ = *i++;
                }else
                    uncheckedEnqueue(first, cr);

                NextClause:;
            }
            ws.shrink(i - j);
        }

        if(initialPropagate && decisionLevel() == 0 && propagate_theories){
            assert(decisionLevel() == 0);
            //propagate any as yet unpropagated literals to each theory

            //Have to check _all_ the theories, even if we haven't eneueued of their literals, in case they have constants to propagate up.
            for(int theoryID = 0; theoryID < theories.size(); theoryID++){
                needsPropagation(theoryID);
            }
            initialPropagate = false;
        }

        //propagate theories;
        while(qhead == trail.size() && confl == CRef_Undef &&
              ((propagate_theories && theory_queue.size()) || unskippable_theory_q.size())){

            while(unskippable_theory_q.size() && (qhead == trail.size()) && confl == CRef_Undef){
                int theoryID = unskippable_theory_q.last();
                if(!propagateTheorySolver(theoryID, confl, theory_conflict)){
                    return confl;
                }else{

                    //only remove theory from propagation queue if it does not conflict
                    //there is a complication here, which is that in certain cases a new theory id may have been pushed into the queue
                    //during theory propagation.
                    assert(unskippable_theory_q.has(theoryID));
                    if(unskippable_theory_q.last() == theoryID){
                        unskippable_theory_q.pop();
                    }else{
                        unskippable_theory_q.remove(theoryID);
                    }
                    assert(!unskippable_theory_q.contains(theoryID));
                }
            }
            while(propagate_theories && theory_queue.size() && (qhead == trail.size())
                  && confl == CRef_Undef){
                int theoryID = theory_queue.last();
                if(!propagateTheorySolver(theoryID, confl, theory_conflict)){
                    return confl;
                }else{
                    //only remove theory from propagation queue if it does not conflict
                    //there is a complication here, which is that in certain cases a new theory id may have been pushed into the queue
                    //during theory propagation.
                    assert(in_theory_queue[theoryID]);

                    if(theory_queue.last() == theoryID){
                        theory_queue.pop();
                    }else{
                        theory_queue.remove(theoryID);
                    }
                    assert(!theory_queue.contains(theoryID));
                    in_theory_queue[theoryID] = false;
                }
            }
        }

    }while(qhead < trail.size());
    propagations += num_props;
    simpDB_props -= num_props;

    assert(confl != CRef_Undef || unskippable_theory_q.size() == 0);
    //one or more theories was not propagated (this is allowed in some settings)
    if(confl == CRef_Undef && theory_queue.size() && decisionLevel() > 0){
        for(int theoryID:theory_queue){
            if(theory_init_prop_trail_pos[theoryID] < 0){
                theory_init_prop_trail_pos[theoryID] = initial_qhead;
            }else{
                assert(theory_init_prop_trail_pos[theoryID] <= initial_qhead);
            }
        }
    }
    return confl;
}

bool Solver::propagateTheorySolver(int theoryID, CRef& confl, vec<Lit>& theory_conflict){
    double start_t = rtime(1);
    theory_conflict.clear();
    //todo: ensure that the bv theory comes first, as otherwise dependent theories may have to be propagated twice...


    if(!theories[theoryID]->propagateTheory(theory_conflict)){
        bool has_conflict = true;
#ifdef DEBUG_CORE
        for(Lit l:theory_conflict)
            assert(value(l) != l_Undef);
#endif
        if(has_conflict && !addConflictClause(theory_conflict, confl)){
            conflicting_heuristic = theories[theoryID]->getConflictingHeuristic();
            if(conflicting_heuristic){
                assert(conflicting_heuristic->getHeuristicIndex() > 0);
            }
            qhead = trail.size();
            stats_theory_conflicts++;
            stats_theory_conflict_time += (rtime(1) - start_t);
            return false;
        }
    }
    stats_theory_prop_time += (rtime(1) - start_t);
    return true;
}

/*_________________________________________________________________________________________________
 |
 |  reduceDB : ()  ->  [void]
 |
 |  Description:
 |    Remove half of the learnt clauses, minus the clauses locked by the current assignment. Locked
 |    clauses are clauses that are reason to some assignment. Binary clauses are never removed.
 |________________________________________________________________________________________________@*/
struct reduceDB_lt {
    ClauseAllocator& ca;

    reduceDB_lt(ClauseAllocator& ca_) :
            ca(ca_){
    }

    bool operator()(CRef x, CRef y){
        return ca[x].size() > 2 && (ca[y].size() == 2 || ca[x].activity() < ca[y].activity());
    }
};

void Solver::reduceDB(){
    int i, j;
    double extra_lim = cla_inc / learnts.size();    // Remove any clause below this activity

    sort(learnts, reduceDB_lt(ca));
    // Don't delete binary or locked clauses. From the rest, delete clauses from the first half
    // and clauses with activity smaller than 'extra_lim':
    for(i = j = 0; i < learnts.size(); i++){
        Clause& c = ca[learnts[i]];
        if(c.size() > 2 && !locked(c) && (i < learnts.size() / 2 || c.activity() < extra_lim)){
            stats_removed_clauses++;
            removeClause(learnts[i]);
        }else
            learnts[j++] = learnts[i];
    }
    learnts.shrink(i - j);
    checkGarbage();
}

void Solver::removeSatisfied(vec<CRef>& cs){
    int i, j;
    for(i = j = 0; i < cs.size(); i++){
        Clause& c = ca[cs[i]];
        if(satisfied(c))
            removeClause(cs[i]);
        else
            cs[j++] = cs[i];
    }
    cs.shrink(i - j);
}

void Solver::rebuildOrderHeap(){
    vec<Var> vs;
    for(Var v = 0; v < nVars(); v++)
        if(decision[v] && value(v) == l_Undef)
            vs.push(v);
    order_heap.build(vs);
}

void Solver::rebuildTheoryOrderHeap(){

    theory_order_heap.clear();
    for(int i = 0; i < decision_heuristics.size(); i++){
        theory_order_heap.insert(decision_heuristics[i]);
    }
}

void Solver::detectPureTheoryLiterals(){
    //Detect any theory literals that are occur only in positive or negative (or neither) polarity
    if(opt_detect_pure_lits){
        stats_pure_lits = 0;
        pure_literal_detections++;
        double startTime = rtime(1);
        lit_counts.growTo(nVars() * 2);
        for(int i = 0; i < lit_counts.size(); i++){
            lit_counts[i].seen = false;
        }
        assert(decisionLevel() == 0);

        //instead of counting the number of occurrences, just check if there are any occurences
        for(Lit l : trail){
            lit_counts[toInt(l)].seen = true;

        }

        for(CRef cr : clauses){
            Clause& c = ca[cr];
            if(!c.derivedClause()){
                for(Lit l : c){
                    lit_counts[toInt(l)].seen = true;
                }
            }
        }
        //the learnt clauses likely don't need to be counted here...
        for(Var v = 0; v < nVars(); v++){

            Lit l = mkLit(v, false);
            if(lit_counts[toInt(l)].seen && !lit_counts[toInt(l)].occurs){
                //stats_pure_lits--;
                lit_counts[toInt(l)].occurs = true;
                for(int n = 0; n < getNTheories(v); n++){
                    int theoryID = getTheoryID(v, n);
                    Lit theoryLit = getTheoryLit(l, n);

                    stats_pure_theory_lits--;
                    if(value(l) == l_Undef){
                        theories[theoryID]->setLiteralOccurs(~theoryLit, true);
                    }
                }
            }
            if(lit_counts[toInt(~l)].seen && !lit_counts[toInt(~l)].occurs){
                //stats_pure_lits--;
                lit_counts[toInt(~l)].occurs = true;

                for(int n = 0; n < getNTheories(v); n++){
                    int theoryID = getTheoryID(v, n);
                    Lit theoryLit = getTheoryLit(l, n);
                    stats_pure_theory_lits--;
                    if(value(l) == l_Undef){
                        theories[theoryID]->setLiteralOccurs(theoryLit, true);
                    }
                }
            }

            if(lit_counts[toInt(l)].occurs && !lit_counts[toInt(l)].seen){
                lit_counts[toInt(l)].occurs = false;

                //stats_pure_lits++;
                if(hasTheory(v)){
                    for(int n = 0; n < getNTheories(v); n++){
                        int theoryID = getTheoryID(v, n);
                        Lit theoryLit = getTheoryLit(l, n);
                        stats_pure_theory_lits++;
                        if(value(l) == l_Undef){
                            theories[theoryID]->setLiteralOccurs(~theoryLit, false);
                        }
                        //setPolarity(v,false);
                    }
                }else{
                    //we can safely assign this now...
                    if(!lit_counts[toInt(~l)].seen){
                        /*	if(decision[v])
						 setDecisionVar(v,false);*/
                    }else{
                        if(value(l) == l_Undef){
                            //uncheckedEnqueue(~l);
                        }
                    }
                }
            }
            if(lit_counts[toInt(~l)].occurs && !lit_counts[toInt(~l)].seen){
                lit_counts[toInt(~l)].occurs = false;

                if(hasTheory(v)){
                    stats_pure_theory_lits++;
                    for(int n = 0; n < getNTheories(v); n++){
                        int theoryID = getTheoryID(v, n);
                        Lit theoryLit = getTheoryLit(l, n);
                        if(value(l) == l_Undef){
                            //level 0 variables are already handled through a separate mechanism in the theory solvers...
                            theories[theoryID]->setLiteralOccurs(theoryLit, false);
                        }
                    }
                    if(!lit_counts[toInt(l)].occurs){
                        //setDecisionVar(v,false);//If v is pure in both polarities, and is a theory var, then just don't assign it at all - it is unconstrained.
                        //This _should_ be safe to do, if the theory semantics make sense... although there are some clause learning schemes that introduce previosuly unused literals that might break with this...
                    }
                }else{
                    //we can safely assign this now...
                    if(!lit_counts[toInt(l)].seen){
                        /*		if(decision[v])
						 setDecisionVar(v,false);*/
                    }else if(value(l) == l_Undef){
                        //uncheckedEnqueue(l);
                    }
                }
            }

        }
        assert(stats_pure_lits <= nVars() * 2);
        assert(stats_pure_lits >= 0);
        stats_pure_lit_time += rtime(1) - startTime;

    }
}

/*_________________________________________________________________________________________________
 |
 |  simplify : [void]  ->  [bool]
 |
 |  Description:
 |    Simplify the clause database according to the current top-level assigment. Currently, the only
 |    thing done here is the removal of satisfied clauses, but more things can be put here.
 |________________________________________________________________________________________________@*/
bool Solver::simplify(){
    assert(decisionLevel() == 0);

    //it is important (for efficiency reasons, not correctness)
    //to run this pure literal detection before propagate, if new theory atoms have been added since the last call to propagate.
    if(ok && initialPropagate && opt_detect_pure_lits){
        detectPureTheoryLiterals();
    }

    bool propagate_theories = !disable_theories;// && opt_propagate_theories_during_fast_simplification; //it is not safe to disable theory props at level 0

    if(!ok || propagate(propagate_theories) != CRef_Undef ||
       !ok) //yes, the second ok check is now required, because propagation of a theory can make the solver unsat at this point...
        return ok = false;

    if(nAssigns() == simpDB_assigns || (simpDB_props >
                                        0)) //if pure literal detection is enabled, then simplify MUST be called at least once after new theory atoms are added.
        return true;

    // Remove satisfied clauses:
    removeSatisfied(learnts);
    if(remove_satisfied){        // Can be turned off.
        removeSatisfied(clauses);

        // Remove all released variables from the trail:
        for(int i = 0; i < released_vars.size(); i++){
            assert(seen[released_vars[i]] == 0);
            seen[released_vars[i]] = 1;
        }

        int i, j;
        for(i = j = 0; i < trail.size(); i++)
            if(seen[var(trail[i])] == 0)
                trail[j++] = trail[i];
        trail.shrink(i - j);
        qhead = trail.size();

        for(int i = 0; i < released_vars.size(); i++)
            seen[released_vars[i]] = 0;

        // Released variables are now ready to be reused:
        append(released_vars, free_vars);
        released_vars.clear();
    }
    checkGarbage();
    rebuildOrderHeap();

    simpDB_assigns = nAssigns();
    simpDB_props = clauses_literals + learnts_literals;   // (shouldn't depend on stats really, but it will do for now)

    detectPureTheoryLiterals();

    return true;
}

void Solver::addClauseSafely(vec<Lit>& ps){
    if(opt_write_learnt_clauses){
        fprintf(opt_write_learnt_clauses, "learnt fact ");
        for(Lit l:ps){
            fprintf(opt_write_learnt_clauses, " %d", dimacs(unmap(l)));
        }
        fprintf(opt_write_learnt_clauses, " 0\n");
        fflush(opt_write_learnt_clauses);
    }

    if(decisionLevel() == 0){
        addClause_(ps, true);
    }else{

        //check if this clause should imply a literal right now.
        //and IF that literal is already implied, then we may still need to adjust the level the literal should be implied at.
        int undef_count = 0;
        int false_count = 0;
        sort(ps);
        Lit p;
        int i, j;
        for(i = j = 0, p = lit_Undef; i < ps.size(); i++){
            assert(var(ps[i]) != var(theoryDecision));
            if(((value(ps[i]) == l_True && level(var(ps[i])) == 0)) || ps[i] == ~p)
                return;
            else if((value(ps[i]) != l_False || level(var(ps[i])) != 0) && ps[i] != p){
                ps[j++] = p = ps[i];
                if(value(ps[i]) == l_Undef)
                    undef_count++;
                else if(value(ps[i]) == l_False){
                    false_count++;
                }
            }
        }
        ps.shrink(i - j);
        if(false_count == ps.size() - 1){
            //this clause is unit under the current assignment.
            //although we _could_ wait until a restart to add this clause, in many cases this will lead to very poor solver behaviour.
            //so we are going to propagate this clause now.

            if(ps.size() == 0){
                ok = false;
                return;
            }else if(ps.size() == 1){
                enqueueLazy(ps[0], 0, CRef_Undef);
            }else{

                bool conflicting = true;
                int nfalse = 0;
                int max_lev = 0;
                int notFalsePos1 = -1;
                int notFalsePos2 = -1;
                for(int j = 0; j < ps.size(); j++){
                    assert(var(ps[j]) < nVars());
                    if(value(ps[j]) == l_False){
                        nfalse++;
                    }else{
                        conflicting = false;
                        if(notFalsePos1 < 0)
                            notFalsePos1 = j;
                        else if(notFalsePos2 < 0){
                            notFalsePos2 = j;
                        }
                    }
                    if(value(ps[j]) != l_Undef){
                        int l = level(var(ps[j]));
                        if(l > max_lev){
                            max_lev = l;
                        }
                    }
                }
                assert(!conflicting);
                if(!conflicting){
                    assert(notFalsePos1 >= 0);
                    if(notFalsePos1 >= 0 && notFalsePos2 >= 0){
                        assert(notFalsePos1 != notFalsePos2);
                        if(notFalsePos1 == 1){
                            std::swap(ps[0], ps[notFalsePos2]);
                        }else{
                            std::swap(ps[0], ps[notFalsePos1]);
                            std::swap(ps[1], ps[notFalsePos2]);
                        }
                    }else{
                        std::swap(ps[0], ps[notFalsePos1]);
                    }
                    assert(value(ps[0]) != l_False);
                    if(notFalsePos2 >= 0){
                        assert(value(ps[1]) != l_False);
                    }
                }

                CRef cr = ca.alloc(ps);
                ca[cr].setDerived(true);
                clauses.push(cr);
                attachClause(cr);
                enqueueLazy(ps[0], max_lev, cr);
            }
        }else{
            //add this clause at the next solver restart
            clauses_to_add.push();
            ps.copyTo(clauses_to_add.last());
        }


    }
}

bool Solver::addConflictClause(vec<Lit>& ps, CRef& confl_out, bool permanent){
    if(opt_write_learnt_clauses){
        fprintf(opt_write_learnt_clauses, "learnt ");
        for(Lit l:ps){
            fprintf(opt_write_learnt_clauses, " %d", dimacs(unmap(l)));
        }
        fprintf(opt_write_learnt_clauses, " 0\n");
        fflush(opt_write_learnt_clauses);
    }

    sort(ps);
    Lit p;
    int i, j;
    for(i = j = 0, p = lit_Undef; i < ps.size(); i++){
        assert(var(ps[i]) != var(theoryDecision));

        if(((value(ps[i]) == l_True && level(var(ps[i])) == 0)) || ps[i] == ~p)
            return true;
        else if((value(ps[i]) != l_False || level(var(ps[i])) != 0) && ps[i] != p)
            ps[j++] = p = ps[i];

    }
    ps.shrink(i - j);
    confl_out = CRef_Undef;
    if(ps.size() == 0){
        ok = false;
        cancelUntil(0);
        return false;
    }else if(ps.size() == 1){
        cancelUntil(0);

        assert(var(ps[0]) < nVars());

        if(!enqueue(ps[0])){
            ok = false;
            return false;
        }

    }else{
        //find the highest level in the conflict (should be the current decision level, but we won't require that)
        bool conflicting = true;
        int nfalse = 0;
        int max_lev = 0;
        bool satisfied = false;
        int notFalsePos1 = -1;
        int notFalsePos2 = -1;
        for(int j = 0; j < ps.size(); j++){
            assert(var(ps[j]) < nVars());
            //assert(value(ps[j])==l_False);
            if(value(ps[j]) == l_False){
                nfalse++;
            }else{
                conflicting = false;
                if(value(ps[j]) == l_True)
                    satisfied = true;
                if(notFalsePos1 < 0)
                    notFalsePos1 = j;
                else if(notFalsePos2 < 0){
                    notFalsePos2 = j;
                }
            }
            if(value(ps[j]) != l_Undef){
                int l = level(var(ps[j]));
                if(l > max_lev){
                    max_lev = l;
                }
            }
        }
        if(!conflicting){
            assert(notFalsePos1 >= 0);
            if(notFalsePos1 >= 0 && notFalsePos2 >= 0){
                assert(notFalsePos1 != notFalsePos2);
                if(notFalsePos1 == 1){
                    std::swap(ps[0], ps[notFalsePos2]);
                }else{
                    std::swap(ps[0], ps[notFalsePos1]);
                    std::swap(ps[1], ps[notFalsePos2]);
                }
            }else{
                std::swap(ps[0], ps[notFalsePos1]);
            }
            assert(value(ps[0]) != l_False);
            if(notFalsePos2 >= 0){
                assert(value(ps[1]) != l_False);
            }
        }

#ifdef DEBUG_CORE
        if(conflicting){
            for(int j = 0; j < ps.size(); j++)
                assert(value(ps[j]) == l_False);
        }
#endif
        if(!permanent && ps.size() > opt_temporary_theory_conflicts){
            if(tmp_clause == CRef_Undef){
                tmp_clause = ca.alloc(ps, false);
                tmp_clause_sz = ps.size();
            }else if(tmp_clause_sz < ps.size()){
                ca[tmp_clause].mark(1);                            //is this needed?
                ca.free(tmp_clause);
                tmp_clause = ca.alloc(ps, false);
                tmp_clause_sz = ps.size();
            }else{
                Clause& c = ca[tmp_clause];
                assert(tmp_clause_sz >= ps.size());
                assert(tmp_clause_sz >= c.size());
                c.grow(tmp_clause_sz - c.size());
                for(int i = 0; i < ps.size(); i++){
                    c[i] = ps[i];
                }
                c.shrink(c.size() - ps.size());
            }

            confl_out = tmp_clause;
        }else{

            CRef cr = ca.alloc(ps, !permanent && !opt_permanent_theory_conflicts);
            ca[cr].setDerived(true);
            if(permanent || opt_permanent_theory_conflicts)
                clauses.push(cr);
            else{
                learnts.push(cr);
                if(--learntsize_adjust_cnt <= 0){
                    learntsize_adjust_confl *= learntsize_adjust_inc;
                    learntsize_adjust_cnt = (int) learntsize_adjust_confl;
                    max_learnts *= learntsize_inc;

                    if(verbosity >= 1)
                        printf("| %9d | %7d %8d %8d | %8d %8d %6.0f | %" PRId64 " removed |\n", (int) conflicts,
                               (int) dec_vars - (trail_lim.size() == 0 ? trail.size() : trail_lim[0]), nClauses(),
                               (int) clauses_literals, (int) max_learnts, nLearnts(),
                               (double) learnts_literals / nLearnts(), stats_removed_clauses);
                }

            }

            attachClause(cr);
            confl_out = cr;
        }

        if(!satisfied && nfalse == ps.size() - 1){
            cancelUntil(max_lev);
            assert(!conflicting);
            assert(value(ps[0]) != l_False);
            assert(value(ps[1]) == l_False);
            if(value(ps[0]) == l_Undef){
                uncheckedEnqueue(ps[0], confl_out);
            }
        }
        if(!conflicting){
            confl_out = CRef_Undef;
        }
        return !conflicting;
    }
    return true;
}

bool Solver::addDelayedClauses(CRef& conflict_out){
    conflict_out = CRef_Undef;
    while(clauses_to_add.size() && ok){
        if(!addConflictClause(clauses_to_add.last(), conflict_out, true)){
            clauses_to_add.pop();
            return false;
        }
        clauses_to_add.pop();
    }
    return true;
}

/*_________________________________________________________________________________________________
 |
 |  search : (nof_conflicts : int) (params : const SearchParams&)  ->  [lbool]
 |
 |  Description:
 |    Search for a model the specified number of conflicts.
 |    NOTE! Use negative value for 'nof_conflicts' indicate infinity.
 |
 |  Output:
 |    'l_True' if a partial assigment that is consistent with respect to the clauseset is found. If
 |    all variables are decision variables, this means that the clause set is satisfiable. 'l_False'
 |    if the clause set is unsatisfiable. 'l_Undef' if the bound on number of conflicts is reached.
 |________________________________________________________________________________________________@*/
lbool Solver::search(int nof_conflicts){
    assert(ok);
    int backtrack_level;
    int conflictC = 0;
    vec<Lit> learnt_clause;
    Heuristic* previous_conflict_heuristic = nullptr;

    Heuristic* last_decision_heuristic = nullptr;
    bool decision_heuristic_changed = false;
    starts++;
    static int decision_iter = 0;
    bool using_theory_decisions = opt_decide_theories && drand(random_seed) < opt_random_theory_freq;
    bool using_theory_order_heap = opt_decide_theories && (opt_theory_order_vsids || opt_theory_order_swapping) &&
                                   drand(random_seed) < opt_random_theory_order_freq;

    bool propagate_theories_during_assumptions = opt_theory_propagate_assumptions;

    if(decisionLevel() == 0 && initialPropagate && opt_detect_pure_lits && !simplify()){
        return l_False;//if using pure literal detection, and the theories haven't been propagated yet, run simpify
    }
    bool last_propagation_was_conflict = false;
    CRef confl = CRef_Undef;
    n_theory_decision_rounds += using_theory_decisions;
    for(;;){
        propagate:

        bool all_assumptions_assigned = decisionLevel() >= assumptions.size();
        bool propagate_theories = (!disable_theories) &&
                                  (propagate_theories_during_assumptions || decisionLevel() == 0 ||
                                   all_assumptions_assigned);
        if(opt_decide_theories_only_prop_decision && propagate_theories && decisionLevel() > 0 &&
           !last_propagation_was_conflict && last_decision_heuristic && !decision_heuristic_changed){
            propagate_theories = false;
        }
        if(!propagate_theories){
            stats_skipped_theory_prop_rounds++;
        }

        confl = propagate(propagate_theories);

        conflict:

        if(!okay() || (confl != CRef_Undef)){
            // CONFLICT
            conflicts++;
            conflictC++;

            if(last_decision_heuristic){
                n_theory_conflicts++;
                consecutive_theory_conflicts++;
                if(opt_theory_conflict_max && consecutive_theory_conflicts >= opt_theory_conflict_max){
                    next_theory_decision = conflicts + consecutive_theory_conflicts;
                    consecutive_theory_conflicts = 0;
                }
            }else{
                consecutive_theory_conflicts = 0;
            }
            if(decisionLevel() == 0)
                return l_False;
            learnt_clause.clear();
            analyze(confl, learnt_clause, backtrack_level);

            int lowest_conflicting_decision_level = decisionLevel();
            if(last_decision_heuristic && (!conflicting_heuristic || conflicting_heuristic->getPriority() <
                                                                     last_decision_heuristic->getPriority())){
                conflicting_heuristic = last_decision_heuristic;
            }
            bool multiple_involved_theories = false;
            bool order_changed = false;
            if(learnt_clause.size() > 1 && last_decision_heuristic && opt_theory_order_swapping){
                //nadel-style theory order swapping
                if(opt_verb >= 3){
                    printf("Last decision: %d, conflicting decision %d\n", last_decision_heuristic->getHeuristicIndex(),
                           conflicting_heuristic->getHeuristicIndex());
                    printf("Old decision order: ");
                    for(Heuristic* h:decision_heuristics){
                        printf("%d, ", h->getHeuristicIndex());
                    }
                    printf("\n");
                }

                assert(decision_heuristics.contains(conflicting_heuristic));
                int start_size = decision_heuristics.size();
                //check to see if the analyzed conflict includes any literals
                //that were either decided or forced by a theory solver
                swapping_involved_theories.clear();
                swapping_uninvolved_pre_theories.clear();
                swapping_uninvolved_post_theories.clear();
                swapping_involved_theory_order.clear();
                //vec<int> involved_theories;
                // if(!opt_theory_order_swapping_last_only) {
                // swapping_involved_theories.insert(last_decision_heuristic->getHeuristicIndex());
                //}
                //note: the last decision heuristic can in some cases be the conflicting heuristic, so these can be the same
                if(opt_theory_order_swapping_prioritize_last_decision &&
                   conflicting_heuristic != last_decision_heuristic)
                    swapping_involved_theories.insert(last_decision_heuristic->getHeuristicIndex());
                swapping_involved_theories.insert(conflicting_heuristic->getHeuristicIndex());

                int max_involved = opt_theory_order_swapping_max_invovled;
                assert(max_involved >= 2);
                if(opt_theory_order_swapping_reset_counts_new_conflict){
                    //find all the involved theories
                    analyzeHeuristicDecisions(confl, swapping_involved_theories, -1,
                                              conflicting_heuristic->getPriority());
                    //check to see if the previous conflict theory is involved in the current conflict
                    if(previous_conflict_heuristic != nullptr &&
                       !swapping_involved_theories.contains(previous_conflict_heuristic->getHeuristicIndex())){
                        //this is a new conflict net, that was not previously involved
                        if(heuristic_swapping_restarts > stats_max_swap_count)
                            stats_max_swap_count = heuristic_swapping_restarts;

                        heuristic_swapping_restarts = 0;//reset the luby series to 0
                        stats_swapping_resets++;

                    }
                }
                stats_swapping_conflict_count++;
                if(opt_theory_order_swapping_luby){
                    max_involved *= luby(restart_inc, heuristic_swapping_restarts++);
                }
                assert(max_involved >= 2);
                if(!opt_theory_order_swapping_reset_counts_new_conflict){
                    analyzeHeuristicDecisions(confl, swapping_involved_theories, max_involved,
                                              conflicting_heuristic->getPriority());
                }else{
                    while(swapping_involved_theories.size() > max_involved){
                        swapping_involved_theories.pop();
                    }
                }
                previous_conflict_heuristic = conflicting_heuristic;
                assert(swapping_involved_theories.has(conflicting_heuristic->getHeuristicIndex()));
                assert(swapping_involved_theories.size() >= 1);
                assert(swapping_involved_theories.size() <= max_involved);
                multiple_involved_theories = swapping_involved_theories.size() > 1;
                if(swapping_involved_theories.size() > 1 || opt_theory_order_swapping_first_on_unit){
                    order_changed = true;
                    //things to try:
                    //1) disable normal restarts in the solver
                    //2) reverse or randomize swapping_involved_theory_order
                    //3) reset the luby series if the _entire_ conflict net set (not just the part up to max involved) doesn't include the previous conflict net

                    //if opt_theory_order_swapping_first_on_unit is set, then on a unit conflict, move the conflict heuristic to the begining (without changing the rest of the decision order).
                    bool unit = (swapping_involved_theories.size() == 1 && opt_theory_order_swapping_first_on_unit);
                    int lowest_involved_position = decision_heuristics.size();

                    if(opt_theory_order_conflict_skip_middle && swapping_involved_theories.size() > 2){
                        //place the conflicting heuristic before the earliest involved heuristic, but leave all the
                        //remaining heuristics in the current order.
                        for(int i = 0; i < decision_heuristics.size(); i++){
                            Heuristic* t = decision_heuristics[i];
                            if(t != conflicting_heuristic && swapping_involved_theories.has(t->getHeuristicIndex())){
                                swapping_involved_theories.clear();
                                swapping_involved_theories.insert(t->getHeuristicIndex());
                                swapping_involved_theories.insert(conflicting_heuristic->getHeuristicIndex());
                                break;
                            }
                        }
                        assert(swapping_involved_theories.size() == 2);
                    }

                    for(int i = 0; i < decision_heuristics.size(); i++){
                        Heuristic* t = decision_heuristics[i];
                        if(t == conflicting_heuristic){

                            if(i < lowest_involved_position){
                                lowest_involved_position = i;
                            }
                            int conflict_level = first_heuristic_decision_level[t->getHeuristicIndex()];
                            if(conflict_level >= 0 && conflict_level < lowest_conflicting_decision_level){
                                lowest_conflicting_decision_level = conflict_level;//use this to backtrack past the earliest decision made by one of the conflicting heuristics.
                            }
                        }else if(swapping_involved_theories.has(t->getHeuristicIndex())){
                            swapping_involved_theory_order.push(t);
                            if(i < lowest_involved_position){
                                lowest_involved_position = i;
                            }
                            int conflict_level = first_heuristic_decision_level[t->getHeuristicIndex()];
                            if(conflict_level >= 0 && conflict_level < lowest_conflicting_decision_level){
                                lowest_conflicting_decision_level = conflict_level;//use this to backtrack past the earliest decision made by one of the conflicting heuristics.
                            }
                        }else if(lowest_involved_position == decision_heuristics.size() && !unit){
                            swapping_uninvolved_pre_theories.push(t);
                        }else{
                            swapping_uninvolved_post_theories.push(t);
                        }
                    }

                    if(opt_verb >= 3){
                        printf("Uninvolved pre theories: ");
                        for(Heuristic* h:swapping_uninvolved_pre_theories){
                            printf("%d, ", h->getHeuristicIndex());
                        }
                        printf("\nInvolved theories: ");
                        printf("%d, ", conflicting_heuristic->getHeuristicIndex());
                        for(Heuristic* h:swapping_involved_theory_order){
                            printf("%d, ", h->getHeuristicIndex());
                        }
                        printf("\nUninvolved post theories: ");
                        for(Heuristic* h:swapping_uninvolved_post_theories){
                            printf("%d, ", h->getHeuristicIndex());
                        }
                        printf("\nConflcit level %d, lowest involved theory level %d\n", backtrack_level,
                               lowest_conflicting_decision_level);
                    }
                    decision_heuristics.clear();
                    swapping_uninvolved_pre_theories.copyTo(decision_heuristics);
                    decision_heuristics.push(conflicting_heuristic);
                    decision_heuristics.extend(swapping_involved_theory_order);
                    decision_heuristics.extend(swapping_uninvolved_post_theories);
                    for(int i = 0; i < decision_heuristics.size(); i++){
                        decision_heuristics[i]->setHeuristicOrder(i);
                    }
                    if(opt_verb >= 3){
                        printf("New decision order: ");
                        for(Heuristic* h:decision_heuristics){
                            printf("%d, ", h->getHeuristicIndex());
                        }
                        printf("\n");
                    }
                }
                assert(decision_heuristics.size() == start_size);

            }
            int conflict_counter_exceeding_count = -1;
            if(conflicting_heuristic && (opt_theory_order_conflict_on_unit || multiple_involved_theories)){
                int64_t conflict_restart_limit = opt_theory_order_conflict_restart;
                if(opt_theory_order_conflict_luby){
                    conflict_restart_limit *= luby(restart_inc, heuristic_conflict_restarts);
                }
                //if multiple theories exceed the conflict count limit, favor the conflicting heuristic
                assert(theory_conflict_counters.size() > conflicting_heuristic->getHeuristicIndex());
                if(++theory_conflict_counters[conflicting_heuristic->getHeuristicIndex()] >= conflict_restart_limit){
                    conflict_counter_exceeding_count = conflicting_heuristic->getHeuristicIndex();
                }

                if(opt_theory_order_conflict_count_analysis > 0){
                    counter_involved_theories.clear();
                    analyzeHeuristicDecisions(confl, counter_involved_theories, -1,
                                              conflicting_heuristic->getPriority());
                    assert(swapping_involved_theories.has(conflicting_heuristic->getHeuristicIndex()));
                    assert(swapping_involved_theories.size() >= 1);
                    for(int i:swapping_involved_theories){
                        if(i != conflicting_heuristic->getHeuristicIndex()){
                            auto count = ++theory_conflict_counters[i];
                            if(opt_theory_order_conflict_count_analysis == 2 && count >= conflict_restart_limit &&
                               conflict_counter_exceeding_count < 0){
                                conflict_counter_exceeding_count = i;
                            }
                        }
                    }

                }
            }

            cancelUntil(backtrack_level);

            if(opt_theory_order_swapping && order_changed){
                //rebuild the theory order queue
                theory_order_heap.build(decision_heuristics);

            }

            //this is now slightly more complicated, if there are multiple lits implied by the super solver in the current decision level:
            //The learnt clause may not be asserting.

            if(learnt_clause.size() == 1){
                uncheckedEnqueue(learnt_clause[0]);
            }else{
                CRef cr = ca.alloc(learnt_clause, true);
                learnts.push(cr);
                attachClause(cr);
                claBumpActivity(ca[cr]);

                if(value(learnt_clause[0]) == l_Undef){
                    uncheckedEnqueue(learnt_clause[0], cr);
                }else{

                    assert(S);
                    if(!S){
                        throw std::runtime_error("Critical error: bad learnt clause.");
                        //throw std::runtime_error("Critical error: bad learnt clause. Aborting\n");
                    }
                    //this is _not_ an asserting clause, its a conflict that must be passed up to the super solver.
                    analyzeFinal(cr, lit_Undef, conflict);

                    if(conflicting_heuristic){
                        theoryBumpActivity(conflicting_heuristic);
                        theoryDecayActivity();
                    }else if(opt_vsids_solver_as_theory){
                        theoryBumpActivity(decisionTheory);
                        theoryDecayActivity();
                    }

                    varDecayActivity();
                    claDecayActivity();
                    conflicting_heuristic = nullptr;
                    return l_False;
                }
            }
            if(conflicting_heuristic){
                theoryBumpActivity(conflicting_heuristic);
                theoryDecayActivity();

            }else if(opt_vsids_solver_as_theory){
                theoryBumpActivity(decisionTheory);
                theoryDecayActivity();
            }
            varDecayActivity();
            claDecayActivity();

            if(conflicting_heuristic && (opt_theory_order_conflict_on_unit || multiple_involved_theories)){
                assert(theory_conflict_counters.size() > conflicting_heuristic->getHeuristicIndex());


                if(opt_theory_order_conflict_restart > 0 &&
                   (!opt_theory_order_conflict_conservative || decision_heuristics[0] != conflicting_heuristic) &&
                   conflict_counter_exceeding_count > -1){
                    //restart the solver
                    heuristic_conflict_restarts++;
                    cancelUntil(0);
                    stats_theory_conflict_counter_restarts++;
                    if(opt_theory_order_conflict_sort_counter){
                        //sort the heuristics by their conflict counter value, putting the highest value first
                        sort(decision_heuristics, [&](Heuristic* a, Heuristic* b) -> bool{
                            return theory_conflict_counters[a->getHeuristicIndex()] >
                                   theory_conflict_counters[b->getHeuristicIndex()];
                        });
                        //assert(decision_heuristics[0]==conflicting_heuristic);
                    }else{

                        if(opt_theory_order_conflict_sort_vsids){
                            //before moving the highest conflict counter to the beginning, sort all heuristics by their activity/priority
                            sort(decision_heuristics, HeuristicActivityOrderLt());
                        }

                        //reorder the conflicting theories to put the decision theory first.
                        int theory_conflict_pos = -1;
                        for(int i = 0; i < decision_heuristics.size(); i++){
                            if(decision_heuristics[i] == conflicting_heuristic){
                                theory_conflict_pos = i;
                                break;
                            }
                        }
                        assert(theory_conflict_pos >= 0);
                        //assert(decision_heuristics[theory_conflict_pos]==conflicting_heuristic);
                        if(theory_conflict_pos > 0){
                            if(opt_theory_order_swapping_preserve_order){
                                //Preserve the previous order of the heuristics, except for moving this one to the begining.
                                Heuristic* h = decision_heuristics[theory_conflict_pos];
                                for(int i = theory_conflict_pos; i >= 1; i--){
                                    assert(i >= 1);
                                    decision_heuristics[i] = decision_heuristics[i - 1];
                                }
                                decision_heuristics[0] = h;
                            }else{
                                std::swap(decision_heuristics[0], decision_heuristics[theory_conflict_pos]);
                            }
                        }
                        assert(decision_heuristics[0] == conflicting_heuristic);
                    }
                    assert(conflict_counter_exceeding_count >= 0);
                    if(opt_theory_order_conflict_count_preserve == 0){
                        //erase all conflict counts when on exceeds the max found
                        for(int i = 0; i < theory_conflict_counters.size(); i++){
                            theory_conflict_counters[i] = 0;
                        }
                    }else if(opt_theory_order_conflict_count_preserve == 1){
                        //only reset the counter of the conflict that triggered the reset
                        theory_conflict_counters[conflict_counter_exceeding_count] = 0;
                    }else if(opt_theory_order_conflict_count_preserve == 2){
                        //only reset the counters involved in the conflict
                        theory_conflict_counters[conflict_counter_exceeding_count] = 0;
                        for(int i:counter_involved_theories){
                            theory_conflict_counters[i] = 0;
                        }
                    }else if(opt_theory_order_conflict_count_preserve == 3){
                        //divide all conflict counter values by 2 on reset
                        for(int i = 0; i < theory_conflict_counters.size(); i++){
                            theory_conflict_counters[i] >>= 1;
                        }
                    }else if(opt_theory_order_conflict_count_preserve == 4){
                        //divide the involved conflict counter values by 2 on reset
                        theory_conflict_counters[conflict_counter_exceeding_count] >>= 1;
                        for(int i:counter_involved_theories){
                            theory_conflict_counters[i] >>= 1;
                        }
                    }
                    for(int i = 0; i < decision_heuristics.size(); i++){
                        decision_heuristics[i]->setHeuristicOrder(i);
                    }
                    theory_order_heap.build(decision_heuristics);

                }
            }

            if(lowest_conflicting_decision_level > 0 && lowest_conflicting_decision_level < decisionLevel()){
                cancelUntil(lowest_conflicting_decision_level - 1);
            }

            if(!propagate_theories_during_assumptions && decisionLevel() > 0 && decisionLevel() < assumptions.size()){
                //improve this in the future!
                //if we backtracked past the assumption level, AND if theory propagate was disabled while assigning assumptions
                //then some theories (notably, the bv theory) currently require us to re-propagate all theory literals of that bv,
                //including (surprisingly) theory literals enqueued at LOWER levels than backtrack_level. This is because those
                //theory literals did not have their effects propagated until after the backtrack level.
                //with better book-keeping, this could probably be avoided... however, it only occurs when backtracking past the
                //assumption level, but not to level 0, which in practice should (?) only happen when the assumptions are in conflict.

                //backtrack to 0, and disable the 'no theory propagation during assumptions'
                //optimization during the current search() call, as a conflict has occured during assumption
                //assignment itself.
                propagate_theories_during_assumptions = true;
                cancelUntil(0); //note: in cases where we have vast nubmers of assumptions,
                //and we expect to always have conflicts during the assumptions, this
                //solution makes the !opt_theory_propagate_assumptions optimization a bad idea as a whole.
                //it might be worth re-exploring the previous, much more complex solution of
                //only backtracking and reasserting the theory solvers in that case... but that strategy
                //had some non-obvious bugs and pitfalls, when combinging the bv theory with the graph theory.
            }

            if(--learntsize_adjust_cnt <= 0){
                learntsize_adjust_confl *= learntsize_adjust_inc;
                learntsize_adjust_cnt = (int) learntsize_adjust_confl;
                max_learnts *= learntsize_inc;

                if(verbosity >= 1)
                    printf("|c%9d | %7d %8d %8d | %8d %8d %6.0f | %" PRId64 " removed |\n", (int) conflicts,
                           (int) dec_vars - (trail_lim.size() == 0 ? trail.size() : trail_lim[0]), nClauses(),
                           (int) clauses_literals, (int) max_learnts, nLearnts(),
                           (double) learnts_literals / nLearnts(), stats_removed_clauses);
            }
            conflicting_heuristic = nullptr;
            last_propagation_was_conflict = true;
        }else{
            last_propagation_was_conflict = false;
            assert(theory_queue.size() == 0 || !propagate_theories);
            assert(!unskippable_theory_q.size());
            if(!addDelayedClauses(confl))
                goto conflict;

            if(opt_subsearch == 0 && decisionLevel() < initial_level){
                return l_Undef;            //give up if we have backtracked past the super solvers decisions
            }else if(opt_subsearch == 2 && S && confl == CRef_Undef){
                if(super_qhead < S->qhead && !propagateTheory(
                        theory_conflict)) //re-enqueue any super solver decisions that we have backtracked past, and keep going
                {
                    if(!addConflictClause(theory_conflict, confl))
                        return l_False;
                    goto conflict;
                }
            }


            // NO CONFLICT
            if((opt_restarts && (nof_conflicts >= 0 && conflictC >= nof_conflicts)) || !withinBudget()){
                // Reached bound on number of conflicts:
                progress_estimate = progressEstimate();
                cancelUntil(initial_level);
                return l_Undef;
            }

            // Simplify the set of problem clauses:
            if(decisionLevel() == 0 && !simplify())
                return l_False;

            if(learnts.size() - nAssigns() >= max_learnts)
                // Reduce the set of learnt clauses:
                reduceDB();
            Heuristic* next_decision_heuristic = nullptr;
            decision_heuristic_changed = false;

            bool assumps_processed = false;
            Lit next = lit_Undef;
            CRef decision_reason = CRef_Undef;  //If a theory makes a decision, then it can supply a 'reason' for that decision.
            //(That reason may be used by some heuristics, but will not be used by conflict analysis.)
            while(decisionLevel() < assumptions.size()){
                // Perform user provided assumption:
                Lit p = assumptions[decisionLevel()];
                assert(p != lit_Undef);
                assert(var(p) < nVars());
                if(value(p) == l_True){
                    // Dummy decision level:
                    newDecisionLevel();
                    assumps_processed = true;
                }else if(value(p) == l_False){
                    analyzeFinal(~p, conflict);
                    return l_False;
                }else{
                    next = p;
                    break;
                }
            }

            if(only_propagate_assumptions && next == lit_Undef && decisionLevel() == assumptions.size())
                return l_True;

            if(next == lit_Undef && assumps_processed){
                //allow theory propagation a chance to be applied to the assumptions, in the case that the last assumption
                //was a dummy assumption, and opt_theory_propagate_assumptions is false
                continue;
            }

            //Note: decision level is now added before theories make their decisions, to allow them to decide multiple literals at once.
            newDecisionLevel();

            /**
             * Give the theory solvers a chance to make decisions
             */
            if(opt_decide_theories && !disable_theories && using_theory_decisions && next == lit_Undef &&
               (opt_theory_conflict_max == 0 || conflicts >= next_theory_decision)){

                int next_var_priority = INT_MIN;

                //remove decided vars from order heap
                while(!order_heap.empty() && value(order_heap.peekMin()) != l_Undef){
                    order_heap.removeMin();
                }
                if(!order_heap.empty()){
                    next_var_priority = priority[order_heap.peekMin()];
                    //printf("Priority var is %d; with priority %d\n",order_heap.peekMin(), next_var_priority);
                    //printf("Vars %d,%d have assignment %d,%d\n", 91569,91580,value(91569),value(91580));
                }

                if(((opt_theory_order_vsids || opt_theory_order_swapping) && using_theory_order_heap)){
                    while(next == lit_Undef && !theory_order_heap.empty() &&
                          theory_order_heap.peekMin()->getPriority() >= next_var_priority){
                        Heuristic* h = theory_order_heap.peekMin();
                        ++decision_iter;
                        //int theoryID = h->getTheoryIndex();
                        if(!heuristicSatisfied(h)){
                            if(opt_vsids_both &&
                               next_var_priority == theory_order_heap.peekMin()->getPriority()){
                                //give the main solver a chance to make a decision, if it has a variable with higher activity
                                if(!order_heap.empty()){
                                    Var v = order_heap.peekMin();
                                    assert(value(v) ==
                                           l_Undef);//because assigned lits should have been removed from queue above.
                                    if(value(v) == l_Undef && activity[v] > h->getActivity()){
                                        stats_solver_preempted_decisions++;
                                        break;
                                    }
                                }

                            }
                            next = h->decideTheory(decision_reason);
                            if(next != lit_Undef){
                                next_decision_heuristic = h;
#ifdef DEBUG_CORE
                                /*     for (int i = 0; i < decision_heuristics.size();i++){
                                    Heuristic * hv = decision_heuristics[i];
                                    if(hv==h){
                                        break;
                                    }else{
                                        CRef ignore;
                                        Lit n2 = hv->decideTheory(ignore);
                                        assert(n2==lit_Undef);

                                    }
                                }*/
#endif
                            }
                        }
                        if(next == lit_Undef){
                            theory_order_heap.removeMin();
                            if(decisionLevel() >
                               0){//this isn't quite right... just because a heuristic cannot suggest anything at level 0, does not mean it cannot suggest something at a later level...
                                theory_decision_trail.push({h, decisionLevel()});
                                /*if(h->getHeuristicIndex()==7){
									printf("%d %d\n",h->getHeuristicIndex(),decisionLevel());
									CRef ignore;
									Lit n2 = h->decideTheory(ignore);
								}*/
                            }
                        }else{
                            assert(var(next) >= 0);
                            assert(var(next) < nVars());
                            assert(value(next) == l_Undef);
                            if(value(next) != l_Undef){
                                throw std::runtime_error("Bad theory decision (already assigned)");
                            }
                        }
                    }
                }else if(opt_randomize_theory_order_all || opt_theory_decision_round_robin){
                    if(opt_randomize_theory_order_all){
                        randomShuffle(random_seed, decision_heuristics);
                    }

                    for(int i = 0; i < decision_heuristics.size() && next == lit_Undef; i++){
                        //j only differs from i if round robin theory decisions are enabled
                        int j = (i + theory_decision_round_robin) % decision_heuristics.size();
                        assert(opt_theory_decision_round_robin || j == i);
                        assert(j >= 0);
                        assert(j <= decision_heuristics.size());
                        Heuristic* t = decision_heuristics[j];
                        if(!heuristicSatisfied(t) && t->getPriority() >= next_var_priority){
                            next = t->decideTheory(decision_reason);
                            next_decision_heuristic = t;
                        }
                    }
                    if(opt_theory_decision_round_robin){
                        theory_decision_round_robin++;
                        theory_decision_round_robin %= decision_heuristics.size();
                    }
                }/*else{


                    for (decision_heuristic_qhead < decision_heuristics.size() && next == lit_Undef; decision_heuristic_qhead++) {
                        //j only differs from i if round robin theory decisions are enabled

                        assert(decision_heuristic_qhead>=0);assert(decision_heuristic_qhead<=decision_heuristics.size());
                        Heuristic * t = decision_heuristics[decision_heuristic_qhead];
                        if (!heuristicSatisfied(t) &&  t->getPriority()>=next_var_priority) {
                            next = t->decideTheory(decision_reason);
                            next_decision_heuristic=t;
                        }
                    }

                }
*/
                if(next != lit_Undef){
                    stats_theory_decisions++;
                    if(theoryDecision != lit_Undef && var(next) == var(theoryDecision)){
                        assigns[var(theoryDecision)] = l_Undef;
                    }

                    if(decision_reason != CRef_Undef){
                        Heuristic* h = getHeuristic(decision_reason);
                        if(h && first_heuristic_decision_level[h->getHeuristicIndex()] < 0){
                            first_heuristic_decision_level[h->getHeuristicIndex()] = decisionLevel();
                        }
                    }
                }else{
                    next_decision_heuristic = nullptr;
                }
            }
            {

                bool has_last_heuristic = last_decision_heuristic;
                decision_heuristic_changed = has_last_heuristic && last_decision_heuristic != next_decision_heuristic;
                if(opt_decide_theories_only_prop_decision && decision_heuristic_changed && !propagate_theories &&
                   has_last_heuristic){
                    assert(trail_lim.size() > 0);
                    cancelUntil(decisionLevel() - 1);
                    goto propagate;
                }else{
                    //update the last decision heuristic only if the above condition was not triggered
                    last_decision_heuristic = next_decision_heuristic;
                }
            }
            if(next == lit_Undef){
                // New variable decision:

                next = pickBranchLit();
                // int p = priority[var(next)];
                if(opt_verb > 2){
                    printf("solver decision\n");
                }
                if(next == lit_Undef){

                    //solve theories if this solver is completely assigned
                    if(!disable_theories){
                        if(!propagate_theories){
                            //if for any reason theories were not propagated on the last run, then attempt to propagate them again
                            confl = propagate(true);
                            if(confl != CRef_Undef){
                                goto conflict;
                            }
                        }
                        for(int i = 0; i < theories.size(); i++){
                            if(opt_subsearch == 3 && track_min_level < initial_level)
                                continue; //Disable attempting to solve sub-solvers if we've backtracked past the super solver's decision level
                            if(!theories[i]->solveTheory(theory_conflict)){
                                if(!addConflictClause(theory_conflict, confl)){
                                    goto conflict;
                                }else{
                                    goto propagate;
                                }
                            }
                            //If propagating one of the sub theories caused this solver to backtrack, then go back to propagation
                            if(qhead < trail.size() || hasNextDecision())
                                goto propagate;
                        }
                    }
                    // Model found:
                    return l_True;
                }
            }

            //last_dec = var(next);
            // Increase decision level and enqueue 'next'
            assert(next != lit_Undef);
            if(next == lit_Undef){
                throw std::runtime_error("Bad decision");
            }
            if(var(next) >= nVars()){
                throw std::runtime_error("Bad decision (unknown variable)");
            }
            if(value(next) != l_Undef){
                throw std::runtime_error("Bad decision (already assigned)");
            }
            if(decisionLevel() > nVars() + assumptions.size() + 1){
                throw std::runtime_error("Bad decision level");
            }

            //if(next!=lit_Error)//lit_Error is used to signify a decision that has no literal in the SAT solver (some theories may support this)
            assert(value(next) == l_Undef);
            //Notice that the decision reason here may be CRef_Undef, and in any case is only used to identify which heuristic suggested this decision (it is _not_ used for conflict analysis)
            enqueue(next,
                    decision_reason);//not unchecked enqueue, because a theory solver _may_ have assigned this literal while making a decision
        }
    }
    //Unreachable
    return l_Undef;
}

double Solver::progressEstimate() const{
    double progress = 0;
    double F = 1.0 / nVars();

    for(int i = 0; i <= decisionLevel(); i++){
        int beg = i == 0 ? 0 : trail_lim[i - 1];
        int end = i == decisionLevel() ? trail.size() : trail_lim[i];
        progress += pow(F, i) * (end - beg);
    }

    return progress / nVars();
}


bool Solver::propagateAssignment(const vec<Lit>& assumps){
    only_propagate_assumptions = true;
    budgetOff();
    assumps.copyTo(assumptions);

    lbool val = solve_();
    only_propagate_assumptions = false;
    return val == l_True;
}

lbool Solver::solveUntilRestart(const vec<Lit>& assumps){
    quit_at_restart = true;
    //budgetOff();
    assumps.copyTo(assumptions);

    lbool val = solve_();
    quit_at_restart = false;
    return val;
}

// NOTE: assumptions passed in member-variable 'assumptions'.
lbool Solver::solve_(){
#ifndef NDEBUG
    if(!shown_warning){
        shown_warning = true;
        fprintf(stderr, "Warning: MonoSAT compiled without -DNDEBUG (will be very slow!).\n");
    }
#endif
#if  defined(DEBUG_MONOSAT) || defined(DEBUG_CORE) || defined(DEBUG_SOLVER) || defined(DEBUG_DGL) || defined(DEBUG_BV) || defined(DEBUG_GRAPH) || defined(DEBUG_THEORY) || defined(DEBUG_GEOMETRY) || defined(DEBUG_PB) || defined(DEBUG_FSM)
    if(!shown_warning){
        shown_warning = true;
        fprintf(stderr, "Warning: MonoSAT was compiled with DEBUG_* (will be very slow!).\n");
    }
#endif
    clearInterrupt();
    cancelUntil(0);
    model.clear();
    conflict.clear();
    if(!ok)
        return l_False;
    if(pbsolver){
        pbsolver->convert();
    }

    solves++;

    max_learnts = nClauses() * learntsize_factor;
    learntsize_adjust_confl = learntsize_adjust_start_confl;
    learntsize_adjust_cnt = (int) learntsize_adjust_confl;
    lbool status = l_Undef;

    if(verbosity >= 1 && !printed_header){
        //on repeated calls, don't print the header again
        printed_header = true;
        printf("============================[ Search Statistics ]==============================\n");
        printf("| Conflicts |          ORIGINAL         |          LEARNT          | Progress |\n");
        printf("|           |    Vars  Clauses Literals |    Limit  Clauses Lit/Cl |          |\n");
        printf("===============================================================================\n");
    }
    initial_level = 0;
    track_min_level = 0;
    //ensure that any theory atoms that were created _after_ the variable was assigned are enqueued in the theory
    //this can be improved
    for(int i = 0; i < qhead; i++){
        Lit p = trail[i];
        for(int n = 0; n < getNTheories(var(p)); n++){
            int theoryID = getTheoryID(p, n);
            Lit theoryLit = getTheoryLit(p, n);
            theories[theoryID]->enqueueTheory(theoryLit);
        }
    }
    // Search:
    int curr_restarts = 0;
    if(quit_at_restart && override_restart_count >= 0){
        curr_restarts = override_restart_count;
    }


    while(status == l_Undef){
        double rest_base = luby_restart ? luby(restart_inc, curr_restarts) : pow(restart_inc, curr_restarts);

        if(opt_rnd_phase){
            for(int i = 0; i < nVars(); i++)
                polarity[i] = irand(random_seed, 1);
        }
        if(opt_decide_theories){
            if(opt_theory_order_conflict_restart && opt_theory_order_conflict_clear_on_restart){
                theory_conflict_counters.clear();
                theory_conflict_counters.growTo(all_decision_heuristics.size());
            }
            if(opt_randomize_theory_order_restart_freq > 0 &&
               drand(random_seed) < opt_randomize_theory_order_restart_freq){
                randomShuffle(random_seed, decision_heuristics);
                for(int i = 0; i < decision_heuristics.size(); i++){
                    decision_heuristics[i]->setHeuristicOrder(i);
                }
                theory_order_heap.build(decision_heuristics);
            }else if(opt_theory_order_initial_sort){
                sort(decision_heuristics, HeuristicActivityOrderLt());
                for(int i = 0; i < decision_heuristics.size(); i++){
                    decision_heuristics[i]->setHeuristicOrder(i);
                }
                theory_order_heap.build(decision_heuristics);
            }
            if(opt_verb >= 3){
                printf("Initial theory order: ");
                for(Heuristic* h:decision_heuristics){
                    printf("%d, ", h->getHeuristicIndex());
                }
                printf("\n");
            }
        }


        status = search(rest_base * restart_first);
        if(verbosity >= 1){
            printf("|r%9d | %7d %8d %8d | %8d %8d %6.0f | %" PRId64 " removed |\n", (int) conflicts,
                   (int) dec_vars - (trail_lim.size() == 0 ? trail.size() : trail_lim[0]), nClauses(),
                   (int) clauses_literals, (int) max_learnts, nLearnts(),
                   (double) learnts_literals / nLearnts(), stats_removed_clauses);
        }

        if(!withinBudget()){

            //printf("Solver is giving up due to budget constraints: ");
            if(conflict_budget >= 0 && conflicts >= conflict_budget){
                //	printf("too many conflicts ");
            }
            if(propagation_budget >= 0 && propagations >= propagation_budget){
                //	printf("too many propagations ");
            }
            if(asynch_interrupt){
                //	printf("resource limit, memory limit, or user interupt triggered ");
            }
            //printf("\n");
            break;
        }
        curr_restarts++;
        if(opt_rnd_restart && status == l_Undef){

            for(int i = 0; i < nVars(); i++){
                activity[i] = drand(random_seed) * 0.00001;
            }
            rebuildOrderHeap();
        }


        if(opt_theory_order_vsids && opt_rnd_theory_restart && status == l_Undef){
            for(int i = 0; i < theories.size(); i++){
                theories[i]->setActivity(drand(random_seed) * 0.00001);
            }
            rebuildTheoryOrderHeap();
        }
        if(quit_at_restart && status == l_Undef){
            override_restart_count = curr_restarts;
            break;//give up
        }else{
            override_restart_count = -1;
        }
    }

    if(status == l_True){
        model.growTo(nVars());
        enqueueAnyUnqueued();//assign any remaining atoms to theories that were trivially satisfied
        // Extend & copy model:

        for(int i = 0; i < nVars(); i++)
            model[i] = value(i);

        for(Lit l:assumptions){
            if(value(l) != l_True){
                throw std::runtime_error("Model is inconsistent with assumptions!");
            }
        }
        if(opt_check_solution){
            if(opt_verb > 0){
                printf("Checking witness...\n");
            }
            ClauseIterator it = clausesBegin();
            bool any_undef_clauses = false;
            while(it != clausesEnd()){
                const Clause& c = *it;
                bool all_false = true;
                bool any_true = false;
                for(Lit l:c){
                    if(value(l) == l_True){
                        any_true = true;
                    }
                    if(value(l) != l_False){
                        all_false = false;
                    }
                }
                if(all_false){
                    throw std::runtime_error("Unsatisfied clause");
                }
                if(!all_false && !any_true){
                    any_undef_clauses = true;
                }
                ++it;
            }
            if(any_undef_clauses){
                printf("Warning: Some clauses are neither SAT nor UNSAT (typically, this is due to variables that have been marked as not decidable)\n");
            }
        }


        if(opt_check_solution && theories.size() && !disable_theories && !only_propagate_assumptions){
            if(opt_verb > 0){
                printf("Checking theory witnesses...\n");
            }
            double check_start = rtime(1);
            for(Theory* t:theories){

                if(!t->check_solved()){
                    throw std::runtime_error("Error! Solution doesn't satisfy theory properties!");
                }
            }
            stats_solution_checking_time += rtime(1) - check_start;
        }
        if(opt_write_learnt_clauses && opt_debug_model){
            //write the model out
            fprintf(opt_write_learnt_clauses, "s SATISFIABLE\n");
            fprintf(opt_write_learnt_clauses, "v ");
            for(int i = 0; i < nVars(); i++){
                Lit l = mkLit(i, false);
                if(value(i) == l_True){
                    fprintf(opt_write_learnt_clauses, " %d", dimacs(unmap(l)));
                }else if(value(i) == l_False){
                    fprintf(opt_write_learnt_clauses, " %d", dimacs(~unmap(l)));
                }else{
                    fprintf(opt_write_learnt_clauses, " %d", dimacs(unmap(l)));//optional
                }
            }
            fprintf(opt_write_learnt_clauses, "\n");
            if(!disable_theories){
                for(Theory* t:theories){
                    std::stringstream ss;
                    t->writeTheoryWitness(ss);
                    fprintf(opt_write_learnt_clauses, "%s", ss.str().c_str());
                }
            }

        }

    }else if(status == l_False && conflict.size() == 0){
        ok = false;
    }else if(status == l_False){
        assert(ok);
    }
    quit_at_restart = false;
    only_propagate_assumptions = false;
    assumptions.clear();
    clearSatisfied();
    return status;
}

bool Solver::solveTheory(vec<Lit>& conflict_out){
    initial_level = decisionLevel();
    track_min_level = initial_level;
    lbool status = l_Undef;
    // Search:
    int curr_restarts = 0;
    conflict.clear();
    conflict_out.clear();
    while(status == l_Undef){
        double rest_base = luby_restart ? luby(restart_inc, curr_restarts) : pow(restart_inc, curr_restarts);
        status = search(rest_base * restart_first);
        if(!withinBudget() || (opt_subsearch == 0 && track_min_level < initial_level))
            break;
        curr_restarts++;
    }
    cancelUntil(track_min_level);
    if(track_min_level < initial_level){
        S->cancelUntil(track_min_level);
    }
    if(!ok){
        return false;
    }
    if(conflict.size()){
        toSuper(conflict.toVec(), conflict_out);
        interpolant.push();
        conflict_out.copyTo(interpolant.last());    //Add this clause into the interpolant vector
        return false;
    }

    return propagateTheory(conflict_out);
}


Var Solver::unmap(Var v){
    if(v >= nVars()){
        throw std::runtime_error("Cannot unmap non-existant variable: " + std::to_string(v));
    }
    if(varRemap){
        return varRemap->unmap(v);
    }else{
        return v;
    }
}

Lit Solver::unmap(Lit l){
    if(var(l) >= nVars()){
        throw std::runtime_error("Cannot unmap non-existant variable: " + std::to_string(var(l)));
    }
    if(varRemap){
        return varRemap->unmap(l);
    }else{
        return l;
    }
}

Var Solver::mapVar(Var v){
    if(varRemap){
        return varRemap->mapVar(*this, v);
    }else{
        return v;
    }
}

int Solver::nMappedVars(){
    if(varRemap){
        return varRemap->nMappedVars();
    }else{
        return nVars();
    }
}

Lit Solver::mapLit(Lit l){
    if(varRemap){
        return varRemap->mapLit(*this, l);
    }else{
        return l;
    }
}


//=================================================================================================
// Garbage Collection methods:

void Solver::relocAll(ClauseAllocator& to){
    if(tmp_clause != CRef_Undef){
        ca[tmp_clause].mark(1);    //is this needed?
        ca.free(tmp_clause);
        tmp_clause = CRef_Undef;
    }
    //Re-allocate the 'theory markers'
    vec<int> marker_theory_tmp;
    // All watchers:
    //
    watches.cleanAll();
    for(int v = 0; v < nVars(); v++){

        for(int s = 0; s < 2; s++){
            Lit p = mkLit(v, s);
            vec<Watcher>& ws = watches[p];
            for(int j = 0; j < ws.size(); j++)
                ca.reloc(ws[j].cref, to);
        }
    }
    // All reasons:
    //
    for(int i = 0; i < trail.size(); i++){
        Var v = var(trail[i]);

        if(ca.isClause(reason(v)) && (ca[reason(v)].reloced() || locked(ca[reason(v)])))
            ca.reloc(vardata[v].reason, to);
    }

    // All learnt:
    //
    for(int i = 0; i < learnts.size(); i++)
        ca.reloc(learnts[i], to);

    // All original:
    //
    for(int i = 0; i < clauses.size(); i++)
        ca.reloc(clauses[i], to);
}

void Solver::garbageCollect(){
    // Initialize the next region to a size corresponding to the estimated utilization degree. This
    // is not precise but should avoid some unnecessary reallocations for the new region:
    ClauseAllocator to(ca.size() - ca.wasted());
    relocAll(to);
    if(verbosity >= 2)
        printf("|  Garbage collection:   %12d bytes => %12d bytes             |\n",
               ca.size() * ClauseAllocator::Unit_Size, to.size() * ClauseAllocator::Unit_Size);
    to.moveTo(ca);
}

void Solver::addLiteralName(Lit l, const std::string& name){
    if(var(l) < 0 || var(l) >= nVars()){
        throw std::invalid_argument("No such variable");
    }

    if(name.size() > 0){

        if(hasLiteral(name)){
            if(getLiteral(name) != l){
                throw std::invalid_argument("All literal names must be unique: " + name);
            }else{
                //this variable already has this name, do nothing
                return;
            }
        }else{
            //check if any chars of name are illegal
            for(char c:name){
                if(!isascii(c) || !isprint(c) || isspace(c) || (c == '~')){
                    throw std::invalid_argument(std::string(
                            "Variable names must consist only of printable, non-whitespace ASCII, and may not include '~'. Invalid character in variable name: ") +
                                                name);
                }
            }
        }

        if(literalNameCount(l) == 0){
            named_literals.push(l);
            litnames.insert({toInt(l), std::vector<std::string>()});
        }
        if(literalNameCount(~l) == 0){
            litnames.insert({toInt(~l), std::vector<std::string>()});
        }

        litnames[toInt(l)].push_back(name);
        assert(!namemap.count(name) || namemap[name] == lit_Undef);
        namemap[name] = l;

        std::string neg_name = "~" + name;
        litnames[toInt(~l)].push_back(neg_name);
        assert(!namemap.count(neg_name) || namemap[neg_name] == lit_Undef);
        namemap[neg_name] = ~l;
    }
}
