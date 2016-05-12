/***************************************************************************************[Solver.cc]
 The MIT License (MIT)

 Copyright (c) 2014, Sam Bayless
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

#include <math.h>

#include "mtl/Sort.h"
#include "core/Solver.h"
#include "core/Config.h"
#include "graph/GraphTheory.h"
#include <unistd.h>
#include "core/Remap.h"
using namespace Monosat;

//=================================================================================================
// Options:
// Collected in Config.h

//=================================================================================================
// Constructor/Destructor:

Solver::Solver() :
		
		// Parameters (user settable):
		//
		verbosity(opt_verb), var_decay(opt_var_decay), clause_decay(opt_clause_decay), theory_decay(opt_var_decay), random_var_freq(
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
				0), stats_pure_theory_lits(0), pure_literal_detections(0), stats_removed_clauses(0), dec_vars(0), clauses_literals(
				0), learnts_literals(0), max_literals(0), tot_literals(0), stats_pure_lit_time(0),  ok(
				true), cla_inc(1), var_inc(1), theory_inc(1), watches(WatcherDeleted(ca)), qhead(0), simpDB_assigns(-1), simpDB_props(
				0), order_heap(VarOrderLt(activity, priority)),theory_order_heap(TheoryOrderLt(theories)), progress_estimate(0), remove_satisfied(true) //lazy_heap( LazyLevelLt(this)),

		// Resource constraints:
		//
				, conflict_budget(-1), propagation_budget(-1) {
			if(opt_vsids_solver_as_theory){
				decisionTheory = new SolverDecisionTheory(*this);
				this->addTheory(decisionTheory);
			}
}

Solver::~Solver() {
	for (Theory * t : theories) {
		delete (t);
	}
}

//=================================================================================================
// Minor methods:

// Creates a new SAT variable in the solver. If 'decision' is cleared, variable will not be
// used as a decision variable (NOTE! This has effects on the meaning of a SATISFIABLE result).
//
Var Solver::newVar(bool sign, bool dvar) {

	int v;
    if (free_vars.size() > 0){
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
	assigns[v]=l_Undef;
	vardata[v]=mkVarData(CRef_Undef, 0);
	int p = 0;
	if (max_decision_var > 0 && v > max_decision_var)
		p = 1;
	if (v < min_decision_var)
		p = 1;
	priority[v]=p;
	theory_vars[v]=TheoryData();
	activity[v] = (rnd_init_act ? drand(random_seed) * 0.00001 : 0);
	seen[v] = 0;
	polarity[v]=opt_init_rnd_phase ? irand(random_seed, 1) : sign;
	//decision.push();//set below

	if (max_decision_var > 0 && v > max_decision_var)
		dvar = false;
	if (v < min_decision_var)
		dvar = false;
	setDecisionVar(v, dvar);

	return v;
}


// Note: at the moment, only unassigned variable will be released (this is to avoid duplicate
// releases of the same variable).
void Solver::releaseVar(Lit l)
{
    if (value(l) == l_Undef && !hasTheory(l)){//dont ever release theory atoms
        released_vars.push(var(l));
    }
    addClause(l);
}


bool Solver::addClause_(vec<Lit>& ps) {
	

	assert(decisionLevel() == 0);
	if (!ok)
		return false;
	resetInitialPropagation();    //Ensure that super solver call propagate on this solver at least once.
	// Check if clause is satisfied and remove false/duplicate literals:
	sort(ps);
	Lit p;
	int i, j;
	for (i = j = 0, p = lit_Undef; i < ps.size(); i++)
		if (value(ps[i]) == l_True || ps[i] == ~p)
			return true;
		else if (value(ps[i]) != l_False && ps[i] != p)
			ps[j++] = p = ps[i];
	ps.shrink(i - j);
	
	if (ps.size() == 0)
		return ok = false;
	else if (ps.size() == 1) {
		uncheckedEnqueue(ps[0]);
		ok = (propagate(false) == CRef_Undef); //do NOT propagate theory solvers here, or else adding unit clauses can become very expensive in some circumstances (such as when constructing the initial CNF for example)
				
		return ok;
	} else {
		CRef cr = ca.alloc(ps, false);
		clauses.push(cr);
		attachClause(cr);
	}
	
	return true;
}

CRef Solver::attachReasonClause(Lit r,vec<Lit> & ps) {
	assert(value(r)==l_True);

	if(opt_write_learnt_clauses){
		if((++opt_n_learnts)==47){
			int a=1;
		}

		fprintf(opt_write_learnt_clauses,"learnt ");
		for (Lit l:ps){
			fprintf(opt_write_learnt_clauses," %d", dimacs(unmap(l)));
		}
		fprintf(opt_write_learnt_clauses," 0\n");
		fflush(opt_write_learnt_clauses);
	}

	//sort(ps);
	Lit p;
	int i, j;
	for (i = j = 0, p = lit_Undef; i < ps.size(); i++)
		if (((value(ps[i]) == l_True && level(var(ps[i])) == 0)) || ps[i] == ~p)
			return CRef_Undef;
		else if ((value(ps[i]) != l_False || level(var(ps[i])) != 0) && ps[i] != p)
			ps[j++] = p = ps[i];
	ps.shrink(i - j);
	
	CRef confl_out = CRef_Undef;
	if (ps.size() == 0) {
		ok = false;
		//cancelUntil(0);
		return CRef_Undef;
	} else if (ps.size() == 1) {
		//cancelUntil(0);
		assert(var(ps[0]) < nVars());

		enqueueLazy(ps[0],0);

		return CRef_Undef;
	} else {

		int nfalse = 0;
		int max_lev = 0;
		//bool satisfied = false;
		int notFalsePos1 = -1;
		int notFalsePos2 = -1;
		for (int j = 0; j < ps.size(); j++) {
			assert(var(ps[j]) < nVars());
			//assert(value(ps[j])==l_False);
			if (value(ps[j]) == l_False) {
				nfalse++;
			} else {

				//if (value(ps[j]) == l_True)
				//	satisfied = true;
				if (notFalsePos1 < 0)
					notFalsePos1 = j;
				else if (notFalsePos2 < 0) {
					notFalsePos2 = j;
				}
			}
			if (value(ps[j]) != l_Undef) {
				int l = level(var(ps[j]));
				if (l > max_lev) {
					max_lev = l;
				}
			}
		}

		if(notFalsePos1<0){
			//all literals in this clause are false - this is the usual case for reason clauses.
			//so just treat the first literal as the non-false one
			notFalsePos1=0;
		}

		assert(notFalsePos1 >= 0);
		if (notFalsePos1 >= 0 && notFalsePos2 >= 0) {
			assert(notFalsePos1 != notFalsePos2);
			if (notFalsePos1 == 1) {
				std::swap(ps[0], ps[notFalsePos2]);
			} else {
				std::swap(ps[0], ps[notFalsePos1]);
				std::swap(ps[1], ps[notFalsePos2]);
			}
		} else {
			std::swap(ps[0], ps[notFalsePos1]);
		}
		assert(value(ps[0])!=l_False);
		if (notFalsePos2 >= 0) {
			assert(value(ps[1])!=l_False);
		}


		CRef cr = ca.alloc(ps);
		ca[cr].setFromTheory(true);
		clauses.push(cr);
		attachClause(cr);
		enqueueLazy(ps[0],max_lev,cr);

		//find the highest level in the conflict (should be the current decision level, but we won't require that)
	/*	if (ps.size() > opt_temporary_theory_reasons) {
			if (tmp_clause == CRef_Undef) {
				tmp_clause = ca.alloc(ps, false);
				tmp_clause_sz = ps.size();
			} else if (tmp_clause_sz < ps.size()) {
				ca[tmp_clause].mark(1); //is this needed?
				ca.free(tmp_clause);
				tmp_clause = ca.alloc(ps, false);
				tmp_clause_sz = ps.size();
			} else {
				Clause & c = ca[tmp_clause];
				assert(tmp_clause_sz >= ps.size());
				assert(tmp_clause_sz >= c.size());
				c.grow(tmp_clause_sz - c.size());
				for (int i = 0; i < ps.size(); i++) {
					c[i] = ps[i];
				}
				c.shrink(c.size() - ps.size());
			}
			
			confl_out = tmp_clause;
		} else {
			CRef cr = ca.alloc(ps, !opt_permanent_theory_conflicts);
			ca[cr].setFromTheory(true);
			if (opt_permanent_theory_conflicts)
				clauses.push(cr);
			else {
				learnts.push(cr);
				if (--learntsize_adjust_cnt <= 0) {
					learntsize_adjust_confl *= learntsize_adjust_inc;
					learntsize_adjust_cnt = (int) learntsize_adjust_confl;
					max_learnts *= learntsize_inc;
					
					if (verbosity >= 1)
						printf("| %9d | %7d %8d %8d | %8d %8d %6.0f | %ld removed |\n", (int) conflicts,
								(int) dec_vars - (trail_lim.size() == 0 ? trail.size() : trail_lim[0]), nClauses(),
								(int) clauses_literals, (int) max_learnts, nLearnts(),
								(double) learnts_literals / nLearnts(), stats_removed_clauses);
				}
				
			}
			attachClause(cr);
			
			confl_out = cr;
		}*/
		return cr;
	}
}

void Solver::attachClause(CRef cr) {
	const Clause& c = ca[cr];
	assert(c.size() > 1);
#ifndef NDEBUG
	for (int p = 0; p <= 1; p++)
		if (value(c[p]) == l_False && level(var(c[p])) == 0) {
			//then c[0] must not have been propagated yet - it must be after qhead. Otherwise this clause will never be enforced
			if (decisionLevel() > 0) {
				assert(false); //because we have already propagated all the level 0 elements
			} else {
				bool found = false;
				for (int i = qhead; i < trail.size(); i++) {
					if (var(trail[i]) == var(c[p])) {
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
	if (c.learnt())
		learnts_literals += c.size();
	else
		clauses_literals += c.size();
}

void Solver::detachClause(CRef cr, bool strict) {
	const Clause& c = ca[cr];
	assert(c.size() > 1);
	
	if (strict) {
		remove(watches[~c[0]], Watcher(cr, c[1]));
		remove(watches[~c[1]], Watcher(cr, c[0]));
	} else {
		// Lazy detaching: (NOTE! Must clean all watcher lists before garbage collecting this clause)
		watches.smudge(~c[0]);
		watches.smudge(~c[1]);
	}
	
	if (c.learnt())
		learnts_literals -= c.size();
	else
		clauses_literals -= c.size();
}

void Solver::removeClause(CRef cr) {
	Clause& c = ca[cr];
	detachClause(cr);
	// Don't leave pointers to free'd memory!
	if (locked(c))
		vardata[var(c[0])].reason = CRef_Undef;
	c.mark(1);
	ca.free(cr);
}

bool Solver::satisfied(const Clause& c) const {
	for (int i = 0; i < c.size(); i++)
		if (value(c[i]) == l_True)
			return true;
	return false;
}

// Revert to the state at given level (keeping all assignment at 'level' but not beyond).
//
void Solver::cancelUntil(int lev) {
	if (decisionLevel() > lev) {

		for (int i = 0; i < theories.size(); i++) {
			if(opt_lazy_backtrack  && theories[i]->supportsLazyBacktracking()){
				//if we _are_ backtracking lazily, then the assumption is that the theory solver will, after backtracking, mostly re-assign the same literals.
				//so instead, we will backtrack the theory lazily, in the future, if it encounters an apparent conflict (and this backtracking may alter or eliminate that conflict.)
			}else{
				theories[i]->backtrackUntil(lev);
			}
		}

		//printf("s: cancel %d\n",lev);
		for (int c = trail.size() - 1; c >= trail_lim[lev]; c--) {
			Var x = var(trail[c]);
			int xlev = level(x);

			if(xlev<=lev){
				to_reenqueue.push(trail[c]);
			}else{
				if(hasTheory(x)){
					theories[getTheoryID(x)]->undecideTheory(getTheoryLit(trail[c]));
				}
				assigns[x] = l_Undef;
				if (phase_saving > 1 || ((phase_saving == 1) && c > trail_lim.last()))
					polarity[x] = sign(trail[c]);
				insertVarOrder(x);
			}
		}
		qhead = trail_lim[lev];
		if (local_qhead > qhead) {
			local_qhead = qhead;
		}
		if (S && super_qhead > S->qhead) {
			super_qhead = S->qhead;
		}
		
		trail.shrink(trail.size() - trail_lim[lev]);
		trail_lim.shrink(trail_lim.size() - lev);

		//remove any lits from the lazy heap that are now unassigned.
/*		while(lazy_heap.size() && value(toLit( lazy_heap.peekMin())) == l_Undef ){
			lazy_heap.pop();
		}*/

		if (decisionLevel() < track_min_level) {
			track_min_level = decisionLevel();
		}
		for (int i : theory_queue) {
			in_theory_queue[i] = false;
		}
		theory_queue.clear();



		//re-enqueue any lazy lits as appropriate
		while(to_reenqueue.size()){
			Lit p = to_reenqueue.last();
			to_reenqueue.pop();
			assert(level(var(p)) <=lev);
			trail.push(p);
			//is this really needed?
			if (hasTheory(p)) {
				int theoryID = getTheoryID(p);
				needsPropagation(theoryID);
				theories[theoryID]->enqueueTheory(getTheoryLit(p));
			}
		}
		//should qhead be adjusted, here? Or do we want to repropagate these literals? (currently, re-propagating these literals).

		//If opt_theory_order_vsids is enabled, need to add back into the theory decision heap theories that may now be able to make decisions again.
		while(theory_decision_trail.size() && theory_decision_trail.last().second>=decisionLevel()){
			int theoryID = theory_decision_trail.last().first;
			assert(!theory_order_heap.inHeap(theoryID));
			theory_order_heap.insert(theoryID);
			theory_decision_trail.pop();
		}
	}
}

void Solver::backtrackUntil(int lev) {
	if (S->trail.size() < super_qhead)
		super_qhead = S->trail.size();
}

//=================================================================================================
// Major methods:

Lit Solver::pickBranchLit() {
	Var next = var_Undef;
	decisions++;
	// Random decision:
	if (drand(random_seed) < random_var_freq && !order_heap.empty()) {
		next = order_heap[irand(random_seed, order_heap.size())];
		if (value(next) == l_Undef && decision[next])
			rnd_decisions++;
	}
	
	// Activity based decision:
	while (next == var_Undef || value(next) != l_Undef || !decision[next])
		if (order_heap.empty()) {
			next = var_Undef;
			break;
		} else
			next = order_heap.removeMin();

	return next == var_Undef ? lit_Undef : mkLit(next, rnd_pol ? drand(random_seed) < 0.5 : polarity[next]);
}

void Solver::instantiateLazyDecision(Lit p,int atLevel, CRef reason){
	assert(value(p) == l_Undef);
	assigns[var(p)] = lbool(!sign(p));
	vardata[var(p)] = mkVarData(reason, atLevel);
	assert(atLevel<=decisionLevel());
	assert(atLevel>0);
	int trail_pos = trail_lim[atLevel-1];
	Lit curDec = trail[trail_pos];
	if(curDec==p)//this lit was already the decision at this level
		return;
	assert(curDec==theoryDecision);
	if(curDec!=theoryDecision){
		exit(6);
	}
	trail[trail_pos]=p;

	if (hasTheory(p)) {
		int theoryID = getTheoryID(p);
		needsPropagation(theoryID);
		theories[theoryID]->enqueueTheory(getTheoryLit(p));
	}
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
void Solver::analyze(CRef confl, vec<Lit>& out_learnt, int& out_btlevel) {
	int pathC = 0;
	CRef original_confl = confl;
	Lit p = lit_Undef;
	assert(confl != CRef_Undef);
	int maxlev = 0;
	{
		Clause & check = ca[confl];
		for(Lit p:check){
			int lev = level(var(p));
			if(lev>maxlev){
				maxlev=lev;
			}
		}
/*		int maxcount = 0;
		for(Lit l:check){
			if(level(var(l))>=maxlev){
				maxcount++;
			}
		}
		if(maxcount==1){
		//	exit(7);
		}*/
	}

	cancelUntil(maxlev);//use of lazily enqueued literals can trigger conflicts at earlier decision levels
	// Generate conflict clause:
	//
#ifndef NDEBUG
	assert(!seen.contains(1));
#endif
	bool possibly_missed_1uip=false;
	to_analyze.clear();
	out_learnt.push(lit_Undef);      // (leave room for the asserting literal)
	int index = trail.size() - 1;
	int stop = trail_lim.last();
	do {
		if (confl != CRef_Undef) {
			assert(!isTheoryCause(confl));
			Clause& c = ca[confl];
			
			if (c.learnt())
				claBumpActivity(c);
			
			for (int j = (p == lit_Undef) ? 0 : 1; j < c.size(); j++) {
				Lit q = c[j];
				if (!seen[var(q)] && level(var(q)) > 0) {
					assert(value(q)==l_False);assert(var(q)!=var(theoryDecision));
					varBumpActivity(var(q));
					seen[var(q)] = 1;
					if(var(q)==10){
						int a=1;
					}
					if (level(var(q)) >= decisionLevel())
						pathC++;
					else
						out_learnt.push(q);
				}
			}
		} else {
			out_learnt.push(~p);
		}
		bool searching=true;
		while(searching){
			searching=false;

			// Select next clause to look at:
			while (index>= stop && (!seen[var(trail[index--])] || (level(var(trail[index+1]))<decisionLevel())));
			assert(index >= -1);
			p = trail[index + 1];
			if(var(p)==10){
				int a=1;
			}
			assert(var(p)!=var(theoryDecision));
			confl = reason(var(p));
			int was_at_level = level(var(p));
			if (isTheoryCause(confl)) {
				//lazily construct the reason for this theory propagation now that we need it
				confl = constructReason(p);
				//for some theories, we may discover while constructing the cause that p is at a lower level than we thought.
			}
			if(level(var(p))<decisionLevel()){
				assert(value(p)==l_True);
				if(was_at_level==decisionLevel()){
					//the level of this variable changed while deriving a reason for it.
					if(level(var(p))>0){
						out_learnt.push(~p);
					}
					pathC--;
					possibly_missed_1uip=true;
				}
				seen[var(p)] = 0;
				searching=pathC>0;
			}
		}

		seen[var(p)] = 0;
		pathC--;
		
	} while (pathC > 0);
	out_learnt[0] = ~p;
	assert(var(p)!=var(theoryDecision));

	if(possibly_missed_1uip){
		//because of literals that were enqueued lazily, at the wrong level, and only discovered during lazy reason construction,
		//in rare circumstances the 1uip may be missed. If that happens, re-start clause learning now that all the relevant reasons have
		//been constructed and all relevant levels are corrected.
		for(Lit l:out_learnt){
			seen[var(l)]=0;
		}
		assert(!seen.contains(1));
		out_learnt.clear();
		analyze(original_confl,out_learnt,out_btlevel);
		return;
	}

#ifndef NDEBUG
	for (Lit p : out_learnt)
		assert(value(p)==l_False);

#endif
	// Simplify conflict clause:
	//
	int i, j;
	out_learnt.copyTo(analyze_toclear);
	if (ccmin_mode == 2) {
		uint32_t abstract_level = 0;
		for (i = 1; i < out_learnt.size(); i++)
			abstract_level |= abstractLevel(var(out_learnt[i])); // (maintain an abstraction of levels involved in conflict)
					
		for (i = j = 1; i < out_learnt.size(); i++)
			if (!ca.isClause(reason(var(out_learnt[i]))) || !litRedundant(out_learnt[i], abstract_level))
				out_learnt[j++] = out_learnt[i];
		
	} else if (ccmin_mode == 1) {
		for (i = j = 1; i < out_learnt.size(); i++) {
			Var x = var(out_learnt[i]);
			
			if (reason(x) == CRef_Undef)
				out_learnt[j++] = out_learnt[i];
			else {
				Clause& c = ca[reason(var(out_learnt[i]))];
				for (int k = 1; k < c.size(); k++)
					if (!seen[var(c[k])] && level(var(c[k])) > 0) {
						out_learnt[j++] = out_learnt[i];
						break;
					}
			}
		}
	} else
		i = j = out_learnt.size();
	
	max_literals += out_learnt.size();
	out_learnt.shrink(i - j);
	tot_literals += out_learnt.size();

	// Find correct backtrack level:
	//
	if (out_learnt.size() == 1)
		out_btlevel = 0;
	else {
		int max_i = 1;
		// Find the first literal assigned at the next-highest level:
		for (int i = 2; i < out_learnt.size(); i++)
			if (level(var(out_learnt[i])) > level(var(out_learnt[max_i])))
				max_i = i;
		// Swap-in this literal at index 1:
		Lit p = out_learnt[max_i];
		out_learnt[max_i] = out_learnt[1];
		out_learnt[1] = p;
		out_btlevel = level(var(p));
	}
#ifndef NDEBUG
	for (Lit p : out_learnt){
		assert(var(p)!=var(theoryDecision));
		assert(value(p)==l_False);
	}
#endif

	for (int j = 0; j < analyze_toclear.size(); j++)
		seen[var(analyze_toclear[j])] = 0;    // ('seen[]' is now cleared)
#ifndef NDEBUG
	assert(!seen.contains(1));
#endif
}

// Check if 'p' can be removed. 'abstract_levels' is used to abort early if the algorithm is
// visiting literals at levels that cannot be removed later.
bool Solver::litRedundant(Lit p, uint32_t abstract_levels) {
	analyze_stack.clear();
	analyze_stack.push(p);
	int top = analyze_toclear.size();
	while (analyze_stack.size() > 0) {
		assert(reason(var(analyze_stack.last())) != CRef_Undef);
		if (!ca.isClause(reason(var(analyze_stack.last())))) {
			//Don't pursue this if it is a subsolver reason
			for (int j = top; j < analyze_toclear.size(); j++)
				seen[var(analyze_toclear[j])] = 0;
			analyze_toclear.shrink(analyze_toclear.size() - top);
			analyze_stack.clear();
			return false;
		}
		
		Clause& c = ca[reason(var(analyze_stack.last()))];
		analyze_stack.pop();
		
		for (int i = 1; i < c.size(); i++) {
			Lit p = c[i];
			if (!seen[var(p)] && level(var(p)) > 0) {
				if (reason(var(p)) != CRef_Undef && (abstractLevel(var(p)) & abstract_levels) != 0) {
					seen[var(p)] = 1;
					analyze_stack.push(p);
					analyze_toclear.push(p);
				} else {
					for (int j = top; j < analyze_toclear.size(); j++)
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
void Solver::analyzeFinal(Lit p, vec<Lit>& out_conflict) {
	out_conflict.clear();
	out_conflict.push(p);
	
	if (decisionLevel() == 0)
		return;
	
	seen[var(p)] = 1;
	
	for (int i = trail.size() - 1; i >= trail_lim[0]; i--) {
		Var x = var(trail[i]);
		if (seen[x]) {
			if (reason(x) == CRef_Undef) {
				assert(level(x) > 0);
				out_conflict.push(~trail[i]);
			} else {
				if (isTheoryCause(reason(x))) {
					constructReason(trail[i]);
				}
				Clause& c = ca[reason(x)];
				assert(var(c[0]) == x);
				for (int j = 1; j < c.size(); j++)
					if (level(var(c[j])) > 0)
						seen[var(c[j])] = 1;
			}
			seen[x] = 0;
		}
	}
	
	seen[var(p)] = 0;
}

void Solver::enqueueLazy(Lit p, int lev, CRef from){
	assert(value(p)!=l_False);
	if(value(p)==l_True && lev < level(var(p))){
		//then the lit was already implied, but needs to be (lazily) moved to an earlier level.
		vardata[var(p)] = mkVarData(from, lev);
	}else if(value(p)==l_Undef && lev<decisionLevel()){
		assert(value(p) == l_Undef);
		assigns[var(p)] = lbool(!sign(p));
		vardata[var(p)] = mkVarData(from, lev);
		trail.push_(p);
		//lazy_heap.insert(toInt(p));
		if (hasTheory(p)) {
			int theoryID = getTheoryID(p);
			needsPropagation(theoryID);
			theories[theoryID]->enqueueTheory(getTheoryLit(p));
		}
	}else if(value(p)==l_Undef){
			uncheckedEnqueue(p, from);
	}else{
		//do nothing
	}
}

void Solver::uncheckedEnqueue(Lit p, CRef from) {
	assert(value(p) == l_Undef);
	assigns[var(p)] = lbool(!sign(p));
	vardata[var(p)] = mkVarData(from, decisionLevel());
	trail.push_(p);
	if (hasTheory(p)) {
		int theoryID = getTheoryID(p);
		needsPropagation(theoryID);
		theories[theoryID]->enqueueTheory(getTheoryLit(p));
	}
}

void Solver::analyzeFinal(CRef confl, Lit skip_lit, vec<Lit>& out_conflict) {
	out_conflict.clear();
	if (decisionLevel() == 0)
		return;
	
	Clause & c = ca[confl];
	
	for (int i = 0; i < c.size(); i++) {
		Var x = var(c[i]);
		if (x == var(skip_lit))
			continue;
		
		assert(x >= 0);
		assert(x < nVars());
		if (level(x) > 0)
			seen[x] = 1;
		assert(value(x)!=l_Undef);
	}
	
	int start = trail.size() - 1;
	int i;
	for (i = start; i >= trail_lim[0]; i--) {
		Var x = var(trail[i]);
		assert(value(x)!=l_Undef);
		assert(x >= 0);
		assert(x < nVars());
		
		if (seen[x]) {
			assert(x != var(skip_lit));
			CRef r = reason(x);
			if (r == CRef_Undef) {
				//Var v = var(trail[i]);
				//int lev = level(v);
				out_conflict.push(~trail[i]);
			} else {
				if (isTheoryCause(r)) {
					r = constructReason(trail[i]);
				}
				Clause& c = ca[r];
				for (int j = 0; j < c.size(); j++) {
					
					if (level(var(c[j])) > 0) {
						seen[var(c[j])] = 1;
					}
				}
			}
			seen[x] = 0;
		}
	}
}

void Solver::buildReason(Lit p, vec<Lit> & reason) {
	Lit local_l = fromSuper(p);    // mkLit(var(p)-super_offset, sign(p));
	analyzeFinal(local_l, reason);
	toSuper(reason, reason);
	interpolant.push();
	reason.copyTo(interpolant.last());    //Add this clause into the interpolant vector
}

void Solver::enqueueTheory(Lit l) {
	
}

//Propagate assignments from the super solver's interface variables to this solver (and, if this solver makes further assignments to the interface, pass those back to the super solver)
bool Solver::propagateTheory(vec<Lit> & conflict_out) {
	//backtrack as needed
	if (S->trail.size() && super_qhead > 0) {
		Lit last_super = S->trail[super_qhead - 1];
		cancelUntil(S->level(var(last_super)));
	} else {
		cancelUntil(0);
	}
	assert(decisionLevel() <= S->decisionLevel());
	
	CRef confl = CRef_Undef;

	while (confl == CRef_Undef && super_qhead < S->qhead) {
		Lit out_l = S->trail[super_qhead++];
		
		if (var(out_l) < min_super || var(out_l) > max_super)
			continue;    //this lit is not on the interface
			
		int lev = S->level(var(out_l));
		assert(decisionLevel() <= lev);
		while (decisionLevel() < lev) {
			
			//We are going to start a new decision level; make sure the current one is fully propagated first
			initial_level = decisionLevel();
			track_min_level = initial_level;
			
			confl = propagate();
			if (confl != CRef_Undef || track_min_level < initial_level) {
				cancelUntil(track_min_level);
				S->cancelUntil(track_min_level);
				goto conflict;
			}
			
			newDecisionLevel();
		}
		
		Lit local_l = fromSuper(out_l);
		
		if (!enqueue(local_l, CRef_Undef)) {
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
	if (confl != CRef_Undef || track_min_level < initial_level) {
		cancelUntil(track_min_level);
		S->cancelUntil(track_min_level);
	}
	conflict:

	if (confl != CRef_Undef) {
		//then we have a conflict which we need to instantiate in S
		conflict_out.clear();
		analyzeFinal(confl, lit_Undef, conflict_out);
		toSuper(conflict_out, conflict_out);
		interpolant.push();
		conflict_out.copyTo(interpolant.last());    //Add this clause into the interpolant vector
		return false;
	} else {
		//find lits to prop
		while (local_qhead < qhead) {
			Lit local_l = trail[local_qhead++];
			if (var(local_l) < min_local || var(local_l) > max_local)
				continue;    //this lit is not on the interface
			Lit out_l = toSuper(local_l);
			lbool v = S->value(out_l);
			
			if (S->value(out_l) != l_True || S->level(var(out_l)) > level(var(local_l))) {
				
				if (S->decisionLevel() > level(var(local_l))) {
					cancelUntil(level(var(local_l)));
					S->cancelUntil(level(var(local_l)));
					
				}
				if (!S->enqueue(out_l, cause_marker)) {    //can probably increment super_qhead here...
					//this is a conflict
					conflict_out.clear();
					analyzeFinal(local_l, conflict_out);
					toSuper(conflict_out, conflict_out);
					interpolant.push();
					conflict_out.copyTo(interpolant.last());			//Add this clause into the interpolant vector
					return false;
				}
			}
		}
		
		while (decisionLevel() < S->decisionLevel())
			newDecisionLevel();
		
		assert(decisionLevel() == S->decisionLevel());
		return true;
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
CRef Solver::propagate(bool propagate_theories) {
	if (qhead == trail.size() && (!initialPropagate || decisionLevel() > 0) && (!propagate_theories || !theory_queue.size())) {//it is possible that the theory solvers need propagation, even if the sat solver has an empty queue.
		return CRef_Undef;
	}
	theoryConflict=-1;
	CRef confl = CRef_Undef;
	int num_props = 0;
	watches.cleanAll();
	if (decisionLevel() == 0 && !propagate_theories) {
		initialPropagate = true;//we will need to propagate this assignment to the theories at some point in the future.
	}
	do {
		
		while (qhead < trail.size()) {
			if (opt_early_theory_prop) {
				//propagate theories;
				while (propagate_theories && theory_queue.size() && confl == CRef_Undef) {
					theory_conflict.clear();
					int theoryID = theory_queue.last();
					theory_queue.pop();
					in_theory_queue[theoryID] = false;
					if (!theories[theoryID]->propagateTheory(theory_conflict)) {
						if (!addConflictClause(theory_conflict, confl)) {
							qhead = trail.size();
							return confl;
						}
					}
				}
			}
			
			Lit p = trail[qhead++];     // 'p' is enqueued fact to propagate.

			vec<Watcher>& ws = watches[p];
			Watcher *i, *j, *end;
			num_props++;
			for (i = j = (Watcher*) ws, end = i + ws.size(); i != end;) {
				// Try to avoid inspecting the clause:
				Lit blocker = i->blocker;
				if (value(blocker) == l_True) {
					*j++ = *i++;
					continue;
				}
				
				// Make sure the false literal is data[1]:
				CRef cr = i->cref;
				Clause& c = ca[cr];
				Lit false_lit = ~p;
				if (c[0] == false_lit)
					c[0] = c[1], c[1] = false_lit;
				assert(c[1] == false_lit);
				i++;
				
				// If 0th watch is true, then clause is already satisfied.
				Lit first = c[0];
				Watcher w = Watcher(cr, first);
				if (first != blocker && value(first) == l_True) {
					*j++ = w;
					continue;
				}
				
				// Look for new watch:
				for (int k = 2; k < c.size(); k++)
					if (value(c[k]) != l_False) {
						c[1] = c[k];
						c[k] = false_lit;
						watches[~c[1]].push(w);
						goto NextClause;
					}
				
				// Did not find watch -- clause is unit under assignment:
				*j++ = w;
				if (value(first) == l_False) {
					confl = cr;
					qhead = trail.size();
					// Copy the remaining watches:
					while (i < end)
						*j++ = *i++;
				} else
					uncheckedEnqueue(first, cr);
				
				NextClause: ;
			}
			ws.shrink(i - j);
		}
		
		if (initialPropagate && decisionLevel() == 0 && propagate_theories) {
			assert(decisionLevel() == 0);
			//propagate any as yet unpropagated literals to each theory
			for (int i = 0; i < qhead; i++) {
				Lit p = trail[i];
				if (hasTheory(p)) {
					int theoryID = getTheoryID(p);
					theories[theoryID]->enqueueTheory(getTheoryLit(p));
					
				}
			}
			//Have to check _all_ the theories, even if we haven't eneueued of their literals, in case they have constants to propagate up.
			for (int theoryID = 0; theoryID < theories.size(); theoryID++) {
				needsPropagation(theoryID);
			}
			initialPropagate = false;
		}
		static int iter = 0;
		if (++iter == 36) {

			int a = 1;
		}
		//printf("iter %d\n",iter);
		//printf("iter %d\n",iter);
		//propagate theories;
		while (propagate_theories && theory_queue.size() && (opt_early_theory_prop || qhead == trail.size())
				&& confl == CRef_Undef) {
			theory_conflict.clear();
			//todo: ensure that the bv theory comes first, as otherwise dependent theories may have to be propagated twice...
			int theoryID = theory_queue.last();
			theory_queue.pop();
			in_theory_queue[theoryID] = false;
			if (!theories[theoryID]->propagateTheory(theory_conflict)) {
				bool has_conflict=true;
#ifndef NDEBUG
				for(Lit l:theory_conflict)
					assert(value(l)!=l_Undef);
#endif
				if (has_conflict && !addConflictClause(theory_conflict, confl)) {
					theoryConflict=theoryID;

					qhead = trail.size();
					stats_theory_conflicts++;

					return confl;
				}
			}
		}
		
		//solve theories if this solver is completely assigned
		/*	for(int i = 0;i<theories.size() && qhead == trail.size() && confl==CRef_Undef && nAssigns()==nVars();i++){
		 if(opt_subsearch==3 && track_min_level<initial_level)
		 continue;//Disable attempting to solve sub-solvers if we've backtracked past the super solver

		 if(!theories[i]->solve(theory_conflict)){
		 if(!addConflictClause(theory_conflict,confl))
		 return confl;

		 }
		 }*/

	} while (qhead < trail.size());
	propagations += num_props;
	simpDB_props -= num_props;
	
	return confl;
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
			ca(ca_) {
	}
	bool operator ()(CRef x, CRef y) {
		return ca[x].size() > 2 && (ca[y].size() == 2 || ca[x].activity() < ca[y].activity());
		/*if(ca[x].size()==2)
		 return false;
		 if(ca[y].size()==2)
		 return true;
		 if(ca[x].fromTheory() && !ca[y].fromTheory())
		 return true;
		 if(ca[y].fromTheory() && ! ca[x].fromTheory())
		 return false;
		 return ca[x].activity() < ca[y].activity();*/
	}
};
void Solver::reduceDB() {
	int i, j;
	double extra_lim = cla_inc / learnts.size();    // Remove any clause below this activity
			
	sort(learnts, reduceDB_lt(ca));
	// Don't delete binary or locked clauses. From the rest, delete clauses from the first half
	// and clauses with activity smaller than 'extra_lim':
	for (i = j = 0; i < learnts.size(); i++) {
		Clause& c = ca[learnts[i]];
		if (c.size() > 2 && !locked(c) && (i < learnts.size() / 2 || c.activity() < extra_lim)) {
			stats_removed_clauses++;
			removeClause(learnts[i]);
		} else
			learnts[j++] = learnts[i];
	}
	learnts.shrink(i - j);
	checkGarbage();
}

void Solver::removeSatisfied(vec<CRef>& cs) {
	int i, j;
	for (i = j = 0; i < cs.size(); i++) {
		Clause& c = ca[cs[i]];
		if (satisfied(c))
			removeClause(cs[i]);
		else
			cs[j++] = cs[i];
	}
	cs.shrink(i - j);
}

void Solver::rebuildOrderHeap() {
	vec<Var> vs;
	for (Var v = 0; v < nVars(); v++)
		if (decision[v] && value(v) == l_Undef)
			vs.push(v);
	order_heap.build(vs);
}
void Solver::rebuildTheoryOrderHeap() {

	theory_order_heap.clear();
	for(int i =0;i<theories.size();i++){
		if(theories[i]->supportsDecisions()){
			theory_order_heap.insert(i);
		}
	}
}

void Solver::detectPureTheoryLiterals(){
	//Detect any theory literals that are occur only in positive or negative (or neither) polarity
	if (opt_detect_pure_lits) {
		//if(pure_literal_detections==0){

		pure_literal_detections++;
		double startTime = rtime(1);
		lit_counts.growTo(nVars() * 2);
		for (int i = 0; i < lit_counts.size(); i++) {
			lit_counts[i].seen = false;
		}
		assert(decisionLevel() == 0);

		//instead of counting the number of occurrences, just check if there are any occurences
		for (Lit l : trail) {
			lit_counts[toInt(l)].seen = true;
		}

		for (CRef cr : clauses) {
			Clause & c = ca[cr];
			for (Lit l : c) {
				lit_counts[toInt(l)].seen = true;
			}
		}
		//the learnt clauses likely don't need to be counted here...
		/*for (CRef cr : learnts) {
			Clause & c = ca[cr];
			for (Lit l : c) {
				lit_counts[toInt(l)].seen = true;
			}
		}*/
		for (Var v = 0; v < nVars(); v++) {

			Lit l = mkLit(v, false);
			if (lit_counts[toInt(l)].seen && !lit_counts[toInt(l)].occurs) {
				stats_pure_lits--;
				lit_counts[toInt(l)].occurs = true;
				if (hasTheory(v)) {
					stats_pure_theory_lits--;
					if(value(l)==l_Undef) {
						theories[getTheoryID(v)]->setLiteralOccurs(~getTheoryLit(l), true);
					}
				}
			}
			if (lit_counts[toInt(~l)].seen && !lit_counts[toInt(~l)].occurs) {
				stats_pure_lits--;
				lit_counts[toInt(~l)].occurs = true;
				if (hasTheory(v)) {
					stats_pure_theory_lits--;
					if(value(l)==l_Undef) {
						theories[getTheoryID(v)]->setLiteralOccurs(~getTheoryLit(~l), true);
					}
				}
			}

			if (lit_counts[toInt(l)].occurs && !lit_counts[toInt(l)].seen) {
				lit_counts[toInt(l)].occurs = false;
				stats_pure_lits++;
				if (hasTheory(v)) {
					stats_pure_theory_lits++;
					if(value(l)==l_Undef) {
						theories[getTheoryID(v)]->setLiteralOccurs(~getTheoryLit(l), false);
					}
					//setPolarity(v,false);
				} else {
					//we can safely assign this now...
					if (!lit_counts[toInt(~l)].seen) {
						/*	if(decision[v])
						 setDecisionVar(v,false);*/
					} else {
						if (value(l) == l_Undef) {
							//uncheckedEnqueue(~l);
						}
					}
				}
			}
			if (lit_counts[toInt(~l)].occurs && !lit_counts[toInt(~l)].seen) {
				lit_counts[toInt(~l)].occurs = false;
				stats_pure_lits++;
				if (hasTheory(v)) {
					stats_pure_theory_lits++;
					if(value(l)==l_Undef) {
						//level 0 variables are already handled through a separate mechanism in the theory solvers...
						theories[getTheoryID(v)]->setLiteralOccurs(~getTheoryLit(~l), false);
					}
					//setPolarity(v,true);
					if (!lit_counts[toInt(l)].occurs) {
						//setDecisionVar(v,false);//If v is pure in both polarities, and is a theory var, then just don't assign it at all - it is unconstrained.
						//This _should_ be safe to do, if the theory semantics make sense... although there are some clause learning schemes that introduce previosuly unused literals that might break with this...
					}
				} else {
					//we can safely assign this now...
					if (!lit_counts[toInt(l)].seen) {
						/*		if(decision[v])
						 setDecisionVar(v,false);*/
					} else if (value(l) == l_Undef) {
						//uncheckedEnqueue(l);
					}
				}
			}

		}
		assert(stats_pure_lits <= nVars() * 2);
		stats_pure_lit_time += rtime(1) - startTime;
		//}
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
bool Solver::simplify() {
	assert(decisionLevel() == 0);
	//it is important (for efficiency reasons, not correctness)
	//to run this pure literal detection before propagate, if new theory atoms have been added since the last call to propagate.
	if(ok && initialPropagate && opt_detect_pure_lits){
		detectPureTheoryLiterals();
	}

	if (!ok || propagate(!disable_theories) != CRef_Undef || !ok) //yes, the second ok check is now required, because propagation of a theory can make the solver unsat at this point...
		return ok = false;
	
	if (nAssigns() == simpDB_assigns || (simpDB_props > 0)) //if pure literal detection is enabled, then simplify MUST be called at least once after new theory atoms are added.
		return true;
	
	// Remove satisfied clauses:
	removeSatisfied(learnts);
	if (remove_satisfied){        // Can be turned off.
		removeSatisfied(clauses);

        // Remove all released variables from the trail:
        for (int i = 0; i < released_vars.size(); i++){
            assert(seen[released_vars[i]] == 0);
            seen[released_vars[i]] = 1;
        }

        int i, j;
        for (i = j = 0; i < trail.size(); i++)
            if (seen[var(trail[i])] == 0)
                trail[j++] = trail[i];
        trail.shrink(i - j);
        //printf("trail.size()= %d, qhead = %d\n", trail.size(), qhead);
        qhead = trail.size();

        for (int i = 0; i < released_vars.size(); i++)
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
void Solver::addClauseSafely(vec<Lit> & ps) {
	if(opt_write_learnt_clauses){
		if(opt_n_learnts++==47){
			int a=1;
		}
		if(ps.size()==2 &&  dimacs(unmap(ps[0])) ==-195621 && dimacs(unmap(ps[1]))== 35361){
			int a=1;
		}
		fprintf(opt_write_learnt_clauses,"learnt fact ");
		for (Lit l:ps){
			fprintf(opt_write_learnt_clauses," %d", dimacs(unmap(l)));
		}
		fprintf(opt_write_learnt_clauses," 0\n");
		fflush(opt_write_learnt_clauses);
	}

	if(decisionLevel()==0){
		addClause(ps);
	}else{
		//clauses_to_add.push();
		//ps.copyTo(clauses_to_add.last());
		//
		//check if this clause should imply a literal right now.
		//and IF that literal is already implied, then we may still need to adjust the level the literal should be implied at.
		int undef_count=0;
		int false_count=0;
		sort(ps);
		Lit p;
		int i, j;
		for (i = j = 0, p = lit_Undef; i < ps.size(); i++){
			assert(var(ps[i])!=var(theoryDecision));
			if (((value(ps[i]) == l_True && level(var(ps[i])) == 0)) || ps[i] == ~p)
				return;
			else if ((value(ps[i]) != l_False || level(var(ps[i])) != 0) && ps[i] != p){
				ps[j++] = p = ps[i];
				if(value(ps[i])==l_Undef)
					undef_count++;
				else if (value(ps[i])==l_False){
					false_count++;
				}
			}
		}
		ps.shrink(i - j);
		if(false_count==ps.size()-1){
			//this clause is unit under the current assignment.
			//although we _could_ wait until a restart to add this clause, in many cases this will lead to very poor solver behaviour.
			//so we are going to propagate this clause now.

			if(ps.size()==0){
				ok=false;
				return;
			}else if (ps.size()==1){
				enqueueLazy(ps[0],0,CRef_Undef);
			}else{

				bool conflicting = true;
				int nfalse = 0;
				int max_lev = 0;
				//bool satisfied = false;
				int notFalsePos1 = -1;
				int notFalsePos2 = -1;
				for (int j = 0; j < ps.size(); j++) {
					assert(var(ps[j]) < nVars());
					//assert(value(ps[j])==l_False);
					if (value(ps[j]) == l_False) {
						nfalse++;
					} else {
						conflicting = false;
						//if (value(ps[j]) == l_True)
						//	satisfied = true;
						if (notFalsePos1 < 0)
							notFalsePos1 = j;
						else if (notFalsePos2 < 0) {
							notFalsePos2 = j;
						}
					}
					if (value(ps[j]) != l_Undef) {
						int l = level(var(ps[j]));
						if (l > max_lev) {
							max_lev = l;
						}
					}
				}
				assert(!conflicting);
				if (!conflicting) {
					assert(notFalsePos1 >= 0);
					if (notFalsePos1 >= 0 && notFalsePos2 >= 0) {
						assert(notFalsePos1 != notFalsePos2);
						if (notFalsePos1 == 1) {
							std::swap(ps[0], ps[notFalsePos2]);
						} else {
							std::swap(ps[0], ps[notFalsePos1]);
							std::swap(ps[1], ps[notFalsePos2]);
						}
					} else {
						std::swap(ps[0], ps[notFalsePos1]);
					}
					assert(value(ps[0])!=l_False);
					if (notFalsePos2 >= 0) {
						assert(value(ps[1])!=l_False);
					}
				}

				CRef cr = ca.alloc(ps);
				ca[cr].setFromTheory(true);
				clauses.push(cr);
				attachClause(cr);
				enqueueLazy(ps[0],max_lev,cr);
			}
		}else{
			//add this clause at the next solver restart
			clauses_to_add.push();
			ps.copyTo(clauses_to_add.last());
		}




	}
}
bool Solver::addConflictClause(vec<Lit> & ps, CRef & confl_out, bool permanent) {
	static int nlearnt=0;
	if(++nlearnt== 111){
		int a=1;
	}
	if(opt_write_learnt_clauses){
		if(++opt_n_learnts==47){
			int a=1;
		}

		fprintf(opt_write_learnt_clauses,"learnt ");
		for (Lit l:ps){
			fprintf(opt_write_learnt_clauses," %d", dimacs(unmap(l)));
		}
		fprintf(opt_write_learnt_clauses," 0\n");
		fflush(opt_write_learnt_clauses);
	}
	//bool any_undef=false;
	sort(ps);
	Lit p;
	int i, j;
	for (i = j = 0, p = lit_Undef; i < ps.size(); i++){
		assert(var(ps[i])!=var(theoryDecision));
		//if(decisionLevel()>0 && value(ps[i])==l_Undef)
		//	any_undef=true;
		if (((value(ps[i]) == l_True && level(var(ps[i])) == 0)) || ps[i] == ~p)
			return true;
		else if ((value(ps[i]) != l_False || level(var(ps[i])) != 0) && ps[i] != p)
			ps[j++] = p = ps[i];

	}
	ps.shrink(i - j);
/*	if(any_undef){
		cancelUntil(0);//this is _not_ a conflict clause.
	}*/
	confl_out = CRef_Undef;
	if (ps.size() == 0) {
		ok = false;
		cancelUntil(0);
		return false;
	} else if (ps.size() == 1) {
		cancelUntil(0);

		assert(var(ps[0]) < nVars());

		if (!enqueue(ps[0])) {
			ok = false;
			return false;
		}

	}else {
		//find the highest level in the conflict (should be the current decision level, but we won't require that)
		bool conflicting = true;
		int nfalse = 0;
		int max_lev = 0;
		bool satisfied = false;
		int notFalsePos1 = -1;
		int notFalsePos2 = -1;
		for (int j = 0; j < ps.size(); j++) {
			assert(var(ps[j]) < nVars());
			//assert(value(ps[j])==l_False);
			if (value(ps[j]) == l_False) {
				nfalse++;
			} else {
				conflicting = false;
				if (value(ps[j]) == l_True)
					satisfied = true;
				if (notFalsePos1 < 0)
					notFalsePos1 = j;
				else if (notFalsePos2 < 0) {
					notFalsePos2 = j;
				}
			}
			if (value(ps[j]) != l_Undef) {
				int l = level(var(ps[j]));
				if (l > max_lev) {
					max_lev = l;
				}
			}
		}
		if (!conflicting) {
			assert(notFalsePos1 >= 0);
			if (notFalsePos1 >= 0 && notFalsePos2 >= 0) {
				assert(notFalsePos1 != notFalsePos2);
				if (notFalsePos1 == 1) {
					std::swap(ps[0], ps[notFalsePos2]);
				} else {
					std::swap(ps[0], ps[notFalsePos1]);
					std::swap(ps[1], ps[notFalsePos2]);
				}
			} else {
				std::swap(ps[0], ps[notFalsePos1]);
			}
			assert(value(ps[0])!=l_False);
			if (notFalsePos2 >= 0) {
				assert(value(ps[1])!=l_False);
			}
		}
		
#ifndef NDEBUG
		if (conflicting) {
			for (int j = 0; j < ps.size(); j++)
				assert(value(ps[j])==l_False);
		}
#endif
		if (!permanent && ps.size() > opt_temporary_theory_conflicts) {
			if (tmp_clause == CRef_Undef) {
				tmp_clause = ca.alloc(ps, false);
				tmp_clause_sz = ps.size();
			} else if (tmp_clause_sz < ps.size()) {
				ca[tmp_clause].mark(1);							//is this needed?
				ca.free(tmp_clause);
				tmp_clause = ca.alloc(ps, false);
				tmp_clause_sz = ps.size();
			} else {
				Clause & c = ca[tmp_clause];
				assert(tmp_clause_sz >= ps.size());
				assert(tmp_clause_sz >= c.size());
				c.grow(tmp_clause_sz - c.size());
				for (int i = 0; i < ps.size(); i++) {
					c[i] = ps[i];
				}
				c.shrink(c.size() - ps.size());
			}
			
			confl_out = tmp_clause;
		} else {
			
			CRef cr = ca.alloc(ps, !permanent && !opt_permanent_theory_conflicts);
			ca[cr].setFromTheory(true);
			if (permanent || opt_permanent_theory_conflicts)
				clauses.push(cr);
			else {
				learnts.push(cr);
				if (--learntsize_adjust_cnt <= 0) {
					learntsize_adjust_confl *= learntsize_adjust_inc;
					learntsize_adjust_cnt = (int) learntsize_adjust_confl;
					max_learnts *= learntsize_inc;
					
					if (verbosity >= 1)
						printf("| %9d | %7d %8d %8d | %8d %8d %6.0f | %" PRId64 " removed |\n", (int) conflicts,
								(int) dec_vars - (trail_lim.size() == 0 ? trail.size() : trail_lim[0]), nClauses(),
								(int) clauses_literals, (int) max_learnts, nLearnts(),
								(double) learnts_literals / nLearnts(), stats_removed_clauses);
				}
				
			}
			
			attachClause(cr);
			confl_out = cr;
		}
		//if (!satisfied )
		//	cancelUntil(max_lev);
		
		if (!satisfied && nfalse == ps.size() - 1) {
			cancelUntil(max_lev);
			assert(!conflicting);
			assert(value(ps[0])!=l_False);
			assert(value(ps[1])==l_False);
			if (value(ps[0]) == l_Undef) {
				uncheckedEnqueue(ps[0], confl_out);
			}
		}
		if(!conflicting){
			confl_out=CRef_Undef;
		}
		return !conflicting;
	}
	return true;
}

bool Solver::addDelayedClauses(CRef & conflict_out) {
	conflict_out = CRef_Undef;
	while (clauses_to_add.size() && ok) {
		if (!addConflictClause(clauses_to_add.last(), conflict_out, true)) {
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
lbool Solver::search(int nof_conflicts) {
	assert(ok);
	int backtrack_level;
	int conflictC = 0;
	vec<Lit> learnt_clause;
	bool last_decision_was_theory=false;
	starts++;
	bool using_theory_decisions= opt_decide_theories && drand(random_seed) < opt_random_theory_freq;
	bool using_theory_vsids= opt_decide_theories && opt_theory_order_vsids && drand(random_seed) < opt_random_theory_vsids_freq;

	if(decisionLevel()==0 && initialPropagate && opt_detect_pure_lits && !simplify()){
		return l_False;//if using pure literal detection, and the theories haven't been propagated yet, run simpify
	}

	n_theory_decision_rounds+=using_theory_decisions;
	for (;;) {
		static int iter = 0;
		if (++iter ==  1564) {
			int a = 1;
		}

		propagate: CRef confl = propagate(!disable_theories);
		conflict: if (!okay() || (confl != CRef_Undef)) {
			// CONFLICT
			conflicts++;
			conflictC++;
			if(last_decision_was_theory){
				n_theory_conflicts++;
				consecutive_theory_conflicts++;
				if(opt_theory_conflict_max && consecutive_theory_conflicts>=opt_theory_conflict_max){
					next_theory_decision=conflicts+consecutive_theory_conflicts;
					consecutive_theory_conflicts=0;
				}
			}else{
				consecutive_theory_conflicts=0;
			}
			if (decisionLevel() == 0)
				return l_False;
			learnt_clause.clear();
			analyze(confl, learnt_clause, backtrack_level);
			
			cancelUntil(backtrack_level);
			
			//this is now slightly more complicated, if there are multiple lits implied by the super solver in the current decision level:
			//The learnt clause may not be asserting.
			
			if (learnt_clause.size() == 1) {
				uncheckedEnqueue(learnt_clause[0]);
			} else {
				CRef cr = ca.alloc(learnt_clause, true);
				learnts.push(cr);
				attachClause(cr);
				claBumpActivity(ca[cr]);

				if (value(learnt_clause[0]) == l_Undef) {
					uncheckedEnqueue(learnt_clause[0], cr);
				} else {

					assert(S);
					if(!S){
						fprintf(stderr,"Critical error: bad learnt clause. Aborting\n");
						fflush(stderr);
						exit(3);
						//throw std::runtime_error("Critical error: bad learnt clause. Aborting\n");
					}
					//this is _not_ an asserting clause, its a conflict that must be passed up to the super solver.
					analyzeFinal(cr, lit_Undef, conflict);

					if(theoryConflict>-1){
						theoryBumpActivity(theoryConflict);
						theoryDecayActivity();
						theoryConflict=-1;
					}else if(opt_vsids_solver_as_theory){
						theoryBumpActivity(decisionTheory->getTheoryIndex());
						theoryDecayActivity();
					}

					varDecayActivity();
					claDecayActivity();
					return l_False;
				}
			}
			if(theoryConflict>-1){
				theoryBumpActivity(theoryConflict);
				theoryDecayActivity();
				theoryConflict=-1;
			}else if(opt_vsids_solver_as_theory){
				theoryBumpActivity(decisionTheory->getTheoryIndex());
				theoryDecayActivity();
			}
			varDecayActivity();
			claDecayActivity();

			if (--learntsize_adjust_cnt <= 0) {
				learntsize_adjust_confl *= learntsize_adjust_inc;
				learntsize_adjust_cnt = (int) learntsize_adjust_confl;
				max_learnts *= learntsize_inc;
				
				if (verbosity >= 1)
					printf("|c%9d | %7d %8d %8d | %8d %8d %6.0f | %" PRId64 " removed |\n", (int) conflicts,
							(int) dec_vars - (trail_lim.size() == 0 ? trail.size() : trail_lim[0]), nClauses(),
							(int) clauses_literals, (int) max_learnts, nLearnts(),
							(double) learnts_literals / nLearnts(), stats_removed_clauses);
			}
			
		} else {
			assert(theory_queue.size() == 0 || disable_theories);
			
			if (!addDelayedClauses(confl))
				goto conflict;
			
			if (opt_subsearch == 0 && decisionLevel() < initial_level) {
				return l_Undef;            //give up if we have backtracked past the super solvers decisions
			} else if (opt_subsearch == 2 && S && confl == CRef_Undef) {
				if (super_qhead < S->qhead && !propagateTheory(theory_conflict)) //re-enqueue any super solver decisions that we have backtracked past, and keep going
						{
					if (!addConflictClause(theory_conflict, confl))
						return l_False;
					goto conflict;
				}
			}

			
			// NO CONFLICT
			if ((nof_conflicts >= 0 && conflictC >= nof_conflicts) || !withinBudget()) {
				// Reached bound on number of conflicts:
				progress_estimate = progressEstimate();
				cancelUntil(initial_level);
				return l_Undef;
			}
			
			// Simplify the set of problem clauses:
			if (decisionLevel() == 0 && !simplify())
				return l_False;
			
			if (learnts.size() - nAssigns() >= max_learnts)
				// Reduce the set of learnt clauses:
				reduceDB();
			last_decision_was_theory=false;
			Lit next = lit_Undef;
			while (decisionLevel() < assumptions.size()) {
				// Perform user provided assumption:
				Lit p = assumptions[decisionLevel()];
				assert(p!=lit_Undef);
				assert(var(p)<nVars());
				if (value(p) == l_True) {
					// Dummy decision level:
					newDecisionLevel();
				} else if (value(p) == l_False) {
					analyzeFinal(~p, conflict);
					return l_False;
				} else {
					next = p;
					break;
				}
			}

			if(only_propagate && next ==lit_Undef && decisionLevel() == assumptions.size())
				return l_True;

			//Note: decision level is now added before theories make their decisions, to allow them to decide multiple literals at once.
			newDecisionLevel();



			if (opt_decide_theories && !disable_theories && using_theory_decisions && next == lit_Undef && (opt_theory_conflict_max==0 || conflicts>=next_theory_decision) ) {

				int next_var_priority=INT_MIN;

				//remove decided vars from order heap
				while(!order_heap.empty() && value(order_heap.peekMin())!=l_Undef){
					order_heap.removeMin();
				}
				if(!order_heap.empty()){
					next_var_priority=priority[ order_heap.peekMin()];
				}
				/**
				 * Give the theory solvers a chance to make decisions
				 */
				if(opt_theory_order_vsids && using_theory_vsids){
					while (next==lit_Undef && !theory_order_heap.empty() && theories[theory_order_heap.peekMin()]->getPriority()>=next_var_priority) {
						int theoryID = theory_order_heap.peekMin();
						if(opt_vsids_both && next_var_priority==theories[theory_order_heap.peekMin()]->getPriority()){
							//give the main solver a chance to make a decision, if it has a variable with higher activity
							if (!order_heap.empty()) {
								Var v = order_heap.peekMin();
								assert(value(v)==l_Undef);//because assigned lits should have been removed from queue above.
								if(value(v)==l_Undef && activity[v]>theories[theoryID]->getActivity()){
									stats_solver_preempted_decisions++;
									break;
								}
							}

						}

						assert(theories[theoryID]->supportsDecisions());
						next = theories[theoryID]->decideTheory();
						if(next==lit_Undef){
							theory_order_heap.removeMin();
							if(decisionLevel()>0)
								theory_decision_trail.push({theoryID,decisionLevel()});
						}
					}
				}else{
					for (int i = 0; i < decidable_theories.size() && next == lit_Undef; i++) {
						if (decidable_theories[i]->getPriority()>=next_var_priority) {
							next = decidable_theories[i]->decideTheory();
						}
					}
				}

				if (next!=lit_Undef){
					stats_theory_decisions++;
					if(theoryDecision !=lit_Undef && var(next)==var(theoryDecision)){
						assigns[var(theoryDecision)]=l_Undef;
					}
					last_decision_was_theory=true;
				}
			}
			
			if (next == lit_Undef) {
				// New variable decision:

				next = pickBranchLit();
				// int p = priority[var(next)];
				if(opt_verb>2){
					printf("solver decision\n");
				}
				if (next == lit_Undef) {
					
					//solve theories if this solver is completely assigned
					if(!disable_theories){
						for (int i = 0; i < theories.size(); i++) {
							if (opt_subsearch == 3 && track_min_level < initial_level)
								continue; //Disable attempting to solve sub-solvers if we've backtracked past the super solver's decision level

							if (!theories[i]->solveTheory(theory_conflict)) {
								if (!addConflictClause(theory_conflict, confl)) {
									goto conflict;
								} else {
									goto propagate;
								}
							}
							//If propagating one of the sub theories caused this solver to backtrack, then go back to propagation
							if (qhead < trail.size() || hasNextDecision())
								goto propagate;
						}
					}
					// Model found:
					return l_True;
				}
			}
			//last_dec = var(next);
			// Increase decision level and enqueue 'next'
			assert(next!=lit_Undef);
			//if(next!=lit_Error)//lit_Error is used to signify a decision that has no literal in the SAT solver (some theories may support this)

			enqueue(next);//not unchecked enqueue, because a theory solver _may_ have assigned this literal while making a decision
		}
	}
	//Unreachable
	return l_Undef;
}

double Solver::progressEstimate() const {
	double progress = 0;
	double F = 1.0 / nVars();
	
	for (int i = 0; i <= decisionLevel(); i++) {
		int beg = i == 0 ? 0 : trail_lim[i - 1];
		int end = i == decisionLevel() ? trail.size() : trail_lim[i];
		progress += pow(F, i) * (end - beg);
	}
	
	return progress / nVars();
}

/*
 Finite subsequences of the Luby-sequence:

 0: 1
 1: 1 1 2
 2: 1 1 2 1 1 2 4
 3: 1 1 2 1 1 2 4 1 1 2 1 1 2 4 8
 ...


 */

static double luby(double y, int x) {
	assert(x>=0);
	// Find the finite subsequence that contains index 'x', and the
	// size of that subsequence:
	int size, seq;
	for (size = 1, seq = 0; size < x + 1; seq++, size = 2 * size + 1);

	assert(size>x);
	while (size - 1 != x) {
		size = (size - 1) >> 1;
		seq--;
		//According to Coverity: size can be zero at this line, leading to a mod by zero...
		//However, since size must be >= x+1 above, and always > x in this loop, and x is positive, this is safe.
		x = x % size;
	}
	
	return pow(y, seq);
}
bool Solver::propagateAssignment(const vec<Lit> & assumps){
	only_propagate=true;
	budgetOff();
	assumps.copyTo(assumptions);

	lbool val = solve_();
	only_propagate=false;
	return val==l_True;
}
// NOTE: assumptions passed in member-variable 'assumptions'.
lbool Solver::solve_() {
	clearInterrupt();
	cancelUntil(0);
	model.clear();
	conflict.clear();
	if (!ok)
		return l_False;
	
	solves++;
	
	max_learnts = nClauses() * learntsize_factor;
	learntsize_adjust_confl = learntsize_adjust_start_confl;
	learntsize_adjust_cnt = (int) learntsize_adjust_confl;
	lbool status = l_Undef;
	
	if (verbosity >= 1 && !printed_header) {
		//on repeated calls, don't print the header again
		printed_header = true;
		printf("============================[ Search Statistics ]==============================\n");
		printf("| Conflicts |          ORIGINAL         |          LEARNT          | Progress |\n");
		printf("|           |    Vars  Clauses Literals |    Limit  Clauses Lit/Cl |          |\n");
		printf("===============================================================================\n");
	}
	initial_level = 0;
	track_min_level = 0;
	// Search:
	int curr_restarts = 0;
	while (status == l_Undef) {
		double rest_base = luby_restart ? luby(restart_inc, curr_restarts) : pow(restart_inc, curr_restarts);

		if (opt_rnd_phase) {
			for (int i = 0; i < nVars(); i++)
				polarity[i] = irand(random_seed, 1);
		}
		if (opt_decide_theories && !opt_theory_order_vsids && opt_randomomize_theory_order) {
			randomShuffle(random_seed, decidable_theories);
		}

		status = search(rest_base * restart_first);
		if (verbosity >= 1) {
			printf("|r%9d | %7d %8d %8d | %8d %8d %6.0f | %" PRId64 " removed |\n", (int) conflicts,
				   (int) dec_vars - (trail_lim.size() == 0 ? trail.size() : trail_lim[0]), nClauses(),
				   (int) clauses_literals, (int) max_learnts, nLearnts(),
				   (double) learnts_literals / nLearnts(), stats_removed_clauses);
		}
		if (!withinBudget())
			break;
		curr_restarts++;
		if (opt_rnd_restart && status == l_Undef) {
			
			for (int i = 0; i < nVars(); i++) {
				activity[i] = drand(random_seed) * 0.00001;
			}
			rebuildOrderHeap();
		}


		if(opt_theory_order_vsids && opt_rnd_theory_restart  && status == l_Undef){
			for(int i = 0;i<theories.size();i++){
				theories[i]->setActivity(  drand(random_seed) * 0.00001);
			}
			rebuildTheoryOrderHeap();
		}
	}
	
	if (status == l_True) {
		// Extend & copy model:
		model.growTo(nVars());
		for (int i = 0; i < nVars(); i++)
			model[i] = value(i);
		
		for(Lit l:assumptions){
			if(value(l)!=l_True){
				throw std::runtime_error("Model is inconsistent with assumptions!");
			}
		}

		if (opt_check_solution && theories.size() && !disable_theories && !only_propagate) {
			if(opt_verb>0){
				printf("Checking witnesses...\n");
			}
			double check_start=rtime(1);
			for(Theory * t:theories){
			
				if (!t->check_solved()) {
					fprintf(stderr, "Error! Solution doesn't satisfy theory properties!\n");
					fflush(stderr);
					exit(4);
				}
			}
			stats_solution_checking_time+=rtime(1)-check_start;
		}
		if(opt_write_learnt_clauses && opt_debug_model){
			//write the model out
			fprintf(opt_write_learnt_clauses,"s SATISFIABLE\n");
			fprintf(opt_write_learnt_clauses,"v ");
			for (int i = 0;i<nVars();i++){
				Lit l = mkLit(i,false);
				if(value(i)==l_True){
					fprintf(opt_write_learnt_clauses," %d",dimacs(unmap(l)));
				}else if(value(i)==l_False){
					fprintf(opt_write_learnt_clauses," %d",dimacs(~unmap(l)));
				}else{
					fprintf(opt_write_learnt_clauses," %d",dimacs(unmap(l)));//optional
				}
			}
			fprintf(opt_write_learnt_clauses,"\n");
			if(!disable_theories){
			for(Theory * t:theories){
				std::stringstream ss;
				t->writeTheoryWitness(ss);
				fprintf(opt_write_learnt_clauses,"%s",ss.str().c_str());
			}
			}

		}
		
	} else if (status == l_False && conflict.size() == 0) {
		ok = false;
	} else if (status == l_False) {
		assert(ok);
	}
	only_propagate=false;
	assumptions.clear();
	return status;
}

bool Solver::solveTheory(vec<Lit> & conflict_out) {
	initial_level = decisionLevel();
	track_min_level = initial_level;
	lbool status = l_Undef;
	// Search:
	int curr_restarts = 0;
	conflict.clear();
	conflict_out.clear();
	while (status == l_Undef) {
		double rest_base = luby_restart ? luby(restart_inc, curr_restarts) : pow(restart_inc, curr_restarts);
		status = search(rest_base * restart_first);
		if (!withinBudget() || (opt_subsearch == 0 && track_min_level < initial_level))
			break;
		curr_restarts++;
	}
	cancelUntil(track_min_level);
	if (track_min_level < initial_level) {
		S->cancelUntil(track_min_level);
	}
	if (!ok) {
		return false;
	}
	if (conflict.size()) {
		toSuper(conflict.toVec(), conflict_out);
		interpolant.push();
		conflict_out.copyTo(interpolant.last());    //Add this clause into the interpolant vector
		return false;
	}
	
	return propagateTheory(conflict_out);
}

//=================================================================================================
// Writing CNF to DIMACS:
// 
// FIXME: this needs to be rewritten completely.

Lit Solver::unmap(Lit l)  {
	if(varRemap){
		return varRemap->unmap(l);
	}else{
		return l;
	}
}

static Var mapVar(Var x, vec<Var>& map, Var& max) {
	if (map.size() <= x || map[x] == -1) {
		map.growTo(x + 1, -1);
		map[x] = max++;
	}
	return map[x];
}

void Solver::toDimacs(FILE* f, Clause& c, vec<Var>& map, Var& max) {
	if (satisfied(c))
		return;
	
	for (int i = 0; i < c.size(); i++)
		if (value(c[i]) != l_False)
			fprintf(f, "%s%d ", sign(c[i]) ? "-" : "", mapVar(var(c[i]), map, max) + 1);
	fprintf(f, "0\n");
}

void Solver::toDimacs(const char *file, const vec<Lit>& assumps) {
	FILE* f = fopen(file, "wr");
	if (f == NULL)
		fprintf(stderr, "could not open file %s\n", file), exit(1);
	toDimacs(f, assumps);
	fclose(f);
}

void Solver::toDimacs(FILE* f, const vec<Lit>& assumps) {
	// Handle case when solver is in contradictory state:
	if (!ok) {
		fprintf(f, "p cnf 1 2\n1 0\n-1 0\n");
		return;
	}
	
	vec<Var> map;
	Var max = 0;
	
	// Cannot use removeClauses here because it is not safe
	// to deallocate them at this point. Could be improved.
	int cnt = 0;
	for (int i = 0; i < clauses.size(); i++)
		if (!satisfied(ca[clauses[i]]))
			cnt++;
	
	for (int i = 0; i < clauses.size(); i++)
		if (!satisfied(ca[clauses[i]])) {
			Clause& c = ca[clauses[i]];
			for (int j = 0; j < c.size(); j++)
				if (value(c[j]) != l_False)
					mapVar(var(c[j]), map, max);
		}
	
	// Assumptions are added as unit clauses:
	cnt += assumptions.size();
	
	fprintf(f, "p cnf %d %d\n", max, cnt);
	
	for (int i = 0; i < assumptions.size(); i++) {
		assert(value(assumptions[i]) != l_False);
		fprintf(f, "%s%d 0\n", sign(assumptions[i]) ? "-" : "", mapVar(var(assumptions[i]), map, max) + 1);
	}
	
	for (int i = 0; i < clauses.size(); i++)
		toDimacs(f, ca[clauses[i]], map, max);
	
	if (verbosity > 0)
		printf("Wrote %d clauses with %d variables.\n", cnt, max);
}

//=================================================================================================
// Garbage Collection methods:

void Solver::relocAll(ClauseAllocator& to) {
	if (tmp_clause != CRef_Undef) {
		ca[tmp_clause].mark(1);    //is this needed?
		ca.free(tmp_clause);
		tmp_clause = CRef_Undef;
	}
	//Re-allocate the 'theory markers'
	vec<int> marker_theory_tmp;
	
	/*	for(int i = 0;i<markers.size();i++){
	 CRef old_cr = markers[i];
	 assert(old_cr!=CRef_Undef);
	 int old_theory = getTheory(old_cr);
	 assert(old_theory>=0);
	 CRef cr=to.makeMarkerReference();
	 assert(cr==markers[i]);//these should be identical in the current implementation
	 markers[i]=cr;
	 int index = CRef_Undef-cr - 1;
	 marker_theory_tmp.growTo(index+1,-1 );
	 marker_theory_tmp[index] = old_theory;

	 assert( marker_theory_tmp[index] == old_theory);

	 assert(markers[ marker_theory_tmp[index]] == cr);
	 }
	 marker_theory_tmp.copyTo(marker_theory);*/
	// All watchers:
	//
	// for (int i = 0; i < watches.size(); i++)
	watches.cleanAll();
	for (int v = 0; v < nVars(); v++){

		for (int s = 0; s < 2; s++) {
			Lit p = mkLit(v, s);
			// printf(" >>> RELOCING: %s%d\n", sign(p)?"-":"", var(p)+1);
			vec<Watcher>& ws = watches[p];
			for (int j = 0; j < ws.size(); j++)
				ca.reloc(ws[j].cref, to);
		}
	}
	// All reasons:
	//
	for (int i = 0; i < trail.size(); i++) {
		Var v = var(trail[i]);
		
		if (ca.isClause(reason(v)) && (ca[reason(v)].reloced() || locked(ca[reason(v)])))
			ca.reloc(vardata[v].reason, to);
	}
	
	// All learnt:
	//
	for (int i = 0; i < learnts.size(); i++)
		ca.reloc(learnts[i], to);
	
	// All original:
	//
	for (int i = 0; i < clauses.size(); i++)
		ca.reloc(clauses[i], to);
}

void Solver::garbageCollect() {
	// Initialize the next region to a size corresponding to the estimated utilization degree. This
	// is not precise but should avoid some unnecessary reallocations for the new region:
	ClauseAllocator to(ca.size() - ca.wasted());
	relocAll(to);
	if (verbosity >= 2)
		printf("|  Garbage collection:   %12d bytes => %12d bytes             |\n",
				ca.size() * ClauseAllocator::Unit_Size, to.size() * ClauseAllocator::Unit_Size);
	to.moveTo(ca);
}

