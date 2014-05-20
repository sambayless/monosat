/****************************************************************************************[Solver.h]
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

#ifndef Minisat_Solver_h
#define Minisat_Solver_h

#include "mtl/Vec.h"
#include "mtl/Heap.h"
#include "mtl/Alg.h"
#include "mtl/Rnd.h"
#include "utils/Options.h"
#include "core/SolverTypes.h"
#include "core/Theory.h"
#include "core/Config.h"
namespace Minisat {

//=================================================================================================
// Solver -- the main class:

class Solver:public Theory {
public:
	friend class Theory;
	friend class GraphTheorySolver;
#ifdef DEBUG_SOLVER
	Solver * dbg_solver;
#endif
    // Constructor/Destructor:
    //
    Solver();
    virtual ~Solver();

    // Problem specification:
    //
    virtual Var     newVar    (bool polarity = true, bool dvar = true); // Add a new variable with parameters specifying variable mode.

    virtual bool    addClause (const vec<Lit>& ps);                     // Add a clause to the solver.
    virtual bool    addEmptyClause();                                   // Add the empty clause, making the solver contradictory.
    virtual bool    addClause (Lit p);                                  // Add a unit clause to the solver.
    virtual bool    addClause (Lit p, Lit q);                           // Add a binary clause to the solver.
    virtual bool    addClause (Lit p, Lit q, Lit r);                    // Add a ternary clause to the solver.
    bool    addClause_(      vec<Lit>& ps);                     // Add a clause to the solver without making superflous internal copy. Will

    void setDecisionPriority(Var v,unsigned int p){
    	priority[v]=p;
    }

   //Theory interface
    void addTheory(Theory*t){
    	theories.push(t);
    	theory_queue.capacity(theories.size());
    	in_theory_queue.push(false);
    	t->setTheoryIndex(theories.size()-1);
    	cancelUntil(0);
    	resetInitialPropagation();
    }



    //Call to force at least one round of propagation to each theory solver at the next solve() call
    void resetInitialPropagation(){
    	if(!initialPropagate){
			initialPropagate=true;//to force propagation to occur at least once to the theory solvers
			if(S){
				S->resetInitialPropagation();
			}
    	}
    }

    void detatchTheory(Theory*t){
    	assert(decisionLevel()==0);
    	assert(theory_queue.size()==0);
    	int i,j =0;
    	for(i = 0;i<theories.size();i++){
    		if(theories[i]==t){
    			int k,l=0;
    			for(k=0;k<markers.size();k++){
    				CRef cr = markers[k];
    			  	int index = CRef_Undef-cr-1;
    			  	int theory = marker_theory[index];
    				if(theory==i){
    					marker_theory[index]=-1;
    				}else{
    					markers[l++]=cr;
    					if(marker_theory[index]>i){
    						marker_theory[index]--;
    					}
    				}
    			}
    			markers.shrink(k-l);
    		}else{
    			theories[j++]=theories[i];
    			in_theory_queue[j++]=in_theory_queue[i];
    		}
    	}
    	theories.shrink(i-j);
    	in_theory_queue.shrink(i-j);

    	cancelUntil(0);
    }

    int getTheoryIndex(){
    	return theory_index;
    }
    void setTheoryIndex(int id){
    	theory_index=id;
    }

    //Generate a new, unique `temporary value' for explaining conflicts
    CRef newReasonMarker(Theory * forTheory){
    	markers.push(ca.makeMarkerReference());


    	int marker_num = CRef_Undef-markers.last()-1;
    	marker_theory.growTo(marker_num+1,-1);

  
    	//this could be done more efficiently
    	for(int i = 0;i<theories.size();i++){
    		if(theories[i]==forTheory){
    			marker_theory[marker_num]=i;
    			break;
    		}
    	}
    	assert(marker_theory[marker_num]>=0);

    	return markers.last();
    }

    void printStats(int detail_level=0){
		double cpu_time = cpuTime();
		double mem_used = memUsedPeak();

	    printf("restarts              : %" PRIu64 "\n", starts);
		printf("conflicts             : %-12" PRIu64 "   (%.0f /sec, %d learnts, %d removed)\n", conflicts   , conflicts   /cpu_time, learnts.size(),stats_removed_clauses);
		printf("decisions             : %-12" PRIu64 "   (%4.2f %% random) (%.0f /sec)\n", decisions, (float)rnd_decisions*100 / (float)decisions, decisions   /cpu_time);
		printf("propagations          : %-12" PRIu64 "   (%.0f /sec)\n", propagations, propagations/cpu_time);
		printf("conflict literals     : %-12" PRIu64 "   (%4.2f %% deleted)\n", tot_literals, (max_literals - tot_literals)*100 / (double)max_literals);
		if(opt_detect_pure_theory_lits){
			printf("pure literals     : %d (%d theory lits) (%d rounds, %f time)\n",stats_pure_lits, stats_pure_theory_lits,pure_literal_detections,stats_pure_lit_time);
		}
		for(int i = 0;i<theories.size();i++){
			theories[i]->printStats(detail_level);
		}
    }

    bool isTheoryCause(CRef cr){
    	return cr != CRef_Undef && !ca.isClause(cr);
    }

    int getTheory(CRef cr){
    	assert(isTheoryCause(cr));
    	// UINT32_MAX-cr - 1;
    	int marker = CRef_Undef-cr-1;
    	return marker_theory[marker];
    }

    bool	 hasTheory(Var v){
    	return theory_vars[v].isTheoryVar;
    }
    bool	hasTheory(Lit l){
    	return theory_vars[var(l)].isTheoryVar;
    }
    int		getTheoryID(Var v){
    	return theory_vars[v].theory-1;
    }
    int		getTheoryID(Lit l){
    	return theory_vars[var(l)].theory-1;
    }
    Var		getTheoryVar(Var v){
    	assert(hasTheory(v));
    	return (Var) theory_vars[v].theory_var;
    }

    //Translate a literal into its corresponding theory literal (if it has a theory literal)
    Lit		getTheoryLit(Lit l){
    	assert(hasTheory(l));
    	return mkLit(getTheoryVar(var(l)),sign(l));
    }

    //Connect a variable in the SAT solver to a variable in a theory.
    virtual void setTheoryVar(Var solverVar, int theory, Var theoryVar){
    	if(hasTheory(solverVar)){
    		fprintf(stderr,"Variable %d is used for multiple atoms. Theory variables may not be re-used! Aborting.\n", solverVar+1);
    		exit(1);
    	}
    	assert(!hasTheory(solverVar));
    	theory_vars[solverVar].theory=theory+1;
    	theory_vars[solverVar].theory_var=theoryVar;
    	all_theory_vars.push(solverVar);
    	assert(hasTheory(solverVar));
    	assert(getTheoryID(solverVar)==theory);
    	assert(getTheoryVar(solverVar)==theoryVar);
    	assert(decisionLevel()==0);
    	if(value(solverVar)!=l_Undef)
    		initialPropagate=true;

    }

    void remapTheoryVar(Var solverVar, Var newTheoryVar){
    	assert(hasTheory(solverVar));
    	theory_vars[solverVar].theory_var=newTheoryVar;
    }

    CRef constructReason(Lit p){
    	CRef cr = reason(var(p));
    	assert(isTheoryCause(cr));
    	assert(!ca.isClause(cr));
    	assert(cr!=CRef_Undef);
    	int t = getTheory(cr);
    	assert(hasTheory(p));
    	theory_reason.clear();
    	theories[t]->buildReason(getTheoryLit(p),theory_reason);
    	/*CRef reason = ca.alloc(theory_reason, !opt_permanent_theory_conflicts);
    	if(opt_permanent_theory_conflicts)
			clauses.push(reason);
		else
			learnts.push(reason);*/
    	assert(theory_reason[0]==p); assert(value(p)==l_True);
#ifdef DEBUG_SOLVER
    	//assert all the other reasons in this cause are earlier on the trail than p...
    	static vec<bool> marks;
    	marks.clear();
    	marks.growTo(nVars());
    	for(int i = 0;i<trail.size() && var(trail[i])!=var(p);i++){
    		marks[var(trail[i])]=true;
    	}
    	for(int i = 1;i<theory_reason.size();i++){
    		assert(marks[var(theory_reason[i])]);
    	}
#endif

    	CRef reason  = attachClauseSafe(theory_reason);
		vardata[var(p)]=mkVarData(reason,level(var(p)));
		return reason;
    }

    //Attach this solver to the next, connecting literals in the super solver to literals in this solver
    //Super_vars and local_vars _must_ be consecutive sets of variables in the super solver.
    void attachTo(Solver * super, vec<Var> & super_vars, vec<Var> & local_vars){
    	assert(decisionLevel()==0);
    	super->addTheory(this);
    	S=super;
    	assert(super_vars.size()==local_vars.size());
    	super_qhead=0;
    	local_qhead=0;
		super_offset=super_vars[0]-local_vars[0];
		min_super=super_vars[0];
		max_super=super_vars.last();
		min_local=local_vars[0];
		max_local=local_vars.last();
		cause_marker=super->newReasonMarker(this);
    }
                                                      // change the passed vector 'ps'.

    // Solving:
    //
    bool    simplify     ();                        // Removes already satisfied clauses.
    virtual bool    solve        (const vec<Lit>& assumps); // Search for a model that respects a given set of assumptions.
    virtual lbool   solveLimited (const vec<Lit>& assumps); // Search for a model that respects a given set of assumptions (With resource constraints).
    virtual bool    solve        ();                        // Search without assumptions.
    virtual bool    solve        (Lit p);                   // Search for a model that respects a single assumption.
    virtual bool    solve        (Lit p, Lit q);            // Search for a model that respects two assumptions.
    virtual bool    solve        (Lit p, Lit q, Lit r);     // Search for a model that respects three assumptions.
    bool    okay         () const;                  // FALSE means solver is in a conflicting state


    void    toDimacs     (FILE* f, const vec<Lit>& assumps);            // Write CNF to file in DIMACS-format.
    void    toDimacs     (const char *file, const vec<Lit>& assumps);
    void    toDimacs     (FILE* f, Clause& c, vec<Var>& map, Var& max);

    // Convenience versions of 'toDimacs()':
    void    toDimacs     (const char* file);
    void    toDimacs     (const char* file, Lit p);
    void    toDimacs     (const char* file, Lit p, Lit q);
    void    toDimacs     (const char* file, Lit p, Lit q, Lit r);
    
    // Variable mode:
    // 
    void    setPolarity    (Var v, bool b); // Declare which polarity the decision heuristic should use for a variable. Requires mode 'polarity_user'.
    void    setDecisionVar (Var v, bool b); // Declare if a variable should be eligible for selection in the decision heuristic.

    // Read state:
    //
    lbool   value      (Var x) const;       // The current value of a variable.
    lbool   value      (Lit p) const;       // The current value of a literal.
    lbool   modelValue (Var x) const;       // The value of a variable in the last model. The last call to solve must have been satisfiable.
    lbool   modelValue (Lit p) const;       // The value of a literal in the last model. The last call to solve must have been satisfiable.
    int     nAssigns   ()      const;       // The current number of assigned literals.
    int     nClauses   ()      const;       // The current number of original clauses.
    int     nLearnts   ()      const;       // The current number of learnt clauses.
    int     nVars      ()      const;       // The current number of variables.
    int		nUnassignedVars()  const;		// The number of variables left to assign. Does not include variables that are not decision vars!
    int     nFreeVars  ()      const;		// The number of non-unit variables

    // Resource contraints:
    //
    void    setConfBudget(int64_t x);
    void    setPropBudget(int64_t x);
    void    budgetOff();
    void    interrupt();          // Trigger a (potentially asynchronous) interruption of the solver.
    void    clearInterrupt();     // Clear interrupt indicator flag.

    // Memory managment:
    //
    virtual void garbageCollect();
    void    checkGarbage(double gf);
    void    checkGarbage();

    // Extra results: (read-only member variable)
    //
    vec<lbool> model;             // If problem is satisfiable, this vector contains the model (if any).
    vec<Lit>   conflict;          // If problem is unsatisfiable (possibly under assumptions),
                                  // this vector represent the final conflict clause expressed in the assumptions.
    vec<vec<Lit> > interpolant;  //This vector represents an interpolant between this module and its super solver ('S'), if it is attached to such a solver and the instance is UNSAT.
    							// Variables in each clause in the interpolant vector are in the super solver's variable space, not the subsolver's.

    vec<Lit>  theory_reason;
    vec<Lit> theory_conflict;
    vec<Theory*> theories;
    vec<Theory*> decidable_theories;
    vec<Var> all_theory_vars;
    struct LitCount{
    	char occurs:1;
    	char seen:1;
    	LitCount():occurs(1),seen(0){ 	}
    };
    vec<LitCount> lit_counts;
    vec<CRef> markers;//a set of special clauses that can be recognized as pointers to theories
    vec<int> marker_theory;
    int theory_index;
    Solver * S;//super solver
    bool initialPropagate;//to force propagation to occur at least once to the theory solvers
    int super_qhead;
    int local_qhead;
    CRef cause_marker;
    int track_min_level;
    int initial_level;
    vec<int> theory_queue;
    vec<bool> in_theory_queue;
    int max_decision_var;
    CRef tmp_clause=CRef_Undef;
    int tmp_clause_sz=0;
    Var max_super;
    Var min_super;
    Var min_local;
    Var max_local;
    int super_offset;



    // Mode of operation:
    //
    int       verbosity;
    double    var_decay;
    double    clause_decay;
    double    random_var_freq;
    double    random_seed;
    bool      luby_restart;
    int       ccmin_mode;         // Controls conflict clause minimization (0=none, 1=basic, 2=deep).
    int       phase_saving;       // Controls the level of phase saving (0=none, 1=limited, 2=full).
    bool      rnd_pol;            // Use random polarities for branching heuristics.
    bool      rnd_init_act;       // Initialize variable activities with a small random value.
    double    garbage_frac;       // The fraction of wasted memory allowed before a garbage collection is triggered.

    int       restart_first;      // The initial restart limit.                                                                (default 100)
    double    restart_inc;        // The factor with which the restart limit is multiplied in each restart.                    (default 1.5)
    double    learntsize_factor;  // The intitial limit for learnt clauses is a factor of the original clauses.                (default 1 / 3)
    double    learntsize_inc;     // The limit for learnt clauses is multiplied with this factor each restart.                 (default 1.1)

    int       learntsize_adjust_start_confl;
    double    learntsize_adjust_inc;

    // Statistics: (read-only member variable)
    //
    uint64_t solves, starts, decisions, rnd_decisions, propagations, conflicts,stats_pure_lits,stats_pure_theory_lits,pure_literal_detections,stats_removed_clauses;
    uint64_t dec_vars, clauses_literals, learnts_literals, max_literals, tot_literals;
    double stats_pure_lit_time;

    Var last_dec;
protected:

    // Helper structures:
    //
    struct VarData { CRef reason; int level; };
    static inline VarData mkVarData(CRef cr, int l){ VarData d = {cr, l}; return d; }

    struct Watcher {
        CRef cref;
        Lit  blocker;
        Watcher(CRef cr, Lit p) : cref(cr), blocker(p) {}
        bool operator==(const Watcher& w) const { return cref == w.cref; }
        bool operator!=(const Watcher& w) const { return cref != w.cref; }
    };

    struct WatcherDeleted
    {
        const ClauseAllocator& ca;
        WatcherDeleted(const ClauseAllocator& _ca) : ca(_ca) {}
        bool operator()(const Watcher& w) const { return ca[w.cref].mark() == 1; }
    };

    struct VarOrderLt {
        const vec<double>&  activity;
        const vec<int> &  priority;
        bool operator () (Var x, Var y) const {
        	if (priority[x] == priority[y])
        		return activity[x] > activity[y];
        	else{
        		return priority[x]>priority[y];
        	}
        }
        VarOrderLt(const vec<double>&  act,const vec<int>&  pri) : activity(act),priority(pri) { }
    };


    struct TheoryData{
    	union{
    		struct{
			unsigned int theory:11;
			unsigned int theory_var:21;
    		};
    		unsigned int isTheoryVar; //true if non-zero - this property is ensured by adding 1 to theory_var
    	};
    	TheoryData():isTheoryVar(0){}
    	TheoryData(unsigned int theory,unsigned int theory_lit):theory(theory),theory_var(theory_lit){

    	}
    };

    // Solver state:
    //
    bool                ok;               // If FALSE, the constraints are already unsatisfiable. No part of the solver state may be used!
    vec<CRef>           clauses;          // List of problem clauses.
    vec<CRef>           learnts;          // List of learnt clauses.
    double              cla_inc;          // Amount to bump next clause with.
    vec<double>         activity;         // A heuristic measurement of the activity of a variable.
    double              var_inc;          // Amount to bump next variable with.
    OccLists<Lit, vec<Watcher>, WatcherDeleted>
                        watches;          // 'watches[lit]' is a list of constraints watching 'lit' (will go there if literal becomes true).
    vec<lbool>          assigns;          // The current assignments.
    vec<char>           polarity;         // The preferred polarity of each variable.
    vec<char>           decision;         // Declares if a variable is eligible for selection in the decision heuristic.
    vec<int>			priority;		  // Static, lexicographic heuristic
    vec<TheoryData>     theory_vars;
    vec<Lit>            trail;            // Assignment stack; stores all assigments made in the order they were made.
    vec<int>            trail_lim;        // Separator indices for different decision levels in 'trail'.
    vec<VarData>        vardata;          // Stores reason and level for each variable.
    int                 qhead;            // Head of queue (as index into the trail -- no more explicit propagation queue in MiniSat).
    int                 simpDB_assigns;   // Number of top-level assignments since last execution of 'simplify()'.
    int64_t             simpDB_props;     // Remaining number of propagations that must be made before next execution of 'simplify()'.
    vec<Lit>            assumptions;      // Current set of assumptions provided to solve by the user.
    Heap<VarOrderLt>    order_heap;       // A priority queue of variables ordered with respect to the variable activity.
    double              progress_estimate;// Set by 'search()'.
    bool                remove_satisfied; // Indicates whether possibly inefficient linear scan for satisfied clauses should be performed in 'simplify'.

    ClauseAllocator     ca;

    // Temporaries (to reduce allocation overhead). Each variable is prefixed by the method in which it is
    // used, exept 'seen' wich is used in several places.
    //
    vec<char>           seen;
    vec<Lit>            analyze_stack;
    vec<Lit>            analyze_toclear;
    vec<Lit>            add_tmp;

    vec<vec<Lit>> 		clauses_to_add;

    double              max_learnts;
    double              learntsize_adjust_confl;
    int                 learntsize_adjust_cnt;

    // Resource contraints:
    //
    int64_t             conflict_budget;    // -1 means no budget.
    int64_t             propagation_budget; // -1 means no budget.
    bool                asynch_interrupt;

    // Main internal methods:
    //
    void     insertVarOrder   (Var x);                                                 // Insert a variable in the decision order priority queue.
    Lit      pickBranchLit    ();                                                      // Return the next decision variable.
public:
    void     newDecisionLevel ();                                                      // Begins a new decision level.
    void     uncheckedEnqueue (Lit p, CRef from = CRef_Undef);                         // Enqueue a literal. Assumes value of literal is undefined.
    bool     enqueue          (Lit p, CRef from = CRef_Undef);                         // Test if fact 'p' contradicts current state, enqueue otherwise.

protected:
    CRef     propagate        (bool propagate_theories=true);                                                      // Perform unit propagation. Returns possibly conflicting clause.
	void 	enqueueTheory(Lit l);
    bool 	propagateTheory(vec<Lit> & conflict);
    bool 	solveTheory(vec<Lit> & conflict_out);


    void 	buildReason(Lit p, vec<Lit> & reason);
    void backtrackUntil(int level);
    //Add a clause to the clause database safely, even if the solver is in the middle of search, propagation, or clause analysis.
    //(In reality, the clause will be added to the database sometime later)
    void 	addClauseSafely(vec<Lit> & ps){
    	clauses_to_add.push();
    	ps.copyTo(clauses_to_add.last());
    }
public:
    void     cancelUntil      (int level);                                             // Backtrack until a certain level.
protected:
    void     analyze          (CRef confl, vec<Lit>& out_learnt, int& out_btlevel);    // (bt = backtrack)
    void 	analyzeFinal(CRef confl, Lit skip_lit, vec<Lit>& out_conflict);
    void     analyzeFinal     (Lit p, vec<Lit>& out_conflict);                         // COULD THIS BE IMPLEMENTED BY THE ORDINARIY "analyze" BY SOME REASONABLE GENERALIZATION?
    bool     litRedundant     (Lit p, uint32_t abstract_levels);                       // (helper method for 'analyze()')
    lbool    search           (int nof_conflicts);                                     // Search for a given number of conflicts.
    lbool    solve_           ();                                                      // Main solve method (assumptions given in 'assumptions').
    void     reduceDB         ();                                                      // Reduce the set of learnt clauses.
    void     removeSatisfied  (vec<CRef>& cs);                                         // Shrink 'cs' to contain only non-satisfied clauses.
    void     rebuildOrderHeap ();

    // Maintaining Variable/Clause activity:
    //
    void     varDecayActivity ();                      // Decay all variables with the specified factor. Implemented by increasing the 'bump' value instead.
    void     varBumpActivity  (Var v, double inc);     // Increase a variable with the current 'bump' value.
    void     varBumpActivity  (Var v);                 // Increase a variable with the current 'bump' value.
    void     claDecayActivity ();                      // Decay all clauses with the specified factor. Implemented by increasing the 'bump' value instead.
    void     claBumpActivity  (Clause& c);             // Increase a clause with the current 'bump' value.

    // Operations on clauses:
    //
    CRef	 attachClauseSafe(vec<Lit> & ps);
    void     attachClause     (CRef cr);               // Attach a clause to watcher lists.
    void     detachClause     (CRef cr, bool strict = false); // Detach a clause to watcher lists.
    void     removeClause     (CRef cr);               // Detach and free a clause.
    bool     locked           (const Clause& c) const; // Returns TRUE if a clause is a reason for some implication in the current state.
    bool     satisfied        (const Clause& c) const; // Returns TRUE if a clause is satisfied in the current state.

    void     relocAll         (ClauseAllocator& to);
public:
    // Misc:
    //
    int      decisionLevel    ()      const; // Gives the current decisionlevel.
    uint32_t abstractLevel    (Var x) const; // Used to represent an abstraction of sets of decision levels.

    CRef     reason           (Var x) const;
    int      level            (Var x) const;
    double   progressEstimate ()      const; // DELETE THIS ?? IT'S NOT VERY USEFUL ...
private:
    bool     withinBudget     ()      const;
    bool 	addConflictClause(vec<Lit> & theory_conflict,CRef & confl_out, bool permanent=false);

    bool	addDelayedClauses(CRef & conflict);
    // Static helpers:
    //
    inline void toSuper(vec<Lit> & from, vec<Lit> & to){

    	to.growTo(from.size());
    	to.shrink(to.size()- from.size());
		for(int i = 0;i<from.size();i++){
			Lit l = from[i];
			assert(local_interface(l));
			to[i]=mkLit(var(l)+super_offset,sign(l));
			assert(super_interface(to[i]));
		}
    }

    inline Var toSuper(Var p){
        	assert(p<nVars());
        	assert(local_interface(p));
        	Var v=  p+super_offset;
        	assert(v<S->nVars());
			assert(super_interface(v));
            return v;
        }

        inline Lit toSuper(Lit p){
        	assert(var(p)<nVars());
        	assert(local_interface(p));

        	Lit l=  mkLit(var(p)+super_offset,sign(p));
        	assert(var(l)<S->nVars());
        	assert(super_interface(l));
        	return l;
        }
        inline Var fromSuper(Var p){
        	return  p-super_offset;
        }

        inline Lit fromSuper(Lit p){
        	return  mkLit(var(p)-super_offset,sign(p));
        }
        inline bool local_interface(Lit p)const{
            	return var(p)>=min_local && var(p) <= max_local;
            }
          inline bool local_interface(Var v)const{
          	return v>=min_local  && v <= max_local;
          }
        inline bool super_interface(Lit p)const{
          	return var(p)>=min_super && var(p) <= max_super;
          }
        inline bool super_interface(Var v)const{
        	return v>=min_super  && v <= max_super;
        }

        inline void dbg_check_propagation(Lit p){
  #ifdef DEBUG_SOLVER
      	  	 if (dbg_solver){

      	  		static vec<Lit> c;
  				c.clear();
  				for(int i = 0;i<trail.size();i++)
  				{
  					c.push(trail[i]);
  				}
  				c.push(~p);
  				bool res = dbg_solver->solve(c);

  				assert(!res);
      	  	 }
  #endif
        }

      inline void dbg_check(const vec<Lit> & clause){
#ifdef DEBUG_SOLVER
    	  	 if (dbg_solver){
    	  		 static bool first = true;
    	  		static vec<Lit> c;
				c.clear();
				for(int i = 0;i<clause.size();i++)
				{
					c.push(~ clause[i]);
				}

				bool res = dbg_solver->solve(c);

				assert(!res);
    	  	 }
#endif
      }

  
};


//=================================================================================================
// Implementation of inline methods:

inline CRef Solver::reason(Var x) const { return vardata[x].reason; }
inline int  Solver::level (Var x) const { return vardata[x].level; }

inline void Solver::insertVarOrder(Var x) {
    if (!order_heap.inHeap(x) && decision[x]) order_heap.insert(x); }

inline void Solver::varDecayActivity() { var_inc *= (1 / var_decay); }
inline void Solver::varBumpActivity(Var v) { varBumpActivity(v, var_inc); }
inline void Solver::varBumpActivity(Var v, double inc) {
    if ( (activity[v] += inc) > 1e100 ) {
        // Rescale:
        for (int i = 0; i < nVars(); i++)
            activity[i] *= 1e-100;
        var_inc *= 1e-100; }

    // Update order_heap with respect to new activity:
    if (order_heap.inHeap(v))
        order_heap.decrease(v); }

inline void Solver::claDecayActivity() { cla_inc *= (1 / clause_decay); }
inline void Solver::claBumpActivity (Clause& c) {
        if (c.learnt() && (c.activity() += cla_inc) > 1e20 ) {
            // Rescale:
            for (int i = 0; i < learnts.size(); i++)
                ca[learnts[i]].activity() *= 1e-20;
            cla_inc *= 1e-20; } }

inline void Solver::checkGarbage(void){ return checkGarbage(garbage_frac); }
inline void Solver::checkGarbage(double gf){
    if (ca.wasted() > ca.size() * gf)
        garbageCollect(); }

// NOTE: enqueue does not set the ok flag! (only public methods do)
inline bool     Solver::enqueue         (Lit p, CRef from)      { return value(p) != l_Undef ? value(p) != l_False : (uncheckedEnqueue(p, from), true); }
inline bool     Solver::addClause       (const vec<Lit>& ps)    { ps.copyTo(add_tmp); return addClause_(add_tmp); }
inline bool     Solver::addEmptyClause  ()                      { add_tmp.clear(); return addClause_(add_tmp); }
inline bool     Solver::addClause       (Lit p)                 { add_tmp.clear(); add_tmp.push(p); return addClause_(add_tmp); }
inline bool     Solver::addClause       (Lit p, Lit q)          { add_tmp.clear(); add_tmp.push(p); add_tmp.push(q); return addClause_(add_tmp); }
inline bool     Solver::addClause       (Lit p, Lit q, Lit r)   { add_tmp.clear(); add_tmp.push(p); add_tmp.push(q); add_tmp.push(r); return addClause_(add_tmp); }
inline bool     Solver::locked          (const Clause& c) const {
	CRef r = reason(var(c[0]));
	bool isClause = ca.isClause( r);
	if(!isClause)
		return false;
	return value(c[0]) == l_True && ca.isClause( reason(var(c[0]))) && ca.lea(reason(var(c[0]))) == &c; }
inline void     Solver::newDecisionLevel()                      { trail_lim.push(trail.size());}

inline int      Solver::decisionLevel ()      const   { return trail_lim.size(); }
inline uint32_t Solver::abstractLevel (Var x) const   { return 1 << (level(x) & 31); }
inline lbool    Solver::value         (Var x) const   { return assigns[x]; }
inline lbool    Solver::value         (Lit p) const   { return assigns[var(p)] ^ sign(p); }
inline lbool    Solver::modelValue    (Var x) const   { return model[x]; }
inline lbool    Solver::modelValue    (Lit p) const   { return model[var(p)] ^ sign(p); }
inline int      Solver::nAssigns      ()      const   { return trail.size(); }
inline int      Solver::nClauses      ()      const   { return clauses.size(); }
inline int      Solver::nLearnts      ()      const   { return learnts.size(); }
inline int      Solver::nVars         ()      const   { return vardata.size(); }
inline int 	    Solver::nUnassignedVars()	  const	  { return (int)dec_vars - trail.size();}
inline int      Solver::nFreeVars     ()      const   { return (int)dec_vars - ((trail_lim.size() == 0) ? trail.size() : trail_lim[0]); }
inline void     Solver::setPolarity   (Var v, bool b) { polarity[v] = b; }
inline void     Solver::setDecisionVar(Var v, bool b) 
{ 
    if      ( b && !decision[v]) dec_vars++;
    else if (!b &&  decision[v]) dec_vars--;

    decision[v] = b;
    insertVarOrder(v);
}
inline void     Solver::setConfBudget(int64_t x){ conflict_budget    = conflicts    + x; }
inline void     Solver::setPropBudget(int64_t x){ propagation_budget = propagations + x; }
inline void     Solver::interrupt(){ asynch_interrupt = true; }
inline void     Solver::clearInterrupt(){ asynch_interrupt = false; }
inline void     Solver::budgetOff(){ conflict_budget = propagation_budget = -1; }
inline bool     Solver::withinBudget() const {
    return !asynch_interrupt &&
           (conflict_budget    < 0 || conflicts < (uint64_t)conflict_budget) &&
           (propagation_budget < 0 || propagations < (uint64_t)propagation_budget); }

// FIXME: after the introduction of asynchronous interrruptions the solve-versions that return a
// pure bool do not give a safe interface. Either interrupts must be possible to turn off here, or
// all calls to solve must return an 'lbool'. I'm not yet sure which I prefer.
inline bool     Solver::solve         ()                    { budgetOff(); assumptions.clear(); return solve_() == l_True; }
inline bool     Solver::solve         (Lit p)               { budgetOff(); assumptions.clear(); assumptions.push(p); return solve_() == l_True; }
inline bool     Solver::solve         (Lit p, Lit q)        { budgetOff(); assumptions.clear(); assumptions.push(p); assumptions.push(q); return solve_() == l_True; }
inline bool     Solver::solve         (Lit p, Lit q, Lit r) { budgetOff(); assumptions.clear(); assumptions.push(p); assumptions.push(q); assumptions.push(r); return solve_() == l_True; }
inline bool     Solver::solve         (const vec<Lit>& assumps){ budgetOff(); assumps.copyTo(assumptions); return solve_() == l_True; }
inline lbool    Solver::solveLimited  (const vec<Lit>& assumps){ assumps.copyTo(assumptions); return solve_(); }
inline bool     Solver::okay          ()      const   { return ok; }

inline void     Solver::toDimacs     (const char* file){ vec<Lit> as; toDimacs(file, as); }
inline void     Solver::toDimacs     (const char* file, Lit p){ vec<Lit> as; as.push(p); toDimacs(file, as); }
inline void     Solver::toDimacs     (const char* file, Lit p, Lit q){ vec<Lit> as; as.push(p); as.push(q); toDimacs(file, as); }
inline void     Solver::toDimacs     (const char* file, Lit p, Lit q, Lit r){ vec<Lit> as; as.push(p); as.push(q); as.push(r); toDimacs(file, as); }


//=================================================================================================
// Debug etc:


//=================================================================================================
}

#endif
