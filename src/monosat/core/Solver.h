/**************************************************************************************************
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

#ifndef Monosat_Solver_h
#define Monosat_Solver_h

#include "monosat/mtl/Vec.h"
#include "monosat/mtl/Heap.h"
#include "monosat/mtl/Alg.h"
#include "monosat/mtl/Rnd.h"
#include "monosat/utils/Options.h"
#include "monosat/core/SolverTypes.h"
#include "monosat/core/Theory.h"
#include "monosat/core/TheorySolver.h"
#include "monosat/core/Config.h"
#include <cinttypes>
#include <map>
#include <string>
#include <sstream>
#include <iostream>
#include <vector>

namespace Monosat {
template<typename Weight>
class GraphTheorySolver;

class FSMTheorySolver;

class DimacsMap;

//=================================================================================================
// Solver -- the main class:
// The MiniSAT Boolean SAT solver, extended to provided basic SMT support.
class Solver : public Theory, public TheorySolver {
public:
    void* _external_data = nullptr;//convenience pointer for external API.
    static bool shown_warning;

    //fix this...
    friend class Theory;

    template<typename Weight> friend
    class GraphTheorySolver;

    friend class FSMTheorySolver;

    friend class LSystemSolver;

    // Constructor/Destructor:
    //
    Solver();

    ~Solver() override;

    // Problem specification:
    //
    Var newVar(bool polarity = true,
               bool dvar = true) override; // Add a new variable with parameters specifying variable mode.
    virtual void releaseVar(
            Lit l);                                  // Make literal true and promise to never refer to variable again.

    bool addClause(const vec<Lit>& ps) override;                     // Add a clause to the solver.
    virtual bool addEmptyClause();                             // Add the empty clause, making the solver contradictory.
    bool addClause(Lit p) override;                                  // Add a unit clause to the solver.
    bool addClause(Lit p, Lit q) override;                           // Add a binary clause to the solver.
    bool addClause(Lit p, Lit q, Lit r) override;                    // Add a ternary clause to the solver.
    bool addClause_(vec<Lit>& ps,
                    bool is_derived_clause = false);           // Add a clause to the solver without making superflous internal copy. Will

    void disableElimination(Var v) override{
        //do nothing
    }

    void setDecisionPriority(Var v, unsigned int p){
        priority[v] = p;
        if(decision[v]){
            if(order_heap.inHeap(v))
                order_heap.decrease(v);
            else
                order_heap.insert(v);
        }
    }

    int getDecisionPriority(Var v) const{
        return priority[v];
    }

    void setTheorySatisfied(Theory* theory) override{
        int theoryID = theory->getTheoryIndex();
        if(!theorySatisfied(theory)){
            //printf("Theory %d sat at %d\n",theoryID, decisionLevel());
            if(trail.size() > 0){
                satisfied_theory_trail_pos[theoryID] = trail.size() - 1;
            }else{
                //this really shouldn't ever happen...
                satisfied_theory_trail_pos[theoryID] = 0;
            }
            post_satisfied_theory_trail_pos[theoryID] = satisfied_theory_trail_pos[theoryID];
            //theory_sat_queue.push(TheorySatisfaction(theoryID,trail.size()));
        }
    }

    bool theorySatisfied(Theory* theory) override{
        int theoryID = theory->getTheoryIndex();
        return satisfied_theory_trail_pos[theoryID] >= 0;
    }

    virtual bool heuristicSatisfied(Heuristic* h){
        int theoryID = h->getTheoryIndex();
        if(theoryID >= 0){
            return satisfied_theory_trail_pos[theoryID] >= 0;
        }else{
            return false;
        }
    }

    void clearSatisfied() override{
        for(Theory* t:theories){
            if(t){
                int theoryID = t->getTheoryIndex();
                satisfied_theory_trail_pos[theoryID] = -1;
                post_satisfied_theory_trail_pos[theoryID] = -1;
                t->clearSatisfied();
            }
        }
    }

    //Theory interface
    void addTheory(Theory* t) override{
        if(t->getName().size() > 0){
            if(theorymap.count(t->getName()) > 0){
                throw std::invalid_argument("All theory names must be unique: " + t->getName());
            }
            theorymap.insert({t->getName(), t});
        }
        satisfied_theory_trail_pos.push(-1);
        post_satisfied_theory_trail_pos.push(-1);
        theories.push(t);
        theory_reprop_trail_pos.push(-1);
        theory_init_prop_trail_pos.push(-1);
        t->setTheoryIndex(theories.size() - 1);
        if(t->supportsDecisions()){
            addHeuristic(t);
        }
        theory_queue.capacity(theories.size());
        in_theory_queue.push(false);

        cancelUntil(0);
        resetInitialPropagation();
    }

    void addHeuristic(Heuristic* t) override{
        if(t->getHeuristicIndex() >= 0){
            assert(t->getHeuristicIndex() < all_decision_heuristics.size());
            assert(all_decision_heuristics[t->getHeuristicIndex()] == t);
            return;
        }
        int heuristic_id = all_decision_heuristics.size();
        assert(heuristic_id > 0);
        t->setHeuristicIndex(heuristic_id);
        all_decision_heuristics.push(t);
        t->setHeuristicOrder(decision_heuristics.size());
        decision_heuristics.push(t);
        theory_order_heap.insert(t);
        theory_conflict_counters.growTo(all_decision_heuristics.size(), 0);
        first_heuristic_decision_level.growTo(all_decision_heuristics.size(), -1);
        t->setActivity(opt_randomize_theory_order_restart_freq > 0 ? drand(random_seed) * 0.00001 : 0);
        t->setPriority(0);
    }

    void activateHeuristic(Heuristic* h) override{
        assert(all_decision_heuristics.contains(h));
        if(h->getHeuristicIndex() > all_decision_heuristics.size() || h->getHeuristicIndex() < 0 ||
           all_decision_heuristics[h->getHeuristicIndex()] != h){
            throw std::runtime_error("Unknown decision heuristic");
        }
        theory_order_heap.update(h);
    }

    void theoryPropagated(Theory* t) override{
        int theoryID = t->getTheoryIndex();
        if(decisionLevel() > 0 && theory_init_prop_trail_pos[theoryID] >= 0 && theory_reprop_trail_pos[theoryID] < 0){
            theory_reprop_trail_pos[theoryID] = qhead;
        }
    }

    //Call to force at least one round of propagation to each theory solver at the next solve() call
    void resetInitialPropagation(){
        if(!initialPropagate){
            initialPropagate = true;                 //to force propagation to occur at least once to the theory solvers
            if(S){
                S->resetInitialPropagation();
            }
        }
    }

    void detatchTheory(Theory* t){
        assert(decisionLevel() == 0);
        assert(theory_queue.size() == 0);
        int i, j = 0;
        for(i = 0; i < theories.size(); i++){
            if(theories[i] == t){
                int k, l = 0;
                for(k = 0; k < markers.size(); k++){
                    CRef cr = markers[k];
                    int index = CRef_Undef - cr - 1;
                    int theory = marker_theory[index];
                    if(theory == i){
                        marker_theory[index] = -1;
                    }else{
                        markers[l++] = cr;
                        if(marker_theory[index] > i){
                            marker_theory[index]--;
                        }
                    }
                }
                markers.shrink(k - l);
            }else{
                theories[j++] = theories[i];
                in_theory_queue[j++] = in_theory_queue[i];
            }
        }
        theories.shrink(i - j);
        in_theory_queue.shrink(i - j);

        cancelUntil(0);
    }

    int getTheoryIndex() const override{
        return theory_index;
    }

    void setTheoryIndex(int id) override{
        theory_index = id;
    }

    //Generate a new, unique `temporary value' for explaining conflicts
    CRef newReasonMarker(Heuristic* forTheory, bool is_decision = false) override{
        markers.push(ca.makeMarkerReference());

        int marker_num = CRef_Undef - markers.last() - 1;
        marker_theory.growTo(marker_num + 1, -1);

        if(is_decision){
            assert(forTheory->getHeuristicIndex() > 0);//heuristic indices must be strictly greater than 0
            marker_theory[marker_num] = -forTheory->getHeuristicIndex();
        }else{
            marker_theory[marker_num] = forTheory->getTheoryIndex();
        }

        //this could be done more efficiently
        /*for (int i = 0; i < theories.size(); i++) {
            if (theories[i] == forTheory) {
                marker_theory[marker_num] = i;
                break;
            }
        }*/
        //assert(marker_theory[marker_num] >= 0);

        return markers.last();
    }

    void printStats(int detail_level = 0) override{

        double cpu_time = cpuTime();
        //double mem_used = memUsedPeak();


        printf("restarts              : %" PRIu64 "\n", starts);
        if(opt_theory_order_conflict_restart && stats_theory_conflict_counter_restarts > 0){
            printf("conflict counter restarts:     %" PRIu64 "\n", stats_theory_conflict_counter_restarts);
        }
        if(opt_theory_order_swapping_reset_counts_new_conflict){
            printf("swapping counter clears/total, max swapping counter value:     %" PRId64 "/%" PRId64 ", %" PRId64 "\n",
                   stats_max_swap_count, stats_swapping_conflict_count, stats_swapping_resets);
        }
        printf("conflicts             : %-12" PRIu64 "   (%.0f /sec, %d learnts (%" PRId64 " theory learnts), %" PRId64 " removed)\n",
               conflicts,
               conflicts / cpu_time, learnts.size(), stats_theory_conflicts, stats_removed_clauses);
        printf("decisions             : %-12" PRIu64 "   (%4.2f %% random) (%.0f /sec)\n", decisions,
               (float) rnd_decisions * 100 / (float) decisions, decisions / cpu_time);
        if(opt_decide_theories){
            printf("Theory decisions: %" PRId64 "\n", stats_theory_decisions);
            printf("Theory decision rounds: %" PRId64 "/%" PRId64 "\n", n_theory_decision_rounds, starts);
        }
        if(opt_vsids_both){
            printf("Sovler pre-empted decisions: %" PRId64 "\n", stats_solver_preempted_decisions);
        }
        printf("propagations          : %-12" PRIu64 "   (%.0f /sec)\n", propagations, propagations / cpu_time);
        printf("conflict literals     : %-12" PRIu64 "   (%4.2f %% deleted)\n", tot_literals,
               (max_literals - tot_literals) * 100 / (double) max_literals);
        if(stats_skipped_theory_prop_rounds > 0){
            printf("theory propagations skipped: %" PRId64 "\n", stats_skipped_theory_prop_rounds);
        }
        if(opt_detect_pure_theory_lits){
            printf("pure literals     : %" PRId64 " (%" PRId64 " theory lits) (%" PRId64 " rounds, %f time)\n",
                   stats_pure_lits,
                   stats_pure_theory_lits, pure_literal_detections, stats_pure_lit_time);
        }
        if(stats_theory_prop_time > 0){
            printf("Time spent in theory propagation: %f\n", stats_theory_prop_time);
        }
        if(stats_theory_conflict_time > 0){
            printf("Time spent in theory conflicts: %f\n", stats_theory_conflict_time);
        }
        if(opt_check_solution){
            printf("Solution double-checking time (disable with -no-check-solution): %f s\n",
                   stats_solution_checking_time);
        }
        for(int i = 0; i < theories.size(); i++){
            theories[i]->printStats(detail_level);
        }
    }

    void writeTheoryWitness(std::ostream& write_to) override{
        if(!ok){
            write_to << "s UNSATISFIABLE\n";
        }else{
            write_to << "s SATISFIABLE\n";
            write_to << "v ";
            if(model.size() >= nVars()){
                for(int v = 0; v < nVars(); v++){
                    if(model[v] == l_True){
                        write_to << " " << (v + 1);
                    }else{
                        write_to << " " << -(v + 1);
                    }
                }
            }
            write_to << " 0\n";

            for(Theory* t : theories){
                t->writeTheoryWitness(write_to);
            }
        }
    }

    bool isTheoryCause(CRef cr) const{
        return cr != CRef_Undef && !ca.isClause(cr);
    }

    int nAssumptions() override{
        return assumptions.size();
    }

    inline Theory* getTheory(CRef cr) const{
        assert(isTheoryCause(cr));
        assert(!isDecisionReason(cr));
        // UINT32_MAX-cr - 1;
        int marker = CRef_Undef - cr - 1;
        assert(marker_theory[marker] >= 0);
        assert(marker_theory[marker] < theories.size());
        return theories[marker_theory[marker]];
    }

    Heuristic* getHeuristic(CRef cr) const{
        assert(isDecisionReason(cr));
        // UINT32_MAX-cr - 1;
        int marker = CRef_Undef - cr - 1;
        int hID = -(marker_theory[marker]);
        assert(hID >= 0);
        assert(hID <= all_decision_heuristics.size());
        return all_decision_heuristics[hID];
    }


    inline bool hasTheory(Var v){
        return v >= 0 && v < theory_vars.size() && theory_vars[v].theories.size() > 0;
    }

    inline bool hasTheory(Lit l){
        return hasTheory(var(l));
    }

    inline int getNTheories(Var v){
        if(!hasTheory(v)){
            return 0;
        }else{
            return theory_vars[v].theories.size();
        }
    }

    inline int getTheoryID(Var v, int theoryN){
        return theory_vars[v].theories[theoryN].theory;
    }

    inline int getTheoryID(Lit l, int theoryN){
        return getTheoryID(var(l), theoryN);
    }

    inline Var getTheoryVar(Var v, int theoryN){
        assert(hasTheory(v));
        //is this check too expensive? Does getTheoryVar appear in any tight loops?
        if(!hasTheory(v)){
            throw std::runtime_error("Literal " + std::to_string(toInt(toLit(v))) + " is not a theory atom");
        }
        return (Var) theory_vars[v].theories[theoryN].theory_var;
    }

    vec<Theory*>& getTheories() override{
        return theories;
    }

    //Translate a literal into its corresponding theory literal (if it has a theory literal)
    Lit getTheoryLit(Lit l, int theoryN){
        assert(hasTheory(l));
        return mkLit(getTheoryVar(var(l), theoryN), sign(l));
    }

    bool theoryHasVar(Var solverVar, Theory* t){
        int theoryID = t->getTheoryIndex();
        return hasTheory(solverVar) && theory_vars[solverVar].theoryMap.has(theoryID) &&
               theory_vars[solverVar].theoryMap[theoryID] != var_Undef;
    }

    Var getTheoryVar(Var v, Theory* t){
        return theory_vars[v].theoryMap[t->getTheoryIndex()];
    }

    Lit getTheoryLit(Lit l, Theory* t){
        return mkLit(getTheoryVar(var(l), t), sign(l));
        /*assert(getNTheories(var(l))<10);
        for(int n = 0;n<getNTheories(var(l));n++){
            if(getTheoryID(l,n)==t->getTheoryIndex()){
                return getTheoryLit(l,n);
            }
        }
        throw std::runtime_error("No such theory lit");*/
    }

    Var newTheoryVar(Var solverVar, int theoryID, Var theoryVar) override{
        while(nVars() <= solverVar)
            newVar();
        Var v = solverVar;
        setTheoryVar(solverVar, theoryID, theoryVar);
        return v;
    }

    //Connect a variable in the SAT solver to a variable in a theory.
    virtual void setTheoryVar(Var solverVar, int theory, Var theoryVar){
        if(solverVar == var(const_true)){
            //handle the True literal specially, so that all the theories can share it...
            return;
        }
        if(theoryHasVar(solverVar, theories[theory]) && getTheoryVar(solverVar, theories[theory]) != theoryVar){
            throw std::runtime_error("Variable used for multiple atoms in one theory.");
        }
        //assert(!hasTheory(solverVar));
        int theoryN = theory_vars[solverVar].theories.size();
        theory_vars[solverVar].theories.push({theory, theoryVar});
        theory_vars[solverVar].theoryMap.insert(theory, theoryVar, var_Undef);
        all_theory_vars.push(solverVar);
        assert(hasTheory(solverVar));
        assert(getTheoryID(solverVar, theoryN) == theory);
        assert(getTheoryVar(solverVar, theoryN) == theoryVar);

        //if (value(solverVar) != l_Undef)
        initialPropagate = true;//if pure literal detection is used, then we NEED to run simplify before
        //the next call to propagate, and initialPropagate is being used to force that.

    }

    const char* getTheoryType() override{
        return "Solver";
    }

    Theory* getTheory(const std::string& name){
        if(theorymap.count(name) > 0){
            return theorymap[name];
        }else{
            return nullptr;
        }
    }

    void preprocess() override{
        for(int i = 0; i < theories.size(); i++){
            theories[i]->preprocess();
        }
    }

    //Lazily construct a reason for a literal propagated from a theory
    CRef constructReason(Lit p){
        assert(value(p) == l_True);
        CRef cr = reason(var(p));
        assert(isTheoryCause(cr));
        assert(!ca.isClause(cr));
        assert(cr != CRef_Undef);
        assert(!isDecisionReason(cr));
        int trail_pos = trail.size();
        Theory* t = getTheory(cr);
        assert(t);
        assert(hasTheory(p));
        theory_reason.clear();

        double start_t = rtime(1);

        t->buildReason(getTheoryLit(p, t), theory_reason, cr);

        stats_theory_conflict_time += (rtime(1) - start_t);
        assert(theory_reason[0] == p);
        assert(value(p) == l_True);

#ifdef DEBUG_CORE
        for(Lit l:theory_reason)
            assert(value(l) != l_Undef);
#endif
#ifdef DEBUG_SOLVER
        //assert all the other reasons in this cause are earlier on the trail than p...
        static vec<bool> marks;
        marks.clear();
        marks.growTo(nVars());
        for(int i = 0;i<trail.size() && var(trail[i])!=var(p);i++) {
            marks[var(trail[i])]=true;
        }
        for(int i = 1;i<theory_reason.size();i++) {
            assert(marks[var(theory_reason[i])]);
        }
#endif

        int lev = decisionLevel();
        CRef reason = attachReasonClause(p, theory_reason);
        vardata[var(p)] = mkVarData(reason, level(var(p)));
        assert(decisionLevel() == lev);//ensure no backtracking happened while adding this clause!
        assert(ok);

        //collect any newly enqueued vars -  we need to analyze them

        for(int i = trail_pos; i < trail.size(); i++){
            to_analyze.push(trail[i]);
        }


        return reason;
    }

    //Attach this solver to the next, connecting literals in the super solver to literals in this solver
    //Super_vars and local_vars _must_ be consecutive sets of variables in the super solver.
    void attachTo(Solver* super, vec<Var>& super_vars, vec<Var>& local_vars){
        assert(decisionLevel() == 0);
        super->addTheory(this);
        S = super;
        assert(super_vars.size() == local_vars.size());
        super_qhead = 0;
        local_qhead = 0;
        super_offset = super_vars[0] - local_vars[0];
        min_super = super_vars[0];
        max_super = super_vars.last();
        min_local = local_vars[0];
        max_local = local_vars.last();
        cause_marker = super->newReasonMarker(this);
    }
    // change the passed vector 'ps'.

    // Solving:
    //
    void detectPureTheoryLiterals(); //if opt_detect_pure_literals, finds pure theory literals.
    bool simplify();                        // Removes already satisfied clauses.
    virtual bool solve(const vec<Lit>& assumps); // Search for a model that respects a given set of assumptions.
    virtual lbool solveLimited(
            const vec<Lit>& assumps); // Search for a model that respects a given set of assumptions (With resource constraints).
    virtual bool solve();                        // Search without assumptions.
    virtual bool solve(Lit p);                   // Search for a model that respects a single assumption.
    virtual bool solve(Lit p, Lit q);            // Search for a model that respects two assumptions.
    virtual bool solve(Lit p, Lit q, Lit r);     // Search for a model that respects three assumptions.
    virtual bool propagateAssignment(
            const vec<Lit>& assumps); //apply unit propagation to the supplied assumptions, and quit without solving
    virtual lbool
    solveUntilRestart(const vec<Lit>& assumps);//attempt to solve the instance, but quit as soon as the solver restarts
    bool okay() const;                  // FALSE means solver is in a conflicting state
    void contradiction(){ //put the solver into a contradictory state
        ok = false;
    }

    Lit True() override{

        if(const_true == lit_Undef){
            cancelUntil(0);
            //try using the first assigned const literal
            /*if (trail.size()>0){
                Lit l = trail[0];
                if(level(var(l))==0){
                    const_true=l;
                }
            }else{*/
            const_true = mkLit(newVar(false, false));
            if(ok){
                addClause(const_true);
                assert(isConstantTrue(const_true));
                addLiteralName(const_true, "True");
            }
        }

        return const_true;
    }

    Lit unmap(Lit l) override;

    Var unmap(Var v);

    Lit mapLit(Lit l);

    Var mapVar(Var v);

    /*
     * The number of externally visible variables (which may be different from the number of internal variables).
     */
    int nMappedVars();

    void setVarMap(DimacsMap* map){
        varRemap = map;
    }

    bool hasLiteral(const std::string& name){
        return getLiteral(name) != lit_Undef;
    }

    int literalNameCount(Lit l){
        if((var(l) < 0) || !litnames.count(toInt(l))){
            return 0;
        }else{
            return litnames[toInt(l)].size();
        }
    }

    bool hasName(Lit l, std::string name){
        return var(l) >= 0 && getLiteral(name) == l;
    }

    //Associate variable v with a name.
    //If name is empty, then remove any existing name associated with v
    //If the variable already has this name, do nothing
    //Assigned the negated polarity the same name prefixed with '~'
    void addLiteralName(Lit l, const std::string& name);

    //Returns the nth name associated with this variable, or an empty string of there is no name.
    const std::string& getLiteralName(Lit l, int nameIndex = 0){
        if(l == lit_Undef || !litnames.count(toInt(l))){
            return empty_name;
        }
        std::vector<std::string>& names = litnames[toInt(l)];
        if(nameIndex < 0 || names.size() < nameIndex){
            return empty_name;
        }else{
            return names[nameIndex];
        }
    }

    //Get the variable associated with this name (if any).
    //Returns Var_Undef if there is no variable assocaited with this name.
    Lit getLiteral(const std::string& name){
        if(!namemap.count(name)){
            return lit_Undef;
        }
        return namemap[name];
    }

    const vec<Lit>& namedLiterals(){
        return named_literals;
    }


    // Variable mode:
    //
    void setPolarity(Var v,
                     bool b); // Declare which polarity the decision heuristic should use for a variable. Requires mode 'polarity_user'.
    bool getPolarity(Var v){
        return polarity[v];
    }

    void setDecisionVar(Var v,
                        bool b) override; // Declare if a variable should be eligible for selection in the decision heuristic.
    bool isDecisionVar(Var v){
        return decision[v];
    }

    // Read state:
    //
    lbool value(Var x) const override;       // The current value of a variable.
    lbool value(Lit p) const override;       // The current value of a literal.
    lbool modelValue(
            Var x) const; // The value of a variable in the last model. The last call to solve must have been satisfiable.
    lbool modelValue(
            Lit p) const; // The value of a literal in the last model. The last call to solve must have been satisfiable.
    bool hasModel() const; //True if the solver has a model
    int nAssigns() const;       // The current number of assigned literals.
    int nClauses() const;       // The current number of original clauses.
    int nLearnts() const;       // The current number of learnt clauses.
    int nVars() const override;       // The current number of variables.
    int
    nUnassignedVars() const;// The number of variables left to assign. Does not include variables that are not decision vars!
    int nFreeVars() const;        // The number of non-unit variables
    bool hasNextDecision() const;

    // Resource contraints:
    //
    void setConfBudget(int64_t x);

    void setPropBudget(int64_t x);

    int64_t getConflictBudget();

    int64_t getPropagationBudget();

    void budgetOff();

    void interrupt();          // Trigger a (potentially asynchronous) interruption of the solver.
    void clearInterrupt();     // Clear interrupt indicator flag.

    // Memory managment:
    //
    virtual void garbageCollect();

    void checkGarbage(double gf);

    void checkGarbage();


    // Extra results: (read-only member variable)
    //
    vec<lbool> model;             // If problem is satisfiable, this vector contains the model (if any).
    LSet conflict;          // If problem is unsatisfiable (possibly under assumptions),
    // this vector represent the final conflict clause expressed in the assumptions.

    vec<vec<Lit>> interpolant; //This vector represents an interpolant between this module and its super solver ('S'), if it is attached to such a solver and the instance is UNSAT.
    // Variables in each clause in the interpolant vector are in the super solver's variable space, not the subsolver's.

    Lit const_true = lit_Undef;
    Lit theoryDecision = lit_Undef;
    vec<Lit> theory_reason;
    vec<Lit> theory_conflict;
    vec<Theory*> theories;
    vec<int> satisfied_theory_trail_pos;
    vec<int> post_satisfied_theory_trail_pos;
    vec<Heuristic*> decision_heuristics;

    int decision_heuristic_qhead = 0;
    vec<int> decision_heuristic_trail_lim;

    vec<Heuristic*> all_decision_heuristics;
    vec<int> first_heuristic_decision_level;//the first decision made by each decision heuristic (other than vsids)
    vec<Heuristic*> decision_heuristic_trail;//trail of decisions made by decision heuristics (other than vsids)
    vec<int> theory_conflict_counters;
    int theory_decision_round_robin = 0;
    Heuristic* decisionTheory = nullptr;//for opt_vsids_solver_as_theory
    vec<Var> all_theory_vars;

    struct LitCount {
        char occurs :1;
        char seen :1;

        LitCount() :
                occurs(0), seen(0){
        }
    };

    vec<LitCount> lit_counts;
    vec<CRef> markers;                        //a set of special clauses that can be recognized as pointers to theories
    vec<int> marker_theory;
    int theory_index = 0;
    Solver* S = nullptr;                                //super solver
    Theory* bvtheory = nullptr;
    PB::PBConstraintSolver* pbsolver = nullptr;
    bool initialPropagate = true;                    //to force propagation to occur at least once to the theory solvers
    int super_qhead = 0;
    int local_qhead = 0;
    CRef cause_marker = CRef_Undef;
    int track_min_level = 0;
    int initial_level = 0;
/*	struct TheorySatisfaction{
		int theoryID=-1;
		int trail_size=-1;
		TheorySatisfaction(int theoryID, int trail_size):theoryID(theoryID),trail_size(trail_size){

		}
	};
	vec<TheorySatisfaction> theory_sat_queue;*/
    vec<int> theory_queue;
    vec<bool> in_theory_queue;
    IntSet<int> unskippable_theory_q;
    vec<int> theory_reprop_trail_pos;
    vec<int> theory_init_prop_trail_pos;
    bool disable_theories = false;
    int min_decision_var = 0;
    int max_decision_var = -1;
    int min_priority_var = 0;
    int max_priority_var = -1;
    CRef tmp_clause = CRef_Undef;
    vec<Lit> tmp_conflict;
    int tmp_clause_sz = 0;
    Var max_super = var_Undef;
    Var min_super = var_Undef;
    Var min_local = var_Undef;
    Var max_local = var_Undef;
    int super_offset = -1;
    DimacsMap* varRemap = nullptr;

    std::map<int, std::vector<std::string>> litnames;  //Optional names associated with each literal
    std::map<std::string, Lit> namemap; //literal name lookup map
    vec<Lit> named_literals;//all literals that have names, in the order they were named

    static const std::string empty_name;


    std::map<std::string, Theory*> theorymap;
    // Mode of operation:
    //
    bool printed_header = false;
    int verbosity;
    double var_decay;
    double clause_decay;
    double random_var_freq;
    double random_seed;
    bool luby_restart;
    int ccmin_mode;         // Controls conflict clause minimization (0=none, 1=basic, 2=deep).
    int phase_saving;       // Controls the level of phase saving (0=none, 1=limited, 2=full).
    bool rnd_pol;            // Use random polarities for branching heuristics.
    bool rnd_init_act;       // Initialize variable activities with a small random value.
    double garbage_frac;       // The fraction of wasted memory allowed before a garbage collection is triggered.

    int restart_first; // The initial restart limit.                                                                (default 100)
    double restart_inc; // The factor with which the restart limit is multiplied in each restart.                    (default 1.5)
    double learntsize_factor; // The intitial limit for learnt clauses is a factor of the original clauses.                (default 1 / 3)
    double learntsize_inc; // The limit for learnt clauses is multiplied with this factor each restart.                 (default 1.1)

    int learntsize_adjust_start_confl;
    double learntsize_adjust_inc;

    // Statistics: (read-only member variable)
    //
    double stats_solution_checking_time = 0;
    uint64_t solves = 0;
    uint64_t starts = 0;
    uint64_t decisions = 0;
    uint64_t rnd_decisions = 0;
    uint64_t propagations = 0;
    uint64_t conflicts = 0;
    int64_t stats_pure_lits = 0;
    uint64_t stats_pure_theory_lits = 0;
    uint64_t pure_literal_detections = 0;
    uint64_t stats_removed_clauses = 0;
    uint64_t dec_vars, clauses_literals, learnts_literals, max_literals, tot_literals;
    uint64_t stats_skipped_theory_prop_rounds = 0;
    uint64_t stats_theory_conflict_counter_restarts = 0;
    uint64_t stats_max_swap_count = 0;
    uint64_t stats_swapping_resets = 0;
    uint64_t stats_swapping_conflict_count = 0;
    uint64_t stats_theory_conflicts = 0;
    double stats_theory_prop_time = 0;
    double stats_theory_conflict_time = 0;

    uint64_t stats_solver_preempted_decisions = 0;
    uint64_t stats_theory_decisions = 0;
    double stats_pure_lit_time = 0;
    uint64_t n_theory_conflicts = 0;
    int consecutive_theory_conflicts = 0;
    uint64_t next_theory_decision = 0;
    uint64_t n_theory_decision_rounds = 0;

    //Var last_dec=var_Undef;
protected:

    // Helper structures:
    //
    struct VarData {
        CRef reason;
        int level;
    };

    static inline VarData mkVarData(CRef cr, int l){
        VarData d = {cr, l};
        return d;
    }

    struct Watcher {
        CRef cref;
        Lit blocker;

        Watcher(CRef cr, Lit p) :
                cref(cr), blocker(p){
        }

        bool operator==(const Watcher& w) const{
            return cref == w.cref;
        }

        bool operator!=(const Watcher& w) const{
            return cref != w.cref;
        }
    };

    struct WatcherDeleted {
        const ClauseAllocator& ca;

        WatcherDeleted(const ClauseAllocator& _ca) :
                ca(_ca){
        }

        bool operator()(const Watcher& w) const{
            return ca[w.cref].mark() == 1;
        }
    };

    struct LazyLevelLt {

        Solver* outer;

        bool operator()(int x, int y) const{
            Var vx = toInt(toLit(x));
            Var vy = toInt(toLit(y));
            return outer->level(vx) > outer->level(vy);
        }

        LazyLevelLt(Solver* outer) :
                outer(outer){
        }
    };

    struct VarOrderLt {
        const vec<double>& activity;
        const vec<int>& priority;

        bool operator()(Var x, Var y) const{
            if(priority[x] == priority[y])
                return activity[x] > activity[y];
            else{
                return priority[x] > priority[y];
            }
        }

        VarOrderLt(const vec<double>& act, const vec<int>& pri) :
                activity(act), priority(pri){
        }
    };

    struct HeuristicOrderLt {

        bool operator()(Heuristic* x, Heuristic* y) const{
            if(x->getPriority() == y->getPriority()){
                if(x->getHeuristicOrder() == y->getHeuristicOrder()){
                    return x->getActivity() > y->getActivity();
                }else{
                    return x->getHeuristicOrder() < y->getHeuristicOrder();//lower ordered heuristic goes earlier
                }
            }else{
                return x->getPriority() > y->getPriority();
            }
        }

        HeuristicOrderLt(){

        }
    };

    struct HeuristicActivityOrderLt {

        bool operator()(Heuristic* x, Heuristic* y) const{
            if(x->getPriority() == y->getPriority()){
                //ignore heuristic order, even if it is set
                return x->getActivity() > y->getActivity();

            }else{
                return x->getPriority() > y->getPriority();
            }
        }

        HeuristicActivityOrderLt(){

        }
    };

    struct TheoryData {
        int theory;
        Var theory_var;
    };
    struct TheoryMap {
        IntMap<int, Var> theoryMap;
        vec<TheoryData> theories;
    };


    /*struct TheoryData {
        union {
            struct {
                unsigned int theory ;
                unsigned int theory_var;//these were previously packed together into a 32bit int, but that has proved insufficient in practice.
            };
            unsigned int:1 isTheoryVar; //true if non-zero - this property is ensured by adding 1 to theory_var
        };
        TheoryData() :
                isTheoryVar(0) {
        }
        TheoryData(unsigned int theory, unsigned int theory_lit) :
                theory(theory), theory_var(theory_lit) {

        }
    };*/

    // Solver state:
    //
    bool ok;            // If FALSE, the constraints are already unsatisfiable. No part of the solver state may be used!
    vec<CRef> clauses;          // List of problem clauses.
    vec<CRef> learnts;          // List of learnt clauses.
    double cla_inc;          // Amount to bump next clause with.
    vec<double> activity;         // A heuristic measurement of the activity of a variable.
    double var_inc;          // Amount to bump next variable with.
    OccLists<Lit, vec<Watcher>, WatcherDeleted> watches; // 'watches[lit]' is a list of constraints watching 'lit' (will go there if literal becomes true).
    Heuristic* conflicting_heuristic = nullptr;

    vec<lbool> assigns;          // The current assignments.
    vec<char> polarity;         // The preferred polarity of each variable.
    vec<char> decision;         // Declares if a variable is eligible for selection in the decision heuristic.
    vec<int> priority;          // Static, lexicographic heuristic. Larger values are higher priority (decided first).

    vec<TheoryMap> theory_vars;
    vec<Lit> to_analyze;
    vec<Lit> to_reenqueue;
    vec<Lit> trail;            // Assignment stack; stores all assigments made in the order they were made.
    vec<int> trail_lim;        // Separator indices for different decision levels in 'trail'.
    vec<int> partially_propagated_levels;//a list of the levels in the sovler for which theory_q was not empty when newDecisionLevel() was called. Used to force theories to be repropagated if theory propagation was skipped, in certain configurations.
    vec<VarData> vardata;          // Stores reason and level for each variable.
    int qhead;            // Head of queue (as index into the trail -- no more explicit propagation queue in MiniSat).
    int simpDB_assigns;   // Number of top-level assignments since last execution of 'simplify()'.
    int64_t simpDB_props;   // Remaining number of propagations that must be made before next execution of 'simplify()'.
    vec<Lit> assumptions;      // Current set of assumptions provided to solve by the user.
    bool only_propagate_assumptions = false; //true if the solver should propagate assumptions and then quit without solving
    bool quit_at_restart = false;//true if the solver should give up as soon as it restarts
    int override_restart_count = -1;
    Heap<Var, VarOrderLt> order_heap;       // A priority queue of variables ordered with respect to the variable activity.
    double theory_inc;
    double theory_decay;

    struct HeuristicToInt {
        int operator()(Heuristic* h) const{return h->getHeuristicIndex();}
    };

    Heap<Heuristic*, HeuristicOrderLt, HeuristicToInt> theory_order_heap;
    vec<std::pair<Heuristic*, int>> theory_decision_trail;
    //Heap<LazyLevelLt> lazy_heap;       // A priority queue of variables to be propagated at earlier levels, lazily.
    double progress_estimate;       // Set by 'search()'.
    bool remove_satisfied; // Indicates whether possibly inefficient linear scan for satisfied clauses should be performed in 'simplify'.

    ClauseAllocator ca;

    vec<Var> released_vars;
    vec<Var> free_vars;

    // Temporaries (to reduce allocation overhead). Each variable is prefixed by the method in which it is
    // used, exept 'seen' wich is used in several places.
    //
    vec<char> seen;
    vec<Lit> analyze_stack;
    vec<Lit> analyze_toclear;
    vec<Lit> add_tmp;

    vec<vec<Lit>> clauses_to_add;

    double max_learnts = 1;
    double learntsize_adjust_confl = 0;
    int learntsize_adjust_cnt = 0;

    // Resource contraints:
    //
    int64_t conflict_budget;    // -1 means no budget.
    int64_t propagation_budget; // -1 means no budget.
    bool asynch_interrupt = false;

    // Main internal methods:
    //
    void insertVarOrder(Var x);                               // Insert a variable in the decision order priority queue.
    Lit pickBranchLit();                                                      // Return the next decision variable.

public:
    void instantiateLazyDecision(Lit l, int atLevel, CRef reason) override;

    Lit theoryDecisionLit(int theoryID) override{
        if(theoryDecision == lit_Undef){
            theoryDecision = mkLit(newVar(true, false));
        }
        return theoryDecision;
    }

    void
    newDecisionLevel() override;                                                      // Begins a new decision level.
    void uncheckedEnqueue(Lit p, CRef from = CRef_Undef);   // Enqueue a literal. Assumes value of literal is undefined.
    bool enqueue(Lit p,
                 CRef from = CRef_Undef) override;       // Test if fact 'p' contradicts current state, enqueue otherwise.
    void enqueueLazy(Lit p, int level, CRef from = CRef_Undef);

    void setBVTheory(Theory* t) override{
        bvtheory = t;
    }

    Theory* getBVTheory() override{
        return bvtheory;
    }

    void setPBSolver(PB::PBConstraintSolver* t){
        pbsolver = t;
    }

    PB::PBConstraintSolver* getPB() override{
        return pbsolver;
    }

    double& getRandomSeed() override{
        return random_seed;
    }


protected:


    CRef propagate(bool propagate_theories = true);    // Perform unit propagation. Returns possibly conflicting clause.
    void enqueueTheory(Lit l) override;

    void
    enqueueAnyUnqueued() override; //in some solving modes, a theory can sometimes delay enqueing some atoms. Check for any such atoms here.
    bool propagateTheory(vec<Lit>& conflict) override;

    bool solveTheory(vec<Lit>& conflict_out) override;

    void buildReason(Lit p, vec<Lit>& reason) override;

    void backtrackUntil(int level) override;
    //Add a clause to the clause database safely, even if the solver is in the middle of search, propagation, or clause analysis.
    //(In reality, the clause may be added to the database sometime later)

public:
    void addClauseSafely(vec<Lit>& ps) override;

    void setTheoriesEnabled(bool enableTheories){
        disable_theories = !enableTheories;
    }

    bool theoriesEnabled() const{
        return !disable_theories;
    }

    // void unsafeUnassign(Lit p);
    void cancelUntil(int level);                                             // Backtrack until a certain level.
    inline void needsPropagation(int theoryID) override{
        if(theories[theoryID]->unskipable()){
            unskippable_theory_q.insert(theoryID);
        }else{
            if(!in_theory_queue[theoryID]){
                in_theory_queue[theoryID] = true;
                theory_queue.push(theoryID);
                assert(theory_queue.size() <= theories.size());

            }
        }
    }

protected:
    void analyze(CRef confl, vec<Lit>& out_learnt, int& out_btlevel);    // (bt = backtrack)
    void analyzeFinal(CRef confl, Lit skip_lit, vec<Lit>& out_conflict);

    void analyzeFinal(Lit p, vec<Lit>& out_conflict);

    void analyzeFinal(CRef confl, Lit skip_lit, LSet& out_conflict){
        analyzeFinal(confl, skip_lit, tmp_conflict);
        out_conflict.clear();
        out_conflict.insertAll(tmp_conflict);
        tmp_conflict.clear();
    }

    void analyzeFinal(Lit p, LSet& out_conflict){
        analyzeFinal(p, tmp_conflict);
        out_conflict.clear();
        out_conflict.insertAll(tmp_conflict);
        tmp_conflict.clear();
    }

    void analyzeHeuristicDecisions(CRef confl, IntSet<int>& conflicting_heuristics, int max_involved_heuristics,
                                   int minimimum_involved_decision_priority);

    bool litRedundant(Lit p, uint32_t abstract_levels);                       // (helper method for 'analyze()')
    lbool search(int nof_conflicts);                                     // Search for a given number of conflicts.
    lbool solve_();                                           // Main solve method (assumptions given in 'assumptions').
    void reduceDB();                                                      // Reduce the set of learnt clauses.
    void removeSatisfied(vec<CRef>& cs);                           // Shrink 'cs' to contain only non-satisfied clauses.
    void rebuildOrderHeap();

    void rebuildTheoryOrderHeap();

    // Maintaining Variable/Clause activity:
    //
    void
    varDecayActivity(); // Decay all variables with the specified factor. Implemented by increasing the 'bump' value instead.
    void varBumpActivity(Var v, double inc);     // Increase a variable with the current 'bump' value.
    void varBumpActivity(Var v);                 // Increase a variable with the current 'bump' value.
    void
    claDecayActivity(); // Decay all clauses with the specified factor. Implemented by increasing the 'bump' value instead.
    void claBumpActivity(Clause& c);             // Increase a clause with the current 'bump' value.

    void theoryBumpActivity(Heuristic* h){
        if(opt_vsids_both){
            theoryBumpActivity(h, var_inc * opt_theory_vsids_balance);
        }else{
            theoryBumpActivity(h, theory_inc);
        }

    }

    void theoryBumpActivity(Heuristic* h, double increase){

        if((h->getActivity() += increase) > 1e100){
            // Rescale:
            for(Heuristic* t:decision_heuristics)
                t->getActivity() *= 1e-100;
            theory_inc *= 1e-100;

            if(opt_vsids_both){
                for(int i = 0; i < nVars(); i++)
                    activity[i] *= 1e-100;
                var_inc *= 1e-100;
            }
        }

        // Update order_heap with respect to new activity:
        if(theory_order_heap.inHeap(h))
            theory_order_heap.decrease(h);
    }

    inline void theoryDecayActivity(){
        if(opt_vsids_both && opt_use_var_decay_for_theory_vsids){
            varDecayActivity();
        }else{
            theory_inc *= (1 / theory_decay);
        }

    }

    // Operations on clauses:
    //
    CRef attachReasonClause(Lit r, vec<Lit>& ps);

    void attachClause(CRef cr);               // Attach a clause to watcher lists.
    void detachClause(CRef cr, bool strict = false); // Detach a clause to watcher lists.
    void removeClause(CRef cr);               // Detach and free a clause.
    bool
    locked(const Clause& c) const; // Returns TRUE if a clause is a reason for some implication in the current state.
    bool satisfied(const Clause& c) const; // Returns TRUE if a clause is satisfied in the current state.

    void relocAll(ClauseAllocator& to);

public:
    // Misc:
    //
    int decisionLevel() const override; // Gives the current decisionlevel.
    uint32_t abstractLevel(Var x) const; // Used to represent an abstraction of sets of decision levels.

    CRef reason(Var x) const override;

    int level(Var x) const override;

    bool isDecisionReason(CRef r) const;

    CRef decisionReason(Var x) const;

    CRef reasonOrDecision(Var x) const;

    double progressEstimate() const; // DELETE THIS ?? IT'S NOT VERY USEFUL ...

    bool isConstant(Var v) const override{
        return value(v) != l_Undef && (level(v) == 0);
    }

    bool isConstantTrue(Lit l) const{
        return value(l) == l_True && (level(var(l)) == 0);
    }

    // Iterate over clauses and top-level assignments:
    ClauseIterator clausesBegin() const;

    ClauseIterator clausesEnd() const;

    TrailIterator trailBegin() const;

    TrailIterator trailEnd() const;

    /**
     * Lit of known solutions (used for debugging only)
     */
    vec<LSet> known_solutions;
private:
    IntSet<int> swapping_involved_theories;
    IntSet<int> counter_involved_theories;
    uint64_t heuristic_swapping_restarts = 0;
    uint64_t heuristic_conflict_restarts = 0;
    vec<Heuristic*> swapping_uninvolved_pre_theories;
    vec<Heuristic*> swapping_uninvolved_post_theories;
    vec<Heuristic*> swapping_involved_theory_order;

    bool withinBudget() const;

    bool addConflictClause(vec<Lit>& theory_conflict, CRef& confl_out, bool permanent = false) override;

    bool addDelayedClauses(CRef& conflict);


    /**
     * Check that this clause doens't violate any known solutions
     * (used for debugging, only)
     */
    void checkClause(vec<Lit>& clause);

    // Static helpers:
    //
    inline void toSuper(const vec<Lit>& from, vec<Lit>& to){

        to.growTo(from.size());
        to.shrink(to.size() - from.size());
        for(int i = 0; i < from.size(); i++){
            Lit l = from[i];
            assert(local_interface(l));
            to[i] = mkLit(var(l) + super_offset, sign(l));
            assert(super_interface(to[i]));
        }
    }

    inline Var toSuper(Var p){
        assert(p < nVars());
        assert(local_interface(p));
        Var v = p + super_offset;
        assert(v < S->nVars());
        assert(super_interface(v));
        return v;
    }

    inline Lit toSuper(Lit p){
        assert(var(p) < nVars());
        assert(local_interface(p));

        Lit l = mkLit(var(p) + super_offset, sign(p));
        assert(var(l) < S->nVars());
        assert(super_interface(l));
        return l;
    }

    inline Var fromSuper(Var p){
        return p - super_offset;
    }

    inline Lit fromSuper(Lit p){
        return mkLit(var(p) - super_offset, sign(p));
    }

    inline bool local_interface(Lit p) const{
        return var(p) >= min_local && var(p) <= max_local;
    }

    inline bool local_interface(Var v) const{
        return v >= min_local && v <= max_local;
    }

    inline bool super_interface(Lit p) const{
        return var(p) >= min_super && var(p) <= max_super;
    }

    inline bool super_interface(Var v) const{
        return v >= min_super && v <= max_super;
    }

    bool propagateTheorySolver(int theoryID, CRef& confl, vec<Lit>& theory_conflict);

    class SolverDecisionTheory : public Theory {
        //This is a stub theory solver, that is only used to conveniently allow the main solver to make Boolean-decisions
        //as part of the theory decision process, if opt_vsids_solver_as_theory is used
        Solver& S;
        int theoryID = 0;
    public:
        const char* getTheoryType() override{
            return "DecisionTheory";
        }

        SolverDecisionTheory(Solver& S) : S(S){

        }

        Lit decideTheory(CRef& decision_reason) override{
            decision_reason = CRef_Undef;
            return S.pickBranchLit();
        }

        bool supportsDecisions() override{
            return true;
        }

        int getTheoryIndex() const override{
            return theoryID;
        }

        void setTheoryIndex(int id) override{
            theoryID = id;
        }

        void backtrackUntil(int untilLevel) override{

        }

        void newDecisionLevel() override{

        }

        void enqueueTheory(Lit p) override{

        }

        bool propagateTheory(vec<Lit>& conflict) override{
            return true;
        }

        bool solveTheory(vec<Lit>& conflict) override{
            return true;
        }
    };

};

//=================================================================================================
// Implementation of inline methods:

inline CRef Solver::reason(Var x) const{
    CRef r = vardata[x].reason;
    if(isDecisionReason(r))
        return CRef_Undef;
    return r;
}

inline CRef Solver::decisionReason(Var x) const{
    CRef r = vardata[x].reason;
    if(isDecisionReason(r))
        return r;
    return CRef_Undef;
}

inline CRef Solver::reasonOrDecision(Var x) const{
    CRef r = vardata[x].reason;
    return r;
}


inline bool Solver::isDecisionReason(CRef cr) const{
    int marker = CRef_Undef - cr - 1;
    if(marker >= 0 && marker < marker_theory.size()){
        return marker_theory[marker] < 0;
    }else{
        return false;
    }
}

inline int Solver::level(Var x) const{
    return vardata[x].level;
}

inline void Solver::insertVarOrder(Var x){
    if(!order_heap.inHeap(x) && decision[x])
        order_heap.insert(x);
}

inline void Solver::varDecayActivity(){
    var_inc *= (1 / var_decay);
}

inline void Solver::varBumpActivity(Var v){
    varBumpActivity(v, var_inc);
}

inline void Solver::varBumpActivity(Var v, double inc){
    if((activity[v] += inc) > 1e100){
        // Rescale:
        for(int i = 0; i < nVars(); i++)
            activity[i] *= 1e-100;
        var_inc *= 1e-100;

        if(opt_vsids_both){
            for(Heuristic* t:decision_heuristics)
                t->getActivity() *= 1e-100;
            theory_inc *= 1e-100;
        }
    }

    // Update order_heap with respect to new activity:
    if(order_heap.inHeap(v))
        order_heap.decrease(v);


}


inline void Solver::claDecayActivity(){
    cla_inc *= (1 / clause_decay);
}

inline void Solver::claBumpActivity(Clause& c){
    if(c.learnt() && (c.activity() += cla_inc) > 1e20){
        // Rescale:
        for(int i = 0; i < learnts.size(); i++)
            ca[learnts[i]].activity() *= 1e-20;
        cla_inc *= 1e-20;
    }
}

inline void Solver::checkGarbage(void){
    return checkGarbage(garbage_frac);
}

inline void Solver::checkGarbage(double gf){
    if(ca.wasted() > ca.size() * gf)
        garbageCollect();
}

// NOTE: enqueue does not set the ok flag! (only public methods do)
inline bool Solver::enqueue(Lit p, CRef from){
    return value(p) != l_Undef ? value(p) != l_False : (uncheckedEnqueue(p, from), true);
}

inline bool Solver::addClause(const vec<Lit>& ps){
    cancelUntil(0);
    ps.copyTo(add_tmp);
    return addClause_(add_tmp);
}

inline bool Solver::addEmptyClause(){
    cancelUntil(0);
    add_tmp.clear();
    return addClause_(add_tmp);
}

inline bool Solver::addClause(Lit p){
    cancelUntil(0);
    add_tmp.clear();
    add_tmp.push(p);
    return addClause_(add_tmp);
}

inline bool Solver::addClause(Lit p, Lit q){
    cancelUntil(0);
    add_tmp.clear();
    add_tmp.push(p);
    add_tmp.push(q);
    return addClause_(add_tmp);
}

inline bool Solver::addClause(Lit p, Lit q, Lit r){
    cancelUntil(0);
    add_tmp.clear();
    add_tmp.push(p);
    add_tmp.push(q);
    add_tmp.push(r);
    return addClause_(add_tmp);
}

inline bool Solver::locked(const Clause& c) const{
    CRef r = reason(var(c[0]));
    bool isClause = ca.isClause(r);
    if(!isClause)
        return false;
    return value(c[0]) == l_True && ca.isClause(reason(var(c[0]))) && ca.lea(reason(var(c[0]))) == &c;
}

inline void Solver::newDecisionLevel(){
    assert(partially_propagated_levels.size() == 0 || partially_propagated_levels.last() < decisionLevel());
    assert(unskippable_theory_q.size() == 0);
    if(theory_queue.size()){
        //some theories were not propagated
        partially_propagated_levels.push(decisionLevel());
    }
    trail_lim.push(trail.size());
}


inline int Solver::decisionLevel() const{
    return trail_lim.size();
}

inline uint32_t Solver::abstractLevel(Var x) const{
    return 1 << (level(x) & 31);
}

inline lbool Solver::value(Var x) const{
    return assigns[x];
}

inline lbool Solver::value(Lit p) const{
    return assigns[var(p)] ^ sign(p);
}

inline lbool Solver::modelValue(Var x) const{
    return model[x];
}

inline lbool Solver::modelValue(Lit p) const{
    return model[var(p)] ^ sign(p);
}

inline bool Solver::hasModel() const{
    return model.size() > 0;
}

inline int Solver::nAssigns() const{
    return trail.size();
}

inline int Solver::nClauses() const{
    return clauses.size();
}

inline int Solver::nLearnts() const{
    return learnts.size();
}

inline int Solver::nVars() const{
    return vardata.size();
}

inline int Solver::nUnassignedVars() const{
    return (int) nVars() - free_vars.size() - trail.size();
}

inline int Solver::nFreeVars() const{
    return (int) dec_vars - ((trail_lim.size() == 0) ? trail.size() : trail_lim[0]);
}

inline bool Solver::hasNextDecision() const{
    return !order_heap.empty();
}

inline void Solver::setPolarity(Var v, bool b){
    polarity[v] = b;
}

inline void Solver::setDecisionVar(Var v, bool b){
    if(b && !decision[v])
        dec_vars++;
    else if(!b && decision[v])
        dec_vars--;

    decision[v] = b;
    insertVarOrder(v);
}

inline void Solver::setConfBudget(int64_t x){
    if(x < 0){
        conflict_budget = -1;
    }else{
        conflict_budget = conflicts + x;
    }
}

inline void Solver::setPropBudget(int64_t x){
    if(x < 0){
        propagation_budget = -1;
    }else{
        propagation_budget = propagations + x;
    }
}

inline int64_t Solver::getConflictBudget(){
    return conflict_budget > -1 ? (conflict_budget - conflicts) : -1;
}

inline int64_t Solver::getPropagationBudget(){
    return propagation_budget > -1 ? (propagation_budget - propagations) : -1;
}

inline void Solver::interrupt(){
    asynch_interrupt = true;
}

inline void Solver::clearInterrupt(){
    asynch_interrupt = false;
}

inline void Solver::budgetOff(){
    conflict_budget = propagation_budget = -1;
}

inline bool Solver::withinBudget() const{
    return !asynch_interrupt && (conflict_budget < 0 || conflicts < (uint64_t) conflict_budget)
           && (propagation_budget < 0 || propagations < (uint64_t) propagation_budget);
}

// FIXME: after the introduction of asynchronous interrruptions the solve-versions that return a
// pure bool do not give a safe interface. Either interrupts must be possible to turn off here, or
// all calls to solve must return an 'lbool'. I'm not yet sure which I prefer.
inline bool Solver::solve(){
    budgetOff();
    assumptions.clear();
    return solve_() == l_True;
}

inline bool Solver::solve(Lit p){
    budgetOff();
    assumptions.clear();
    assumptions.push(p);
    return solve_() == l_True;
}

inline bool Solver::solve(Lit p, Lit q){
    budgetOff();
    assumptions.clear();
    assumptions.push(p);
    assumptions.push(q);
    return solve_() == l_True;
}

inline bool Solver::solve(Lit p, Lit q, Lit r){
    budgetOff();
    assumptions.clear();
    assumptions.push(p);
    assumptions.push(q);
    assumptions.push(r);
    return solve_() == l_True;
}

inline bool Solver::solve(const vec<Lit>& assumps){
    budgetOff();
    assumps.copyTo(assumptions);
    return solve_() == l_True;
}

inline lbool Solver::solveLimited(const vec<Lit>& assumps){
    assumps.copyTo(assumptions);
    return solve_();
}

inline bool Solver::okay() const{
    return ok;
}

inline ClauseIterator Solver::clausesBegin() const{return ClauseIterator(ca, &clauses[0]);}

inline ClauseIterator Solver::clausesEnd() const{return ClauseIterator(ca, &clauses[clauses.size()]);}

inline TrailIterator Solver::trailBegin() const{return TrailIterator(&trail[0]);}

inline TrailIterator Solver::trailEnd() const{
    return TrailIterator(&trail[decisionLevel() == 0 ? trail.size() : trail_lim[0]]);
}

//=================================================================================================
// Debug etc:

//=================================================================================================
}

#endif
