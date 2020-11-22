/**************************************************************************************************
 The MIT License (MIT)

 Copyright (c) 2014, Sam Bayless

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

#ifndef FSM_THEORY_H_
#define FSM_THEORY_H_

#include "monosat/utils/System.h"
#include "monosat/core/Theory.h"
#include "monosat/dgl/DynamicGraph.h"
#include "monosat/fsm/DynamicFSM.h"
#include "monosat/fsm/FSMAcceptDetector.h"
#include "monosat/fsm/FSMGeneratesDetector.h"
#include "monosat/fsm/FSMTransducesDetector.h"
#include "monosat/fsm/FSMGeneratorAcceptorDetector.h"
#include "monosat/core/SolverTypes.h"
#include "monosat/mtl/Map.h"


#include "monosat/utils/System.h"
#include "monosat/core/Solver.h"

#include <vector>

#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <sstream>

using namespace dgl;
namespace Monosat {


class FSMTheorySolver;

class FSMTheorySolver : public Theory {
public:
    struct Transition {
        Var v = var_Undef;
        Var outerVar = var_Undef;
        int from = -1;
        int to = -1;
        int inputchar;
        int outputchar;
    };

    struct Assignment {
        int isEdge :1;
        int assign :1;
        int isSat  :1;
        int edgeID:29;
        Var var;

        Assignment(bool isEdge, bool assign, int edgeID, Var v) : isEdge(isEdge), assign(assign), isSat(false),
                                                                  edgeID(edgeID), var(v){

        }

        Assignment(bool isEdge, bool assign, bool issat, int edgeID, Var v) : isEdge(isEdge), assign(assign),
                                                                              isSat(issat), edgeID(edgeID), var(v){

        }
    };

    double rnd_seed;
    vec<vec<int>>* strings = nullptr;
private:
    Solver* S;

public:
    int id;

    const char* getTheoryType() override{
        return "FSM";
    }

    vec<lbool> assigns;

    vec<int> n_in_alphabets;//number of transition labels. Transition labels start from 0 (which is the non-consuming epsilon transition) and go to n_labels-1.
    vec<int> n_out_alphabets;

    vec<vec<vec<Transition>>> edge_labels;
    vec<DynamicFSM*> g_unders;
    vec<DynamicFSM*> g_overs;

    vec<FSMAcceptDetector*> accepts;
    vec<FSMGeneratesDetector*> generates;
    vec<FSMTransducesDetector*> transduces;

    // If an fsm ID is frozen, no new nodes or transition can be added to it
    vec<bool> frozenFSM;

    /**
     * The cutgraph is (optionally) used for conflict analysis by some graph theories.
     * It has two edges for every edge in the real graph (with indices edgeID*2 and edgeID*2+1).
     * If edge ID is assigned to FALSE, then edge ID*2 is assigned to enabled in the cutgraph, and
     * edge ID*2+1 is disabled.
     * Otherwise, if edge ID is unassigned or true, then edge ID*2 is disabled in the cutgraph, and
     * edge ID*2+1 is enabled.
     */



    vec<Assignment> trail;
    vec<int> trail_lim;

    vec<FSMGeneratorAcceptorDetector*> gen_accept_lit_map;

    struct StringAcceptConstraint {
        int fsmID;
        int strID;
        int sourceState;
        int acceptState;
        Var outer_var;
    };

    vec<StringAcceptConstraint> unimplemented_string_accept_constraints;

    struct GeneratorAcceptConstraint {
        int generatorFsmID;
        int acceptorFsmID;
        int strID;// not currently used (must be -1)
        int generatorSourceState;
        int generatorAcceptState;
        int acceptorSourceState;
        int acceptorAcceptState;
        bool generator_is_deterministic;
        Var outer_var;
    };

    vec<GeneratorAcceptConstraint> unimplemented_generator_accept_constraints;

    struct TransduceConstraint {
        int fsmID;
        int source;
        int dest;
        int strID;
        int strID2;
        Var outer_var;
    };

    vec<TransduceConstraint> unimplemented_transduce_constraints;

    struct GenerateStringConstraint {
        int fsmID;
        int source;
        int strID;
        Var outer_var;
    };

    vec<GenerateStringConstraint> unimplemented_generate_string_constraints;


    template<typename... TTypes>
    class tuple_hash {
    private:
        typedef std::tuple<TTypes...> Tuple;

        template<int N>
        size_t operator()(const Tuple& value) const{return 0;}

        template<int N, typename THead, typename... TTail>
        size_t operator()(const Tuple& value) const{
            constexpr int Index = N - sizeof...(TTail) - 1;
            return hash_combine(std::hash<THead>()(std::get<Index>(value)), operator()<N, TTail...>(value));
        }

        // adapted from boost::hash_combine
        static std::size_t hash_combine(const std::size_t& h, const std::size_t& v){
            return h ^ (v + 0x9e3779b9 + (h << 6) + (h >> 2));
        }

    public:
        size_t operator()(const Tuple& value) const{
            return operator()<sizeof...(TTypes), TTypes...>(value);
        }
    };

    Map<std::tuple<int, int, int, int, bool>, FSMGeneratorAcceptorDetector*, tuple_hash<int, int, int, int, bool>> gen_accept_map;
public:
    vec<FSMDetector*> detectors;
    vec<FSMAcceptDetector*> reach_detectors;


    vec<int> marker_map;


    bool requiresPropagation = true;

    vec<char> seen;
    vec<int> to_visit;
    vec<Lit> tmp_clause;
    //Data about local theory variables, and how they connect to the sat solver's variables
    struct VarData {
        int isEdge :1;
        int occursPositive :1;
        int occursNegative :1;
        int isSatisfied:1;
        //the detector this variable belongs to, or its edge number, if it is an edge variable
        int detector_edge:28;
        int input;
        int output;
        int fsmID;
        Var solverVar;
    };

    vec<VarData> vars;
    int theory_index = 0;
public:

    double mctime = 0;
    double reachtime = 0;
    double unreachtime = 0;
    double pathtime = 0;
    double propagationtime = 0;
    int64_t stats_propagations = 0;
    int64_t propagations = 0;
    int64_t stats_num_conflicts = 0;
    int64_t stats_decisions = 0;
    int64_t stats_num_reasons = 0;

    double reachupdatetime = 0;
    double unreachupdatetime = 0;
    double stats_initial_propagation_time = 0;
    double stats_decision_time = 0;
    double stats_reason_initial_time = 0;
    double stats_reason_time = 0;
    int64_t num_learnt_paths = 0;
    int64_t learnt_path_clause_length = 0;
    int64_t num_learnt_cuts = 0;
    int64_t learnt_cut_clause_length = 0;
    int64_t stats_pure_skipped = 0;
    int64_t stats_mc_calls = 0;
    int64_t stats_propagations_skipped = 0;


    FSMTheorySolver(Solver* S_, int _id = -1) :
            S(S_), id(_id){
        strings = new vec<vec<int>>();
        rnd_seed = opt_random_seed;
        S->addTheory(this);
    }

    ~FSMTheorySolver() override{
        if(this->strings){
            //delete(this->strings);//why is this problematic?
        }
        for(DynamicFSM* f:g_unders){
            delete (f);
        }
        for(DynamicFSM* f:g_overs){
            delete (f);
        }
    }

    Solver* getSolver(){
        return S;
    }

    void printStats(int detailLevel) override{
        printf("FSM %d stats:\n", getGraphID());
        if(stats_decisions > 0){
            printf("FSM decisions: %" PRId64 "\n", stats_decisions);
        }
        if(detailLevel > 0){
            for(FSMDetector* d : detectors)
                d->printStats();
        }


        fflush(stdout);
    }

    void writeTheoryWitness(std::ostream& write_to) override{

        for(FSMDetector* d : detectors){
            write_to << "Graph " << this->getGraphID() << ", detector " << d->getID() << ":\n";
            d->printSolution(write_to);
        }
    }

    inline int getTheoryIndex() const override{
        return theory_index;
    }

    inline void setTheoryIndex(int id) override{
        theory_index = id;
    }

    inline int getGraphID(){
        return id;
    }

    inline int inAlphabet(int fsmID) const{
        return n_in_alphabets[fsmID];
    }

    inline int outAlphabet(int fsmID) const{
        return n_out_alphabets[fsmID];
    }

    void setAlphabets(int fsmID, int inputs, int outputs){
        assert(g_unders[fsmID]);
        DynamicFSM& g_under = *g_unders[fsmID];
        DynamicFSM& g_over = *g_overs[fsmID];
        inputs += 1;
        outputs += 1;
        assert(inputs >= n_in_alphabets[fsmID]);
        while(inputs > inAlphabet(fsmID)){
            n_in_alphabets[fsmID]++;
            g_under.addInCharacter();
            g_over.addInCharacter();
        }
        while(outputs > outAlphabet(fsmID)){
            n_out_alphabets[fsmID]++;
            g_under.addOutCharacter();
            g_over.addOutCharacter();
        }
        for(int i = 0; i < edge_labels[fsmID].size(); i++){
            edge_labels[fsmID][i].growTo(inAlphabet(fsmID) * outAlphabet(fsmID));
        }
    }

    int newString(vec<int>& str){
        strings->push();
        str.copyTo(strings->last());
        return strings->size() - 1;
    }

    int nStrings() const{
        return strings->size();
    }

    inline bool isEdgeVar(Var v) const{
        assert(v < vars.size());
        return vars[v].isEdge;
    }

    inline int getEdgeID(Var v) const{
        assert(isEdgeVar(v));
        return vars[v].detector_edge;
    }

    inline int getFsmID(Var v) const{

        return vars[v].fsmID;
    }

    inline Transition& getTransition(int fsmID, int edgeID, int input, int output){
        return edge_labels[fsmID][edgeID][input + output * inAlphabet(fsmID)];
    }

    inline int getInput(Var v) const{
        assert(isEdgeVar(v));
        return vars[v].input;
    }

    inline int getOutput(Var v) const{
        assert(isEdgeVar(v));
        return vars[v].output;
    }

    inline int getDetector(Var v) const{
        assert(!isEdgeVar(v));
        return vars[v].detector_edge;
    }

    inline Var getTransitionVar(int fsmID, int edgeID, int inputChar, int outputChar){
        assert(inputChar >= 0);
        assert(outputChar >= 0);
        Var v = getTransition(fsmID, edgeID, inputChar, outputChar).v;
        assert(v < vars.size());
        return v;
    }

    void makeEqual(Lit l1, Lit l2){
        Lit o1 = toSolver(l1);
        Lit o2 = toSolver(l2);
        S->addClause(~o1, o2);
        S->addClause(o1, ~o2);
    }

    void makeEqualInSolver(Lit l1, Lit l2){
        S->addClause(~l1, l2);
        S->addClause(l1, ~l2);
    }

    void addClause(Lit l1){
        Lit o1 = toSolver(l1);
        S->addClause(o1);
    }

    void addClause(Lit l1, Lit l2){
        Lit o1 = toSolver(l1);
        Lit o2 = toSolver(l2);
        S->addClause(o1, o2);
    }

    void addClause(Lit l1, Lit l2, Lit l3){
        Lit o1 = toSolver(l1);
        Lit o2 = toSolver(l2);
        Lit o3 = toSolver(l3);
        S->addClause(o1, o2, o3);
    }

    void addClause(vec<Lit>& c){
        tmp_clause.clear();
        c.copyTo(tmp_clause);
        toSolver(tmp_clause);
        S->addClause(tmp_clause);
    }

    void addClauseSafely(vec<Lit>& c){
        tmp_clause.clear();
        c.copyTo(tmp_clause);
        toSolver(tmp_clause);

        S->addClauseSafely(tmp_clause);
    }

    Var newAuxVar(int forDetector = -1, bool connectToTheory = false){
        Var s = S->newVar();
        return newVar(s, forDetector, -1, -1, -1, false, connectToTheory);
    }

    Var newVar(Var solverVar, int detector_edge, int fsmID = -1, int label = -1, int output = -1, bool isEdge = false,
               bool connectToTheory = true){
        while(S->nVars() <= solverVar)
            S->newVar();
        Var v = vars.size();

        vars.push();
        vars[v].isEdge = isEdge;
        vars[v].fsmID = fsmID;
        vars[v].detector_edge = detector_edge;
        vars[v].solverVar = solverVar;
        vars[v].input = label;
        vars[v].output = output;
        vars[v].isSatisfied = 0;
        vars[v].occursPositive = 0;
        vars[v].occursNegative = 0;
        assigns.push(l_Undef);
        if(connectToTheory){
            S->setTheoryVar(solverVar, getTheoryIndex(), v);
            assert(toSolver(v) == solverVar);
        }
        if(isEdge){
            assert(label > -1);
            assert(detector_edge > -1);
        }
        if(!isEdge && detector_edge >= 0)
            detectors[detector_edge]->addVar(v);
        return v;
    }

    inline int level(Var v){
        return S->level(toSolver(v));
    }

    inline int decisionLevel(){
        return trail_lim.size();
    }

    inline int nVars() const{
        return vars.size();
    }

    inline Var toSolver(Var v){
        assert(v < vars.size());

        return vars[v].solverVar;
    }

    inline Lit toSolver(Lit l){

        return mkLit(vars[var(l)].solverVar, sign(l));
    }

    void toSolver(vec<Lit>& c){
        for(int i = 0; i < c.size(); i++){
            c[i] = toSolver(c[i]);
        }
    }

    inline lbool value(Var v){
        if(assigns[v] != l_Undef)
            assert(S->value(toSolver(v)) == assigns[v]);

        return assigns[v]; //S->value(toSolver(v));
    }

    inline lbool value(Lit l){
        if(assigns[var(l)] != l_Undef){
            assert(S->value(toSolver(l)) == (assigns[var(l)] ^ sign(l)));
        }
        return assigns[var(l)] ^ sign(l);
    }

    inline lbool dbg_value(Var v){
        return S->value(toSolver(v));
    }

    inline lbool dbg_value(Lit l){
        return S->value(toSolver(l));
    }

    inline bool enqueue(Lit l, CRef reason){
        assert(assigns[var(l)] == l_Undef);

        Lit sl = toSolver(l);
        if(S->enqueue(sl, reason)){
            enqueueTheory(l);
            return true;
        }else{
            return false;
        }
    }

    bool isFrozen(int fsmID){
        frozenFSM.growTo(fsmID + 1, false);
        return frozenFSM[fsmID];
    }

    void freezeFSM(int fsmID){
        frozenFSM.growTo(fsmID + 1, false);
        frozenFSM[fsmID] = true;
    }

    int newNode(int fsmID){
        if(isFrozen(fsmID)){
            throw std::runtime_error("Cannot add nodes to fsm "
                                     + std::to_string(fsmID) +
                                     " after constraints implemented");
        }
        assert(g_overs[fsmID]);
        g_overs[fsmID]->addNode();

        seen.growTo(nNodes(fsmID));

        return g_unders[fsmID]->addNode();
    }

    void newNodes(int fsmID, int n){
        for(int i = 0; i < n; i++)
            newNode(fsmID);
    }

    int nNodes(int fsmID){
        return g_overs[fsmID]->nodes();
    }

    bool isNode(int fsmID, int n){
        return n >= 0 && n < nNodes(fsmID);
    }

    void backtrackUntil(int level) override{

        bool changed = false;
        //need to remove and add edges in the two graphs accordingly.
        if(trail_lim.size() > level){

            int stop = trail_lim[level];
            for(int i = trail.size() - 1; i >= trail_lim[level]; i--){

                Assignment& e = trail[i];
                assert(assigns[e.var] != l_Undef);
                if(e.isSat){
                    Lit l = mkLit(e.var, !e.assign);
                    vars[var(l)].isSatisfied = 0;
                    detectors[getDetector(e.var)]->setSatisfied(l, false);

                }else if(e.isEdge){
                    assert(dbg_value(e.var) == value(e.var));
                    int fsmID = getFsmID(e.var);
                    int edgeID = getEdgeID(e.var); //e.var-min_edge_var;
                    int input = getInput(e.var);
                    int output = getOutput(e.var);
                    assert(edgeID == e.edgeID);

                    if(e.assign){
                        g_unders[fsmID]->disableTransition(edgeID, input, output);

                    }else{
                        g_overs[fsmID]->enableTransition(edgeID, input, output);

                    }
                }else{
                    //This is a reachability literal
                    detectors[getDetector(e.var)]->unassign(mkLit(e.var, !e.assign));
                }
                assigns[e.var] = l_Undef;
                changed = true;
            }
            trail.shrink(trail.size() - stop);
            trail_lim.shrink(trail_lim.size() - level);
            assert(trail_lim.size() == level);

            if(changed){
                requiresPropagation = true;
            }

            for(FSMDetector* d : detectors){
                d->backtrack(level);
            }
        }

        assert(dbg_graphsUpToDate());

    };

    bool supportsDecisions() override{
        return true;
    }

    Lit decideTheory(CRef& decision_reason) override{
        decision_reason = CRef_Undef;
        if(!opt_decide_theories)
            return lit_Undef;
        double start = rtime(1);

        for(int i = 0; i < detectors.size(); i++){
            FSMDetector* r = detectors[i];
            Lit l = r->decide(decisionLevel());
            if(l != lit_Undef){
                stats_decisions++;
                r->stats_decisions++;
                stats_decision_time += rtime(1) - start;
                return toSolver(l);
            }
        }
        stats_decision_time += rtime(1) - start;
        return lit_Undef;
    }

    void backtrackUntil(Lit p){
        //need to remove and add edges in the two graphs accordingly.
        assert(value(p) == l_True);
        int i = trail.size() - 1;
        for(; i >= 0; i--){
            Assignment e = trail[i];
            if(var(p) == e.var){
                assert(sign(p) != e.assign);
                break;
            }
            if(e.isSat){
                Lit l = mkLit(e.var, !e.assign);
                vars[var(l)].isSatisfied = 0;
                detectors[getDetector(e.var)]->setSatisfied(l, false);
            }else if(e.isEdge){
                int fsmID = getFsmID(e.var);
                int edgeID = getEdgeID(e.var); //e.var-min_edge_var;
                int input = getInput(e.var);
                int output = getOutput(e.var);
                assert(assigns[e.var] != l_Undef);
                assigns[e.var] = l_Undef;
                if(e.assign){
                    g_unders[fsmID]->disableTransition(edgeID, input, output);
                }else{
                    g_overs[fsmID]->enableTransition(edgeID, input, output);
                }
            }else{

                assigns[e.var] = l_Undef;
                detectors[getDetector(e.var)]->unassign(mkLit(e.var, !e.assign));
            }
        }

        trail.shrink(trail.size() - (i + 1));

        requiresPropagation = true;//is this really needed?

        for(FSMDetector* d : detectors){
            d->backtrack(this->decisionLevel());
        }

    };

    void newDecisionLevel() override{
        trail_lim.push(trail.size());
    };

    void buildReason(Lit p, vec<Lit>& reason) override{
        CRef marker = S->reason(var(toSolver(p)));
        assert(marker != CRef_Undef);
        int pos = CRef_Undef - marker;
        int d = marker_map[pos];
        //double initial_start = rtime(1);
        double start = rtime(1);
        backtrackUntil(p);

        assert(d < detectors.size());
        detectors[d]->buildReason(p, reason, marker);
        toSolver(reason);
        double finish = rtime(1);
        stats_reason_time += finish - start;
        stats_num_reasons++;
        //stats_reason_initial_time+=start-initial_start;

    }


    bool dbg_graphsUpToDate(){

        return true;
    }

    bool isFsmDefined(int fsmID){
        return fsmID >= 0 && fsmID < g_overs.size() && g_overs[fsmID] != nullptr;
    }

private:
    void ensureStateDefined(int fsmId, int state){
        if(!isFsmDefined(fsmId)){
            throw std::runtime_error("Undefined FSM with ID " + std::to_string(fsmId));
        }
        while(nNodes(fsmId) <= state){
            newNode(fsmId);
        }
    }

public:
    void implementConstraints(){
        if(!S->okay())
            return;

        // before implementing any constraints, ensure all states are defined
        for(auto& c:unimplemented_string_accept_constraints){
            ensureStateDefined(c.fsmID, c.sourceState);
            ensureStateDefined(c.fsmID, c.acceptState);
        }

        for(auto& c:unimplemented_generator_accept_constraints){
            ensureStateDefined(c.generatorFsmID, c.generatorSourceState);
            ensureStateDefined(c.generatorFsmID, c.generatorAcceptState);

            ensureStateDefined(c.acceptorFsmID, c.acceptorSourceState);
            ensureStateDefined(c.acceptorFsmID, c.acceptorAcceptState);
        }

        for(auto& c:unimplemented_transduce_constraints){
            ensureStateDefined(c.fsmID, c.source);
            ensureStateDefined(c.fsmID, c.dest);
        }

        for(auto& c:unimplemented_generate_string_constraints){
            ensureStateDefined(c.fsmID, c.source);
        }

        // implement constraints

        for(auto& c:unimplemented_string_accept_constraints){
            implementStringAcceptPrivate(c);
        }
        unimplemented_string_accept_constraints.clear();

        for(auto& c:unimplemented_generator_accept_constraints){
            implementComposeAccept(c);
        }
        unimplemented_generator_accept_constraints.clear();

        for(auto& c:unimplemented_transduce_constraints){
            implementTransduce(c);
        }
        unimplemented_transduce_constraints.clear();

        for(auto& c:unimplemented_generate_string_constraints){
            implementGenerateLit(c);
        }
        unimplemented_generate_string_constraints.clear();
    }

    void preprocess() override{
        implementConstraints();
        for(int i = 0; i < detectors.size(); i++){
            detectors[i]->preprocess();
        }
    }

    void setLiteralOccurs(Lit l, bool occurs) override{
        if(isEdgeVar(var(l))){
            //don't do anything
        }else{
            if(!sign(l)){
                vars[var(l)].occursPositive = occurs;
            }else{
                vars[var(l)].occursNegative = occurs;
            }
            //this is a graph property detector var

            detectors[getDetector(var(l))]->setOccurs(l, occurs);

        }

    }

    void enqueueTheory(Lit l) override{
        Var v = var(l);

        int lev = level(v);

        assert(decisionLevel() <= lev);

        while(lev > trail_lim.size()){
            newDecisionLevel();
        }

        if(assigns[var(l)] != l_Undef){
            return;            //this is already enqueued.
        }
        assert(assigns[var(l)] == l_Undef);
        assigns[var(l)] = sign(l) ? l_False : l_True;
        requiresPropagation = true;

        if(isEdgeVar(var(l))){

            //this is an edge assignment
            int fsmID = getFsmID(var(l));
            int edgeID = getEdgeID(var(l)); //v-min_edge_var;
            int input = getInput(var(l));
            int output = getOutput(var(l));
            assert(getTransition(fsmID, edgeID, input, output).v == var(l));

            trail.push(Assignment(true, !sign(l), edgeID, v));

            if(!sign(l)){
                g_unders[fsmID]->enableTransition(edgeID, input, output);
            }else{
                g_overs[fsmID]->disableTransition(edgeID, input, output);
            }

        }else{

            trail.push(Assignment(false, !sign(l), -1, v));
            //this is an assignment to a non-edge atom. (eg, a reachability assertion)
            detectors[getDetector(var(l))]->assign(l);

        }

    }

    //mark an atom as satisfied in the theory, so it doesn't need to be tracked in the future
    void enqueueSat(Lit l){
        assert(value(l) == l_True);//the literal must be assigned true
        if(opt_detect_satisfied_predicates){
            if(!vars[var(l)].isSatisfied){

                while(decisionLevel() < S->decisionLevel()){
                    newDecisionLevel();
                }

                vars[var(l)].isSatisfied = 1;
                int d = getDetector(var(l));
                detectors[d]->setSatisfied(l, true);
            }
        }
    }

    bool isSatisfied(Lit l){
        assert(!(vars[var(l)].isSatisfied && value(l) == l_Undef));
        // if the var is marked satisfied, it SHOULD be the case
        // that it is assigned true.
        return vars[var(l)].isSatisfied;
    }

    bool literalOccurs(Lit l){
        if(!sign(l)){
            return vars[var(l)].occursPositive || (value(l) == l_False);//false, not true, here
        }else{
            return vars[var(l)].occursNegative || (value(l) == l_False);
        }
    }

    bool litIsRelevant(Lit l){
        return literalOccurs(~l) && !isSatisfied(l);
    }

    bool propagateTheory(vec<Lit>& conflict) override{
        return propagateTheory(conflict, false);
    }

    bool propagateTheory(vec<Lit>& conflict, bool force_propagation){

        stats_propagations++;

        if(!requiresPropagation){
            stats_propagations_skipped++;
            assert(dbg_graphsUpToDate());
            return true;
        }

        propagations++;

        if(propagations > 1 && (!force_propagation && (propagations % opt_fsm_prop_skip != 0))){
            stats_propagations_skipped++;
            return true;
        }

        bool any_change = false;
        double startproptime = rtime(1);


        conflict.clear();
        //Can probably speed this up alot by a) constant propagating reaches that I care about at level 0, and b) Removing all detectors for nodes that appear only in the opposite polarity (or not at all) in the cnf.
        //That second one especially.

        //At level 0, need to propagate constant reaches/source nodes/edges...

        assert(dbg_graphsUpToDate());

        for(int d = 0; d < detectors.size(); d++){
            assert(conflict.size() == 0);
            bool r = detectors[d]->propagate(conflict);
            if(!r){
                stats_num_conflicts++;
                toSolver(conflict);
                propagationtime += rtime(1) - startproptime;
                return false;
            }else{
                assert(conflict.size() == 0);
            }
        }


        requiresPropagation = false;
        for(int i = 0; i < nFsms(); i++){
            if(g_unders[i]){
                g_unders[i]->clearChanged();
                g_overs[i]->clearChanged();

                g_unders[i]->clearHistory();
                g_overs[i]->clearHistory();
            }
        }


        double elapsed = rtime(1) - startproptime;
        propagationtime += elapsed;

        return true;
    };

    bool solveTheory(vec<Lit>& conflict) override{
        requiresPropagation = true;        //Just to be on the safe side... but this shouldn't really be required.
        bool ret = propagateTheory(conflict, true);
        //Under normal conditions, this should _always_ hold
        // (as propagateTheory should have been called and checked by the parent solver before getting to this point).
        assert(ret || opt_fsm_prop_skip > 1);
        return ret;
    };

    void drawFull(int from, int to){

    }

    void drawFSM(int fsmID, int source, int dest, bool over){
        if(over){
            g_overs[fsmID]->draw(source, dest);
        }else{
            g_unders[fsmID]->draw(source, dest);
        }
    }

    bool check_solved() override{

        if(unimplemented_generate_string_constraints.size() > 0 ||
           unimplemented_transduce_constraints.size() > 0 ||
           unimplemented_generator_accept_constraints.size() > 0 ||
           unimplemented_string_accept_constraints.size() > 0){
            return false;
        }

        for(int fsmID = 0; fsmID < nFsms(); fsmID++){
            if(!g_unders[fsmID])
                continue;
            DynamicFSM& g_under = *g_unders[fsmID];
            DynamicFSM& g_over = *g_overs[fsmID];
            for(int edgeID = 0; edgeID < edge_labels[fsmID].size(); edgeID++){
                for(int input = 0; input < inAlphabet(fsmID); input++){
                    for(int output = 0; output < outAlphabet(fsmID); output++){
                        if(getTransition(fsmID, edgeID, input, output).v < 0)
                            continue;
                        Var v = getTransition(fsmID, edgeID, input, output).v;
                        if(v == var_Undef)
                            continue;
                        lbool val = value(v);
                        if(val == l_Undef){
                            return false;
                        }

                        if(val == l_True){
                            if(!g_under.transitionEnabled(edgeID, input, output)){
                                return false;
                            }
                            if(!g_over.transitionEnabled(edgeID, input, output)){
                                return false;
                            }
                        }else{
                            if(g_under.transitionEnabled(edgeID, input, output)){
                                return false;
                            }
                            if(g_over.transitionEnabled(edgeID, input, output)){
                                return false;
                            }
                        }
                    }
                }
            }
        }
        for(int i = 0; i < detectors.size(); i++){
            if(!detectors[i]->checkSatisfied()){
                return false;
            }
        }
        return true;
    }

    bool dbg_solved(){

        return true;
    }

    void drawCurrent(){

    }

    int nFsms() const{
        return g_overs.size();
    }

    int nEdges(int fsmID){
        return edge_labels[fsmID].size();
    }

    CRef newReasonMarker(int detectorID, bool is_decision = false){
        CRef reasonMarker = S->newReasonMarker(this, is_decision);
        int mnum = CRef_Undef - reasonMarker;
        marker_map.growTo(mnum + 1);
        marker_map[mnum] = detectorID;
        return reasonMarker;
    }

    bool hasFSM(int fsmID){
        return fsmID >= 0 && fsmID < g_unders.size() && g_unders[fsmID] != nullptr;
    }

    int newFSM(int fsmID = -1){
        if(fsmID < 0){
            fsmID = g_unders.size();
        }
        g_unders.growTo(fsmID + 1, nullptr);
        g_overs.growTo(fsmID + 1, nullptr);
        if(g_unders[fsmID]){
            throw std::runtime_error("Redefined fsm " + std::to_string(fsmID));
        }
        edge_labels.growTo(fsmID + 1);
        g_unders[fsmID] = new DynamicFSM(fsmID);
        g_overs[fsmID] = new DynamicFSM(fsmID);
        n_in_alphabets.growTo(fsmID + 1,
                              1);//number of transition labels. Transition labels start from 0 (which is the non-consuming epsilon transition) and go to n_labels-1.
        n_out_alphabets.growTo(fsmID + 1, 1);
        return fsmID;
    }

    Lit newTransition(int fsmID, int from, int to, int input, int output, Var outerVar = var_Undef){
        if(isFrozen(fsmID)){
            throw std::runtime_error("Cannot add transition to fsm "
                                     + std::to_string(fsmID) +
                                     " after constraints implemented");
        }
        assert(outerVar != var_Undef);
        assert(input >= 0);
        assert(outerVar != var_Undef);
        if(from == to && input == 0){
            //don't add this transition; self-looping e transitions have no effect.
            return lit_Undef;
        }
        assert(g_unders[fsmID]);
        DynamicFSM& g_under = *g_unders[fsmID];
        DynamicFSM& g_over = *g_overs[fsmID];
        int edgeID = -1;
        if(g_under.states() > from && g_under.states() > to && (edgeID = g_under.getEdge(from, to)) > -1){
            if(g_over.inAlphabet() > input && g_over.outAlphabet() > output &&
               g_over.transitionEnabled(edgeID, input, output)){
                //we already have this transition implemented
                Var ov = getTransition(fsmID, edgeID, input, output).outerVar;
                assert(ov != var_Undef);
                assert(getTransition(fsmID, edgeID, input, output).from == from);
                assert(getTransition(fsmID, edgeID, input, output).to == to);
                makeEqualInSolver(mkLit(outerVar), mkLit(ov));
                return lit_Undef;
            }

        }else{


        }


        g_under.addTransition(from, to, edgeID, input, output, false);
        edgeID = g_over.addTransition(from, to, edgeID, input, output, true);
        Var v = newVar(outerVar, edgeID, fsmID, input, output, true);
        assert(input < inAlphabet(fsmID));
        assert(output < outAlphabet(fsmID));


        edge_labels[fsmID].growTo(edgeID + 1);
        edge_labels[fsmID][edgeID].growTo(inAlphabet(fsmID) * outAlphabet(fsmID));


        getTransition(fsmID, edgeID, input, output).v = v;
        getTransition(fsmID, edgeID, input, output).outerVar = outerVar;
        getTransition(fsmID, edgeID, input, output).from = from;
        getTransition(fsmID, edgeID, input, output).to = to;

        return mkLit(v, false);
    }

    void printSolution() override{

        for(auto* d : detectors){
            assert(d);
            d->printSolution();
        }
    }

    void setStrings(vec<vec<int>>* strings){
        if(this->strings){
            delete (this->strings);
        }
        this->strings = strings;
    }

    void addAcceptLit(int fsmID, int source, int reach, int strID, Var outer_var){
        unimplemented_string_accept_constraints.push({fsmID, strID, source, reach, outer_var});
    }

private:
    void implementStringAcceptPrivate(const StringAcceptConstraint& constraint){
        assert(g_unders[constraint.fsmID]);
        DynamicFSM& g_under = *g_unders[constraint.fsmID];
        DynamicFSM& g_over = *g_overs[constraint.fsmID];
        for(int i = 0; i < (*strings)[constraint.strID].size(); i++){
            int l = (*strings)[constraint.strID][i];
            if(l < 0 || l >= g_overs[constraint.fsmID]->inAlphabet()){
                throw std::runtime_error("String has letter " + std::to_string(l) + " out of range for fsm " +
                                         std::to_string(constraint.fsmID));
            }
        }
        while(constraint.sourceState >= g_over.nodes() || constraint.acceptState >= g_over.nodes()){
            g_over.addNode();
            g_under.addNode();
        }

        accepts.growTo(constraint.sourceState + 1);
        if(!accepts[constraint.sourceState]){
            accepts[constraint.sourceState] = new FSMAcceptDetector(detectors.size(), this,
                                                                    g_under, g_over, constraint.sourceState, *strings,
                                                                    drand(rnd_seed));
            detectors.push(accepts[constraint.sourceState]);
        }
        accepts[constraint.sourceState]->addAcceptLit(constraint.acceptState, constraint.strID, constraint.outer_var);
        freezeFSM(constraint.fsmID);
    }

public:
    void addGenerateLit(int fsmID, int source, int strID, Var outer_var){
        unimplemented_generate_string_constraints.push({
                                                               fsmID, source, strID, outer_var
                                                       });
    }

private:
    void implementGenerateLit(const GenerateStringConstraint& constraint){
        assert(g_unders[constraint.fsmID]);
        DynamicFSM& g_under = *g_unders[constraint.fsmID];
        DynamicFSM& g_over = *g_overs[constraint.fsmID];
        generates.growTo(constraint.source + 1);
        if(!generates[constraint.source]){
            generates[constraint.source] = new FSMGeneratesDetector(detectors.size(), this, g_under, g_over,
                                                                    constraint.source, *strings, drand(rnd_seed));
            detectors.push(generates[constraint.source]);
        }
        generates[constraint.source]->addGeneratesLit(constraint.strID, constraint.outer_var);
        freezeFSM(constraint.fsmID);
    }

public:

    void addTransduceLit(int fsmID, int source, int dest, int strID, int strID2, Var outer_var){
        unimplemented_transduce_constraints.push({fsmID, source, dest, strID, strID2, outer_var});
    }

private:
    void implementTransduce(const TransduceConstraint& constraint){
        assert(g_unders[constraint.fsmID]);
        DynamicFSM& g_under = *g_unders[constraint.fsmID];
        DynamicFSM& g_over = *g_overs[constraint.fsmID];
        transduces.growTo(constraint.source + 1);
        if(!transduces[constraint.source]){
            transduces[constraint.source] =
                    new FSMTransducesDetector(detectors.size(), this, g_under, g_over, constraint.source, *strings,
                                              drand(rnd_seed));
            detectors.push(transduces[constraint.source]);
        }
        transduces[constraint.source]->addTransducesLit(constraint.dest, constraint.strID, constraint.strID2,
                                                        constraint.outer_var);
        freezeFSM(constraint.fsmID);
    }

public:

    void addComposeAcceptLit(int generatorFsmId, int acceptorFsmId, int generatorSourceState,
                             int generatorAcceptState, int acceptorSourceState, int acceptorAcceptState, int strID,
                             Var reachVar, bool generator_is_deterministic = false){
        if(strID >= 0){
            throw std::invalid_argument("String inputs are not yet supported in compositions (strId must be -1)");
        }
        unimplemented_generator_accept_constraints.push({generatorFsmId, acceptorFsmId, strID, generatorSourceState,
                                                         generatorAcceptState, acceptorSourceState, acceptorAcceptState,
                                                         generator_is_deterministic, reachVar});
    }

private:
    void implementComposeAccept(const GeneratorAcceptConstraint& constraint){
        //for now, only linear generator/acceptor compositions are supported
        if(constraint.strID >= 0){
            throw std::invalid_argument("String inputs are not yet supported in compositions (strId must be -1)");
        }

        if(!g_overs[constraint.generatorFsmID] || !g_overs[constraint.acceptorFsmID]){
            throw std::invalid_argument("Undefined FSM ID");

        }
        while(constraint.generatorSourceState >= g_overs[constraint.generatorFsmID]->nodes()
              || constraint.generatorAcceptState >= g_overs[constraint.generatorFsmID]->nodes()){
            g_overs[constraint.generatorFsmID]->addNode();
            g_unders[constraint.generatorFsmID]->addNode();
        }

        while(constraint.acceptorSourceState >= g_overs[constraint.acceptorFsmID]->nodes()
              || constraint.acceptorAcceptState >= g_overs[constraint.acceptorFsmID]->nodes()){
            g_overs[constraint.acceptorFsmID]->addNode();
            g_unders[constraint.acceptorFsmID]->addNode();
        }

        if(constraint.generatorSourceState >= g_overs[constraint.generatorFsmID]->states()
           || constraint.acceptorSourceState >= g_overs[constraint.acceptorFsmID]->states() ||
           constraint.generatorAcceptState >= g_overs[constraint.generatorFsmID]->states()
           || constraint.acceptorAcceptState >= g_overs[constraint.acceptorFsmID]->states()){
            throw std::invalid_argument("Undefined fsm state in formula");

        }
        if(!g_overs[constraint.generatorFsmID]->isGenerator()){
            throw std::invalid_argument("Only compositions of linear generators with FSAs are supported");

        }
        if(!g_overs[constraint.generatorFsmID]->isLinear()){
            throw std::invalid_argument("Only compositions of linear generators with FSAs are supported");

        }
        if(!g_overs[constraint.acceptorFsmID]->isAcceptor()){
            throw std::invalid_argument("Only compositions of linear generators with FSAs are supported");

        }
        if(g_overs[constraint.generatorFsmID]->out_alphabet != g_overs[constraint.acceptorFsmID]->in_alphabet){
            throw std::invalid_argument(
                    "Size of output alphabet of first fsm (was " +
                    std::to_string(g_overs[constraint.generatorFsmID]->out_alphabet) +
                    ") must match size of input alphabet of second fsm (was " +
                    std::to_string(g_overs[constraint.acceptorFsmID]->in_alphabet) + ")");
        }

        bool generator_is_deterministic = constraint.generator_is_deterministic;

        FSMGeneratorAcceptorDetector* d = nullptr;
        auto key = std::tuple<int, int, int, int, bool>(constraint.generatorFsmID, constraint.acceptorFsmID,
                                                        constraint.generatorSourceState, constraint.acceptorSourceState,
                                                        generator_is_deterministic);
        if(gen_accept_map.has(key)){
            d = gen_accept_map[key];
        }else{

            d = new FSMGeneratorAcceptorDetector(detectors.size(), this, *g_unders[constraint.generatorFsmID],
                                                 *g_overs[constraint.generatorFsmID],
                                                 *g_unders[constraint.acceptorFsmID],
                                                 *g_overs[constraint.acceptorFsmID], constraint.generatorSourceState,
                                                 constraint.acceptorSourceState,
                                                 drand(rnd_seed));
            detectors.push(d);
            if(generator_is_deterministic){
                d->setGeneratorDeterministic(true);
            }
            gen_accept_map.insert(key, d);
        }

        d->addAcceptLit(constraint.generatorAcceptState, constraint.acceptorAcceptState, constraint.outer_var);
        gen_accept_lit_map.growTo(constraint.outer_var + 1, nullptr);
        gen_accept_lit_map[constraint.outer_var] = d;

        freezeFSM(constraint.generatorFsmID);
        freezeFSM(constraint.acceptorFsmID);
    }

public:

    void addComposeAcceptSuffixFSM(Var composeAcceptVar, int suffix_fsmID){
        assert(composeAcceptVar < gen_accept_lit_map.size());
        assert(gen_accept_lit_map[composeAcceptVar]);
        gen_accept_lit_map[composeAcceptVar]->setSuffixGenerator(this->g_overs[suffix_fsmID]);
    }

    void addComposeAcceptSuffixLit(Var composeAcceptVar, int startSuffixState, int acceptSuffixState, Var suffixVar){
        assert(composeAcceptVar < gen_accept_lit_map.size());
        assert(gen_accept_lit_map[composeAcceptVar]);
        Var v = this->S->getTheoryVar(composeAcceptVar, this);
        gen_accept_lit_map[composeAcceptVar]->addSuffixLit(mkLit(v), startSuffixState, acceptSuffixState, suffixVar);
    }

};

};

#endif
