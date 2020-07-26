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
#ifndef FSM_GENERATORACCEPTDETECTOR_H_
#define FSM_GENERATORACCEPTDETECTOR_H_

#include "monosat/utils/System.h"

#include "monosat/dgl/DynamicGraph.h"

#include "monosat/fsm/DynamicFSM.h"

#include "monosat/core/SolverTypes.h"
#include "monosat/mtl/Map.h"
#include "monosat/mtl/Bitset.h"

#include "monosat/utils/System.h"
#include "monosat/fsm/FSMDetector.h"
#include "monosat/fsm/alg/NFALinearGeneratorAcceptor.h"
#include "monosat/graph/GraphTheory.h"

using namespace dgl;
namespace Monosat {

class FSMTheorySolver;

class FSMGeneratorAcceptorDetector : public FSMDetector {
public:
    FSMTheorySolver* outer;
    DynamicFSM& g_under;
    DynamicFSM& g_over;
    DynamicFSM& acceptor_under;
    DynamicFSM& acceptor_over;

    DynamicFSM* suffix_fsm = nullptr;

    int first_destination = -1;
    int gen_source;
    int accept_source;
    double rnd_seed;

    struct AcceptStatus {
        FSMGeneratorAcceptorDetector& detector;
        bool polarity;

        void accepts(int string, int state, int edgeID, int label, bool accepts);

        AcceptStatus(FSMGeneratorAcceptorDetector& _outer, bool _polarity) :
                detector(_outer), polarity(_polarity){
        }
    };

    AcceptStatus* underReachStatus = nullptr;
    AcceptStatus* overReachStatus = nullptr;

    NFALinearGeneratorAcceptor<AcceptStatus>* underapprox_detector = nullptr;
    NFALinearGeneratorAcceptor<AcceptStatus>* overapprox_detector = nullptr;
    NFALinearGeneratorAcceptor<AcceptStatus>* inverted_overapprox_detector = nullptr;


    CRef underprop_marker;
    CRef overprop_marker;
    CRef forcededge_marker;
    CRef forced_positive_nondet_generator_transition_marker;
    CRef deterministic_forcededge_marker;
    CRef chokepoint_acceptor_transition_marker;
    struct Change {
        Lit l;
        int u;
        int str;
    };
    vec<Change> changed;

    Var first_var = var_Undef;

    struct SuffixLit {
        Lit l;
        Lit accept_l;
        int suffix_from;
        int suffix_to;

        SuffixLit(Lit l, Lit accept_l, int _suffix_from, int _suffix_to) : l(l), accept_l(accept_l),
                                                                           suffix_from(_suffix_from),
                                                                           suffix_to(_suffix_to){

        }
    };

    struct AcceptLit {
        Lit l = lit_Undef;
        int gen_to = -1;
        int accept_to = -1;
        vec<SuffixLit>* suffixLits = nullptr;

        AcceptLit(){

        }

        AcceptLit(Lit l, int gen_to, int accept_to) : l(l), gen_to(gen_to), accept_to(accept_to){

        }

        AcceptLit(Lit l, int gen_to, int accept_to, vec<SuffixLit>* suffixLits) : l(l), gen_to(gen_to),
                                                                                  accept_to(accept_to),
                                                                                  suffixLits(suffixLits){

        }
    };

    vec<int> accept_lit_map;
    vec<AcceptLit> all_accept_lits;
    vec<Lit> all_lits;

    bool relp = false;
    bool reln = false;
    bool accp = false;
    bool accn = false;
    vec<SuffixLit> suffix_lit_map;
    vec<SuffixLit> all_suffix_lits;
    //stats

    int64_t stats_full_updates = 0;
    int64_t stats_fast_updates = 0;
    int64_t stats_fast_failed_updates = 0;
    int64_t stats_skip_deletes = 0;
    int64_t stats_skipped_updates = 0;
    int64_t stats_num_skipable_deletions = 0;
    int64_t stats_learnt_components = 0;
    int64_t stats_learnt_components_sz = 0;
    int64_t stats_forced_edges = 0;
    double mod_percentage = 0.2;
    int64_t stats_pure_skipped = 0;
    int64_t stats_shrink_removed = 0;
    double stats_full_update_time = 0;
    double stats_fast_update_time = 0;

    int64_t stats_chokepoint_edge_propagations = 0;
    int64_t stats_forced_edge_propagations = 0;
    int64_t stats_nondet_gen_edge_propagations = 0;
    int64_t stats_det_gen_edge_propagations = 0;

    Map<Var, Lit> forcedVars;
    vec<Var> lit_backward_map;
    vec<bool> seen_chars;

    GraphTheorySolver<int64_t>* graph = nullptr;//for if we reduce the nfa to a graph
    vec<vec<int>> nodes;

    void printStats() override{

        if(graph){
            graph->printStats(1);
        }else{
            FSMDetector::printStats();
            if(opt_detect_pure_theory_lits)
                printf("\tPropagations skipped by pure literal detection: %" PRId64 "\n", stats_pure_skipped);
            if(opt_shrink_theory_conflicts){
                printf("\t%" PRId64 " lits removed by shrinking conflicts\n", stats_shrink_removed);
            }
            if(opt_learn_unreachable_component){
                printf("\t%" PRId64 " components learned, average component size: %f\n", stats_learnt_components,
                       stats_learnt_components_sz / (float) stats_learnt_components);
            }
            if(opt_fsm_forced_edge_prop){
                printf("Forced edge assignments: %" PRId64 "\n", stats_forced_edge_propagations);
            }
            if(opt_fsm_chokepoint_prop){
                printf("Chokepoint edge assignments: %" PRId64 "\n", stats_chokepoint_edge_propagations);
            }
            if(opt_fsm_edge_prop){
                printf("Generator edge assignments (deterministic,nondeterministic): %" PRId64 ",%" PRId64 "\n",
                       stats_det_gen_edge_propagations, stats_nondet_gen_edge_propagations);
            }

        }
    }

    inline void setForcedVar(Var edgeVar, Lit forcedBy){
        if(forcedVars.has(edgeVar)){
            forcedVars.remove(edgeVar);
        }

        forcedVars.insert(edgeVar, forcedBy);
    }

    inline Lit getForcedVar(Var edgeVar){
        return forcedVars[edgeVar];
    }

    inline int indexOf(Var v) const{
        int index = v - first_var;
        assert(index < accept_lit_map.size());
        return index;
    }

    int getGeneratorFinal(Var reachVar){
        assert(reachVar >= first_var);
        int index = indexOf(reachVar);
        assert(all_accept_lits[accept_lit_map[index]].gen_to >= 0);
        return all_accept_lits[accept_lit_map[index]].gen_to;
    }

    int getAcceptorFinal(Var reachVar){
        assert(reachVar >= first_var);
        int index = indexOf(reachVar);

        assert(all_accept_lits[accept_lit_map[index]].accept_to >= 0);
        return all_accept_lits[accept_lit_map[index]].accept_to;
    }

    void setSuffixGenerator(DynamicFSM* suffix_fsm){
        this->suffix_fsm = suffix_fsm;
    }

    Lit decide(int level) override;

    bool propagate(vec<Lit>& conflict) override;

    void buildAcceptReason(int genFinal, int acceptFinal, vec<Lit>& conflict);

    void buildNonAcceptReason(int genFinal, int acceptFinal, vec<Lit>& conflict);

    void buildNonSuffixAcceptReason(int genFinal, int acceptFinal, vec<int>& acceptStartStates, SuffixLit& suffixLit,
                                    vec<Lit>& conflict);

    void buildForcedEdgeReason(int genFinal, int acceptFinal, int forcedEdge, int forcedLabel, vec<Lit>& conflict);

    void
    buildDeterministicForcedEdgeReason(int genFinal, int acceptFinal, int forcedGeneratorEdge, int forcedGeneratorLabel,
                                       vec<Lit>& conflict);

    void
    buildAcceptorChokepointEdgeReason(int genFinal, int acceptFinal, int forcedAcceptorEdge, int forcedAcceptorLabel,
                                      vec<Lit>& conflict);

    void
    buildForcedNondetEdgeReason(int genFinal, int acceptFinal, int forcedEdge, int forcedLabel, vec<Lit>& conflict);

    void preprocess() override;

    void buildReason(Lit p, vec<Lit>& reason, CRef marker) override;

    bool checkSatisfied() override;

    void printSolution(std::ostream& write_to) override;

    bool checkNegatedPolarity() override;

    void addAcceptLit(int state, int strID, Var reach_var);

    void addSuffixLit(Lit acceptorLit, int suffixStartState, int suffixAcceptState, Var outer_reach_var);

    Var getDetectorVar(int gen_final, int accept_final){
        return lit_backward_map[gen_final + accept_final * g_over.states()];
    }

    FSMGeneratorAcceptorDetector(int _detectorID, FSMTheorySolver* _outer, DynamicFSM& g_under, DynamicFSM& g_over,
                                 DynamicFSM& acceptor_under, DynamicFSM& acceptor_over,
                                 int gen_source, int acceptor_source, double seed = 1);

    const char* getName() override{
        return "NFA Generator Acceptor Detector";
    }

    void setGeneratorDeterministic(bool is_deterministic){
        generator_is_deterministic = is_deterministic;
    }

private:

    struct Transition {
        int edgeID;
        int in;
        int out;
    };

    void constructAllPaths(bool add_graph_symbols = false);

    void stepGeneratorForward(vec<Transition>& store, vec<bool>& store_seen, int& cur_gen_state);

    vec<int> next;
    vec<int> cur;

    vec<bool> next_seen;
    vec<bool> cur_seen;
    bool generator_is_deterministic = false;
    vec<int> gen_cur;
    vec<int> gen_next;
    vec<bool> gen_next_seen;
    vec<bool> gen_cur_seen;
    vec<int> chars;
    vec<bool> tmp_seen_chars;
    vec<vec<Bitset>> prefixTables;
    vec<vec<ForcedTransition>> forced_edges;
    vec<ChokepointTransition> chokepoint_edges;
    vec<int> pre_accepting_states;//states of the fsm that lead to the accepting state in one transition, in the current generator
    vec<NFATransition> tmp_path;
    int last_prefix_update = -1;
    vec<int> dfs_q;

    bool isAttractor(int acceptorState);

    //note: forced edge and forced_label refer to the generator fsm, while ignoredEdge and  ignoredLabel refer to the acceptor fsm
    bool find_gen_path(int gen_final, int accept_final, int forcedEdge, int forcedLabel, vec<NFATransition>& path,
                       bool invertAcceptance = false, bool all_paths = false, bool under_approx = false,
                       int ignoredEdge = -1, int ignoredLabel = -1);

    bool stepGenerator(int final, int forcedEdge, int forcedLabel, vec<int>& store, vec<bool>& store_seen,
                       int& cur_gen_state, vec<NFATransition>* path = nullptr);

    bool buildSuffixCut(int gen_final, int accept_final, vec<Lit>& cut, bool accepting_state_is_attractor,
                        bool invertAcceptance, int ignoreGeneratorEdge = -1, int ignoreGeneratorLabel = -1,
                        int forcedAcceptorEdge = -1, int forcedAcceptorLabel = -1);

    bool
    stepGeneratorBackward(int final, vec<Bitset>& prefixTable, vec<Lit>& cut, vec<int>& store, vec<bool>& store_seen,
                          vec<NFATransition>* path = nullptr, int ignoreGeneratorEdge = -1,
                          int ignoreGeneratorLabel = -1);

    void updatePrefixTable(int gen_final, int accept_final);

    //True if the acceptor fsm, starting from any of 'acceptor_start_states', and starting from ALL suffix_start_states (!), accepts g_suffix
    bool acceptsSuffixes(DynamicFSM& acceptor, int acceptor_accept_state, vec<SuffixLit>& suffix_lits,
                         vec<int>& acceptor_start_states, int& failedSuffixLit);

    bool acceptsSuffix(DynamicFSM& g_suffix, DynamicFSM& acceptor, int acceptor_accept_state, int suffix_start_state,
                       int suffix_accept_state, vec<int>& acceptor_start_states);

    vec<Bitset>& getPrefixTable(int gen_final, int accept_final){
        updatePrefixTable(gen_final, accept_final);
        Var v = getDetectorVar(gen_final, accept_final);
        assert(v != var_Undef);
        int index = v - first_var;
        return prefixTables[index];
    }

    bool isAttractor(DynamicFSM& accept, int acceptorState){
        if(acceptorState < 0){
            return true;
        }else{
            //should really fix this to work correctly for epsilon transitions between multiple acceptor states...
            for(int i = 0; i < accept.nIncident(acceptorState); i++){
                int edgeID = accept.incident(acceptorState, i).id;
                int to = accept.incident(acceptorState, i).node;
                if(to == acceptorState){
                    for(int c = 1; c < accept.inAlphabet(); c++){
                        if(!accept.transitionEnabled(edgeID, c, -1)){
                            return false;
                        }
                    }
                    return true;
                }
            }
        }
        return false;
    }

};
};
#endif
