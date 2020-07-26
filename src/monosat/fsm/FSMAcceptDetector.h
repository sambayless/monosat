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
#ifndef FSM_ACCEPTDETECTOR_H_
#define FSM_ACCEPTDETECTOR_H_

#include "monosat/utils/System.h"

#include "monosat/dgl/DynamicGraph.h"

#include "monosat/fsm/DynamicFSM.h"

#include "monosat/core/SolverTypes.h"
#include "monosat/mtl/Map.h"
#include "monosat/mtl/Heap.h"

#include "monosat/utils/System.h"
#include "monosat/fsm/FSMDetector.h"
#include "monosat/fsm/alg/NFAAccept.h"
#include "monosat/fsm/alg/NFAAcceptor.h"

using namespace dgl;
namespace Monosat {

class FSMTheorySolver;

class FSMAcceptDetector : public FSMDetector {
public:
    FSMTheorySolver* outer;
    DynamicFSM& g_under;
    DynamicFSM& g_over;


    int source;
    int first_destination = -1;
    vec<vec<int>>& strings;
    double rnd_seed;

    struct AcceptStatus {
        FSMAcceptDetector& detector;
        bool polarity;

        void accepts(int string, int state, int edgeID, int label, bool accepts);

        AcceptStatus(FSMAcceptDetector& _outer, bool _polarity) :
                detector(_outer), polarity(_polarity){
        }
    };

    AcceptStatus* underReachStatus = nullptr;
    AcceptStatus* overReachStatus = nullptr;

    NFAAcceptor* underapprox_detector;
    NFAAcceptor* overapprox_detector;


    CRef underprop_marker;
    CRef overprop_marker;

    struct Change {
        Lit l;
        int u;
        int str;
        bool polarity;
    };
    vec<Change> changed;

    vec<vec<Lit>> accept_lits;
    Var first_var = var_Undef;
    struct AcceptLit {
        int str;
        int to;
    };
    vec<AcceptLit> accept_lit_map;
    vec<Lit> all_lits;
    vec<Lit> to_decide;
    vec<NFATransition> decision_path;
    int last_decision_status = -1;
    //stats

    int64_t stats_full_updates = 0;
    int64_t stats_fast_updates = 0;
    int64_t stats_fast_failed_updates = 0;
    int64_t stats_skip_deletes = 0;
    int64_t stats_skipped_updates = 0;
    int64_t stats_num_skipable_deletions = 0;
    int64_t stats_learnt_components = 0;
    int64_t stats_learnt_components_sz = 0;
    double mod_percentage = 0.2;
    int64_t stats_pure_skipped = 0;
    int64_t stats_shrink_removed = 0;
    double stats_full_update_time = 0;
    double stats_fast_update_time = 0;
    int64_t stats_symmetry_conflicts = 0;
    vec<int> accepting_state;//maintains a list of all accepting states, which are not considered during symmetry breaking.


    void printStats() override{
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
        if(opt_fsm_symmetry_breaking){
            printf("Symmetry breaking conflicts: %" PRId64 "\n", stats_symmetry_conflicts);
        }
    }

    void assign(Lit l) override{
        FSMDetector::assign(l);
        int index = indexOf(var(l));
        if(index >= 0 && index < accept_lit_map.size() && accept_lit_map[index].to != -1){
            int node = accept_lit_map[index].to;
            int str = accept_lit_map[index].str;

            changed.push({l, node, str});
            if(opt_theory_internal_vsids_fsm && !order_heap.inHeap(var(l))){
                order_heap.insert(var(l));
            }
        }
    }

    void unassign(Lit l) override{
        FSMDetector::unassign(l);
        int index = indexOf(var(l));
        if(index >= 0 && index < accept_lit_map.size() && accept_lit_map[index].to != -1){
            int node = accept_lit_map[index].to;
            int str = accept_lit_map[index].str;
            changed.push({l, node, str});
        }
    }

    inline int indexOf(Var v) const{
        int index = v - first_var;
        assert(index < accept_lit_map.size());
        return index;
    }

    int getState(Var reachVar){
        assert(reachVar >= first_var);
        int index = indexOf(reachVar);

        assert(accept_lit_map[index].to >= 0);
        return accept_lit_map[index].to;
    }

    int getString(Var reachVar){
        assert(reachVar >= first_var);
        int index = indexOf(reachVar);

        assert(accept_lit_map[index].to >= 0);
        return accept_lit_map[index].str;
    }

    void backtrack(int level) override{
        to_decide.clear();
        last_decision_status = -1;

    }

    bool propagate(vec<Lit>& conflict) override;

    Lit decide(int level) override;

    void buildAcceptReason(int node, int str, vec<Lit>& conflict);

    void buildNonAcceptReason(int node, int str, vec<Lit>& conflict);

    void preprocess() override{
        accepting_state.growTo(g_over.states());
    }

    void buildReason(Lit p, vec<Lit>& reason, CRef marker) override;

    bool checkSatisfied() override;

    void printSolution(std::ostream& write_to) override;

    void addAcceptLit(int state, int strID, Var reach_var);

    bool checkSymmetryConstraints(vec<Lit>& conflict);

    bool checkSymmetryConstraintsPopCount(vec<Lit>& conflict);

    Bitset _bitvec1;
    Bitset _bitvec2;

    bool checkSymmetryConstraintsBitVec(vec<Lit>& conflict);

    void learnClauseSymmetryConflict(vec<Lit>& conflict, int a, int b);

    FSMAcceptDetector(int _detectorID, FSMTheorySolver* _outer, DynamicFSM& g_under, DynamicFSM& g_over,
                      int _source, vec<vec<int>>& strs, double seed = 1);

    ~FSMAcceptDetector() override{

    }

    const char* getName() override{
        return "NFA Accepts Detector";
    }

    double var_inc = 1;
    double var_decay = 1;
    vec<double> activity;

    struct AcceptOrderLt {
        FSMAcceptDetector& outer;
        const vec<double>& activity;

        bool operator()(Var x, Var y) const{
            int indexX = outer.indexOf(x);
            int indexY = outer.indexOf(y);
            return activity[indexX] > activity[indexY];

        }

        AcceptOrderLt(FSMAcceptDetector& outer, const vec<double>& act) : outer(outer),
                                                                          activity(act){
        }
    };

    //Local vsids implementation for edges...
    Heap<int, AcceptOrderLt> order_heap;


    inline void insertAcceptOrder(Lit acceptLit){
        if(opt_theory_internal_vsids_fsm){
            int index = indexOf(var(acceptLit));
            if(!order_heap.inHeap(var(acceptLit)) && activity[index] > 0)
                order_heap.insert(var(acceptLit));
        }
    }

    inline void bumpConflict(vec<Lit>& conflict){

        if(opt_theory_internal_vsids_fsm){
            for(Lit l:conflict){
                int index = indexOf(var(l));
                if(index >= 0 && index < accept_lit_map.size() && accept_lit_map[index].to != -1){

                    edgeBumpActivity(var(l));
                    edgeDecayActivity();
                    if(!order_heap.inHeap(var(l))){
                        order_heap.insert(var(l));
                    }
                }
            }
        }

    }

    inline void bumpAcceptLit(Lit conflict){

        if(opt_theory_internal_vsids_fsm){
            edgeBumpActivity(var(conflict));
            edgeDecayActivity();
        }

    }

    inline void edgeDecayActivity(){
        var_inc *= (1 / var_decay);
    }

    inline void edgeBumpActivity(Var v){
        edgeBumpActivity(v, var_inc);
    }

    inline void edgeBumpActivity(Var v, double inc){
        if((activity[indexOf(v)] += inc) > 1e100){
            // Rescale:
            for(int i = 0; i < activity.size(); i++)
                activity[i] *= 1e-100;
            var_inc *= 1e-100;
        }

        // Update order_heap with respect to new activity:
        if(order_heap.inHeap(v))
            order_heap.decrease(v);
    }

};
};
#endif
