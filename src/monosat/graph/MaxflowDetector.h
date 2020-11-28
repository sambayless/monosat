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
#ifndef MAXFLOWDETECTOR_H_
#define MAXFLOWDETECTOR_H_

#include "monosat/utils/System.h"
#include "monosat/dgl/KohliTorr.h"
#include "GraphTheoryTypes.h"
#include "monosat/dgl/DynamicGraph.h"
#include "monosat/dgl/MaxFlow.h"
#include "monosat/dgl/EdmondsKarp.h"
#include "monosat/core/SolverTypes.h"
#include "monosat/mtl/Map.h"
#include "monosat/mtl/Deque.h"
#include "monosat/utils/System.h"
#include "Detector.h"
#include "monosat/bv/BVTheorySolver.h"
#include "monosat/dgl/AcyclicFlow.h"
#include "monosat/core/Config.h"
#include <vector>

using namespace dgl;
namespace Monosat {
template<typename Weight>
class GraphTheorySolver;

template<typename Weight = int, typename Graph = DynamicGraph<Weight>>
class MaxflowDetector : public Detector, public EdgeDecider<Weight>, public DynamicGraphAlgorithm {
public:
    GraphTheorySolver<Weight>* outer;
    std::vector<Weight> capacities;

    Graph& g_under;
    Graph& g_over;

    int source;
    int target;
    double rnd_seed;
    CRef underprop_marker;
    CRef overprop_marker;

    BVTheorySolver<Weight>* bvTheory = nullptr;

    class FlowOp : public GraphTheorySolver<Weight>::GraphTheoryOp {
        MaxflowDetector* outer;
        int bvID;
    public:
        FlowOp(BVTheorySolver<Weight>& theory, MaxflowDetector* outer, int bvID)
                : GraphTheorySolver<Weight>::GraphTheoryOp(theory, outer->outer), outer(outer), bvID(bvID){

        }

        int getBV() override{
            return bvID;
        }

        bool propagate(bool& changed, vec<Lit>& conflict) override{
            return true;
        }

        void updateApprox(Var ignore_bv, Weight& under_new, Weight& over_new,
                          typename BVTheorySolver<Weight>::Cause& under_cause_new,
                          typename BVTheorySolver<Weight>::Cause& over_cause_new) override{
            ;
        }

        void analyzeReason(bool compareOver, Comparison op, Weight to, vec<Lit>& conflict) override;

        bool checkSolved() override{
            return true;
        }
    };


    MaxFlow<Weight>* underapprox_detector = nullptr;
    MaxFlow<Weight>* overapprox_detector = nullptr;
    MaxFlow<Weight>* underapprox_conflict_detector = nullptr;
    MaxFlow<Weight>* overapprox_conflict_detector = nullptr;
    AcyclicFlow<Weight>* acyclic_flow = nullptr;
    std::vector<Weight> refined_flow;
    std::vector<Weight> refined_flow_model;

    int last_decision_status = -1;
    int last_decision_q_pos = 0;
    int alg_id = -1;

    int64_t stats_decision_calculations = 0;
    double stats_total_prop_time = 0;
    double stats_flow_calc_time = 0;
    double stats_flow_recalc_time = 0;
    double stats_redecide_time = 0;
    int64_t stats_heuristic_recomputes = 0;
    Lit last_decision_lit = lit_Undef;

    //vec<Lit> decisions;
    vec<bool> is_potential_decision;//not used for bitvector edges
    vec<bool> in_decision_q;
    vec<int> potential_decisions;

    struct DecisionS {
        //int path_decision:1;
        int edgeID:32;

        DecisionS() : edgeID(-1){}

        DecisionS(const int edgeID) : edgeID(edgeID){}

        DecisionS(int edgeID, bool pathDecision) : edgeID(edgeID){}

        operator int() const{return edgeID;}
    };

    Deque<DecisionS> potential_decisions_q;
    vec<Lit> to_decide;
    std::vector<int> q;

    Graph learn_graph;
    vec<int> back_edges;
    int learngraph_history_qhead = 0;
    int learngraph_history_clears = -1;
    bool overIsEdgeSet = false;
    MaxFlow<Weight>* learn_cut = nullptr;
    //int current_decision_edge=-1;
    //vec<Lit>  reach_lits;
    Var first_reach_var;
    vec<int> reach_lit_map;
    vec<int> force_reason;

    vec<int> priority_decisions;

    struct DistLit {
        Lit l = lit_Undef;
        Weight max_flow = (Weight) -1;
        BitVector<Weight> bv;
        bool inclusive = false;//If inclusive true, l is true iff the maximum flow is >= max_flow; else, l is true iff the maximum flow is > max_flow.
        FlowOp* op = nullptr;
    };
    vec<DistLit> flow_lits;

    int n_satisfied_lits = 0;

    std::vector<MaxFlowEdge> cut;

    vec<MaxFlowEdge> tmp_cut;
    vec<int> visit;
    vec<bool> seen;
    vec<bool> seen_path;

    struct FlowListener {
        virtual void edgeFlowChange(int edgeID, const Weight& flow) = 0;
    } * flowListener = nullptr;

    void backtrack(int level) override{
        to_decide.clear();
        last_decision_status = -1;
    }

    void collectChangedEdges();

    void collectDisabledEdges();

    void updateHistory() override{
        collectDisabledEdges();
    }

    bool propagate(vec<Lit>& conflict) override{
        Lit ignore = lit_Undef;
        return propagate(conflict, false, ignore);
    }


    bool propagate(vec<Lit>& conflict, bool backtrackOnly, Lit& conflictLit) override;

    void analyzeMaxFlowLEQ(Weight flow, vec<Lit>& conflict, bool force_maxflow = false);

    void analyzeMaxFlowGEQ(Weight flow, vec<Lit>& conflict);

    void buildMaxFlowTooHighReason(Weight flow, vec<Lit>& conflict);

    Lit findFirstReasonTooHigh(Weight flow);

    Lit findFirstReasonTooLow(Weight flow);

    void buildMaxFlowTooLowReason(Weight flow, vec<Lit>& conflict, bool force_maxflow = false);

    void buildForcedEdgeReason(int reach_node, int forced_edge_id, vec<Lit>& conflict);

    void buildReason(Lit p, vec<Lit>& reason, CRef marker) override;

    bool checkSatisfied() override;

    bool decideEdgeWeight(int edgeID, Weight& store, DetectorComparison& op) override;

    void undecideEdgeWeight(int edgeID) override;

    void undecide(Lit l) override;

    void debug_decidable(Var v) override;

    void assignBV(int bvID) override;

    void unassignBV(int bvID) override;

    void assign(Lit l) override{
        Detector::assign(l);
        if(default_heuristic){
            outer->activateHeuristic(default_heuristic);
        }
    }

    void setSatisfied(Lit l, bool isSatisfied) override;

    bool detectorIsSatisfied() override;

    Lit decide(CRef& decision_reason) override;

    bool supportsEdgeDecisions() override{
        return true;
    }

    void attachFlowListener(FlowListener* l){
        assert(!flowListener);
        flowListener = l;
    }

    void suggestDecision(Lit l) override;

    void printStats() override{
        Detector::printStats();
        printf("\tTotal Detector Propagation Time: %fs\n", stats_total_prop_time);
        if(mincutalg == MinCutAlg::ALG_KOHLI_TORR){
            KohliTorr<Weight>* kt = (KohliTorr<Weight>*) overapprox_detector;
            printf(
                    "\tInit Time %f, Decision flow calculations: %" PRId64 ", (redecide: %f s) flow_calc %f s, flow_discovery %f s, (%" PRId64 ") (maxflow %f,flow assignment %f),  inits: %" PRId64 ",re-inits %" PRId64 "\n",
                    kt->stats_init_time, stats_decision_calculations, stats_redecide_time, stats_flow_calc_time,
                    stats_flow_recalc_time,
                    kt->stats_flow_calcs, kt->stats_flow_time, kt->stats_calc_time, kt->stats_inits, kt->stats_reinits);
        }else
            printf("\tDecision flow calculations: %" PRId64 "\n", stats_decision_calculations);
        if(n_stats_priority_decisions > 0){
            printf("\tPriority decisions: %" PRId64 "\n", n_stats_priority_decisions);
        }
        if(opt_theory_internal_vsids){
            printf("\tVsids decisions: %" PRId64 "\n", n_stats_vsids_decisions);
        }

    }

    Weight computeUnderApprox(Weight& computed_under_weight){
        if(computed_under_weight > -1){
            return computed_under_weight;
        }
        if(underapprox_detector && (!opt_detect_pure_theory_lits || unassigned_positives > 0)){
            double startdreachtime = rtime(2);
            stats_under_updates++;

            computed_under_weight = underapprox_detector->maxFlow();
            assert(computed_under_weight == underapprox_conflict_detector->maxFlow());
            double reachUpdateElapsed = rtime(2) - startdreachtime;
            stats_under_update_time += reachUpdateElapsed;
            return computed_under_weight;
        }
        return 0;
    }

    Weight computeOverApprox(Weight& computed_over_weight){
        if(computed_over_weight > -1){
            return computed_over_weight;
        }
        if(overapprox_detector && (!opt_detect_pure_theory_lits || unassigned_negatives > 0)){
            double startunreachtime = rtime(2);
            stats_over_updates++;
            computed_over_weight = overapprox_detector->maxFlow();
            assert(computed_over_weight == overapprox_conflict_detector->maxFlow());
            double unreachUpdateElapsed = rtime(2) - startunreachtime;
            stats_over_update_time += unreachUpdateElapsed;
            return computed_over_weight;
        }
        return 0;
    }

    //Lit decideByPath(int level);
    void dbg_decisions();

    void printSolution(std::ostream& write_to) override;

    Weight getModel_Maxflow(){
        underapprox_detector->update();
        return underapprox_detector->maxFlow();
    }

    Weight getModel_EdgeFlow(int edgeID){

        return underapprox_detector->getEdgeFlow(edgeID);
    }

    Weight getModel_AcyclicEdgeFlow(int edgeID){
        if(!acyclic_flow){
            acyclic_flow = new AcyclicFlow<Weight>(g_under);
        }
        if(refined_flow_model.size() == 0){
            underapprox_detector->update();
            refined_flow_model.resize(g_under.edges());
            for(int i = 0; i < g_under.edges(); i++){
                if(g_under.hasEdge(i) && g_under.edgeEnabled(i)){
                    refined_flow_model[i] = underapprox_detector->getEdgeFlow(i);
                }else{
                    refined_flow_model[i] = 0;
                }
            }
            acyclic_flow->getAcyclicFlow(source, target, refined_flow_model);
        }
        return refined_flow_model[edgeID];

    }

    void buildModel() override{
        underapprox_detector->update();
        refined_flow_model.clear();
    }

    void setFlowBV(const BitVector<Weight>& bv);

    Lit addFlowLit(Weight max_flow, Var reach_var, bool inclusive);

    Lit addMaxFlowGEQ_BV(const BitVector<Weight>& bv, Var v, bool inclusive);

    MaxflowDetector(int _detectorID, GraphTheorySolver<Weight>* _outer,
                    Graph& _g, Graph& _antig, int _source, int _target, double seed = 1,
                    bool overIsEdgeSet = false); //:Detector(_detectorID),outer(_outer),within(-1),source(_source),rnd_seed(seed),positive_reach_detector(NULL),negative_reach_detector(NULL),positive_path_detector(NULL),positiveReachStatus(NULL),negativeReachStatus(NULL){}
    ~MaxflowDetector() override{
        if(underapprox_conflict_detector && underapprox_conflict_detector != underapprox_detector)
            delete underapprox_conflict_detector;
        if(overapprox_conflict_detector && overapprox_conflict_detector != overapprox_detector)
            delete overapprox_conflict_detector;

        if(underapprox_detector)
            delete underapprox_detector;
        if(overapprox_detector)
            delete overapprox_detector;

        if(learn_cut)
            delete learn_cut;
        if(acyclic_flow)
            delete acyclic_flow;
    }

    std::string getName() override{
        return "Max-flow Detector";
    }

    bool preprocessed = false;

    void preprocess() override{
        preprocessed = true;
        activity.growTo(g_under.edges());
        seen.growTo(outer->nNodes());
        is_potential_decision.growTo(g_under.edges(), false);

        in_decision_q.growTo(g_under.edges(), false);
        int max_decision_priority = -1;
        for(int i = 0; i < flow_lits.size(); i++){
            Lit l = flow_lits[i].l;
            int priority = outer->getSolver()->getDecisionPriority(var(outer->toSolver(l)));
            if(priority > max_decision_priority){
                max_decision_priority = priority;
            }
        }
        default_heuristic->setPriority(max_decision_priority);
    }

private:

    void buildLearnGraph(){
        if(learn_graph.nodes() < g_under.nodes()){
            while(learn_graph.nodes() < g_under.nodes())
                learn_graph.addNode();
            back_edges.growTo(g_under.edges(), -1);

            for(auto& e : g_under.getEdges()){
                if(back_edges[e.id] == -1){
                    learn_graph.addEdge(e.from, e.to, -1, 0);
                }
            }

            for(auto& e : g_under.getEdges()){
                if(back_edges[e.id] == -1){
                    back_edges[e.id] = learn_graph.addEdge(e.to, e.from, -1, 0);
                }
            }

            learn_graph.invalidate();

            learn_graph.clearHistory(true);
        }
    }

    double var_inc = 1;
    double var_decay = 1;
    vec<double> activity;

    struct EdgeOrderLt {
        const vec<double>& activity;

        bool operator()(Var x, Var y) const{
            return activity[x] > activity[y];
        }

        EdgeOrderLt(const vec<double>& act) :
                activity(act){
        }
    };

    //Local vsids implementation for edges...
    Heap<int, EdgeOrderLt> order_heap;

    inline void insertEdgeOrder(int edgeID){
        if(opt_theory_internal_vsids){
            if(!order_heap.inHeap(edgeID) && activity[edgeID] > 0)
                order_heap.insert(edgeID);
        }
    }

    inline void bumpConflictEdges(vec<Lit>& conflict){
        if(opt_theory_prioritize_conflicts){
            for(Lit l:conflict){
                if(outer->isEdgeVar(var(l))){
                    int edgeID = outer->getEdgeID(var(l));
                    priority_decisions.push(edgeID);
                }
            }
        }
        if(opt_theory_internal_vsids){
            for(Lit l:conflict){
                if(outer->isEdgeVar(var(l))){
                    int edgeID = outer->getEdgeID(var(l));
                    edgeBumpActivity(edgeID);
                }
            }
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
        if((activity[v] += inc) > 1e100){
            // Rescale:
            for(int i = 0; i < g_over.edges(); i++)
                activity[i] *= 1e-100;
            var_inc *= 1e-100;
        }

        // Update order_heap with respect to new activity:
        if(order_heap.inHeap(v))
            order_heap.decrease(v);
    }


};

};

#endif /* DistanceDetector_H_ */
