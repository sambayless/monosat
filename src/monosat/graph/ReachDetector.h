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
#ifndef REACHDETECTOR_H_
#define REACHDETECTOR_H_

#include "monosat/utils/System.h"
#include "GraphTheoryTypes.h"
#include "monosat/dgl/DynamicGraph.h"
#include "monosat/dgl/Reach.h"
#include "monosat/dgl/Dijkstra.h"
#include "monosat/dgl/BFS.h"
#include "monosat/dgl/DFS.h"

#include "monosat/core/SolverTypes.h"
#include "monosat/mtl/Map.h"
#include "monosat/dgl/MaxFlow.h"

#include "monosat/dgl/EdmondsKarp.h"
#include "monosat/dgl/EdmondsKarpAdj.h"
#include "monosat/dgl/Chokepoint.h"
#include "WeightedDijkstra.h"

#include "monosat/utils/System.h"
#include "Detector.h"
#include "monosat/graph/MaxflowDetector.h"

using namespace dgl;
namespace Monosat {
template<typename Weight>
class GraphTheorySolver;

template<typename Weight, typename Graph = DynamicGraph<Weight>>
class ReachDetector : public Detector {
public:
    GraphTheorySolver<Weight>* outer;
    Graph& g_under;
    Graph& g_over;
    Graph& cutGraph;
    DynamicGraph<Weight> localCutGraph;//todo: this is not really used any more, remove
    int within;
    int source;
    double rnd_seed;
    int constraintsBuiltOver = -1;
    int constraintsBuiltUnder = -1;

    CRef underprop_marker;
    CRef overprop_marker;
    CRef forced_edge_marker;

    Reach* underapprox_detector = nullptr;
    Reach* overapprox_reach_detector = nullptr;
    Reach* underapprox_path_detector = nullptr;
    Reach* overapprox_path_detector = nullptr;
    Reach* cutgraph_detector = nullptr;
    Reach* underapprox_fast_detector = nullptr;
    Distance<int>* negative_distance_detector = nullptr;
    vec<bool> original_reach_lits;
    vec<Lit> reach_lits;
    vec<Lit> cnf_reach_lits;
    Var first_reach_var;

    vec<int> reach_lit_map;
    vec<int> force_reason;
    vec<Heuristic*> all_reach_heuristics;
    vec<Heuristic*> reach_heuristics;

    /*
     struct DistLit{
     Lit l;
     int min_distance;

     };*/
    vec<vec<Lit>> dist_lits;

    std::vector<ForceReason> forced_edges;
    struct Change {
        Var v;
        int u;
        bool polarity;
    };
    vec<bool> is_changed_under;
    vec<bool> is_changed_over;
    vec<Change> changed;

    vec<int> to_visit;
    vec<char> seen;
    vec<Lit> extra_conflict;
    vec<int> removed_edges;
    //stats

    int stats_full_updates = 0;
    int stats_fast_updates = 0;
    int stats_fast_failed_updates = 0;
    int stats_skip_deletes = 0;
    int stats_skipped_updates = 0;
    int stats_num_skipable_deletions = 0;
    int stats_learnt_components = 0;
    int stats_learnt_components_sz = 0;
    int64_t stats_heuristic_recomputes = 0;
    double mod_percentage = 0.2;
    int stats_pure_skipped = 0;
    int stats_shrink_removed = 0;
    double stats_full_update_time = 0;
    double stats_fast_update_time = 0;
    Heuristic* conflictingHeuristic = nullptr;

    Heuristic* getConflictingHeuristic() override{
        return conflictingHeuristic;
    }


    void printStats() override{
        //printf("Reach detector\n");
        Detector::printStats();
        if(opt_detect_pure_theory_lits)
            printf("\tPropagations skipped by pure literal detection: %d\n", stats_pure_skipped);
        if(opt_shrink_theory_conflicts){
            printf("\t%d lits removed by shrinking conflicts\n", stats_shrink_removed);
        }
        if(opt_learn_unreachable_component){
            printf("\t%d components learned, average component size: %f\n", stats_learnt_components,
                   stats_learnt_components_sz / (float) stats_learnt_components);
        }
        if(opt_decide_theories){
            printf("\t%" PRId64 " heuristic path recomputations\n", stats_heuristic_recomputes);
        }
        if(overapprox_reach_detector){
            printf("\t\tOverapproxReach: ");
            overapprox_reach_detector->printStats();
        }
        if(overapprox_path_detector){
            printf("\t\tOverapproxPath: ");
            overapprox_path_detector->printStats();
        }
        if(underapprox_path_detector){
            printf("\t\tUnderapproxReach: ");
            underapprox_path_detector->printStats();
        }

        if(underapprox_fast_detector){
            printf("\t\tUnderapproxFast: ");
            underapprox_fast_detector->printStats();
        }
    }

    struct ReachStatus {
        ReachDetector& detector;
        bool polarity;

        void setReachable(int u, bool reachable);

        bool isReachable(int u) const{
            return false;
        }

        void setMininumDistance(int u, bool reachable, Weight distance);

        ReachStatus(ReachDetector& _outer, bool _polarity) :
                detector(_outer), polarity(_polarity){
        }
    };

    ReachStatus* positiveReachStatus = nullptr;
    ReachStatus* negativeReachStatus = nullptr;
    MaxFlow<Weight>* conflict_flow = nullptr;
    std::vector<MaxFlow<Weight>*> conflict_flows;

    Reach* chokepoint_detector = nullptr;

    std::vector<MaxFlowEdge> cut;

    //This class can stand in as a reach algorithm if we have encoded reachability directly into the cnf.
    class CNFReachability : public Reach {
        ReachDetector& detector;
        bool over_approx;
        int num_updates = 0;

    public:
        CNFReachability(ReachDetector& _outer, bool over_approx) :
                detector(_outer), over_approx(over_approx){

        }

        int numUpdates() const override{
            return num_updates;
        }

        void update() override{
            num_updates++;
        }

        void setSource(int s) override{
            assert(s == detector.source);
        }

        int getSource() override{
            return detector.source;
        }

        bool connected_unsafe(int t) override{
            return connected(t);
        }

        bool connected_unchecked(int t) override{
            return connected(t);
        }

        bool connected(int t) override{
            assert(detector.reach_lits[t] != lit_Undef);
            if(over_approx)
                return detector.outer->value(detector.reach_lits[t]) != l_False;
            else
                return detector.outer->value(detector.reach_lits[t]) == l_True;
        }

        int previous(int node) override{
            assert(false);
            throw std::runtime_error("Unsupported operation"); //not supported
        }

        int incomingEdge(int node) override{
            assert(false);
            throw std::runtime_error("Unsupported operation");
        }

    };

    void unassign(Lit l) override{
        Detector::unassign(l);
        int index = var(l) - first_reach_var;
        if(index >= 0 && index < reach_lit_map.size() && reach_lit_map[index] != -1){
            int node = reach_lit_map[index];

            if(!is_changed_under[node]){
                changed.push({var(l), node, true});
                is_changed_under[node] = true;
            }
            if(!is_changed_over[node]){
                changed.push({var(l), node, false});
                is_changed_over[node] = true;
            }
        }
    }

    void activateHeuristic() override{
        for(Heuristic* h:all_reach_heuristics){
            outer->activateHeuristic(h);
        }
    }

    void assign(Lit l) override{
        Detector::assign(l);
        int index = var(l) - first_reach_var;
        if(index >= 0 && index < reach_lit_map.size()){
            int to = reach_lit_map[index];
            if(to >= 0 && to < reach_heuristics.size() && reach_heuristics[to]){
                outer->activateHeuristic(reach_heuristics[to]);
            }
        }
    }

    int getNode(Var reachVar){
        assert(reachVar >= first_reach_var);
        int index = reachVar - first_reach_var;
        assert(index < reach_lit_map.size());
        assert(reach_lit_map[index] >= 0);
        return reach_lit_map[index];
    }

    void backtrack(int level) override{

    }

    /*	Lit getLit(int node){

     return reach_lits[node];

     }*/
    void attachSubHeuristic(Heuristic* h, int to);

    void buildSATConstraints(bool onlyUnderApprox = false, int within_steps = -1);

    bool propagate(vec<Lit>& conflict) override;

    void buildReachReason(int node, vec<Lit>& conflict);

    void buildNonReachReason(int node, vec<Lit>& conflict, bool force_maxflow = false);

    void buildForcedEdgeReason(int reach_node, int forced_edge_id, vec<Lit>& conflict);

    void buildReason(Lit p, vec<Lit>& reason, CRef marker) override;

    bool checkSatisfied() override;

    void printSolution(std::ostream& write_to) override;

    void addLit(int from, int to, Var reach_var);

    Lit decide(CRef& decision_reason) override;

    void preprocess() override;

    void dbg_sync_reachability();

    ReachDetector(int _detectorID, GraphTheorySolver<Weight>* _outer, Graph& g_under, Graph& g_over, Graph& cutGraph,
                  int _source, double seed = 1);

    ~ReachDetector() override{

        if(chokepoint_detector)
            delete chokepoint_detector;
        if(cutgraph_detector)
            delete cutgraph_detector;

        if(positiveReachStatus)
            delete positiveReachStatus;
        if(negativeReachStatus)
            delete negativeReachStatus;

        if(underapprox_path_detector && underapprox_path_detector != underapprox_detector)
            delete underapprox_path_detector;

        if(overapprox_path_detector && overapprox_path_detector != overapprox_reach_detector)
            delete overapprox_path_detector;

        if(underapprox_fast_detector && underapprox_fast_detector != underapprox_path_detector){
            delete underapprox_fast_detector;
        }

        if(negative_distance_detector && negative_distance_detector != overapprox_path_detector)
            delete negative_distance_detector;

        if(underapprox_detector)
            delete underapprox_detector;

        if(overapprox_reach_detector)
            delete overapprox_reach_detector;


        if(conflict_flow)
            delete conflict_flow;

        for(auto* c: conflict_flows){
            if(c)
                delete c;
        }
    }

    std::string getName() override{
        return "Reachability Detector";
    }

    bool dbg_cut(std::vector<MaxFlowEdge>& cut, Graph& graph, int source, int node){
#ifdef DEBUG_GRAPH

        DynamicGraph<int>  t;
        for (int i = 0; i < graph.nodes(); i++)
            t.addNode();

        for (int id = 0; id < graph.edges(); id++) {
            t.addEdge(graph.getEdge(id).from, graph.getEdge(id).to, id);
            if (id % 2 == 0) {
                bool incut = false;
                for (int i = 0; i < cut.size(); i++) {
                    if (cut[i].id == id) {
                        incut = true;
                        break;
                    }
                }
                if (incut) {
                    t.setEdgeWeight(id, 0);
                    t.disableEdge(id);
                } else
                    t.setEdgeWeight(id, 1);
            } else {
                t.setEdgeWeight(id, 0xFFFF);
            }

        }
        EdmondsKarpAdj<int> check(t,  source, node);
        std::vector<MaxFlowEdge> check_cut;
        int flow = check.minCut(check_cut);
        assert(flow < 0xFFFF);
        if (flow < 0xFFFF) {
            throw std::runtime_error("bad cut");
        }

#endif
        return true;
    }

    int getReachNode(Lit reachLit);

    bool isConnected(int node, bool overapprox);

    bool isConnected(Lit reachLit, bool overapprox);

    //Return the path (in terms of nodes)
    bool getModel_Path(int node, std::vector<int>& store_path);

    bool getModel_PathByEdgeLit(int node, std::vector<Lit>& store_path);
};
};
#endif /* REACHDETECTOR_H_ */
