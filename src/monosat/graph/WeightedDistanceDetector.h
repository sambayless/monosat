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

#ifndef WEIGHTED_DISTANCEETECTOR_H_
#define WEIGHTED_DISTANCEETECTOR_H_

#include "monosat/utils/System.h"

#include "monosat/graph/GraphTheoryTypes.h"
#include "monosat/dgl/Graph.h"
#include "monosat/dgl/Reach.h"
#include "monosat/dgl/Distance.h"
#include "monosat/dgl/Dijkstra.h"
#include "monosat/dgl/BFS.h"
#include "monosat/dgl/MaxFlow.h"
#include "monosat/core/SolverTypes.h"
#include "monosat/mtl/Map.h"
#include "monosat/graph/WeightedDijkstra.h"
#include <gmpxx.h>
#include "monosat/utils/System.h"
#include "monosat/graph/Detector.h"
#include "monosat/bv/BVTheorySolver.h"
#include <vector>

using namespace dgl;
namespace Monosat {
template<typename Weight>
class GraphTheorySolver;

template<typename Weight, typename Graph = DynamicGraph<Weight>>
class WeightedDistanceDetector : public Detector {
public:
    GraphTheorySolver<Weight>* outer;

    Graph& g_under;
    Graph& g_over;
    Graph& cutGraph;
    //int within;
    int source;
    double rnd_seed;
    int constraintsBuilt;

    CRef weighted_underprop_marker = CRef_Undef;
    CRef weighted_overprop_marker = CRef_Undef;
    CRef weighted_underprop_bv_marker = CRef_Undef;
    CRef weighted_overprop_bv_marker = CRef_Undef;


    Distance<Weight>* underapprox_weighted_distance_detector = nullptr;
    Distance<Weight>* overapprox_weighted_distance_detector = nullptr;
    Distance<Weight>* underapprox_weighted_path_detector = nullptr;


    //vec<Lit>  reach_lits;
    Var first_reach_var = var_Undef;
    enum DistLitType {
        None, WeightedConstLit, WeightedBVLit
    };
    struct ReachLit {
        int to;
        int within;
        DistLitType type;
    };
    vec<ReachLit> reach_lit_map;
    vec<int> force_reason;

    bool has_weighted_shortest_paths_overapprox = false;
    vec<int> unweighted_over_approx_shortest_paths;
    vec<Weight> over_approx_shortest_paths;
    MaxFlow<Weight>* conflict_flow = nullptr;

    vec<int> to_visit;
    vec<char> seen;

    int max_unweighted_distance;

    int64_t stats_pure_skipped = 0;
    int64_t stats_distance_gt_reasons = 0;
    int64_t stats_distance_leq_reasons = 0;
    int64_t stats_unweighted_gt_reasons = 0;
    int64_t stats_unweighted_leq_reasons = 0;
    int64_t stats_gt_unweighted_edges_skipped = 0;
    int64_t stats_gt_weighted_edges_skipped = 0;

    BVTheorySolver<Weight>* bvTheory = nullptr;

    class DistanceOp : public GraphTheorySolver<Weight>::GraphTheoryOp {
        WeightedDistanceDetector* outer;
        int bvID;
        int to;
        bool strictCompare;
        Lit comparisonLit;
    public:
        DistanceOp(BVTheorySolver<Weight>& theory, WeightedDistanceDetector* outer, int bvID, int to,
                   bool strictCompare, Lit comparisonLit) : GraphTheorySolver<Weight>::GraphTheoryOp(theory,
                                                                                                     outer->outer),
                                                            outer(outer), bvID(bvID), to(to),
                                                            strictCompare(strictCompare), comparisonLit(comparisonLit){

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


    struct WeightedDistLit {
        Lit l;
        int u;
        Weight min_distance;
        bool strictComparison;

        bool operator<(WeightedDistLit& b) const{
            return min_distance < b.min_distance;
        }
    };

    vec<WeightedDistLit> weighted_dist_lits;


    struct WeightedDistBVLit {
        Lit l;
        int u;
        BitVector <Weight> bv;
        bool strictComparison;
        DistanceOp* op;
    };
    vec<WeightedDistBVLit> weighted_dist_bv_lits;

    struct Change {
        //Var v;
        int u;
        //int min_distance;
    };
    vec<Change> changed;
    vec<bool> is_changed_under;
    vec<bool> is_changed_over;
    vec<Var> tmp_nodes;

    std::vector<double> rnd_weight;

    WeightedDijkstra<Weight, Graph, double>* rnd_path;

    struct ReachStatus {
        WeightedDistanceDetector& detector;
        bool polarity;

        void setReachable(int u, bool reachable);

        bool isReachable(int u) const{
            return false;
        }

        void setMininumDistance(int u, bool reachable, int distance);

        ReachStatus(WeightedDistanceDetector& _outer, bool _polarity) :
                detector(_outer), polarity(_polarity){
        }
    };

    struct DistanceStatus {
        WeightedDistanceDetector& detector;
        bool polarity;

        void setReachable(int u, bool reachable);

        bool isReachable(int u) const{
            return false;
        }

        void setMininumDistance(int u, bool reachable, Weight& distance);

        DistanceStatus(WeightedDistanceDetector& _outer, bool _polarity) :
                detector(_outer), polarity(_polarity){
        }
    };

    std::vector<MaxFlowEdge> cut;

    DistanceStatus* positiveDistanceStatus;
    DistanceStatus* negativeDistanceStatus;

    DistLitType getLitType(Lit reachLit){
        return getLitType(var(reachLit));
    }

    DistLitType getLitType(Var reachVar){
        assert(reachVar >= first_reach_var);
        int index = reachVar - first_reach_var;
        assert(index < reach_lit_map.size());
        assert(reach_lit_map[index].to >= 0);
        return reach_lit_map[index].type;
    }

    WeightedDistLit& getDistLit(Var v){
        assert(v >= first_reach_var);
        int index = v - first_reach_var;
        assert(index < reach_lit_map.size());
        assert(reach_lit_map[index].type == WeightedConstLit);
        assert(reach_lit_map[index].to >= 0);
        int w_index = reach_lit_map[index].within;
        assert(w_index >= 0);
        assert(w_index < weighted_dist_lits.size());
        assert(var(weighted_dist_lits[w_index].l) == v);
        return weighted_dist_lits[w_index];
    }

    WeightedDistBVLit& getBVDistLit(Var v){
        assert(v >= first_reach_var);
        int index = v - first_reach_var;
        assert(index < reach_lit_map.size());
        assert(reach_lit_map[index].type == WeightedBVLit);
        assert(reach_lit_map[index].to >= 0);
        int w_index = reach_lit_map[index].within;
        assert(w_index >= 0);
        assert(w_index < weighted_dist_bv_lits.size());
        assert(var(weighted_dist_bv_lits[w_index].l) == v);
        return weighted_dist_bv_lits[w_index];
    }

    int getNode(Var reachVar){
        assert(reachVar >= first_reach_var);
        int index = reachVar - first_reach_var;
        assert(index < reach_lit_map.size());
        assert(reach_lit_map[index].to >= 0);
        return reach_lit_map[index].to;
    }

    void printStats() override{
        //printf("Distance detector\n");
        Detector::printStats();
        if(opt_verb > 0){
            if(opt_detect_pure_theory_lits)
                printf("\tPropagations skipped by pure literal detection: %" PRId64 "\n", stats_pure_skipped);
            printf("\tUnweighted Reasons (leq,gt): %" PRId64 ",%" PRId64 "\n", stats_unweighted_leq_reasons,
                   stats_unweighted_gt_reasons);
            printf("\tWeighted Reasons (leq,gt): %" PRId64 ",%" PRId64 "\n", stats_distance_leq_reasons,
                   stats_distance_gt_reasons);
            printf("\tConflict Edges Skipped (unweighted %" PRId64 ", weighted %" PRId64 ")\n",
                   stats_gt_unweighted_edges_skipped,
                   stats_gt_weighted_edges_skipped);
        }
    }

    void printSolution(std::ostream& write_to) override;

    void unassign(Lit l) override{
        Detector::unassign(l);
    }

    void preprocess() override;

    bool propagate(vec<Lit>& conflict) override;

    void buildDistanceLEQReason(int to, Weight& min_distance, vec<Lit>& conflict, bool strictComparison = false);

    void buildDistanceGTReason(int to, Weight& min_distance, vec<Lit>& conflict, bool strictComparison = true);

    void analyzeDistanceLEQReason(int to, Weight& min_distance, vec<Lit>& conflict, bool strictComparison = true);

    void analyzeDistanceGTReason(int to, Weight& min_distance, vec<Lit>& conflict, bool strictComparison = true);

    void buildReason(Lit p, vec<Lit>& reason, CRef marker) override;

    bool checkSatisfied() override;

    Lit decide(CRef& decision_reason) override;

    void updateShortestPaths();

    void addWeightedShortestPathLit(int from, int to, Var reach_var, Weight within_distance, bool strictComparison);

    void
    addWeightedShortestPathBVLit(int from, int to, Var reach_var, const BitVector <Weight>& bv, bool strictComparison);

    bool getModel_Path(int node, std::vector<int>& store_path);

    bool getModel_PathByEdgeLit(int node, std::vector<Lit>& store_path);

    WeightedDistanceDetector(int _detectorID, GraphTheorySolver<Weight>* _outer,
                             Graph& _g, Graph& _antig, Graph& cutGraph, int _source,
                             double seed = 1);//:Detector(_detectorID),outer(_outer),within(-1),source(_source),rnd_seed(seed),positive_reach_detector(NULL),negative_reach_detector(NULL),positive_path_detector(NULL),positiveReachStatus(NULL),negativeReachStatus(NULL){}
    ~WeightedDistanceDetector() override{
        if(positiveDistanceStatus)
            delete positiveDistanceStatus;
        if(negativeDistanceStatus)
            delete negativeDistanceStatus;


        if(underapprox_weighted_path_detector &&
           underapprox_weighted_path_detector != underapprox_weighted_distance_detector)
            delete underapprox_weighted_path_detector;


        if(underapprox_weighted_distance_detector)
            delete underapprox_weighted_distance_detector;

        if(overapprox_weighted_distance_detector)
            delete overapprox_weighted_distance_detector;


        if(rnd_path)
            delete rnd_path;

        if(conflict_flow)
            delete conflict_flow;
    }

    std::string getName() override{
        return "Weighted Distance Detector";
    }

private:
    void buildUnweightedSATConstraints(bool onlyUnderApprox, int distance = -1);
};

};
#endif /* DistanceDetector_H_ */
