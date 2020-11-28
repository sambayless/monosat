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

#ifndef DISTANCEETECTOR_H_
#define DISTANCEETECTOR_H_

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
class DistanceDetector : public Detector {
public:
    GraphTheorySolver<Weight>* outer;

    Graph& g_under;
    Graph& g_over;
    Graph& cutGraph;
    //int within;
    int source;
    double rnd_seed;
    int constraintsBuilt;
    CRef unweighted_underprop_marker = CRef_Undef;
    CRef unweighted_overprop_marker = CRef_Undef;


    Distance<int>* underapprox_unweighted_distance_detector = nullptr;
    Distance<int>* overapprox_unweighted_distance_detector = nullptr;


    Reach* underapprox_path_detector = nullptr;

    Var first_reach_var = var_Undef;

    struct ReachLit {
        int to;
        int within;

    };
    vec<ReachLit> reach_lit_map;
    bool has_unweighted_shortest_paths_overapprox = false;

    vec<int> unweighted_over_approx_shortest_paths;
    MaxFlow<Weight>* conflict_flow = nullptr;


    int max_unweighted_distance = 0;

    int64_t stats_pure_skipped = 0;
    int64_t stats_distance_gt_reasons = 0;
    int64_t stats_distance_leq_reasons = 0;
    int64_t stats_unweighted_gt_reasons = 0;
    int64_t stats_unweighted_leq_reasons = 0;
    int64_t stats_gt_unweighted_edges_skipped = 0;
    int64_t stats_gt_weighted_edges_skipped = 0;


    vec<vec<Lit>> unweighted_sat_lits;

    struct UnweightedDistLit {
        Lit l;
        int min_unweighted_distance;

    };

    vec<vec<UnweightedDistLit>> unweighted_dist_lits;


    struct Change {
        int u;
        bool polarity;
    };
    vec<Change> changed;
    vec<bool> is_changed_under;
    vec<bool> is_changed_over;
    std::vector<double> rnd_weight;

    WeightedDijkstra<Weight, Graph, double>* rnd_path;


    struct ReachStatus {
        DistanceDetector& detector;
        bool polarity;

        void setReachable(int u, bool reachable);

        bool isReachable(int u) const{
            return false;
        }

        void setMininumDistance(int u, bool reachable, int distance);

        ReachStatus(DistanceDetector& _outer, bool _polarity) :
                detector(_outer), polarity(_polarity){
        }
    };

    struct DistanceStatus {
        DistanceDetector& detector;
        bool polarity;

        void setReachable(int u, bool reachable);

        bool isReachable(int u) const{
            return false;
        }

        void setMininumDistance(int u, bool reachable, Weight& distance);

        DistanceStatus(DistanceDetector& _outer, bool _polarity) :
                detector(_outer), polarity(_polarity){
        }
    };

    std::vector<MaxFlowEdge> cut;

    ReachStatus* positiveReachStatus;
    ReachStatus* negativeReachStatus;

    vec<int> to_visit;
    vec<char> seen;

    int getNode(Var reachVar){
        assert(reachVar >= first_reach_var);
        int index = reachVar - first_reach_var;
        assert(index < reach_lit_map.size());
        assert(reach_lit_map[index].to >= 0);
        return reach_lit_map[index].to;
    }

    int getMaximumDistance(Var reachVar){
        assert(reachVar >= first_reach_var);

        int index = reachVar - first_reach_var;
        assert(index < reach_lit_map.size());

        assert(reach_lit_map[index].within >= 0);
        return reach_lit_map[index].within;
    }

    void printStats() override{
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
        int index = var(l) - first_reach_var;

        //at the moment, change in assignments are only tracked this way for unweighted lits:
        if(index >= 0 && index < reach_lit_map.size() && reach_lit_map[index].to != -1){
            int node = reach_lit_map[index].to;
            if(!is_changed_under[node]){
                changed.push({node, true});
                is_changed_under[node] = true;
            }
            if(!is_changed_over[node]){
                changed.push({node, false});
                is_changed_over[node] = true;
            }
        }
    }

    void preprocess() override;

    bool propagate(vec<Lit>& conflict) override;

    void buildUnweightedDistanceLEQReason(int node, vec<Lit>& conflict);

    void buildUnweightedDistanceGTReason(int node, int within_steps, vec<Lit>& conflict);

    void buildReason(Lit p, vec<Lit>& reason, CRef marker) override;

    bool checkSatisfied() override;

    Lit decide(CRef& decision_reason) override;

    void updateShortestPaths();

    void addUnweightedShortestPathLit(int from, int to, Var reach_var, int within_steps = -1);

    bool getModel_Path(int node, std::vector<int>& store_path);

    bool getModel_PathByEdgeLit(int node, std::vector<Lit>& store_path);

    DistanceDetector(int _detectorID, GraphTheorySolver<Weight>* _outer,
                     Graph& g_under, Graph& g_over, Graph& cutGraph, int _source,
                     double seed = 1);//:Detector(_detectorID),outer(_outer),within(-1),source(_source),rnd_seed(seed),positive_reach_detector(NULL),negative_reach_detector(NULL),positive_path_detector(NULL),positiveReachStatus(NULL),negativeReachStatus(NULL){}
    ~DistanceDetector() override{

        if(overapprox_unweighted_distance_detector)
            delete overapprox_unweighted_distance_detector;

        if(underapprox_path_detector && underapprox_path_detector != underapprox_unweighted_distance_detector)
            delete underapprox_path_detector;


        if(positiveReachStatus)
            delete positiveReachStatus;
        if(negativeReachStatus)
            delete negativeReachStatus;

        if(underapprox_unweighted_distance_detector)
            delete underapprox_unweighted_distance_detector;


        if(rnd_path)
            delete rnd_path;

        if(conflict_flow)
            delete conflict_flow;
    }

    std::string getName() override{
        return "Shortest Path Detector";
    }

private:
    void buildUnweightedSATConstraints(bool onlyUnderApprox, int distance = -1);
};

};
#endif /* DistanceDetector_H_ */
