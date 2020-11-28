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

#ifndef STEINER_DETECTOR_H_
#define STEINER_DETECTOR_H_

#include "monosat/utils/System.h"
#include "GraphTheoryTypes.h"
#include "monosat/dgl/Graph.h"
#include "monosat/dgl/MinimumSpanningTree.h"
#include "monosat/dgl/SteinerTree.h"
#include "monosat/core/SolverTypes.h"
#include "monosat/mtl/Map.h"

#include "monosat/dgl/alg/DisjointSets.h"
#include "monosat/dgl/DynamicNodes.h"
#include "monosat/utils/System.h"
#include "Detector.h"

using namespace dgl;
namespace Monosat {
template<typename Weight>
class GraphTheorySolver;

template<typename Weight = int, typename Graph=DynamicGraph<Weight>>
class SteinerDetector : public Detector {
public:
    GraphTheorySolver<Weight>* outer;

    Graph& g_under;
    Graph& g_over;

    DynamicNodes underTerminalSet;
    DynamicNodes overTerminalSet;
    double rnd_seed;
    CRef underprop_marker;
    CRef overprop_marker;
    CRef underprop_edge_marker;
    CRef overprop_edge_marker;

    SteinerTree<Weight>* underapprox_detector = nullptr;
    SteinerTree<Weight>* overapprox_detector = nullptr;
    SteinerTree<Weight>* underapprox_conflict_detector = nullptr;
    SteinerTree<Weight>* overapprox_conflict_detector = nullptr;

    vec<Var> terminal_map;
    vec<int> terminal_var_map;

    struct WeightLit {
        Lit l;
        Weight min_weight;

        WeightLit() :
                l(lit_Undef), min_weight(-1){
        }
    };

    bool checked_unique;
    bool all_unique;
    vec<WeightLit> weight_lits;

    Var first_reach_var;
    struct ChangedWeight {
        Lit l;
        Weight weight;
    };
    vec<ChangedWeight> changed_weights;
    vec<Var> tmp_nodes;
    vec<bool> seen;
    vec<bool> black;
    vec<int> ancestors;

    vec<Lit> tmp_conflict;
    vec<int> visit;
    DisjointSets sets;

    struct SteinerStatus {
        SteinerDetector& detector;
        bool polarity;

        void setMinimumSteinerTree(Weight& weight);

        SteinerStatus(SteinerDetector& _outer, bool _polarity) :
                detector(_outer), polarity(_polarity){
        }
    };

    SteinerStatus* negativeStatus;
    SteinerStatus* positiveStatus;

    bool propagate(vec<Lit>& conflict) override;

    void buildMinWeightTooSmallReason(Weight& weight, vec<Lit>& conflict);

    void buildMinWeightTooLargeReason(Weight& weight, vec<Lit>& conflict);

    void buildReason(Lit p, vec<Lit>& reason, CRef marker) override;

    bool checkSatisfied() override;

    Lit decide(CRef& decision_reason) override;

    void addWeightLit(Weight& min_weight, Var weight_var);

    void addTerminalNode(int node, Var theoryVar);

    SteinerDetector(int _detectorID, GraphTheorySolver<Weight>* _outer, Graph& _g,
                    Graph& _antig,
                    double seed = 1); //:Detector(_detectorID),outer(_outer),within(-1),source(_source),rnd_seed(seed),positive_reach_detector(NULL),negative_reach_detector(NULL),positive_path_detector(NULL),positiveReachStatus(NULL),negativeReachStatus(NULL){}

    ~SteinerDetector() override{
        delete underapprox_detector;
        delete overapprox_detector;
    }

    void assign(Lit l) override{
        Detector::assign(l);
        if(var(l) < terminal_var_map.size()){
            int node = terminal_var_map[var(l)];
            if(node >= 0){
                if(sign(l)){
                    underTerminalSet.setNodeEnabled(node, false);
                    overTerminalSet.setNodeEnabled(node, false);
                }else{
                    underTerminalSet.setNodeEnabled(node, true);
                    overTerminalSet.setNodeEnabled(node, true);
                }
            }
        }
    }

    void unassign(Lit l) override{
        Detector::unassign(l);
        if(var(l) < terminal_var_map.size()){
            int node = terminal_var_map[var(l)];
            if(node >= 0){
                underTerminalSet.setNodeEnabled(node, false);
                overTerminalSet.setNodeEnabled(node, true);
            }
        }
    }

    std::string getName() override{
        return "Steiner Detector";
    }

};
};
#endif /* DistanceDetector_H_ */
