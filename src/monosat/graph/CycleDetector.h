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
#ifndef CYCLE_DETECTOR_H_
#define CYCLE_DETECTOR_H_

#include "monosat/utils/System.h"

#include "GraphTheoryTypes.h"
#include "monosat/dgl/DynamicGraph.h"
#include "monosat/dgl/ConnectedComponents.h"
#include "monosat/dgl/Cycle.h"
#include "monosat/core/SolverTypes.h"
#include "monosat/mtl/Map.h"

#include "monosat/utils/System.h"
#include "Detector.h"

using namespace dgl;
namespace Monosat {
template<typename Weight>
class GraphTheorySolver;

template<typename Weight, typename Graph = DynamicGraph<Weight>>
class CycleDetector : public Detector {
public:
    GraphTheorySolver<Weight>* outer;
    //int within;
    Graph& g_under;
    Graph& g_over;

    double rnd_seed;
    CRef directed_cycle_marker;
    CRef no_directed_cycle_marker;

    CRef undirected_cycle_marker;
    CRef no_undirected_cycle_marker;
    Cycle* underapprox_directed_cycle_detector = nullptr;
    Cycle* overapprox_directed_cycle_detector = nullptr;
    Cycle* underapprox_undirected_cycle_detector = nullptr;
    Cycle* overapprox_undirected_cycle_detector = nullptr;
    //Reach *  positive_path_detector;

    //vec<Lit>  reach_lits;

    vec<Var> tmp_nodes;
    vec<bool> seen;
    vec<bool> black;
    vec<int> ancestors;

    vec<bool> edge_in_clause;
    vec<int> visit;

    Lit undirected_acyclic_lit = lit_Undef;
    Lit directed_acyclic_lit = lit_Undef;

    bool propagate(vec<Lit>& conflict) override;

    void buildNoUndirectedCycleReason(vec<Lit>& conflict);

    void buildNoDirectedCycleReason(vec<Lit>& conflict);

    void buildUndirectedCycleReason(vec<Lit>& conflict);

    void buildDirectedCycleReason(vec<Lit>& conflict);

    //void buildForcedMinWeightReason(int reach_node, int forced_edge_id,vec<Lit> & conflict);
    void buildReason(Lit p, vec<Lit>& reason, CRef marker) override;

    bool checkSatisfied() override;

    Lit decide(CRef& decision_reason) override;

    Lit addAcyclicLit(bool directed, Var v = var_Undef);

    bool hasAcyclicLit(bool directed){
        if(directed){
            return directed_acyclic_lit != lit_Undef;
        }else{
            return undirected_acyclic_lit != lit_Undef;
        }
    }

    CycleDetector(int _detectorID, GraphTheorySolver<Weight>* _outer, Graph& _g, Graph& _antig,
                  bool detect_directed_cycles = true,
                  double seed = 1); //:Detector(_detectorID),outer(_outer),within(-1),source(_source),rnd_seed(seed),positive_reach_detector(NULL),negative_reach_detector(NULL),positive_path_detector(NULL),positiveReachStatus(NULL),negativeReachStatus(NULL){}
    ~CycleDetector() override{
        if(overapprox_undirected_cycle_detector &&
           overapprox_undirected_cycle_detector != overapprox_directed_cycle_detector)
            delete overapprox_undirected_cycle_detector;
        if(underapprox_undirected_cycle_detector &&
           underapprox_undirected_cycle_detector != underapprox_directed_cycle_detector)
            delete underapprox_undirected_cycle_detector;
        if(underapprox_directed_cycle_detector)
            delete underapprox_directed_cycle_detector;
        if(overapprox_directed_cycle_detector)
            delete overapprox_directed_cycle_detector;

    }

    std::string getName() override{
        return "Cycle Detector";
    }
};
};
#endif /* DistanceDetector_H_ */
