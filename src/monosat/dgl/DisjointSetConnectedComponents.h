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

#ifndef DISJOINTSETSCONNECTEDCOMPONENTS_H_
#define DISJOINTSETSCONNECTEDCOMPONENTS_H_

#include <vector>
#include "monosat/dgl/alg/Heap.h"
#include "Graph.h"
#include "monosat/core/Config.h"
#include "ConnectedComponents.h"
#include "monosat/dgl/alg/DisjointSets.h"
#include <limits>

namespace dgl {
template<typename Weight, class Status = ConnectedComponents::NullConnectedComponentsStatus>
class DisjointSetsConnectedComponents : public ConnectedComponents {
public:

    Graph<Weight>& g;
    Status& status;
    int last_modification;

    int last_addition;
    int last_deletion;
    int history_qhead;

    int last_history_clear;
    bool hasParents = false;
    int INF;
    DisjointSets sets;

    std::vector<int> q;
    std::vector<int> check;
    const int reportPolarity;
    struct ConnectCheck {
        int u;
        int v;
    };
    std::vector<ConnectCheck> connectChecks;

    //stats

    int stats_full_updates = 0;
    int stats_fast_updates = 0;
    int stats_fast_failed_updates = 0;
    int stats_skip_deletes = 0;
    int stats_skipped_updates = 0;
    int stats_num_skipable_deletions = 0;
    double mod_percentage = 0.2;

    double stats_full_update_time = 0;
    double stats_fast_update_time = 0;

public:
    DisjointSetsConnectedComponents(Graph<Weight>& graph, Status& _status, int _reportPolarity = 0) :
            g(graph), status(_status), last_modification(-1), last_addition(-1), last_deletion(-1), history_qhead(0),
            last_history_clear(
                    0), INF(0), reportPolarity(_reportPolarity){

    }

    DisjointSetsConnectedComponents(Graph<Weight>& graph, int _reportPolarity = 0) :
            g(graph), status(nullConnectedComponentsStatus), last_modification(-1), last_addition(-1), last_deletion(
            -1), history_qhead(0), last_history_clear(0), INF(0), reportPolarity(_reportPolarity){

    }

    void setNodes(int n){
        q.reserve(n);
        check.reserve(n);

        INF = std::numeric_limits<int>::max();
        sets.AddElements(n);

    }

    void addConnectedCheck(int u, int v){
        connectChecks.push_back({u, v});
    }

    void update(){
        if(last_modification > 0 && g.getCurrentHistory() == last_modification){
            stats_skipped_updates++;
            return;
        }
        stats_full_updates++;

        if(last_deletion == g.nDeletions()){
            stats_num_skipable_deletions++;
        }
        hasParents = false;
        sets.Reset();
        setNodes(g.nodes());

        for(int i = 0; i < g.edges(); i++){
            if(g.edgeEnabled(i)){
                int u = g.getEdge(i).from;
                int v = g.getEdge(i).to;
                sets.UnionElements(u, v);
            }
        }

        status.setComponents(sets.NumSets());

        for(auto c : connectChecks){
            int u = c.u;
            int v = c.v;
            bool connected = sets.FindSet(u) == sets.FindSet(v);
            if(reportPolarity >= 0 && connected){
                status.setConnected(u, v, true);
            }else if(reportPolarity <= 0 && !connected){
                status.setConnected(u, v, false);
            }
        }

        last_modification = g.getCurrentHistory();
        last_deletion = g.nDeletions();
        last_addition = g.nAdditions();

        history_qhead = g.historySize();
        last_history_clear = g.nHistoryClears();
    }

    bool connected(int from, int to){
        update();
        return sets.FindSet(from) == sets.FindSet(to);
    }

    int numComponents(){
        update();
        return sets.NumSets();
    }

    int getComponent(int node){
        update();
        return sets.FindSet(node);
    }

    int getElement(int set){
        update();
        return sets.GetElement(set);
    }

    bool dbg_uptodate(){
        return true;
    };

};
};

#endif /* DISJOINTSETCONNECTEDCOMPONENTS_H_ */
