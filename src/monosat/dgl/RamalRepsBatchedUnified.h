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

#ifndef RAMAL_REPS_BATCHED_UNIFIED_H_
#define RAMAL_REPS_BATCHED_UNIFIED_H_

#include <monosat/dgl/alg/Heap.h>
#include <monosat/dgl/alg/Rnd.h>
#include <monosat/dgl/Dijkstra.h>
#include <monosat/dgl/Distance.h>
#include <monosat/dgl/Graph.h>
#include <monosat/dgl/Reach.h>
//#include "monosat/core/Config.h"
//#include <algorithm>
#include <cassert>
#include <cstdio>
//#include <exception>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
/*#ifndef NDEBUG
#define DEBUG_RAMAL
#define DEBUG_RAMAL2
#endif*/
namespace dgl {
template<typename Weight = int, typename Graph = DynamicGraph<Weight>, class Status = typename Distance<Weight>::NullStatus>
class RamalRepsBatchedUnified : public Distance<Weight>, public DynamicGraphAlgorithm {
public:
    static bool ever_warned_about_zero_weights;
    Graph& g;
    std::vector<Weight>& weights;
    std::vector<Weight> local_weights;
    Status& status;
    int reportPolarity;
    bool reportDistance;
    uint64_t stats_updates = 0;
    uint64_t stats_resets = 0;
    uint64_t stats_all_updates = 0;

    int last_modification = -1;
    int last_addition = -1;
    int last_deletion = -1;
    int history_qhead = 0;

    int last_history_clear = 0;

    int source;
    Weight INF;

    std::vector<Weight> old_dist;
    std::vector<int> changed;
    std::vector<bool> node_changed;
    std::vector<Weight> dist;
    std::vector<int> prev;

    struct DistCmp {
        std::vector<Weight>& _dist;

        bool operator()(int a, int b) const{
            return _dist[a] < _dist[b];
        }

        DistCmp(std::vector<Weight>& d) :
                _dist(d){
        };
    };

    alg::Heap<DistCmp> q_batch;


    std::vector<int> edgeInShortestPathGraph; //edgeInShortestPathGraph(v,u) is true if there exists any shortest path from s to u that includes the edge (v,u)
    std::vector<int> delta; //delta records, for each node u, the number of incoming edges (v ,u) for which edgeInShortestPathGraph (v,u) is true
    std::vector<int> changeset;
    std::vector<bool> in_changeset;
    int alg_id = -1;

    struct LocalDistanceStatus {
        RamalRepsBatchedUnified& outer;

        void setReachable(int u, bool reachable){

        }

        bool isReachable(int u) const{
            return false;
        }

        void setMininumDistance(int u, bool reachable, Weight& distance){
            if(reachable){
                if(outer.dist[u] != distance){
                    outer.dist[u] = distance;
                    if(!outer.node_changed[u]){
                        outer.node_changed[u] = true;
                        outer.changed.push_back(u);
                    }
                }
            }else{
                if(outer.dist[u] != outer.INF){
                    outer.dist[u] = outer.INF;
                    if(!outer.node_changed[u]){
                        outer.node_changed[u] = true;
                        outer.changed.push_back(u);
                    }
                }
            }
        }

        LocalDistanceStatus(RamalRepsBatchedUnified& _outer) :
                outer(_outer){
        }
    } local_distance_status;

    Dijkstra<Weight, Graph, Status> dijkstras;
    bool has_zero_weights = false;
public:

    int64_t stats_full_updates = 0;
    int64_t stats_fast_updates = 0;
    int64_t stats_fast_failed_updates = 0;
    int64_t stats_skip_deletes = 0;
    int64_t stats_skipped_updates = 0;
    int64_t stats_num_skipable_deletions = 0;
    double mod_percentage = 0;

    double stats_full_update_time = 0;
    double stats_fast_update_time = 0;

    RamalRepsBatchedUnified(int s, Graph& graph, Status& status, int reportPolarity = 0,
                            bool reportDistance = false) :
            g(graph), weights(g.getWeights()), status(status), reportPolarity(reportPolarity),
            reportDistance(reportDistance), last_modification(
            -1), last_addition(-1), last_deletion(-1), history_qhead(0), last_history_clear(0), source(s), INF(
            0), q_batch(DistCmp(dist)), local_distance_status(*this), dijkstras(s, graph, status, reportPolarity){

        mod_percentage = 0.2;
        alg_id = g.addDynamicAlgorithm(this);
    }

    RamalRepsBatchedUnified(int s, Graph& graph, int reportPolarity = 0, bool reportDistance = false) :
            g(graph), weights(g.getWeights()), status(Distance<Weight>::nullStatus), reportPolarity(reportPolarity),
            reportDistance(reportDistance), last_modification(
            -1), last_addition(-1), last_deletion(-1), history_qhead(0), last_history_clear(0), source(s), INF(
            0), q_batch(DistCmp(dist)), local_distance_status(*this), dijkstras(s, graph, status, reportPolarity){

        mod_percentage = 0.2;
        alg_id = g.addDynamicAlgorithm(this);
    }

    std::string getName() override{
        return "RamalRepsBatchedUnified(" + std::to_string(getSource()) + ")";
    }

    void setSource(int s) override{
        source = s;
        last_modification = -1;
        last_addition = -1;
        last_deletion = -1;
    }

    int getSource() override{
        return source;
    }

    std::vector<int>& getChanged(){
        return changed;
    }

    void clearChanged(){
        changed.clear();
    }

    void drawFull(){

    }

    void dbg_delta(){
#ifdef DEBUG_RAMAL
        dbg_delta_lite();
        assert(delta.size() == g.nodes());

        for (int i = 0; i < g.nEdgeIDs(); i++) {
            if (!g.edgeEnabled(i)) {
                assert(!edgeInShortestPathGraph[i]);
            }
        }

        std::vector<int> dbg_delta;
        std::vector<Weight> dbg_dist;
        dbg_dist.resize(g.nodes(), INF);
        dbg_delta.resize(g.nodes());
        dbg_dist[getSource()] = 0;
        struct DistCmp {
            std::vector<Weight> & _dist;
            bool operator()(int a, int b) const {
                return _dist[a] < _dist[b];
            }
            DistCmp(std::vector<Weight> & d) :
                    _dist(d) {
            }
            ;
        };
        alg::Heap<DistCmp> q(dbg_dist);

        q.insert(getSource());

        while (q.size()) {
            int u = q.removeMin();
            if (dbg_dist[u] == INF)
                break;
            dbg_delta[u] = 0;

            for (int i = 0; i < g.nIncoming(u); i++) {
                if (!g.edgeEnabled(g.incoming(u, i).id))
                    continue;

                int edgeID = g.incoming(u, i).id;
                int v = g.getEdge(edgeID).from;
                Weight alt = dbg_dist[v] + weights[edgeID];
                assert(alt >= dbg_dist[u]);
                /*	if (alt==dbg_dist[u]){
                 dbg_delta[u]++;
                 }*/
            }

            for (int i = 0; i < g.nIncident(u); i++) {
                if (!g.edgeEnabled(g.incident(u, i).id))
                    continue;

                int edgeID = g.incident(u, i).id;
                int v = g.getEdge(edgeID).to;
                Weight alt = dbg_dist[u] + weights[edgeID];
                if (alt < dbg_dist[v]) {

                    dbg_dist[v] = alt;

                    if (!q.inHeap(v))
                        q.insert(v);
                    else
                        q.decrease(v);
                }/*else if (alt==dbg_dist[v]){
				 dbg_delta[v]++;
				 }*/
            }
        }
        Dijkstra<Weight> d2(source,g);
        for (int u = 0; u < g.nodes(); u++) {
            Weight & d = dist[u];
            if(d2.connected(u)) {
                assert(d2.distance(u) == dist[u]);
            }else{
                assert(dist[u]>=INF);
            }
            assert(dbg_dist[u] == dist[u]);
            for (int i = 0; i < g.nIncoming(u); i++) {
                if (!g.edgeEnabled(g.incoming(u, i).id))
                    continue;

                int edgeID = g.incoming(u, i).id;
                int v = g.getEdge(edgeID).from;

                Weight alt = dbg_dist[v] + weights[edgeID];
                assert(alt >= dbg_dist[u]);
                if (alt == dbg_dist[u]) {
                    dbg_delta[u]++; //is this right?
                    assert(edgeInShortestPathGraph[edgeID]);
                } else {
                    assert(!edgeInShortestPathGraph[edgeID]);
                }

            }

        }
        for (int u = 0; u < g.nodes(); u++) {
            Weight& d = dist[u];
            Weight& d_expect = dbg_dist[u];
            assert(d == dbg_dist[u]);

        }
        for (int u = 0; u < g.nodes(); u++) {
            int d = delta[u];
            int d_expect = dbg_delta[u];
            assert(d == dbg_delta[u]);

        }
        dbg_delta_lite();
#endif
    }

    //Called when an edge is enabled
    void AddEdge(int edgeID){
        assert(g.edgeEnabled(edgeID));
        if(edgeInShortestPathGraph[edgeID])
            return;
        int ru = g.getEdge(edgeID).from;
        int rv = g.getEdge(edgeID).to;
        Weight& rdv = dist[rv];
        Weight& rdu = dist[ru];

        Weight& weight = weights[edgeID];
        if(dist[rv] < (dist[ru] + weight))
            return;
        else if(dist[rv] == (dist[ru] + weight)){
            assert(!edgeInShortestPathGraph[edgeID]);
            edgeInShortestPathGraph[edgeID] = true;
            delta[rv]++; //we have found an alternative shortest path to v
            return;
        }
        edgeInShortestPathGraph[edgeID] = true;
        dist[rv] = dist[ru] + weight;
        assert(delta[rv] >= 0);
        delta[rv]++;
        q_batch.update(rv);
        //maintain delta invariant
        for(int j = 0; j < g.nIncoming(rv); j++){
            auto& e = g.incoming(rv, j);
            int adjID = e.id;

            if(g.edgeEnabled(adjID)){
                if(edgeInShortestPathGraph[adjID]){
                    assert(g.getEdge(adjID).to == rv);
                    int v = g.getEdge(adjID).from;
                    Weight& w = weights[adjID]; //assume a weight of one for now
                    Weight alt = dist[v] + w;

                    if(dist[rv] < alt || alt >= INF){
                        edgeInShortestPathGraph[adjID] = false;
                        delta[rv]--;
                    }
                }
            }

        }
        assert(delta[rv] >= 1);
        assert(edgeInShortestPathGraph[edgeID]);
    }

    void dbg_delta_lite(){
#ifdef DEBUG_RAMAL
        for (int u = 0; u < g.nodes(); u++) {
            int del = delta[u];
            Weight d = dist[u];
            int num_in = 0;
            for (int i = 0; i < g.nIncoming(u); i++) {
                auto & e = g.incoming(u, i);
                int adjID = e.id;
                int from = g.getEdge(adjID).from;

                Weight dfrom = dist[from];
                if (edgeInShortestPathGraph[adjID])
                    num_in++;
            }
            assert(del == num_in);
        }
#endif

    }

    //Called if an edge weight is decreased
    void DecreaseWeight(int edgeID){

        assert(g.edgeEnabled(edgeID));
        int ru = g.getEdge(edgeID).from;
        int rv = g.getEdge(edgeID).to;
        Weight& rdv = dist[rv];
        Weight& rdu = dist[ru];

        Weight& weight = weights[edgeID];
        assert(weight > 0);
        if(dist[rv] < (dist[ru] + weight)){
            assert(!edgeInShortestPathGraph[edgeID]);
            return;
        }else if(dist[rv] == (dist[ru] + weight)){
            if(!edgeInShortestPathGraph[edgeID]){
                edgeInShortestPathGraph[edgeID] = true;
                delta[rv]++; //we have found an alternative shortest path to v
            }
            return;
        }
        //decreasing this edge weight has decreased the shortest path length to rv

        assert(dist[ru] + weight < dist[rv]);
        dist[rv] = dist[ru] + weight;


        if(!edgeInShortestPathGraph[edgeID]){
            edgeInShortestPathGraph[edgeID] = true;
            assert(delta[rv] >= 0);
            delta[rv]++;
        }else{
            //this edge was _already_ the shortest path to the ru.
            assert(delta[rv] > 0);
        }
        q_batch.update(rv);

        //maintain delta invariant
        for(int j = 0; j < g.nIncoming(rv); j++){
            auto& e = g.incoming(rv, j);
            int adjID = e.id;

            if(g.edgeEnabled(adjID)){
                if(edgeInShortestPathGraph[adjID]){
                    assert(g.getEdge(adjID).to == rv);
                    int v = g.getEdge(adjID).from;
                    Weight& w = weights[adjID]; //assume a weight of one for now
                    Weight alt = dist[v] + w;
                    assert(w > 0);
                    if(dist[rv] < alt || alt >= INF){
                        edgeInShortestPathGraph[adjID] = false;
                        delta[rv]--;
                    }
                }
            }

        }
        assert(delta[rv] >= 1);
        assert(edgeInShortestPathGraph[edgeID]);
    }

    //Called if the weight of an edge is increased
    void IncreaseWeight(int edgeID){

        assert(g.edgeEnabled(edgeID)); //the edge must be enabled, in order to have its weight decreased
        //First, check if this edge is actually in the shortest path graph. If it isn't, then increasing its weight has no effect
        if(!edgeInShortestPathGraph[edgeID]){
            return;
        }

        int ru = g.getEdge(edgeID).from;
        int rv = g.getEdge(edgeID).to;
        assert(weights[edgeID] > 0);
        assert(delta[rv] > 0);
        delta[rv]--;
        edgeInShortestPathGraph[edgeID] = false;                        //remove this edge from the shortest path graph
        if(delta[rv] > 0){
            //there was another path to rv of the same previous length,
            //so increasing this edge doesn't effect the total path length to rv,
            //but does remove this edge from the shortest path graph
            return; //the shortest path hasn't changed in length, because there was an alternate route of the same length to this node.
        }
        //we increased this edge weight.
        //this edge MAY or MAY NOT be in the shortest path after that change.

        if(!in_changeset[rv]){
            changeset.push_back(rv);
            in_changeset[rv] = true;
        }

    }

    //Called if an edge is removed
    void RemoveEdge(int edgeID){

        assert(!g.edgeEnabled(edgeID));
        //First, check if this edge is actually in the shortest path graph
        if(!edgeInShortestPathGraph[edgeID])
            return;
        edgeInShortestPathGraph[edgeID] = false;                        //remove this edge from the shortest path graph

        int ru = g.getEdge(edgeID).from;
        int rv = g.getEdge(edgeID).to;
        assert(delta[rv] > 0);
        delta[rv]--;
        if(delta[rv] > 0)
            return; //the shortest path hasn't changed in length, because there was an alternate route of the same length to this node.


        if(!in_changeset[rv]){
            changeset.push_back(rv);
            in_changeset[rv] = true;
        }


    }

    int64_t num_updates = 0;

    int numUpdates() const override{
        return num_updates;
    }

    void update() override{

        if(g.outfile()){
            fprintf(g.outfile(), "r %d %d %d %d %d\n", getSource(), last_modification, g.getCurrentHistory(),
                    g.changed(), g.historySize());
        }

        stats_all_updates++;
        if(last_modification > 0 && g.getCurrentHistory() == last_modification)
            return;
        stats_updates++;
        if(last_modification <= 0 || g.changed()){
            stats_resets++;
            Weight oldInf = INF;

            INF = 1;
            has_zero_weights = false;
            for(Weight& w : weights){
                if(w <= 0){
                    //Note: in the future, we could implement the DFMN algorithm (Maintaining Shortest Paths in Digraphs with Arbitrary Arc Weights: An Experimental Study), which does support negative length weights, but is slower than RR.
                    //throw std::invalid_argument("Ramalingham-Reps doesn't support zero-weight edges (select a different distance algorithm, such as dijkstra)");
                    //for the moment: the _first_ time <0 weights are detected, simply fallback on dijkstra's, permanently.

                    has_zero_weights = true;
                }
                INF += w;
            }
            if(INF < (INF * 100)){
                INF = INF * 100;
            }
            if(INF < oldInf){
                INF = oldInf;
            }
            if(INF != oldInf){
                for(int i = 0; i < dist.size(); i++){
                    if(dist[i] == oldInf){
                        dist[i] = INF;
                    }
                }
            }
            while(local_weights.size() < weights.size()){
                local_weights.push_back(weights[local_weights.size()]);
            }
            dist.resize(g.nodes(), INF);
            dist[getSource()] = 0;
            delta.resize(g.nodes());
            node_changed.resize(g.nodes(), true);
            in_changeset.resize(g.nodes(), false);
            for(int i = 0; i < g.nodes(); i++){
                if((dist[i] >= INF && reportPolarity <= 0) || (dist[i] < INF && reportPolarity >= 0)){
                    node_changed[i] = true;
                    changed.push_back(i);                            //On the first round, report status of all nodes.
                }
            }
        }else{
            dbg_delta_lite();
        }
        edgeInShortestPathGraph.resize(g.nEdgeIDs());

        if(has_zero_weights){
            if(!ever_warned_about_zero_weights){
                ever_warned_about_zero_weights = true;
                fprintf(stderr,
                        "Warning: Ramalingham-Reps doesn't support zero-weight edges; falling back on Dijkstra's (which is much slower)\n");
            }
            dijkstras.update();
        }else{

            if(last_history_clear != g.nHistoryClears()){
                if(!g.changed() && last_history_clear >= 0 && last_history_clear == g.nHistoryClears() - 1 &&
                   history_qhead == g.getPreviousHistorySize()){
                    //no information was lost in the history clear
                    history_qhead = 0;
                    last_history_clear = g.nHistoryClears();
                }else{
                    history_qhead = g.historySize();
                    last_history_clear = g.nHistoryClears();
                    for(int edgeid = 0; edgeid < g.edges(); edgeid++){
                        if(g.edgeEnabled(edgeid)){
                            if(weights[edgeid] < local_weights[edgeid]){
                                local_weights[edgeid] = weights[edgeid];
                                if(weights[edgeid] == 0){
                                    assert(!has_zero_weights);
                                    if(!ever_warned_about_zero_weights){
                                        ever_warned_about_zero_weights = true;
                                        fprintf(stderr,
                                                "Warning: Ramalingham-Reps doesn't support zero-weight edges; falling back on Dijkstra's (which is much slower)\n");
                                    }
                                    has_zero_weights = true;
                                    dijkstras.update();
                                    break;
                                }
                                DecreaseWeight(edgeid);
                            }else if(weights[edgeid] > local_weights[edgeid]){
                                local_weights[edgeid] = weights[edgeid];
                                IncreaseWeight(edgeid);
                            }else{
                                AddEdge(edgeid);
                            }
                        }else{
                            local_weights[edgeid] = weights[edgeid];
                            RemoveEdge(edgeid);
                        }
                    }
                }
            }
            for(int i = history_qhead; i < g.historySize(); i++){
                int edgeid = g.getChange(i).id;

                if(g.getChange(i).weight_increase || g.getChange(i).weight_decrease){
                    if(g.getChange(i).weight_decrease && g.edgeEnabled(edgeid)){
                        local_weights[edgeid] = weights[edgeid];
                        if(weights[edgeid] == 0){
                            assert(!has_zero_weights);
                            if(!ever_warned_about_zero_weights){
                                ever_warned_about_zero_weights = true;
                                fprintf(stderr,
                                        "Warning: Ramalingham-Reps doesn't support zero-weight edges; falling back on Dijkstra's (which is much slower)\n");
                            }
                            has_zero_weights = true;
                            dijkstras.update();
                            break;
                        }
                        DecreaseWeight(
                                edgeid);//need to run this EVEN IF local weights==weights, because the DecreaseWeight
                        //code runs on weights, not local weights, and so it might be out of sync with local weights
                    }else if(g.getChange(i).weight_increase && g.edgeEnabled(edgeid)){
                        local_weights[edgeid] = weights[edgeid];
                        IncreaseWeight(
                                edgeid);//need to run this EVEN IF local weights==weights, because the DecreaseWeight
                        //code runs on weights, not local weights, and so it might be out of sync with local weights
                    }
                }else if(g.getChange(i).addition && g.edgeEnabled(edgeid)){
                    assert(!g.getChange(i).weight_increase);
                    assert(!g.getChange(i).weight_decrease);
                    local_weights[edgeid] = weights[edgeid];
                    AddEdge(edgeid);
                }else if(!g.getChange(i).addition && !g.edgeEnabled(edgeid)){
                    assert(!g.getChange(i).weight_increase);
                    assert(!g.getChange(i).weight_decrease);
                    local_weights[edgeid] = weights[edgeid];
                    RemoveEdge(edgeid);
                }
            }
        }
        if(!has_zero_weights){
            processChanged();
            processQ();

            dbg_delta_lite();
        }
        //for(int i = 0;i<g.nodes();i++){
        //	int u=i;
        if(reportPolarity > -2 && !has_zero_weights){
            for(int u : changed){
                //int u = changed[i];
                node_changed[u] = false;
                //CANNOT clear the change flag here, because if we backtrack and then immediately re-propagate an edge before calling theoryPropagate, the change to this node may be missed.
                if(reportPolarity <= 0 && dist[u] >= INF){
                    status.setReachable(u, false);
                    status.setMininumDistance(u, dist[u] < INF, dist[u]);
                }else if(reportPolarity >= 0 && dist[u] < INF){
                    status.setReachable(u, true);
                    status.setMininumDistance(u, dist[u] < INF, dist[u]);
                }
            }
            changed.clear();
        }
        assert(dbg_uptodate());

        num_updates++;
        last_modification = g.getCurrentHistory();
        last_deletion = g.nDeletions();
        last_addition = g.nAdditions();
        g.updateAlgorithmHistory(this, alg_id, history_qhead);
        history_qhead = g.historySize();
        last_history_clear = g.nHistoryClears();

#ifdef DEBUG_RAMAL
        for(int edgeID = 0;edgeID<weights.size();edgeID++){
          assert(has_zero_weights || (weights[edgeID]==local_weights[edgeID]) || !g.edgeEnabled(edgeID) );
      }
#endif
    }

    void processChanged(){
        //find all effected nodes whose shortest path lengths may now be increased (or that may have become unreachable)
        for(int i = 0; i < changeset.size(); i++){
            int u = changeset[i];
            dist[u] = INF;
            for(int i = 0; i < g.nIncident(u); i++){
                auto& e = g.incident(u, i);
                int adjID = e.id;
                if(g.edgeEnabled(adjID)){
                    if(edgeInShortestPathGraph[adjID]){
                        edgeInShortestPathGraph[adjID] = false;
                        assert(g.getEdge(adjID).from == u);
                        int s = g.getEdge(adjID).to;
                        assert(delta[s] > 0);
                        delta[s]--;
                        if(delta[s] == 0){
                            if(!in_changeset[s]){
                                changeset.push_back(s);
                                in_changeset[s] = true;
                            }
                        }
                    }
                }
            }
        }

        for(int i = 0; i < changeset.size(); i++){
            int u = changeset[i];
            assert(in_changeset[u]);
            in_changeset[u] = false;
            assert(dist[u] == INF);
            int shortest_edge = -1;
            for(int j = 0; j < g.nIncoming(u); j++){
                auto& e = g.incoming(u, j);
                int adjID = e.id;

                if(g.edgeEnabled(adjID)){
                    assert(g.getEdge(adjID).to == u);
                    int v = g.getEdge(adjID).from;
                    Weight& w = weights[adjID]; //assume a weight of one for now
                    Weight alt = dist[v] + w;

                    if(dist[u] > alt){
                        assert(alt < INF);
                        dist[u] = alt;
                        shortest_edge = adjID;
                    }
                }

            }


            if(dist[u] != INF){
                assert (shortest_edge >= 0);

                if(!edgeInShortestPathGraph[shortest_edge]){
                    edgeInShortestPathGraph[shortest_edge] = true;
                    delta[u]++;
                }

                //q.insert(u);
                //dbg_Q_add(q,u);
                q_batch.update(u);

                if(!reportDistance && reportPolarity >= 0){
                    if(!node_changed[u]){
                        node_changed[u] = true;
                        changed.push_back(u);
                    }
                }
            }else if(reportPolarity <= 0){
                //have to mark this change even if we are reporting distanec, as u has not been added to the queue.
                if(!node_changed[u]){
                    node_changed[u] = true;
                    changed.push_back(u);
                }
            }
        }
        changeset.clear();
#ifndef NDEBUG
        for(int i = 0; i < in_changeset.size(); i++){
            assert(!in_changeset[i]);
        }
#endif
    }

    void processQ(){
        while(q_batch.size() > 0){
            int u = q_batch.removeMin();
            if(reportDistance){
                if(dist[u] != INF){
                    if(reportPolarity >= 0){
                        if(!node_changed[u]){
                            node_changed[u] = true;
                            changed.push_back(u);
                        }
                    }
                }else if(reportPolarity <= 0){
                    if(!node_changed[u]){
                        node_changed[u] = true;
                        changed.push_back(u);
                    }
                }
            }
            delta[u] = 0;
            for(int i = 0; i < g.nIncoming(u); i++){
                auto& e = g.incoming(u, i);
                int adjID = e.id;
                if(g.edgeEnabled(adjID)){
                    assert(g.getEdge(adjID).to == u);
                    int v = g.getEdge(adjID).from;
                    Weight& dv = dist[v];
                    Weight& du = dist[u];
                    //Weight & w = ;
                    Weight alt = dist[v] + weights[adjID];
                    if(alt > INF){
                        alt = INF;
                    }
                    if(alt < INF && dist[u] == alt){
                        edgeInShortestPathGraph[adjID] = true;
                        delta[u]++;
                    }else if(dist[u] < alt || alt == INF){
                        edgeInShortestPathGraph[adjID] = false;
                    }else if(dist[u] > alt){
                        //assert(false);//not clear if this should ever happen, because it means that dist[u] is incorrect at this point, but that might still be corrected below
                    }else{
                        assert(!edgeInShortestPathGraph[adjID]);
                    }
                }else{
                    edgeInShortestPathGraph[adjID] = false;    //need to add this, because we may have disabled multiple edges at once.
                }
            }

            for(int i = 0; i < g.nIncident(u); i++){
                auto& e = g.incident(u, i);
                int adjID = e.id;
                if(g.edgeEnabled(adjID)){
                    assert(g.getEdge(adjID).from == u);
                    int s = g.getEdge(adjID).to;
                    Weight& w = weights[adjID];
                    Weight alt = dist[u] + w;
                    if(alt > INF){
                        alt = INF;
                    }
                    if(dist[s] > alt){
                        if(reportPolarity >= 0 && dist[s] >= 0){
                            //This check is needed (in addition to the above), because even if we are NOT reporting distances, it is possible for a node that was previously not reachable
                            //to become reachable here. This is ONLY possible because we are batching multiple edge incs/decs at once (otherwise it would be impossible for removing an edge to decrease the distance to a node).
                            if(!node_changed[s]){
                                node_changed[s] = true;
                                changed.push_back(s);
                            }
                        }

                        dist[s] = alt;
                        q_batch.update(s);
                    }else if(alt < INF && dist[s] == alt){
                        if(!edgeInShortestPathGraph[adjID]){
                            edgeInShortestPathGraph[adjID] = true;
                            delta[s]++;
                        }
                    }else if(alt == INF && edgeInShortestPathGraph[adjID]){
                        edgeInShortestPathGraph[adjID] = false;
                        delta[s]--;
                    }else if(dist[s] < alt){
                        //assert(!edgeInShortestPathGraph[adjID] || ( dist[s] >=dist[u] && q_batch.inHeap(s)));
                    }
                }else{
                    assert(!edgeInShortestPathGraph[adjID]);
                }
            }


        }
    }


    void printStats() override{
        printf("Updates: %" PRId64 " (+%" PRId64 " skipped), %" PRId64 " restarts\n", stats_updates,
               stats_all_updates - stats_updates, stats_resets);
    }

    void updateHistory() override{
        update();
    }

    bool dbg_path(int to){
#ifdef DEBUG_RAMAL
        /*	assert(connected(to));
         if(to == source){
         return true;
         }
         int p = previous(to);

         if(p<0){
         return false;
         }
         if(p==to){
         return false;
         }

         return dbg_path(p);*/

#endif
        return true;
    }

    bool dbg_manual_uptodate() override{

        if(last_modification < 0)
            return true;
        update();
        dbg_delta();
        Dijkstra<Weight, Graph> d(source, g);

        for(int i = 0; i < g.nodes(); i++){
            Weight dis = dist[i];
            bool c = i < dist.size() && dist[i] < INF;
            if(!c)
                dis = this->unreachable();

            Weight dbgdist = d.distance(i);

            if(dis != dbgdist){
                assert(false);
                throw std::logic_error("Internal error in Ramal Reps");
            }
        }
//#endif

        return true;
    }

    bool dbg_uptodate(){
#ifdef DEBUG_RAMAL2
        if(last_modification<0 || has_zero_weights)
         return true;
         dbg_delta();
         Dijkstra<Weight> d(source,g);

         for(int i = 0;i<g.nodes();i++){
         Weight dis = dist[i];
         bool c = i<dist.size() && dist[i]<INF;
         if(!c)
         dis = this->unreachable();

         Weight dbgdist = d.distance(i);

         if(dis!=dbgdist){
         assert(false);
         throw std::logic_error("Internal error in Ramal Reps");
         }
         }
//#endif
#endif
        return true;
    }

    bool connected_unsafe(int t) override{
        //dbg_uptodate();
        if(has_zero_weights){
            return dijkstras.connected_unsafe(t);
        }
        return t < dist.size() && dist[t] < INF;
    }

    bool connected_unchecked(int t) override{
        assert(last_modification == g.getCurrentHistory());
        return connected_unsafe(t);
    }

    bool connected(int t) override{

        update();

        assert(dbg_uptodate());
        if(has_zero_weights){
            return dijkstras.connected(t);
        }
        return dist[t] < INF;
    }

    Weight& distance(int t) override{

        update();
        if(has_zero_weights){
            return dijkstras.distance(t);
        }
        if(connected_unsafe(t))
            return dist[t];
        else
            return this->unreachable();
    }

    Weight& distance_unsafe(int t) override{
        if(has_zero_weights){
            return dijkstras.distance_unsafe(t);
        }
        if(connected_unsafe(t))
            return dist[t];
        else
            return this->unreachable();
    }

    int incomingEdge(int t) override{
        if(has_zero_weights){
            return dijkstras.incomingEdge(t);
        }
        if(!connected_unsafe(t)){
            return -1;
        }
        if(t == source)
            return -1;

        assert(dist[t] >= 0);
        assert(dist[t] != INF);
        assert(delta[t] > 0);
        int prev = -1;
        int prev_edgeID = -1;

        //it should be possible to maintain an explicit list of all the edges in the shortest path tree,
        //or perhaps at least one such edge for each node, and avoid this search, at the cost of more storage and slightly more expensive
        //edge updates
        for(int i = 0; i < g.nIncoming(t); i++){

            int edgeID = g.incoming(t, i).id;
            if(edgeInShortestPathGraph[edgeID]){
                assert(g.edgeEnabled(edgeID));
                int from = g.incoming(t, i).node;
                assert(connected_unsafe(from));
                assert(dist[from] >= 0);
                assert(dist[from] != INF);
                prev = from;
                prev_edgeID = edgeID;
                break;
            }
        }
        assert(prev != -1);
        //assert(min_prev_dist<dist[t]);
        assert(prev_edgeID != -1);
        return prev_edgeID;
    }

    int previous(int t) override{
        if(has_zero_weights){
            return dijkstras.previous(t);
        }
        if(!connected_unsafe(t)){
            return -1;
        }
        if(t == source)
            return -1;

        assert(dist[t] >= 0);
        assert(dist[t] != INF);
        assert(delta[t] > 0);
        int prev = -1;

        for(int i = 0; i < g.nIncoming(t); i++){

            int edgeID = g.incoming(t, i).id;
            if(edgeInShortestPathGraph[edgeID]){
                assert(g.edgeEnabled(edgeID));
                int from = g.incoming(t, i).node;
                assert(connected_unsafe(from));
                assert(dist[from] >= 0);
                assert(dist[from] != INF);
                prev = from;
                break;
            }
        }
        assert(prev != -1);


        return prev;
    }

    //Select a randomly chosen, but still shortest path, previous node or edge (falling back on selecting a deterministic edge if this functionality is not supported)
    int randomPrevious(int t, double& seed) override{
        if(has_zero_weights){
            return dijkstras.previous(t);
        }
        if(!connected_unsafe(t)){
            return -1;
        }
        if(t == source)
            return -1;

        assert(dist[t] >= 0);
        assert(dist[t] != INF);
        assert(delta[t] > 0);
        int prev = -1;

        int start = dgl::alg::irand(seed, g.nIncoming(t));
        assert(start >= 0);
        assert(start < g.nIncoming(t));
        for(int j = 0; j < g.nIncoming(t); j++){
            int i = (j + start) % g.nIncoming(t);
            assert(i >= 0);
            assert(i < g.nIncoming(t));
            int edgeID = g.incoming(t, i).id;
            if(edgeInShortestPathGraph[edgeID]){
                assert(g.edgeEnabled(edgeID));
                int from = g.incoming(t, i).node;
                assert(connected_unsafe(from));
                assert(dist[from] >= 0);
                assert(dist[from] != INF);
                prev = from;
                break;
            }
        }
        assert(prev != -1);


        return prev;
    }

    int randomIncomingEdge(int t, double& seed) override{
        if(has_zero_weights){
            return dijkstras.incomingEdge(t);
        }
        if(!connected_unsafe(t)){
            return -1;
        }
        if(t == source)
            return -1;

        assert(dist[t] >= 0);
        assert(dist[t] != INF);
        assert(delta[t] > 0);
        int prev = -1;
        int prev_edgeID = -1;

        //it should be possible to maintain an explicit list of all the edges in the shortest path tree,
        //or perhaps at least one such edge for each node, and avoid this search, at the cost of more storage and slightly more expensive
        //edge updates
        int start = dgl::alg::irand(seed, g.nIncoming(t));
        assert(start >= 0);
        assert(start < g.nIncoming(t));
        for(int j = 0; j < g.nIncoming(t); j++){
            int i = (j + start) % g.nIncoming(t);
            assert(i >= 0);
            assert(i < g.nIncoming(t));

            int edgeID = g.incoming(t, i).id;
            if(edgeInShortestPathGraph[edgeID]){
                assert(g.edgeEnabled(edgeID));
                int from = g.incoming(t, i).node;
                assert(connected_unsafe(from));
                assert(dist[from] >= 0);
                assert(dist[from] != INF);
                prev = from;
                prev_edgeID = edgeID;
                break;
            }
        }
        assert(prev != -1);

        assert(prev_edgeID != -1);
        return prev_edgeID;
    }

};

template<typename Weight, typename Graph = DynamicGraph<Weight>, class Status= typename Distance<Weight>::NullStatus>
class UnweightedRamalRepsBatchedUnified : public Distance<int>, public DynamicGraphAlgorithm {
public:
    Graph& g;
    Status& status;
    int reportPolarity;
    bool reportDistance;
    int last_modification = -1;
    int last_addition = -1;
    int last_deletion = -1;
    int history_qhead = 0;
    uint64_t stats_updates = 0;
    uint64_t stats_resets = 0;
    uint64_t stats_all_updates = 0;
    int last_history_clear = 0;

    int source;
    int INF = -1;

    int maxDistance = -1;

    std::vector<int> old_dist;
    std::vector<int> changed;
    std::vector<bool> node_changed;
    std::vector<int> dist;
    std::vector<int> prev;

    struct DistCmp {
        std::vector<int>& _dist;

        bool operator()(int a, int b) const{
            return _dist[a] < _dist[b];
        }

        DistCmp(std::vector<int>& d) :
                _dist(d){
        };
    };

    //Heap<DistCmp> q;
    //std::vector<bool> in_queue_inc;
    std::vector<bool> in_queue;
    std::vector<bool> in_queue2;
    //std::vector<bool> in_queue_inc2;

    std::vector<int> q_batched;
    std::vector<int> q_batched2;
    std::vector<int> edgeInShortestPathGraph; //edgeInShortestPathGraph(v,u) is true if there exists any shortest path from s to u that includes the edge (v,u)
    std::vector<int> delta; //delta records, for each node u, the number of incoming edges (v ,u) for which edgeInShortestPathGraph (v,u) is true
    std::vector<int> changeset;
    std::vector<bool> in_changeset;
    int alg_id = -1;
public:

    int64_t stats_full_updates = 0;
    int64_t stats_fast_updates = 0;
    int64_t stats_fast_failed_updates = 0;
    int64_t stats_skip_deletes = 0;
    int64_t stats_skipped_updates = 0;
    int64_t stats_num_skipable_deletions = 0;
    double mod_percentage = 0;

    double stats_full_update_time = 0;
    double stats_fast_update_time = 0;

    UnweightedRamalRepsBatchedUnified(int s, Graph& graph, Status& status, int reportPolarity = 0,
                                      bool reportDistance = true) :
            g(graph), status(status), reportPolarity(reportPolarity), reportDistance(reportDistance), last_modification(
            -1), last_addition(-1), last_deletion(-1), history_qhead(0), last_history_clear(0), source(s), INF(
            -1){
        maxDistance = -1;
        mod_percentage = 0.2;
        alg_id = g.addDynamicAlgorithm(this);
    }

    UnweightedRamalRepsBatchedUnified(int s, Graph& graph, int reportPolarity = 0,
                                      bool reportDistance = true) :
            g(graph), status(Distance<Weight>::nullStatus), reportPolarity(reportPolarity),
            reportDistance(reportDistance), last_modification(
            -1), last_addition(-1), last_deletion(-1), history_qhead(0), last_history_clear(0), source(s), INF(
            -1){
        maxDistance = -1;
        mod_percentage = 0.2;
        alg_id = g.addDynamicAlgorithm(this);
    }

    std::string getName() override{
        return "UnweightedRamalRepsBatchedUnified(" + std::to_string(getSource()) + ")";
    }

    void printStats() override{
        printf("Updates: %" PRId64 " (+%" PRId64 " skipped), %" PRId64 " restarts\n", stats_updates,
               stats_all_updates - stats_updates, stats_resets);
    }

    //Dijkstra(const Dijkstra& d):g(d.g), last_modification(-1),last_addition(-1),last_deletion(-1),history_qhead(0),last_history_clear(0),source(d.source),INF(0),q(DistCmp(dist)),stats_full_updates(0),stats_fast_updates(0),stats_skip_deletes(0),stats_skipped_updates(0),stats_full_update_time(0),stats_fast_update_time(0){marked=false;};
    void setMaxDistance(int& _maxDistance) override{
        if(_maxDistance != maxDistance){
            last_modification = -1;        //force the next update to recompute from scratch
            if(_maxDistance < 0){
                maxDistance = INF;
            }else
                maxDistance = _maxDistance;
        }
    }

    void setSource(int s) override{
        source = s;
        last_modification = -1;
        last_addition = -1;
        last_deletion = -1;
    }

    int getSource() override{
        return source;
    }

    std::vector<int>& getChanged(){
        return changed;
    }

    void clearChanged(){
        changed.clear();
    }

    void drawFull(){

    }

    void dbg_delta(){
#ifdef DEBUG_RAMAL
        //g.drawFull();
        dbg_delta_lite();
        assert(delta.size() == g.nodes());

        for (int i = 0; i < g.nEdgeIDs(); i++) {
            if (!g.edgeEnabled(i)) {
                assert(!edgeInShortestPathGraph[i]);
                if (edgeInShortestPathGraph[i]) {
                    throw std::runtime_error("");
                }
            }
        }

        std::vector<int> dbg_delta;
        std::vector<int> dbg_dist;
        dbg_dist.resize(g.nodes(), INF);
        dbg_delta.resize(g.nodes());
        dbg_dist[getSource()] = 0;

        struct DistCmp {
            std::vector<int> & _dist;
            bool operator()(int a, int b) const {
                return _dist[a] < _dist[b];
            }
            DistCmp(std::vector<int> & d) :
                    _dist(d) {
            }
            ;
        };
        alg::Heap<DistCmp> q(dbg_dist);

        q.insert(getSource());

        while (q.size()) {
            int u = q.removeMin();
            if (dbg_dist[u] == INF)
                break;
            dbg_delta[u] = 0;

            for (int i = 0; i < g.nIncoming(u); i++) {
                if (!g.edgeEnabled(g.incoming(u, i).id))
                    continue;

                int edgeID = g.incoming(u, i).id;
                int v = g.getEdge(edgeID).from;
                int alt = dbg_dist[v] + 1;
                if (maxDistance >= 0 && alt > maxDistance)
                    alt = INF;
                assert(alt >= dbg_dist[u]);
                /*	if (alt==dbg_dist[u]){
                 dbg_delta[u]++;
                 }*/
            }

            for (int i = 0; i < g.nIncident(u); i++) {
                if (!g.edgeEnabled(g.incident(u, i).id))
                    continue;

                int edgeID = g.incident(u, i).id;
                int v = g.getEdge(edgeID).to;
                int alt = dbg_dist[u] + 1;
                if (maxDistance >= 0 && alt > maxDistance)
                    alt = INF;
                if (alt < dbg_dist[v]) {

                    dbg_dist[v] = alt;

                    if (!q.inHeap(v))
                        q.insert(v);
                    else
                        q.decrease(v);
                }
            }
        }

        for (int u = 0; u < g.nodes(); u++) {
            int d = dist[u];

            int db = dbg_dist[u];
            assert(dbg_dist[u] == dist[u]);

            for (int i = 0; i < g.nIncoming(u); i++) {
                if (!g.edgeEnabled(g.incoming(u, i).id))
                    continue;

                int edgeID = g.incoming(u, i).id;
                int v = g.getEdge(edgeID).from;
                int alt = dbg_dist[v] + 1;
                int du = dbg_dist[u];
                if (maxDistance >= 0 && alt > maxDistance)
                    alt = INF;
                assert(alt >= dbg_dist[u]);
                if (alt == dbg_dist[u] && alt < INF) {
                    dbg_delta[u]++;
                    assert(edgeInShortestPathGraph[edgeID]);
                } else if (alt < INF) {
                    assert(!edgeInShortestPathGraph[edgeID]);
                }
            }
        }

        for (int u = 0; u < g.nodes(); u++) {
            int d = dist[u];
            int d_expect = dbg_dist[u];
            assert(d == dbg_dist[u]);
        }
        for (int u = 0; u < g.nodes(); u++) {
            int du = dist[u];
            if (dist[u] < INF) {
                int d = delta[u];
                int d_expect = dbg_delta[u];
                assert(d == dbg_delta[u]);
            }
        }
        dbg_delta_lite();
#endif
    }

    void AddEdge(int edgeID){
        dbg_delta_lite();
        assert(g.edgeEnabled(edgeID));
        if(edgeInShortestPathGraph[edgeID])
            return;
        int ru = g.getEdge(edgeID).from;
        int rv = g.getEdge(edgeID).to;

        int rdv = dist[rv];
        int rdu = dist[ru];

        int weight = 1;
        int altw = dist[ru] + weight;
        if(altw > maxDistance)
            altw = INF;
        if(dist[rv] < altw)
            return;
        else if(dist[rv] == altw && dist[rv] < INF){
            assert(!edgeInShortestPathGraph[edgeID]);
            edgeInShortestPathGraph[edgeID] = true;
            delta[rv]++;        //we have found an alternative shortest path to v
            return;
        }else if(altw == INF){
            return;        //don't do anything
        }
        assert(altw < INF);
        edgeInShortestPathGraph[edgeID] = true;
        assert(delta[rv] >= 0);
        delta[rv]++;
        dist[rv] = altw;

        if(!in_queue[rv]){
            q_batched.push_back(rv);
            in_queue[rv] = true;
        }

        //maintain delta invariant
        for(int j = 0; j < g.nIncoming(rv); j++){
            auto& e = g.incoming(rv, j);
            int adjID = e.id;
            if(g.edgeEnabled(adjID)){
                if(edgeInShortestPathGraph[adjID]){
                    assert(g.getEdge(adjID).to == rv);
                    int v = g.getEdge(adjID).from;
                    Weight w = 1;
                    Weight alt = dist[v] + w;
                    if(dist[rv] < alt || alt >= INF){
                        edgeInShortestPathGraph[adjID] = false;
                        delta[rv]--;
                    }
                }
            }

        }
        assert(delta[rv] >= 1);
        assert(edgeInShortestPathGraph[edgeID]);
    }

    void dbg_not_seen_q(std::vector<int>& q, int u, int from){
#ifdef DEBUG_RAMAL
        bool found = false;
        for (int i = from; i < q.size(); i++) {
            if (q[i] == u) {
                found = true;
                break;
            }
        }
        assert(found);
#endif
    }

    void dbg_Q_add(std::vector<int>& q, int u){
#ifdef DEBUG_RAMAL
        //assert(!in_queue[u]);
        for (int v : q) {
            assert(u != v);
            assert(dist[v] <= dist[u]);
            //assert(in_queue[v]);
        }

#endif
    }

    void dbg_Q_order(std::vector<int>& _q){
#ifdef DEBUG_RAMAL

        for (int i = 1; i < _q.size(); i++) {
            int v = _q[i];
            int u = _q[i - 1];

            if (&_q == &q_batched) {
                assert(in_queue[v]);
                assert(in_queue[u]);
                if (!(in_queue2[u] || in_queue2[v])) {

                    assert(dist[u] <= dist[v]);
                }
            }else if (&_q == &q_batched2) {
                assert(in_queue2[v]);
                assert(in_queue2[u]);
                assert(dist[u] <= dist[v]);
            }
        }

#endif
    }

    void dbg_delta_lite(){
#ifdef DEBUG_RAMAL
        /*	for (int u = 0; u < g.nodes(); u++) {
            int del = delta[u];
            int d = dist[u];
            int num_in = 0;
            for (int i = 0; i < g.nIncoming(u); i++) {
                auto & e = g.incoming(u, i);
                int adjID = e.id;
                int from = g.getEdge(adjID).from;

                int dfrom = dist[from];
                if (edgeInShortestPathGraph[adjID]) {

                    assert(dist[u] < INF);
                    num_in++;
                }
            }
            assert(del == num_in);
        }*/
#endif

    }

    void RemoveEdge(int edgeID){
        dbg_delta_lite();
        assert(!g.edgeEnabled(edgeID));
        //First, check if this edge is actually in the shortest path graph
        if(!edgeInShortestPathGraph[edgeID])
            return;
        edgeInShortestPathGraph[edgeID] = false;                        //remove this edge from the shortest path graph

        int ru = g.getEdge(edgeID).from;
        int rv = g.getEdge(edgeID).to;

        assert(delta[rv] > 0);
        delta[rv]--;
        if(delta[rv] > 0)
            return; //the shortest path hasn't changed in length, because there was an alternate route of the same length to this node.


        if(!in_changeset[rv]){
            changeset.push_back(rv);
            in_changeset[rv] = true;
        }


    }

    int64_t num_updates = 0;

    int numUpdates() const override{
        return num_updates;
    }

    void update() override{

        if(g.outfile()){
            fprintf(g.outfile(), "r %d %d %d %d %d\n", getSource(), last_modification, g.getCurrentHistory(),
                    g.changed(), g.historySize());
        }
        stats_all_updates++;
        if(last_modification > 0 && g.getCurrentHistory() == last_modification){
            return;
        }
        stats_updates++;
        if(last_modification <= 0 || g.changed()){//Note for the future: there is probably room to improve this further.
            stats_resets++;
            stats_full_updates++;
            int oldInf = INF;
            INF = g.nodes() + 1;
            if(INF < oldInf){
                INF = oldInf;
            }
            if(INF != oldInf){
                for(int i = 0; i < dist.size(); i++){
                    if(dist[i] == oldInf){
                        dist[i] = INF;
                    }
                }
            }
            if(maxDistance < 0 || maxDistance >= oldInf){
                maxDistance = INF;
            }
            last_history_clear = -1;
            dist.resize(g.nodes(), INF);
            dist[getSource()] = 0;
            delta.resize(g.nodes());
            edgeInShortestPathGraph.resize(g.nEdgeIDs());
            node_changed.resize(g.nodes());
            changed.clear();

            in_queue.resize(g.nodes());
            in_queue.resize(g.nodes());
            in_queue2.resize(g.nodes());
            in_queue2.resize(g.nodes());
            in_changeset.resize(g.nodes(), false);
            for(int i = 0; i < g.nodes(); i++){
                if((dist[i] >= INF && reportPolarity <= 0) || (dist[i] < INF && reportPolarity >= 0)){
                    node_changed[i] = true;
                    changed.push_back(i);                            //On the first round, report status of all nodes.
                }
            }
        }

        if(last_history_clear != g.nHistoryClears()){
            if(!g.changed() && last_history_clear >= 0 && last_history_clear == g.nHistoryClears() - 1 &&
               history_qhead == g.getPreviousHistorySize()){
                //no information was lost in the history clear
                history_qhead = 0;
                last_history_clear = g.nHistoryClears();
            }else{
                history_qhead = g.historySize();
                last_history_clear = g.nHistoryClears();
                for(int edgeid = 0; edgeid < g.edges(); edgeid++){
                    if(g.edgeEnabled(edgeid)){
                        AddEdge(edgeid);
                    }else{
                        RemoveEdge(edgeid);
                    }
                }
            }
        }

        for(int i = history_qhead; i < g.historySize(); i++){
            int edgeid = g.getChange(i).id;
            if(g.getChange(i).addition && g.edgeEnabled(edgeid)){
                AddEdge(edgeid);
            }else if(!g.getChange(i).addition && !g.edgeEnabled(edgeid)){
                RemoveEdge(edgeid);
            }
        }
        processChanged();
        processQ();

        dbg_delta_lite();
        //for(int i = 0;i<g.nodes();i++){
        for(int u : changed){
            //int u=i;
            //int u = changed[i];
            node_changed[u] = false;

            if(reportPolarity <= 0 && dist[u] >= INF){
                status.setReachable(u, false);
                status.setMininumDistance(u, dist[u] < INF, dist[u]);
            }else if(reportPolarity >= 0 && dist[u] < INF){
                status.setReachable(u, true);
                status.setMininumDistance(u, dist[u] < INF, dist[u]);
            }
        }
        changed.clear();
        //}
        dbg_delta();
        num_updates++;
        last_modification = g.getCurrentHistory();
        last_deletion = g.nDeletions();
        last_addition = g.nAdditions();

        history_qhead = g.historySize();
        g.updateAlgorithmHistory(this, alg_id, history_qhead);
        last_history_clear = g.nHistoryClears();
        assert(dbg_uptodate());

    }

    void processChanged(){
        //find all effected nodes whose shortest path lengths may now be increased (or that may have become unreachable)
        for(int i = 0; i < changeset.size(); i++){
            int u = changeset[i];
            assert(in_changeset[u]);
            dist[u] = INF;
            for(int i = 0; i < g.nIncident(u); i++){
                auto& e = g.incident(u, i);
                int adjID = e.id;
                if(g.edgeEnabled(adjID)){
                    if(edgeInShortestPathGraph[adjID]){
                        edgeInShortestPathGraph[adjID] = false;
                        assert(g.getEdge(adjID).from == u);
                        int s = g.getEdge(adjID).to;

                        assert(delta[s] > 0);
                        delta[s]--;
                        if(delta[s] == 0){
                            if(!in_changeset[s]){
                                in_changeset[s] = true;
                                changeset.push_back(s);
                            }
                        }
                    }
                }
            }
        }

        for(int i = 0; i < changeset.size(); i++){
            int u = changeset[i];
            assert(in_changeset[u]);
            in_changeset[u] = false;
            int shortest_edge = -1;
            assert(dist[u] == INF);
            //for(auto & e:g.inverted_adjacency[u]){
            for(int i = 0; i < g.nIncoming(u); i++){
                auto& e = g.incoming(u, i);
                int adjID = e.id;

                if(g.edgeEnabled(adjID)){
                    assert(g.getEdge(adjID).to == u);
                    int v = g.getEdge(adjID).from;
                    int w = 1; //assume a weight of one for now
                    int alt = dist[v] + w;
                    if(alt > maxDistance)
                        alt = INF;

                    if(dist[u] > alt){
                        assert(alt < INF);
                        dist[u] = alt;
                        shortest_edge = adjID;
                    }
                }

            }

            if(dist[u] < INF){
                assert (shortest_edge >= 0);

                if(!edgeInShortestPathGraph[shortest_edge]){
                    edgeInShortestPathGraph[shortest_edge] = true;
                    delta[u]++;
                }
                //q.insert(u);
                //dbg_Q_add(q,u);
                if(!in_queue[u]){
                    q_batched.push_back(u);
                    in_queue[u] = true;
                }
                if(reportPolarity >= 0){
                    if(!node_changed[u]){
                        node_changed[u] = true;
                        changed.push_back(u);
                    }
                }
            }else if(reportPolarity <= 0){
                //call this even if we are reporting distance, because u hasn't been placed in the queue!
                if(!node_changed[u]){
                    node_changed[u] = true;
                    changed.push_back(u);
                }
            }

        }
        changeset.clear();
#ifndef NDEBUG
        for(int i = 0; i < in_changeset.size(); i++){
            assert(!in_changeset[i]);
        }
#endif
    }

    void processQ(){
        //breaking this up into a 'q' that is really a vector, and a second q_batched2, also just a vector that will store the elements
        //added to it afterward
        //the idea here is to use vectors instead of heaps.
        //the values in the original q might actually have differing values, and so do require sorting.
        //this is why the line " else if (dist[q[i]] < dist[q_batched2[j]]) " checks whether the next smallest element is in q or q_batched2.
        //however, because the graph is unweighted, it is safe to visit the elements of q_batched2 in the order they are seen as well, storing them in a simple vector.

        assert(!q_batched2.size());
#ifdef DEBUG_RAMAL
        for(int n:q_batched){
            assert(in_queue[n]);
        }
        for(int i = 0;i<in_queue2.size();i++){
            assert(!in_queue2[i]);
        }
#endif
        std::sort(q_batched.begin(), q_batched.end(), DistCmp(dist));
        int i = 0, j = 0;
        while(i < q_batched.size() || j < q_batched2.size()){
            int u;
            if(i == q_batched.size()){
                assert(j < q_batched2.size());
                u = q_batched2[j++];
            }else if(j == q_batched2.size()){
                assert(i < q_batched.size());
                u = q_batched[i++];
                if(in_queue2[u]){
                    continue;
                }
            }else if(dist[q_batched[i]] < dist[q_batched2[j]]){
                u = q_batched[i++];
                if(in_queue2[u]){
                    continue;
                }
            }else{
                assert(dist[q_batched2[j]] <= dist[q_batched[i]]);
                u = q_batched2[j++];
            }

            dbg_Q_order(q_batched);
            dbg_Q_order(q_batched2);
            if(!node_changed[u]){
                node_changed[u] = true;
                changed.push_back(u);
            }
            if(dist[u] >= INF){
                int a = 1;
            }
            delta[u] = 0;
            //assert(dist[u] < INF);
            for(int k = 0; k < g.nIncoming(u); k++){
                auto& e = g.incoming(u, k);
                int adjID = e.id;
                if(g.edgeEnabled(adjID)){

                    assert(g.getEdge(adjID).to == u);
                    int v = g.getEdge(adjID).from;
                    int w = 1;        //assume a weight of one for now
                    int du = dist[u];
                    int dv = dist[v];
                    int alt = dist[v] + w;
                    if(alt > maxDistance)
                        alt = INF;
                    if(alt >= INF){
                        edgeInShortestPathGraph[adjID] = false;
                    }else if(dist[u] == alt){
                        assert(alt < INF);
                        edgeInShortestPathGraph[adjID] = true;
                        delta[u]++;
                    }else if(dist[u] < alt){
                        //This doesn't hold for us, because we are allowing multiple edges to be added at once.
                        //assert(dist[u]<(dist[v]+w));

                        edgeInShortestPathGraph[adjID] = false;
                    }else{
                        //don't do anything. This will get corrected in a future call to AddEdge.
                        //assert(false);

                    }
                }else{
                    edgeInShortestPathGraph[adjID] = false;    //need to add this, because we may have disabled multiple edges at once.
                }
            }

            for(int k = 0; k < g.nIncident(u); k++){
                auto& e = g.incident(u, k);
                int adjID = e.id;
                if(g.edgeEnabled(adjID)){
                    assert(g.getEdge(adjID).from == u);
                    int s = g.getEdge(adjID).to;
                    int w = 1;                            //assume a weight of one for now
                    int du = dist[u];
                    int ds = dist[s];
                    int alt = dist[u] + w;
                    if(alt > maxDistance)
                        alt = INF;
                    if(dist[s] > alt){
                        dist[s] = alt;
                        assert(alt < INF);
                        if(reportPolarity >= 0 && dist[s] >= 0){
                            //This check is needed (in addition to the above), because even if we are NOT reporting distances, it is possible for a node that was previously not reachable
                            //to become reachable here. This is ONLY possible because we are batching multiple edge incs/decs at once (otherwise it would be impossible for removing an edge to decrease the distance to a node).
                            if(!node_changed[s]){
                                node_changed[s] = true;
                                changed.push_back(s);
                            }
                        }
                        if(!in_queue2[s]){
                            dbg_Q_add(q_batched2, s);
                            q_batched2.push_back(s);
                            in_queue2[s] = true;
                        }
                    }else if(dist[s] == alt && !edgeInShortestPathGraph[adjID] && alt < INF){
                        assert(alt < INF);
                        edgeInShortestPathGraph[adjID] = true;
                        delta[s]++;
                    }
                }else{
                    assert(!edgeInShortestPathGraph[adjID]);
                }
            }
        }
        for(int j:q_batched){
            assert(in_queue[j]);
            in_queue[j] = false;
        }
        q_batched.clear();
        for(int j:q_batched2){
            assert(in_queue2[j]);
            in_queue2[j] = false;
        }
        q_batched2.clear();
#ifdef DEBUG_RAMAL

        for(int i = 0;i<in_queue2.size();i++){
            assert(!in_queue2[i]);
        }

        for(int i = 0;i<in_queue.size();i++){
            assert(!in_queue[i]);
        }

#endif
    }

    void updateHistory() override{
        update();
    }

    bool dbg_path(int to){
#ifdef DEBUG_RAMAL
        assert(connected(to));
        if(to == source) {
            return true;
        }
        int p = previous(to);

        if(p<0) {
            return false;
        }
        if(p==to) {
            return false;
        }

        return dbg_path(p);

#endif
        return true;
    }

    bool dbg_manual_uptodate() override{

        if(last_modification < 0)
            return true;
        update();
        dbg_delta();
        UnweightedDijkstra<Weight, Graph> d(source, g);

        for(int i = 0; i < g.nodes(); i++){
            Weight dis = dist[i];
            bool c = i < dist.size() && dist[i] < INF;
            if(!c)
                dis = this->unreachable();

            Weight dbgdist = d.distance(i);

            if(dis != dbgdist){
                assert(false);
                throw std::logic_error("Internal error in Ramal Reps");
            }
        }
//#endif

        return true;
    }

    bool dbg_uptodate(){
//#ifdef DEBUG_GRAPH
#ifdef DEBUG_RAMAL2
        if(last_modification<0)
        return true;
        dbg_delta();
        UnweightedDijkstra<Weight> d(source,g);

        for(int i = 0;i<g.nodes();i++) {
            int dis = dist[i];
            bool c = i<dist.size() && dist[i]<INF;
            if(!c)
            dis = this->unreachable();
            if(maxDistance>=0 && dis>maxDistance) {
                dis=-1;
            }
            int dbgdist = d.distance(i);
            if(maxDistance>=0 && dbgdist>maxDistance) {
                dbgdist=-1;
            }
            if(dis!=dbgdist) {
                assert(false);
                throw std::runtime_error("");
            }
            if(d.connected(i) && d.distance(i)<maxDistance){

                int dd =d.dist[i];
                int mdis = dist[i];
                if(! (mdis< INF)){
                    assert(false);
                    throw std::runtime_error("");
                }
                if(dd!=mdis) {
                    assert(false);
                    throw std::runtime_error("");
                }
            }else{
                if(dist[i]<maxDistance){
                    assert(false);
                    throw std::runtime_error("");
                }
            }

        }
//#endif
#endif
        return true;
    }

    bool connected_unsafe(int t) override{
        dbg_uptodate();
        return t < dist.size() && dist[t] < INF;
    }

    bool connected_unchecked(int t) override{
        assert(last_modification == g.getCurrentHistory());
        return connected_unsafe(t);
    }

    bool connected(int t) override{

        update();


        assert(dbg_uptodate());

        return dist[t] < INF;
    }

    int& distance(int t) override{

        update();

        if(connected_unsafe(t))
            return dist[t];
        else
            return this->unreachable();
    }

    int& distance_unsafe(int t) override{
        if(connected_unsafe(t))
            return dist[t];
        else
            return this->unreachable();
    }

    int incomingEdge(int t) override{

        if(!connected_unsafe(t)){
            return -1;
        }
        if(t == source)
            return -1;

        int d = dist[t];
        assert(delta[t] > 0);
        assert(d >= 0);
        assert(d != INF);
        int prev = -1;
        int prev_edgeID = -1;


        for(int i = 0; i < g.nIncoming(t); i++){

            int edgeID = g.incoming(t, i).id;
            if(edgeInShortestPathGraph[edgeID]){
                assert(g.edgeEnabled(edgeID));
                int from = g.incoming(t, i).node;
                assert(connected_unsafe(from));
                int from_dist = dist[from];
                assert(from_dist >= 0);
                assert(from_dist != INF);
                prev = from;
                prev_edgeID = edgeID;
                break;
            }
        }
        assert(prev != -1);


        return prev_edgeID;
    }

    int previous(int t) override{

        if(!connected_unsafe(t)){
            return -1;
        }
        if(t == source)
            return -1;

        int d = dist[t];

        assert(d >= 0);
        assert(d != INF);
        assert(delta[t] > 0);
        int prev = -1;

        for(int i = 0; i < g.nIncoming(t); i++){

            int edgeID = g.incoming(t, i).id;
            if(edgeInShortestPathGraph[edgeID]){
                assert(g.edgeEnabled(edgeID));
                int from = g.incoming(t, i).node;
                assert(connected_unsafe(from));
                int from_dist = dist[from];
                assert(from_dist >= 0);
                assert(from_dist != INF);
                prev = from;
                break;
            }
        }
        assert(prev != -1);


        return prev;
    }

    int randomIncomingEdge(int t, double& seed) override{
        if(!connected_unsafe(t)){
            return -1;
        }
        if(t == source)
            return -1;

        int d = dist[t];
        assert(delta[t] > 0);
        assert(d >= 0);
        assert(d != INF);
        int prev = -1;
        int prev_edgeID = -1;


        int start = dgl::alg::irand(seed, g.nIncoming(t));
        assert(start >= 0);
        assert(start < g.nIncoming(t));
        for(int j = 0; j < g.nIncoming(t); j++){
            int i = (j + start) % g.nIncoming(t);
            assert(i >= 0);
            assert(i < g.nIncoming(t));

            int edgeID = g.incoming(t, i).id;
            if(edgeInShortestPathGraph[edgeID]){
                assert(g.edgeEnabled(edgeID));
                int from = g.incoming(t, i).node;
                assert(connected_unsafe(from));
                int from_dist = dist[from];
                assert(from_dist >= 0);
                assert(from_dist != INF);
                prev = from;
                prev_edgeID = edgeID;
                break;
            }
        }
        assert(prev != -1);


        return prev_edgeID;
    }

    int randomPrevious(int t, double& seed) override{
        if(!connected_unsafe(t)){
            return -1;
        }
        if(t == source)
            return -1;

        int d = dist[t];

        assert(d >= 0);
        assert(d != INF);
        assert(delta[t] > 0);
        int prev = -1;

        int start = dgl::alg::irand(seed, g.nIncoming(t));
        assert(start >= 0);
        assert(start < g.nIncoming(t));
        for(int j = 0; j < g.nIncoming(t); j++){
            int i = (j + start) % g.nIncoming(t);
            assert(i >= 0);
            assert(i < g.nIncoming(t));

            int edgeID = g.incoming(t, i).id;
            if(edgeInShortestPathGraph[edgeID]){
                assert(g.edgeEnabled(edgeID));
                int from = g.incoming(t, i).node;
                assert(connected_unsafe(from));
                int from_dist = dist[from];
                assert(from_dist >= 0);
                assert(from_dist != INF);
                prev = from;
                break;
            }
        }
        assert(prev != -1);


        return prev;
    }
};

template<typename Weight, typename Graph, class Status>
bool RamalRepsBatchedUnified<Weight, Graph, Status>::ever_warned_about_zero_weights = false;
};
#endif
