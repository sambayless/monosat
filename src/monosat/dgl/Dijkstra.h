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
 NONinf()RINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
 OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 **************************************************************************************************/

#ifndef DIJKSTRA_H_
#define DIJKSTRA_H_

#include <vector>
#include "monosat/dgl/alg/Heap.h"
#include "Graph.h"
#include "DynamicGraph.h"
#include "Reach.h"
#include <monosat/dgl/alg/Rnd.h>
#include "Distance.h"
#include "monosat/core/Config.h"
#include <limits>

namespace dgl {

template<typename Weight = int64_t, typename Graph = DynamicGraph<Weight>, class Status = typename Distance<Weight>::NullStatus, bool undirected = false>
class Dijkstra : public Distance<Weight> {
    using Distance<Weight>::inf;
    using Distance<Weight>::unreachable;
public:
    Graph& g;

    Status& status;
    int reportPolarity;

    int last_modification = -1;
    int last_addition = 0;
    int last_deletion = 0;
    int last_edge_inc = 0;
    int last_edge_dec = 0;
    int history_qhead = 0;

    int last_history_clear = 0;

    int source;
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

    alg::Heap<DistCmp> q;

public:

    int stats_full_updates = 0;
    int stats_fast_updates = 0;
    int stats_fast_failed_updates = 0;
    int stats_skip_deletes = 0;
    int stats_skipped_updates = 0;
    int stats_num_skipable_deletions = 0;
    double mod_percentage;

    double stats_full_update_time = 0;
    double stats_fast_update_time = 0;

    Dijkstra(int s, Graph& graph, Status& status, int reportPolarity = 0) :
            g(graph), status(status), reportPolarity(reportPolarity), source(s), q(DistCmp(dist)){

        mod_percentage = 0.2;

    }

    Dijkstra(int s, Graph& graph, int reportPolarity = 0) :
            g(graph), status(Distance<Weight>::nullStatus), reportPolarity(reportPolarity), source(s), q(DistCmp(dist)){

        mod_percentage = 0.2;

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

    void drawFull(){

    }

    int64_t num_updates = 0;

    int numUpdates() const override{
        return num_updates;
    }

    void update() override{
        if(last_modification > 0 && g.getCurrentHistory() == last_modification)
            return;

        if(last_addition == g.nAdditions() && last_edge_inc == g.lastEdgeIncrease() &&
           last_edge_dec == g.lastEdgeDecrease() && last_modification > 0){
            //if none of the deletions were to edges that were the previous edge of some shortest path, then we don't need to do anything
            if(last_history_clear != g.nHistoryClears()){
                history_qhead = 0;
                last_history_clear = g.nHistoryClears();
            }
            bool need_recompute = false;
            //ok, now check if any of the added edges allow for a decrease in distance.
            for(int i = history_qhead; i < g.historySize(); i++){
                assert(!g.getChange(i).addition);
                int edgeid = g.getChange(i).id;
                int u = g.getEdge(edgeid).from;
                int v = g.getEdge(edgeid).to;
                if(incomingEdge(u) == edgeid || incomingEdge(v) == edgeid){
                    history_qhead = i - 1;
                    need_recompute = true;
                    //this deletion matters, so we need to recompute.
                    break;
                }
            }
            if(!need_recompute){
                //none of these deletions touched any shortest paths, so we can ignore them.

                last_modification = g.getCurrentHistory();
                last_deletion = g.nDeletions();
                last_addition = g.nAdditions();
                last_edge_inc = g.lastEdgeIncrease();
                last_edge_dec = g.lastEdgeDecrease();
                history_qhead = g.historySize();
                last_history_clear = g.nHistoryClears();

                assert(dbg_uptodate());
                stats_skip_deletes++;
                return;
            }
        }

        stats_full_updates++;
        if(dist.size() != g.nodes()){

            dist.resize(g.nodes());
            prev.resize(g.nodes());
        }

        q.clear();
        for(int i = 0; i < g.nodes(); i++){
            dist[i] = inf();
            prev[i] = -1;
        }

        dist[source] = 0;
        q.insert(source);
        while(q.size()){
            int u = q.peekMin();
            if(dist[u] == inf())
                break;
            q.removeMin();
            for(int i = 0; i < g.nIncident(u, undirected); i++){
                if(!g.edgeEnabled(g.incident(u, i, undirected).id))
                    continue;
                int edgeID = g.incident(u, i, undirected).id;
                int v = g.incident(u, i, undirected).node;
                Weight alt = dist[u] + g.getWeight(edgeID);
                if(alt < dist[v]){
                    dist[v] = alt;
                    prev[v] = edgeID;
                    if(!q.inHeap(v))
                        q.insert(v);
                    else
                        q.decrease(v);
                }
            }
        }

        assert(dbg_uptodate());
        for(int u = 0; u < g.nodes(); u++){
            if(reportPolarity <= 0 && dist[u] >= inf()){
                status.setReachable(u, false);
                status.setMininumDistance(u, dist[u] < inf(), dist[u]);
            }else if(reportPolarity >= 0 && dist[u] < inf()){
                status.setReachable(u, true);
                status.setMininumDistance(u, dist[u] < inf(), dist[u]);
            }
        }
        num_updates++;
        last_modification = g.getCurrentHistory();
        last_deletion = g.nDeletions();
        last_addition = g.nAdditions();
        last_edge_inc = g.lastEdgeIncrease();
        last_edge_dec = g.lastEdgeDecrease();
        history_qhead = g.historySize();
        last_history_clear = g.nHistoryClears();

    }

    bool dbg_path(int to){
#ifdef DEBUG_DIJKSTRA
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

    bool dbg_uptodate(){
#ifdef DEBUG_DIJKSTRA
        if(last_modification<=0)
        return true;
#endif
        return true;
    }

    bool connected_unsafe(int t) override{
        return t < dist.size() && dist[t] < inf();
    }

    bool connected_unchecked(int t) override{
        assert(last_modification == g.getCurrentHistory());
        return connected_unsafe(t);
    }

    bool connected(int t) override{
        if(last_modification != g.getCurrentHistory())
            update();

        assert(dbg_uptodate());

        return dist[t] < inf();
    }

    Weight& distance(int t) override{
        if(last_modification != g.getCurrentHistory())
            update();
        if(connected_unsafe(t))
            return dist[t];
        return this->unreachable();
    }

    Weight& distance_unsafe(int t) override{
        if(connected_unsafe(t))
            return dist[t];
        else
            return this->unreachable();
    }

    int incomingEdge(int t) override{
        assert(t >= 0 && t < prev.size());
        assert(prev[t] >= -1);
        return prev[t];
    }

    int previous(int t) override{
        if(prev[t] < 0)
            return -1;
        if(undirected && g.getEdge(incomingEdge(t)).from == t){
            return g.getEdge(incomingEdge(t)).to;
        }
        assert(g.getEdge(incomingEdge(t)).to == t);
        return g.getEdge(incomingEdge(t)).from;
    }

    //Select a randomly chosen, but still shortest path, previous node or edge (falling back on selecting a deterministic edge if this functionality is not supported)
    int randomPrevious(int t, double& seed) override{
        assert(t >= 0 && t < prev.size());
        assert(prev[t] >= -1);
        int prev = -1;


        int start = dgl::alg::irand(seed, g.nIncoming(t));
        assert(start >= 0);
        assert(start < g.nIncoming(t));
        for(int j = 0; j < g.nIncoming(t); j++){
            int i = (j + start) % g.nIncoming(t);
            assert(i >= 0);
            assert(i < g.nIncoming(t));
            int edgeID = g.incoming(t, i).id;
            int from = g.incoming(t, i).node;
            if(g.edgeEnabled(edgeID)){
                if(connected_unsafe(from)){
                    Weight alt = dist[from] + g.getWeight(edgeID);
                    assert(alt >= dist[t]);
                    if(alt == dist[t]){
                        prev = from;
                        break;
                    }
                }
            }
        }
        assert(prev != -1);
        return prev;
    }

    int randomIncomingEdge(int t, double& seed) override{
        assert(t >= 0 && t < prev.size());
        assert(prev[t] >= -1);

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
            int from = g.incoming(t, i).node;
            if(g.edgeEnabled(edgeID)){
                if(connected_unsafe(from)){
                    Weight alt = dist[from] + g.getWeight(edgeID);
                    assert(alt >= dist[t]);
                    if(alt == dist[t]){
                        prev = from;
                        prev_edgeID = edgeID;
                        break;
                    }
                }
            }
        }
        assert(prev != -1);
        assert(prev_edgeID != -1);
        return prev_edgeID;
    }

};

template<typename Weight, typename Graph = DynamicGraph<Weight>, class Status =typename Distance<int>::NullStatus, bool undirected = false>
class UnweightedDijkstra : public Distance<int> {
    using Distance<int>::inf;
    using Distance<int>::unreachable;
public:
    Graph& g;
    Status& status;
    int reportPolarity;

    int last_modification = -1;
    int last_addition = -1;
    int last_deletion = -1;
    int history_qhead = 0;

    int last_history_clear = 0;

    int source;

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

    alg::Heap<DistCmp> q;

public:

    int stats_full_updates = 0;
    int stats_fast_updates = 0;
    int stats_fast_failed_updates = 0;
    int stats_skip_deletes = 0;
    int stats_skipped_updates = 0;
    int stats_num_skipable_deletions = 0;
    double mod_percentage;

    double stats_full_update_time = 0;
    double stats_fast_update_time = 0;

    UnweightedDijkstra(int s, Graph& graph, Status& status, int reportPolarity = 0) :
            g(graph), status(status), reportPolarity(reportPolarity), source(s), q(DistCmp(dist)){

        mod_percentage = 0.2;

    }

    UnweightedDijkstra(int s, Graph& graph, int reportPolarity = 0) :
            g(graph), status(Distance<int>::nullStatus), reportPolarity(reportPolarity), source(s), q(DistCmp(dist)){

        mod_percentage = 0.2;

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

    void drawFull(){

    }

    int64_t num_updates = 0;

    int numUpdates() const override{
        return num_updates;
    }

    void update() override{
        if(last_modification > 0 && g.getCurrentHistory() == last_modification)
            return;

        if(last_addition == g.nAdditions() && last_modification > 0){
            //if none of the deletions were to edges that were the previous edge of some shortest path, then we don't need to do anything
            if(last_history_clear != g.nHistoryClears()){
                history_qhead = 0;
                last_history_clear = g.nHistoryClears();
            }
            bool need_recompute = false;
            //ok, now check if any of the added edges allow for a decrease in distance.
            for(int i = history_qhead; i < g.historySize(); i++){
                assert(!g.getChange(i).addition);
                int edgeid = g.getChange(i).id;
                int u = g.getEdge(edgeid).from;
                int v = g.getEdge(edgeid).to;
                if(incomingEdge(u) == edgeid || incomingEdge(v) == edgeid){
                    history_qhead = i - 1;
                    need_recompute = true;
                    //this deletion matters, so we need to recompute.
                    break;
                }
            }
            if(!need_recompute){
                //none of these deletions touched any shortest paths, so we can ignore them.

                last_modification = g.getCurrentHistory();
                last_deletion = g.nDeletions();
                last_addition = g.nAdditions();

                history_qhead = g.historySize();
                last_history_clear = g.nHistoryClears();

                assert(dbg_uptodate());
                stats_skip_deletes++;
                return;
            }
        }

        stats_full_updates++;

        dist.resize(g.nodes());
        prev.resize(g.nodes());
        q.clear();
        for(int i = 0; i < g.nodes(); i++){
            //old_dist[i]=last_modification > 0 ? dist[i]:inf();//this won't work properly if we added nodes...
            dist[i] = inf();
            prev[i] = -1;
        }

        dist[source] = 0;
        q.insert(source);
        while(q.size()){
            int u = q.peekMin();
            if(dist[u] == inf())
                break;
            q.removeMin();
            for(int i = 0; i < g.nIncident(u, undirected); i++){
                if(!g.edgeEnabled(g.incident(u, i, undirected).id))
                    continue;
                int edgeID = g.incident(u, i, undirected).id;
                int v = g.incident(u, i, undirected).node;
                int alt = dist[u] + 1;
                if(alt < dist[v]){
                    dist[v] = alt;
                    prev[v] = edgeID;
                    if(!q.inHeap(v))
                        q.insert(v);
                    else
                        q.decrease(v);
                }
            }
        }

        assert(dbg_uptodate());
        for(int u = 0; u < g.nodes(); u++){
            if(reportPolarity <= 0 && dist[u] >= inf()){
                status.setReachable(u, false);
                status.setMininumDistance(u, dist[u] < inf(), dist[u]);
            }else if(reportPolarity >= 0 && dist[u] < inf()){
                status.setReachable(u, true);
                status.setMininumDistance(u, dist[u] < inf(), dist[u]);
            }
        }
        num_updates++;
        last_modification = g.getCurrentHistory();
        last_deletion = g.nDeletions();
        last_addition = g.nAdditions();

        history_qhead = g.historySize();
        last_history_clear = g.nHistoryClears();

    }

    bool dbg_path(int to){
#ifdef DEBUG_DIJKSTRA
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

    bool dbg_uptodate(){
#ifdef DEBUG_DIJKSTRA
        if(last_modification<=0)
        return true;
        /*	Graph gdbg;
         for(int i = 0;i<g.nodes();i++){
         gdbg.addNode();
         }

         for(int i = 0;i<g.nodes();i++){
         for(int j = 0;j<g.nIncident(u);j++){
         int u = g.incident(i,j);
         gdbg.addEdge(i,u);
         }
         }*/

        /*	Dijkstra d(source,g);
         d.update();
         for(int i = 0;i<g.nodes();i++){
         int distance = dist[i];
         int dbgdist = d.dist[i];
         assert(distance==dbgdist);
         }*/
#endif
        return true;
    }

    bool connected_unsafe(int t) override{
        return t < dist.size() && dist[t] < inf();
    }

    bool connected_unchecked(int t) override{
        assert(last_modification == g.getCurrentHistory());
        return connected_unsafe(t);
    }

    bool connected(int t) override{
        if(last_modification != g.getCurrentHistory())
            update();

        assert(dbg_uptodate());

        return dist[t] < inf();
    }

    int& distance(int t) override{
        if(last_modification != g.getCurrentHistory())
            update();
        if(connected_unsafe(t))
            return dist[t];
        return this->unreachable();
    }

    int& distance_unsafe(int t) override{
        if(connected_unsafe(t))
            return dist[t];
        else
            return this->unreachable(); //return inf();
    }

    int incomingEdge(int t) override{
        assert(t >= 0 && t < prev.size());
        assert(prev[t] >= -1);
        return prev[t];
    }

    int previous(int t) override{
        if(prev[t] < 0)
            return -1;
        if(undirected && g.getEdge(incomingEdge(t)).from == t){
            return g.getEdge(incomingEdge(t)).to;
        }
        assert(g.getEdge(incomingEdge(t)).to == t);
        return g.getEdge(incomingEdge(t)).from;
    }

    //Select a randomly chosen, but still shortest path, previous node or edge (falling back on selecting a deterministic edge if this functionality is not supported)
    int randomPrevious(int t, double& seed) override{
        assert(t >= 0 && t < prev.size());
        assert(prev[t] >= -1);

        int prev = -1;
        int start = dgl::alg::irand(seed, g.nIncoming(t));
        assert(start >= 0);
        assert(start < g.nIncoming(t));
        for(int j = 0; j < g.nIncoming(t); j++){
            int i = (j + start) % g.nIncoming(t);
            assert(i >= 0);
            assert(i < g.nIncoming(t));
            int edgeID = g.incoming(t, i).id;
            int from = g.incoming(t, i).node;
            if(g.edgeEnabled(edgeID)){
                if(connected_unsafe(from)){
                    Weight alt = dist[from] + 1;
                    assert(alt >= dist[t]);
                    if(alt == dist[t]){
                        prev = from;
                        break;
                    }
                }
            }
        }
        assert(prev != -1);
        //assert(min_prev_dist<dist[t]);

        return prev;
    }

    int randomIncomingEdge(int t, double& seed) override{
        assert(t >= 0 && t < prev.size());
        assert(prev[t] >= -1);

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
            if(g.edgeEnabled(edgeID)){
                int from = g.incoming(t, i).node;
                if(connected_unsafe(from)){
                    Weight alt = dist[from] + 1;
                    assert(alt >= dist[t]);
                    if(alt == dist[t]){
                        prev = from;
                        prev_edgeID = edgeID;
                        break;
                    }
                }
            }
        }
        assert(prev != -1);
        assert(prev_edgeID != -1);
        return prev_edgeID;
    }
};
};
#endif /* DIJKSTRA_H_ */
