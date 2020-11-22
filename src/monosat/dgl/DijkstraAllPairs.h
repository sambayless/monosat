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

#ifndef DIJKSTRA_ALLPAIRS_H_
#define DIJKSTRA_ALLPAIRS_H_

#include <vector>
#include "monosat/dgl/alg/Heap.h"
#include "Graph.h"
#include "monosat/core/Config.h"
#include "AllPairs.h"


namespace dgl {
template<typename Weight, typename Graph = DynamicGraph<Weight>, class Status = AllPairs::NullStatus>
class DijkstraAllPairs : public AllPairs {
public:

    Graph& g;
    Status& status;
    int last_modification;
    int last_addition;
    int last_deletion;
    int history_qhead;

    int last_history_clear;

    std::vector<int> sources;
    int INF;

    std::vector<int> check;
    const int reportPolarity;

    std::vector<std::vector<int>> dist;
    std::vector<std::vector<int>> prev;

    std::vector<int>* dist_ptr = nullptr;

    struct DistCmp {
        std::vector<int>** _dist;

        bool operator()(int a, int b) const{
            return (**_dist)[a] < (**_dist)[b];
        }

        DistCmp(std::vector<int>** d) :
                _dist(d){
        };
    };

    alg::Heap<DistCmp> q;

public:
    int stats_full_updates;
    int stats_fast_updates;
    int stats_fast_failed_updates;
    int stats_skip_deletes;
    int stats_skipped_updates;
    int stats_num_skipable_deletions;
    double mod_percentage;

    double stats_full_update_time;
    double stats_fast_update_time;

    DijkstraAllPairs(Graph& graph, Status& _status, int _reportPolarity = 0) :
            g(graph), status(_status), last_modification(-1), last_addition(-1), last_deletion(-1), history_qhead(0),
            last_history_clear(
                    0), INF(0), reportPolarity(0), q(DistCmp(&dist_ptr)){

        mod_percentage = 0.2;
        stats_full_updates = 0;
        stats_fast_updates = 0;
        stats_skip_deletes = 0;
        stats_skipped_updates = 0;
        stats_full_update_time = 0;
        stats_fast_update_time = 0;
        stats_num_skipable_deletions = 0;
        stats_fast_failed_updates = 0;
    }

    void addSource(int s) override{
        assert(!std::count(sources.begin(), sources.end(), s));
        sources.push_back(s);

        last_modification = -1;
        last_addition = -1;
        last_deletion = -1;
    }

    void setNodes(int n){

        check.reserve(n);

        INF = g.nodes() + 1;
        if(dist.size() < n){
            dist.resize(n);

            for(int i = 0; i < dist.size(); i++)
                dist[i].resize(n);

            prev.resize(n);
            for(int i = 0; i < dist.size(); i++)
                prev[i].resize(n);
        }
    }

    int64_t num_updates = 0;

    int numUpdates() const override{
        return num_updates;
    }

    void update() override{
        stats_full_updates++;

        if(last_modification > 0 && g.getCurrentHistory() == last_modification){
            stats_skipped_updates++;
            return;
        }

        if(last_deletion == g.nDeletions()){
            stats_num_skipable_deletions++;
        }
        INF = g.nodes() + 1;
        stats_full_updates++;

        setNodes(g.nodes());
        q.clear();
        for(int i = 0; i < g.nodes(); i++){
            for(int j = 0; j < g.nodes(); j++){
                prev[i][j] = -1;
                dist[i][j] = INF;
            }
            dist[i][i] = 0;
        }

        for(int s = 0; s < sources.size(); s++){
            int source = sources[s];
            q.clear();
            dist_ptr = &dist[source];
            dist[source][source] = 0;
            q.insert(source);
            while(q.size()){
                int u = q.peekMin();
                if(dist[source][u] == INF)
                    break;

                q.removeMin();
                for(int i = 0; i < g.nIncident(u); i++){
                    if(!g.edgeEnabled(g.incident(u, i).id))
                        continue;
                    int v = g.incident(u, i).node;
                    int alt = dist[source][u] + 1;
                    if(alt < dist[source][v]){
                        dist[source][v] = alt;
                        prev[source][v] = u;
                        if(!q.inHeap(v))
                            q.insert(v);
                        else
                            q.decrease(v);
                    }
                }
            }
        }

        for(int i = 0; i < sources.size(); i++){
            int s = sources[i];
            for(int u = 0; u < g.nodes(); u++){
                if(dist[s][u] >= INF && reportPolarity < 1){
                    status.setReachable(s, u, false);
                    status.setMininumDistance(s, u, false, INF);
                }else if(dist[s][u] < INF && reportPolarity > -1){
                    status.setReachable(s, u, true);
                    status.setMininumDistance(s, u, true, dist[s][u]);
                }
            }
        }
        assert(dbg_uptodate());
        num_updates++;
        last_modification = g.getCurrentHistory();
        last_deletion = g.nDeletions();
        last_addition = g.nAdditions();

        history_qhead = g.historySize();
        last_history_clear = g.nHistoryClears();

    }

    void getPath(int from, int to, std::vector<int>& path) override{
        update();
        assert(dist[from][to] < INF);
        int d = dist[from][to];
        path.clear();
        path.push_back(to);
        while(prev[from][to] != -1){
            int p = prev[from][to];
            path.push_back(p);
            to = p;
        }

        for(int i = 0; i < path.size(); i++){
            if(i >= path.size() - i - 1)
                break;
            int t = path[i];
            int e = path[path.size() - i - 1];
            path[i] = e;
            path[path.size() - i - 1] = t;
        }
    }

    bool dbg_path(int from, int to){

        return true;
    }

    void drawFull(){
    }

    bool dbg_uptodate(){

        return true;
    }

    bool connected_unsafe(int from, int t) override{
        return dist[from][t] < INF;
    }

    bool connected_unchecked(int from, int t) override{
        assert(last_modification == g.getCurrentHistory());
        return connected_unsafe(from, t);
    }

    bool connected(int from, int t) override{
        if(last_modification != g.getCurrentHistory())
            update();

        assert(dbg_uptodate());

        return dist[from][t] < INF;
    }

    int distance(int from, int t) override{
        if(connected(from, t))
            return dist[from][t];
        else
            return INF;
    }

    int distance_unsafe(int from, int t) override{
        if(connected_unsafe(from, t))
            return dist[from][t];
        else
            return INF;
    }
};
};
#endif
