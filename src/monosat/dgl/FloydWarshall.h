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

#ifndef FLOYD_WARSHALL_H_
#define FLOYD_WARSHALL_H_

#include <vector>
#include "monosat/dgl/alg/Heap.h"
#include "Graph.h"
#include "monosat/core/Config.h"
#include "AllPairs.h"
#include "monosat/mtl/Sort.h"

namespace dgl {

template<typename Weight, typename Graph = DynamicGraph<Weight>, class Status = AllPairs::NullStatus>
class FloydWarshall : public AllPairs {
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

    std::vector<int> order;
//	std::vector<int> q;
//	std::vector<int> check;
    const int reportPolarity;

    std::vector<std::vector<int>> dist;
    std::vector<std::vector<int>> next;

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

    FloydWarshall(Graph& graph, Status& _status = AllPairs::nullStatus, int _reportPolarity = 0) :
            g(graph), status(_status), last_modification(-1), last_addition(-1), last_deletion(-1), history_qhead(0),
            last_history_clear(
                    0), INF(0), reportPolarity(0){

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
        //q.reserve(n);
        //check.reserve(n);
        order.clear();
        for(int i = 0; i < n; i++)
            order.push_back(i);
        INF = g.nodes() + 1;
        if(dist.size() < n){
            dist.resize(n);
            for(int i = 0; i < dist.size(); i++)
                dist[i].resize(n);
        }
        if(next.size() < n){
            next.resize(n);
            for(int i = 0; i < next.size(); i++)
                next[i].resize(n);
        }
    }

    struct lt_key {
        std::vector<int>& _dist;

        bool operator()(int a, int b) const{
            return _dist[a] < _dist[b];
        }

        lt_key(std::vector<int>& d) :
                _dist(d){
        };
    };

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

        setNodes(g.nodes());

        for(int i = 0; i < g.nodes(); i++){
            for(int j = 0; j < g.nodes(); j++){
                next[i][j] = -1;
                dist[i][j] = INF;
            }
            dist[i][i] = 0;
        }

        for(int i = 0; i < g.edges(); i++){
            if(g.hasEdge(i) && g.edgeEnabled(i)){
                int u = g.getEdge(i).from;
                int v = g.getEdge(i).to;
                if(u != v)
                    dist[u][v] = 1;
            }
        }

        //for(int l = 0;l<sources.size();l++){
        //	int k = sources[l];
        for(int k = 0; k < g.nodes(); k++){
            for(int i = 0; i < g.nodes(); i++){
                for(int j = 0; j < g.nodes(); j++){
                    int d = dist[i][k] + dist[k][j];
                    if(d < dist[i][j]){
                        dist[i][j] = d;
                        next[i][j] = k;
                    }
                }
            }
        }
        for(int i = 0; i < sources.size(); i++){
            int s = sources[i];
            //sort(order,lt_key(dist[s]));//disabled because it is NOT required
            for(int j = 0; j < order.size(); j++){
                int u = j;        // order[j];
                /*		if(j<order.size()-1){
                 assert(dist[s][u]<=dist[s][order[j+1]]);
                 }*/
                //Wrong. This is only required if we are returning learnt clauses that include other reachability lits.
                //it is crucial to return the nodes in order of distance, so that they are enqueued in the correct order in the solver.
                if(dist[s][u] >= INF && reportPolarity < 1){
                    status.setReachable(s, u, false);
                    status.setMininumDistance(s, u, false, INF);
                }else if(dist[s][u] < INF && reportPolarity > -1){
                    status.setReachable(s, u, true);
                    status.setMininumDistance(s, u, true, dist[s][u]);
                }
            }
            /*	for(int u = 0;u<g.nodes();u++){
             if(dist[s][u]>=INF && reportPolarity<1){
             //status.setReachable(s,u,false);
             //status.setMininumDistance(s,u,false,INF);
             }else if(dist[s][u]<INF && reportPolarity>-1){
             //status.setReachable(s,u,true);
             //status.setMininumDistance(s,u,true,dist[s][u]);
             }
             }*/
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
        path.push_back(from);
        getPath_private(from, to, path);
        assert(path.back() != to);
        path.push_back(to);
    }

    void getPath_private(int from, int to, std::vector<int>& path){
        assert(dist[from][to] < INF);
        int intermediate = next[from][to];
        if(intermediate > -1){
            getPath_private(from, intermediate, path);
            path.push_back(intermediate);
            getPath_private(intermediate, to, path);

        }

    }

    bool dbg_path(int from, int to){

        return true;
    }

    void drawFull(){
        /*printf("digraph{\n");
         for(int i = 0;i< g.nodes();i++){

         if(seen[i]){
         printf("n%d [fillcolor=blue style=filled]\n", i);
         }else{
         printf("n%d \n", i);
         }


         }

         for(int i = 0;i< g.nodes();i++){
         for(int j =0;j<g.nIncident(u);j++){
         int id  =g.incident(i,j).id;
         int u =  g.incident(i,j).node;
         const char * s = "black";
         if( g.edgeEnabled(id))
         s="blue";
         else
         s="red";



         printf("n%d -> n%d [label=\"v%d\",color=\"%s\"]\n", i,u, id, s);
         }
         }

         printf("}\n");*/
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
    /*	int previous(int t){
     assert(t>=0 && t<next.size());
     assert(next[t]>=-1 && next[t]<next.size());
     return next[t];
     }*/

};
};
#endif
