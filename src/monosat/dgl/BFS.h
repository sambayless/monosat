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

#ifndef BFSReachability_H_
#define BFSReachability_H_

#include <vector>
#include "monosat/dgl/alg/Heap.h"
#include "Graph.h"
#include "DynamicGraph.h"
#include "monosat/core/Config.h"
#include "Reach.h"
#include "Distance.h"

namespace dgl {

template<typename Weight, typename Graph = DynamicGraph<Weight>, class Status = Reach::NullStatus, bool undirected = false>
class BFSReachability : public Reach {
public:

    Graph& g;
    Status& status;
    int last_modification;
    int last_addition;
    int last_deletion;
    int history_qhead;

    int last_history_clear;

    int source;
    int INF;

    bool opt_skip_deletions = false;
    bool opt_skip_additions = false;
    bool opt_inc_graph = false;
    int opt_dec_graph = 0;

    std::vector<int> q;
    std::vector<int> check;
    const int reportPolarity;

    std::vector<char> seen;

    std::vector<int> prev;

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

    BFSReachability(int s, Graph& graph, Status& _status = Reach::nullStatus, int _reportPolarity = 0) :
            g(graph), status(_status), last_modification(-1), last_addition(-1), last_deletion(-1), history_qhead(0),
            last_history_clear(
                    0), source(s), INF(0), reportPolarity(_reportPolarity){

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

    void setSource(int s) override{
        source = s;
        last_modification = -1;
        last_addition = -1;
        last_deletion = -1;
    }

    int getSource() override{
        return source;
    }


    void setNodes(int n){
        q.reserve(n);
        check.reserve(n);
        seen.resize(n);
        prev.resize(n);
        INF = g.nodes() + 1;
    }

    inline void add_update(int to, bool update){
        q.clear();
        q.push_back(to);

        for(int i = 0; i < q.size(); i++){
            int u = q[i];
            assert(seen[u]);
            if(update)
                status.setReachable(u, seen[u]);

            for(int i = 0; i < g.nIncident(u, undirected); i++){
                if(!g.edgeEnabled(g.incident(u, i, undirected).id))
                    continue;
                int v = g.incident(u, i, undirected).node;
                int edgeID = g.incident(u, i, undirected).id;
                if(!seen[v]){
                    seen[v] = 1;
                    prev[v] = edgeID;
                    q.push_back(v);
                }
            }
        }
    }

    inline void delete_update(int to){
        q.clear();
        q.push_back(to);
        seen[to] = false;
        prev[to] = -1;
        check.clear();
        check.push_back(to);

        for(int i = 0; i < q.size(); i++){
            int u = q[i];
            assert(!seen[u]);

            for(int i = 0; i < g.nIncident(u, undirected); i++){
                int v = g.incident(u, i, undirected).node;
                if(seen[v] && previous(v) == u){
                    seen[v] = 0;
                    prev[v] = -1;
                    check.push_back(v);
                    q.push_back(v);
                }
            }
        }

        q.clear();
        for(int i = 0; i < check.size(); i++){
            int u = check[i];
            if(!seen[u]){
                if(!undirected){
                    for(int i = 0; i < g.nIncoming(u); i++){
                        if(g.edgeEnabled(g.incoming(u, i).id)){
                            int from = g.incoming(u, i).node;
                            int edgeID = g.incoming(u, i).id;
                            int to = u;
                            if(seen[from]){
                                seen[to] = 1;
                                prev[to] = edgeID;
                                add_update(to, false);
                                break;
                            }
                        }
                    }
                }else{
                    for(int i = 0; i < g.nIncident(u, undirected); i++){

                        if(g.edgeEnabled(g.incident(u, i, undirected).id)){
                            int from = g.incident(u, i, undirected).node;
                            int edgeID = g.incident(u, i, undirected).id;
                            assert(from != u);
                            int to = u;
                            if(seen[from]){
                                seen[to] = 1;
                                prev[to] = edgeID;
                                add_update(to, false);
                                break;
                            }
                        }
                    }
                }
            }
        }

    }

    bool update_additions(){

        if(g.nHistoryClears() != last_history_clear){
            last_history_clear = g.nHistoryClears();
            history_qhead = 0;
        }

        assert(INF > g.nodes());
        assert(seen.size() >= g.nodes());
        q.clear();

        for(int i = history_qhead; i < g.historySize(); i++){
            int edgeid = g.getChange(i).id;
            int from = g.getEdge(edgeid).from;
            int to = g.getEdge(edgeid).to;
            if(g.getChange(i).addition){
                //incrementally add edge

                if(seen[from] && !seen[to]){
                    seen[to] = 1;
                    prev[to] = edgeid;
                    add_update(to, true);
                }else if(undirected){
                    if(seen[to] && !seen[from]){
                        seen[from] = 1;
                        prev[from] = edgeid;
                        add_update(from, true);
                    }
                }

            }else if(!undirected
                     && (to == source || !seen[from] || (seen[to] && seen[previous(to)] && previous(to) != from))){
                //then deleting this edge has no impact on BFSReachability, so don't need to do anything
            }else if(undirected
                     && ((to == source || !seen[from] || (seen[to] && seen[previous(to)] && previous(to) != from))
                         && (from == source || !seen[to]
                             || (seen[from] && seen[previous(from)] && previous(from) != to))

                     )){
                //then deleting this edge has no impact on BFSReachability, so don't need to do anything
            }else{
                stats_fast_failed_updates++;

                return false;
            }

        }

        stats_fast_updates++;


        assert(dbg_uptodate());

        last_modification = g.getCurrentHistory();
        last_deletion = g.nDeletions();
        last_addition = g.nAdditions();

        history_qhead = g.historySize();
        last_history_clear = g.nHistoryClears();

        return true;
    }

    bool incrementalUpdate(){

        if(g.nHistoryClears() != last_history_clear){
            last_history_clear = g.nHistoryClears();
            history_qhead = 0;
        }

        assert(INF > g.nodes());
        assert(seen.size() >= g.nodes());
        q.clear();

        for(int i = history_qhead; i < g.historySize(); i++){
            int edgeid = g.getChange(i).id;
            int from = g.getEdge(edgeid).from;
            int to = g.getEdge(edgeid).to;

            if(g.getChange(i).addition && g.edgeEnabled(g.getChange(i).id)){
                //incrementally add edge

                if(seen[from] && !seen[to]){
                    seen[to] = 1;
                    prev[to] = edgeid;
                    add_update(to, true);
                }else if(undirected){
                    if(seen[to] && !seen[from]){
                        seen[from] = 1;
                        prev[from] = edgeid;
                        add_update(from, true);
                    }
                }
            }else if(!g.getChange(i).addition && !g.edgeEnabled(g.getChange(i).id)){

                if(!undirected
                   && (to == source || !seen[from] || (seen[to] && seen[previous(to)] && previous(to) != from))){
                    //then deleting this edge has no impact on BFSReachability, so don't need to do anything
                }else if(undirected
                         && ((to == source || !seen[from] || (seen[to] && seen[previous(to)] && previous(to) != from))
                             && (from == source || !seen[to]
                                 || (seen[from] && seen[previous(from)] && previous(from) != to))

                         )){
                    //then deleting this edge has no impact on BFSReachability, so don't need to do anything
                }else{
                    delete_update(to);
                }
            }
        }

        for(int u = 0; u < g.nodes(); u++){
            status.setReachable(u, seen[u]);
        }

        stats_fast_updates++;


        assert(dbg_uptodate());

        last_modification = g.getCurrentHistory();
        last_deletion = g.nDeletions();
        last_addition = g.nAdditions();

        history_qhead = g.historySize();
        last_history_clear = g.nHistoryClears();

        return true;
    }

    int64_t num_updates = 0;

    int numUpdates() const override{
        return num_updates;
    }

    void update() override{
        if(last_modification > 0 && g.getCurrentHistory() == last_modification){
            stats_skipped_updates++;
            assert(dbg_uptodate());
            return;
        }

        if(last_modification > 0 && last_deletion == g.nDeletions()){
            stats_num_skipable_deletions++;
            if(opt_skip_deletions && reportPolarity < 1){
                assert(dbg_uptodate());
                return;            //I don't trust the correctness of these shortcuts
            }
        }
        if(last_modification > 0 && last_addition == g.nAdditions()){
            if(opt_skip_additions && reportPolarity > -1){
                assert(dbg_uptodate());
                return;            //I don't trust the correctness of these shortcuts
            }
        }

        setNodes(g.nodes());

        if(g.nHistoryClears() != last_history_clear){
            last_history_clear = g.nHistoryClears();
            history_qhead = 0;
        }else if(opt_inc_graph && last_modification > 0 && (g.nHistoryClears() <= (last_history_clear +
                                                                                   1))){
            if(opt_dec_graph == 2){
                if(incrementalUpdate())
                    return;
            }else{
                if(opt_dec_graph == 1 && last_deletion < g.nDeletions()){

                    //scan through the deletions and check if any of them matter..
                    bool safe = true;
                    for(int i = history_qhead; i < g.historySize(); i++){
                        int edgeid = g.getChange(i).id;
                        int from = g.getEdge(edgeid).from;
                        int to = g.getEdge(edgeid).to;
                        if(g.getChange(i).addition){
                            //safe
                        }else if(!undirected
                                 && (!seen[from] || (seen[to] && seen[previous(to)] && previous(to) != from))){
                            //then deleting this edge has no impact on BFSReachability, so don't need to do anything
                        }else if(undirected
                                 && ((!seen[from] || (seen[to] && seen[previous(to)] && previous(to) != from))
                                     && (!seen[to] || (seen[from] && seen[previous(from)] && previous(from) != to)))){
                            //then deleting this edge has no impact on BFSReachability, so don't need to do anything
                        }else{
                            safe = false;
                            break;
                        }
                    }
                    if(safe){
                        last_deletion = g.nDeletions();
                    }

                }

                if(last_deletion == g.nDeletions()){
                    if(update_additions())
                        return;
                }
            }
        }

        stats_full_updates++;

        q.clear();
        for(int i = 0; i < g.nodes(); i++){
            seen[i] = 0;
            prev[i] = -1;
        }

        seen[source] = 1;
        q.push_back(source);
        for(int i = 0; i < q.size(); i++){
            int u = q[i];
            assert(seen[u]);
            if(reportPolarity == 1)
                status.setReachable(u, true);

            for(int i = 0; i < g.nIncident(u, undirected); i++){
                if(!g.edgeEnabled(g.incident(u, i, undirected).id))
                    continue;
                int v = g.incident(u, i, undirected).node;
                int edgeID = g.incident(u, i, undirected).id;
                if(!seen[v]){
                    seen[v] = 1;
                    prev[v] = edgeID;
                    q.push_back(v);
                }
            }
        }

        if(reportPolarity < 1){
            for(int u = 0; u < g.nodes(); u++){
                if(!seen[u]){
                    status.setReachable(u, false);
                }else if(reportPolarity == 0){
                    status.setReachable(u, true);
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

        ;
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

    void drawFull(){
        printf("digraph{\n");
        for(int i = 0; i < g.nodes(); i++){

            if(seen[i]){
                printf("n%d [fillcolor=blue style=filled]\n", i);
            }else{
                printf("n%d \n", i);
            }

        }

        for(int i = 0; i < g.nodes(); i++){
            for(int j = 0; j < g.nIncident(i, undirected); j++){
                int id = g.incident(i, j, undirected).id;
                int u = g.incident(i, j, undirected).node;
                const char* s = "black";
                if(g.edgeEnabled(id))
                    s = "blue";
                else
                    s = "red";

                printf("n%d -> n%d [label=\"v%d\",color=\"%s\"]\n", i, u, id, s);
            }
        }

        printf("}\n");
    }

    bool dbg_uptodate(){
        return true;
    }

    bool connected_unsafe(int t) override{
        return t < seen.size() && seen[t];
    }

    bool connected_unchecked(int t) override{
        assert(last_modification == g.getCurrentHistory());
        return connected_unsafe(t);
    }

    bool connected(int t) override{
        if(last_modification != g.getCurrentHistory())
            update();

        assert(dbg_uptodate());

        return seen[t];
    }

    int distance(int t){
        if(connected(t))
            return 1;
        else
            return -1;
    }

    int distance_unsafe(int t){
        if(connected_unsafe(t))
            return 1;
        else
            return -1;
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
};

/**
 * Detect connectivity within a number of steps in unweighted, directed graphs
 */
template<typename Weight, typename Graph = DynamicGraph<Weight>, class Status = Distance<int>::NullStatus, bool undirected = false>
class UnweightedBFS : public Distance<int> {
    using Distance<int>::inf;
    using Distance<int>::unreachable;
public:

    Graph& g;
    Status& status;
    int last_modification;
    int last_addition;
    int last_deletion;
    int history_qhead;

    int last_history_clear;

    int source;

    int maxDistance;

    std::vector<int> q;
    std::vector<int> check;
    const int reportPolarity;

    std::vector<int> dist;

    std::vector<int> prev;

    //stats

    int stats_full_updates;
    int stats_fast_updates;
    int stats_fast_failed_updates;
    int stats_skip_deletes;
    int stats_skipped_updates;
    int stats_num_skipable_deletions;
    double mod_percentage;

    double stats_full_update_time;
    double stats_fast_update_time;

public:

    UnweightedBFS(int s, Graph& graph, int _reportPolarity = 0) :
            g(graph), status(Distance<int>::nullStatus), last_modification(-1), last_addition(-1), last_deletion(-1),
            history_qhead(
                    0), last_history_clear(0), source(s), reportPolarity(_reportPolarity){
        maxDistance = -1;
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

    UnweightedBFS(int s, Graph& graph, Status& _status, int _reportPolarity = 0) :
            g(graph), status(_status), last_modification(-1), last_addition(-1), last_deletion(-1), history_qhead(0),
            last_history_clear(
                    0), source(s), reportPolarity(_reportPolarity){
        maxDistance = -1;
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

    void setMaxDistance(int& _maxDistance) override{
        if(_maxDistance < 0){
            maxDistance = inf();
        }else
            maxDistance = _maxDistance;
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

    void setNodes(int n){
        q.reserve(n);
        dist.resize(n);
        prev.resize(n);
        if(maxDistance < 0)
            maxDistance = inf();
    }

    int64_t num_updates = 0;

    int numUpdates() const override{
        return num_updates;
    }

    void update() override{
        if(last_modification > 0 && g.getCurrentHistory() == last_modification){
            stats_skipped_updates++;
            return;
        }
        stats_full_updates++;

        if(last_deletion == g.nDeletions()){
            stats_num_skipable_deletions++;
        }

        setNodes(g.nodes());

        if(g.nHistoryClears() != last_history_clear){
            last_history_clear = g.nHistoryClears();
            history_qhead = 0;
        }

        q.clear();
        for(int i = 0; i < g.nodes(); i++){
            dist[i] = inf();
            prev[i] = -1;
        }

        dist[source] = 0;
        q.push_back(source);
        for(int i = 0; i < q.size(); i++){
            int u = q[i];
            assert(dist[u] < inf());
            if(reportPolarity >= 0)
                status.setMininumDistance(u, true, dist[u]);
            int d = dist[u];
            for(int i = 0; i < g.nIncident(u, undirected); i++){
                if(!g.edgeEnabled(g.incident(u, i, undirected).id))
                    continue;
                int edgeID = g.incident(u, i, undirected).id;
                int v = g.incident(u, i, undirected).node;
                int dv = dist[v];
                int alt = d + 1;
                if(alt > maxDistance)
                    alt = inf();            //Abort BFS early
                if(dist[v] > alt){
                    dist[v] = alt;
                    prev[v] = edgeID;
                    q.push_back(v);
                }
            }
        }

        if(reportPolarity <= 0){
            for(int u = 0; u < g.nodes(); u++){
                if(dist[u] >= inf()){
                    status.setMininumDistance(u, dist[u] < inf(), dist[u]);
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

        ;
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

    void drawFull(){
        printf("digraph{\n");
        for(int i = 0; i < g.nodes(); i++){

            if(dist[i] < inf()){
                printf("n%d [label=\"n%d %d \" fillcolor=blue style=filled]\n", i, i, dist[i]);
            }else{
                printf("n%d \n", i);
            }

        }

        for(int i = 0; i < g.nodes(); i++){
            for(int j = 0; j < g.nIncident(i, undirected); j++){
                int id = g.incident(i, j).id;
                int u = g.incident(i, j).node;
                const char* s = "black";
                if(g.edgeEnabled(id))
                    s = "blue";
                else
                    s = "red";

                printf("n%d -> n%d [label=\"v%d\",color=\"%s\"]\n", i, u, id, s);
            }
        }

        printf("}\n");
    }

    bool dbg_uptodate(){
#ifdef DEBUG_DIJKSTRA
        if(last_modification<=0)
        return true;
        UnweightedDijkstra<Reach::NullStatus, undirected> d(source,g);
        d.update();
        //drawFull();

        for(int i = 0;i<g.nodes();i++) {
            int distance = dist[i];
            int dbgdist = d.dist[i];
            if(dbgdist>maxDistance)
            dbgdist= inf();
            assert(distance==dbgdist);
        }
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
        if(connected(t)){
            return dist[t];
        }else{
            return inf();
        }
    }

    int& distance_unsafe(int t) override{
        if(connected_unsafe(t))
            return dist[t];
        else
            return inf();
    }

    int incomingEdge(int t) override{

        assert(t >= 0 && t < prev.size());
        assert(prev[t] >= -1);
        return prev[t];
    }

    int previous(int t) override{

        if(incomingEdge(t) < 0)
            return -1;
        if(undirected && g.getEdge(incomingEdge(t)).from == t){
            return g.getEdge(incomingEdge(t)).to;
        }
        assert(g.getEdge(incomingEdge(t)).to == t);
        return g.getEdge(incomingEdge(t)).from;
    }

};

};
#endif
