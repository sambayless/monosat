//
// Created by sam on 20/03/17.
//

#ifndef MONOSAT_DSTARLITE_H
#define MONOSAT_DSTARLITE_H


#include <vector>
#include "monosat/dgl/alg/Heap.h"
#include "DynamicGraph.h"
#include "monosat/core/Config.h"
#include "Reach.h"

namespace dgl {

template<typename Weight,class Status = Reach::NullStatus, bool undirected = false>
class DStarLite: public Reach {
public:

    DynamicGraph<Weight> & g;
    Status & status;
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
    alg::Heap<DistCmp> U;

    std::vector<int> q;
    std::vector<int> check;
    const int reportPolarity;

    //std::vector<char> old_seen;
    std::vector<char> seen;
//	std::vector<int> changed;

    std::vector<int> prev;

    long stats_full_updates=0;
    long stats_fast_updates=0;
    long stats_fast_failed_updates=0;
    long stats_skip_deletes=0;
    long stats_skipped_updates=0;
    long stats_num_skipable_deletions=0;


    double stats_full_update_time=0;
    double stats_fast_update_time=0;
    std::vector<bool> removed_from_heap;
    std::vector<int> rhs;
    int destination=0;
    int k_m = 0;


    int H(int f, int t){
        return 0;
    }
    std::vector<int> g_vec;
    int G(int s){
        return g_vec[s];
    }
    int edgeCost(int edgeID){
        return g.edgeEnabled(edgeID) ? 1:INF;
    }
    struct Key{
        int a=0;
        int b=0;
        bool operator<(const Key& rhs) const {
            return (rhs.a < this->a) || (rhs.a == this->a && rhs.b < this->b);
        }
        bool operator>(const Key& rhs) const {
            return (rhs.a > this->a) || (rhs.a == this->a && rhs.b > this->b);
        }
    };

    Key calculateKey(int s){
        int r_out = std::min(G(s), rhs[s]);
        return Key{r_out + H(source, s) + k_m, r_out};
    }
    std::vector<Key> priorities;

    void init(){
        rhs.clear();
        rhs.resize(g.nodes());
        for(int i = 0;i<rhs.size();i++){
            rhs[i]=INF;
        }
        rhs[destination]=0;
        node_on_path.clear();
        node_on_path.resize(g.nodes());
        prev.clear();
        prev.resize(g.nodes(),-1);
        prev_edge.clear();
        prev_edge.resize(g.nodes(),-1);

        stats_full_updates++;
        int oldInf = INF;
        INF = g.nodes() + 1;
        if(INF!=oldInf){
            for(int i = 0;i<dist.size();i++){
                if(dist[i]==oldInf){
                    dist[i]=INF;
                }
            }
        }
        last_history_clear=-1;
        dist.resize(g.nodes(), INF);
        dist[getSource()] = 0;

        removed_from_heap.clear();
        removed_from_heap.resize(g.nodes(),false);

        priorities.resize(g.nodes());
        for(int i = 0;i<priorities.size();i++){
            priorities[i].a = INF;
            priorities[i].b = INF;
        }
        priorities[destination] =  calculateKey(destination);
        U.clear();
        U.insert(destination);


    }
    void updateVertex(int u){
        if (u!= destination){
            int m = INF;
            for (int i = 0; i < g.nIncident(u, undirected); i++) {
                if (!g.edgeEnabled(g.incident(u, i, undirected).id))
                    continue;
                int edgeID = g.incident(u, i, undirected).id;
                int v = g.incident(u, i, undirected).node;
                if (v!=u){
                    int n = edgeCost(edgeID)+g(v);
                    if (n<m){
                        m = n;
                    }
                }
            }
            rhs[u] = m;
        }

        if(G(u)!=rhs[u]){
            priorities[u]= calculateKey(u);
            U.update(u);//is this always a decrease or an increase?
        }else{
            removed_from_heap[u]=true;//u should be removed from U
        }
    }

    void computeShortestPath(int s_start){
        while(true){
            int u = U.peekMin();
            while (U.size() && removed_from_heap[u]){
                removed_from_heap[u]=false;
                U.removeMin();
                if(U.size()) {
                    u = U.peekMin();
                }else{
                    u=-1;
                }
            }
            if(u<0)
                break;
            Key k_old = priorities[u];
            if(k_old<priorities[source] || rhs[source]!=G(source)){
                Key k_new = calculateKey(u);
                if(k_old<k_new){
                    priorities[u]=k_new;
                    U.insert(u);
                }else if (G(u)>rhs[u]){
                    g_vec[u]=rhs[u];
                    for (int i = 0; i < g.nIncoming(u, undirected); i++) {
                        if (!g.edgeEnabled(g.incoming(u, i, undirected).id))
                            continue;
                        int edgeID = g.incoming(u, i, undirected).id;
                        int v = g.incoming(u, i, undirected).node;
                        if (v!=u){
                            updateVertex(v);
                        }
                    }
                }else{
                    g_vec[u]=INF;
                    updateVertex(u);
                    for (int i = 0; i < g.nIncoming(u, undirected); i++) {
                        if (!g.edgeEnabled(g.incoming(u, i, undirected).id))
                            continue;
                        int edgeID = g.incoming(u, i, undirected).id;
                        int v = g.incoming(u, i, undirected).node;
                        if (v!=u){
                            updateVertex(v);
                        }
                    }
                }
            }

        }
    }

public:


    DStarLite(int s, DynamicGraph<Weight> & graph, Status & _status= Reach::nullStatus, int _reportPolarity = 0) :
            g(graph), status(_status), last_modification(-1), last_addition(-1), last_deletion(-1), history_qhead(0), last_history_clear(
            0), source(s), INF(0), reportPolarity(_reportPolarity),U(DistCmp(dist)) {


    }
    DStarLite(DynamicGraph<Weight> & graph,Status & _status= Reach::nullStatus, int _reportPolarity = 0):			g(graph), status(_status), last_modification(-1), last_addition(-1), last_deletion(-1), history_qhead(0), last_history_clear(
            0), source(0), INF(0), reportPolarity(_reportPolarity),U(DistCmp(dist)) {

    }
    //Connectivity(const Connectivity& d):g(d.g), last_modification(-1),last_addition(-1),last_deletion(-1),history_qhead(0),last_history_clear(0),source(d.source),INF(0),mod_percentage(0.2),stats_full_updates(0),stats_fast_updates(0),stats_skip_deletes(0),stats_skipped_updates(0),stats_full_update_time(0),stats_fast_update_time(0){marked=false;};

    void setSource(int s) {
        source = s;
        last_modification = -1;
        last_addition = -1;
        last_deletion = -1;
    }
    int getSource() {
        return source;
    }

    std::vector<bool> node_on_path;

    std::vector<int> prev_edge;
    long num_updates = 0;
    int numUpdates() const {
        return num_updates;
    }
    void update() {
        static int iteration = 0;
        int local_it = ++iteration;
        if (last_modification > 0 && g.modifications == last_modification){
            return;
        }
        int s_start=source;
        int s_last=s_start;
        if (last_modification <= 0 || g.changed()) {//Note for the future: there is probably room to improve this further.
            init();
            prev[source]=-1;
            computeShortestPath(s_start);
        }

        if (last_history_clear != g.historyclears) {
            if (!g.changed() && last_history_clear>=0 && last_history_clear == g.historyclears-1 && history_qhead==g.getPreviousHistorySize()){
                //no information was lost in the history clear
                history_qhead = 0;
                last_history_clear = g.historyclears;
            }else {
                history_qhead = g.historySize();
                last_history_clear = g.historyclears;
                for (int edgeid = 0; edgeid < g.edges(); edgeid++) {
                    if (g.edgeEnabled(edgeid)) {

                    } else {

                    }
                }
            }
        }

        while(s_start!=destination) {
            int m = INF;
            int lowest_suc = -1;
            for (int i = 0; i < g.nIncident(s_start, undirected); i++) {
                if (!g.edgeEnabled(g.incident(s_start, i, undirected).id))
                    continue;
                int edgeID = g.incident(s_start, i, undirected).id;
                int v = g.incident(s_start, i, undirected).node;
                if (v != s_start) {
                    int n = edgeCost(edgeID) + G(v);
                    if (n < m) {
                        m = n;
                        lowest_suc = v;
                    }
                }
            }
            assert(lowest_suc >= 0);
            s_start = lowest_suc;
            k_m = k_m + H(s_last, s_start);
            s_last = s_start;
            for (int i = history_qhead; i < g.historySize(); i++) {
                int edgeid = g.getChange(i).id;
                if (g.getChange(i).addition && g.edgeEnabled(edgeid)) {
                    int v = g.getEdge(edgeid).to;
                    updateVertex(v);
                } else if (!g.getChange(i).addition && !g.edgeEnabled(edgeid)) {
                    int v = g.getEdge(edgeid).to;
                    updateVertex(v);
                }
            }
            computeShortestPath(s_start);
        }
        num_updates++;
        last_modification = g.modifications;
        last_deletion = g.deletions;
        last_addition = g.additions;

        history_qhead = g.historySize();

        last_history_clear = g.historyclears;
        assert(dbg_uptodate());
    }

    bool dbg_path(int to) {
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
    void drawFull() {
        printf("digraph{\n");
        for (int i = 0; i < g.nodes(); i++) {

            if (seen[i]) {
                printf("n%d [fillcolor=blue style=filled]\n", i);
            } else {
                printf("n%d \n", i);
            }

        }

        for (int i = 0; i < g.nodes(); i++) {
            for (int j = 0; j < g.nIncident(i); j++) {
                int id = g.incident(i, j).id;
                int u = g.incident(i, j).node;
                const char * s = "black";
                if (g.edgeEnabled(id))
                    s = "blue";
                else
                    s = "red";

                printf("n%d -> n%d [label=\"v%d\",color=\"%s\"]\n", i, u, id, s);
            }
        }

        printf("}\n");
    }
    bool dbg_uptodate() {

        return true;
    }

    bool connected_unsafe(int t) {
        return t < seen.size() && seen[t];
    }
    bool connected_unchecked(int t) {
        assert(last_modification == g.modifications);
        return connected_unsafe(t);
    }
    bool connected(int t) {
        if (last_modification != g.modifications)
            update();

        assert(dbg_uptodate());

        return seen[t];
    }
    int distance(int t) {
        if (connected(t))
            return 1;
        else
            return -1;
    }
    int distance_unsafe(int t) {
        if (connected_unsafe(t))
            return 1;
        else
            return -1;
    }
    int incomingEdge(int t) {
        assert(t >= 0 && t < prev.size());
        assert(prev[t] >= -1);
        return prev[t];
    }
    int previous(int t) {
        if (prev[t] < 0)
            return -1;
        if (undirected && g.getEdge(incomingEdge(t)).from == t) {
            return g.getEdge(incomingEdge(t)).to;
        }
        assert(g.getEdge(incomingEdge(t)).to == t);
        return g.getEdge(incomingEdge(t)).from;
    }

};
}
;

#endif //MONOSAT_DSTARLITE_H
