
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


#ifndef PKTOPOLOGICALSORT_H_
#define PKTOPOLOGICALSORT_H_

//Algorithm PK from "A Dynamic Topological Sort Algorithm for Directed Acyclic Graphs", 2006


#include <vector>
#include "monosat/dgl/alg/Heap.h"
#include "Graph.h"
#include "monosat/core/Config.h"
#include "Reach.h"
#include <algorithm>
#include "TarjansSCC.h"
#include "DFSCycle.h"

namespace dgl {
template<typename Weight, typename Graph = DynamicGraph<Weight>>
class PKToplogicalSort : public Cycle, public DynamicGraphAlgorithm {
public:

    Graph& g;
    DFSCycle<Weight, Graph, true, false> dfs_cycle;

    std::vector<bool> is_strict_scc;
    std::vector<int> strict_sccs;
    int num_strict_sccs = 0;

    int last_modification = -1;
    int last_addition = 0;
    int last_deletion = 0;
    int history_qhead = 0;

    int last_history_clear = 0;
    int alg_id = -1;

    int INF;
    std::vector<bool> edge_in_cycle;
    std::vector<bool> in_cycle;
    std::vector<int> cycle;
    std::vector<int> store_cycle;
    std::vector<int> prev;

    const int reportPolarity;

    //If true, force the (analyzed) graph to remain a DAG.
    //This is only valid if it is known that any discovered cycles will be eliminated.
    //In practice, what will happen is that as soon as a cycle is created on edge addition, that edge addition will be cancelled (but the cycle will be stored in "cycle" for analysis)
    bool force_dag = false;

    std::vector<int> ignore;
    std::vector<int> ord;
    std::vector<bool> visited;
    std::vector<bool> tmp_mark;
    std::vector<int> L;
    std::vector<int> R;
    std::vector<int> l_xy_B;

    std::vector<int> l_xy_F;
    //bool might_not_have_cycle=false;
    bool has_cycle = false;
    bool cycleComputed = false;
    bool has_topo = false;
    int nextOrd = 0;

    int lower_bound = 0;
    int upper_bound = 0;
    std::vector<bool> edge_enabled;

    struct Ord_LT {
        std::vector<int>& ord;

        bool operator()(int a, int b){
            return ord[a] < ord[b];
        }

        Ord_LT(std::vector<int>& ord) : ord(ord){

        }
    } ord_lt;

public:

    void forceDAG() override{
        //If true, force the (analyzed) graph to remain a DAG.
        //This is only valid if it is known that any discovered cycles will be eliminated.
        //In practice, what will happen is that as soon as a cycle is created on edge addition, that edge addition will be cancelled (but the cycle will be stored in "cycle" for analysis)
        force_dag = true;
    }

    PKToplogicalSort(Graph& graph, int _reportPolarity = 0) :
            g(graph), dfs_cycle(g), INF(0), reportPolarity(_reportPolarity), ord_lt(ord){
        alg_id = g.addDynamicAlgorithm(this);

    }

    std::string getName() override{
        return "PKTopologicalSort";
    }

    void setNodes(int n){
        in_cycle.resize(n);
        edge_in_cycle.resize(g.edges());
        visited.resize(n);
        ord.resize(n);//how should this be initialised?
        is_strict_scc.resize(n);
        tmp_mark.resize(n);
        edge_enabled.clear();
        edge_enabled.resize(g.edges());
        INF = g.nodes() + 1;
    }

private:

    bool inSCC(int node){
        assert(has_cycle);

        return in_cycle[node];


    }

    void invalidateSCC(int node){
        assert(has_cycle);
        if(cycle.size()){
#ifdef DEBUG_DGL
            for(int edgeID:cycle){
                int from = g.getEdge(edgeID).from;
                int to = g.getEdge(edgeID).to;
                assert(in_cycle[from]);
                assert(in_cycle[to]);
                assert(edge_in_cycle[edgeID]);
            }

#endif
            if(in_cycle[node]){

                num_strict_sccs = 0;
                for(int edgeID:cycle){
                    int from = g.getEdge(edgeID).from;
                    int to = g.getEdge(edgeID).to;
                    in_cycle[from] = false;
                    in_cycle[to] = false;
                    edge_in_cycle[edgeID] = false;
                }
                cycle.clear();
                cycleComputed = false;
            }
        }else{
            num_strict_sccs = 0;
/*
			int sccID =scc.getComponentUnsafe(node);
			if (is_strict_scc[sccID]){
				is_strict_scc[sccID]=false;
				num_strict_sccs--;
				if(num_strict_sccs==0){
					cycleComputed=false;
				}
			}
			assert(num_strict_sccs>=0);*/
        }
    }

    void updateSCCs(){
        if(has_cycle && !cycleComputed){
            for(int edgeID:cycle){
                int from = g.getEdge(edgeID).from;
                int to = g.getEdge(edgeID).to;
                in_cycle[from] = false;
                in_cycle[to] = false;
                edge_in_cycle[edgeID] = false;
            }
            cycle.clear();
            //search for a cycle
            if(dfs_cycle.hasDirectedCycle()){
                cycle = dfs_cycle.getDirectedCycle();
                for(int edgeID:cycle){
                    int from = g.getEdge(edgeID).from;
                    int to = g.getEdge(edgeID).to;
                    in_cycle[from] = true;
                    in_cycle[to] = true;
                    edge_in_cycle[edgeID] = true;
                }
            }
            has_cycle = cycle.size() > 0;
            num_strict_sccs = has_cycle ? 1 : 0;
            cycleComputed = has_cycle;
            if(has_cycle){
                checkCycle();
            }
        }
        if(has_cycle){
            assert(cycle.size() > 0);
            assert(num_strict_sccs > 0);
            assert(cycleComputed);
        }
    }

    void removeEdge(int edgeID){
        if(!has_cycle){
            //don't need to do anything at all.
            return;
        }else{
            checkCycle();
            //int from = g.getEdge(edgeID).from;
            int to = g.getEdge(edgeID).to;
            if(!cycleComputed){
                updateSCCs();
            }
            if(edge_in_cycle[edgeID]){

                //then we have just broken an existing cycle, and need to check if the graph is still cyclic.
                invalidateSCC(to);
                if(force_dag){
                    has_cycle = false;
                }else{
                    //what to do in this case? for now, first: check if tarjan's scc still has any sccs that are not broken.
                    if(num_strict_sccs > 0){
                        //then we still have at least one scc, somewhere. mark the current cycle invalid in scc.
                    }else{
                        updateSCCs();

                    }
                }
            }else{
                //there is still at least one cycle in this graph.
                return;
            }
        }
    }

    //PK Algorithm:

    void addEdge(int edgeID){
#ifdef DEBUG_DGL
        for(int i = 0;i<visited.size();i++){
            assert(!visited[i]);
        }
#endif
        //once the PK algorithm has a cycle, it is in an invalid state.
        if(has_cycle)
            return;
        int from = g.getEdge(edgeID).from;
        int to = g.getEdge(edgeID).to;
        if (from == to){
            // this edge is a self-loop
            has_cycle = true;
            cycle.clear();
            cycle.push_back(edgeID);
            checkCycle();
            return;
        }

        if(!has_topo){
            topologicalSort();//re-initialized topological sort
            if(has_cycle){
                checkCycle();
                return;
            }
        }

        lower_bound = ord[to];
        upper_bound = ord[from];
        if(lower_bound < upper_bound){
            dfs_forward(to);

            if(!has_cycle){
                dfs_backward(from, has_cycle ? to : -1, edgeID);
                reorder();
#ifdef DEBUG_DGL
                for(int i = 0;i<visited.size();i++){
                    assert(!visited[i]);
                }
#endif
            }else{
                cycle.push_back(edgeID);
                in_cycle[from] = true;
                in_cycle[to] = true;
                edge_in_cycle[edgeID] = true;
                assert(visited[to]);
                checkCycle();

                for(int i = 0; i < l_xy_B.size(); i++){
                    visited[l_xy_B[i]] = false;
                }
                for(int i = 0; i < l_xy_F.size(); i++){
                    visited[l_xy_F[i]] = false;
                }
#ifdef DEBUG_DGL
                for(int i = 0;i<visited.size();i++){
                    assert(!visited[i]);
                }
#endif
                l_xy_F.clear();
                l_xy_B.clear();
                /*if(force_dag){
                    has_cycle=false;//ignore the edge that caused a cycle to be introduced, because it is guaranteed that any edges creating cycles will be removed.
                }else{*/
                has_topo = false;//topological sort is no longer valid.
                //}
            }
        }

    }

    void checkCycle(){
#ifdef DEBUG_DGL
        std::vector<bool> seen;
        seen.clear();
        seen.resize(g.nodes(),false);
        //first, check that the in_cycle vector is correct
        for(int edgeID:cycle){
            int f = g.getEdge(edgeID).from;
            int t = g.getEdge(edgeID).to;
            assert(in_cycle[f]);
            assert(in_cycle[t]);
            assert(edge_in_cycle[edgeID]);
            seen[f]=true;
            seen[t]=true;
        }
        for(int i = 0;i<in_cycle.size();i++){
            if(in_cycle[i]){
                assert(seen[i]);
            }
        }

        std::vector<int> count;
        count.clear();
        count.resize(g.nodes(),0);
        for(int edgeID:cycle){
            int f = g.getEdge(edgeID).from;
            int t = g.getEdge(edgeID).to;
            assert(in_cycle[f]);
            assert(in_cycle[t]);
            count[f]= count[f]-1;
            count[t]= count[t] +1;
        }

        for(int i:count){
            if(i!=0){
                throw std::runtime_error("Error in PKTopologicalSort");
            }
            assert(i==0);
        }

#endif
    }

    bool dfs_forward(int n){
        assert(!has_cycle);
        visited[n] = true;
        l_xy_F.push_back(n);
        for(int i = 0; !has_cycle && i < g.nIncident(n); i++){
            int edgeID = g.incident(n, i).id;
            if(edge_enabled[edgeID]){
                int w = g.incident(n, i).node;
                if(ord[w] == upper_bound){
                    //cycle detected.
                    has_cycle = true;
                    //assert(cycleComputed==false);//cycle will be computed next time it is relevant.
                    cycleComputed = true;
                    num_strict_sccs = 1;
                    assert(cycle.size() == 0);
                    cycle.push_back(edgeID);
                    in_cycle[n] = true;
                    in_cycle[w] = true;
                    edge_in_cycle[edgeID] = true;
                    return false;
                }
                if(!visited[w] && ord[w] < upper_bound){
                    if(!dfs_forward(w)){
                        cycle.push_back(edgeID);
                        in_cycle[n] = true;
                        in_cycle[w] = true;
                        edge_in_cycle[edgeID] = true;
                        return false;
                    }
                }
            }
        }
        return true;
    }


    void dfs_backward(int n, int cycleNode = -1, int skipEdge = -1){

        visited[n] = true;
        l_xy_B.push_back(n);
        for(int i = 0; i < g.nIncoming(n); i++){
            int edgeID = g.incoming(n, i).id;
            if(edgeID != skipEdge && edge_enabled[edgeID]){
                int w = g.incoming(n, i).node;
                if(!visited[w] && w == cycleNode){
                    visited[w] = true;
                    l_xy_B.push_back(w);
                    if(has_cycle){

                        cycle.push_back(edgeID);
                        in_cycle[n] = true;
                        in_cycle[w] = true;
                        edge_in_cycle[edgeID] = true;
                    }
                }else if(!visited[w] && lower_bound < ord[w]){
                    if(has_cycle){
                        cycle.push_back(edgeID);
                        edge_in_cycle[edgeID] = true;
                        in_cycle[n] = true;
                        in_cycle[w] = true;
                    }
                    dfs_backward(w);
                }
            }
        }
    }

    void merge(std::vector<int>& left, std::vector<int>& right, std::vector<int>& store){
        store.clear();
        int lpos = 0;
        int rpos = 0;
        while(lpos < left.size() && rpos < right.size()){
            if(left[lpos] <= right[rpos]){
                store.push_back(left[lpos++]);
            }else{
                store.push_back(right[rpos++]);
            }
        }
        while(lpos < left.size()){
            store.push_back(left[lpos++]);
        }
        while(rpos < right.size()){
            store.push_back(right[rpos++]);
        }
    }

    void reorder(){
        assert(!has_cycle);
        //sort l_xy arrays with respect to their current ordinal
        std::sort(l_xy_B.begin(), l_xy_B.end(), ord_lt);
        std::sort(l_xy_F.begin(), l_xy_F.end(), ord_lt);


        L.clear();
        R.clear();
        for(int i = 0; i < l_xy_B.size(); i++){
            int w = l_xy_B[i];
            l_xy_B[i] = ord[w];
            visited[w] = false;
            L.push_back(w);
        }
        for(int i = 0; i < l_xy_F.size(); i++){
            int w = l_xy_F[i];
            l_xy_F[i] = ord[w];
            visited[w] = false;
            L.push_back(w);
        }

        merge(l_xy_B, l_xy_F, R);
        for(int i = 0; i < L.size(); i++){
            ord[L[i]] = R[i];
        }
        l_xy_F.clear();
        l_xy_B.clear();
    }

    void dbg_check_topo(){
#ifdef DEBUG_DGL
        if(!has_topo){
            return;
        }
        //assert ords are unique
        static std::vector<bool> seen;
        seen.clear();
        seen.resize(g.nodes());
        for(int i = 0;i<ord.size();i++){
            int o = ord[i];
            if(seen[o]){
                assert(false);
            }
            seen[o]=true;
        }

        for(int edgeID = 0;edgeID<g.edges();edgeID++){
            if(g.edgeEnabled(edgeID)){
                int from = g.getEdge(edgeID).from;
                int to = g.getEdge(edgeID).to;

                int ordF = ord[from];
                int ordT = ord[to];
                assert(ordF<ordT);
            }
        }
#endif
    }

public:
    void update() override{
        if(last_modification > 0 && g.getCurrentHistory() == last_modification){
            stats_skipped_updates++;
            return;
        }
        if(last_modification <= 0 || g.nHistoryClears() != last_history_clear || g.changed()){
            setNodes(g.nodes());
            cycleComputed = false;
            has_cycle = false;
            has_topo = false;
            cycle.clear();
            history_qhead = g.historySize();
            stats_history_clears++;
            for(int edgeID = 0; edgeID < g.edges(); edgeID++){
                edge_enabled[edgeID] = g.edgeEnabled(edgeID);
            }
            if(!topologicalSort()){
                updateSCCs();
            }

        }

        for(int i = history_qhead; i < g.historySize(); i++){
            int edgeID = g.getChange(i).id;

            if(g.getChange(i).addition && g.edgeEnabled(edgeID) && !edge_enabled[edgeID]){
                edge_enabled[edgeID] = true;
                addEdge(edgeID);
            }else if(!g.getChange(i).addition && !g.edgeEnabled(edgeID) && edge_enabled[edgeID]){
                edge_enabled[edgeID] = false;
                removeEdge(edgeID);
            }
        }
        dbg_check_topo();
        last_modification = g.getCurrentHistory();
        last_deletion = g.nDeletions();
        last_addition = g.nAdditions();

        history_qhead = g.historySize();
        g.updateAlgorithmHistory(this, alg_id, history_qhead);
        last_history_clear = g.nHistoryClears();

    }

    void updateHistory() override{
        update();
    }

private:
    int cycle_node = -1;

    bool sortVisit(int node){

        if(tmp_mark[node]){
            has_cycle = true;
            cycle_node = node;
            tmp_mark[node] = false;
            //assert(cycleComputed==false);//cycle will be computed next time it is relevant.
            cycleComputed = true;
            num_strict_sccs = 1;
            assert(cycle.size() == 0);


            return false;
        }
        if(!visited[node]){
            tmp_mark[node] = true;
            for(int j = 0; j < g.nIncident(node); j++){
                int edgeID = g.incident(node, j).id;
                if(g.edgeEnabled(edgeID)){
                    if(!sortVisit(g.incident(node, j).node)){
                        assert(has_cycle);
                        tmp_mark[node] = false;
                        if(cycle_node >= 0){
                            if(node == cycle_node)
                                cycle_node = -1;
                            cycle.push_back(edgeID);
                            edge_in_cycle[edgeID] = true;
                            in_cycle[node] = true;
                            in_cycle[g.incident(node, j).node] = true;
                        }
                        return false;
                    }
                }
            }
            visited[node] = true;
            tmp_mark[node] = false;
            ord[node] = --nextOrd;
        }
        return true;
    }

    //topologically sorts the graph, or returns false if it has a directed cycle.
    bool topologicalSort(){
        if(has_topo)
            return true;
        has_topo = true;
        nextOrd = g.nodes();
        //from wikipedia's pseudocode: http://en.wikipedia.org/wiki/Topological_sorting
        L.clear();
#ifdef DEBUG_DGL
        for(int i = 0;i<tmp_mark.size();i++){
            assert(!tmp_mark[i]);
            assert(!visited[i]);
        }
#endif
        cycle_node = -1;
        for(int n = 0; n < g.nodes(); n++){
            if(!visited[n]){
                if(!sortVisit(n)){
                    assert(has_cycle);
                    has_topo = false;
                    for(int i = 0; i < visited.size(); i++){
                        visited[i] = false;
                    }
                    return false;
                }
            }
        }
        assert(cycle_node == -1);
        for(int i = 0; i < visited.size(); i++){
            visited[i] = false;
        }
        dbg_check_topo();

        return has_topo;
    }

public:
    bool hasDirectedCycle() override{
        update();
        return has_cycle;
    }

    //get _any_ directed cycle from this graph (must be cyclic)
    std::vector<int>& getDirectedCycle() override{
        update();
        checkCycle();


        return cycle;

        /*	if(!hasDirectedCycle()){
                assert(cycle.size()==0);
                return cycle;
            }

            store_cycle.clear();
            for(int id:strict_sccs){
                int element = scc.getElement(id);
                //in the future, it might be a good idea to put this into the cached 'cycle'...
                store_cycle.clear();

                if(store_cycle.size()>1){
                    return store_cycle;
                }
                store_cycle.clear();
                if(cycle.size()>1){
                    for(int n:cycle){
                        in_cycle[n]=true;
                    }
                    return cycle;
                }
            }
            return store_cycle;*/
    }

    bool hasUndirectedCycle() override{
        assert(false);
        return false;//not implemented
    }

    std::vector<int>& getUndirectedCycle() override{
        assert(false);
        return ignore;
    }


};
};


#endif /* PKTOPOLOGICALSORT_H_ */
