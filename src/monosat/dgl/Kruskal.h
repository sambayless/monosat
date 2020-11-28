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
#ifndef KRUSKAL_H_
#define KRUSKAL_H_

#include <vector>
#include "monosat/dgl/alg/Heap.h"
#include "monosat/mtl/Sort.h"
#include "Graph.h"
#include "monosat/core/Config.h"
#include "MinimumSpanningTree.h"
#include "monosat/dgl/alg/DisjointSets.h"
#include <limits>
#include <algorithm>
#include <stdexcept>

namespace dgl {
template<class Status, typename Weight = int>
class Kruskal : public MinimumSpanningTree<Weight> {
public:

    Graph<Weight>& g;

    Status& status;
    int last_modification;
    Weight min_weight;
    int last_addition;
    int last_deletion;
    int history_qhead;

    int last_history_clear;
    bool hasParents;
    Weight INF;
    DisjointSets sets;
    std::vector<int> mst;
    std::vector<int> q;
    std::vector<int> check;
    const int reportPolarity;

    //std::vector<char> old_seen;
    std::vector<bool> in_tree;;
    std::vector<int> parents;
    std::vector<int> parent_edges;

    struct EdgeLt {
        const std::vector<Weight>& edge_weights;

        bool operator()(int x, int y) const{
            return edge_weights[x] < edge_weights[y];
        }

        EdgeLt(const std::vector<Weight>& _edge_weights) :
                edge_weights(_edge_weights){
        }
    };

    alg::Heap<EdgeLt> edge_heap;
    std::vector<int> edge_list;

    std::vector<int> prev;

public:

    int stats_full_updates = 0;
    int stats_fast_updates = 0;
    int stats_fast_failed_updates = 0;
    int stats_skip_deletes = 0;
    int stats_skipped_updates = 0;
    int stats_num_skipable_deletions = 0;
    double mod_percentage = 0;

    double stats_full_update_time = 0;
    double stats_fast_update_time = 0;

    void dbg_printSpanningTree(bool showWeights = true){
#ifdef DEBUG_DGL

#ifdef DEBUG_DGL
        printf("graph{\n");

        for (int i = 0; i < g.adjacency_list.size(); i++) {
            for (int j = 0; j < g.adjacency_list[i].size(); j++) {
                int id = g.adjacency_list[i][j].id;

                int u = g.adjacency_list[i][j].node;
                const char * s = "black";
                if(in_tree[id]){
                    s="green";
                    assert(g.edgeEnabled(id));
                }else if (g.edgeEnabled(id))
                    s = "red";
                else
                    s = "blue";
                if(showWeights){
                    std::stringstream ss;
                    ss<<g.getWeight(id);
                    printf("n%d -- n%d [label=\"v%d w=%s\",color=\"%s\"]\n", i,u, id,ss.str().c_str(), s);
                }else{
                    printf("n%d -- n%d [label=\"v%d\",color=\"%s\"]\n", i, u, id, s);
                }
            }
        }
        printf("}\n");
#endif


#endif

    }

    Kruskal(Graph<Weight>& graph, Status& _status =
    MinimumSpanningTree<Weight>::nullStatus, int _reportPolarity = 0) :
            g(graph), status(_status), last_modification(-1), last_addition(-1), last_deletion(-1), history_qhead(
            0), last_history_clear(0), INF(0), reportPolarity(_reportPolarity), edge_heap(EdgeLt(g.getWeights())){

        mod_percentage = 0.2;

        min_weight = -1;
        hasParents = false;
    }

    void setNodes(int n){
        q.reserve(n);
        check.reserve(n);
        in_tree.resize(g.nEdgeIDs());

        //INF=std::numeric_limits<int>::max();
        sets.AddElements(n);
        parents.resize(n);
        parent_edges.resize(n);
    }

    int64_t num_updates = 0;

    int numUpdates() const{
        return num_updates;
    }

    void update(){

        if(g.outfile()){
            fprintf(g.outfile(), "m\n");
            fflush(g.outfile());
        }

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
        if(last_modification <= 0 || g.changed() || last_history_clear != g.nHistoryClears()){
            INF = 1; //g.nodes()+1;

            for(auto& w : g.getWeights())
                INF += w;
        }
        setNodes(g.nodes());
        min_weight = 0;

        mst.clear();

        if(edge_list.size() < g.nEdgeIDs()){
            edge_list.clear();
            for(int i = 0; i < g.nEdgeIDs(); i++){
                edge_list.push_back(i);

            }
            std::sort(edge_list.begin(), edge_list.end(), EdgeLt(g.getWeights()));
        }
        for(int i = 0; i < in_tree.size(); i++)
            in_tree[i] = false;
        for(int i = 0; i < edge_list.size(); i++){
            int edge_id = edge_list[i];
            if(!g.edgeEnabled(edge_id))
                continue;
            int u = g.getEdge(edge_id).from;
            int v = g.getEdge(edge_id).to;

            int set1 = sets.FindSet(u);
            int set2 = sets.FindSet(v);
            if(set1 != set2){
                assert(g.edgeEnabled(edge_id));
                in_tree[edge_id] = true;
                mst.push_back(edge_id);
                min_weight += g.getWeight(edge_id);
                sets.UnionSets(set1, set2);
                assert(sets.FindSet(u) == sets.FindSet(v));
            }
        }

        status.setMinimumSpanningTree(sets.NumSets() > 1 ? INF : min_weight, sets.NumSets() <= 1);

        for(int i = 0; i < in_tree.size(); i++){
            //Note: for the tree edge detector, polarity is effectively reversed.
            if(reportPolarity < 1 && (!g.edgeEnabled(i) || in_tree[i])){
                status.inMinimumSpanningTree(i, true);
            }else if(reportPolarity > -1 && (g.edgeEnabled(i) && !in_tree[i])){
                status.inMinimumSpanningTree(i, false);
            }
        }

        num_updates++;
        last_modification = g.getCurrentHistory();
        last_deletion = g.nDeletions();
        last_addition = g.nAdditions();

        history_qhead = g.historySize();
        last_history_clear = g.nHistoryClears();

        assert(dbg_uptodate());
    }

    std::vector<int>& getSpanningTree(){
        update();
        return mst;
    }

    int getParent(int node){
        update();
        //kruskals doesn't actually give us the parents, normally. need to construct it on demand, here.
        if(!hasParents)
            buildParents();
        return parents[node];
    }

    int getParentEdge(int node){
        if(getParent(node) != -1)
            return parent_edges[node];
        else
            return -1;
    }

    bool edgeInTree(int edgeid){
        update();
        return in_tree[edgeid];
    }

    bool dbg_mst(){

        return true;
    }

    Weight& weight(){

        update();

        assert(dbg_uptodate());
        if(sets.NumSets() > 1)
            return INF;
        return min_weight;
    }

    Weight& forestWeight(){
        update();
        assert(dbg_uptodate());
        return min_weight;
    }

    int numComponents(){
        update();
        return sets.NumSets();
    }

    int getComponent(int node){
        update();
        return sets.FindSet(node);
    }

    int getRoot(int component = 0){
        update();
        if(!hasParents)
            buildParents();
        int u = 0;
        if(sets.NumSets() > 1)
            u = sets.GetElement(component);

        while(int p = getParent(u) != -1){
            u = p;
        }
        assert(getParent(u) == -1);
        return u;
    }

    bool dbg_uptodate(){
#ifdef DEBUG_DGL
        Weight sumweight = 0;
        in_tree.resize(g.nEdgeIDs());
        for (int i = 0; i < g.edges(); i++) {
            if (in_tree[i]) {
                sumweight += g.getWeight(i);
            }
        }
        assert(sumweight == min_weight || min_weight == INF);

#endif
        return true;
    };
private:

    void buildParents(){
        hasParents = true;
        for(int i = 0; i < parents.size(); i++){
            parents[i] = -1;
            parent_edges[i] = -1;
        }
        q.clear();
        for(int i = 0; i < sets.NumSets(); i++){
            int root = sets.GetElement(i);                //chose an element arbitrarily from the nth set.
            int set = sets.FindSet(root);
            q.push_back(root);
            assert(parents[root] == -1);
            while(q.size()){
                int u = q.back();
                assert(sets.FindSet(u) == set);
                q.pop_back();
                for(int j = 0; j < g.nIncident(u, true); j++){
                    int edge = g.incident(u, j, true).id;
                    int to = g.incident(u, j, true).node;
                    if(edge < 0)
                        continue;
                    if(edge >= in_tree.size())
                        throw std::logic_error("Error in Kruskal");
                    if(in_tree[edge]){

                        if(parents[to] == -1 && to != root){
                            assert(to != root);
                            //int st = sets.FindSet(to);
                            assert(sets.FindSet(to) == set);
                            parents[to] = u;
                            parent_edges[to] = edge;
                            q.push_back(to);
                        }
                    }
                }
            }
            assert(parents[root] == -1);
        }
    }
};
};
#endif
