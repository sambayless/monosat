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

#ifndef DINICS_LINKCUT_H
#define DINICS_LINKCUT_H

#include "Graph.h"
#include "MaxFlow.h"
#include <vector>
#include "monosat/core/Config.h"
#include "EdmondsKarpAdj.h"
#include "monosat/dgl/alg/LinkCutCost.h"
#include <algorithm>
#include <climits>
//Implementation of the Link/Cut Tree (aka Dynamic Tree) max-flow algorithm from Sleator and Tarjan, 1983).
namespace dgl {
template<typename Weight>
class DinitzLinkCut : public MaxFlow<Weight> {

public:

    std::vector<Weight> F;

    struct LocalEdge {
        int from;
        int id;
        bool backward = false;

        LocalEdge(int from = -1, int id = -1, bool backward = false) :
                from(from), id(id), backward(backward){

        }
    };

    bool dinics_recursive = false;
    Weight curflow;
    int last_modification;
    int last_deletion;
    int last_addition;

    int history_qhead;
    int last_history_clear;
    std::vector<int> dist;
    std::vector<int> pos;    //position in the combined forward and backward adjacency list of each node in the DFS.
    std::vector<bool> changed;
    Graph<Weight>& g;

    int source = -1;
    int sink = -1;
    Weight INF;
    int src = -1;
    int dst = -1;
    struct ParentEdge {
        bool backward;
        int edgeID;
    };
    std::vector<ParentEdge> parentEdge;

    LinkCutCost<Weight> forest;
    std::vector<int> Q;
    struct Link {
        int u;
        int v;
        bool backward :1;
        int edgeID :31;
    };
    std::vector<Link> toLink;


    int64_t stats_augmenting_rounds = 0;
    int64_t stats_rounds = 0;
    int64_t stats_backtracks = 0;
    int64_t stats_avoided_backtracks = 0;

    void printStats(){
        printf("Dinics Link Cut:\n");
        printf("Rounds: %" PRId64 ", Augmenting Rounds: %" PRId64 "\n", stats_rounds, stats_augmenting_rounds);
        printf("Backtracks %" PRId64 " (%" PRId64 " avoided)\n", stats_backtracks, stats_avoided_backtracks);
    }

public:
    DinitzLinkCut(Graph<Weight>& _g, int source = -1, int sink = -1) :
            g(_g), source(source), sink(sink), INF(0xF0F0F0){
        curflow = 0;
        last_modification = -1;
        last_deletion = -1;
        last_addition = -1;

        history_qhead = -1;
        last_history_clear = -1;

    }

    int getSource() const{
        return source;
    }

    int getSink() const{
        return sink;
    }

    void dbg_print_graph(int from, int to){

#ifdef DEBUG_DGL
        /*	return;

        printf("digraph{\n");
        for (int i = 0; i < g.nodes(); i++) {
            if (i == from) {
                printf("n%d [label=\"From\", style=filled, fillcolor=blue]\n", i);
            } else if (i == to) {
                printf("n%d [label=\"To\", style=filled, fillcolor=red]\n", i);
            } else if (dist[i] == -1) {
                printf("n%d [style=filled, fillcolor=gray]\n", i);
            } else
                printf("n%d\n", i);
        }

        for (int i = 0; i < g.edges(); i++) {
            if (g.edgeEnabled(i)) {
                auto  e = g.getEdge(i);
                const char * s = "black";

                bool link = false;
                bool backward = false;
                for (auto & e : toLink) {
                    if (e.edgeID == i) {
                        link = true;
                        backward = e.backward;
                        if (backward) {
                            s = "orange";
                        } else {
                            s = "green";
                        }
                        break;
                    }
                }
                int hasParent = false;
                bool backwardParent = false;
                if (parentEdge[e.from].edgeID == i) {
                    hasParent = true;
                    s = "red";
                } else if (parentEdge[e.to].edgeID == i) {
                    s = "blue";
                    //hasParent=true;
                    backwardParent = true;
                }

                //assert(!(hasParent && link));

                if (dist[e.to] == dist[e.from] + 1) {
                    s = "blue";
                }

                std::cout << "n" << e.from << " -> n" << e.to << " [label=\"" << i << ": " << F[i] << "/" << g.getWeight(i)
                        << "\" color=\"" << s << "\"]\n";
            }
        }

        printf("}\n");*/
#endif

    }

    bool buildLevelGraph(int src, int dst){

        dist.clear();
        dist.resize(g.nodes(), -1);
        dist[src] = 0;
        Q.push_back(src);
        //Build the level graph using a simple BFS
        for(int i = 0; i < Q.size(); i++){
            int u = Q[i];
            for(int j = 0; j < g.nIncident(u); j++){
                int edgeID = g.incident(u, j).id;
                if(!g.edgeEnabled(edgeID))
                    continue;
                int v = g.incident(u, j).node;
                if(dist[v] < 0 && F[edgeID] < g.getWeight(edgeID)){
                    dist[v] = dist[u] + 1;
                    Q.push_back(v);
                }
            }
            for(int j = 0; j < g.nIncoming(u); j++){
                int edgeID = g.incoming(u, j).id;
                if(!g.edgeEnabled(edgeID))
                    continue;
                int v = g.incoming(u, j).node;
                //this is a backward edge, so it has capacity exactly if the forward edge has flow
                if(dist[v] < 0 && F[edgeID]){
                    dist[v] = dist[u] + 1;
                    Q.push_back(v);
                }
            }
        }
        Q.clear();

        return dist[dst] >= 0;
    }

    bool dbg_hasLink(int u){
#ifdef DEBUG_DGL
        for (auto e : toLink) {
            if (e.u == u)
                return true;
        }
#endif
        return false;
    }

    bool dbg_isLinkRoot(int v){
#ifdef DEBUG_DGL
        int u = src;
        int i = 0;
        while (true) {
            u = forest.findRoot(u);
            if (i < toLink.size()) {
                auto & e = toLink[i++];
                assert(e.u == u);
                u = e.v;
            } else {
                break;
            }
        }
        /*if(u!=v){
         dbg_print_graph(src,dst);
         }*/
        assert(u == v);
        return u == v;
#endif
        return true;
    }

    Weight findAugmentingPath(int src, int dst){
        dbg_print_graph(src, dst);
        toLink.clear();
        Weight f = 0;
        int u = forest.findRoot(src);
        while(true){
            bool found = true;
            bool foundPath = false;
            dbg_print_graph(src, dst);

            while(found){
                found = false;
                assert(dbg_isLinkRoot(u));
                if(u == dst){
                    foundPath = true;
                    break;
                }else{
                    //extend the path
                    for(; pos[u] < g.nIncident(u); pos[u]++){
                        auto& edge = g.incident(u, pos[u]);
                        int edgeID = edge.id;
                        int v = edge.node;
                        if((dist[v] != dist[u] + 1) || !g.edgeEnabled(edgeID))
                            continue;
                        if(F[edgeID] < g.getWeight(edgeID)){

                            assert(!dbg_hasLink(u));
                            toLink.push_back({u, v, false, edgeID});
                            u = forest.findRoot(v);
                            found = true;
                            break;
                        }
                    }
                    if(!found){
                        for(; pos[u] - g.nIncident(u) < g.nIncoming(u); pos[u]++){
                            //auto & edge = g.inverted_adjacency[u][pos[u]-g.nIncident(u)];
                            auto& edge = g.incoming(u, pos[u] - g.nIncident(u));
                            int edgeID = edge.id;
                            int v = edge.node;
                            if((dist[v] != dist[u] + 1) || !g.edgeEnabled(edgeID))
                                continue;

                            //these are backwards edges, which have capacity exactly if the forward edge has non-zero flow
                            if(F[edgeID] != 0){

                                assert(!dbg_hasLink(u));
                                toLink.push_back({u, v, true, edgeID});
                                u = forest.findRoot(v);
                                found = true;
                                break;
                            }
                        }
                    }
                }

            }
            dbg_isLinkRoot(u);
            dbg_print_graph(src, dst);

            assert(!found);
            if(foundPath){
                stats_augmenting_rounds++;
                dbg_print_graph(src, dst);
                //ok, link the path together now that we know there is a path.
                for(int i = 0; i < toLink.size(); i++){
                    auto& e = toLink[i];
                    int u = e.u;
                    int v = e.v;
                    int edgeID = e.edgeID;
                    bool backward = e.backward;
                    assert(edgeID >= 0);
                    assert(forest.findRoot(src) == u);
                    assert(g.getWeight(edgeID) >= F[edgeID]);
                    assert(F[edgeID] >= 0);
                    if(backward){
                        forest.link(u, v, F[edgeID]);
                        parentEdge[u].edgeID = edgeID;
                        parentEdge[u].backward = true;
                    }else{
                        forest.link(u, v, g.getWeight(edgeID) - F[edgeID]);
                        parentEdge[u].backward = false;
                        parentEdge[u].edgeID = edgeID;
                    }
                    assert(forest.findRoot(src) == forest.findRoot(v));
                }
                toLink.clear();
                dbg_print_graph(src, dst);
                assert(forest.findRoot(src) == dst);

                //found s-t augmenting path.
                //Find the minimum capacity of any edge from src to dst on that path
                Weight c = forest.minCost(src);
                //subtract that capacity from all edges on the path
                f += c;
                forest.updateCostOfPathToRoot(src, -c);
                //delete edges with no remaining capacity
                Weight minC = 0;
                while(minC == 0){
                    int u = forest.ancecstorFindMin(src);
                    assert(forest.getCost(u) == minC);
                    int edgeID = parentEdge[u].edgeID;
                    assert(edgeID >= 0);
                    if(!parentEdge[u].backward){
                        assert(F[edgeID] + c <= g.getWeight(edgeID));
                        F[edgeID] += c;
                    }else{
                        assert(c <= F[edgeID]);
                        F[edgeID] -= c;
                    }
                    forest.cut(u);
                    parentEdge[u].edgeID = -1;
                    minC = forest.minCost(src);
                }
                toLink.clear();
                u = forest.findRoot(src);
                dbg_print_graph(src, dst);

            }else{
                //couldn't find any s-t augmenting path.
                //remove the vertex from the graph
                if(u == src){
                    //done - no paths remain.
                    //Clean up the state of the tree:
                    for(int u = 0; u < g.nodes(); u++){
                        if(parentEdge[u].edgeID >= 0){
                            assert(!forest.isRoot(u));
                            Weight c = forest.getCost(u);

                            int edgeID = parentEdge[u].edgeID;
                            assert(c <= g.getWeight(edgeID));
                            if(!parentEdge[u].backward){
                                F[edgeID] = g.getWeight(edgeID) - c;
                            }else{
                                F[edgeID] = c;
                            }

                            forest.cut(u);
                            parentEdge[u].edgeID = -1;
                        }
                        forest.undeleteNode(u);
                    }

                    toLink.clear();

                    break;
                }else{
                    dist[u] = -1;                         //prevents u from being explored again.
                    if(toLink.size() && u == toLink.back().v){
                        stats_avoided_backtracks++;
                        u = toLink.back().u;
                        toLink.pop_back();
                    }else{
                        stats_backtracks++;
                        forest.removeNode(u);
                        for(int i = 0; i < g.nIncident(u); i++){
                            auto& edge = g.incident(u, i, false);
                            int edgeID = edge.id;
                            int v = edge.node;
                            if(!g.edgeEnabled(edgeID))
                                continue;
                            if(parentEdge[v].edgeID == edgeID){
                                //need to remember the remaining flow on this edge...
                                assert(F[edgeID] > 0);                         //else this edge wouldn't be in the tree
                                Weight residual_capacity = forest.getCost(v);
                                F[edgeID] = residual_capacity;
                                assert(F[edgeID] >= 0);
                                forest.cut(v);                         //this is a backward edge into u
                                parentEdge[v].edgeID = -1;
                            }
                        }
                        for(int i = 0; i < g.nIncoming(u); i++){
                            auto& edge = g.incoming(u, i, false);
                            int edgeID = edge.id;
                            int v = edge.node;
                            if(!g.edgeEnabled(edgeID))
                                continue;
                            if(parentEdge[v].edgeID == edgeID){
                                Weight residual_capacity = forest.getCost(v);
                                F[edgeID] = g.getWeight(edgeID) - residual_capacity;
                                assert(F[edgeID] >= 0);
                                forest.cut(v);
                                parentEdge[v].edgeID = -1;
                            }
                        }
                        if(toLink.size()){
                            u = forest.findRoot(toLink.back().v);                         //should this be u or v?
                        }else{
                            u = forest.findRoot(src);
                        }
                    }
                }
            }
        }

        return f;
    }

    Weight findAugmentingPath_pointers(int src, int dst){
        Weight f = 0;
        while(true){

            bool found = true;
            bool foundPath = false;

            int u = forest.findRoot(src);
            while(found){
                found = false;

                if(u == dst){
                    foundPath = true;
                    break;
                }else{
                    //extend the path

                    for(; pos[u] < g.nIncident(u); pos[u]++){
                        auto& edge = g.incident(u, pos[u]);
                        int edgeID = edge.id;
                        int v = edge.node;
                        if((dist[v] != dist[u] + 1) || !g.edgeEnabled(edgeID))
                            continue;

                        if(F[edgeID] < g.getWeight(edgeID)){
                            forest.link(u, v, g.getWeight(edgeID) - F[edgeID]);
                            parentEdge[u].backward = false;
                            parentEdge[u].edgeID = edgeID;
                            assert(forest.findRoot(u) == v);
                            u = forest.findRoot(v);
                            found = true;
                            break;
                        }
                    }
                    if(!found){
                        for(; pos[u] - g.nIncident(u) < g.nIncoming(u); pos[u]++){

                            auto& edge = g.incoming(u, pos[u] - g.nIncident(u));
                            int edgeID = edge.id;
                            int v = edge.node;
                            if((dist[v] != dist[u] + 1) || !g.edgeEnabled(edgeID))
                                continue;

                            //these are backwards edges, which have capacity exactly if the forward edge has non-zero flow
                            if(F[edgeID]){
                                forest.link(u, v, F[edgeID]);
                                parentEdge[u].backward = true;
                                parentEdge[u].edgeID = edgeID;
                                assert(forest.findRoot(u) == v);
                                u = forest.findRoot(v);
                                found = true;
                                break;
                            }
                        }
                    }
                }

            }

            assert(!found);
            if(foundPath){
                stats_augmenting_rounds++;
                //found s-t augmenting path.
                //Find the minimum capacity of any edge from src to dst on that path
                int c = forest.minCost(src);
                //subtract that capacity from all edges on the path
                f += c;
                forest.updateCostOfPathToRoot(src, -c);
                //delete edges with no remaining capacity
                Weight minC = 0;
                while(minC == 0){

                    int u = forest.ancecstorFindMin(src);
                    assert(forest.getCost(u) == minC);
                    int edgeID = parentEdge[u].edgeID;
                    if(!parentEdge[u].backward){
                        assert(F[edgeID] + c <= g.getWeight(edgeID));
                        F[edgeID] += c;
                    }else{

                        assert(c <= F[edgeID]);
                        F[edgeID] -= c;
                    }
                    forest.cut(u);
                    parentEdge[u].edgeID = -1;
                    minC = forest.minCost(src);
                }

            }else{
                //couldn't find any s-t augmenting path.
                //remove the vertex from the graph
                if(u == src){
                    //done - no paths remain.
                    //Clean up the state of the tree:
                    for(int u = 0; u < g.nodes(); u++){
                        if(parentEdge[u].edgeID >= 0){
                            Weight c = forest.getCost(u);
                            int edgeID = parentEdge[u].edgeID;
                            if(!parentEdge[u].backward){
                                F[edgeID] = g.getWeight(edgeID) - c;
                            }else{

                                F[edgeID] = c;
                            }

                            forest.cut(u);
                            parentEdge[u].edgeID = -1;
                        }
                    }
                    break;
                }else{
                    dist[u] = -1;
                    for(int i = 0; i < g.nIncident(u); i++){
                        auto& edge = g.incident(u, i);
                        int edgeID = edge.id;
                        int v = edge.node;
                        if(!g.edgeEnabled(edgeID))
                            continue;
                        if(parentEdge[v].edgeID == edgeID){
                            //need to remember the remaining flow on this edge...
                            assert(F[edgeID] > 0);                     //else this edge wouldn't be in the tree
                            Weight residual_capacity = forest.getCost(v);
                            F[edgeID] = residual_capacity;
                            assert(F[edgeID] >= 0);
                            forest.cut(v);                     //this is a backward edge into u
                            parentEdge[v].edgeID = -1;
                        }
                    }
                    for(int i = 0; i < g.nIncoming(u); i++){
                        auto& edge = g.incoming(u, i);
                        int edgeID = edge.id;
                        int v = edge.node;
                        if(!g.edgeEnabled(edgeID))
                            continue;
                        if(parentEdge[v].edgeID == edgeID){
                            Weight residual_capacity = forest.getCost(v);
                            F[edgeID] = g.getWeight(edgeID) - residual_capacity;
                            assert(F[edgeID] >= 0);
                            forest.cut(v);
                            parentEdge[v].edgeID = -1;
                        }
                    }
                }
            }

        }

        return f;
    }

    Weight findAugmentingPath_original(int src, int dst){

        Weight f = 0;
        while(true){
            stats_augmenting_rounds++;
            bool found = true;
            bool foundPath = false;

            int u = forest.findRoot(src);
            while(found){
                found = false;

                if(u == dst){
                    foundPath = true;
                    break;
                }else{
                    //extend the path

                    for(int i = 0; i < g.nIncident(u); i++){
                        auto& edge = g.incident(u, i);
                        int edgeID = edge.id;
                        int v = edge.node;
                        if((dist[v] != dist[u] + 1) || !g.edgeEnabled(edgeID))
                            continue;

                        if(F[edgeID] < g.getWeight(edgeID)){
                            forest.link(u, v, g.getWeight(edgeID) - F[edgeID]);

                            parentEdge[u].backward = false;
                            parentEdge[u].edgeID = edgeID;
                            assert(forest.findRoot(u) == v);
                            found = true;
                            break;
                        }
                    }
                    if(!found){
                        for(int i = 0; i < g.nIncoming(u); i++){
                            auto& edge = g.incoming(u, i);
                            int edgeID = edge.id;
                            int v = edge.node;
                            if((dist[v] != dist[u] + 1) || !g.edgeEnabled(edgeID))
                                continue;

                            //these are backwards edges, which have capacity exactly if the forward edge has non-zero flow
                            if(F[edgeID]){
                                forest.link(u, v, F[edgeID]);
                                parentEdge[u].backward = true;
                                parentEdge[u].edgeID = edgeID;
                                assert(forest.findRoot(u) == v);
                                found = true;
                                break;
                            }
                        }
                    }
                }
                u = forest.findRoot(src);
            }
            assert(!found);
            if(foundPath){
                //found s-t augmenting path.
                //Find the minimum capacity of any edge from src to dst on that path
                Weight c = forest.minCost(src);
                //subtract that capacity from all edges on the path
                f += c;
                forest.updateCostOfPathToRoot(src, -c);
                //delete edges with no remaining capacity
                Weight minC = 0;
                while(minC == 0){

                    int u = forest.ancecstorFindMin(src);
                    assert(forest.getCost(u) == minC);
                    int edgeID = parentEdge[u].edgeID;
                    if(!parentEdge[u].backward){
                        assert(F[edgeID] + c <= g.getWeight(edgeID));
                        F[edgeID] += c;
                    }else{
                        assert(c <= F[edgeID]);
                        F[edgeID] -= c;
                    }
                    forest.cut(u);
                    parentEdge[u].edgeID = -1;
                    minC = forest.minCost(src);
                }

            }else{
                //couldn't find any s-t augmenting path.
                //remove the vertex from the graph
                if(u == src){
                    //done - no paths remain.
                    //Clean up the state of the tree:
                    for(int u = 0; u < g.nodes(); u++){
                        if(parentEdge[u].edgeID >= 0){
                            Weight c = forest.getCost(u);
                            int edgeID = parentEdge[u].edgeID;
                            if(!parentEdge[u].edgeID){
                                // if(tree_edges[edgeID].in_tree){
                                F[edgeID] = g.getWeight(edgeID) - c;
                            }else{
                                //assert(tree_edges[edgeID].in_tree_backward);
                                F[edgeID] = c;
                            }

                            forest.cut(u);
                            parentEdge[u].edgeID = -1;
                        }
                    }

                    break;
                }else{
                    dist[u] = -1;

                    for(int i = 0; i < g.nIncident(u); i++){
                        auto& edge = g.incident(u, i);
                        int edgeID = edge.id;
                        int v = edge.node;
                        if(!g.edgeEnabled(edgeID))
                            continue;
                        if(parentEdge[v].edgeID == edgeID){
                            //need to remember the remaining flow on this edge...
                            assert(F[edgeID] > 0);                        //else this edge wouldn't be in the tree
                            Weight residual_capacity = forest.getCost(v);
                            F[edgeID] = residual_capacity;
                            assert(F[edgeID] >= 0);
                            forest.cut(v);                        //this is a backward edge into u
                            parentEdge[v].edgeID = -1;
                        }
                    }
                    for(int i = 0; i < g.nIncoming(u); i++){
                        auto& edge = g.incoming(u, i);
                        int edgeID = edge.id;
                        int v = edge.node;
                        if(!g.edgeEnabled(edgeID))
                            continue;
                        if(parentEdge[v].edgeID == edgeID){
                            Weight residual_capacity = forest.getCost(v);
                            F[edgeID] = g.getWeight(edgeID) - residual_capacity;
                            assert(F[edgeID] >= 0);
                            forest.cut(v);
                            parentEdge[v].edgeID = -1;
                        }
                    }
                }
            }

        }

        return f;
    }

    int64_t num_updates = 0;

    int numUpdates() const{
        return num_updates;
    }

    const Weight update(){
        return maxFlow(source, sink);
    }

    std::vector<int> changed_edges;

    std::vector<int>& getChangedEdges(){
        return changed_edges;
    }

    void clearChangedEdges(){
        for(int edgeID : changed_edges){
            assert(changed[edgeID]);
            changed[edgeID] = false;
        }
        changed_edges.clear();
    }

private:

    void markChanged(int edgeID){
        if(!changed[edgeID]){
            changed[edgeID] = true;
            changed_edges.push_back(edgeID);
        }
    }

public:
    void setSource(int s){
        if(source == s){
            return;
        }
        source = s;
        last_modification = g.getCurrentHistory() - 1;
    }

    void setSink(int t){
        if(sink == t){
            return;
        }
        sink = t;
        last_modification = g.getCurrentHistory() - 1;
    }

    const Weight maxFlow(int s, int t){
        Weight f = 0;
        if(g.outfile()){
            fprintf(g.outfile(), "f %d %d\n", s, t);
            fflush(g.outfile());
        }

        if(last_modification > 0 && g.getCurrentHistory() == last_modification){
            return curflow;
        }
        src = s;
        dst = t;
        F.clear();
        F.resize(g.nEdgeIDs());
        changed.resize(g.nEdgeIDs());
        dist.clear();
        dist.resize(g.nodes());
        f = 0;
        forest.reset();
        while(forest.nNodes() < g.nodes())
            forest.addNode();
        pos.clear();
        pos.resize(g.nodes());

        parentEdge.clear();
        parentEdge.resize(g.nodes(), {false, -1});

        while(buildLevelGraph(s, t)){
            stats_rounds++;

            dbg_print_graph(s, t);
            for(int i = 0; i < pos.size(); i++)
                pos[i] = 0;

            for(int i = 0; i < parentEdge.size(); i++){
                parentEdge[i] = {false, -1};
            }

            if(dinics_recursive){
                Weight delta = findAugmentingPath_original(src, dst);
                f += delta;
            }else{
                Weight delta = findAugmentingPath(src, dst);
                f += delta;
            }
            dbg_print_graph(s, t);

        }

#ifdef DEBUG_MAXFLOW
        Weight expected_flow = ek.maxFlow(s,t);
#endif

#ifdef DEBUG_MAXFLOW
        assert(f==expected_flow);
#endif

        curflow = f;
        num_updates++;
        last_modification = g.getCurrentHistory();
        last_deletion = g.nDeletions();
        last_addition = g.nAdditions();

        history_qhead = g.historySize();
        last_history_clear = g.nHistoryClears();
        return f;
    }

    std::vector<bool> seen;
    std::vector<bool> visited;

    const Weight minCut(std::vector<MaxFlowEdge>& cut){
        return minCut(source, sink, cut);
    }

    const Weight minCut(int s, int t, std::vector<dgl::MaxFlowEdge>& cut){
        const Weight f = maxFlow(s, t);
        //ok, now find the cut
        Q.clear();
        Q.push_back(s);
        seen.clear();
        seen.resize(g.nodes());
        seen[s] = true;
        cut.clear();
        //explore the residual graph
        for(int j = 0; j < Q.size(); j++){
            int u = Q[j];

            for(int i = 0; i < g.nIncident(u); i++){
                if(!g.edgeEnabled(g.incident(u, i).id))
                    continue;
                int v = g.incident(u, i).node;
                int id = g.incident(u, i).id;
                if(g.getWeight(id) - F[id] == 0){
                    cut.push_back(MaxFlowEdge{u, v, id});                        //potential element of the cut
                }else if(!seen[v]){
                    Q.push_back(v);
                    seen[v] = true;
                }
            }
            for(int i = 0; i < g.nIncoming(u); i++){
                if(!g.edgeEnabled(g.incoming(u, i).id))
                    continue;
                int v = g.incoming(u, i).node;
                int id = g.incoming(u, i).id;
                if(F[id] == 0){

                }else if(!seen[v]){
                    Q.push_back(v);
                    seen[v] = true;
                }
            }
        }
        //Now keep only the edges from a seen vertex to an unseen vertex
        int i, j = 0;
        for(i = 0; i < cut.size(); i++){
            if(!seen[cut[i].v] && seen[cut[i].u]){
                cut[j++] = cut[i];
            }
        }
        cut.resize(j);
#ifdef DEBUG_DGL
        Weight dbg_sum = 0;
        for (int i = 0; i < cut.size(); i++) {
            int id = cut[i].id;
            assert(F[id] == g.getWeight(id));
            dbg_sum += F[id];
        }
        assert(dbg_sum == f);
#endif
        return f;
    }

    const Weight getEdgeCapacity(int id){
        assert(g.edgeEnabled(id));
        return g.getWeight(id);
    }

    const Weight getEdgeFlow(int id){
        assert(g.edgeEnabled(id));
        return F[id];                        // reserve(id);
    }

    const Weight getEdgeResidualCapacity(int id){
        assert(g.edgeEnabled(id));
        return g.getWeight(id) - F[id];                        // reserve(id);
    }
};
};
#endif

