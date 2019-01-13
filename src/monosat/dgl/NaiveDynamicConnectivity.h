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
#ifndef NAIVE_DYNAMICCONNECT_H_
#define NAIVE_DYNAMICCONNECT_H_

#include "DynamicConnectivityImpl.h"
#include "monosat/dgl/alg/DisjointSets.h"
#include <vector>
#include "monosat/mtl/Sort.h"
#include <cmath>

namespace dgl {
class NaiveDynamicConnectivity : public DynamicConnectivityImpl {

    struct Edge {
        int edgeID;
        int from;
        int to;
        bool enabled;
    };
    DisjointSets sets;
    std::vector<Edge> edges;
    int nodes;

private:

    bool needsRebuild;

    void rebuild(){
        needsRebuild = false;
        sets.Reset();

        sets.AddElements(nodes);

        for(int i = 0; i < edges.size(); i++){
            if(edges[i].enabled){
                insert(i);
            }
        }
        assert(!needsRebuild);
    }

    void insert(int edgeID){
        if(!needsRebuild){
            sets.UnionElements(edges[edgeID].from, edges[edgeID].to);
        }
    }

    void cut(int edgeID){
        if(sets.FindSet(edges[edgeID].from) == sets.FindSet(edges[edgeID].to)){
            needsRebuild = true;
        }else{
            //do nothing
        }
    }

public:

    NaiveDynamicConnectivity(){
        needsRebuild = false;
        nodes = 0;
    }

    bool connected(int u, int v) override{
        if(needsRebuild){
            rebuild();
        }
        return sets.FindSet(u) == sets.FindSet(v);
    }

    int numComponents() override{
        return sets.NumSets();
    }

    void dbg_print() override{
#ifdef DEBUG_DGL
        std::vector<bool> seen;
        seen.resize(nodes);
        for (int n = 0; n < nodes; n++) {
            if (!seen[n]) {
                seen[n] = true;

                std::vector<int> node_set;
                node_set.push_back(n);
                for (int j = 0; j < nodes; j++) {
                    if (j != n) {
                        if (sets.FindSet(j) == sets.FindSet(n)) {
                            node_set.push_back(j);
                        }
                    }
                }
                std::sort(node_set.begin(), node_set.end());
                if (node_set.size() > 1) {
                    for (int n : node_set) {
                        seen[n] = true;
                        printf("%d,", n);
                    }
                    printf("\n");
                }
                /*	for(int i =0;i<e.size();i++){
                 printf("(%d->%d)", edges[e[i]].from,edges[e[i]].to );
                 }
                 printf("\n");*/
            }
        }
#endif
    }

    void addNode() override{
        nodes++;
        edges.push_back({});
        sets.AddElements(1);
    }

    void addEdge(int from, int to, int edgeID) override{
        if(edges.size() <= edgeID)
            edges.resize(edgeID + 1);
        edges[edgeID].from = from;
        edges[edgeID].to = to;
        edges[edgeID].enabled = false;
    }

    bool edgeEnabled(int edgeid) const override{
        return edges[edgeid].enabled;
    }

    bool setEdgeEnabled(int from, int to, int edgeid, bool enabled) override{
        if(enabled && !edges[edgeid].enabled){
            edges[edgeid].enabled = true;
            insert(edgeid);
            return true;
        }else if(!enabled && edges[edgeid].enabled){
            edges[edgeid].enabled = false;
            cut(edgeid);
            return true;
        }
        return false;
    }

};
};
#endif /* DYNAMICCONNECT_H_ */
