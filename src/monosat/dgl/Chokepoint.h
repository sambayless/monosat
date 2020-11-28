/*
 * Chokepoint.h
 *
 *  Created on: 2013-08-04
 *      Author: sam
 */
#ifndef CHOKEPOINT_H_
#define CHOKEPOINT_H_
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

#include "Graph.h"
#include "DynamicGraph.h"
#include "Reach.h"

namespace dgl {

struct ForceReason {
    int edge_id;
    int node;
};

template<typename Weight, typename Graph = DynamicGraph<Weight>, class EdgeStatus= typename Reach::NullStatus>
class Chokepoint {

    Graph& g;
    EdgeStatus& status;
    int source;
    std::vector<int> current;
    std::vector<int> prev;
    std::vector<int> queue;
    std::vector<int> dist;
    std::vector<bool> in_queue;

    static const int SOURCE = -3;
    static const int UNDEF = -2;
    static const int EMPTY = -1;
public:
    Chokepoint(EdgeStatus& _status, Graph& _graph, int _source) :
            g(_graph), status(_status), source(_source){

    }

    void update(){

        in_queue.clear();
        queue.clear();

        current.resize(g.nodes());
        prev.resize(g.nodes());
        dist.resize(g.nodes());
        in_queue.resize(g.nodes());

        for(int i = 0; i < g.nodes(); i++){
            prev[i] = UNDEF;
            current[i] = UNDEF;
            dist[i] = UNDEF;
        }

        prev[source] = SOURCE;
        current[source] = SOURCE;
        dist[source] = 0;

        in_queue[source] = 1;
        queue.push_back(source);

        for(int i = 0; i < queue.size(); i++){
            int u = queue[i];

            int old_dist = dist[u];

            assert(in_queue[u]);

            if(u == source){
                //a better way to model this would be to give the source node an extra incoming edge
            }else{
                in_queue[u] = false;
                current[u] = UNDEF;
                prev[u] = UNDEF;

                for(int j = 0; j < g.nIncoming(u); j++){
                    int from = g.incoming(u, j).node;
                    int id = g.incoming(u, j).id;
                    if(!g.edgeEnabled(id))
                        continue;

                    //From does not (yet) have a path from source (if it later gets one, then u will be re-added to the queue and re-processed
                    if(prev[from] == UNDEF)
                        continue;

                    if(prev[u] == UNDEF){
                        prev[u] = from;
                        if(status(id)){
                            //then this edge is a candidate edge
                            current[u] = id;
                        }else{
                            current[u] = UNDEF;
                        }
                        dist[u] = dist[from] + 1;

                        continue;
                    }

                    //compare the incoming list to the previous list.
                    //First, make the distances equal
                    while(dist[u] > dist[from] + 1){
                        prev[u] = prev[prev[u]];
                        dist[u]--;
                        assert(dist[i] == dist[prev[u]] + 1);
                    }

                    while(dist[u] >= 0 && dist[from] > dist[u]){
                        from = prev[from];
                    }

                    //now, backtrack in both linked lists until they have the same current element
                    while(current[u] != current[from]){
                        prev[u] = prev[prev[u]];
                        dist[u]--;
                        assert(dist[u] >= 0);
                        from = prev[from];
                        assert(dist[i] == dist[prev[u]] + 1);
                    }

                    assert(dist[u] >= 0);

                }
            }
            //ok, if we changed anything, update all outgoing.edges()

            if(dist[u] != old_dist || u == source){
                assert(dist[u] >= 0);

                for(int j = 0; j < g.nIncident(u); j++){
                    int to = g.incident(u, j).node;
                    int id = g.incident(u, j).id;
                    if(!g.edgeEnabled(id))
                        continue;

                    if(!in_queue[to]){
                        queue.push_back(to);
                        in_queue[to] = true;
                    }

                }

            }

        }

    }

    void collectForcedEdges(std::vector<ForceReason>& forced_ids){
        update();
        for(int i = 0; i < g.nodes(); i++){
            if(status.mustReach(i)){
                int u = i;
                while(prev[u] >= 0){
                    if(current[u] != UNDEF){
                        forced_ids.push_back({current[u], i});
                        current[u] = UNDEF;
                    }
                    int v = prev[u];
                    prev[u] = UNDEF;
                    u = v;
                }
            }
        }
    }

};

};
#endif

