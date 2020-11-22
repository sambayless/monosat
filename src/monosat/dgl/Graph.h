/**************************************************************************************************
 The MIT License (MIT)

 Copyright (c) 2018, Sam Bayless

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

#ifndef GRAPH_H_
#define GRAPH_H_

#include <vector>
#include <algorithm>
#include <string>
#include <cassert>
#include <sstream>
#include <cstdint>
#include <sstream>
#include <cstdio>

namespace dgl {
/**
 * Optional callback interface for dynamic graph algorithms.
 */
class DynamicGraphAlgorithm {
public:
    virtual std::string getName() = 0;

    virtual ~DynamicGraphAlgorithm(){}

    virtual void updateHistory() = 0;
};

/**
 * Abstract interface for dynamic graphs.
 * See DynamicGraph.h for a concrete implementation.
 */
template<typename Weight>
class Graph {
public:
    struct Edge {
        int node;
        int id;
    };

    struct FullEdge {
        int from;
        int to;
        int id;

        FullEdge() :
                from(-1), to(-1), id(-1){
        }

        FullEdge(int from, int to, int id) :
                from(from), to(to), id(id){
        }
    };

    struct EdgeChange {
        bool addition;
        bool deletion;
        bool weight_increase;
        bool weight_decrease;
        int id;
        int mod;
        int prev_mod;
    };


    Graph(){
    }

    virtual ~Graph(){

    }

    virtual FILE* outfile() = 0;

    virtual void addNodes(int n) = 0;

    //Returns true iff the edge exists and is a self loop
    virtual bool selfLoop(int edgeID) = 0;

    //SLOW!
    virtual bool hasEdge(int from, int to) const = 0;

    //SLOW! Returns -1 if there is no edge
    virtual int getEdge(int from, int to) const = 0;

    virtual bool hasEdgeUndirected(int from, int to) const = 0;

    virtual int addNode() = 0;

    //true iff the edge's current assignment will never be altered.
    virtual bool isConstant(int edgeID) const = 0;

    virtual void makeEdgeAssignmentConstant(int edgeID) = 0;

    virtual bool edgeEnabled(int edgeID) const = 0;

    virtual bool isEdge(int edgeID) const = 0;

    virtual bool hasEdge(int edgeID) const = 0;

    //Instead of actually adding and removing edges, tag each edge with an 'enabled/disabled' label, and just expect reading algorithms to check and respect that label.
    virtual int addEdge(int from, int to, int id = -1, Weight weight = 1) = 0;

    virtual int nEdgeIDs() = 0;

    virtual int nodes() const = 0;

    virtual int edges() const = 0;

    virtual int nIncident(int node, bool undirected = false) = 0;

    virtual int nDirectedEdges(int node, bool incoming) = 0;

    virtual Edge& directedEdges(int node, int i, bool is_incoming) = 0;

    virtual int nIncoming(int node, bool undirected = false) = 0;

    virtual Edge& incident(int node, int i, bool undirected = false) = 0;

    virtual Edge& incoming(int node, int i, bool undirected = false) = 0;

    virtual std::vector<FullEdge>& getEdges() = 0;

    virtual std::vector<Weight>& getWeights() = 0;

    virtual Weight getWeight(int edgeID) = 0;

    virtual FullEdge& getEdge(int id) = 0;

    virtual void setEdgeEnabled(int id, bool enable) = 0;

    virtual void enableEdge(int id) = 0;

    virtual void disableEdge(int id) = 0;

    virtual void enableEdge(int from, int to, int id) = 0;

    virtual bool undoEnableEdge(int id) = 0;

    virtual void disableEdge(int from, int to, int id) = 0;

    virtual bool undoDisableEdge(int id) = 0;

    virtual Weight getEdgeWeight(int edgeID) = 0;

    virtual void setEdgeWeight(int id, const Weight& w) = 0;

    virtual void drawFull(bool showWeights = false, bool force_draw = false) = 0;

    virtual bool rewindHistory(int steps) = 0;

    /**
     * Returns a unique identifier for this algorithm.
     */
    virtual int addDynamicAlgorithm(DynamicGraphAlgorithm* alg) = 0;

    virtual void updateAlgorithmHistory(DynamicGraphAlgorithm* alg, int algorithmID, int historyPos) = 0;

    virtual EdgeChange& getChange(int64_t historyPos) = 0;

    virtual int historySize() = 0;

    virtual int nHistoryClears() const = 0;

    virtual int getCurrentHistory() const = 0;

    virtual int nDeletions() const = 0;

    virtual int nAdditions() const = 0;

    virtual int lastEdgeIncrease() const = 0;

    virtual int lastEdgeDecrease() const = 0;

    virtual void clearHistory(bool forceClear = false) = 0;

    //force a new modification
    virtual void invalidate() = 0;

    virtual int64_t getPreviousHistorySize() const = 0;

    virtual void markChanged() = 0;

    virtual bool changed() = 0;

    virtual void clearChanged() = 0;

    virtual void clear() = 0;


};

};
#endif /* DYNAMICGRAPH_H_ */
