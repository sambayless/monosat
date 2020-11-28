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

#ifndef DYNAMIC_BACK_GRAPH_H_
#define DYNAMIC_BACK_GRAPH_H_

#include <vector>
#include <algorithm>
#include <cassert>
#include <sstream>
#include <cstdint>
#include <sstream>
#include <cstdio>
#include "Graph.h"
#include "DynamicGraph.h"

namespace dgl {


/**
 * A thin wrapper around a dynamic graph, which reverses all the edges in the original graph.
 */
template<typename Weight>
class DynamicBackGraph final : public Graph<Weight> {
public:
    typedef typename Graph<Weight>::Edge Edge;
    typedef typename Graph<Weight>::FullEdge FullEdge;
    typedef typename Graph<Weight>::EdgeChange EdgeChange;
private:
    DynamicGraph<Weight>& base;
    std::vector<FullEdge> all_back_edges;

    void updateEdges(){
        for(int id = all_back_edges.size(); id < base.nEdgeIDs(); id++){
            auto& forward_edge = base.getEdge(id);
            all_back_edges.push_back({forward_edge.to, forward_edge.from, forward_edge.id});
            assert(all_back_edges[id].id == id);
        }
        assert(nEdgeIDs() == base.nEdgeIDs());
        assert(all_back_edges.size() == base.nEdgeIDs());
    }

public:

    DynamicBackGraph(DynamicGraph<Weight>& base) : base(base){

    }

    ~DynamicBackGraph(){

    }

    FILE* outfile() override{
        return base.outfile();
    };

    void addNodes(int n) override{base.addNodes(n);};

    //Returns true iff the edge exists and is a self loop
    bool selfLoop(int edgeID) override{return base.selfLoop(edgeID);};

    //SLOW!
    bool hasEdge(int from, int to) const override{return base.hasEdge(to, from);};

    //SLOW! Returns -1 if there is no edge
    int getEdge(int from, int to) const override{return base.getEdge(to, from);};

    bool hasEdgeUndirected(int from, int to) const override{return base.hasEdgeUndirected(to, from);};

    int addNode() override{return base.addNode();};

    //true iff the edge's current assignment will never be altered.
    bool isConstant(int edgeID) const override{return base.isConstant(edgeID);};

    void makeEdgeAssignmentConstant(int edgeID) override{base.makeEdgeAssignmentConstant(edgeID);};

    bool edgeEnabled(int edgeID) const override{return base.edgeEnabled(edgeID);};

    bool isEdge(int edgeID) const override{return base.isEdge(edgeID);};

    bool hasEdge(int edgeID) const override{return base.hasEdge(edgeID);};

    int addEdge(int from, int to, int id, Weight weight = 1) override{
        updateEdges();
        int edgeID = base.addEdge(to, from, id, weight);
        assert(all_back_edges.size() == base.nEdgeIDs() - 1);
        if(all_back_edges.size() <= id)
            all_back_edges.resize(id + 1);
        all_back_edges[id] = {from, to, id}; //,weight};
        return edgeID;
    };

    int nEdgeIDs() override{return base.nEdgeIDs();};

    int nodes() const override{return base.nodes();};

    int edges() const override{return base.edges();};

    int nIncident(int node, bool undirected = false) override{return base.nIncoming(node, undirected);};

    int nDirectedEdges(int node, bool incoming) override{return base.nDirectedEdges(node, !incoming);};


    Edge& directedEdges(int node, int i, bool is_incoming) override{
        return base.directedEdges(node, i, !is_incoming);
    };

    int nIncoming(int node, bool undirected = false) override{
        return base.nIncident(node, undirected);
    };

    Edge& incident(int node, int i, bool undirected = false) override{
        return base.incoming(node, i, undirected);
    };

    Edge& incoming(int node, int i, bool undirected = false) override{
        return base.incident(node, i, undirected);
    };

    std::vector<FullEdge>& getEdges() override{
        updateEdges();
        return all_back_edges;
    };

    std::vector<Weight>& getWeights() override{return base.getWeights();};

    Weight getWeight(int edgeID) override{return base.getWeight(edgeID);};

    FullEdge& getEdge(int id) override{
        updateEdges();
        assert(all_back_edges.size() == base.nEdgeIDs());
        assert(all_back_edges.size() > id);
        assert(all_back_edges[id].id == id);
        return all_back_edges[id];
    };

    void setEdgeEnabled(int id, bool enable) override{base.setEdgeEnabled(id, enable);};

    void enableEdge(int id) override{base.enableEdge(id);};

    void disableEdge(int id) override{base.disableEdge(id);};

    void enableEdge(int from, int to, int id) override{base.enableEdge(to, from, id);};

    bool undoEnableEdge(int id) override{return base.undoEnableEdge(id);};

    void disableEdge(int from, int to, int id) override{return base.disableEdge(to, from, id);};

    bool undoDisableEdge(int id) override{return base.undoDisableEdge(id);};

    Weight getEdgeWeight(int edgeID) override{return base.getEdgeWeight(edgeID);};

    void setEdgeWeight(int id, const Weight& w) override{base.setEdgeWeight(id, w);};

    void drawFull(bool showWeights = false, bool force_draw = false) override{base.drawFull(showWeights, force_draw);};

    bool rewindHistory(int steps) override{return base.rewindHistory(steps);};

    /**
     * Returns a unique identifier for this algorithm.
     */
    int addDynamicAlgorithm(DynamicGraphAlgorithm* alg) override{
        if(alg == nullptr){
            throw std::runtime_error("Internal error in DynamicBackGraph (null dynamic graph algorithm)");
        }
        return base.addDynamicAlgorithm(alg);
    };

    void updateAlgorithmHistory(DynamicGraphAlgorithm* alg, int algorithmID, int historyPos) override{
        if(algorithmID < 0 || alg == nullptr){
            if(alg == nullptr){
                throw std::runtime_error("Internal error in DynamicBackGraph (invalid algorithm ID): algorithm " +
                                         std::to_string(algorithmID) + ", null algorithm");
            }else{
                throw std::runtime_error("Internal error in DynamicBackGraph (invalid algorithm ID): algorithm " +
                                         std::to_string(algorithmID) + ", " + alg->getName());
            }
        }
        base.updateAlgorithmHistory(alg, algorithmID, historyPos);
    };

    EdgeChange& getChange(int64_t historyPos) override{
        return base.getChange(historyPos);
    };

    int historySize() override{return base.historySize();};

    int nHistoryClears() const override{
        return base.nHistoryClears();
    }

    int getCurrentHistory() const override{return base.getCurrentHistory();};

    int nDeletions() const override{
        return base.nDeletions();
    };

    int nAdditions() const override{
        return base.nAdditions();
    }

    int lastEdgeIncrease() const override{
        return base.lastEdgeIncrease();
    }

    int lastEdgeDecrease() const override{
        return base.lastEdgeDecrease();
    }

    void clearHistory(bool forceClear = false) override{base.clearHistory(forceClear);};

    void invalidate() override{base.invalidate();};

    int64_t getPreviousHistorySize() const override{return base.getPreviousHistorySize();};

    void markChanged() override{base.markChanged();};

    bool changed() override{return base.changed();};

    void clearChanged() override{base.clearChanged();};

    void clear() override{base.clear();};
};

};
#endif /* DYNAMIC_BACK_GRAPH_H_ */
