#ifndef MAX_FLOW_H
#define MAX_FLOW_H


#include <iostream>
#include <climits>
#include <vector>

namespace dgl{
struct MaxFlowEdge {
	int u;
	int v;
	int id;
};

template<typename Weight=int>
class MaxFlow{
public:


	virtual ~MaxFlow(){};

    virtual void setCapacity(int u, int w, Weight c)=0;

    virtual void setAllEdgeCapacities(Weight c)=0;
    virtual Weight maxFlow(int s, int t)=0;

    virtual void printStats(){

    }


    virtual Weight minCut(int s, int t, std::vector<MaxFlowEdge> & cut)=0;
    virtual Weight getEdgeFlow(int edgeID)=0;
    virtual Weight getEdgeCapacity(int id)=0;

    virtual  Weight getEdgeResidualCapacity(int id)=0;
};
};
#endif

