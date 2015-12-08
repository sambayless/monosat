/****************************************************************************************[Solver.h]
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
#ifndef MAX_FLOW_H
#define MAX_FLOW_H

#include <iostream>
#include <climits>
#include <vector>

namespace dgl {
struct MaxFlowEdge {
	int u;
	int v;
	int id;
};

template<typename Weight = int>
class MaxFlow {
public:
	
	virtual ~MaxFlow() {
	}
	;

/*	virtual void setCapacity(int u, int w, Weight c)=0;

	virtual void setAllEdgeCapacities(Weight c)=0;*/
	const Weight maxFlow() {
		return update();
	}
	virtual const Weight maxFlow(int from, int to) {
		setSource(from);
		setSink(to);
		return maxFlow();
	}
	virtual const Weight update()=0;
	virtual int numUpdates() const=0;
	virtual void setSource(int s)=0;
	virtual void setSink(int t)=0;

	virtual std::vector<int> & getChangedEdges()=0;
	virtual void clearChangedEdges() {
		
	}
	virtual std::vector<int> & getChangedPartition(){
		throw std::runtime_error("Not implemented");
	}
	virtual void clearChangedPartition() {

	}

	virtual void printStats() {
		
	}
	virtual int getSource() const=0;
	virtual int getSink() const=0;
	virtual const Weight minCut(std::vector<MaxFlowEdge> & cut)=0;
	//Return a valid cut that may or may not be the mincut.
	virtual const Weight approxMinCut(std::vector<MaxFlowEdge> & cut) {
		return minCut(cut);
	}
	
	virtual const Weight getEdgeFlow(int edgeID)=0;
	virtual const Weight getEdgeCapacity(int id)=0;

	virtual const Weight getEdgeResidualCapacity(int id)=0;

	virtual const bool isOnCut(int edgeID){
		throw std::runtime_error("Not implemented");
	}
	virtual const bool inSourcePartition(int node){
		throw std::runtime_error("Not implemented");
	}
};
}
;
#endif

