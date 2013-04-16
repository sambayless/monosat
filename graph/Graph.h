/*
 * Graph.h
 *
 *  Created on: 2013-04-13
 *      Author: sam
 */

#ifndef GRAPH_H_
#define GRAPH_H_

#include "mtl/Vec.h"
#include "mtl/Heap.h"
#include "mtl/Alg.h"
#include "utils/Options.h"
#include "core/SolverTypes.h"
#include "core/Solver.h"

namespace Minisat {


class Graph {
public:
	const int Undef_Node =-1;
private:
	Solver * S;
	int nodes;

public:



	Graph(Solver * S_):S(S_){

	}
    virtual ~Graph();


public:
	virtual int newNode(){
		return nodes++;
	}
	int nNodes(){
		return nodes;
	}
	bool isNode(int n){
		return n>=0 && n<nodes;
	}
	virtual Lit newEdge(int from,int to)=0;

	virtual void reachesAny(int from, vec<Lit> properties_out,int within_steps=nNodes())=0;
};
};

#endif /* GRAPH_H_ */
