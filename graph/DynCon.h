/*
 * DynCon.h
 *
 *  Created on: 2013-05-28
 *      Author: sam
 */

#ifndef DYNCON_H_
#define DYNCON_H_
#include "mtl/Vec.h"
using namespace Minisat;
class dDynamicGraph{
	struct edge{
		int level;
		bool in_tree;
		int l;
		int r;
	};
	LinkCut forest;
	vec<LinkCut::Node*> nodes;
	vec<vec<edge> > adjacency;
	void addNode(){
		nodes.push(forest.addNode());
		adjacency.push();
	}
	void addEdge(int  l, int  r){
		LinkCut::Node * ln = nodes[l];
		LinkCut::Node * rn = nodes[r];

		adjacency[l].push();
		adjacency[l].last().l = l;
		adjacency[l].last().r = r;



		if(!forest.connected(ln,rn)){
			forest.link(ln,rn);
		}
	}
	void removeEdge(int l, int r){
		edge e = adjacency[l][r];

	}
};


#endif /* DYNCON_H_ */
