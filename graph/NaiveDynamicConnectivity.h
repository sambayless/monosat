/*
 * DynamicConnect.h
 *
 *  Created on: Mar 26, 2014
 *      Author: sam
 */

#ifndef NAIVE_DYNAMICCONNECT_H_
#define NAIVE_DYNAMICCONNECT_H_
#include "DynamicConnectivityImpl.h"
#include "DisjointSets.h"
#include "mtl/Vec.h"
#include <cmath>
using namespace Minisat;
class NaiveDynamicConnectivity:public DynamicConnectivityImpl{



struct Edge{
	int edgeID;
	int from;
	int to;
	bool enabled;
};
DisjointSets sets;
vec<Edge> edges;
int nodes;

private:

bool needsRebuild;

void rebuild(){
	needsRebuild=false;
	sets.Reset();

	sets.AddElements(nodes);

	for(int i =0;i<edges.size();i++){
		if(edges[i].enabled){
			insert(i);
		}
	}
	assert(!needsRebuild);
}

void insert(int edgeID){
	if(!needsRebuild){
		sets.UnionElements(edges[edgeID].from,edges[edgeID].to);
	}
}

void cut(int edgeID){
	if(sets.FindSet(edges[edgeID].from)==sets.FindSet(edges[edgeID].to)){
		needsRebuild=true;
	}else{
		//do nothing
	}
}


public:

NaiveDynamicConnectivity(){
	needsRebuild=false;
	nodes=0;
}

bool connected(int u, int v){
	if(needsRebuild){
		rebuild();
	}
	return sets.FindSet(u)==sets.FindSet(v);
}

int numComponents(){
	return sets.NumSets();
}

void addNode(){
	nodes++;
	edges.push();
	sets.AddElements(1);
}

void addEdge(int edgeID, int from, int to){
	edges.growTo(edgeID+1);
	edges[edgeID].from=from;
	edges[edgeID].to=to;
	edges[edgeID].enabled=false;
}

bool edgeEnabled(int edgeid)const{
	return edges[edgeid].enabled;
}

void setEdgeEnabled(int edgeid, bool enabled){
	if(enabled && ! edges[edgeid].enabled){
		edges[edgeid].enabled=true;
		insert(edgeid);
	}else if(!enabled && edges[edgeid].enabled){
		edges[edgeid].enabled=false;
		cut(edgeid);
	}
}

};

#endif /* DYNAMICCONNECT_H_ */
