/*
 * DynamicConnect.h
 *
 *  Created on: Mar 26, 2014
 *      Author: sam
 */

#ifndef THORUPDYNAMICCONNECT_H_
#define THORUPDYNAMICCONNECT_H_

#include "EulerTree.h"
#include <cmath>
class ThorupDynamicConnectivity{



struct Edge{
	int edgeID;
	int from;
	int to;
	bool in_forest;
	bool enabled;
	int level;

};

vec<Edge> edges;
int nodes;

vec<EulerTree> et;

private:

int levels;

void insert(int edgeID){
	Edge & e = edges[edgeID];
	e.level=levels-1;
	if (!et.last().connected(e.from,e.to)){
		et.last().link(e.from,e.to,edgeID);
		e.in_forest=true;
	}else{
		e.in_forest=false;
	}
}

void cut(int edgeID){
	Edge & e = edges[edgeID];
	if(!e.in_forest){
		return;//this edge is not in the top level forest (and hence also not in _any_ level forest), so we don't need to do anything special to remove it
	}
	int u = e.from;
	int v = e.to;
	//find a replacement edge to connect u and v, if one exists
	for(int i = e.level;i<levels;i++){

		EulerTree::EulerVertex * tu = et[i].getVertex(u);
		EulerTree::EulerVertex * tv = et[i].getVertex(v);

		/*et[i].cutEdge(edgeID);

		assert(tu->getSize() + tv->getSize()<= 1<<i);//invariant from paper...
		assert(!et[i].connected(tu,tv));//these must be disconnected now that we have cut them

		if(tv->getSize()>tu->getSize()){
			std::swap(tv,tu);
		}
		assert(tv->getSize()<=tu->getSize());

		//ok, now explore tv and see if we can find a replacement edge at level i to connect tv to tu.
		if(tv->has_incident_edges){
			//need to traverse edges in graphI (but not in forest I) here.
			//If no incident edges actually exit, then clear the flag
		}
		if(tv->subtree_has_incident_edges){
			//ok, explore the subtree.
			EulerTree::EulerVertex * n = et[i].getNext(tv);
			while(n!=tv){
				if(n->has_incident_edges){
					//explore
				}
				n= et[i].getNext(tv);
			}
		}*/
		//also need to clear the subtree_has_incident_edges flag, if no incident edges existed...


	}

	e.level=levels;
}

void setEdgeLevel(int edgeID, int level){
	assert(level==edges[edgeID].level);
	//et[level].setHasIncidentEdges(edges[edgeID].from,true);
	//et[level].setHasIncidentEdges(edges[edgeID].to,true);
}

public:

bool connected(int u, int v){
	return et.last().connected(u,v);
}

int numComponents(){
	return et.last().numComponents();
}

void addNode(){
	nodes++;

	levels = (int)(floor(log(nodes)/log(2))+1);
}

void addEdge(int edgeID, int from, int to){
	edges.growTo(edgeID+1);
	edges.last().from=from;
	edges.last().to=to;
	edges.last().level=levels-1;
	setEdgeLevel(edgeID,edges.last().level);
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
