/*
 * DynamicConnect.h
 *
 *  Created on: Mar 26, 2014
 *      Author: sam
 */

#ifndef THORUPDYNAMICCONNECT_H_
#define THORUPDYNAMICCONNECT_H_

#include "DynamicConnectivityImpl.h"
#include "EulerTree.h"
#include <cmath>
#ifndef NDEBUG
#include "NaiveDynamicConnectivity.h"
#endif

class ThorupDynamicConnectivity:public DynamicConnectivityImpl{

#ifndef NDEBUG
	NaiveDynamicConnectivity dbg;
#endif

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
//This is a list of _ALL_ incident edges to each node, except the ones that are actually in the forest, regardless of their level.
vec<vec<int> > incident_edges;
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
		incident_edges[e.from].push(edgeID);
		incident_edges[e.to].push(edgeID);
	}

}

void dbg_checkGraph(){
#ifndef NDEBUG
	for(int i = 0;i<nodes;i++){
		for(int j = 0;j<nodes;j++){
			assert(connected(i,j)==dbg.connected(i,j));
		}
	}

#endif
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

		et[i].cut(edgeID);

		assert(et[i].getFullTreeSize(tu) + et[i].getFullTreeSize(tv) <= 1<<i);//invariant from paper...
		assert(!et[i].connected(tu,tv));//these must be disconnected now that we have cut them

		if(et[i].getFullTreeSize(tv)>et[i].getFullTreeSize(tu)){
			std::swap(tv,tu);
			std::swap(u,v);
		}
		assert(et[i].getFullTreeSize(tv)<=et[i].getFullTreeSize(tu));

		//iterate through all the nodes in this tree in the forest at level i, and see if any have an incident edge at or above this level
		for(auto it = et[i].begin(tv->index);it != et[i].end();++it){
			int v = *it;
			int j,k;
			//Note that because I am only storing one list of incident edges per vertex - rather than one per vertex per level -
			//this loop may do more work than it needs to. In sparsely connected graphs, I am hoping the improved space/locality still makes this a win.
			for(k = 0;k<incident_edges[v].size();k++){
				assert(incident_edges[v][k]!=edgeID);//edges that are in the forest cannot be in the list of incident edges
				Edge & e = edges[incident_edges[v][k]];

				if(e.level==i){
					int otherNode = e.from==v?e.to:e.from;
					if(et[i].connected(tu->index,otherNode)){
						//this is a replacement edge to keep the two components connected
						et[i].link(tv->index,tu->index,e.edgeID);
						//now link all the higher level trees
						for(int h = i+1;j<et.size();j++){
							et[j].cut(edgeID);
							assert(!et[j].connected(u,v));
							et[j].link(u,v,e.edgeID);
						}
						return;
					}else{
						e.level--;
						assert(e.level>=0);
					}
				}
			}
		}

	}

	e.level=levels;

}

void setEdgeLevel(int edgeID, int level){
	assert(level==edges[edgeID].level);
	//et[level].setHasIncidentEdges(edges[edgeID].from,true);
	//et[level].setHasIncidentEdges(edges[edgeID].to,true);
}

public:
ThorupDynamicConnectivity():nodes(0),levels(0){

}
bool connected(int u, int v){
	assert(et.last().connected(u,v)==dbg.connected(u,v));
	return et.last().connected(u,v);
}

int numComponents(){
	assert(et.last().numComponents()==dbg.numComponents());
	return et.last().numComponents();
}

void addNode(){
	nodes++;
	edges.push();
	incident_edges.push();
	levels = (int)(floor(log(nodes)/log(2))+1);
	et.growTo(levels);
	for(EulerTree & t:et){
		while(t.nVertices()<nodes)
			t.createVertex();
		assert(t.nVertices()==nodes);
	}
#ifndef NDEBUG
	dbg.addNode();
#endif
}

void addEdge(int edgeID, int from, int to){
	edges.growTo(edgeID+1);
	edges[edgeID].from=from;
	edges[edgeID].to=to;
	edges[edgeID].level=levels-1;
	setEdgeLevel(edgeID,edges[edgeID].level);
#ifndef NDEBUG
	dbg.addEdge(edgeID,from,to);
#endif
}

bool edgeEnabled(int edgeid)const{
	return edges[edgeid].enabled;
}

void setEdgeEnabled(int edgeID, bool enabled){
	if(enabled && ! edges[edgeID].enabled){
		edges[edgeID].enabled=true;
		insert(edgeID);
	}else if(!enabled && edges[edgeID].enabled){
		edges[edgeID].enabled=false;
		cut(edgeID);
	}
#ifndef NDEBUG
	dbg.setEdgeEnabled(edgeID,enabled);
	dbg_checkGraph();
#endif
}

};

#endif /* DYNAMICCONNECT_H_ */
