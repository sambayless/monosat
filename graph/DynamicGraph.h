/*
 * DynamicGraph.h
 *
 *  Created on: 2013-07-15
 *      Author: sam
 */

#ifndef DYNAMICGRAPH_H_
#define DYNAMICGRAPH_H_
#include "mtl/Vec.h"
using namespace Minisat;

class DynamicGraph{
public:
	int nodes;
	int modifications;
	int additions;
	int deletions;
	int addlistclears;
	int dellistclears;
	vec<vec<int> > adjacency;//adj list
	struct EdgeChange{
		int u;//from
		int v;//top
		int mod;
	};
	vec<EdgeChange> addition_list;
	vec<EdgeChange> deletion_list;


	DynamicGraph():nodes(0),modifications(0),additions(0),deletions(0),addlistclears(0),dellistclears(0){}
	void addNodes(int n){
		for(int i = 0;i<n;i++)
			addNode();
	}

	bool hasEdge(int from, int to){
		for(int i = 0;i< adjacency[from].size();i++){
			if(adjacency[from][i]==to){
				return true;
			}
		}
		return false;
	}

	int addNode(){

		adjacency.push();//adj list

	/*	g.push(); // full matrix
		for(int i = 0;i<g.size();g++){
			g[i].growTo(nodes+1);
		}*/
		modifications++;
		additions=modifications;
		deletions=modifications;
		return nodes++;
	}
	void addEdge(int from, int to){
		adjacency[from].push(to);

		modifications++;
		additions=modifications;
		addition_list.push({from,to,modifications});
	}
	//Removes _all_ edges (from, to)
	void removeEdge(int from, int to){
		vec<int>& adj= adjacency[from];
		int i,j = 0;
		for(i = 0;i<adj.size();i++){
			if(adj[i]==to){

			}else{
				adj[j++]=adj[i];
			}
		}
		adj.shrink(i-j);
		modifications++;
		deletions=modifications;
		deletion_list.push({from,to,modifications});
	}

	void clearChangeSets(){
		addition_list.clear();
		deletion_list.clear();
		addlistclears++;
		dellistclears++;
	}
};


#endif /* DYNAMICGRAPH_H_ */
