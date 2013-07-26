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
	int edges;
	int modifications;
	int additions;
	int deletions;
	int historyclears;

	vec<vec<int> > adjacency;//adj list
	//vec<vec<int> > inverted_adjacency;//adj list
	struct EdgeChange{
		bool addition;
		int u;//from
		int v;//top
		int mod;
	};
	vec<EdgeChange> history;

	DynamicGraph():nodes(0),edges(0),modifications(0),additions(0),deletions(0),historyclears(0){}
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
		//inverted_adjacency.push();
		modifications++;
		additions=modifications;
		deletions=modifications;
		clearHistory();
		return nodes++;
	}
	void addEdge(int from, int to){
		assert(from<nodes);
		assert(to<nodes);
		edges++;
		adjacency[from].push(to);
		//inverted_adjacency[to].push(from);
		modifications++;
		additions=modifications;
		history.push({true,from,to,modifications});
	}
	//Removes _all_ edges (from, to)
	void removeEdge(int from, int to){
		{
			vec<int>& adj= adjacency[from];
			int i,j = 0;
			for(i = 0;i<adj.size();i++){
				if(adj[i]==to){
					edges--;
				}else{
					adj[j++]=adj[i];
				}
			}
			adj.shrink(i-j);
		}
	/*	{
			vec<int>& inv_adj= inverted_adjacency[to];
			int i,j = 0;
			for(i = 0;i<inv_adj.size();i++){
				if(inv_adj[i]==from){

				}else{
					inv_adj[j++]=inv_adj[i];
				}
			}
			inv_adj.shrink(i-j);
		}*/
		modifications++;
		deletions=modifications;
		history.push({false,from,to,modifications});
	}

	void clearHistory(){
		if(history.size()>1000){
			history.clear();
			historyclears++;
		}
	}
};


#endif /* DYNAMICGRAPH_H_ */
