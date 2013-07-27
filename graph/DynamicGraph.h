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
	int next_id;

	struct Edge{
		int to;
		int id;
	};
	vec<vec<Edge> > adjacency;//adj list
	//vec<vec<int> > inverted_adjacency;//adj list
	struct EdgeChange{
		bool addition;
		int u;//from
		int v;//top
		int id;
		int mod;
	};
	vec<EdgeChange> history;
	vec<char> edge_status;

	DynamicGraph():nodes(0),edges(0),modifications(0),additions(0),deletions(0),historyclears(0),next_id(0){}
	void addNodes(int n){
		for(int i = 0;i<n;i++)
			addNode();
	}

	bool hasEdge(int from, int to){
		for(int i = 0;i< adjacency[from].size();i++){
			if(adjacency[from][i].to==to && edgeEnabled(adjacency[from][i].id)){
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

	bool edgeEnabled(int edgeID){
		assert(edgeID<edge_status.size());
		return edge_status[edgeID];
	}

	//Instead of actually adding and removing edges, tag each edge with an 'enabled/disabled' label, and just expect reading algorithms to check and respect that label.
	void addEdge(int from, int to, int id=-1){
		assert(from<nodes);
		assert(to<nodes);
		if(id<0){
			id = next_id++;
		}else{
			if(id>=next_id){
				next_id=id+1;
			}
		}
		edges++;
		adjacency[from].push({to,id});
		edge_status.growTo(id+1);
		//inverted_adjacency[to].push(from);
		modifications++;
		additions=modifications;
		history.push({true,from,to,id,modifications});
		enableEdge(from,to,id);//default to enabled
	}

	void enableEdge(int from, int to, int id){
		assert(id>=0);
		assert(id<edge_status.size());
		if(edge_status[id]!=true){
			edge_status[id]=true;
			modifications++;
			additions=modifications;
			history.push({true,from,to,id,modifications});
		}
	}

	void disableEdge(int from, int to, int id){
		assert(id>=0);
		assert(id<edge_status.size());
		if(edge_status[id]!=false){
			edge_status[id]=false;
			modifications++;
			additions=modifications;
			history.push({false,from,to,id,modifications});
		}
	}

	//Removes _all_ edges (from, to)
	/*void removeEdge(int from, int to, int id){
		assert(id>=0);
		assert(id<edge_status.size());
		{
			vec<Edge>& adj= adjacency[from];
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

		modifications++;
		deletions=modifications;
		history.push({false,from,to,id,modifications});
	}*/

	void clearHistory(){
		if(history.size()>1000){
			history.clear();
			historyclears++;
		}
	}
};


#endif /* DYNAMICGRAPH_H_ */
