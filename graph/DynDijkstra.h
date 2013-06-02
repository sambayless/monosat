/*
 * Dijkstra.h
 *
 *  Created on: 2013-05-28
 *      Author: sam
 */

#ifndef DIJKSTRA_H_
#define DIJKSTRA_H_

#include "mtl/Vec.h"
#include "mtl/Heap.h"
using namespace Minisat;




class DynamicGraph{
public:
	int nodes;
	int modifications;
	vec<vec<int> > adjacency;//adj list


	DynamicGraph():nodes(0),modifications(0){}
	void addNodes(int n){
		for(int i = 0;i<n;i++)
			addNode();
	}
	int addNode(){

		adjacency.push();//adj list

	/*	g.push(); // full matrix
		for(int i = 0;i<g.size();g++){
			g[i].growTo(nodes+1);
		}*/
		modifications++;
		return nodes++;
	}
	void addEdge(int from, int to){
		adjacency[from].push(to);

		modifications++;
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
	}
};

class Dijkstra{
public:
	DynamicGraph & g;
	int last_modification;
	int source;
	int INF;
	bool marked;
	vec<int> dist;
	vec<int> prev;
	struct DistCmp{
		vec<int> & _dist;
		 bool operator()(int a, int b)const{
			return _dist[a]<_dist[b];
		}
		 DistCmp(vec<int> & d):_dist(d){};
	};
	Heap<DistCmp> q;
public:

	Dijkstra(int s,DynamicGraph & graph):g(graph), last_modification(-1),source(s),INF(0),marked(false),q(DistCmp(dist)){	}
	Dijkstra(const Dijkstra& d):g(d.g), last_modification(-1),source(d.source),INF(0),marked(false),q(DistCmp(dist)){};


	void setSource(int s){
		source = s;
		last_modification=-1;
	}
	int getSource(){
		return source;
	}
	void update( ){
		last_modification=g.modifications;
		INF=g.nodes+1;
		dist.growTo(g.nodes);
		prev.growTo(g.nodes);
		q.clear();
		for(int i = 0;i<g.nodes;i++){
			dist[i]=i==source?0 :INF;
			prev[i]=-1;
			q.insert(i);
		}
		while(q.size()){
			int u = q.removeMin();
			if(dist[u]==INF)
				break;
			for(int i = 0;i<g.adjacency[u].size();i++){
				int v = g.adjacency[u][i];
				int alt = dist[u]+ 1;
				if(alt<dist[v]){
					dist[v]=alt;
					prev[v]=u;
					q.decrease(v);
				}

			}
		}
	}
	bool connected_unsafe(int t){
		return dist[t]<INF;
	}
	bool connected(int t){
		if(last_modification!=g.modifications)
			update();
		return dist[t]<INF;
	}
	int distance(int t){
		if(last_modification!=g.modifications)
					update();
		return dist[t];
	}
	int previous(int t){
		return prev[t];
	}

};

#endif /* DIJKSTRA_H_ */
