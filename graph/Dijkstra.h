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
#include "DynamicGraph.h"
using namespace Minisat;

/*
class GraphListener{
	void addEdge(int u, int v, int mod);
	void removeEdge(int u, int v, int mod);
};*/


class Dijkstra{
public:
	DynamicGraph & g;
	int last_modification;
	int last_addition;
	int last_deletion;
	int addition_qhead;
	int deletion_qhead;
	int lastaddlist;
	int lastdellist;
	int source;
	int INF;

	bool marked;

	vec<int> old_dist;
	vec<int> changed;

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

	Dijkstra(int s,DynamicGraph & graph):g(graph), last_modification(-1),last_addition(-1),last_deletion(-1),addition_qhead(0),deletion_qhead(0),lastaddlist(0),lastdellist(0),source(s),INF(0),marked(false),q(DistCmp(dist)){	}
	Dijkstra(const Dijkstra& d):g(d.g), last_modification(-1),last_addition(-1),last_deletion(-1),addition_qhead(0),deletion_qhead(0),lastaddlist(0),lastdellist(0),source(d.source),INF(0),marked(false),q(DistCmp(dist)){};


	void setSource(int s){
		source = s;
		last_modification=-1;
		last_addition=-1;
		last_deletion=-1;
	}
	int getSource(){
		return source;
	}

	void updateFast(){
		/*for(int i = 0;i<g.nodes;i++)
					changed.push(i);*/
		assert(last_deletion==g.deletions);
		last_modification=g.modifications;
		last_addition=g.additions;
		INF=g.nodes+1;
		dist.growTo(g.nodes);
		prev.growTo(g.nodes);
		q.clear();
		if(lastaddlist!=g.addlistclears){
			addition_qhead=0;
			lastaddlist=g.addlistclears;
		}
		//ok, now check if any of the added edges allow for a decrease in distance.
		for (int i = addition_qhead;i<g.addition_list.size();i++){
			int u=g.addition_list[i].u;
			int v=g.addition_list[i].v;
			int alt = dist[u]+1 ;
			if(alt< dist[v]){

				if(dist[v]>=INF){
					//this was changed
					changed.push(v);
				}

				dist[v]=alt;
				prev[v]=u;

				if(!q.inHeap(v))
					q.insert(v);
			}

			/*
			 *
			 *
			//Is this altered code still correct?
			if(dist[v]>=INF){
				//this was changed
				changed.push(v);

				dist[v]=alt;
				prev[v]=u;

				if(!q.inHeap(v))
					q.insert(v);
			}
			 *
			 */

		}
		addition_qhead=g.addition_list.size();

		while(q.size()){
			int u = q.removeMin();
			if(dist[u]==INF)
				break;
			for(int i = 0;i<g.adjacency[u].size();i++){
				int v = g.adjacency[u][i];
				int alt = dist[u]+ 1;
				if(alt<dist[v]){
					if(dist[v]>=INF){
						//this was changed
						changed.push(v);
					}
					dist[v]=alt;
					prev[v]=u;
					if(!q.inHeap(v))
						q.insert(v);
					else
						q.decrease(v);
				}

			}
		}
	}
	vec<int> & getChanged(){
		return changed;
	}
	void clearChanged(){
		changed.clear();
	}
	void update( ){
		//changed.clear();
		//for(int i = 0;i<g.nodes;i++)
		//	changed.push(i);
		/*if (last_addition==g.additions && last_modification>0){
			if(lastdellist!=g.dellistclears){
				deletion_qhead=0;
				lastdellist=g.dellistclears;
			}
			//ok, now check if any of the added edges allow for a decrease in distance.
			for (int i = deletion_qhead;i<g.deletion_list.size();i++){
				int u=g.deletion_list[i].u;
				int v=g.deletion_list[i].v;
				if(prev[v]==u){
					deletion_qhead = i-1;
					//this deletion matters, so we need to recompute.
					break;
				}
			}
			//none of these deletions touched any shortest paths, so we can ignore them.
			deletion_qhead=g.deletion_list.size();
			last_deletion = g.deletions;
		}*/

		if(last_deletion==g.deletions && last_modification>0  ){
			//we can use a faster, simple dijkstra update method
			updateFast();
			return;
		}
		//if none of the deletions were to 'previous' edges, then we don't need to do anything

		//for(int i = 0;i<g.nodes;i++)
		//	changed.push(i);

		INF=g.nodes+1;
		dist.growTo(g.nodes);
		prev.growTo(g.nodes);
		old_dist.growTo(g.nodes);
		q.clear();
		for(int i = 0;i<g.nodes;i++){
			old_dist[i]=last_modification > 0 ? dist[i]:INF;//this won't work properly if we added nodes...
			dist[i]=i==source?0 :INF;
			prev[i]=-1;
			q.insert(i);
		}
		while(q.size()){
			int u = q.peakMin();
			if(dist[u]==INF)
				break;
			if(old_dist[u]>=INF){
				changed.push(u);
			}
			q.removeMin();
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
		while(q.size()){
			//iterate through the unreached nodes and check which ones were previously reached
			int u = q.removeMin();
			if(last_modification <=0  || (old_dist[u]<INF && dist[u]>=INF)){
				changed.push(u);
			}
		}
	

	last_modification=g.modifications;
		last_deletion = g.deletions;
		last_addition=g.additions;

		addition_qhead=g.addition_list.size();
		lastaddlist=g.addlistclears;

		deletion_qhead=g.deletion_list.size();
		lastdellist=g.dellistclears;


	}

	bool dbg_uptodate(){
	/*	DynamicGraph gdbg;
		for(int i = 0;i<g.nodes;i++){
			gdbg.addNode();
		}

		for(int i = 0;i<g.adjacency.size();i++){
			for(int j = 0;j<g.adjacency[i].size();j++){
				int u = g.adjacency[i][j];
				gdbg.addEdge(i,u);
			}
		}*/

		Dijkstra d(source,g);
		d.update();
		for(int i = 0;i<g.nodes;i++){
			int distance = dist[i];
			int dbgdist = d.dist[i];
			assert(distance==dbgdist);
		}
		return true;
	}

	bool connected_unsafe(int t){
		return t<dist.size() && dist[t]<INF;
	}
	bool connected(int t){
		if(last_modification!=g.modifications)
			update();

		assert(dbg_uptodate());

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
