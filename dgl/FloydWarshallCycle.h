/****************************************************************************************[Solver.h]
The MIT License (MIT)

Copyright (c) 2014, Sam Bayless

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/

#ifndef FLOYD_WARSHALL_CYCLE_H_
#define FLOYD_WARSHALL_CYCLE_H_

#include "graph/DynamicGraph.h"

namespace dgl{
class FloydWarshallCycle{
	DynamicGraph & g;
	int last_modification;
	int last_addition;
	int last_deletion;
	int addition_qhead;
	int deletion_qhead;
	int lastaddlist;
	int lastdellist;
	std::vector<std::vector<int > > dist;
	std::vector<std::vector<int > > next;
	std::vector<int> dist_changed;
	int INF;
	int cycle_start;

public:
	Cycle():cycle_start(-1), last_modification(-1),last_addition(-1),last_deletion(-1),addition_qhead(0),deletion_qhead(0),lastaddlist(0),lastdellist(0){

	}

private:

	void floyd_warshal(){
		INF = g.nodes()+1;
		dist.resize(g.nodes());
		next.resize(g.nodes());
		for(int i = 0;i<dist.size();i++){
			dist[i].resize(g.nodes());
			for(int j = 0;j<dist[i].size();j++)
				dist[i][j]=INF;
			next[i].resize(g.nodes());
			for(int j = 0;j<next[i].size();j++)
				next[i][j]=-1;

		}

		for(int u = 0;u<g.nodes();u++){
			for(int j = 0;j<g.nIncident(u);j++){
				int v = g.adjacency[j];
				dist[u][v]=1;
			}
		}

		//Apply floyd-warshal
		for(int i = 0;i<g.nodes();i++)
			for(int j = 0;j<g.nodes();j++)
				for(int k = 0;k<g.nodes();k++){
					if (dist[j][i] + dist[i][k] < dist[j][k]){
						dist[j][k] = dist[j][i] + dist[i][k];
						next[j][k]=i;
					}
				}

		//ok, check for cycles
		//As per //http://stackoverflow.com/questions/3911626/find-cycle-of-shortest-length-in-a-directed-graph-with-positive-weights
		//If any vertex has a path to itself < INF, then it is part of a cycle. The minimum cycle is the lowest of these.
		cycle_start = -1;
		int min_cycle_length = INF;
		for(int i = 0;i<g.nodes();i++){
			if(dist[i][i]<min_cycle_length){
				min_cycle_length=dist[i][i];
				cycle_start=i;
			}
		}

	}


	void update( ){

		if(last_deletion==g.deletions && last_modification>0  ){
			//only additions have been made since last time.
			//lets quickly run through those additions; so long as none of them decrease distance, we don't need to do anything more.

			assert(last_deletion==g.deletions);
					last_modification=g.modifications;
					last_addition=g.additions;
					INF=g.nodes()+1;
					dist.resize(g.nodes());

					if(lastaddlist!=g.addlistclears){
						addition_qhead=0;
						lastaddlist=g.addlistclears;
					}
					dist_changed.clear();
					//ok, now check if any of the added edges allow for a decrease in distance.
					for (int i = addition_qhead;i<g.addition_list.size();i++){
						int u=g.addition_list[i].u;
						int v=g.addition_list[i].v;

						if( dist[u][v]>=INF){
							dist_changed.push_back(u);
						}
					}

					//else if(dist_changed.size()> g.nodes()/3){
					if(dist_changed.size())
						floyd_warshal();
					/*}else{
						//try dijkstras for each changed node, instead of running full floyd_warshal.
						for(int i = 0;i<dist_changed.size();i++){
							dijkstras(dist_changed[i]);
						}
					}*/

			return;
		}

		floyd_warshal();

		last_modification=g.modifications;
		last_deletion = g.deletions;
		last_addition=g.additions;

		addition_qhead=g.addition_list.size();
		lastaddlist=g.addlistclears;

		deletion_qhead=g.deletion_list.size();
		lastdellist=g.dellistclears;
	}

	bool hasCycle(){
		if(last_modification!=g.modifications)
				update();
		return cycle_start>-1;
	}
private:
	//From wikipedia; is there a non-recursive version of this?
	void getPath(int i, int j, std::vector<int> & path_out){
		int v = next[i][j];
		if(v>0 && v != cycle_start){
			path_out.push_back(v);
			getPath(i,v,path_out);
			getPath(v,j,path_out);
		}
	}
public:

	//Returns a minimum length cycle, if one exists
	void getCycle(std::vector<int> & cycle_out){
		cycle_out.clear();
		if(hasCycle()){
			cycle_out.push_back(cycle_start);
			getPath(cycle_start,cycle_start,cycle_out);
		}
	}
};
};

#endif /* CYCLE_H_ */
