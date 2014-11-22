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

#ifndef TARJANSSCC_H_
#define TARJANSSCC_H_


#include <vector>
#include "alg/Heap.h"
#include "DynamicGraph.h"
#include "core/Config.h"
#include "ConnectedComponents.h"
#include "alg/DisjointSets.h"
#include <limits>
#include <algorithm>

namespace dgl{
template<class Status=ConnectedComponents::NullConnectedComponentsStatus>
class TarjansSCC:public ConnectedComponents{
public:

	DynamicGraph & g;
	Status &  status;
	int last_modification;

	int last_addition;
	int last_deletion;
	int history_qhead;

	int last_history_clear;
	int INF;

	struct Component{
		int id=-1;//the unique id of the component
		int next=-1;//pointer to the next node in the component
		Component(){}
		Component(int id, int next):id(id),next(next){}
	};

	std::vector<Component> scc; //the components of each node
	struct SCC{
		int sz;//size of the component
		int element;//pointer to an arbitrary node in the scc.
	};
	std::vector<SCC> scc_set;

	std::vector<int> q;

	std::vector<int> indices;
	std::vector<int> lowlink;
	std::vector<bool> in_q;


	int stats_full_updates=0;
	int stats_fast_updates=0;
	int stats_fast_failed_updates=0;
	int stats_skip_deletes=0;
	int stats_skipped_updates=0;
	int stats_num_skipable_deletions=0;
	double mod_percentage=0.2;

	double stats_full_update_time=0;
	double stats_fast_update_time=0;

public:
	TarjansSCC(DynamicGraph & graph):g(graph), status(nullConnectedComponentsStatus), last_modification(-1),last_addition(-1),last_deletion(-1),history_qhead(0),last_history_clear(0),INF(0){

	}

	TarjansSCC(DynamicGraph & graph, Status & _status ):g(graph), status(_status), last_modification(-1),last_addition(-1),last_deletion(-1),history_qhead(0),last_history_clear(0),INF(0){

	}

	void setNodes(int n){
		q.reserve(n);
		in_q.resize(n);
		indices.resize(n,-1);
		lowlink.resize(n,-1);
		scc.resize(n);
		INF=std::numeric_limits<int>::max();
	}


	void strongConnect(int node, int & index, std::vector<int> * scc_out =nullptr){
		//Following wikipedia's pseudocode, from http://en.wikipedia.org/wiki/Tarjan%27s_strongly_connected_components_algorithm
		assert(indices[node]<0);
		assert(lowlink[node]<0);
		indices[node]=index;
		lowlink[node]=index;
		index++;
		assert(!in_q[node]);
		q.push_back(node);
		in_q[node]=true;

		//for(auto & edge:g.adjacency[node]){
		for(int i = 0;i<g.nIncident(node);i++){
			auto & edge = g.incident(node,i);
			if(g.edgeEnabled(edge.id)){
				int to = edge.node;
				assert(in_q[node]==(std::count( q.begin(),q.end(),node)));
				if(indices[to]<0){
					strongConnect(to,index,scc_out);
					lowlink[node]=std::min(lowlink[node],lowlink[to]);
				}else if(in_q[to]){
					lowlink[node]= std::min(lowlink[node],indices[to]);
				}
			}
		}

		// If v is a root node, pop the stack and generate an SCC
		if(lowlink[node]==indices[node]){
			int sccID = scc_set.size();
			int sz = q.size();
			assert(q.size());
			int first = q.back();
			scc_set.push_back({sz,first});
			if(scc_out){
				scc_out->clear();
			}
			do{
				int n = q.back();

				q.pop_back();
				assert(in_q[n]);
				in_q[n]=false;
				int next;
				if(q.size() && n!=node){
					next = q.back();
				}else{
					next = first;
				}
				scc[n]={sccID,next};
				if(scc_out){
					scc_out->push_back(n);
				}
#ifndef NDEBUG
				{
					UnweightedDijkstra<> d(node,g);
					assert(d.connected(n));
					UnweightedDijkstra<> d2(n,g);
					assert(d2.connected(node));
				}
#endif
				if(n==node)
					break;
			}while(q.size());

		}
	}

	void update( ){
		static int iteration = 0;
		int local_it = ++iteration ;

		if(last_modification>0 && g.modifications==last_modification){
			stats_skipped_updates++;
			return;
		}
		stats_full_updates++;
		
		if(last_deletion==g.deletions){
			stats_num_skipable_deletions++;
		}

		setNodes(g.nodes());
		scc_set.clear();
		q.clear();
		for(int i = 0;i<g.nodes();i++){
			indices[i]=-1;
		}
		int index =0;
		for(int i = 0;i<g.nodes();i++){
			if(indices[i]<0){
				strongConnect(i,index);
			}
		}

		status.setComponents(scc_set.size());
		last_modification=g.modifications;
		last_deletion = g.deletions;
		last_addition=g.additions;

		history_qhead=g.history.size();
		last_history_clear=g.historyclears;
		;
	}

	bool connected(int from, int to){
		update();
		return scc[from].id == scc[to].id;
	}

	int numComponents(){
		 update();
		 return scc_set.size();
	}
	int getComponent(int node){
		update();
		return scc[node].id;
	}
	int getElement(int set){
		update();
		assert(set<scc_set.size());
		return scc_set[set].element;
	}

	bool dbg_uptodate(){
		return true;
	};

	//Compute the SCC for a single node.
	//This is faster for one-shot scc computations from a single source than doing a full update().
	static void getSCC(int node, DynamicGraph & graph, std::vector<int>&scc ){
		TarjansSCC<> s(graph);
		s.setNodes(graph.nodes());
		int index=0;
		s.strongConnect(node,index,&scc);
		assert(scc.size());
		assert(std::count( scc.begin(),scc.end(), node));
	}

};

};


#endif /* TARJANSSCC_H_ */
