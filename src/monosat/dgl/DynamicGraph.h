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

#ifndef DYNAMICGRAPH_H_
#define DYNAMICGRAPH_H_
#include <vector>
#include <algorithm>
#include <cassert>
#include <sstream>
#include <cstdint>
#include <sstream>
#include <cstdio>

namespace dgl {

class DynamicGraphAlgorithm{
public:
	virtual ~DynamicGraphAlgorithm(){}

	virtual void updateHistory()=0;
};

/**
 * A dynamic graph.
 * It supports efficiently recomputing graph properties as edges are added and removed ('enabled' and 'disabled').
 *
 * DynamicGraph expects all edges (and all nodes) that it will ever use to be added in advance, after which
 * those (already declared) edges can be efficiently enabled or disabled.
 *
 * Although it also allows adding new edges (and new nodes) dynamically, this library is not optimized for that use case;
 * adding new edges or nodes (as opposed to enabling or disabling existing edges) will typically cause all properties to be
 * recomputed from scratch.
 *
 * Most algorithms in the library are optimized for moderate sized, sparsely connected graphs (<10,000 edges/nodes).
 */
template<typename Weight>
class DynamicGraph {
	
	std::vector<bool> edge_status;
	std::vector<bool> edge_status_const;
	std::vector<Weight> weights;
	int num_nodes=0;
	int num_edges=0;
	int next_id=0;
	bool is_changed=false;
	std::vector<DynamicGraphAlgorithm*> dynamic_algs;
	std::vector<int> dynamic_history_pos;

	int64_t history_offset=0;

public:
	bool disable_history_clears=false;
	int dynamic_history_clears=0;

	bool adaptive_history_clear = false;
	int64_t historyClearInterval = 1000;
	int modifications=0;
	int additions=0;
	int deletions=0;
	int edge_increases = 0;
	int edge_decreases = 0;
	int64_t previous_history_size=0;
	int64_t historyclears=0;
	int64_t skipped_historyclears=0;
	struct Edge {
		int node;
		int id;
	};
	std::vector<std::vector<Edge> > adjacency_list;
	std::vector<std::vector<Edge> > inverted_adjacency_list;
	std::vector<std::vector<Edge> > adjacency_undirected_list;
public:
	struct FullEdge {
		int from;
		int to;
		int id;
		//int weight;
		FullEdge() :
				from(-1), to(-1), id(-1) {
		} //,weight(1){}
		FullEdge(int from, int to, int id) :
				from(from), to(to), id(id) {
		} //,weight(weight){}
	};

private:
#ifdef DEBUG_DGL
public:
#endif
	std::vector<FullEdge> all_edges;
public:
	struct EdgeChange {
		bool addition;
		bool deletion;
		bool weight_increase;
		bool weight_decrease;
		int id;
		int mod;
		int prev_mod;
	};
private:
	std::vector<EdgeChange> history;
public:
	//Logfile information if recording is enabled.
	FILE * outfile=nullptr;

	DynamicGraph() {
	}

	~DynamicGraph() {

	}

	void addNodes(int n) {
		for (int i = 0; i < n; i++)
			addNode();
	}
	//Returns true iff the edge exists and is a self loop
	inline bool selfLoop(int edgeID)  {
		return hasEdge(edgeID) && getEdge(edgeID).from == getEdge(edgeID).to;
	}
	//SLOW!
	bool hasEdge(int from, int to) const {
		for (int i = 0; i < adjacency_list[from].size(); i++) {
			if (adjacency_list[from][i].node == to && edgeEnabled(adjacency_list[from][i].id)) {
				return true;
			}
		}
		return false;
	}
	//SLOW! Returns -1 if there is no edge
	int getEdge(int from, int to) const {
		for (int i = 0; i < adjacency_list[from].size(); i++) {
			if (adjacency_list[from][i].node == to && edgeEnabled(adjacency_list[from][i].id)) {
				return adjacency_list[from][i].id;
			}
		}
		return -1;
	}
	bool hasEdgeUndirected(int from, int to) const {
		for (int i = 0; i < adjacency_undirected_list[from].size(); i++) {
			if (adjacency_undirected_list[from][i].node == to && edgeEnabled(adjacency_undirected_list[from][i].id)) {
				return true;
			}
		}
		return false;
	}

	int addNode() {

		adjacency_list.push_back( { }); //adj list
		adjacency_undirected_list.push_back( { });
		inverted_adjacency_list.push_back( { });
		modifications++;
		additions = modifications;
		deletions = modifications;
		edge_increases = modifications;
		edge_decreases = modifications;
		markChanged();
		clearHistory(true);

		if (outfile) {
			fprintf(outfile, "node %d\n", num_nodes);
			fflush(outfile);
		}

		return num_nodes++;
	}
	//true iff the edge's current assignment will never be altered.
	bool isConstant(int edgeID) const{
		return edge_status_const[edgeID];
	}

	void makeEdgeAssignmentConstant(int edgeID){
		edge_status_const[edgeID]=true;
	}

	bool edgeEnabled(int edgeID) const {
		assert(edgeID < edge_status.size());
		return edge_status[edgeID];
	}
	bool isEdge(int edgeID) const {
		return edgeID < all_edges.size() && all_edges[edgeID].id == edgeID;
	}
	bool hasEdge(int edgeID) const {
		return isEdge(edgeID);
	}
	//Instead of actually adding and removing edges, tag each edge with an 'enabled/disabled' label, and just expect reading algorithms to check and respect that label.
	int addEdge(int from, int to, int id = -1, Weight weight=1){
		assert(from < num_nodes);
		assert(to < num_nodes);
		assert(from >= 0);
		assert(to >= 0);
		if (id < 0) {
			id = next_id++;
		} else {
			if (id >= next_id) {
				next_id = id + 1;
			}
		}

		num_edges = next_id;
		adjacency_list[from].push_back( { to, id });
		adjacency_undirected_list[from].push_back( { to, id });
		adjacency_undirected_list[to].push_back( { from, id });
		if (edge_status.size() <= id)
			edge_status.resize(id + 1);

		if(edge_status_const.size()<=id){
			edge_status_const.resize(id+1,false);
		}

		inverted_adjacency_list[to].push_back( { from, id });
		if (all_edges.size() <= id)
			all_edges.resize(id + 1);
		all_edges[id]= {from,to,id}; //,weight};
		if(weights.size()<=id)
			weights.resize(id+1,0);
		weights[id]=weight;

		modifications++;
		additions = modifications;
		edge_increases = modifications;
		markChanged();


		if (outfile) {
			fprintf(outfile, "edge %d %d %d %d\n", from, to, 1, id + 1);
			std::stringstream ss;
			ss<<weight;
			fprintf(outfile, "edge_weight %d %s\n", id + 1, ss.str().c_str());
			fflush(outfile);
		}

//		history.push_back({true,id,modifications});
		enableEdge(from, to, id);		//default to enabled

		return id;
	}
	int nEdgeIDs() {
		assert(num_edges == all_edges.size());
		return num_edges;		//all_edges.size();
	}
	inline int nodes() const {
		return num_nodes;
	}
	inline int edges() const {
		return num_edges;
	}

	inline int nIncident(int node, bool undirected = false) {
		assert(node >= 0);
		assert(node < nodes());
		if (undirected) {
			return adjacency_undirected_list[node].size();
		} else {
			return adjacency_list[node].size();
		}
	}

	inline int nDirectedEdges(int node, bool incoming) {
		assert(node >= 0);
		assert(node < nodes());
		if (incoming) {
			return nIncoming(node, false);
		} else {
			return nIncident(node, false);
		}
	}
	inline Edge & directedEdges(int node, int i, bool is_incoming) {
		if (is_incoming) {
			return incoming(node, i, false);
		} else {
			return incident(node, i, false);
		}
	}

	inline int nIncoming(int node, bool undirected = false) {
		assert(node >= 0);
		assert(node < nodes());
		if (undirected) {
			return adjacency_undirected_list[node].size();
		} else {
			return inverted_adjacency_list[node].size();
		}
	}

	inline Edge & incident(int node, int i, bool undirected = false) {
		assert(node >= 0);
		assert(node < nodes());
		assert(i < nIncident(node, undirected));
		if (undirected) {
			return adjacency_undirected_list[node][i];
		} else {
			return adjacency_list[node][i];
		}
	}
	inline Edge & incoming(int node, int i, bool undirected = false) {
		assert(node >= 0);
		assert(node < nodes());
		assert(i < nIncoming(node, undirected));
		if (undirected) {
			return adjacency_undirected_list[node][i];
		} else {
			return inverted_adjacency_list[node][i];
		}
	}
	std::vector<FullEdge> & getEdges(){
		return all_edges;
	}

	std::vector<Weight> & getWeights(){
		return weights;
	 }
/*	 Weight getWeight(int edgeID){
		 return weights[edgeID];
	 //return all_edges[edgeID].weight;
	 }*/
	 Weight  getWeight(int edgeID){
		 return weights[edgeID];
	 //return all_edges[edgeID].weight;
	 }
	FullEdge & getEdge(int id)  {
		return all_edges[id];
	}
	void setEdgeEnabled(int id, bool enable){
		if(enable){
			enableEdge(id);
		}else{
			disableEdge(id);
		}
	}

	void enableEdge(int id) {
		enableEdge(all_edges[id].from, all_edges[id].to, id);
	}
	void disableEdge(int id) {
		disableEdge(all_edges[id].from, all_edges[id].to, id);
	}
	void enableEdge(int from, int to, int id) {
		assert(id >= 0);
		assert(id < edge_status.size());
		assert(isEdge(id));
		if (!edge_status[id] ) {
			edge_status[id] = true;
			//edge_status.setStatus(id,true);

			modifications++;
			additions = modifications;
			history.push_back( { true,false,false,false, id, modifications, additions });

			if (outfile) {

				fprintf(outfile, "%d\n", id + 1);
				fflush(outfile);
			}

		}
	}

	bool undoEnableEdge(int id) {
		assert(id >= 0);
		assert(id < edge_status.size());
		assert(isEdge(id));
		if (!history.size())
			return false;

		if (history.back().addition && history.back().id == id && history.back().mod == modifications) {
			//edge_status.setStatus(id,false);
			edge_status[id] = false;

			if (outfile) {
				if (id + 1 == 89486) {
					int a = 1;
				}
				fprintf(outfile, "-%d\n", id + 1);
				fflush(outfile);
			}

			modifications--;
			additions = history.back().prev_mod;
			history.pop_back();
			return true;
		}
		return false;
	}

	void disableEdge(int from, int to, int id) {
		assert(id >= 0);
		assert(id < edge_status.size());
		assert(isEdge(id));
		if (edge_status[id]) {
			//edge_status.setStatus(id,false);
			edge_status[id] = false;

			if (outfile) {
				fprintf(outfile, "-%d\n", id + 1);
				fflush(outfile);
			}

			modifications++;

			history.push_back( { false,true,false,false, id, modifications, deletions });
			deletions = modifications;
		}
	}

	bool undoDisableEdge(int id) {
		assert(id >= 0);
		assert(id < edge_status.size());
		assert(isEdge(id));
		if (!history.size())
			return false;

		if (!history.back().addition && history.back().id == id && history.back().mod == modifications) {
			//edge_status.setStatus(id,true);
			edge_status[id] = true;

			if (outfile) {

				fprintf(outfile, "%d\n", id + 1);
				fflush(outfile);
			}

			modifications--;
			deletions = history.back().prev_mod;
			history.pop_back();
			return true;
		}
		return false;
	}
	inline Weight  getEdgeWeight(int edgeID){
		return getWeight(edgeID);
	}
	void setEdgeWeight(int id,const Weight & w) {
			assert(id >= 0);
			assert(id < edge_status.size());
			assert(isEdge(id));
			if(w==getWeight(id)){
				return;
			}

			modifications++;
			if(w>getWeight(id)){
				history.push_back( {false,false, true,false, id, modifications, additions });
				edge_increases = modifications;
			}else{
				assert(w<getWeight(id));
				history.push_back( {false,false, false, true, id, modifications, additions });
				edge_decreases = modifications;
			}
			weights[id]=w;

			if (outfile) {
				std::stringstream ss;
				ss<<w;
				fprintf(outfile, "edge_weight %d %s\n", id + 1, ss.str().c_str());
				fflush(outfile);
			}


		}


	void drawFull(bool showWeights = false, bool force_draw=false) {
#ifndef DEBUG_DGL
		if(!force_draw)
			return;
#endif
//#ifdef DEBUG_DGL
		printf("digraph{\n");
		for (int i = 0; i < num_nodes; i++) {
			printf("n%d\n", i);
		}

		for (int i = 0; i < adjacency_list.size(); i++) {
			for (int j = 0; j < adjacency_list[i].size(); j++) {
				int id = adjacency_list[i][j].id;
				int u = adjacency_list[i][j].node;
				const char * s = "black";
				if (edgeEnabled(id))
					s = "red";
				else
					s = "blue";
				if(showWeights){
					std::stringstream ss;
					ss<<getWeight(id);
					printf("n%d -> n%d [label=\"v%d w=%s\",color=\"%s\"]\n", i,u, id,ss.str().c_str(), s);
				}else{
				printf("n%d -> n%d [label=\"v%d\",color=\"%s\"]\n", i, u, id, s);
				}
			}
		}
		printf("}\n");
//#endif

	}

	bool rewindHistory(int steps) {

		int cur_modifications = modifications;
		for (int i = 0; i < steps; i++) {
			EdgeChange & e = history.back();
			if (e.addition) {
				if (!undoEnableEdge(e.id)) {
					return false;
				}
			} else if(e.deletion) {
				if (!undoDisableEdge(e.id)) {
					return false;
				}
			}

		}
		assert(modifications == cur_modifications - steps);
		return true;
	}
	/**
	 * Returns a unique identifier for this algorithm.
	 */
	int addDynamicAlgorithm(DynamicGraphAlgorithm*alg){
		dynamic_algs.push_back(alg);
		dynamic_history_pos.push_back(0);
		//n_dynamic_algs_updtodate+= (historySize()==0);
		return dynamic_algs.size()-1;
	}

	void updateAlgorithmHistory(DynamicGraphAlgorithm * alg, int algorithmID, int historyPos){
		assert(dynamic_algs[algorithmID]==alg);//sanity check
		//bool was_uptodate = dynamic_history_pos[algorithmID]==historySize();
		dynamic_history_pos[algorithmID]=historyPos;
/*		if(!was_uptodate && historyPos ==historySize()){
			dynamic_history_pos[algorithmID]++;
		}*/
	}

	EdgeChange & getChange(int64_t historyPos){
		assert(historyPos-history_offset>=0);
		assert(historyPos-history_offset<history.size());
		return history[historyPos-history_offset];
	}

	int historySize(){
		return history.size() + history_offset;
	}

	int getCurrentHistory() {
		return modifications;
	}



	void clearHistory(bool forceClear = false) {
		//int64_t expect=std::max(1000,historyClearInterval*edges());
		//check whether we can do a cheap history cleanup (without resetting all the dynamic algorithms)
		if(disable_history_clears)
			return;

		if (history.size()
				&& (forceClear
						|| (history.size()
								>= (std::min((int64_t)history.max_size(), (adaptive_history_clear ?
										std::max(1000L, historyClearInterval * edges()) : historyClearInterval)))))) {//){


			if(!forceClear && dynamic_history_clears>0){
				int n_uptodate=0;
				for(int algorithmID = 0;algorithmID<dynamic_algs.size();algorithmID++){
					if (dynamic_history_pos[algorithmID]!=historySize()){
						if(dynamic_history_clears==2){
							dynamic_algs[algorithmID]->updateHistory();
							if (dynamic_history_pos[algorithmID]==historySize()){
								n_uptodate++;
							}
						}
					}else{
						n_uptodate++;
					}
				}

				if(n_uptodate==dynamic_algs.size()){
					//we can skip this history clear.
					history_offset=historySize();
					skipped_historyclears++;
					history.clear();
					assert(history_offset == historySize());
					return;
				}
			}
			previous_history_size=history.size();
			history_offset=0;
			history.clear();
			historyclears++;

			if (outfile) {
				fprintf(outfile, "clearHistory\n");
				fflush(outfile);
			}

		}
	}
	//force a new modification
	void invalidate() {
		modifications++;
		additions = modifications;
		modifications++;
		deletions = modifications;
		is_changed = true;

		if (outfile) {
			fprintf(outfile, "invalidate\n");
			fflush(outfile);
		}

	}
	int64_t getPreviousHistorySize()const{
		return previous_history_size;
	}
	void markChanged() {
		is_changed = true;

		if (outfile) {
			fprintf(outfile, "markChanged\n");
			fflush(outfile);
		}

	}
	bool changed() {
		return is_changed;
	}
	
	void clearChanged() {
		is_changed = false;

		if (outfile) {
			fprintf(outfile, "clearChanged\n");
			fflush(outfile);
		}

	}

	void clear(){
		edge_status.clear();
		num_nodes=0;
		num_edges=0;
		next_id=0;


		adjacency_list.clear();
		inverted_adjacency_list.clear();
		adjacency_undirected_list.clear();
		all_edges.clear();
		history.clear();
		invalidate();
		clearHistory(true);
	}
	void copyTo(DynamicGraph & to){
		to.clear();

		to.num_nodes = num_nodes;
		to.num_edges = num_edges;
		to.next_id = next_id;
		to.edge_status = edge_status;
		to.historyClearInterval=historyClearInterval;
		to.adjacency_list = adjacency_list;
		to.adjacency_undirected_list=adjacency_undirected_list;
		to.all_edges =all_edges;
		to.inverted_adjacency_list=inverted_adjacency_list;


	}

	void crc(){
		unsigned long prime = 31;
		unsigned long result = 1;

		for (int i = 0;i<edges();i++){
			if (edgeEnabled(i)){
				result = prime * result + i;
			}
		}
		printf("Graph CRC: %ld",result);
	}


};

}
;
#endif /* DYNAMICGRAPH_H_ */
