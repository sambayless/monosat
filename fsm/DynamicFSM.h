/*
 * DynamicFSM.h
 *
 *  Created on: Dec 15, 2014
 *      Author: sam
 */

#ifndef DYNAMICFSM_H_
#define DYNAMICFSM_H_

#include <vector>
#include "mtl/Vec.h"
#include "mtl/Bitset.h"
#include <algorithm>
#include <cassert>

#include "dgl/DynamicGraph.h"
using namespace dgl;
namespace Monosat {



class DynamicFSM{
	DynamicGraph g;
	std::vector<bool> edge_status;
	bool has_epsilon=true;
	bool is_changed = true;
	vec<Bitset> transitions;
public:
	bool adaptive_history_clear = false;
	long historyClearInterval = 1000;
	int modifications=0;
	int additions=0;
	int deletions=0;
	int n_labels =0;
	long historyclears=0;
	struct EdgeChange {
		bool addition;

		int id;
		int label;
		int mod;
		int prev_mod;

	};
	std::vector<EdgeChange> history;

public:
	DynamicFSM() {

	}

	~DynamicFSM() {

	}

	bool useEmoves()const{
		return has_epsilon;
	}

	bool emove(int edgeID)const{
		return useEmoves() && transitions[edgeID][0];
	}

	unsigned int nLabels()const{
		return n_labels;
	}
	void addLabel(){
		n_labels++;
	}
	bool transitionEnabled(int edgeID, int label)const{
		return transitions[edgeID][label];
	}


	void enableTransition(int edgeID, int label) {
		assert(edgeID >= 0);
		assert(edgeID < g.edges());
		assert(isEdge(edgeID));
		if (transitions[edgeID][label]!= true) {
			transitions[edgeID].set(label);
			//edge_status.setStatus(id,true);
			modifications++;
			additions = modifications;
			history.push_back( { true, edgeID,label, modifications, additions });
		}
	}
	void disableTransition(int edgeID, int label) {
		assert(edgeID >= 0);
		assert(edgeID < g.edges());
		assert(isEdge(edgeID));
		if (edge_status[edgeID] != false) {
			edge_status[edgeID] = false;
			modifications++;
			history.push_back( { false, edgeID,label, modifications, deletions });
			deletions = modifications;
		}
	}

	int addNode() {

		g.addNode();
		modifications++;
		additions = modifications;
		deletions = modifications;
		markChanged();
		clearHistory(true);

		return g.nodes();
	}

	bool edgeEnabled(int edgeID) const {
		return g.edgeEnabled(edgeID);
	}
	bool isEdge(int edgeID) const {
		return g.isEdge(edgeID);
	}
	bool hasEdge(int edgeID) const {
		return isEdge(edgeID);
	}
	//Instead of actually adding and removing edges, tag each edge with an 'enabled/disabled' label, and just expect reading algorithms to check and respect that label.
	int addEdge(int from, int to, int nid = -1) { //, int weight=1
		int id = g.addEdge(from,to,nid);

		modifications++;
		additions = modifications;
		markChanged();

		if(transitions.size()<=id){
			transitions.growTo(id+1);
		}

		return id;
	}

	int nEdgeIDs() {
		return g.nEdgeIDs();
	}
	inline int states() const {
		return g.nodes();
	}
	inline int nodes() const {
		return g.nodes();
	}
	inline int edges() const {
		return g.edges();
	}

	inline int nIncident(int node, bool undirected = false) {
		return g.nIncident(node,undirected);
	}

	inline int nDirectedEdges(int node, bool incoming) {
		return g.nDirectedEdges(node,incoming);
	}
	inline DynamicGraph::Edge & directedEdges(int node, int i, bool is_incoming) {
		return g.directedEdges(node,i,is_incoming);
	}

	inline int nIncoming(int node, bool undirected = false) {
		return g.nIncoming(node,undirected);
	}

	inline DynamicGraph::Edge & incident(int node, int i, bool undirected = false) {
		return g.incident(node,i,undirected);
	}
	inline DynamicGraph::Edge & incoming(int node, int i, bool undirected = false) {
		return g.incoming(node,i,undirected);
	}

	DynamicGraph::FullEdge getEdge(int id) const {
		return g.getEdge(id);
	}

	int getCurrentHistory() {
		return modifications;
	}

	void clearHistory(bool forceClear = false) {
		//long expect=std::max(1000,historyClearInterval*edges());
		if (history.size()
				&& (forceClear
						|| (history.size()
								> (adaptive_history_clear ?
										std::max(1000L, historyClearInterval * edges()) : historyClearInterval)))) {//){
			history.clear();
			historyclears++;

		}
		g.clearHistory();
	}
	//force a new modification
	void invalidate() {
		modifications++;
		additions = modifications;
		modifications++;
		deletions = modifications;
		is_changed = true;

	}

	void markChanged() {
		is_changed = true;

	}
	bool changed() {
		return is_changed;
	}

	void clearChanged() {
		is_changed = false;
		g.clearChanged();
	}
};

}
;



#endif /* DYNAMICFSM_H_ */
