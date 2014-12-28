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
	//std::vector<Bitset> edge_status;
	bool has_epsilon=true;
	bool is_changed = true;
public:
	vec<Bitset> transitions;

	bool adaptive_history_clear = false;
	long historyClearInterval = 1000;
	int modifications=0;
	int additions=0;
	int deletions=0;
	int in_alphabet =1;
	int out_alphabet=1;
	long historyclears=0;
	struct EdgeChange {
		bool addition;

		int id;
		int input;
		int output;
		int mod;
		int prev_mod;

	};
	std::vector<EdgeChange> history;

public:
	DynamicFSM() {

	}

	~DynamicFSM() {

	}

	void setEmovesEnabled(bool enabled){
		has_epsilon=enabled;
	}

	bool emovesEnabled()const{
		return has_epsilon;
	}

/*
	bool emove(int edgeID)const{
		return emovesEnabled() && transitions[edgeID][0];
	}
*/

	unsigned int inAlphabet()const{
		return in_alphabet;
	}
	unsigned int outAlphabet()const{
		return out_alphabet;
	}
	void addInCharacter(){
		in_alphabet++;
	}
	void addOutCharacter(){
		out_alphabet++;
	}
	bool transitionEnabled(int edgeID, int input, int output)const{
		assert(input<inAlphabet());
		assert(output<outAlphabet());
		int pos = input +output*inAlphabet();
		return transitions[edgeID][pos];
	}

	int addTransition(int from, int to,int edgeID, int input,int output, bool defaultEnabled=true){
		assert(input<inAlphabet());
		assert(output<outAlphabet());
		while(from>=g.nodes() || to>=g.nodes())
			g.addNode();
		if(edgeID==-1){
			edgeID = g.addEdge(from, to, edgeID);
		}
		transitions.growTo(edgeID+1);
		transitions[edgeID].growTo(inAlphabet()*outAlphabet());
		int pos = input +output*inAlphabet();
		if(defaultEnabled)
			transitions[edgeID].set(pos);
		return edgeID;
	}

	void enableTransition(int edgeID, int input,int output) {
		assert(edgeID >= 0);
		assert(edgeID < g.edges());
		assert(isEdge(edgeID));
		int pos = input +output*inAlphabet();
		if (!transitions[edgeID][pos]) {
			transitions[edgeID].set(pos);
			//edge_status.setStatus(id,true);
			modifications++;
			additions = modifications;
			history.push_back( { true, edgeID,input,output, modifications, additions });
		}
	}
	void disableTransition(int edgeID, int input,int output) {
		assert(edgeID >= 0);
		assert(edgeID < g.edges());
		assert(isEdge(edgeID));
		int pos = input +output*inAlphabet();
		if (transitions[edgeID][pos]) {
			transitions[edgeID].clear(pos);
			modifications++;
			history.push_back( { false, edgeID,input,output, modifications, deletions });
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

	bool hasEdge(int from, int to)const{

		return g.hasEdge(from,to);
	}
	int getEdge(int from,int to)const{
		return g.getEdge(from,to);
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

	void draw(int source=-1, int dest=-1){

		printf("digraph{\n");
		if(source>=0){
			printf("start->%d\n",source);
		}
		if(dest>=0){
			printf("%d [shape=doublecircle]\n",dest);
		}
		for(int i = 0;i<transitions.size();i++){
			bool any_enabled=false;
			for(int l= 0;l<transitions[i].size();l++){
				if(transitions[i][l]){
					any_enabled=true;
					break;
				}
			}
			if (any_enabled){
				printf("%d->%d [label=\"", g.getEdge(i).from,g.getEdge(i).to);

				for(int in = 0;in<inAlphabet();in++){
					for(int out = 0;out<outAlphabet();out++){
						int pos = in + inAlphabet()*out;
						if(transitions[i][pos]){
							if(out==0){
								if(in==0){
									printf("{},");
								}else{
									printf("%c:,",'a'+in-1);
								}
							}else{
								if(in==0){
									printf(":%c,",'a'+out-1);
								}else{
									printf("%c:%c,",'a'+in-1,'a'+out-1);
								}
							}

						}
					}
				}


				printf("\"]\n");
			}
		}


		printf("}\n");

	}

};

}
;



#endif /* DYNAMICFSM_H_ */
