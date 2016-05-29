/*
 * NFAReach.h
 *
 *  Created on: Dec 16, 2014
 *      Author: sam
 */

#ifndef NFA_GRAPH_DAWG_ACCEPT_H_
#define NFA_GRAPH_DAWG_ACCEPT_H_


#include <dgl/DynamicGraph.h>
#include <fsm/alg/NFATypes.h>
#include <fsm/alg/NFAAcceptor.h>
#include <fsm/DynamicFSM.h>
#include <mtl/Vec.h>
#include <mtl/Alg.h>
//#include "mtl/Bitset.h"
#include <cassert>
#include <vector>
#include <dgl/RamalReps.h>
#include <mtl/Bitset.h>
#include <string>
#include <iostream>
#include <stdexcept>
using namespace Monosat;

template<class Status>
class NFAGraphDawgAccept:public NFAAcceptor{
	DynamicFSM & f;
	Status & status;
	DynamicGraph<int> g;

	int last_modification=-1;

	int last_addition=-1;
	int last_deletion=-1;
	int history_qhead=0;
	int last_history_clear=0;

	long stats_full_updates=0;
	long stats_fast_updates=0;
	long stats_fast_failed_updates=0;
	long stats_skip_deletes=0;
	long stats_skipped_updates=0;
	long stats_skipped_string_updates=0;
	long stats_num_skipable_deletions=0;
	double mod_percentage=0;

	double stats_full_update_time=0;
	double stats_fast_update_time=0;

	int source;

	Dawg * root_dawg=nullptr;


	bool hasUsed=false;

	void setDawgAccepted(int node,int acceptingState,bool reachable){
		Dawg * d;
		//look up the
		status.accepts(d, acceptingState,-1,-1, reachable);
	}

	//the graph being checked is set up as a prefix tree; each accepting node corresponds to one dawg node (many accepting graph nodes to each dawg),


	struct NFADawgReachStatus {
		NFAGraphDawgAccept & outer;

		void setReachable(int u, bool reachable){
			for(std::pair<int,int> track: outer.tracked_nodes[u]){
				outer.setDawgAccepted(track.first,track.second, reachable);
			}
		}
		bool isReachable(int u) const {
			return false;
		}

		void setMininumDistance(int u, bool reachable, int distance){

		}

		NFADawgReachStatus(NFAGraphDawgAccept & _outer) :
				outer(_outer) {
		}
	};
	UnweightedRamalReps<int,NFADawgReachStatus> * rr=nullptr;

public:
	NFAGraphDawgAccept(DynamicFSM & f,int source, Dawg * root_dawg,Status & status=fsmNullStatus, bool trackUsedTransitions=false):f(f),status(status),source(source),root_dawg(root_dawg){


	}


private:

	NFADawgReachStatus *positiveReachStatus = nullptr;
	vec<vec<vec<int>>> edge_map;
	//vec<vec<int>> states;

	vec<NFATransition> reverse_edge_map;

	struct ToCheck{
		int str;
		int accepting_state;
	};
	vec<ToCheck> acceptancesToCheck;
	//graph states for each unrolled string fsm are held in a tree structure


	struct UnrolledStep{
		UnrolledStep * parent=nullptr;
		vec<UnrolledStep*> children;//alphabet size of these
		int l=-1;
		int depth =0;
		vec<int> states;
	};
	vec<UnrolledStep*> string_last_nodes;
	UnrolledStep * root=nullptr;
	//int startNode =-1;

	void addStringAcceptanceCheck(int str, int acceptingState){
		string_last_nodes.growTo(strings.size(),nullptr);
		if (!root){
			edge_map.growTo(f.edges());
			//map from edge ids in f
			for(int edgeID = 0;edgeID<f.edges();edgeID++){
				if(f.hasEdge(edgeID)){
					edge_map[edgeID].growTo(f.inAlphabet());

				}
			}

			//startNode = g.addNode();
			root = new UnrolledStep();
			root->children.growTo(f.inAlphabet(),nullptr);

			for(int i = 0;i<f.states();i++){
				root->states.push(g.addNode());
			}
			if(f.hasEmoves()){
				for(int s = 0;s<f.states();s++){
					for(int j = 0;j<f.nIncoming(s);j++){
						int edgeID = f.incoming(s,j).id;
						int from = f.incoming(s,j).node;
						//if(f.transitionEnabled(edgeID,0,0)){
							int newEdgeID = g.addEdge(root->states[from],root->states[s]);
							edge_map[edgeID][0].push(newEdgeID);
							reverse_edge_map.growTo(newEdgeID+1);
							reverse_edge_map[newEdgeID].edgeID = edgeID;
							reverse_edge_map[newEdgeID].input = 0;
							reverse_edge_map[newEdgeID].output=0;
						//}
						if (!f.transitionEnabled(edgeID, 0, 0)){
							g.disableEdge(newEdgeID);
						}
					}
				}
			}
			positiveReachStatus = new NFADawgReachStatus(*this);
			rr = new UnweightedRamalReps<int,NFADawgReachStatus> (root->states[source], g,*positiveReachStatus,0,false);
		}



		vec<int> & string = strings[str];


		//build nodes

		UnrolledStep * parent = root;
		if(!string_last_nodes[str]) {
			for (int p = 1; p <= string.size(); p++) {
				int l = string[p - 1];

				if (parent->children[l]) {
					parent = parent->children[l];
					continue;
				}
				UnrolledStep *child = new UnrolledStep();
				child->l = l;
				child->parent = parent;
				child->children.growTo(f.inAlphabet(), nullptr);
				child->depth = parent->depth + 1;
				parent->children[l] = child;


				for (int s = 0; s < f.states(); s++) {
					child->states.push(g.addNode());

					for (int j = 0; j < f.nIncoming(s); j++) {
						int edgeID = f.incoming(s, j).id;
						int from = f.incoming(s, j).node;
						//if (f.transitionEnabled(edgeID, l, 0)) {
							int newEdgeID = g.addEdge(parent->states[from], child->states[s]);
							edge_map[edgeID][l].push(newEdgeID);
							reverse_edge_map.growTo(newEdgeID+1);
							reverse_edge_map[newEdgeID].edgeID = edgeID;
							reverse_edge_map[newEdgeID].input = l;
							reverse_edge_map[newEdgeID].output=0;
						if (!f.transitionEnabled(edgeID, l, 0)){
							g.disableEdge(newEdgeID);
						}

						//}
					}
				}
				if (f.hasEmoves()) {
					for (int s = 0; s < f.states(); s++) {
						for (int j = 0; j < f.nIncoming(s); j++) {
							int edgeID = f.incoming(s, j).id;
							int from = f.incoming(s, j).node;
							//if (f.transitionEnabled(edgeID, 0, 0)) {
								int newEdgeID = g.addEdge(child->states[from], child->states[s]);
								edge_map[edgeID][0].push(newEdgeID);
								reverse_edge_map.growTo(newEdgeID+1);
								reverse_edge_map[newEdgeID].edgeID = edgeID;
								reverse_edge_map[newEdgeID].input = 0;
								reverse_edge_map[newEdgeID].output=0;
							if (!f.transitionEnabled(edgeID, 0, 0)){
								g.disableEdge(newEdgeID);
							}
							//}
						}
					}
				}
				parent = child;
			}
		}else{
			parent=string_last_nodes[str];
		}
		tracked_nodes.growTo(g.nodes());
		string_last_nodes[str]=parent;
		bool alreadyTracked=false;
		for(std::pair<int,int>  track: tracked_nodes[parent->states[acceptingState]]){
			if(track.first==str && track.second==acceptingState){
				alreadyTracked=true;
				break;
			}
		}
		if(!alreadyTracked) {
			acceptancesToCheck.push({str, parent->states[acceptingState]});
		}
		tracked_nodes[parent->states[acceptingState]].push({str,acceptingState});
		rr->update();
		status.accepts(str, acceptingState,-1,-1, rr->connected(parent->states[acceptingState]));

	}

public:

	void update(){

		if (last_modification > 0 && f.modifications == last_modification){
			return;
		}

		if (last_history_clear != f.historyclears) {
			history_qhead = f.historySize();;
			last_history_clear = f.historyclears;
			for (int edgeid = 0; edgeid < f.edges(); edgeid++) {
				for(int l = 0;l<f.inAlphabet();l++){
					if (f.transitionEnabled(edgeid,l,0)) {
						for (int graphEdgeID:edge_map[edgeid][l]) {
							g.enableEdge(graphEdgeID);
						}
					}else {
						for(int graphEdgeID:edge_map[edgeid][l]){
							g.disableEdge(graphEdgeID);
						}
					}
				}
			}
		}

		for (int i = history_qhead; i < f.historySize(); i++) {
			DynamicFSM::EdgeChange & c = f.history[i];
			assert(c.output==0);
			if(c.addition && f.transitionEnabled(c.id,c.input,c.output)){
				for(int edgeID:edge_map[c.id][c.input]){
					g.enableEdge(edgeID);
				}
			}else if (!c.addition && !f.transitionEnabled(c.id,c.input,c.output)){
				for(int edgeID:edge_map[c.id][c.input]){
					g.disableEdge(edgeID);
				}
			}
		}
		rr->update();


		last_modification = f.modifications;
		last_deletion = f.deletions;
		last_addition = f.additions;

		history_qhead = f.historySize();

		last_history_clear = f.historyclears;
	}

	int numUpdates() const {
		if(!rr)
			return 0;
		return rr->numUpdates();
	}



public:





	//inefficient!
	//If state is -1, then this is true if any state accepts the string.
	bool acceptsDawg(Dawg * dawg, int state){

		update();
		int s =string_last_nodes[string]->states[state];
		return rr->connected(s);

	}
	bool getPath(int string, int state, vec<NFATransition> &path) {
		update();
		int s =string_last_nodes[string]->states[state];
		if(! rr->connected(s)){
			return false;
		}

		//need to map the edges of the unrolled graph to transitions in the fsm
		while(s != root->states[source]){
			int edgeID = rr->incomingEdge(s);
			path.push(reverse_edge_map[edgeID]);
			int p = rr->previous(s);
			s = p;
		}
		reverse(path);
		return true;
	}
	vec<Bitset> used_transition;
	bool getAbstractPath(int string, int state, vec<NFATransition> &path, bool reversed) {
		update();
		int s =string_last_nodes[string]->states[state];
		if(! rr->connected(s)){
			return false;
		}
		used_transition.growTo(f.edges());
		//need to map the edges of the unrolled graph to transitions in the fsm
		while(s != root->states[source]){
			int edgeID = rr->incomingEdge(s);

			NFATransition & transition = reverse_edge_map[edgeID];
			used_transition[transition.edgeID].growTo(f.inAlphabet());
			int l = transition.input;
			if(!used_transition[transition.edgeID][l]) {
				path.push(transition);
				used_transition[transition.edgeID].set(transition.input);
			}
			int p = rr->previous(s);
			s = p;
		}

		if(!reversed)
			reverse(path);

		for(NFATransition & t:path){
			used_transition[t.edgeID].clear(t.input);
		}
		return true;
	}


};


#endif /* NFAREACH_H_ */
