/*
 * NFAReach.h
 *
 *  Created on: Dec 16, 2014
 *      Author: sam
 */

#ifndef NFA_GRAPH_DAWG_ACCEPT_H_
#define NFA_GRAPH_DAWG_ACCEPT_H_

#include "monosat/fsm/alg/Dawg.h"
#include "monosat/dgl/DynamicGraph.h"
#include "monosat/fsm/alg/NFATypes.h"
#include "monosat/fsm/alg/NFAAcceptor.h"
#include "monosat/fsm/DynamicFSM.h"
#include "monosat/mtl/Vec.h"
#include "monosat/mtl/Alg.h"
//#include "mtl/Bitset.h"
#include <cassert>
#include <vector>
#include "monosat/dgl/RamalReps.h"
#include "monosat/mtl/Bitset.h"
#include <string>
#include <iostream>
#include <stdexcept>
using namespace Monosat;

template<class Status>
class NFAGraphDawgAccept{
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

	FSMDawg * root_dawg=nullptr;


	bool hasUsed=false;
	vec<FSMDawg*> trackedDawgs;
	void setDawgAccepted(int dawgID,int acceptingState,bool reachable){
		assert(dawgID>=0);
		FSMDawg * d = trackedDawgs[dawgID];
		//look up the
		status.accepts(d, acceptingState,-1,-1, reachable);
	}

	//the graph being checked is set up as a prefix tree; each accepting node corresponds to one dawg node (many accepting graph nodes to each dawg),

	//one for every node in the constructed dynamic graph
	/*struct AcceptanceNode{
		//FSMDawg * d =nullptr;
	};*/
	struct TrackedNode{
		int dawg_id=-1;
		int state=-1;
	};
	vec<TrackedNode> tracked_node_dawg_ids;

	struct NFADawgReachStatus {
		NFAGraphDawgAccept & outer;

		void setReachable(int u, bool reachable){
			/*for(std::pair<int,int> track: outer.tracked_nodes[u]){
				outer.setDawgAccepted(track.first,track.second, reachable);
			}*/
			if(outer.tracked_node_dawg_ids[u].dawg_id>=0){
				outer.setDawgAccepted(outer.tracked_node_dawg_ids[u].dawg_id,outer.tracked_node_dawg_ids[u].state,reachable);
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
	NFAGraphDawgAccept(DynamicFSM & f,int source, Status & status, bool trackUsedTransitions=false):f(f),status(status),source(source){


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
	vec<UnrolledStep*> dawg_last_nodes;
	UnrolledStep * root=nullptr;
	//int startNode =-1;


	void addDawgAcceptanceCheck(FSMDawg * dawg, int acceptingState){
		//dawg_last_nodes.growTo(strings.size(),nullptr);
		int dawgID = dawg->id;
		assert(dawgID>=0);
		trackedDawgs.growTo(dawgID+1,nullptr);
		assert(trackedDawgs[dawgID]==nullptr);

		trackedDawgs[dawgID] = dawg;

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


		//build nodes
		vec<UnrolledStep*> last_nodes;
		UnrolledStep * parent = root;
		//if(!dawg_last_nodes[str]) {
		for(int l = 0;l<dawg->transitions.size();l++) {
			if(dawg->transitions[l])
				buildStep(dawg->transitions[l],l, parent, last_nodes);
		}

		tracked_node_dawg_ids.growTo(g.nodes());
		for(UnrolledStep * last: last_nodes) {
			int acceptingNode = last->states[acceptingState];
			tracked_node_dawg_ids[acceptingNode].dawg_id = dawgID;
			tracked_node_dawg_ids[acceptingNode].state = acceptingState;
		}


			//for (int p = 1; p <= string.size(); p++) {

			//}
		/*}else{
			parent=dawg_last_nodes[str];
		}*/
		//tracked_nodes.growTo(g.nodes());
		//dawg_last_nodes[str]=parent;
		/*bool alreadyTracked=false;
		for(std::pair<int,int>  track: tracked_nodes[parent->states[acceptingState]]){
			if(track.first==str && track.second==acceptingState){
				alreadyTracked=true;
				break;
			}
		}
		if(!alreadyTracked) {
			acceptancesToCheck.push({str, parent->states[acceptingState]});
		}
		tracked_nodes[parent->states[acceptingState]].push({str,acceptingState});*/
		rr->update();
		status.accepts(dawg, acceptingState,-1,-1, rr->connected(parent->states[acceptingState]));

	}
	void buildStep(FSMDawg * dawg, int l, UnrolledStep * parent, vec<UnrolledStep*> & last_nodes){
		//int l = string[p - 1];
		if(!dawg)
			return;
		if (parent->children[l]) {
			parent = parent->children[l];
			return;
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
		bool any_transitions=false;
		for(int l = 0;l<dawg->transitions.size();l++) {
			if(dawg->transitions[l]) {
				any_transitions=true;
				buildStep(dawg->transitions[l], l, child,last_nodes);
			}
		}
		if(!any_transitions){
			last_nodes.push(child);
		}

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



	void setTrackDawgAcceptance(FSMDawg * d, int state,bool trackPositiveAcceptance, bool trackNegativeAcceptance){
		addDawgAcceptanceCheck(d,state);

	}



	//inefficient!
	//If state is -1, then this is true if any state accepts the string.
	bool acceptsDawg(FSMDawg * dawg, int state){
		int dawgID = dawg->id;
		assert(dawgID>=0);
		update();
		int s =dawg_last_nodes[dawg->id]->states[state];
		return rr->connected(s);

	}
	bool getPath(FSMDawg * dawg, int state, vec<NFATransition> &path) {
		update();
		int s =dawg_last_nodes[dawg->id]->states[state];
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
	bool getAbstractPath(FSMDawg * dawg, int state, vec<NFATransition> &path, bool reversed) {
		update();
		int dawgID = dawg->id;
		int s =dawg_last_nodes[dawgID]->states[state];
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
