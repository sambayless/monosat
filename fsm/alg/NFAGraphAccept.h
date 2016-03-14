/*
 * NFAReach.h
 *
 *  Created on: Dec 16, 2014
 *      Author: sam
 */

#ifndef NFA_GRAPH_ACCEPT_H_
#define NFA_GRAPH_ACCEPT_H_


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
#include <fsm/alg/radix_tree/radix_tree.hpp>
#include <string>
#include <iostream>
#include <stdexcept>
using namespace Monosat;

class rtentry {
public:
	int str_id;
	vec<vec<int>> * strings;
	int       start;
	int 	  end;

	//rtentry() : str_id(-1), prefix_len(0) { }
	rtentry(int strID=-1,vec<vec<int>> * strings=nullptr, int start=0,int _end=-1) : str_id(strID),strings(strings), start(start) {
		if (strings) {
			if (_end < 0 || _end > (*strings)[strID].size()) {
				end = (*strings)[strID].size();
			} else {
				end = _end;
			}
		}
	}
	rtentry(const  rtentry & from):str_id(from.str_id),strings(from.strings),start(from.start),end(from.end){

	}
	int operator[] (int n) const {
		return (*strings)[str_id][n+start];
	}

	bool operator== (const rtentry &rhs) const {
		if(str_id == rhs.str_id && start==rhs.start && end ==rhs.end) {
			return true;
		}else if (size()==rhs.size()){
			int i = 0;
			int j = 0;
			for(;i<size() && j<rhs.size();i++,j++){
				if((*this)[i]<rhs[j])
					return false;
				if(rhs[j]<(*this)[i])
					return false;
			}
			return true;
		}else{
			return false;
		}

	}

	bool operator!= (const rtentry &rhs) const {
		return !((*this)==rhs);
	}
	int size()const{
		return end-start;
	}
	bool operator< (const rtentry &rhs) const {
		if (str_id == rhs.str_id)
			return start < rhs.start || end<rhs.end;
		else{
			int i = 0;
			int j = 0;
			for(;i<size() && j<rhs.size();i++,j++){
				if((*this)[i]<rhs[j])
					return true;
				if(rhs[j]<(*this)[i])
					return true;
			}
			return (i == size()) && (j != rhs.size());
		}

	}
	void print()const{
		printf("str %d [%d,%d]: ",str_id,start,end);
		for(int i = 0;i<size();i++){
			printf("%d", (*this)[i]);
		}
		printf("\n");
	}
};

template<>
inline rtentry radix_substr<rtentry>(const rtentry &entry, int begin, int num)
{
	return rtentry(entry.str_id,entry.strings,entry.start+begin,entry.start+begin+num);
}

template<>
inline rtentry radix_join<rtentry>(const rtentry &entry1, const rtentry &entry2)
{

	if (entry1.end==entry2.start) {
		//only implementing this join for the case where entry1 is a prefix of the underlying string of entry2
#ifndef NDEBUG
		vec<int> &str1 = (*entry1.strings)[entry1.str_id];
		vec<int> &str2 = (*entry2.strings)[entry2.str_id];
		assert(entry1.size() + entry2.size() <= str2.size());
		assert(entry1.end == entry2.start);
		for (int j = entry1.size(); j >= 0; j--) {
			int pos = entry2.start - j - 1;
			assert(pos >= 0);
			int a = str2[pos];
			assert(a == str1[j]);
		}
#endif

		return rtentry(entry2.str_id, entry2.strings, entry2.start - entry1.size(), entry1.size() + entry2.size());
	}else if(entry2.end == entry1.start){
		//only implementing this join for the case where entry1 is a prefix of the underlying string of entry2
#ifndef NDEBUG
		vec<int> &str1 = (*entry1.strings)[entry1.str_id];
		vec<int> &str2 = (*entry2.strings)[entry2.str_id];
		assert(entry1.size() + entry2.size() <= str2.size());
		assert(entry2.end == entry1.start);
		for (int j = entry2.size(); j >= 0; j--) {
			int pos = entry1.start - j - 1;
			assert(pos >= 0);
			int a = str1[pos];
			assert(a == str2[j]);
		}
#endif

		return rtentry(entry1.str_id, entry1.strings, entry1.start - entry2.size(), entry1.size() + entry2.size());
	}else{
		throw std::runtime_error("Radix tree join error");
	}
}

template<>
inline int radix_length<rtentry>(const rtentry &entry)
{
	return entry.size();
}


template<class Status=FSMNullStatus>
class NFAGraphAccept:public NFAAcceptor{
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
	vec<vec<int>> & strings;
	vec<bool> tracked_strings;
	vec<vec<std::pair<int,int>>> tracked_nodes;

	bool checkUsed=true;
	bool hasUsed=false;

	void setStringAccepted(int strID,int acceptingState,bool reachable){
		status.accepts(strID, acceptingState,-1,-1, reachable);
	}

	struct NFAReachStatus {
		NFAGraphAccept & outer;

		void setReachable(int u, bool reachable){
			for(std::pair<int,int> track: outer.tracked_nodes[u]){
				outer.setStringAccepted(track.first,track.second, reachable);
			}
		}
		bool isReachable(int u) const {
			return false;
		}

		void setMininumDistance(int u, bool reachable, int distance){

		}

		NFAReachStatus(NFAGraphAccept & _outer) :
				outer(_outer) {
		}
	};
	UnweightedRamalReps<int,NFAReachStatus> * rr=nullptr;
	//tx_tool::tx * trie=nullptr;
	//cedar::da<int> trie;


	radix_tree<rtentry, int> trie;
public:
	NFAGraphAccept(DynamicFSM & f,int source, vec<vec<int>> & strings,Status & status=fsmNullStatus, bool trackUsedTransitions=false):f(f),status(status),source(source),strings(strings),checkUsed(trackUsedTransitions){


	}


private:

	NFAReachStatus *positiveReachStatus = nullptr;
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
			positiveReachStatus = new NFAReachStatus(*this);
			rr = new UnweightedRamalReps<int,NFAReachStatus> (root->states[source], g,*positiveReachStatus,0,false);
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




public:
	int last_n_strings = 0;
	std::vector<std::vector<int>> t_strings;

	void updateTrie(){
/*		if (last_n_strings>strings.size() ){
			while(t_strings.size()<strings.size()){
				int s = t_strings.size();
				t_strings.push_back();
				for(int i:strings[s]){
					t_strings[s].push_back[i];
				}
				trie[t_strings[s]]=s;
			}
		}*/
	}
	void setTrackStringAcceptance(int str, int state,bool trackPositiveAcceptance, bool trackNegativeAcceptance){
	/*	tracked_strings.growTo(strings.size(),false);
		if(tracked_strings[str])
			return;*/
		addStringAcceptanceCheck(str,state);
		/*vec<int> & string = strings[str];
		rtentry key(str,&strings);
		printf("string %d: ",str);
		for(int i:strings[str]){
			printf("%d",i);
		}
		printf("\n");
		std::vector<typename radix_tree<rtentry, int>::iterator> vect;

		//First check if this key is a prefix of (one or more) existing keys
		trie.prefix_match(key, vect);
		if(vect.size()) {
			const rtentry & prefix_containing_key = vect[0]->first;
			prefix_containing_key.print();

			//the new string is a prefix of an existing string
			int pref_strID = prefix_containing_key.str_id;
			vec<int> & pref_str = strings[pref_strID];
			assert(pref_str.size()>=string.size());
			int i = 0;
			for(;i<pref_str.size();i++){
				if(pref_str[i]!= string[i])
					break;
			}
			//matching prefix length is i



		}else{
			//next check if an existing key is a prefix of this key (and find the longest such key)
			auto it = trie.longest_match(key);
			if(it != trie.end()){
				//found at least one element that is a prefix of this key
				const rtentry & longest_prefix = it->first;
				longest_prefix.print();

			}else{
				trie.greedy_match(key, vect);
				if(vect.size()) {
					const rtentry & best_match = vect[0]->first;
					best_match.print();
					int pref_strID = prefix_containing_key.str_id;
					vec<int> & pref_str = strings[pref_strID];

					int i = 0;
					for(;i<pref_str.size() && i<string.size();i++){
						if(pref_str[i]!= string[i])
							break;
					}

					//longest matching prefix is i

				}
			}
		}

		trie[key] = str;*/
		//tracked_strings[str]=true;


	}



	//inefficient!
	//If state is -1, then this is true if any state accepts the string.
	bool acceptsString(int string, int state){

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


};


#endif /* NFAREACH_H_ */
