/*
 * NFAReach.h
 *
 *  Created on: Dec 16, 2014
 *      Author: sam
 */

#ifndef NFAREACH_H_
#define NFAREACH_H_

#include <dgl/DynamicGraph.h>
#include <fsm/alg/NFATypes.h>
#include <fsm/DynamicFSM.h>
#include <mtl/Vec.h>
//#include "mtl/Bitset.h"
#include <cassert>
#include <vector>

using namespace Monosat;


template<class Status=FSMNullStatus>
class NFAAccept{
	DynamicFSM & g;
	Status & status;
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



	vec<int> next;
	vec<int> accepts;

	vec<bool> next_seen;
	vec<bool> cur_seen;

	int source;
	vec<vec<int>> & strings;
	vec<vec<bool>>  usedTransitions;
	bool checkUsed=true;
	bool hasUsed=false;

	vec<vec<int>> states_to_track_positive;//for each string, which states that might be accepted that we should track (for acceptance)
	vec<vec<int>> states_to_track_negative;//for each string, which states that might be accepted that we should track (for rejection)

	vec<vec<bool>> states_were_accepting;
	//vec<vec<bool>> states_were_rejecting;
	vec<int> n_trackingString;
	bool hasAcceptanceStates=false;

	int n_track_positive=0;
	int n_track_negative=0;
public:
	NFAAccept(DynamicFSM & f,int source, vec<vec<int>> & strings,Status & status=fsmNullStatus, bool trackUsedTransitions=false):g(f),status(status),source(source),strings(strings),checkUsed(trackUsedTransitions){

		buildStringTrackers();
	}

private:
	void buildStringTrackers(){
		states_to_track_positive.growTo(strings.size());
		states_to_track_negative.growTo(strings.size());
		n_trackingString.growTo(strings.size(),g.states()*2);
		states_were_accepting.growTo(strings.size());
		//states_were_rejecting.growTo(strings.size());

		for(int i = 0;i<strings.size();i++){
			states_to_track_positive[i].growTo(g.states(),true);
			states_to_track_negative[i].growTo(g.states(),true);

			states_were_accepting[i].growTo(g.states());
			//states_were_rejecting[i].growTo(f.states());
		}

		n_track_positive = strings.size();
		n_track_negative=strings.size();
	}
public:
	bool isTrackedState(int str, int state, bool positive){
		if(positive)
			return states_to_track_positive[str][state];
		else
			return states_to_track_negative[str][state];
	}

	void setTrackStringAcceptance(int str, int state,bool trackPositiveAcceptance, bool trackNegativeAcceptance){
		buildStringTrackers();
		if(states_to_track_positive[str][state]!=trackPositiveAcceptance){
			if(trackPositiveAcceptance){
				n_track_positive++;
				n_trackingString[str]++;
			}else{
				n_track_positive--;
				n_trackingString[str]--;
			}
			states_to_track_positive[str][state]=trackPositiveAcceptance;
		}
		if(states_to_track_negative[str][state]!=trackNegativeAcceptance){
			if(trackNegativeAcceptance){
				n_track_negative++;
				n_trackingString[str]++;
			}else{
				n_track_negative--;
				n_trackingString[str]--;
			}
			states_to_track_negative[str][state]=trackNegativeAcceptance;
		}
	}

private:

	void markUsed(int edgeID, int labelIn, int labelOut){
		assert(labelOut==0);
		if(checkUsed){
			hasUsed=true;
			usedTransitions[edgeID][labelIn]=true;
		}
	}

	void find_accepts(int str){

		for(int s:accepts){
			assert(cur_seen);
			cur_seen[s]=false;
		}
		accepts.clear();
		assert(next.size()==0);
		cur_seen[source]=true;
		accepts.push(source);

		vec<int> & string = strings[str];

		//initial emove pass:
		if(g.emovesEnabled()){
			for(int i = 0;i<accepts.size();i++){
				int s = accepts[i];
				for(int j = 0;j<g.nIncident(s);j++){
					//now check if the label is active
					int edgeID= g.incident(s,j).id;
					int to = g.incident(s,j).node;
					if(!cur_seen[to] && g.transitionEnabled(edgeID,0,0)){
						cur_seen[to]=true;
						accepts.push(to);
						markUsed(edgeID,0,0);
					}

				}
			}
		}

		for(int l:string)
		{
			assert(l>0);
			for(int i = 0;i<accepts.size();i++){
				int s = accepts[i];
				for(int j = 0;j<g.nIncident(s);j++){
					//now check if the label is active
					int edgeID= g.incident(s,j).id;
					int to = g.incident(s,j).node;
					if(!cur_seen[to] && g.transitionEnabled(edgeID,0,0)){
						cur_seen[to]=true;
						accepts.push(to);
						markUsed(edgeID,0,0);
						//status.reaches(str,to,edgeID,0);
					}

					if (!next_seen[to] && g.transitionEnabled(edgeID,l,0)){
						//status.reaches(str,to,edgeID,l);
						next_seen[to]=true;
						next.push(to);
						markUsed(edgeID,l,0);
					}
				}
			}

			next.swap(accepts);
			next_seen.swap(cur_seen);

			for(int s:next){
				assert(next_seen[s]);
				next_seen[s]=false;
			}
			next.clear();
		}

		//final emove pass:
		if(g.emovesEnabled()){
			for(int i = 0;i<accepts.size();i++){
				int s = accepts[i];
				for(int j = 0;j<g.nIncident(s);j++){
					//now check if the label is active
					int edgeID= g.incident(s,j).id;
					int to = g.incident(s,j).node;
					if(!cur_seen[to] && g.transitionEnabled(edgeID,0,0)){
						cur_seen[to]=true;
						accepts.push(to);
						markUsed(edgeID,0,0);
					}

				}
			}
		}


	}

	bool path_rec(int s, int dest,int string,int str_pos,int emove_count, vec<NFATransition> & path){
		if(str_pos==strings[string].size() && (s==dest || dest<0) ){
			return true;
		}
		if (emove_count>=g.states()){
			return false;//this is not a great way to solve the problem of avoiding infinite e-move cycles...
		}



		for(int j = 0;j<g.nIncident(s);j++){
			//now check if the label is active
			int edgeID= g.incident(s,j).id;
			int to = g.incident(s,j).node;
			if( g.transitionEnabled(edgeID,0,0)){
				path.push({edgeID,0,0});
				if(path_rec(to,dest,string,str_pos,emove_count+1,path)){//str_pos is NOT incremented!
					return true;
				}else{
					path.pop();
				}
			}
			if(str_pos< strings[string].size()){
				int l = strings[string][str_pos];
				if (g.transitionEnabled(edgeID,l,0)){
					path.push({edgeID,l,0});
					if(path_rec(to,dest,string,str_pos+1,0,path)){//str_pos is incremented
						return true;
					}else{
						path.pop();
					}
				}
			}
		}
		return false;
	}

public:

	void update(){

		if (last_modification > 0 && g.modifications == last_modification) {
			stats_skipped_updates++;
			return;
		}
		//This shouldn't normally happen
		if(n_track_positive==0 && n_track_negative==0){
			stats_skipped_updates++;
			return;
		}

		static int iteration = 0;
		int local_it = ++iteration;
		stats_full_updates++;

		if (last_deletion == g.deletions) {
			stats_num_skipable_deletions++;
		}



		bool edgesAdded=false;
		bool edgesRemoved=false;
		bool usedTransitionRemoved=!hasUsed;

		if (last_modification <= 0 || g.changed() || last_history_clear != g.historyclears) {
			next_seen.clear();
			next_seen.growTo(g.states());
			cur_seen.clear();
			cur_seen.growTo(g.states());
			hasUsed=false;
			if(checkUsed){
				usedTransitions.growTo(g.nEdgeIDs() );
				for(int i =0;i<g.nEdgeIDs();i++){
					usedTransitions[i].growTo(g.inAlphabet()+1);

				}
			}
			hasAcceptanceStates=false;
			for(int i = 0;i<strings.size();i++){
				for(int j = 0;j<states_were_accepting[i].size();j++)
					states_were_accepting[i][j]=false;
				//for(int j = 0;j<states_were_rejecting[i].size();j++)
				//	states_were_rejecting[i][j]=false;
			}
			edgesAdded=true;
			edgesRemoved=true;
			usedTransitionRemoved=true;
		}

		//it should be safe to skip updates, under the following conditions
		//1) if none of the removed transitions were traversed by the (relevant) accepting traces.
		//2) if all of the changes are added edges, AND all of the properties we are checking are asserted positive
		//3) if all of the changes are removed edges, AND all of the properties we are checking are asserted false.
		//4) If no transitions removed edges, AND all strings that were previously accepted
		//5) if no transitions added edges, AND all strings were previously rejected

		//These checks can be done per-string, or as a whole.
		for (int i = history_qhead; i < g.history.size(); i++) {
			DynamicFSM::EdgeChange & c = g.history[i];
			if(c.addition){
				edgesAdded=true;
			}else{
				edgesRemoved=true;
				if(hasUsed){
					//if we are keeping track of which transitions were used, we can apply
					usedTransitionRemoved|=usedTransitions[c.id][c.input];
				}
			}

			if(usedTransitionRemoved && edgesAdded && edgesRemoved)
				break;
		}

/*
 	 //not sure if these are safe, and they might be superseeded by the (slightly less efficient) check below, anyhow...
		if(n_track_positive==0 && !edgesRemoved){
			//we can skip this update, because no edges were removed, and we don't care if any strings were made positive
			stats_skipped_updates++;
			return;
		}
		if(n_track_negative==0 && !edgesAdded){
			//we can skip this update, because no edges were added, and we don't care if any strings were made negative
			stats_skipped_updates++;
			return;
		}*/






		if(hasUsed){
			//clear used
			hasUsed=false;
		/*	for(int i = 0;i<usedTransitions.size();i++){
				for(int l = 0;l<usedTransitions[i].size()){
					usedTransitions[i][l]=false;
				}
			}*/
		}

		//for(int i = 0;i<g.states();i++)
		for (int str = 0;str<strings.size();str++){
			if(n_trackingString[str]==0)
				continue;
			bool all_satisfied=false;
			if(hasAcceptanceStates){
				all_satisfied=true;
				if(edgesAdded){
					for(int i = 0;i<g.states();i++){
						//If we are tracking the acceptance of this state, and it didn't previously accept this string, and we enabled an edge
						//then we need to recompute acceptance
						if(states_to_track_positive[str][i] && !states_were_accepting[str][i]){
							all_satisfied=false;
							break;
						}

						//If we are tracking the rejection of this state, and it previously accepted this string, and we removed an edge
						//then we need to recompute acceptance
						if(states_to_track_negative[str][i] && !states_were_accepting[str][i]){
							all_satisfied=false;
							break;
						}
					}
				}
				if(edgesRemoved && all_satisfied){
					for(int i = 0;i<g.states();i++){

						if(states_to_track_positive[str][i] && states_were_accepting[str][i]){
							all_satisfied=false;
							break;
						}

						if(states_to_track_negative[str][i] && states_were_accepting[str][i]){
							all_satisfied=false;
							break;
						}
					}
				}

			}
			if(!all_satisfied){

				find_accepts(str);
				for(int i = 0;i<g.states();i++){
					states_were_accepting[str][i]=false;
				}
				for(int s:accepts){
					status.accepts(str,s,-1,-1,true);
					states_were_accepting[str][s]=true;
				}

				//improve this:
				for(int s = 0;s<g.states();s++){
					if (!cur_seen[s]){
						status.accepts(str,s,-1,-1,false);
					}
				}
			}else{
				for(int i = 0;i<g.states();i++){
					status.accepts(str,i,-1,-1,states_were_accepting[str][i]);
				}
			}
		}

		hasAcceptanceStates=true;

		last_modification = g.modifications;
		last_deletion = g.deletions;
		last_addition = g.additions;

		history_qhead = g.history.size();
		last_history_clear = g.historyclears;
	}




public:
	void run(int str){
		next_seen.growTo(g.states());
		cur_seen.growTo(g.states());
		find_accepts(str);
	}
	//If state is -1, then this is true if any state accepts the string.
	bool accepting( int state){
		if(state<0){
			return accepts.size();
		}
		return cur_seen[state];
	}

	//inefficient!
	//If state is -1, then this is true if any state accepts the string.
	bool acceptsString(int string, int state){

		run(string);
		return accepting(state);
	}
	bool getPath(int string, int state, vec<NFATransition> & path){
		return path_rec(source,state,string,0,0,path);
	}



};


#endif /* NFAREACH_H_ */
