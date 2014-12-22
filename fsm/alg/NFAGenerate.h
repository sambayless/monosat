/*
 * NFAReach.h
 *
 *  Created on: Dec 16, 2014
 *      Author: sam
 */

#ifndef NFA_GENERATE_H_
#define NFA_GENERATE_H_
#include "../DynamicFSM.h"
#include "mtl/Bitset.h"
#include "mtl/Vec.h"
#include "NFATypes.h"

using namespace Monosat;



template<class Status=FSMNullStatus>
class NFAGenerate{
	DynamicFSM & g;
	Status & status;
	int last_modification=-1;

	int last_addition=-1;
	int last_deletion=-1;
	int history_qhead=0;
	int last_history_clear=0;

	int stats_full_updates=0;
	int stats_fast_updates=0;
	int stats_fast_failed_updates=0;
	int stats_skip_deletes=0;
	int stats_skipped_updates=0;
	int stats_num_skipable_deletions=0;
	double mod_percentage=0;

	double stats_full_update_time=0;
	double stats_fast_update_time=0;


	vec<int> next;
	vec<int> accepts;
	vec<bool> next_seen;
	vec<bool> cur_seen;

	int source;
	vec<vec<int>> & strings;



public:
	NFAGenerate(DynamicFSM & f,int source, vec<vec<int>> & strings,Status & status=fsmNullStatus):g(f),status(status),source(source),strings(strings){

	}

private:
	struct UsedTransition{
		int edge=-1;
		int label=-1;
	};

	vec<UsedTransition> used_transitions;

	bool unique_path_rec(int s,int string,int str_pos,int emove_count, vec<NFATransition> & path){
		if(str_pos==strings[string].size()){
			return true;
		}
		if (emove_count>=g.states()){
			return false;//this is not a great way to solve the problem of avoiding infinite e-move cycles...
		}

		int l = strings[string][str_pos];
		for(int j = 0;j<g.nIncident(s);j++){
			//now check if the label is active
			int edgeID= g.incident(s,j).id;
			int to = g.incident(s,j).node;
			if( g.transitionEnabled(edgeID,0,0)){
				bool set_transition=false;
				if(used_transitions[s].edge==-1){
					used_transitions[s].edge=edgeID;
					used_transitions[s].label=0;
					set_transition=true;
				//	used_transitions[s].depth = path.size();
				}

				if(used_transitions[s].label==0 && used_transitions[s].edge==edgeID){

					path.push({edgeID,0,0});
					if(unique_path_rec(to,string,str_pos,emove_count+1,path)){//str_pos is NOT incremented!
						if(set_transition){
							used_transitions[s].edge=-1;
							used_transitions[s].label=-1;
						}
						return true;
					}else{
						if(set_transition){
							used_transitions[s].edge=-1;
							used_transitions[s].label=-1;
						}
						path.pop();
					}
				}
			}
			if(str_pos< strings[string].size()){
				if (g.transitionEnabled(edgeID,0,l)){
					bool set_transition=false;
					if(used_transitions[s].edge==-1){
						used_transitions[s].edge=edgeID;
						used_transitions[s].label=l;
						set_transition=true;
						//used_transitions[s].depth = path.size();
					}
					if(used_transitions[s].label==l && used_transitions[s].edge==edgeID){
						path.push({edgeID,0,l});
						if(unique_path_rec(to,string,str_pos+1,0,path)){//str_pos is incremented
							if(set_transition){
								used_transitions[s].edge=-1;
								used_transitions[s].label=-1;
							}
							return true;
						}else{
							if(set_transition){
								used_transitions[s].edge=-1;
								used_transitions[s].label=-1;
							}
							path.pop();
						}
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
		static int iteration = 0;
		int local_it = ++iteration;
		stats_full_updates++;

		if (last_deletion == g.deletions) {
			stats_num_skipable_deletions++;
		}

		if (last_modification <= 0 || g.changed() || last_history_clear != g.historyclears) {
			next_seen.clear();
			next_seen.growTo(g.states());
			cur_seen.clear();
			cur_seen.growTo(g.states());
		}
		for(int i = 0;i<strings.size();i++){
			if(generatesString(i)){
				status.generates(i,true);
			}else{
				status.generates(i,false);
			}
		}
/*
		//first, apply e-moves

		//for(int i = 0;i<g.states();i++)
		for (int str = 0;str<strings.size();str++){

			for(int s:accepts){
				status.generates(str,s,-1,,true);
			}
			//improve this:
			for(int s = 0;s<g.states();s++){
				if (!cur_seen[s]){
					status.generates(str,s,-1,,false);
				}
			}
		}*/



		last_modification = g.modifications;
		last_deletion = g.deletions;
		last_addition = g.additions;

		history_qhead = g.history.size();
		last_history_clear = g.historyclears;
	}




public:



	//inefficient!
	bool generatesString(int string){
		static vec<NFATransition> ignore;
		ignore.clear();
		return getPath(string, ignore);
	}

	bool getPath(int string, vec<NFATransition> & path){
		used_transitions.growTo(g.states());
		return unique_path_rec(source,string,0,0,path);
	}



};


#endif /* NFAREACH_H_ */
