/*
 * NFAReach.h
 *
 *  Created on: Dec 16, 2014
 *      Author: sam
 */

#ifndef NFAREACH_H_
#define NFAREACH_H_
#include "../DynamicFSM.h"
#include "mtl/Bitset.h"
#include "mtl/Vec.h"
#include "NFATypes.h"
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
	NFAAccept(DynamicFSM & f,int source, vec<vec<int>> & strings,Status & status=fsmNullStatus):g(f),status(status),source(source),strings(strings){

	}

private:
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
						//status.reaches(str,to,edgeID,0);
					}

					if (!next_seen[to] && g.transitionEnabled(edgeID,l,0)){
						//status.reaches(str,to,edgeID,l);
						next_seen[to]=true;
						next.push(to);
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


		//first, apply e-moves

		//for(int i = 0;i<g.states();i++)
		for (int str = 0;str<strings.size();str++){
			find_accepts(str);
			for(int s:accepts){
				status.accepts(str,s,-1,-1,true);
			}
			//improve this:
			for(int s = 0;s<g.states();s++){
				if (!cur_seen[s]){
					status.accepts(str,s,-1,-1,false);
				}
			}
		}



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
