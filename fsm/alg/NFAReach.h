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
using namespace Monosat;
struct FSMNullStatus {
	void reaches(int string, int state,int edgeID,int label) {

	}
};
static FSMNullStatus fsmNullStatus;
template<class Status=FSMNullStatus>
class NFAReach{
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
	vec<int> reaches;
	vec<bool> next_seen;
	vec<bool> cur_seen;

	int source;
	vec<vec<int>> & strings;
public:
	NFAReach(DynamicFSM & f,int source, vec<vec<int>> & strings,Status & status=fsmNullStatus):g(f),status(status),source(source),strings(strings){

	}

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
		for(int s:reaches){
			assert(cur_seen);
			cur_seen[s]=false;
		}
		reaches.clear();
		assert(next.size()==0);

		cur_seen[source]=true;
		reaches.push(source);
		//first, apply e-moves

		//for(int i = 0;i<g.states();i++)
		for (int str = 0;str<strings.size();str++){
			vec<int> & string = strings[str];
			for(int l:string)
			{
				assert(l>0);
				for(int s:reaches){
					for(int j = 0;j<g.nIncident(s);j++){
						//now check if the label is active
						int edgeID= g.incident(s,j).id;
						int to = g.incident(s,j).node;
						if(!cur_seen[to] && g.emove(edgeID)){
							cur_seen[to]=true;
							reaches.push(to);
							status.reaches(str,to,edgeID,0);
						}

						if (!next_seen[to] && g.transitionEnabled(edgeID,l)){
							status.reaches(str,to,edgeID,l);
							next_seen[to]=true;
							next.push(to);
						}
					}
				}

				next.swap(reaches);
				next_seen.swap(cur_seen);

				for(int s:next){
					assert(next_seen[s]);
					next_seen[s]=false;
				}
				next.clear();
			}
		}


		last_modification = g.modifications;
		last_deletion = g.deletions;
		last_addition = g.additions;

		history_qhead = g.history.size();
		last_history_clear = g.historyclears;
	}

	bool reachable(int state){
		update();
		return reaches[state];
	}

};


#endif /* NFAREACH_H_ */
