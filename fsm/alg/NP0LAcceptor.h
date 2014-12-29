/*
 * NFAReach.h
 *
 *  Created on: Dec 16, 2014
 *      Author: sam
 */

#ifndef NPOL_ACCEPT_H_
#define NPOL_ACCEPT_H_
#include "../LSystem.h"
#include "mtl/Bitset.h"
#include "mtl/Vec.h"
#include "NFATypes.h"
#include <fst/fstlib.h>
using namespace Monosat;



class NP0LAccept{
	LSystem & g;

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
	NP0LAccept(LSystem & f, vec<vec<int>> & strings):g(f),strings(strings){
		assert(f.strictlyProducing);
		fst::StdVectorFst t;
		const auto one = fst::StdArc::Weight::One();
		auto  start = t.AddState();
		t.SetStart(start);
		t.SetFinal(start,one);

		for(int c = 0;c<f.nCharacters();c++){
			int  cState = t.AddState();
			//Add an arc that takes epsilon and outputs this character.
			t.AddArc(start, fst::StdArc(0, c+1, one, cState));


			//In the future: combine rules into a prefix tree first.
			for(int rID : f.getRules(c)){
				vec<int> & rule = f.getRule(rID);
				int from = cState;
				for(int i = 0;i<rule.size();i++){
					int o = rule[i];
					int next;
					if(i<rule.size()-1){
						next = t.AddState();
					}else{
						next=start;
					}

					//Add an arc that takes this rule character, and outputs epsilon
					t.AddArc(from, fst::StdArc(o+1,0, one, next));
					from = next;
				}
			}
		}

		t.Write("/home/sam/uncombined.fst");

		fst::StdVectorFst * c = new fst::StdVectorFst();
		int depth = 5;
		printf("composing %d...\n",0);
		fst::Compose(t,t,c);
		for(int i = 1;i<depth;i++){
			printf("composing %d...\n",i);
			fst::StdVectorFst * c_out = new fst::StdVectorFst();
			fst::Compose(*c,t,c_out);
			delete(c);
			c=c_out;
			long narcs = fst::CountArcs(*c);
			printf("done composing %d, %d states, %d transitions...\n",i,c->NumStates(),narcs);
		}
		c->Write("/home/sam/combined.fst");
		delete(c);
	}




public:

	void update(){


	}




public:
	void run(int str){

	}
	//If state is -1, then this is true if any state accepts the string.
	bool accepting( int str){
		return false;
	}

	//inefficient!
	//If state is -1, then this is true if any state accepts the string.
	bool acceptsString(int string){

		run(string);
		return false;
	}
	bool getPath(int string, vec<NFATransition> & path){
		return false;
	}



};


#endif
