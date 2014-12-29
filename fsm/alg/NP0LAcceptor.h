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
#include "../DynamicFSM.h"
//#include <fst/fstlib.h>
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
	int atom = 1;
	vec<vec<int>> & strings;
	vec<vec<int>>  stringset;
	vec<vec<Bitset>> suffixTables;
	vec<vec<int>> toChecks;

	DynamicFSM acceptor;
	struct RuleTransition{

		int edgeID=-1;
		int inChar;
		int outChar;
	};
	vec<RuleTransition> ruleMap;
public:
	NP0LAccept(LSystem & f, vec<vec<int>> & strings):g(f),strings(strings){
		assert(f.strictlyProducing);


		//bool buildSuffixTable(int startState, int finalState, vec<int> & string, vec<Bitset> & table){
		for(int i = 0;i<=g.nCharacters();i++){
			acceptor.addInCharacter();
			acceptor.addOutCharacter();

		}
		int start = acceptor.addState();
		for(int c = 0;c<f.nCharacters();c++){
			int  cState = acceptor.addState();
			//Add an arc that takes epsilon and outputs this character.
			acceptor.addTransition(start,cState,-1,0,c+1);
			//In the future: combine rules into a prefix tree first.

			for(int rID : f.getRules(c)){
				vec<int> & rule = f.getRule(rID);
				int from = cState;
				int firstRule = -1;
				for(int i = 0;i<rule.size();i++){
					int o = rule[i];
					int next;
					if(i<rule.size()-1){
						next = acceptor.addState();
					}else{
						next=start;
					}

					//Add an arc that takes this rule character, and outputs epsilon
					int edgeID = acceptor.addTransition(from,next,-1,o+1,0);
					if(firstRule<0){
						firstRule=edgeID;
						ruleMap.growTo(rID+1);
						ruleMap[rID].edgeID = edgeID;
						ruleMap[rID].inChar=o+1;
						ruleMap[rID].outChar=0;
					}
					from = next;
				}
			}
		}


	}


private:
	void setRuleEnabled(int ruleID, bool enable){
		assert(ruleMap[ruleID].edgeID>=0);
		if(enable){
			acceptor.enableTransition(ruleMap[ruleID].edgeID,ruleMap[ruleID].inChar,ruleMap[ruleID].outChar);
		}else{
			acceptor.disableTransition(ruleMap[ruleID].edgeID,ruleMap[ruleID].inChar,ruleMap[ruleID].outChar);
		}
	}

	bool ruleEnabled(int ruleID){
		assert(ruleMap[ruleID].edgeID>=0);
		return acceptor.transitionEnabled(ruleMap[ruleID].edgeID,ruleMap[ruleID].inChar,ruleMap[ruleID].outChar);
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
			for(int i = 0;i<g.nRules();i++){
				if(g.hasRule(i)){
					setRuleEnabled(i,g.ruleEnabled(i));
				}
			}
		}else{
			for (int i = history_qhead; i < g.history.size(); i++) {
				int edgeid = g.history[i].id;

				if (g.history[i].addition && g.ruleEnabled(edgeid) && !ruleEnabled(edgeid)) {

					setRuleEnabled(edgeid,true);

				} else if (!g.history[i].addition && !g.ruleEnabled(edgeid) && ruleEnabled(edgeid)) {

					setRuleEnabled(edgeid,false);

				}
			}
		}



		last_modification = g.modifications;
		last_deletion = g.deletions;
		last_addition = g.additions;

		history_qhead = g.history.size();
		last_history_clear = g.historyclears;

	}

private:

	bool path_rec(int s, int dest,vec<int> & string,int str_pos,int emove_count,int depth,vec<Bitset> & suffixTable, vec<int> & path){
		if(str_pos==string.size() && (s==dest || dest<0) ){
			//string accepted by lsystem.
			if(path.size()==0){
				return false;
			}else if(path.size()==1 && path[0]==atom){
				return true;
			}else
				return accepts_rec(-1,depth+1);

		}
		if (emove_count>=acceptor.states()){
			return false;//this is not a great way to solve the problem of avoiding infinite e-move cycles...
		}
		//if(!suffixTable[str_pos][s])
		//	return false;
		int l = string[str_pos];
		for(int j = 0;j<acceptor.nIncident(s);j++){
			//now check if the label is active
			int edgeID= acceptor.incident(s,j).id;
			int to = acceptor.incident(s,j).node;
			for(int o = 0;o<acceptor.outAlphabet();o++){
				if(acceptor.transitionEnabled(edgeID,0,o)){
					//assert(suffixTable[str_pos][to]);
					if(o>0)
						path.push(o);

					if(path_rec(to,dest,string,str_pos,emove_count+1,depth,suffixTable,path)){//str_pos is NOT incremented!
						return true;
					}else if (o>0){
						path.pop();
					}
				}
			}

			if(str_pos< string.size()){// && suffixTable[str_pos+1][to]){
				for(int o = 0;o<acceptor.outAlphabet();o++){
					if (acceptor.transitionEnabled(edgeID,l,o)){
						if(o>0)
							path.push(o);

						if(path_rec(to,dest,string,str_pos+1,0,depth,suffixTable,path)){//str_pos is incremented
							return true;
						}else if (o>0){
							path.pop();
						}
					}
				}
			}
		}
		return false;
	}

	bool accepts_rec(int str,int depth){
		if(depth>5)
			return false;
		vec<int> & string = depth==0 ? strings[str] :stringset[depth];
		stringset.growTo(depth+2);
		suffixTables.growTo(depth+1);
		vec<Bitset> & suffixTable = suffixTables[depth];
		//build suffix table of states that can reach the final state from the nth suffix of the string
		acceptor.buildSuffixTable(0,0,string,suffixTable);

		//now find all paths through the nfa using a dfs, but filtering the search using the suffix table so that all paths explored are valid paths.
		toChecks.growTo(depth+1);
		vec<int> & toCheck = toChecks[depth];
		int str_pos = 0;
		stringset[depth+1].clear();
		return path_rec(0,0,string,0,0,depth,suffixTable,stringset[depth+1]);

	}


public:


	//inefficient!
	//If state is -1, then this is true if any state accepts the string.
	bool acceptsString(int string){
		update();

		return accepts_rec(string,0);
	}



};


#endif
