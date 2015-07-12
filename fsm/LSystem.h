/*
 * DynamicFSM.h
 *
 *  Created on: Dec 15, 2014
 *      Author: sam
 */

#ifndef DYNAMICLSYSTEM_H_
#define DYNAMICLSYSTEM_H_

#include <vector>
#include "mtl/Vec.h"
#include "mtl/Bitset.h"
#include <algorithm>
#include <cassert>

#include "dgl/DynamicGraph.h"
using namespace dgl;
namespace Monosat {


//This models a context-free, non-deterministic lsystem
class LSystem{

	bool is_changed = true;
public:
	vec<bool> enabled_rules;
	int characters=0;
	vec<bool> terminal_character;
	struct Rule{
		int predessor=-1;
		vec<int> production;
	};

	vec<Rule> rules;
	vec<vec<int> > all_rules;
	bool producing=true;
	bool strictlyProducing=true;
	int nTerminalCharacters=0;
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
	LSystem() {

	}

	~LSystem() {

	}

	bool isProducing(){
		return producing;
	}

	//True if every character is either terminal, or only has rules that produce 2 or more characters.
	bool isStrictlyProducing(){
		return strictlyProducing;
	}

	void addCharacter(int chars=1){
		characters +=chars;
		nTerminalCharacters+=chars;
		enabled_rules.growTo(characters);
		all_rules.growTo(characters);
		terminal_character.growTo(characters,true);

	}
	bool hasRule(int rule){
		return rule>=0 && rule<rules.size() && rules[rule].predessor>=0;
	}

	int getPredessor(int rule){
		assert(hasRule(rule));
		return rules[rule].predessor;
	}

	vec<int> & getRule(int rule){
		assert(hasRule(rule));
		return rules[rule].production;
	}

	vec<int> & getRules(int character){
		return all_rules[character];
	}
	int nCharacters()const{
		return characters;
	}
	int nRules()const{
		return rules.size();
	}
	//Add a production rule to a character
	int addRule(int predecessor,vec<int> & production,  bool defaultEnabled=false){
		int ruleID = rules.size();
		assert(predecessor<characters);
		rules.growTo(ruleID+1);
		assert(rules[ruleID].predessor==-1);
		rules[ruleID].predessor=predecessor;
		production.copyTo(rules[ruleID].production);

		if(production.size()!= 1 || production[0]!=predecessor){
			if(terminal_character[predecessor]){
				nTerminalCharacters--;
				terminal_character[predecessor]=false;
			}
		}
		if(production.size()==0){
			producing=false;
			strictlyProducing=false;
		}
		if(production.size()== 1 && production[0]!=predecessor){
			strictlyProducing=false;
		}

		all_rules[predecessor].push(ruleID);

		enabled_rules.growTo(ruleID+1);
		enabled_rules[ruleID]=defaultEnabled;
		return ruleID;
	}

	bool isTerminal(int character){
		return terminal_character[character];
	}

	unsigned int inAlphabet()const{
		return in_alphabet;
	}

	void addInCharacter(){
		in_alphabet++;
	}

	bool ruleEnabled(int ruleID)const{
		return enabled_rules[ruleID];
	}

	void enableRule(int ruleID) {

		if (!enabled_rules[ruleID]) {
			enabled_rules[ruleID]=true;
			//edge_status.setStatus(id,true);
			modifications++;
			additions = modifications;
			history.push_back( { true, ruleID, modifications, additions });
		}
	}

	void disableRule(int ruleID) {
		if (enabled_rules[ruleID]) {
			enabled_rules[ruleID]=false;
			modifications++;
			history.push_back( { false, ruleID, modifications, deletions });
			deletions = modifications;
		}
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
										std::max(1000L, historyClearInterval * nRules()) : historyClearInterval)))) {//){
			history.clear();
			historyclears++;

		}

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

	}

	void draw(int source=-1, int dest=-1){


	}

};

}
;



#endif /* DYNAMICLSYSTEM_H_ */
