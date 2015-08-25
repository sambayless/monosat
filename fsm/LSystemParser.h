/****************************************************************************************[Solver.h]
 The MIT License (MIT)

 Copyright (c) 2014, Sam Bayless

 Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
 associated documentation files (the "Software"), to deal in the Software without restriction,
 including without limitation the rights to use, copy, modify, merge, publish, distribute,
 sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all copies or
 substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
 NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
 OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 **************************************************************************************************/

#ifndef LSYSTEM_PARSER_H_
#define LSYSTEM_PARSER_H_

#include <stdio.h>

#include "utils/ParseUtils.h"
#include "core/SolverTypes.h"
#include "fsm/LSystemTheory.h"
#include <algorithm>
#include "core/Config.h"

#include "core/Dimacs.h"

#include <set>
#include <string>
#include <sstream>
namespace Monosat {


//=================================================================================================
// GRAPH Parser:
template<class B, class Solver>
class LSystemParser: public Parser<B, Solver> {

	vec<LSystemSolver*> lsystems;
	vec<int> alphabets;

	struct Transition{
		int lsystem;
		int from;
		Var edgeVar;
		vec<int> production;
	};

	vec<vec<Transition> > rules;
	vec<bool> created_strings;
	vec<vec<int>> strings;
	vec<int> stringLabels;

	struct Accepts{
		int lsystem;
		int atom;
		int from;
		int to;
		int strID;
		Var reachVar;

	};
	vec<vec<Accepts> > produces;


	vec<Lit> lits;
	int count = 0;

	void readLSystem(B& in,  Solver& S) {
		if (opt_ignore_theories) {
			skipLine(in);
			return;
		}

		int lsystemID = parseInt(in);  //id of the lsystem
		int chars = parseInt(in);
		int nrules= parseInt(in);
		if (lsystemID < 0 ) {
			parse_errorf("lsystem id must be >=0, was %d\n", lsystemID);
		}


		lsystems.growTo(lsystemID + 1);
		if(lsystems[lsystemID]){
			parse_errorf("lsystem id %d declared twice!\n", lsystemID);
		}
		lsystems[lsystemID]= new LSystemSolver(&S);
		S.addTheory(lsystems[lsystemID]);
		rules.growTo(lsystemID+1);
		produces.growTo(lsystemID+1);

		alphabets.growTo(lsystemID+1);

	}
	
	void readString(B& in, Solver & S){
		if (opt_ignore_theories) {
			skipLine(in);
			return;
		}

		int strID = parseInt(in);
		strings.growTo(strID+1);
		created_strings.growTo(strID+1);
		stringLabels.growTo(strID+1);
		if(strID<0 || created_strings[strID]){
			parse_errorf("Bad string id %d\n", strID);
		}

		created_strings[strID]=true;

		//allow zero-length strings
		if (isEof(in) || *in == '\n')
			return;
		while(int i = parseInt(in)){
			if (i<=0){
				parse_errorf("lsystem strings must contain only positive (non-zero) integers, found %d\n", i);
			}
			strings[strID].push(i);
			stringLabels[strID]= std::max(stringLabels[strID],i+1);
			skipWhitespace(in);
			if (isEof(in) || *in == '\n')
				break;

		}

	}

	void readRule(B& in, Solver& S) {
		if (opt_ignore_theories) {
			skipLine(in);
			return;
		}
		
		++in;
		
		int lsystemID = parseInt(in);
		int from = parseInt(in);
		int edgeVar = parseInt(in) - 1;


		if (lsystemID < 0 || lsystemID >= lsystems.size()) {
			parse_errorf("Undeclared lsystem identifier %d for edge %d\n", lsystemID, edgeVar);
		}

		if (edgeVar < 0) {
			parse_errorf("Edge variables must be >=0, was %d\n", edgeVar);
		}


		while (edgeVar >= S.nVars())
			S.newVar();

		alphabets[lsystemID]=std::max(alphabets[lsystemID],from);

		rules[lsystemID].push();
		rules[lsystemID].last().lsystem =lsystemID;
		rules[lsystemID].last().from =from;
		rules[lsystemID].last().edgeVar =edgeVar;

		//allow zero-length strings
		if (isEof(in) || *in == '\n')
			return;
		while(true){
			int i = parseInt(in);
			rules[lsystemID].last().production.push(i);
			alphabets[lsystemID]=std::max(alphabets[lsystemID],i);
			skipWhitespace(in);
			if (isEof(in) || *in == '\n')
				break;
		}
	}

	void readProduces(B& in, Solver& S) {
		if (opt_ignore_theories) {
			skipLine(in);
			return;
		}

		++in;
		
		int lsystemID = parseInt(in);
		int atom = parseInt(in);
		int strID = parseInt(in);

		int reachVar = parseInt(in) - 1;

		//now read in the string
		produces[lsystemID].push();

		produces[lsystemID].last().lsystem=lsystemID;
		produces[lsystemID].last().atom = atom;
		produces[lsystemID].last().strID = strID;
		produces[lsystemID].last().reachVar = reachVar;

		if (lsystemID < 0 || lsystemID >= lsystems.size()) {
			parse_errorf("Undeclared lsystem identifier %d for edge %d\n", lsystemID, reachVar);
		}



		if (reachVar < 0) {
			parse_errorf("Edge variables must be >=0, was %d\n", reachVar);
		}
		
		while (reachVar >= S.nVars())
			S.newVar();


	}


public:
	LSystemParser():Parser<B, Solver>("LSystem"){
		
	}

	bool parseLine(B& in, Solver& S) {

		skipWhitespace(in);
		if (*in == EOF)
			return false;

		if (match(in, "lsystem")) {
			skipWhitespace(in);
			readLSystem(in,S);
			skipWhitespace(in);

			return true;
		} else if (match(in, "production")) {
			count++;
			readRule(in, S);
			return true;
		}else if (match(in,"lstr")){
			readString(in,S);
			return true;
		}else if (match(in, "produces")) {
			readProduces(in, S);
			return true;
		}
		return false;
	}
	

	void implementConstraints(Solver & S) {
		
		for(int i = 0;i<lsystems.size();i++){
			if(lsystems[i]){

				lsystems[i]->setAlphabet(alphabets[i]+1);

				lsystems[i]->setStrings(&strings);


				for (auto &t:rules[i]){
					lsystems[i]->newRule(t.from,t.production,t.edgeVar);
				}

				for(auto & a: produces[i]){

					if (a.strID<0 || !created_strings[a.strID]){
						parse_errorf("String ID must be a non-negative integer, was %d\n", a.strID);
					}


					lsystems[i]->addProducesLit(a.atom, a.strID,a.reachVar);
				}

			}
		}




	}

	
};

//=================================================================================================
}
;

#endif /* GRAPH_PARSER_H_ */
