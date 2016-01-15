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

#ifndef FSM_PARSER_H_
#define FSM_PARSER_H_

#include <stdio.h>

#include "utils/ParseUtils.h"
#include "core/SolverTypes.h"
#include "fsm/FSMTheory.h"
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
class FSMParser: public Parser<B, Solver> {
	using Parser<B, Solver>::mapVar;
	vec<int> fsmIDs;

	vec<int> inAlphabets;
	vec<int> outAlphabets;
	struct Transition{
		int fsm;
		int from;
		int to;
		int input;
		int output;
		Var edgeVar;
	};
	vec<bool> hasEpsilonTransitions;
	vec<vec<Transition> > transitions;
	vec<bool> created_strings;
	vec<int> stringIDMap;

	vec<vec<int>> * strings;
	vec<int> stringLabels;
	struct Accepts{
		int fsm;
		int from;
		int to;
		int strID;
		Var reachVar;

	};
	vec<vec<Accepts> > accepts;

	struct ComposeAccepts{
		int fsmID1;
		int fsmID2;
		int from1;
		int from2;
		int to1;
		int to2;
		int strID;
		Var reachVar;

	};
	vec<ComposeAccepts>  compose_accepts;
	struct Generates{
		int fsm;
		int from;
		int strID;
		Var reachVar;
	};

	vec<vec<Generates> > generates;

	struct Transduces{
			int fsm;
			int from;
			int to;
			int strID;
			int strID2;
			Var reachVar;
		};

	vec<vec<Transduces> > transduces;


	vec<Lit> lits;
	int count = 0;

	void readFSM(B& in,  Solver& S) {
		if (opt_ignore_theories) {
			skipLine(in);
			return;
		}

		int fsmID = parseInt(in);  //id of the fsm
		int n_labels = parseInt(in);
		bool hasEpsilon= parseInt(in)>0;
		if (fsmID < 0 ) {
			parse_errorf("FSM id must be >=0, was %d\n", fsmID);
		}
		if (n_labels<0){
			parse_errorf("Number of transition labels must be >=0, was %d\n", n_labels);
		}

		fsmIDs.growTo(fsmID + 1,-1);
		if(fsmIDs[fsmID]>=0){
			parse_errorf("FSM id %d declared twice!\n", fsmID);
		}
		//fsms[fsmID]= new FSMTheorySolver(&S);
		fsmIDs[fsmID]=fsmID;
		//S.addTheory(fsms[fsmID]);
		transitions.growTo(fsmID+1);
		accepts.growTo(fsmID+1);
		generates.growTo(fsmID+1);
		transduces.growTo(fsmID+1);
		hasEpsilonTransitions.growTo(fsmID+1);
		inAlphabets.growTo(fsmID+1,0);
		outAlphabets.growTo(fsmID+1,0);
		hasEpsilonTransitions[fsmID]=false;
	}
	
	void readString(B& in, Solver & S){
		if (opt_ignore_theories) {
			skipLine(in);
			return;
		}

		int strID = parseInt(in);
		//strings->growTo(strID+1);

		stringIDMap.growTo(strID+1,-1);
		stringIDMap[strID]=strings->size();
		strID = strings->size();
		strings->push();

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
				parse_errorf("FSM strings must contain only positive (non-zero) integers, found %d\n", i);
			}
			strings->[strID].push(i);
			stringLabels[strID]= std::max(stringLabels[strID],i+1);
			skipWhitespace(in);
			if (isEof(in) || *in == '\n')
				break;

		}

	}

	void readTransition(B& in, Solver& S) {
		if (opt_ignore_theories) {
			skipLine(in);
			return;
		}
		
		++in;
		
		int fsmID = parseInt(in);
		int from = parseInt(in);
		int to = parseInt(in);
		int input = parseInt(in);
		int output = parseInt(in);
		int edgeVar = parseInt(in) - 1;
		
		if (fsmID < 0 || fsmID >= fsmIDs.size()) {
			parse_errorf("Undeclared fsm identifier %d for edge %d\n", fsmID, edgeVar);
		}
		if (input<0){
			parse_errorf("Transition inputs  must be >=0, was %d\n", input);
		}
		if (output<0){
				parse_errorf("Transition outputs  must be >=0, was %d\n", output);
			}
		if (edgeVar < 0) {
			parse_errorf("Edge variables must be >=0, was %d\n", edgeVar);
		}

		if (input==0){
			hasEpsilonTransitions[fsmID]=true;
		}

		edgeVar= mapVar(S,edgeVar);
		inAlphabets[fsmID]=std::max(inAlphabets[fsmID],input+1);
		outAlphabets[fsmID]=std::max(outAlphabets[fsmID],output+1);
		transitions[fsmID].push({fsmID,from,to,input,output,edgeVar});
	}

	void readAccepts(B& in, Solver& S) {
		if (opt_ignore_theories) {
			skipLine(in);
			return;
		}

		++in;
		
		int fsmID = parseInt(in);
		int from = parseInt(in);
		int to = parseInt(in);
		int strID = parseInt(in);
		int reachVar = parseInt(in) - 1;
		reachVar= mapVar(S,reachVar);
		if (fsmID < 0 || fsmID >= fsmIDs.size()) {
			parse_errorf("Undefined fsm %d",fsmID);
		}
		//now read in the string
		accepts[fsmID].push();

		accepts[fsmID].last().fsm=fsmID;
		accepts[fsmID].last().from=from;
		accepts[fsmID].last().to=to;
		accepts[fsmID].last().strID = strID;
		accepts[fsmID].last().reachVar = reachVar;

		if (fsmID < 0 || fsmID >= fsmIDs.size()) {
			parse_errorf("Undeclared fsm identifier %d for edge %d\n", fsmID, reachVar);
		}

		if (from<0){
				parse_errorf("Source state must be a node id (a non-negative integer), was %d\n", from);
			}
		if (to<0){
			parse_errorf("Accepting state must be a node id (a non-negative integer), was %d\n", to);
		}
		if (reachVar < 0) {
			parse_errorf("Edge variables must be >=0, was %d\n", reachVar);
		}



	}
	
	void readCompositionAccepts(B& in, Solver& S) {
		if (opt_ignore_theories) {
			skipLine(in);
			return;
		}

		++in;

		int fsmID1 = parseInt(in);
		int fsmID2 = parseInt(in);
		int from1 = parseInt(in);
		int to1 = parseInt(in);
		int from2 = parseInt(in);
		int to2 = parseInt(in);
		int strID = parseInt(in);
		int reachVar = parseInt(in) - 1;
		reachVar= mapVar(S,reachVar);
		//now read in the string
		compose_accepts.push();

		compose_accepts.last().fsmID1=fsmID1;
		compose_accepts.last().fsmID2=fsmID2;
		compose_accepts.last().from1=from1;
		compose_accepts.last().from2=from2;
		compose_accepts.last().to1=to1;
		compose_accepts.last().to2=to2;
		compose_accepts.last().strID = strID;
		compose_accepts.last().reachVar = reachVar;


		if (from1<0){
				parse_errorf("Source state must be a node id (a non-negative integer), was %d\n", from1);
			}
		if (to1<0){
			parse_errorf("Accepting state must be a node id (a non-negative integer), was %d\n", to1);
		}

		if (fsmID2 < 0 || fsmID2 >= fsmIDs.size()) {
			parse_errorf("Undeclared fsm identifier %d for edge %d\n", fsmID2, reachVar);
		}

		if (from2<0){
				parse_errorf("Source state must be a node id (a non-negative integer), was %d\n", from2);
			}
		if (to2<0){
			parse_errorf("Accepting state must be a node id (a non-negative integer), was %d\n", to2);
		}

		if (reachVar < 0) {
			parse_errorf("Edge variables must be >=0, was %d\n", reachVar);
		}



	}

	void readGenerates(B& in, Solver& S) {
		if (opt_ignore_theories) {
			skipLine(in);
			return;
		}

		++in;

		int fsmID = parseInt(in);
		int from = parseInt(in);

		int strID = parseInt(in);
		int reachVar = parseInt(in) - 1;
		reachVar= mapVar(S,reachVar);
		//now read in the string
		generates[fsmID].push();

		generates[fsmID].last().fsm=fsmID;
		generates[fsmID].last().from=from;

		generates[fsmID].last().strID = strID;
		generates[fsmID].last().reachVar = reachVar;

		if (fsmID < 0 || fsmID >= fsmIDs.size()) {
			parse_errorf("Undeclared fsm identifier %d for edge %d\n", fsmID, reachVar);
		}

		if (from<0){
				parse_errorf("Source state must be a node id (a non-negative integer), was %d\n", from);
			}

		if (reachVar < 0) {
			parse_errorf("Edge variables must be >=0, was %d\n", reachVar);
		}




	}

	void readTransduces(B& in, Solver& S) {
		if (opt_ignore_theories) {
			skipLine(in);
			return;
		}

		++in;

		int fsmID = parseInt(in);
		int from = parseInt(in);
		int to = parseInt(in);
		int strID = parseInt(in);
		int strID2 = parseInt(in);
		int reachVar = parseInt(in) - 1;
		reachVar= mapVar(S,reachVar);
		//now read in the string
		transduces[fsmID].push();

		transduces[fsmID].last().fsm=fsmID;
		transduces[fsmID].last().from=from;
		transduces[fsmID].last().to=to;

		transduces[fsmID].last().strID = strID;
		transduces[fsmID].last().strID2 = strID2;
		transduces[fsmID].last().reachVar = reachVar;

		if (fsmID < 0 || fsmID >= fsmIDs.size()) {
			parse_errorf("Undeclared fsm identifier %d for edge %d\n", fsmID, reachVar);
		}

		if (from<0){
				parse_errorf("Source state must be a node id (a non-negative integer), was %d\n", from);
			}
		if (to<0){
				parse_errorf("Accepting state must be a node id (a non-negative integer), was %d\n", to);
			}
		if (reachVar < 0) {
			parse_errorf("Edge variables must be >=0, was %d\n", reachVar);
		}



	}

public:
	FSMParser():Parser<B, Solver>("Finite State Machine"){
		strings = new vec<vec<int>>();
	}

	bool parseLine(B& in, Solver& S) {

		skipWhitespace(in);
		if (*in == EOF)
			return false;

		if (match(in, "fsm")) {
			skipWhitespace(in);
			readFSM(in,S);
			skipWhitespace(in);

			return true;
		} else if (match(in, "transition")) {
			count++;
			readTransition(in, S);
			return true;
		}else if (match(in,"str")){
			readString(in,S);
			return true;
		}else if (match(in,"accepts_composition")){//previously, was 'composition_accepts'; changed to avoid conflicting with comments
			readCompositionAccepts(in,S);
			return true;
		}else if (match(in, "accepts")) {
			readAccepts(in, S);
			return true;
		}else if (match(in,"generates")){
			readGenerates(in,S);
			return true;
		}else if (match(in,"transduces")){
			readTransduces(in,S);
			return true;
		}
		return false;
	}
	

	void implementConstraints(Solver & S) {
		FSMTheorySolver * theory=nullptr;

		for(int i = 0;i<fsmIDs.size();i++){
			int fsmID = fsmIDs[i];
			if(fsmID<0)
				continue;

			if(!theory){
				theory = new FSMTheorySolver(&S);
				S.addTheory(theory);
				theory->setStrings(strings);
			}
				theory->newFSM(fsmID);


				theory->setAlphabets(fsmID,inAlphabets[i],outAlphabets[i]);

				for (auto &t:transitions[i]){
					theory->newTransition(fsmID,t.from,t.to,t.input,t.output,t.edgeVar);
				}

				for(auto & a: accepts[i]){
					if (a.strID<0 ){
							parse_errorf("String ID must be a non-negative integer, was %d\n", a.strID);
						}
					stringIDMap.growTo(a.strID+1,-1);
					a.strID = stringIDMap[a.strID];

					if (a.strID<0  || a.strID>=created_strings.size() || !created_strings[a.strID]){
						parse_errorf("String ID must be a non-negative integer, was %d\n", a.strID);
					}

					if(a.from<0 || a.from>=theory->nNodes(fsmID)){
							parse_errorf("%d is not a valid state\n", a.from);
						}
					if(a.to<0 || a.to>=theory->nNodes(fsmID)){
						parse_errorf("%d is not a valid state\n", a.to);
					}

					theory->addAcceptLit(fsmID,a.from, a.to,a.strID,a.reachVar);
				}

				for(auto & a: generates[i]){
					if (a.strID<0 ){
							parse_errorf("String ID must be a non-negative integer, was %d\n", a.strID);
						}
					stringIDMap.growTo(a.strID+1,-1);
					a.strID = stringIDMap[a.strID];
					if (a.strID<0  || a.strID>=created_strings.size() || !created_strings[a.strID]){
						parse_errorf("String ID must be a non-negative integer, was %d\n", a.strID);
					}
					if(a.from<0 || a.from>=theory->nNodes(fsmID)){
							parse_errorf("%d is not a valid state\n", a.from);
						}

					theory->addGenerateLit(fsmID,a.from, a.strID,a.reachVar);
				}

				for(auto & a: transduces[i]){
					if (a.strID<0 ){
							parse_errorf("String ID must be a non-negative integer, was %d\n", a.strID);
						}
					if (a.strID2<0 ){
							parse_errorf("String ID must be a non-negative integer, was %d\n", a.strID2);
						}
					stringIDMap.growTo(a.strID+1,-1);
					a.strID = stringIDMap[a.strID];
					stringIDMap.growTo(a.strID2+1,-1);
					a.strID2 = stringIDMap[a.strID2];
					if (a.strID<0 || a.strID>=created_strings.size() || !created_strings[a.strID]){
						parse_errorf("String ID must be a non-negative integer, was %d\n", a.strID);
					}
					if(a.from<0 || a.from>=theory->nNodes(fsmID)){
							parse_errorf("%d is not a valid state\n", a.from);
						}
					if(a.to<0 || a.to>=theory->nNodes(fsmID)){
									parse_errorf("%d is not a valid state\n", a.from);
								}
					theory->addTransduceLit(fsmID,a.from,a.to, a.strID,a.strID2,a.reachVar);
				}



		}
		for(auto & c: compose_accepts){
			if(!theory){
				parse_errorf("No fsms declared!\n");

			}
			if (c.strID<0 ){
				parse_errorf("String ID must be a non-negative integer, was %d\n", c.strID);
			}
			stringIDMap.growTo(c.strID+1,-1);
			c.strID = stringIDMap[c.strID];
			if (c.strID<0 || c.strID>=created_strings.size() || !created_strings[c.strID]){
				parse_errorf("String ID must be a non-negative integer, was %d\n", c.strID);
			}
			theory->addComposeAcceptLit(c.fsmID1,c.fsmID2,c.from1,c.to1,c.from2,c.to2, c.strID,c.reachVar);
		}


	}

	
};

//=================================================================================================
}
;

#endif /* GRAPH_PARSER_H_ */
