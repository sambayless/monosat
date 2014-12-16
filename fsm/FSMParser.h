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

	vec<FSMTheorySolver*> fsms;

	struct Transition{
		int fsm;
		int from;
		int to;
		int label;
		Var edgeVar;
	};
	vec<vec<Transition> > transitions;

	struct Accepts{
		int fsm;
		int state;
		Var acceptVar;
		vec<int> str;
	};
	vec<vec<Accepts> > accepts;
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
			printf("PARSE ERROR! FSM id must be >=0, was %d\n", fsmID), exit(1);
		}
		if (n_labels<0){
			printf("PARSE ERROR! Number of transition labels must be >=0, was %d\n", n_labels), exit(1);
		}

		fsms.growTo(fsmID + 1);
		if(fsms[fsmID]){
			printf("PARSE ERROR! FSM id %d declared twice!\n", fsmID), exit(1);
		}
		fsms[fsmID]= new FSMTheorySolver();
		transitions.growTo(fsmID+1);
		accepts.growTo(fsmID+1);
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
		int label = parseInt(in);
		int edgeVar = parseInt(in) - 1;
		
		if (fsmID < 0 || fsmID >= fsms.size()) {
			printf("PARSE ERROR! Undeclared fsm identifier %d for edge %d\n", graphID, edgeVar), exit(1);
		}
		if (label<0){
			printf("PARSE ERROR! Transition labels  must be >=0, was %d\n", label), exit(1);
		}

		if (edgeVar < 0) {
			printf("PARSE ERROR! Edge variables must be >=0, was %d\n", edgeVar), exit(1);
		}

		while (edgeVar >= S.nVars())
			S.newVar();
		
		transitions[fsmID].push({fsm,from,to,label,edgeVar});
	}

	void readAccepts(B& in, Solver& S) {
		if (opt_ignore_theories) {
			skipLine(in);
			return;
		}

		++in;
		
		int fsmID = parseInt(in);
		int acceptingState = parseInt(in);

		int reachVar = parseInt(in) - 1;

		//now read in the string
		accepts.push();

		accepts.back().fsm=fsmID;
		accepts.back().state=acceptingState;
		accepts.back().reachVar = reachVar;
		while(int i = parseInt(in)){
			if (i<=0){
				printf("PARSE ERROR! FSM strings must contain only positive (non-zero) integers, found %d\n", i), exit(1);
			}
			accepts.back().str.push(i);
		}

		if (fsmID < 0 || fsmID >= fsms.size()) {
			printf("PARSE ERROR! Undeclared fsm identifier %d for edge %d\n", graphID, reachVar), exit(1);
		}
		if (acceptingState<0){
			printf("PARSE ERROR! Accepting state must be a node id (a non-negative integer), was %d\n", acceptingState), exit(1);
		}
		if (reachVar < 0) {
			printf("PARSE ERROR! Edge variables must be >=0, was %d\n", reachVar), exit(1);
		}
		
		while (reachVar >= S.nVars())
			S.newVar();


	}
	
public:
	FSMParser(){
		
	}

	bool parseLine(B& in, Solver& S) {

		skipWhitespace(in);
		if (*in == EOF)
			return false;
		else if (*in == 'c') {
		if (match(in, "fsm")) {
			skipWhitespace(in);
			readFSM(in,S);
			skipWhitespace(in);

			return true;
		} else if (match(in, "transition")) {
			count++;
			readTransition(in, S);
			return true;
		} else if (match(in, "accepts")) {
			readAccepts(in, S);
			return true;
		}
		return false;
	}
	
	void implementConstraints(Solver & S) {
		



		for (int i = 0; i < fsms.size(); i++) {
			if (fsms[i])
				fsms[i]->implementConstraints();
		}

	}

	
};

//=================================================================================================
}
;

#endif /* GRAPH_PARSER_H_ */
