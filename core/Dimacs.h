/****************************************************************************************[Dimacs.h]
Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
Copyright (c) 2007-2010, Niklas Sorensson

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

#ifndef Minisat_Dimacs_h
#define Minisat_Dimacs_h

#include <stdio.h>

#include "utils/ParseUtils.h"
#include "core/SolverTypes.h"
#include "mtl/Vec.h"
namespace Minisat {

template<class B, class Solver>
class Parser{
public:
	Parser(){}
	virtual ~Parser(){}
	virtual bool parseLine(B& in, Solver& S)=0;
	virtual void implementConstraints(Solver & S)=0;
};

//The MiniSAT DIMACS Parser, converted to be extensible...
template<class B, class Solver>
class Dimacs{
	vec<Parser<char*,Solver>*> parsers;
public:
	Dimacs(){

	}

	virtual ~Dimacs(){}

private:

	void readClause(char * in, Solver& S, vec<Lit>& lits) {
		int     parsed_lit, var;
		lits.clear();
		for (;;){
			parsed_lit = parseInt(in);
			if (parsed_lit == 0) break;
			var = abs(parsed_lit)-1;
			while (var >= S.nVars()) S.newVar();
			lits.push( (parsed_lit > 0) ? mkLit(var) : ~mkLit(var) );
		}
	}

	virtual bool parseLine(const char * line, Solver& S){
		for(auto * p:parsers){
			char * ln = (char*)line;//intentionally discard const qualifier
			if(p->parseLine(ln,S)){
				return true;
			}
		}
		return false;
	}

	bool readLine(vec<char> & linebuf, B& in){
		linebuf.clear();
		  for (;;){
			if (isEof(in))
				return false;
			else if (*in == '\n') {
				linebuf.push(*in);
				++in;
				break;
			}else{
				linebuf.push(*in);
				++in;
			}
		  }
		  linebuf.push(0);
		  return true;
	}

	void parse(B& in, Solver& S) {
		vec<Lit> lits;
		int vars    = 0;
		int clauses = 0;
		int cnt     = 0;

		vec<char> linebuf;
		for (;;){
			skipWhitespace(in);
			if (*in == EOF)
				break;
			readLine(linebuf,in);
			if (parseLine(linebuf.begin(),S)){
				//do nothing
			}else if (linebuf[0] == 'p'){
				char * b =linebuf.begin();
				if (eagerMatch(b, "p cnf")){
					vars    = parseInt(b);
					clauses = parseInt(b);
				}else{
					printf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
				}
			}else if (linebuf[0] == 'c'){
				//comment line
				//skipLine(in);
			}else{
				//if nothing else works, attempt to parse this line as a clause.
				cnt++;
				readClause(linebuf.begin(), S, lits);
				S.addClause_(lits);
			}
		}
		if (vars != S.nVars())
			fprintf(stderr, "WARNING! DIMACS header mismatch: wrong number of variables.\n");
		if (cnt  != clauses)
			fprintf(stderr, "WARNING! DIMACS header mismatch: wrong number of clauses.\n");

		for(auto * p:parsers){
			p->implementConstraints(S);
		}
	}
public:
	void addParser(Parser<char*,Solver> * parser ){
		parsers.push(parser);
	}
	// Inserts problem into solver.
	void parse_DIMACS(gzFile input_stream, Solver& S) {
		StreamBuffer in(input_stream);
		parse(in, S);
	}

};
};
#endif
