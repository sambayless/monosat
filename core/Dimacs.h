/****************************************************************************************[Solver.h]

 MonoSAT -- Copyright (c) 2014, Sam Bayless
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
#include <string>
#include <algorithm>
namespace Monosat {

template<class B, class Solver>
class Parser {
public:
	Parser() {
	}
	virtual ~Parser() {
	}
	virtual bool parseLine(B& in, Solver& S)=0;
	virtual void implementConstraints(Solver & S)=0;
};

//A simple parser to allow for named variables
template<class B, class Solver>
class SymbolParser: public Parser<B, Solver> {
	vec<std::pair<int, std::string> >  symbols;
	std::string symbol;
public:
	SymbolParser(){

	}

	bool parseLine(B& in, Solver& S) {
		if (*in == EOF)
			return false;
		else if (*in == 'c') {
			if (match(in, "c var ")) {

				//this is a variable symbol map
				skipWhitespace(in);
				int v = parseInt(in);
				if (v <= 0) {
					//printf("PARSE ERROR! Variables must be positive: %c\n", *in), exit(1);
					v = -v;
				}

				v--; //subtract one to get the variable id

				symbol.clear();
				skipWhitespace(in);
				while (*in != '\n' && !isWhitespace(*in)) {
					symbol.push_back(*in);
					++in;
				}
				if (symbol.size() == 0) {
					printf("PARSE ERROR! Empty symbol: %c\n", *in), exit(1);
				}
				/*   		if(symbols && used_symbols.count(symbol)){
				 printf("PARSE ERROR! Duplicated symbol: %c\n", *symbol.c_str()), exit(1);
				 }
				 used_symbols.insert(symbol);*/

				symbols.push();
				symbols.last().first = v;
				symbols.last().second = symbol;
				return true;
			} else {
				//just a comment
				return false;
			}
		}
		return false;
	}
	void implementConstraints(Solver & S){
		//clear any unmapped symbols (have to do this _before_ implementing constraints, which may introduce new variables)
		int i, j = 0;
		for (i = 0; i < symbols.size(); i++) {
			std::pair<int, std::string> p = symbols[i];
			Var v = p.first;
			if (v <= S.nVars()) {
				//keep this symbol
				symbols[j++] = symbols[i];
			}
		}
		symbols.shrink(i - j);
	}
	vec<std::pair<int, std::string> > &  getSymbols(){
		return symbols;
	}
};

//The MiniSAT DIMACS Parser, converted to be extensible...
template<class B, class Solver>
class Dimacs {
	vec<Parser<char*, Solver>*> parsers;
public:
	Dimacs() {
		
	}
	
	virtual ~Dimacs() {
	}

private:
	
	void readClause(char * in, Solver& S, vec<Lit>& lits) {
		int parsed_lit, var;
		lits.clear();
		for (;;) {
			parsed_lit = parseInt(in);
			if (parsed_lit == 0)
				break;
			var = abs(parsed_lit) - 1;
			while (var >= S.nVars())
				S.newVar();
			lits.push((parsed_lit > 0) ? mkLit(var) : ~mkLit(var));
		}
	}
	
	virtual bool parseLine(const char * line, Solver& S) {
		for (auto * p : parsers) {
			char * ln = (char*) line; //intentionally discard const qualifier
			if (p->parseLine(ln, S)) {
				return true;
			}
		}
		return false;
	}
	
	bool readLine(vec<char> & linebuf, B& in) {
		linebuf.clear();
		for (;;) {
			if (isEof(in))
				return false;
			else if (*in == '\n') {
				linebuf.push(*in);
				++in;
				break;
			} else {
				linebuf.push(*in);
				++in;
			}
		}
		linebuf.push(0);
		return true;
	}
	
	void parse(B& in, Solver& S) {
		vec<Lit> lits;
		int vars = 0;
		int clauses = 0;
		int cnt = 0;
		
		vec<char> linebuf;
		for (;;) {
			skipWhitespace(in);
			if (*in == EOF)
				break;
			readLine(linebuf, in);
			if (parseLine(linebuf.begin(), S)) {
				//do nothing
			} else if (linebuf[0] == 'p') {
				char * b = linebuf.begin();
				if (eagerMatch(b, "p cnf")) {
					vars = parseInt(b);
					clauses = parseInt(b);
				} else {
					printf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
				}
			} else if (linebuf[0] == 'c') {
				//comment line
				//skipLine(in);
			} else {
				//if nothing else works, attempt to parse this line as a clause.
				cnt++;
				readClause(linebuf.begin(), S, lits);
				S.addClause_(lits);
			}
		}
		//Disabling this for now, as it is always triggered when there are theory atoms...
		/*if (vars != S.nVars())
			fprintf(stderr, "WARNING! DIMACS header mismatch: wrong number of variables.\n");
		if (cnt != clauses)
			fprintf(stderr, "WARNING! DIMACS header mismatch: wrong number of clauses.\n");*/
		
		for (auto * p : parsers) {
			p->implementConstraints(S);
		}
	}
public:
	void addParser(Parser<char*, Solver> * parser) {
		parsers.push(parser);
	}
	// Inserts problem into solver.
	void parse_DIMACS(gzFile input_stream, Solver& S) {
		StreamBuffer in(input_stream);
		parse(in, S);
	}
	
};
}
;
#endif
