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

#ifndef Assumption_Parser_h
#define Assumption_Parser_h
#include <unordered_map>
#include <stdio.h>
#include <string>
#include "utils/ParseUtils.h"
#include "core/SolverTypes.h"
#include "mtl/Vec.h"
#include "mtl/Sort.h"
using namespace std;
namespace Monosat {

//=================================================================================================
// Assumption Parser:
template<class B>
static Lit readLit(B& in) {
	int parsed_lit, var;
	parsed_lit = parseInt(in);
	if (parsed_lit == 0) {
		return lit_Error;
	}
	var = abs(parsed_lit) - 1;
	
	return ((parsed_lit > 0) ? mkLit(var) : ~mkLit(var));
}

template<class Solver>
static void parse_Assumptions(gzFile input_stream, vec<Lit>& assumptions, Solver& S,
		vec<std::pair<int, std::string> > * symbols = NULL) {
	StreamBuffer in(input_stream);
	std::string symbol;
	std::unordered_map<std::string, int> symbol_table;
	if (symbols) {
		for (int i = 0; i < symbols->size(); i++) {
			int var = (*symbols)[i].first;
			string symbol = (*symbols)[i].second;
			symbol_table[symbol] = var;
		}
	}
	for (;;) {
		skipWhitespace(in);
		if (*in == EOF)
			break;
		else if (*in == 'c' || *in == 'p')
			skipLine(in);
		else {
			bool neg = false;
			if (*in == '-') {
				neg = true;
				++in;
			}
			if (isNumber(*in)) {
				Lit l = readLit(in);
				if (neg) {
					l = ~l;
				}
				assumptions.push(l);
			} else {
				//look this up in the symbol table
				symbol.clear();
				while (*in != '\n' && !isWhitespace(*in)) {
					symbol.push_back(*in);
					++in;
				}
				if (symbol.size() > 0) {
					if (symbol_table.count(symbol) == 0) {
						parse_errorf("Unknown symbol: %s\n", symbol.c_str()), exit(3);
					} else {
						int v = symbol_table[symbol];
						Lit l = mkLit(v, neg);
						assumptions.push(l);
					}
				} else {
					parse_errorf("Empty symbol!\n"), exit(3);
				}
				
			}
		}
	}
	sort(assumptions);
	
	Lit prev = lit_Undef;
	for (int i = 0; i < assumptions.size(); i++) {
		Lit l = assumptions[i];
		if (var(l) < 0) {
			parse_errorf("Bad assumption %d: %d is not a valid variable id\n", dimacs(l), var(l) + 1);
		} else if (var(l) >= S.nVars()) {
			parse_errorf("Bad assumption %d: %d is larger than the largest variable in the GNF\n", dimacs(l),
					var(l) + 1);
		} else if (l == ~prev) {
			fprintf(stderr,
					"Warning: Bad assumption %d: Variable %d assumed multiple times with conflicting polarity, this is unsolveable\n",
					dimacs(l), var(l) + 1);
			
		} else if (l == prev) {
			fprintf(stderr,
					"Warning: Variable %d assumed multiple times with the same polarity, this is probably an error\n",
					var(l) + 1);
		}
		prev = l;
	}
//=================================================================================================
}
}
;

#endif
