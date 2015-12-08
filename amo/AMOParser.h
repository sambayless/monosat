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

#ifndef AMO_PARSER_H_
#define AMO_PARSER_H_

#include <stdio.h>

#include "utils/ParseUtils.h"
#include "core/SolverTypes.h"
#include "amo/AMOTheory.h"
#include "core/Config.h"
#include "core/Dimacs.h"
#include <set>
#include <string>
#include <sstream>
namespace Monosat {


//=================================================================================================
// GRAPH Parser:
template<class B, class Solver>
class AMOParser: public Parser<B, Solver> {
	using Parser<B, Solver>::mapVar;

	vec<Var> vars;
	
public:
	AMOParser():Parser<B, Solver> ("At-Most-One") {
		
	}
	bool parseLine(B& in, Solver& S) {
		
		skipWhitespace(in);
		if (*in == EOF)
			return false;
		if (match(in,"amo")){
			//read in amo constraints the same way as clauses, except that all ints are interpreted as variables, not lits, for now:
			vars.clear();
			for (;;) {
				int parsed_var = parseInt(in);
				if (parsed_var == 0)
					break;
				if(parsed_var<0){
					fprintf(stderr,"at-most-one constraints currently only accept variables (eg, must be >0), not literals, but found %d\n",parsed_var);
					exit(1);
				}
				Var var = parsed_var-1;
				var= mapVar(S,var);
				vars.push(var);
			}
			if(vars.size()>1){
				//else this constraint has no effect
				AMOTheory * amo = new AMOTheory(&S);
				for (Var v : vars)
					amo->addVar(v);

			}
			return true;
		}

		return false;
	}
	
	void implementConstraints(Solver & S) {

	}

	
};

//=================================================================================================
}
;

#endif /* GRAPH_PARSER_H_ */
