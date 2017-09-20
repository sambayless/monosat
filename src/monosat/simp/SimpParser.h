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

#ifndef SIMP_PARSER_H_
#define SIMP_PARSER_H_

#include <stdio.h>

#include "monosat/utils/ParseUtils.h"
#include "monosat/core/SolverTypes.h"
#include "monosat/amo/AMOTheory.h"
#include "monosat/core/Config.h"
#include "monosat/core/Dimacs.h"
#include <set>
#include <string>
#include <sstream>
namespace Monosat {


//=================================================================================================
// GRAPH Parser:
template<class B, class Solver>
class SimpParser: public Parser<B, Solver> {
	using Parser<B, Solver>::mapVar;

	vec<Var> vars;
	
public:
	SimpParser():Parser<B, Solver> ("Preprocessor") {
		
	}
	bool parseLine(B& in, Solver& S) {
		
		skipWhitespace(in);
		if (*in == EOF)
			return false;
		if (match(in,"disable_pre")){
			if(!opt_respect_preprocessor_gnf_directives)
				return true;
			S.disablePreprocessing();
			return true;
		}else if (match(in,"disallow_simp")){
			if(!opt_respect_preprocessor_gnf_directives)
				return true;
			int parsed_int = parseInt(in);
			int var = abs(parsed_int)-1;
			var = mapVar(S,var);

			if(S.isEliminated(var)){
				fprintf(stderr,"Warning: Var %d has already been eliminated by the pre-processor\n", var+1);
				return false;
			}else {
				S.setFrozen(var, true);
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
