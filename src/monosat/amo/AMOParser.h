/**************************************************************************************************
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

#include <cstdio>

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
class AMOParser : public Parser<B, Solver> {
    using Parser<B, Solver>::mapVar;

    vec<Lit> lits;

public:
    AMOParser() : Parser<B, Solver>("At-Most-One"){

    }

    bool parseLine(B& in, Solver& S){

        skipWhitespace(in);
        if(*in == EOF)
            return false;
        if(match(in, "amo")){
            //read in amo constraints the same way as clauses, except that all ints are interpreted as variables, not lits, for now:
            lits.clear();
            for(;;){
                int parsed_var = parseInt(in);
                if(parsed_var == 0)
                    break;
                bool sign = false;
                if(parsed_var < 0){
                    sign = true;
                    parsed_var = -parsed_var;
                }
                Var var = parsed_var - 1;
                var = mapVar(S, var);
                Lit l = mkLit(var, sign);
                lits.push(l);
            }
            if(lits.size() > 1){
                //else this constraint has no effect
                AMOTheory* amo = new AMOTheory(&S);
                for(Lit l : lits)
                    amo->addLit(l);

            }
            return true;
        }

        return false;
    }

    void implementConstraints(Solver& S){

    }


};

//=================================================================================================
};

#endif /* GRAPH_PARSER_H_ */
