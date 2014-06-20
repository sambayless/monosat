/*
 * graph_parser.h
 *
 *  Created on: 2013-06-08
 *      Author: sam
 */

#ifndef PB_PARSER_H_
#define PB_PARSER_H_


#include <stdio.h>

#include "utils/ParseUtils.h"
#include "core/SolverTypes.h"
#include "pb/PbTheory.h"
#include "core/Config.h"
#include <set>
namespace Minisat {


template<class B, class Solver>
static void parse_PB_main(B& in, Solver& S, vec<std::pair<int,std::string> > * symbols=NULL ) {
	PbTheory * theory = new PbTheory(&S);
    S.addTheory(theory);
	vec<Lit> lits;
	vec<int> weights;
    int vars    = 0;
    int clauses = 0;
    int cnt     = 0;
    vec<bool> hasSymbol;
    //std::set<std::string> used_symbols;
    int line =0;

    for (;;){
    	line++;
        skipWhitespace(in);
        if (*in == EOF) break;
        else if (*in == '*'){
        	skipLine(in);

        }else if (*in == 'm'){
        	fprintf(stderr,"PARSE ERROR! Minimization statements are not supported!\n");
        	fflush(stderr);
            exit(1);
        }else{
            cnt++;
            //a line is a sequence of pairs (integer, varID) followed by an operater and an integer and a semicolon


		   lits.clear();
		   weights.clear();
		   for (;;){
			   skipWhitespace(in);
			   if(*in=='>' || *in=='<' || *in=='='){
			        int rhs= 0;


			        PbTheory::PbType op;
			        if(*in=='>'){
			        	++in;
			        	if(*in=='='){
			        		op=PbTheory::PbType::GE;
			        		++in;
			        	}else{
			        		op=PbTheory::PbType::GT;
			        	}
			        }else if(*in=='<'){
			        	++in;
			        	if(*in=='='){
			        		op=PbTheory::PbType::LE;
			        		++in;
			        	}else{
			        		op=PbTheory::PbType::LT;
			        	}
			        }else{
			        	assert(*in=='=');
			        	++in;
			        	op=PbTheory::PbType::EQ;
			        }

			        skipWhitespace(in);
			        rhs = parseInt(in);

			        if(lits.size()==0){
			        	  fprintf(stderr,"Parse ERROR at line %d! Clause cannot be empty\n",line);
						   exit(1);
			        }

			        assert(lits.size()==weights.size());
			        theory->addConstraint(lits,weights,rhs,lit_Undef,op);
				   break;
			   }else{

				   if(*in=='+')
					   ++in;//skip this
				   int parsed_weight = parseInt(in);
				   weights.push(parsed_weight);

				   skipWhitespace(in);
				   //A variableID is the letter 'x' followed by an integer
				   if(*in=='x'){
					   ++in;
				   }else{
					   fprintf(stderr,"Parse ERROR at line %d! Expected variable id\n",line);
					   exit(1);
				   }

				   int parsed_ID = parseInt(in);
				   if (parsed_ID <= 0){
					   fprintf(stderr,"Parse ERROR  at line %d! variable ids must be positive integers of the form 'xN'.\n",line);
						exit(1);
				   }
				   Var v = parsed_ID - 1;

				   while (v >= S.nVars()) S.newVar();
				   lits.push(mkLit(v) );

				   if(symbols){
					   hasSymbol.growTo(v+1);
					   if(!hasSymbol[v]){
						   hasSymbol[v]=true;

							stringstream ss;
							ss<<"x" << parsed_ID;

							symbols->push();
							symbols->last().first=v;
							symbols->last().second=ss.str();
					   }
				   }

			   }
		   }
		   skipWhitespace(in);
		   if(*in !=';'){
			   fprintf(stderr,"Parse ERROR  at line %d! Line must end with semicolon\n",line);
				exit(1);
		   }
		   ++in;
        }

    }
    theory->implementConstraints();

}

// Inserts problem into solver.
//
template<class Solver>
static void parse_PB(gzFile input_stream, Solver& S, vec<std::pair<int,std::string> > * symbols=NULL) {
    StreamBuffer in(input_stream);
    parse_PB_main(in, S,symbols);

}

//=================================================================================================
}



#endif /* GRAPH_PARSER_H_ */
