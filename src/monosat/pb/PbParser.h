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

#ifndef PB_PARSER_H_
#define PB_PARSER_H_

#include <cstdio>

#include "monosat/utils/ParseUtils.h"
#include "monosat/core/SolverTypes.h"
#include "monosat/pb/PbTheory.h"
#include "monosat/core/Config.h"
#include <sstream>
#include <set>
#include "PbSolver.h"
#include "monosat/pb/Pb.h"

namespace Monosat {


template<class B, class Solver>
class PBParser : public Parser<B, Solver> {
    using Parser<B, Solver>::mapVar;
    using Parser<B, Solver>::mapBV;
    PbTheory* pbtheory = nullptr;
    PB::PbSolver* pbsolver = nullptr;
    Solver& S;
public:
    PBParser(Solver& S) : Parser<B, Solver>("PB"), S(S){

    }

    void init(){
        if(opt_pb_theory){
            pbtheory = new PbTheory(&S);

        }else{
            if(!S.getPB()){
                pbsolver = new PB::PbSolver(S);
                S.setPBSolver(pbsolver);
            }else{
                pbsolver = (PB::PbSolver*) S.getPB();//fix this!
            }
        }
    }

    vec<Lit> lits;
    vec<int> weights;
    int vars = 0;
    int clauses = 0;
    vec<PB::Int> coefs;

    bool parseLine(B& in, Solver& S) override{

        skipWhitespace(in);
        if(*in == EOF)
            return false;
        else if(*in == 'c'){
            //just a comment
            return false;
        }else if(match(in, "pb")){
            init();
            readPB(in);
            return true;
        }
        return false;
    }

    void readPB(B& in){
        if(opt_ignore_theories){
            skipLine(in);
            return;
        }
        skipWhitespace(in);
        //pb constraints are in this form:
        //pb lt rhs <size> lit1 lit2 ... [0 | <n_weights> weight1 weight2 ...]
        PbTheory::PbType op = PbTheory::PbType::EQ;
        //read the operator:

        if(*in == '<'){
            ++in;
            if(*in == '='){
                ++in;
                op = PbTheory::PbType::LE;
            }else{
                op = PbTheory::PbType::LT;
            }
        }else if(*in == '>'){
            ++in;
            if(*in == '='){
                ++in;
                op = PbTheory::PbType::GE;
            }else{
                op = PbTheory::PbType::GT;
            }
        }else if(*in == '!'){
            ++in;
            if(*in != '='){
                parse_errorf("PARSE ERROR! Unexpected char: %c\n", *in);
            }
            ++in;
            op = PbTheory::PbType::NE;
        }else if(*in == '='){
            ++in;
            if(*in == '=')
                ++in; //second '=' is optional
            op = PbTheory::PbType::EQ;
        }else{
            parse_errorf("PARSE ERROR! Bad PB constraint\n");
        }
        if(op == PbTheory::PbType::NE && !opt_pb_theory){
            parse_errorf("PARSE ERROR! Minisat+ doesnt support dis-equalities (!=)\n");
        }

        lits.clear();
        weights.clear();
        int rhs = parseInt(in);
        int size = parseInt(in);
        if(size <= 0){
            parse_errorf("PARSE ERROR! Empty PB clause\n");
        }

        for(int i = 0; i < size; i++){
            int parsed_lit = parseInt(in);
            if(parsed_lit == 0)
                break;
            int v = abs(parsed_lit) - 1;
            v = mapVar(S, v);
            lits.push((parsed_lit > 0) ? mkLit(v) : ~mkLit(v));
        }

        int wsize = parseInt(in);
        if(wsize != 0 && wsize != size){
            parse_errorf("PARSE ERROR! Number of weights must either be the same as the size of the clause, or 0.\n");
        }
        for(int i = 0; i < wsize; i++){
            int parsed_weight = parseInt(in);
            weights.push(parsed_weight);
        }
        if(wsize == 0){
            for(int i = 0; i < size; i++)
                weights.push(1);
        }
        skipWhitespace(in);

        if(opt_pb_theory){
            assert(lits.size() == weights.size());
            pbtheory->addConstraint(lits, weights, rhs, lit_Undef, op,
                                    PbTheory::ConstraintSide::Both);
        }else{

            coefs.clear();
            for(int w:weights){
                coefs.push(PB::Int(w));
            }
            PB::Ineq ineq_op;
            if(op == PbTheory::PbType::EQ){
                ineq_op = PB::Ineq::EQ;
            }else if(op == PbTheory::PbType::LE){
                ineq_op = PB::Ineq::LEQ;
            }else if(op == PbTheory::PbType::LT){
                ineq_op = PB::Ineq::LT;
            }else if(op == PbTheory::PbType::GE){
                ineq_op = PB::Ineq::GEQ;
            }else if(op == PbTheory::PbType::GT){
                ineq_op = PB::Ineq::GT;
            }else if(op == PbTheory::PbType::NE){
                parse_errorf("PARSE ERROR! Minisat+ doesnt support dis-equalities (!=)\n");
            }

            pbsolver->addConstr(lits, coefs, PB::Int(rhs), ineq_op);
        }
    }

    void implementConstraints(Solver& S) override{
        if(pbtheory){
            pbtheory->implementConstraints();
        }
        if(pbsolver){
            pbsolver->convert();
        }
    }
};
//=================================================================================================
}

#endif /* GRAPH_PARSER_H_ */
