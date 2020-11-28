/**************************************************************************************************

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
#include "monosat/core/Config.h"
#include "monosat/utils/ParseUtils.h"
#include "monosat/core/SolverTypes.h"
#include "monosat/core/Optimize.h"
#include "monosat/mtl/Vec.h"
#include <string>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <iostream>
#include <sstream>
#include <cstdarg>
#include <set>
#include "monosat/core/Remap.h"

namespace Monosat {

template<class B, class Solver>
class Parser {
    const char* parser_name;
protected:
    DimacsMap* dimacsParser = nullptr;
    BVMap* bvmap = nullptr;
public:

    Parser(const char* parser_name) : parser_name(parser_name){
    }

    virtual ~Parser(){
    }

    virtual bool parseLine(B& in, Solver& S) = 0;

    virtual void implementConstraints(Solver& S) = 0;

    const char* getParserName() const{
        return parser_name;
    }

    Var mapVar(Solver& S, Var v){
        return dimacsParser->mapVar(S, v);
    }

    Lit mapLit(Solver& S, Lit l){
        Var m = mapVar(S, var(l));
        return mkLit(m, sign(l));
    }

    void setDimacs(DimacsMap* dimacs){
        this->dimacsParser = dimacs;
    }

    void setBVMap(BVMap* bvmap){
        this->bvmap = bvmap;
    }

    bool inVarMap(Var externalVar){
        return dimacsParser->inVarMap(externalVar);
    }

    Var unmap(Var internalVar){
        return dimacsParser->unmap(internalVar);
    }

    void addVarToMap(Var v, Var map_to){
        dimacsParser->addVarToMap(v, map_to);
    }

    bool hasMappedVar(Var internalVar){
        return dimacsParser->hasMappedVar(internalVar);
    }

    Var getVarFromExternalVar(Var externalVar){
        return dimacsParser->getVarFromExternalVar(externalVar);
    }

    int dimacs(Var v){
        return dimacsParser->dimacs(v);
    }

    int dimacs(Lit l){
        return dimacsParser->dimacs(l);
    }


    inline int mapBV(Solver& S, int bv){
        return bvmap->mapBV(S, bv);
    }

    inline bool inBVMap(int externalBV){
        return bvmap->inBVMap(externalBV);
    }

    inline void addBVToMap(int bvID, int map_to){
        bvmap->addBVToMap(bvID, map_to);
    }


};


//A simple parser to allow for named variables
template<class B, class Solver>
class SymbolParser : public Parser<B, Solver> {
    using Parser<B, Solver>::mapVar;
    std::stringstream ss;

public:
    SymbolParser() : Parser<B, Solver>("Symbol"){

    }

    bool parseLine(B& in, Solver& S){
        if(*in == EOF)
            return false;
        else if(*in == 's'){
            if(match(in, "symbol")){ //used to use "c var" for symbols

                //this is a variable symbol map
                skipWhitespace(in);
                bool sign = false;
                int externalVar = parseInt(in);
                if(externalVar <= 0){
                    externalVar = -externalVar;
                    sign = true;
                }

                externalVar--; //subtract one to get the variable id

                ss.str(std::string());
                skipWhitespace(in);

                while(!isEof(in) && *in != '\n' && !isWhitespace(*in)){
                    ss << ((char) *in);
                    ++in;
                }

                std::string name_str = ss.str();

                if(name_str.size() == 0){
                    parse_errorf("Empty symbol: %c\n", *in);
                }

                Var internalVar = mapVar(S, externalVar);
                Lit internalLit = mkLit(internalVar, sign);
                S.addLiteralName(internalLit, name_str);

                return true;
            }else{
                //just a comment
                return false;
            }
        }
        return false;
    }

    void implementConstraints(Solver& S){

    }

};


//The MiniSAT DIMACS Parser, converted to be extensible...
template<class B, class Solver>
class Dimacs : public DimacsMap, public BVMap {
    vec<Parser<char*, Solver>*> parsers;


public:
    vec<Objective> objectives;
    vec<Lit> assumptions;

    Dimacs() : DimacsMap(opt_remap_vars), BVMap(opt_remap_vars){

    }

    ~Dimacs() override{
    }

    vec<Parser<char*, Solver>*>& getParsers(){
        return parsers;
    };


private:


    void readClause(B& in, Solver& S, vec<Lit>& lits){
        int parsed_lit, var;
        lits.clear();
        for(;;){
            parsed_lit = parseInt(in);
            if(parsed_lit == 0)
                break;
            var = abs(parsed_lit) - 1;
            var = mapVar(S, var);
            lits.push((parsed_lit > 0) ? mkLit(var) : ~mkLit(var));
        }
    }

    virtual bool parseLine(const char* line, int line_number, Solver& S){
        for(auto* p : parsers){
            char* ln = (char*) line; //intentionally discard const qualifier
            try{
                if(p->parseLine(ln, S)){
                    if(opt_parser_immediate_mode){
                        p->implementConstraints(S);
                    }
                    return true;
                }
            }catch(const parse_error& e){
                std::stringstream ss;
                ss << "PARSE ERROR in " << p->getParserName() << " parser at line " << line_number << ": " << line
                   << "; ";
                ss << e.what();
                throw parse_error(ss.str());
            }catch(const std::exception& e){
                std::stringstream ss;
                ss << "PARSE ERROR in " << p->getParserName() << " parser at line " << line_number << ": " << line
                   << "; ";
                ss << e.what();
                throw parse_error(ss.str());
            }catch(...){
                std::stringstream ss;
                ss << "PARSE ERROR in " << p->getParserName() << " parser at line " << line_number << ": " << line
                   << "; (Unknown Error)";
                throw parse_error(ss.str());
            }
        }
        return false;
    }

    bool readLine(vec<char>& linebuf, B& in){
        linebuf.clear();
        for(;;){
            if(isEof(in))
                return false;
            else if(*in == '\n'){
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

    int vars = 0;
    int clauses = 0;
    int clause_count = 0;
    int line_num = 0;
    int solves = 0;

    bool parse_(B& in, Solver& S){
        vec<Lit> lits;
        if(opt_remap_vars){
            S.setVarMap(this);
        }
        S.cancelUntil(0);
        objectives.clear();
        assumptions.clear();
        bool solve = false;
        vec<char> linebuf;
        try{
            while(!solve){
                skipWhitespace(in);
                if(*in == EOF)
                    break;
                line_num++;//this will merge line counts if there are multiple blank lines...
                //Typically, 99% of lines are either comments or clauses, and so it makes a lot of sense to handle these first, and before reading the whole line into a buffer.
                if(*in == '-' || (*in >= '0' && *in <= '9')){
                    //this is a clause
                    clause_count++;
                    readClause(in, S, lits);
                    S.addClause_(lits);
                    continue;
                }
                if(*in == 'c'){
                    skipLine(in);
                    continue;//comment
                }
                readLine(linebuf, in);
                char* b = linebuf.begin();
                if(match(b, "solve")){
                    int parsed_lit, var;
                    lits.clear();
                    for(;;){
                        while(*b == ' ')
                            ++b;
                        if(*b == '\n')
                            break;
                        parsed_lit = parseInt(b);
                        if(parsed_lit == 0)
                            break;
                        var = abs(parsed_lit) - 1;
                        var = mapVar(S, var);
                        assumptions.push((parsed_lit > 0) ? mkLit(var) : ~mkLit(var));
                    }
                    solves++;
                    solve = true;
                }else if(match(b, "dbg_solution")){
                    //add a known solution that
                    //all learnt clauses should be checked against
                    int parsed_lit, var;
                    lits.clear();
                    S.known_solutions.push();
                    for(;;){
                        while(*b == ' ')
                            ++b;
                        if(*b == '\n')
                            break;
                        parsed_lit = parseInt(b);
                        if(parsed_lit == 0)
                            break;
                        var = abs(parsed_lit) - 1;
                        var = mapVar(S, var);
                        S.known_solutions.last().push((parsed_lit > 0) ? mkLit(var) : ~mkLit(var));
                    }
                }else if(match(b, "priority")){
                    int parsed_int = parseInt(b);
                    int var = abs(parsed_int) - 1;
                    var = mapVar(S, var);
                    int priority = parseInt(b);
                    S.setDecisionPriority(var, priority);
                }else if(match(b, "decision")){
                    int parsed_int = parseInt(b);
                    int var = abs(parsed_int) - 1;
                    var = mapVar(S, var);
                    bool decision = parseInt(b);
                    S.setDecisionVar(var, decision);
                }else if(match(b, "clear_opt")){
                    objectives.clear();
                }else if(match(b, "minimize_core")){
                    int parsed_lit, var;
                    lits.clear();
                    for(;;){
                        while(*b == ' ')
                            ++b;
                        if(*b == '\n')
                            break;
                        parsed_lit = parseInt(b);
                        if(parsed_lit == 0)
                            break;
                        var = abs(parsed_lit) - 1;
                        var = mapVar(S, var);
                        //assumptions.push((parsed_lit > 0) ? mkLit(var) : ~mkLit(var));
                    }
                    //do nothing
                }else if(match(b, "minimize bv")){
                    //fprintf(stderr,"minimize statements not yet supported\n");
                    skipWhitespace(b);
                    int bvID = parseInt(b);
                    assert(bvID >= 0);

                    objectives.push(Objective(bvID, false));
                }else if(match(b, "maximize bv")){
                    //fprintf(stderr,"minimize statements not yet supported\n");
                    skipWhitespace(b);
                    int bvID = parseInt(b);
                    assert(bvID >= 0);

                    objectives.push(Objective(bvID, true));
                }else if(match(b, "maximize lits")){
                    //fprintf(stderr,"minimize statements not yet supported\n");
                    skipWhitespace(b);
                    objectives.push();
                    objectives.last().type = Objective::Type::PB;
                    objectives.last().maximize = true;
                    int n_lits = parseInt(b);
                    for(int i = 0; i < n_lits; i++){
                        int parsed_lit = parseInt(b);
                        if(parsed_lit == 0){
                            parse_errorf("Bad literal: 0\n");
                        }
                        Var var = abs(parsed_lit) - 1;
                        var = mapVar(S, var);
                        S.setFrozen(var, true);
                        objectives.last().pb_lits.push((parsed_lit > 0) ? mkLit(var) : ~mkLit(var));
                    }

                    for(int i = 0; i < n_lits; i++){
                        skipWhitespaceNoNewLines(b);
                        if(*b == '\n')
                            break;
                        int weight = parseInt(b);
                        objectives.last().pb_weights.push(weight);
                    }

                    if(objectives.last().pb_weights.size() > objectives.last().pb_lits.size()){
                        objectives.last().pb_weights.shrink(
                                objectives.last().pb_weights.size() - objectives.last().pb_lits.size());
                    }
                    objectives.last().pb_weights.growTo(objectives.last().pb_lits.size(), 1);
                }else if(match(b, "minimize lits")){
                    //fprintf(stderr,"minimize statements not yet supported\n");
                    skipWhitespace(b);
                    objectives.push();
                    objectives.last().type = Objective::Type::PB;
                    objectives.last().maximize = false;
                    int n_lits = parseInt(b);
                    for(int i = 0; i < n_lits; i++){
                        int parsed_lit = parseInt(b);
                        if(parsed_lit == 0){
                            parse_errorf("Bad literal: 0\n");
                        }
                        Var var = abs(parsed_lit) - 1;
                        var = mapVar(S, var);
                        S.setFrozen(var, true);
                        objectives.last().pb_lits.push((parsed_lit > 0) ? mkLit(var) : ~mkLit(var));
                    }

                    for(int i = 0; i < n_lits; i++){
                        skipWhitespaceNoNewLines(b);
                        if(*b == '\n')
                            break;
                        int weight = parseInt(b);
                        objectives.last().pb_weights.push(weight);
                    }

                    if(objectives.last().pb_weights.size() > objectives.last().pb_lits.size()){
                        objectives.last().pb_weights.shrink(
                                objectives.last().pb_weights.size() - objectives.last().pb_lits.size());
                    }
                    objectives.last().pb_weights.growTo(objectives.last().pb_lits.size(), 1);
                }else if(parseLine(b, line_num, S)){
                    //do nothing
                }else if(linebuf[0] == 'p'){

                    if(eagerMatch(b, "p cnf")){
                        vars = parseInt(b);
                        clauses = parseInt(b);
                    }else{
                        parse_errorf("Unexpected char: %c\n", *in);
                    }
                }else{
                    //if nothing else works, attempt to parse this line as a clause.
                    parse_errorf("Bad line at %d: %s", line_num, linebuf.begin());
                }
            }
            if(solve){
                //continue reading any blank/comment lines
                while(*in != EOF){
                    skipWhitespace(in);
                    if(*in == EOF)
                        break;
                    else if(*in == 'c'){
                        skipLine(in);
                        //continue
                    }else{
                        break;
                    }
                }
            }

            //Disabling this for now, as it is always triggered when there are theory atoms...
            /*if (vars != S.nVars())
                fprintf(stderr, "WARNING! DIMACS header mismatch: wrong number of variables.\n");
            if (cnt != clauses)
                fprintf(stderr, "WARNING! DIMACS header mismatch: wrong number of clauses.\n");*/
            for(auto* p : parsers){
                try{
                    p->implementConstraints(S);
                }catch(const std::exception& e){
                    std::stringstream ss;
                    ss << "PARSE ERROR in " << p->getParserName() << " parser ; ";
                    ss << e.what();
                    throw parse_error(ss.str());
                }catch(...){
                    std::stringstream ss;
                    ss << "PARSE ERROR in " << p->getParserName() << " parser ; (Unknown Error)";
                    throw parse_error(ss.str());
                }
            }
            for(int i = 0; i < objectives.size(); i++){
                if(objectives[i].isBV()){
                    int bvID = objectives[i].bvID;
                    bvID = this->mapBV(S, bvID);
                    objectives[i].bvID = bvID;
                }else{
                    //lits are already remapped
                }
            }

        }catch(const parse_error& e){
            std::stringstream ss;
            ss << "PARSE ERROR in DIMACS parser at line " << line_num << "; ";
            ss << e.what();
            throw parse_error(ss.str());
        }catch(const std::exception& e){
            std::stringstream ss;
            ss << "PARSE ERROR in DIMACS parser at line " << line_num << "; ";
            ss << e.what();
            throw parse_error(ss.str());
        }catch(...){
            std::stringstream ss;
            ss << "PARSE ERROR in DIMACS parser at line " << line_num << "; (Unknown Error)";
            throw parse_error(ss.str());
        }
        return solve;
    }

public:
    void addParser(Parser<char*, Solver>* parser){
        parser->setDimacs(this);
        parser->setBVMap(this);
        parsers.push(parser);
    }

/*	// Inserts problem into solver.
	bool parse_DIMACS(gzFile input_stream, Solver& S) {
		StreamBuffer in(input_stream);
		return parse(in, S);
	}*/
    bool parse(StreamBuffer& in, Solver& S){
        return parse_(in, S);
    }
};
};
#endif
