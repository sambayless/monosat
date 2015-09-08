
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


#ifndef BVPARSER_H_
#define BVPARSER_H_


#include <stdio.h>

#include "utils/ParseUtils.h"
#include "core/SolverTypes.h"
#include "bv/BVTheorySolver.h"

#include "core/Config.h"

#include "core/Dimacs.h"
#include <gmpxx.h>
#include <set>
#include <string>
#include <sstream>
namespace Monosat {


template<class B, class Solver>
class BVParser: public Parser<B, Solver> {
public:
	BVTheorySolver<long>* theory=nullptr;
private:
	vec<std::pair<int, std::string> >  symbols;
	std::string symbol;

	struct BV{
		int id=-1;
		int width=0;
		long constval=-1;
		vec<Var> vector;
	};
	vec<BV> bvs;



	vec<int> weights;
	vec<Lit> lits;
	int count = 0;
	vec<char> tmp;
	struct Compare{
		int bvID;
		long w;
		Comparison c;
		Var var;
	};
	vec<Compare> compares;
	struct CompareBV{
		int bvID;
		int compareID;
		Comparison c;
		Var var;
	};
	vec<CompareBV> comparebvs;

/*	struct AddConst{
		int resultID;
		int aBV;
		long b;
	};

	vec<AddConst> addconsts;*/
	struct AddBV{
		int resultID;
		int aBV;
		int bBV;

	};
	vec<AddBV> addbvs;

	struct IteBV{
		Lit condition;
		int thenId;
		int elseId;
		int resultId;
	};
	vec<IteBV> itebvs;

	void readConstBV(B& in,  Solver& S) {
		//bv id width l0 l1 l2 ...

		int id = parseInt(in);
		int width = parseInt(in);

		bvs.growTo(id + 1);

		if(bvs[id].id!=-1){

			parse_errorf("Re-defined bitvector %d\n", id);

		}
		bvs[id].id = id;
		bvs[id].width=width;
		bvs[id].constval=parseLong(in);
	}

	void readBV(B& in,  Solver& S) {
		//bv id width l0 l1 l2 ...

		int id = parseInt(in);
		int width = parseInt(in);

		bvs.growTo(id + 1);

		if(bvs[id].id!=-1){

			parse_errorf("Re-defined bitvector %d\n", id);

		}
		bvs[id].id = id;
		bvs[id].width=width;
		for(int i =0;i<width;i++){
			int v = parseInt(in) - 1;
			while (v >= S.nVars())
				S.newVar();
			bvs[id].vector.push(v);
		}
	}


	void readAddBV(B& in, Solver& S) {
		if (opt_ignore_theories) {
			skipLine(in);
			return;
		}
		skipWhitespace(in);
	/*	if(! match(in,"bv")){
			parse_errorf("Result of addition must be bitvector\n");
		}*/
		int resultID = parseInt(in);
		skipWhitespace(in);
		//bool arg1_is_bv=match(in,"bv");
		long arg1 = parseLong(in);
		skipWhitespace(in);
		//bool arg2_is_bv=match(in,"bv");
		long arg2 = parseLong(in);

		//if(arg1_is_bv && arg2_is_bv){
		if(arg2<arg1){

			std::swap(arg1,arg2);
		}

		addbvs.push();
		addbvs.last().resultID = resultID;
		addbvs.last().aBV =  (int) arg1;
		addbvs.last().bBV = (int) arg2;
		return;
	/*	}else if(!(arg1_is_bv || arg2_is_bv)){
			bvconstants.push({resultID,arg1+arg2});
			return;
		}else if (arg2_is_bv && ! arg1_is_bv){
			std::swap(arg1_is_bv,arg2_is_bv);
			std::swap(arg1,arg2);
		}
		addconsts.push();
		addconsts.last().resultID = resultID;
		addconsts.last().aBV = (int) arg1;
		addconsts.last().b = arg2;*/
	}

	void readSymbol(B& in, Solver& S){
		//this is a variable symbol map
		skipWhitespace(in);
		int v = parseInt(in);

		symbol.clear();
		skipWhitespace(in);
		while (*in != '\n' && !isWhitespace(*in)) {
			symbol.push_back(*in);
			++in;
		}
		if (symbol.size() == 0) {
			parse_errorf("Empty symbol: %c\n", *in);
		}
		/*   		if(symbols && used_symbols.count(symbol)){
		 parse_errorf("Duplicated symbol: %c\n", *symbol.c_str());
		 }
		 used_symbols.insert(symbol);*/

		symbols.push();
		symbols.last().first = v;
		symbols.last().second = symbol;
	}

	void readIteBV(B& in, Solver& S){
		//"bv ite %d %d %d %d\n"%(dimacs(condition_lit),aID,bID,resultID))
		int parsed_lit = parseInt(in);
		if (parsed_lit == 0)
			parse_errorf("If argument to bv If-Then-Else must be a valid dimacs literal (was 0)\n");
		int var = abs(parsed_lit) - 1;

		int thenId = parseInt(in);
		int elseId = parseInt(in);

		int resultId = parseInt(in);

		Lit l = mkLit(var,false);
		itebvs.push();
		itebvs.last().condition=l;
		itebvs.last().thenId=thenId;
		itebvs.last().elseId=elseId;
		itebvs.last().resultId=resultId;

	}

	void readCompareBV(B& in, Solver& S,Comparison c) {
		if (opt_ignore_theories) {
			skipLine(in);
			return;
		}
		//bv_lt bvID var weight
		skipWhitespace(in);
		int v = parseInt(in) - 1;
		while (v >= S.nVars())
			S.newVar();
		skipWhitespace(in);
		//bool arg1_is_bv=match(in,"bv");
		long arg1 = parseInt(in);
		skipWhitespace(in);
		//bool arg2_is_bv=match(in,"bv");
		long arg2 = parseInt(in);


/*
		if (arg2_is_bv && ! arg1_is_bv){
			c=~c;
			std::swap(arg1_is_bv,arg2_is_bv);
			std::swap(arg1,arg2);
		}*/
		//if(arg2_is_bv){
			comparebvs.push();
			comparebvs.last().bvID = (int) arg1;
			comparebvs.last().compareID = (int) arg2;
			comparebvs.last().c = c;
			comparebvs.last().var = v;
		/*}else{
			compares.push();
			compares.last().bvID =(int)  arg1;
			compares.last().w = arg2;
			compares.last().c = c;
			compares.last().var = v;
		}*/
	}

	void readCompare(B& in, Solver& S,Comparison c) {
			if (opt_ignore_theories) {
				skipLine(in);
				return;
			}
			//bv_lt bvID var weight
			skipWhitespace(in);
			int v = parseInt(in) - 1;
			while (v >= S.nVars())
				S.newVar();
			skipWhitespace(in);
			//bool arg1_is_bv=match(in,"bv");
			long arg1 = parseInt(in);
			skipWhitespace(in);
			//bool arg2_is_bv=match(in,"bv");
			long arg2 = parseLong(in);
			compares.push();
			compares.last().bvID =  arg1;
			compares.last().w = arg2;
			compares.last().c = c;
			compares.last().var = v;

		}

public:
	BVParser():Parser<B, Solver>("BitVector"){

	}
	bool parseLine(B& in, Solver& S) {

		skipWhitespace(in);
		if (*in == EOF)
			return false;
		else if (match(in, "c bv")) {
			//this is a bitvector symbol definition
			readSymbol(in,S);
			return true;
		}else if (*in == 'c') {
			//just a comment
			return false;
		} else if (match(in, "bv")) {
			skipWhitespace(in);
			if (match(in,"const")){
				skipWhitespace(in);
				if (match(in, "<=")) {
					readCompare(in, S,Comparison::leq);
					return true;
				}else if (match(in, "<")) {

					readCompare(in, S,Comparison::lt);
					return true;
				}else if (match(in, ">=")) {

					readCompare(in, S,Comparison::geq);
					return true;
				}  else if (match(in, ">")) {

					readCompare(in, S,Comparison::gt);
					return true;
				}else
					readConstBV(in,S);
				return true;
			}else if (match(in, "+")) {

				readAddBV(in, S);
				return true;
			}else if (match(in, "<=")) {

				readCompareBV(in, S,Comparison::leq);
				return true;
			}else if (match(in, "<")) {

				readCompareBV(in, S,Comparison::lt);
				return true;
			}else if (match(in, ">=")) {

				readCompareBV(in, S,Comparison::geq);
				return true;
			}  else if (match(in, ">")) {

				readCompareBV(in, S,Comparison::gt);
				return true;
			}else if (match(in,"ite")){
				readIteBV(in,S);
				return true;
			}else{
				readBV(in,S);
				return true;
			}

			return false;
		}
		return false;
	}

	void implementConstraints(Solver & S) {
		if(bvs.size()){
			theory = new BVTheorySolver<long>(&S);


			for (auto & bv:bvs){
				if(bv.id>-1){
					if(bv.constval>=0){
						theory->newBitvector(bv.id,bv.width,bv.constval);
					}else{
						theory->newBitvector(bv.id,bv.vector);
					}
				}
			}


			/*for(auto & c:addconsts){
				int newBV = theory->newBitvector(-1,theory->getBV(c.aBV).width()).getID();
				theory->makeConst(newBV,c.b);
				addbvs.push();
				addbvs.last().resultID = c.resultID;
				addbvs.last().aBV =  c.aBV;
				addbvs.last().bBV = newBV;
			}*/

			for(auto & c:compares){
				if(!theory->hasBV(c.bvID)){
					parse_errorf("Undefined bitvector ID %d",c.bvID);
				}
				theory->newComparison(c.c,c.bvID,c.w,c.var);
			}

			for(auto & c:comparebvs){
				if(!theory->hasBV(c.bvID)){
					parse_errorf("Undefined bitvector ID %d",c.bvID);
				}
				if(!theory->hasBV(c.compareID)){
					parse_errorf("Undefined bitvector ID %d",c.compareID);
				}

				theory->newComparisonBV(c.c,c.bvID,c.compareID,c.var);
			}

			for(auto & c:addbvs){
				if(!theory->hasBV(c.aBV)){
					parse_errorf("Undefined bitvector ID %d",c.aBV);
				}
				if(!theory->hasBV(c.bBV)){
					parse_errorf("Undefined bitvector ID %d",c.bBV);
				}
				if(!theory->hasBV(c.resultID)){
					parse_errorf("Undefined bitvector ID %d",c.resultID);
				}
				theory->newAdditionBV(c.resultID,c.aBV,c.bBV);
			}

			for (auto & c:itebvs){
				if(!theory->hasBV(c.thenId)){
					parse_errorf("Undefined bitvector ID %d",c.thenId);
				}
				if(!theory->hasBV(c.elseId)){
					parse_errorf("Undefined bitvector ID %d",c.elseId);
				}
				if(!theory->hasBV(c.resultId)){
					parse_errorf("Undefined bitvector ID %d",c.resultId);
				}
				theory->newConditionalBV(c.condition,c.thenId,c.elseId, c.resultId);
			}

			for (int i = 0; i < symbols.size(); i++) {
				int bvID = symbols[i].first;
				string & s = symbols[i].second;
				if (theory->hasBV(bvID)){
					theory->setSymbol(bvID, s.c_str());
				}else{
					fprintf(stderr,"Unmatched bv symbol definition for %d : %s\n",bvID,s.c_str());
				}
			}

		}else if (addbvs.size() || comparebvs.size() || compares.size()){

			parse_errorf("Undefined bitvector\n");

		}


	}


};

}
;




#endif /* BVPARSER_H_ */
