
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


#ifndef COMPARISONPARSER_H_
#define COMPARISONPARSER_H_


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
class ComparisonParser: public Parser<B, Solver> {
public:
	BVTheorySolver<long>* theory=nullptr;
private:
	struct BV{
		int id=-1;
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

	void readBV(B& in,  Solver& S) {
		//bv id width l0 l1 l2 ...
		++in;
		int id = parseInt(in);
		int width = parseInt(in);

		bvs.growTo(id + 1);

		if(bvs[id].id!=-1){

			printf("PARSE ERROR! Re-defined bitvector %d\n", id), exit(1);

		}
		bvs[id].id = id;
		for(int i =0;i<width;i++){
			int v = parseInt(in) - 1;
			while (v >= S.nVars())
				S.newVar();
			bvs[id].vector.push(v);
		}
	}

	void readCompareBV(B& in, Solver& S,Comparison c) {
		if (opt_ignore_theories) {
			skipLine(in);
			return;
		}
		//bv_lt bvID var weight
		++in;

		int bvID = parseInt(in);
		int v = parseInt(in) - 1;

		while (v >= S.nVars())
			S.newVar();

		int compareID = parseInt(in);

		comparebvs.push();
		comparebvs.last().bvID = bvID;
		comparebvs.last().compareID = compareID;
		comparebvs.last().c = c;
		comparebvs.last().var = v;
	}

	void readCompare(B& in, Solver& S,Comparison c) {
		if (opt_ignore_theories) {
			skipLine(in);
			return;
		}
		//bv_lt bvID var weight
		++in;

		int bvID = parseInt(in);
		int v = parseInt(in) - 1;

		while (v >= S.nVars())
			S.newVar();

		long weight = parseLong(in);

		compares.push();
		compares.last().bvID = bvID;
		compares.last().w = weight;
		compares.last().c = c;
		compares.last().var = v;
	}

public:
	ComparisonParser(){

	}
	bool parseLine(B& in, Solver& S) {

		skipWhitespace(in);
		if (*in == EOF)
			return false;
		else if (*in == 'c') {
			//just a comment
			return false;
		} else if (match(in, "bv")) {
			if (match(in, "_lt_bv")) {
				count++;
				readCompareBV(in, S,Comparison::lt);
				return true;
			} else if (match(in, "_leq_bv")) {
				count++;
				readCompareBV(in, S,Comparison::leq);
				return true;
			} else if (match(in, "_gt_bv")) {
				count++;
				readCompareBV(in, S,Comparison::gt);
				return true;
			} else if (match(in, "_geq_bv")) {
				count++;
				readCompareBV(in, S,Comparison::geq);
				return true;
			}else if (match(in, "_lt")) {
				count++;
				readCompare(in, S,Comparison::lt);
				return true;
			} else if (match(in, "_leq")) {
				count++;
				readCompare(in, S,Comparison::leq);
				return true;
			} else if (match(in, "_gt")) {
				count++;
				readCompare(in, S,Comparison::gt);
				return true;
			} else if (match(in, "_geq")) {
				count++;
				readCompare(in, S,Comparison::geq);
				return true;
			}else{
				readBV(in,S);
			}

			return true;
		}
		return false;
	}

	void implementConstraints(Solver & S) {
		if(compares.size() || bvs.size()){
			theory = new BVTheorySolver<long>(&S);


			for (auto & bv:bvs){
				if(bv.id>-1){
					theory->newBitvector(bv.id,bv.vector);
				}
			}

			for(auto & c:compares){
				theory->newComparison(c.c,c.bvID,c.w,c.var);
			}

			for(auto & c:comparebvs){
				theory->newComparisonBV(c.c,c.bvID,c.compareID,c.var);
			}

		}


	}


};

}
;




#endif /* COMPARISONPARSER_H_ */
