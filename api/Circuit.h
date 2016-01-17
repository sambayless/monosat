/****************************************************************************************[Solver.h]
 The MIT License (MIT)

 Copyright (c) 2016, Sam Bayless

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

#ifndef CIRCUIT_H_
#define CIRCUIT_H_

#include "core/SolverTypes.h"
#include "mtl/Vec.h"
#include <list>

namespace Monosat{

//Barebones helper methods for expressing combinatorial logic in CNF.
template<class Solver>
class Circuit{
	Solver & S;
	Lit lit_True;

	bool isConst(Lit l){
		return isConstTrue(l) || isConstFalse(l);
	}

	bool isConstTrue(Lit l){
		return S.value(l)==l_True && S.level(var(l))==0;
	}
	bool isConstFalse(Lit l){
		return S.value(l)==l_False && S.level(var(l))==0;
	}
	vec<Lit> tmp;
	vec<Lit> store;
	template<typename... Args>
	void collect(vec<Lit> & store, Lit a, Args... args ){
		store.push(a);
		collect(store,args...);
	}
	void collect(vec<Lit> & store, Lit a){
		store.push(a);
	}


	Lit bin_op(vec<Lit> & store,Lit (Circuit::*f)(Lit,Lit)){
		int n = store.size();
		while(n>1){
			int p = 0;
			for(int i = 0;i<n;i+=2){
				if(i+1<n){
					store[p++] = (this->**&f)(store[i],store[i+1]);
				}else{
					store[p++] = store[i];
				}
			}
			n=p;
		}
		Lit a = store[0];
		store.clear();
		return a;
	}
public:

	Circuit(Solver & S):S(S){
		lit_True = mkLit(S.newVar());
		S.addClause(lit_True);
	}


	Lit True(){
		return lit_True;
	}

	Lit False(){
		return ~lit_True;
	}

	Lit And(Lit a)
	{
		return a;
	}

	Lit And(Lit a, Lit b){
		//special case these
		if(isConst(a) || isConst(b)){
			if (isConstTrue(a)){
				return b;
			}else if (isConstFalse(a)){
				return a;
			}else if (isConstTrue(b)){
				return a;
			}else if (isConstFalse(b)){
				return b;
			}
		}

		Lit out = mkLit(S.newVar());
		S.addClause(a,~out);
		S.addClause(b,~out);
		S.addClause(~a,~b,out);
		return out;
	}
	Lit And(const std::list<Lit> & vals){

		tmp.clear();
		for (Lit l:vals)
			tmp.push(l);
		return bin_op(tmp,&Circuit::And);

	}
	Lit And(const vec<Lit> & vals){

		tmp.clear();
		for (Lit l:vals)
			tmp.push(l);
		return bin_op(tmp,&Circuit::And);
	}

	template<typename... Args>
	Lit And(Lit a, Lit b, Args... args)
	{
		store.clear();
		store.push(a);
		collect(store,b,args...);

		return And(store);
	}



	Lit Or(Lit a)
	{
		return a;
	}

	Lit Or(Lit a, Lit b){
		Lit out = mkLit(S.newVar());
		S.addClause(~a,out);
		S.addClause(~b,out);
		S.addClause(a,b,~out);
		return out;
	}

	Lit Or(const std::list<Lit> & vals){

		tmp.clear();
		for (Lit l:vals)
			tmp.push(l);
		return bin_op(tmp,&Circuit::Or);

	}
	Lit Or(const vec<Lit> & vals){
		tmp.clear();
		for (Lit l:vals)
			tmp.push(l);
		return bin_op(tmp,&Circuit::Or);
	}

	template<typename... Args>
	Lit Or(Lit a, Lit b, Args... args)
	{
		store.clear();
		store.push(a);
		collect(store,b,args...);

		return Or(store);
	}



	Lit Nor(Lit a)
	{
		return ~a;
	}

	Lit Nor(Lit a, Lit b){
		return ~Or(a,b);
	}

	Lit Nor(const std::list<Lit> & vals){
		return ~Or(vals);
	}
	Lit Nor(const vec<Lit> & vals){
		return ~Or(vals);
	}

	template<typename... Args>
	Lit Nor(Lit a, Lit b, Args... args)	{
		store.clear();
		store.push(a);
		collect(store,b,args...);

		return Nor(store);
	}


	Lit Nand(Lit a)
	{
		return ~a;
	}

	Lit Nand(Lit a, Lit b){
		return ~And(a,b);
	}
	Lit Nand(const std::list<Lit> & vals){
		return ~And(vals);

	}
	Lit Nand(const vec<Lit> & vals){
		return ~And(vals);
	}

	template<typename... Args>
	Lit Nand(Lit a, Lit b, Args... args)
	{
		store.clear();
		store.push(a);
		collect(store,b,args...);
		return Nand(store);
	}

	Lit Xor(Lit a)
	{
		return ~a;
	}

	Lit Xor(Lit a, Lit b){
		if(isConst(a) || isConst(b)){
			if (isConstTrue(a)){
				return ~b;
			}else if (isConstFalse(a)){
				return b;
			}else if (isConstTrue(b)){
				return ~a;
			}else if (isConstFalse(b)){
				return a;
			}
		}
		//return Or(And(a, ~b), And(~a,b));
		Lit out = mkLit(S.newVar());
		S.addClause(a,b,out);
		S.addClause(~a,b,~out);
		S.addClause(a,~b,~out);
		S.addClause(~a,~b,out);
		return out;
	}

	Lit Xor(const std::list<Lit> & vals){
		tmp.clear();
		for (Lit l:vals)
			tmp.push(l);
		return bin_op(tmp,&Circuit::Xor);

	}
	Lit Xor(const vec<Lit> & vals){
		tmp.clear();
		for (Lit l:vals)
			tmp.push(l);
		return bin_op(tmp,&Circuit::Xor);
	}

	template<typename... Args>
	Lit Xor(Lit a, Lit b, Args... args)
	{
		store.clear();
		store.push(a);
		collect(store,b,args...);
		return Xor(store);
	}


	Lit Xnor(Lit a)
	{
		return a;
	}
	Lit Xnor(Lit a, Lit b){
		return ~Xor(a,b);
	}
	Lit Xnor(const std::list<Lit> & vals){
		return ~Xor(vals);

	}
	Lit Xnor(const vec<Lit> & vals){
		return ~Xor(vals);
	}

	template<typename... Args>
	Lit Xnor(Lit a, Lit b, Args... args)
	{
		store.clear();
		store.push(a);
		collect(store,b,args...);
		return Xnor(store);
	}





	Lit Ite(Lit cond, Lit thn, Lit els){
		Lit l = ~And(cond,~thn);
		Lit r = ~And(~cond,~els);
		return And(l,r);
	}


};

};

#endif /* CIRCUIT_H_ */
