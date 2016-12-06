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

#include "monosat/core/SolverTypes.h"
#include "monosat/mtl/Vec.h"
#include <list>

namespace Monosat{

//Barebones helper methods for expressing combinatorial logic in CNF.
template<class Solver>
class Circuit{
	Solver & S;
	Lit lit_True=lit_Undef;

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
	vec<Lit> tmp2;
	vec<Lit> clause;
	vec<Lit> store;
	template<typename... Args>
	void collect(vec<Lit> & store, Lit a, Args... args ){
		store.push(a);
		collect(store,args...);
	}
	void collect(vec<Lit> & store, Lit a){
		store.push(a);
	}

	//Note: a vector of size zero will always return lit_True
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
		Lit a;
		if(store.size()){
			a = store[0];
			store.clear();
		}else{
			a = lit_True;
		}
		return a;
	}
public:

	Circuit(Solver & S):S(S){
		lit_True = mkLit(S.newVar());
		S.addClause(lit_True);
	}

	Solver & getSolver(){
		return S;
	}
	Lit getTrue(){
		return lit_True;
	}

	Lit getFalse(){
		return ~lit_True;
	}

	Lit And(Lit a)
	{
		return a;
	}

	Lit And_(Lit a, Lit b, Lit out){
		//special case these
		if(isConst(a) || isConst(b)){
			if(isConstTrue(a) && isConstTrue(b)){
				if(out!=lit_Undef) {
					Assert(out);
				}
				return getTrue();
			}else if (isConstFalse(a) || isConstFalse(b)){
				if(out!=lit_Undef) {
					Assert(~out);
				}
				return getFalse();
			}

			if (isConstFalse(a)){
				if(out!=lit_Undef) {
					Assert(~out);
				}
				return getFalse();
			}else if (isConstFalse(b)){
				if(out!=lit_Undef) {
					Assert(~out);
				}
				return getFalse();
			}else if (isConstTrue(a)){
				if(out!=lit_Undef) {
					AssertEqual(b,out);
				}
				return b;
			}else if (isConstTrue(b)){
				if(out!=lit_Undef) {
					AssertEqual(a,out);
				}
				return a;
			}
		}
		if(out==lit_Undef){
			out = mkLit(S.newVar());
		}
		S.addClause(a,~out);
		S.addClause(b,~out);
		S.addClause(~a,~b,out);
		return out;
	}

	Lit And_(const std::list<Lit> & vals, Lit out){

		tmp2.clear();
		for (Lit l:vals)
			tmp2.push(l);
		return And_(tmp2,out);

	}
	Lit And_(const vec<Lit> & vals, Lit out){
		assert(tmp.size()==0);
		for(Lit l:vals){
			if (isConstFalse(l)){
				if(out!=lit_Undef){
					Assert(~out);
				}
				return getFalse();
			}else if (isConstTrue(l)){
				//leave literal out
			}else{
				tmp.push(l);
			}
		}
		//all arguments are constant true
		if(tmp.size()==0){
			if(out!=lit_Undef){
				Assert(out);
			}
			return getTrue();
		}else if (tmp.size()==1){
			if(out!=lit_Undef){
				AssertEqual(tmp[0],out);
			}
			return tmp[0];
		}
		if(out==lit_Undef){
			out = mkLit(S.newVar());
		}
		for(Lit l:tmp){
			S.addClause(l,~out);
		}
		for(int i = 0;i<tmp.size();i++){
			tmp[i]=~tmp[i];
		}
		tmp.push(out);
		S.addClause(tmp);
		tmp.clear();
		return out;
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
		tmp2.clear();
		for (Lit l:vals)
			tmp2.push(l);
		return And(tmp2);

	}
	Lit And(const vec<Lit> & vals){
		return And_(vals,lit_Undef);
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
		if(isConstTrue(a)){
			return getTrue();
		}else if (isConstTrue(b)){
			return getTrue();
		}else if(isConstFalse(a)){
			return b;
		}else if (isConstFalse(b)){
			return a;
		}

		Lit out = mkLit(S.newVar());
		S.addClause(~a,out);
		S.addClause(~b,out);
		S.addClause(a,b,~out);
		return out;
	}
	Lit Or_(Lit a, Lit b, Lit out){

		//special case these
		if(isConst(a) || isConst(b)){

			if(isConstFalse(a) && isConstFalse(b)){
				if(out!=lit_Undef) {
					Assert(~out);
				}
				return getFalse();
			}else if (isConstTrue(a) || isConstTrue(b)){
				if(out!=lit_Undef) {
					Assert(out);
				}
				return getTrue();
			}

			if (isConstTrue(a)){
				if(out!=lit_Undef) {
					Assert(out);
				}
				return getTrue();
			}else if (isConstTrue(b)){
				if(out!=lit_Undef) {
					Assert(out);
				}
				return getTrue();
			}else if (isConstFalse(a)){
				if(out!=lit_Undef) {
					AssertEqual(b,out);
				}
				return b;
			}else if (isConstFalse(b)){
				if(out!=lit_Undef) {
					AssertEqual(a,out);
				}
				return a;
			}
		}
		if(out==lit_Undef){
			out = mkLit(S.newVar());
		}
		S.addClause(~a,out);
		S.addClause(~b,out);
		S.addClause(a,b,~out);
		return out;
	}
	Lit Or_(const vec<Lit> & vals, Lit out){
		assert(tmp.size()==0);
		for(Lit l:vals){
			if (isConstFalse(l)){
				//leave literal out
			}else if (isConstTrue(l)){
				if(out!=lit_Undef){
					Assert(out);
				}
				return getTrue();
			}else{
				tmp.push(l);
			}
		}
		//all arguments are constant true
		if(tmp.size()==0){
			if(out!=lit_Undef){
				Assert(~out);
			}
			return getFalse();//or should this be true?
		}else if (tmp.size()==1){
			if(out!=lit_Undef){
				AssertEqual(tmp[0],out);
			}
			return tmp[0];
		}
		if(out==lit_Undef){
			out = mkLit(S.newVar());
		}
		for(Lit l:tmp){
			S.addClause(~l,out);
		}

		tmp.push(~out);
		S.addClause(tmp);
		tmp.clear();
		return out;
	}

	Lit Or(const std::list<Lit> & vals){

		tmp2.clear();
		for (Lit l:vals)
			tmp2.push(l);
		return Or(tmp2);

	}
	Lit Or(const vec<Lit> & vals){
		return Or_(vals,lit_Undef);
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
		S.addClause(a,b,~out);
		S.addClause(~a,b,out);
		S.addClause(a,~b,out);
		S.addClause(~a,~b,~out);
		return out;
	}
	Lit Xor_(Lit a, Lit b, Lit out){
		if(isConst(a) || isConst(b)){
			if (isConstTrue(a)){
				if(out!=lit_Undef){
					AssertEqual(~b,out);
				}
				return ~b;
			}else if (isConstFalse(a)){
				if(out!=lit_Undef){
					AssertEqual(b,out);
				}
				return b;
			}else if (isConstTrue(b)){
				if(out!=lit_Undef){
					AssertEqual(~a,out);
				}
				return ~a;
			}else if (isConstFalse(b)){
				if(out!=lit_Undef){
					AssertEqual(a,out);
				}
				return a;
			}
		}
		//return Or(And(a, ~b), And(~a,b));
		if(out!=lit_Undef){
			out = mkLit(S.newVar());
		}
		S.addClause(a,b,~out);
		S.addClause(~a,b,out);
		S.addClause(a,~b,out);
		S.addClause(~a,~b,~out);
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
	Lit Xnor_(Lit a, Lit b, Lit out){
		return ~Xnor(a,b,out);
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

	Lit Implies(Lit a, Lit b){
		return Or(~a, b);
	}

	Lit Ite(Lit cond, Lit thn, Lit els){
		Lit l = ~And(cond,~thn);
		Lit r = ~And(~cond,~els);
		return And(l,r);
	}
	Lit Ite_(Lit cond, Lit thn, Lit els,Lit out){
		Lit l = ~And(cond,~thn);
		Lit r = ~And(~cond,~els);
		return And_(l,r,out);
	}

	Lit HalfAdder(Lit a, Lit b, Lit & carry_out){
		carry_out = And(a,b);
		return Xor(a,b);
	}


	Lit HalfAdder_(Lit a, Lit b, Lit & carry_out, Lit & out){
		carry_out = And_(a,b,carry_out);
		return Xor_(a,b,out);
	}

	void Add(vec<Lit> & a, vec<Lit> & b,vec<Lit> & store_out, Lit & carry_out){
		assert(b.size()==a.size());
		store_out.clear();
		Lit carry=lit_Undef;
		Lit sum = HalfAdder(a[0],b[0],carry);
		store_out.push(sum);
		for(int i = 1;i<a.size();i++){
			Lit x1 = Xor(a[i],b[i]);
			store_out.push(Xor(x1, carry));
			Lit c1 = And(x1,carry);
			carry = Or(And(a[i],b[i]),c1);
		}
		carry_out= carry;
	}
	void Add_(vec<Lit> & a, vec<Lit> & b,vec<Lit> & store_out, Lit & carry_out){
		assert(b.size()==a.size());
		assert(store_out.size()==a.size());
		Lit carry=lit_Undef;
		Lit sum = HalfAdder_(a[0],b[0],carry,store_out[0]);
		for(int i = 1;i<a.size();i++){
			Lit x1 = Xor(a[i],b[i]);
			Xor_(x1, carry,store_out[i]);
			Lit c1 = And(x1,carry);
			carry = Or(And(a[i],b[i]),c1);
		}
		if(carry_out!=lit_Undef){
			AssertEqual(carry_out,carry);
		}else {
			carry_out = carry;
		}
	}

	void Assert(Lit l){
		S.addClause(l);
	}

	void AssertOr(Lit a)
	{
		S.addClause(a);
	}

	void AssertOr(Lit a, Lit b){
		S.addClause(a,b);
	}

	void AssertOr(const std::list<Lit> & vals){

		tmp.clear();
		for (Lit l:vals)
			tmp.push(l);
		S.addClause(tmp);

	}

	void AssertOr(const vec<Lit> & vals){
		tmp.clear();
		for (Lit l:vals)
			tmp.push(l);
		S.addClause(tmp);
	}

	template<typename... Args>
	void AssertOr(Lit a, Lit b, Args... args)
	{
		store.clear();
		store.push(a);
		collect(store,b,args...);
		S.addClause(store);
	}

	void AssertNand(Lit a)
	{
		S.addClause(~a);
	}

	void AssertNand(Lit a, Lit b){
		S.addClause(~a,~b);
	}

	void AssertNand(const std::list<Lit> & vals){

		tmp.clear();
		for (Lit l:vals)
			tmp.push(~l);
		S.addClause(tmp);

	}

	void AssertNand(const vec<Lit> & vals){
		tmp.clear();
		for (Lit l:vals)
			tmp.push(~l);
		S.addClause(tmp);
	}

	template<typename... Args>
	void AssertNand(Lit a, Lit b, Args... args)
	{
		store.clear();
		store.push(a);
		collect(store,b,args...);

		tmp.clear();
		for (Lit l:store)
			tmp.push(~l);
		S.addClause(tmp);
	}

	void AssertAnd(Lit a)
	{
		S.addClause(a);
	}

	void AssertAnd(Lit a, Lit b){
		S.addClause(a);
		S.addClause(b);
	}

	void AssertAnd(const std::list<Lit> & vals){
		for (Lit l:vals)
			S.addClause(l);
	}

	void AssertAnd(const vec<Lit> & vals){
		for (Lit l:vals)
			S.addClause(l);
	}

	template<typename... Args>
	void AssertAnd(Lit a, Lit b, Args... args)
	{
		store.clear();
		store.push(a);
		collect(store,b,args...);
		for (Lit l:store)
			S.addClause(l);
	}

	void AssertNor(Lit a)
	{
		S.addClause(~a);
	}

	void AssertNor(Lit a, Lit b){
		S.addClause(~a);
		S.addClause(~b);
	}

	void AssertNor(const std::list<Lit> & vals){
		for (Lit l:vals)
			S.addClause(~l);
	}

	void AssertNor(const vec<Lit> & vals){
		for (Lit l:vals)
			S.addClause(~l);
	}

	template<typename... Args>
	void AssertNor(Lit a, Lit b, Args... args)
	{
		store.clear();
		store.push(a);
		collect(store,b,args...);
		for (Lit l:store)
			S.addClause(~l);
	}

	void AssertXor(Lit a)
	{
		S.addClause(a);//is this correct?
	}

	void AssertXor(Lit a, Lit b){
		S.addClause(a, b);
		S.addClause(~a, ~b);
	}

	void AssertXor(const std::list<Lit> & vals){
		tmp.clear();
		for(Lit l:vals)
			tmp.push(l);
		AssertXor(tmp);
	}

	void AssertXor(const vec<Lit> & vals){
		if(vals.size()==0)
			Assert(getFalse());
		else if (vals.size()==1){
			Assert(vals[0]);
		}else if (vals.size()==2){
			S.addClause(vals[0], vals[1]);
			S.addClause(~vals[0], ~vals[1]);
		}else{
			//this can probably be done better...
			Assert(Xor(vals));
		}
	}

	template<typename... Args>
	void AssertXor(Lit a, Lit b, Args... args)
	{
		store.clear();
		store.push(a);
		collect(store,b,args...);
		AssertXor(store);
	}

	void AssertXnor(Lit a)
	{
		S.addClause(a);//is this correct?
	}

	void AssertXnor(Lit a, Lit b){
		S.addClause(~a, b);
		S.addClause(a, ~b);
	}

	void AssertXnor(const std::list<Lit> & vals){
		tmp.clear();
		for(Lit l:vals)
			tmp.push(l);
		AssertXnor(tmp);
	}

	void AssertXnor(const vec<Lit> & vals){
		if(vals.size()==0)
			Assert(getFalse());
		else if (vals.size()==1){
			Assert(vals[0]);
		}else if (vals.size()==2){
			S.addClause(~vals[0], vals[1]);
			S.addClause(vals[0], ~vals[1]);
		}else{
			//this can probably be done better...
			Assert(Xnor(vals));
		}
	}

	template<typename... Args>
	void AssertXnor(Lit a, Lit b, Args... args)
	{
		store.clear();
		store.push(a);
		collect(store,b,args...);
		AssertXnor(store);
	}

	void AssertImplies(Lit a, Lit b){
		S.addClause(~a, b);
	}

	void AssertEqual(Lit a, Lit b){
		AssertXnor(a,b);
	}

	void AssertEqual(const std::list<Lit> & vals){
		tmp.clear();
		for(Lit l:vals)
			tmp.push(l);
		AssertXnor(tmp);
	}

	void AssertEqual(const vec<Lit> & vals){
		if(vals.size()==0)
			Assert(getFalse());
		else if (vals.size()==1){
			Assert(vals[0]);
		}else if (vals.size()==2){
			S.addClause(~vals[0], vals[1]);
			S.addClause(vals[0], ~vals[1]);
		}else{
			//this can probably be done better...
			Assert(Xnor(vals));
		}
	}

	template<typename... Args>
	void AssertEqual(Lit a, Lit b, Args... args)
	{
		store.clear();
		store.push(a);
		collect(store,b,args...);
		AssertXnor(store);
	}
};
};

#endif /* CIRCUIT_H_ */
