/**************************************************************************************************
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
#include <stdio.h>

namespace Monosat {

//Helper methods for expressing combinatorial logic in CNF.
//Does not perform any kind of circuit rewriting, simplification, or memoization.
template<class Solver>
class Circuit {
    Solver& S;
    Lit lit_True = lit_Undef;

    bool isConst(Lit l){
        return isConstTrue(l) || isConstFalse(l);
    }

    bool isConstTrue(Lit l){
        return S.value(l) == l_True && S.level(var(l)) == 0;
    }

    bool isConstFalse(Lit l){
        return S.value(l) == l_False && S.level(var(l)) == 0;
    }

    vec<Lit> tmp;
    vec<Lit> tmp2;
    vec<Lit> clause;
    vec<Lit> store;
    FILE* outfile = nullptr;

    inline int dimacs(Solver& S, Lit internalLit){
        Lit l = S.unmap(internalLit);
        return sign(l) ? -(var(l) + 1) : (var(l) + 1);
    }

    bool _addClause(Lit a){
        if(outfile){
            fprintf(outfile, "%d 0\n ", dimacs(S, a));
            fflush(outfile);
        }
        return S.addClause(a);
    }

    bool _addClause(Lit a, Lit b){
        if(outfile){
            fprintf(outfile, "%d %d 0\n ", dimacs(S, a), dimacs(S, b));
            fflush(outfile);
        }
        return S.addClause(a, b);
    }

    bool _addClause(Lit a, Lit b, Lit c){
        if(outfile){
            fprintf(outfile, "%d %d %d 0\n ", dimacs(S, a), dimacs(S, b), dimacs(S, c));
            fflush(outfile);
        }
        return S.addClause(a, b, c);
    }

    bool _addClause(vec<Lit>& clause){
        if(outfile){
            for(Lit l:clause){
                fprintf(outfile, "%d ", dimacs(S, l));
            }
            fprintf(outfile, "0\n");
            fflush(outfile);
        }
        return S.addClause(clause);
    }

    template<typename... Args>
    void collect(vec<Lit>& store, Lit a, Args... args){
        store.push(a);
        collect(store, args...);
    }

    void collect(vec<Lit>& store, Lit a){
        store.push(a);
    }

    //Note: a vector of size zero will always return lit_True
    Lit bin_op(vec<Lit>& store, Lit (Circuit::*f)(Lit, Lit)){
        int n = store.size();
        while(n > 1){
            int p = 0;
            for(int i = 0; i < n; i += 2){
                if(i + 1 < n){
                    store[p++] = (this->**&f)(store[i], store[i + 1]);
                }else{
                    store[p++] = store[i];
                }
            }
            n = p;
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
    /**
     * Specify a file to write constructed clauses to (in GNF format)
     * @param file
     */
    void setOutputFile(FILE* file){
        this->outfile = file;
    }

    Circuit(Solver& S) : S(S){
        lit_True = S.True();
        //_addClause(lit_True);
    }

    Lit newLit(bool decisionLit = true){
        return mkLit(S.newVar(false, decisionLit));
    }

    Solver& getSolver(){
        return S;
    }

    Lit getTrue(){
        return lit_True;
    }

    Lit getFalse(){
        return ~lit_True;
    }

    Lit And(Lit a){
        return a;
    }

    Lit And_(Lit a, Lit b, Lit out){
        //special case these
        if(isConst(a) || isConst(b)){
            if(isConstTrue(a) && isConstTrue(b)){
                if(out != lit_Undef){
                    Assert(out);
                    return out;
                }
                return getTrue();
            }else if(isConstFalse(a) || isConstFalse(b)){
                if(out != lit_Undef){
                    Assert(~out);
                    return out;
                }
                return getFalse();
            }

            if(isConstFalse(a)){
                if(out != lit_Undef){
                    Assert(~out);
                    return out;
                }
                return getFalse();
            }else if(isConstFalse(b)){
                if(out != lit_Undef){
                    Assert(~out);
                    return out;
                }
                return getFalse();
            }else if(isConstTrue(a)){
                if(out != lit_Undef){
                    AssertEqual(b, out);
                    return out;
                }
                return b;
            }else if(isConstTrue(b)){
                if(out != lit_Undef){
                    AssertEqual(a, out);
                    return out;
                }
                return a;
            }
        }
        if(out == lit_Undef){
            out = mkLit(S.newVar());
        }
        _addClause(a, ~out);
        _addClause(b, ~out);
        _addClause(~a, ~b, out);
        return out;
    }

    Lit And_(const std::list<Lit>& vals, Lit out){

        tmp2.clear();
        for(Lit l:vals)
            tmp2.push(l);
        return And_(tmp2, out);

    }

    Lit And_(const vec<Lit>& vals, Lit out){
        tmp.clear();
        for(Lit l:vals){
            if(isConstFalse(l)){
                if(out != lit_Undef){
                    Assert(~out);
                }
                return getFalse();
            }else if(isConstTrue(l)){
                //leave literal out
            }else{
                tmp.push(l);
            }
        }
        //all arguments are constant true
        if(tmp.size() == 0){
            if(out != lit_Undef){
                Assert(out);
            }
            return getTrue();
        }else if(tmp.size() == 1){
            if(out != lit_Undef){
                AssertEqual(tmp[0], out);
                return out;
            }
            return tmp[0];
        }
        if(out == lit_Undef){
            out = mkLit(S.newVar());
        }
        for(Lit l:tmp){
            _addClause(l, ~out);
        }
        for(int i = 0; i < tmp.size(); i++){
            tmp[i] = ~tmp[i];
        }
        tmp.push(out);
        _addClause(tmp);
        tmp.clear();
        return out;
    }

    void AssertImpliesAnd_(Lit implies, const vec<Lit>& vals, Lit out = lit_Undef){
        tmp.clear();
        for(Lit l:vals){
            if(isConstFalse(l)){
                if(out != lit_Undef){
                    AssertImplies(implies, ~out);
                }else{
                    Assert(~implies);
                }
                return;
            }else if(isConstTrue(l)){
                //leave literal out
            }else{
                tmp.push(l);
            }
        }

        if(tmp.size() == 0){
            if(out == lit_Undef){
                //do nothing
            }else{
                Assert(out);//out holds always
            }
        }else{
            if(out == lit_Undef){
                for(Lit l:tmp){
                    _addClause(l, ~implies);
                }
            }else{
                for(Lit l:tmp){
                    _addClause(l, ~out, ~implies);
                }
            }

            for(int i = 0; i < tmp.size(); i++){
                tmp[i] = ~tmp[i];
            }
            if(out == lit_Undef){
                tmp.push(out);
            }
            tmp.push(~implies);
            _addClause(tmp);
            tmp.clear();
        }
    }

    Lit And(Lit a, Lit b){
        if(a == b){
            return a;
        }else if(var(a) == var(b)){
            assert(sign(a) != sign(b));
            return getFalse();
        }
        //special case these
        if(isConst(a) || isConst(b)){
            if(isConstTrue(a)){
                return b;
            }else if(isConstFalse(a)){
                return a;
            }else if(isConstTrue(b)){
                return a;
            }else if(isConstFalse(b)){
                return b;
            }
        }

        Lit out = mkLit(S.newVar());
        _addClause(a, ~out);
        _addClause(b, ~out);
        _addClause(~a, ~b, out);
        return out;
    }

    Lit And(const std::list<Lit>& vals){
        tmp2.clear();
        for(Lit l:vals)
            tmp2.push(l);
        return And(tmp2);

    }

    Lit And(const vec<Lit>& vals){
        return And_(vals, lit_Undef);
    }

    template<typename... Args>
    Lit And(Lit a, Lit b, Args... args){
        store.clear();
        store.push(a);
        collect(store, b, args...);

        return And(store);
    }


    Lit Or(Lit a){
        return a;
    }

    Lit Or(Lit a, Lit b){
        if(a == b){
            return a;
        }else if(var(a) == var(b)){
            assert(sign(a) != sign(b));
            return getTrue();
        }
        if(isConstTrue(a)){
            return getTrue();
        }else if(isConstTrue(b)){
            return getTrue();
        }else if(isConstFalse(a)){
            return b;
        }else if(isConstFalse(b)){
            return a;
        }

        Lit out = mkLit(S.newVar());
        _addClause(~a, out);
        _addClause(~b, out);
        _addClause(a, b, ~out);
        return out;
    }

    Lit Or_(Lit a, Lit b, Lit out){

        //special case these
        if(isConst(a) || isConst(b)){

            if(isConstFalse(a) && isConstFalse(b)){
                if(out != lit_Undef){
                    Assert(~out);
                }
                return getFalse();
            }else if(isConstTrue(a) || isConstTrue(b)){
                if(out != lit_Undef){
                    Assert(out);
                }
                return getTrue();
            }

            if(isConstTrue(a)){
                if(out != lit_Undef){
                    Assert(out);
                }
                return getTrue();
            }else if(isConstTrue(b)){
                if(out != lit_Undef){
                    Assert(out);
                }
                return getTrue();
            }else if(isConstFalse(a)){
                if(out != lit_Undef){
                    AssertEqual(b, out);
                }
                return b;
            }else if(isConstFalse(b)){
                if(out != lit_Undef){
                    AssertEqual(a, out);
                }
                return a;
            }
        }
        if(out == lit_Undef){
            out = mkLit(S.newVar());
        }
        _addClause(~a, out);
        _addClause(~b, out);
        _addClause(a, b, ~out);
        return out;
    }

    Lit Or_(const vec<Lit>& vals, Lit out){
        tmp.clear();
        for(Lit l:vals){
            if(isConstFalse(l)){
                //leave literal out
            }else if(isConstTrue(l)){
                if(out != lit_Undef){
                    Assert(out);
                }
                return getTrue();
            }else{
                tmp.push(l);
            }
        }
        //all arguments are constant false
        if(tmp.size() == 0){
            if(out != lit_Undef){
                Assert(~out);
            }
            return getFalse();//or should this be true?
        }else if(tmp.size() == 1){
            if(out != lit_Undef){
                AssertEqual(tmp[0], out);
                return out;
            }
            return tmp[0];
        }
        if(out == lit_Undef){
            out = mkLit(S.newVar());
        }
        for(Lit l:tmp){
            _addClause(~l, out);
        }

        tmp.push(~out);
        _addClause(tmp);
        tmp.clear();
        return out;
    }

    //If this gate is true, then all of vals must be true.
    //But if this gate is false, vals may be true or false.
    Lit ImpliesAnd(const vec<Lit>& vals, Lit out = lit_Undef){
        if(out == lit_Undef){
            out = mkLit(S.newVar());
        }
        for(Lit l:vals){
            AssertImplies(out, l);
        }
        return out;
    }

    //If this gate is true, then at least one of vals must be true.
    //But if this gate is false, vals may be true or false.
    Lit ImpliesOr(const vec<Lit>& vals, Lit out = lit_Undef){
        if(out == lit_Undef){
            out = mkLit(S.newVar());
        }
        AssertImpliesOr_(out, vals, lit_Undef);
        return out;
    }

    //This is an OR condition that holds only if implies is true
    void AssertImpliesOr_(Lit implies, const vec<Lit>& vals, Lit out = lit_Undef){
        tmp.clear();
        for(Lit l:vals){
            if(isConstFalse(l)){
                //leave literal out
            }else if(isConstTrue(l)){
                if(out != lit_Undef){
                    AssertImplies(implies, out);
                }
                //no effect on implies, otherwise
                return;
            }else{
                tmp.push(l);
            }
        }

        if(vals.size() == 0){
            //all literals in vals are false (or vals is empty)
            if(out != lit_Undef){
                AssertImplies(out, ~implies);
            }else{
                Assert(~implies);
            }
            return;
        }
        if(out != lit_Undef){
            //is this correct?
            tmp.push(~out);
            tmp.push(
                    ~implies);//if implies is true, then at least one of the elements of tmp must be true, or out must be false.
            _addClause(tmp);
            tmp.clear();
        }else{
            tmp.push(~implies);//if implies is true, then at least one of the elements of tmp must be true.
            _addClause(tmp);
            tmp.clear();
        }


    }

    Lit Or(const std::list<Lit>& vals){

        tmp2.clear();
        for(Lit l:vals)
            tmp2.push(l);
        return Or(tmp2);

    }

    Lit Or(const vec<Lit>& vals){
        return Or_(vals, lit_Undef);
    }

    template<typename... Args>
    Lit Or(Lit a, Lit b, Args... args){
        store.clear();
        store.push(a);
        collect(store, b, args...);

        return Or(store);
    }


    Lit Nor(Lit a){
        return ~a;
    }

    Lit Nor(Lit a, Lit b){
        return ~Or(a, b);
    }

    Lit Nor(const std::list<Lit>& vals){
        return ~Or(vals);
    }

    Lit Nor(const vec<Lit>& vals){
        return ~Or(vals);
    }

    template<typename... Args>
    Lit Nor(Lit a, Lit b, Args... args){
        store.clear();
        store.push(a);
        collect(store, b, args...);

        return Nor(store);
    }


    Lit Nand(Lit a){
        return ~a;
    }

    Lit Nand(Lit a, Lit b){
        return ~And(a, b);
    }

    Lit Nand(const std::list<Lit>& vals){
        return ~And(vals);

    }

    Lit Nand(const vec<Lit>& vals){
        return ~And(vals);
    }

    template<typename... Args>
    Lit Nand(Lit a, Lit b, Args... args){
        store.clear();
        store.push(a);
        collect(store, b, args...);
        return Nand(store);
    }

    Lit Xor(Lit a){
        return ~a;
    }

    Lit Xor(Lit a, Lit b){
        if(isConst(a) || isConst(b)){
            if(isConstTrue(a)){
                return ~b;
            }else if(isConstFalse(a)){
                return b;
            }else if(isConstTrue(b)){
                return ~a;
            }else if(isConstFalse(b)){
                return a;
            }
        }
        //return Or(And(a, ~b), And(~a,b));
        Lit out = mkLit(S.newVar());
        _addClause(a, b, ~out);
        _addClause(~a, b, out);
        _addClause(a, ~b, out);
        _addClause(~a, ~b, ~out);
        return out;
    }

    Lit Xor_(Lit a, Lit b, Lit out){
        if(isConst(a) || isConst(b)){
            if(isConstTrue(a)){
                if(out != lit_Undef){
                    AssertEqual(~b, out);
                }
                return ~b;
            }else if(isConstFalse(a)){
                if(out != lit_Undef){
                    AssertEqual(b, out);
                }
                return b;
            }else if(isConstTrue(b)){
                if(out != lit_Undef){
                    AssertEqual(~a, out);
                }
                return ~a;
            }else if(isConstFalse(b)){
                if(out != lit_Undef){
                    AssertEqual(a, out);
                }
                return a;
            }
        }
        //return Or(And(a, ~b), And(~a,b));
        if(out != lit_Undef){
            out = mkLit(S.newVar());
        }
        _addClause(a, b, ~out);
        _addClause(~a, b, out);
        _addClause(a, ~b, out);
        _addClause(~a, ~b, ~out);
        return out;
    }

    Lit Xor(const std::list<Lit>& vals){
        tmp.clear();
        for(Lit l:vals)
            tmp.push(l);
        return bin_op(tmp, &Circuit::Xor);

    }

    Lit Xor(const vec<Lit>& vals){
        tmp.clear();
        for(Lit l:vals)
            tmp.push(l);
        return bin_op(tmp, &Circuit::Xor);
    }

    template<typename... Args>
    Lit Xor(Lit a, Lit b, Args... args){
        store.clear();
        store.push(a);
        collect(store, b, args...);
        return Xor(store);
    }


    Lit Xnor(Lit a){
        return a;
    }

    Lit Xnor(Lit a, Lit b){
        return ~Xor(a, b);
    }

    Lit Xnor_(Lit a, Lit b, Lit out){
        return ~Xnor(a, b, out);
    }

    Lit Xnor(const std::list<Lit>& vals){
        return ~Xor(vals);

    }

    Lit Xnor(const vec<Lit>& vals){
        return ~Xor(vals);
    }

    template<typename... Args>
    Lit Xnor(Lit a, Lit b, Args... args){
        store.clear();
        store.push(a);
        collect(store, b, args...);
        return Xnor(store);
    }

    Lit Implies(Lit a, Lit b){
        return Or(~a, b);
    }


    Lit Implies_(Lit a, Lit b, Lit out){
        return Or_(~a, b, out);
    }


    Lit Ite(Lit cond, Lit thn, Lit els){
        Lit l = ~And(cond, ~thn);
        Lit r = ~And(~cond, ~els);
        return And(l, r);
    }

    Lit Ite_(Lit cond, Lit thn, Lit els, Lit out){
        Lit l = ~And(cond, ~thn);
        Lit r = ~And(~cond, ~els);
        return And_(l, r, out);
    }

    Lit HalfAdder(Lit a, Lit b, Lit& carry_out){
        carry_out = And(a, b);
        return Xor(a, b);
    }

    Lit FullAdder(Lit a, Lit b, Lit c_in, Lit& carry_out){
        Lit a_xor_b = Xor(a, b);
        Lit a_and_b = And(a, b);
        carry_out = Or(a_and_b, And(c_in, a_xor_b));
        return Xor(a_xor_b, c_in);

    }

    Lit HalfAdder_(Lit a, Lit b, Lit& carry_out, Lit& out){
        carry_out = And_(a, b, carry_out);
        return Xor_(a, b, out);
    }

    //a-b
    Lit HalfSubtractor(Lit a, Lit b, Lit& borrow_out){
        borrow_out = And(~a, b);
        return Xor(a, b);
    }

    //a-b
    Lit HalfSubtractor_(Lit a, Lit b, Lit& borrow_out, Lit& out){
        borrow_out = And_(~a, b, borrow_out);
        return Xor_(a, b, out);
    }


    void Add(vec<Lit>& a, vec<Lit>& b, vec<Lit>& store_out, Lit& carry_out){
        assert(b.size() == a.size());
        store_out.clear();
        Lit carry = lit_Undef;
        Lit sum = HalfAdder(a[0], b[0], carry);
        store_out.push(sum);
        for(int i = 1; i < a.size(); i++){
            Lit x1 = Xor(a[i], b[i]);
            store_out.push(Xor(x1, carry));
            Lit c1 = And(x1, carry);
            carry = Or(And(a[i], b[i]), c1);
        }
        carry_out = carry;
    }

    void Add_(vec<Lit>& a, vec<Lit>& b, vec<Lit>& store_out, Lit& carry_out){
        assert(b.size() == a.size());
        assert(store_out.size() == a.size());
        Lit carry = lit_Undef;
        Lit sum = HalfAdder_(a[0], b[0], carry, store_out[0]);
        for(int i = 1; i < a.size(); i++){
            Lit x1 = Xor(a[i], b[i]);
            Xor_(x1, carry, store_out[i]);
            Lit c1 = And(x1, carry);
            carry = Or(And(a[i], b[i]), c1);
        }
        if(carry_out != lit_Undef){
            AssertEqual(carry_out, carry);
        }else{
            carry_out = carry;
        }
    }

    void Subtract(vec<Lit>& a, vec<Lit>& b, vec<Lit>& store_out, Lit& borrow_out){
        assert(b.size() == a.size());
        store_out.clear();
        Lit borrow = lit_Undef;
        Lit dif = HalfSubtractor(a[0], b[0], borrow);
        store_out.push(dif);
        for(int i = 1; i < a.size(); i++){
            Lit x1 = Xor(a[i], b[i]);
            store_out.push(Xor(x1, borrow));
            Lit c1 = And(~x1, borrow);
            borrow = Or(And(~a[i], b[i]), c1); //check this!
        }
        borrow_out = borrow;
    }

    void Subtract_(vec<Lit>& a, vec<Lit>& b, vec<Lit>& store_out, Lit& borrow_out){
        assert(b.size() == a.size());
        assert(store_out.size() == a.size());
        Lit borrow = lit_Undef;
        Lit dif = HalfSubtractor_(a[0], b[0], borrow, store_out[0]);
        for(int i = 1; i < a.size(); i++){
            Lit x1 = Xor(a[i], b[i]);
            Xor_(x1, borrow, store_out[i]);
            Lit c1 = And(~x1, borrow);
            borrow = Or(And(~a[i], b[i]), c1); //check this!
        }
        if(borrow_out != lit_Undef){
            AssertEqual(borrow_out, borrow);
        }else{
            borrow_out = borrow;
        }
    }

    //perform two's complement negation (invert bits and add 1)
    void Negate(vec<Lit>& a, vec<Lit>& store_out){
        store_out.clear();
        Lit carry = lit_Undef;
        Lit sum = HalfAdder(~a[0], getTrue(), carry);
        store_out.push(sum);
        for(int i = 1; i < a.size(); i++){
            Lit x1 = Xor(~a[i], carry);
            store_out.push(x1);
            carry = And(~a[i], carry);

        }
    }

    //perform two's complement negation (invert bits and add 1)
    void Negate_(vec<Lit>& a, vec<Lit>& store_out){
        store_out.clear();
        assert(store_out.size() == a.size());
        Lit carry = lit_Undef;
        Lit sum = HalfAdder(~a[0], getTrue(), carry);
        store_out.push(sum);
        for(int i = 1; i < a.size(); i++){
            Xor_(~a[i], carry, store_out[i]);
            carry = And(~a[i], carry);
        }
    }

protected:
    Lit mBlock(Lit m0, Lit m1, Lit q0, Lit q1, Lit c_in, Lit& c_out){
        Lit a = And(q1, m0);
        Lit b = And(q0, m1);
        Lit res = FullAdder(a, b, c_in, c_out);
        return res;
    }

    Lit nBlock(Lit m0, Lit m1, Lit q0, Lit c_in, Lit& c_out){
        Lit a = And(q0, m0);
        Lit res = FullAdder(a, m1, c_in, c_out);
        return res;
    }

public:
    void Multiply(vec<Lit>& m, vec<Lit>& q, vec<Lit>& store_out){
        assert(m.size() == q.size());
        store_out.clear();

        Lit res = And(m[0], q[0]);
        store_out.push(res);

        if(m.size() == 1){
            return;
        }

        vec<Lit> cur_results;
        vec<Lit> cur_passthrough;

        //first row
        Lit prev_c_out = getFalse();
        for(int i = 1; i <= m.size(); i++){
            Lit m0 = m[i - 1];
            Lit m1 = i < m.size() ? m[i] : getFalse();
            Lit q0 = q[0];
            Lit q1 = q[1];
            Lit c_in = prev_c_out;
            Lit c_out;
            Lit res = mBlock(m0, m1, q0, q1, c_in, c_out);
            prev_c_out = c_out;
            cur_passthrough.push(m0);
            cur_results.push(res);
            if(i == 1){
                store_out.push(res);
            }
        }
        cur_results.push(prev_c_out);
        vec<Lit> prev_results;
        vec<Lit> prev_passthrough;
        //remaining rows
        for(int j = 2; j < q.size(); j++){
            assert(cur_passthrough.size() == cur_results.size() - 1);
            cur_results.copyTo(prev_results);
            cur_results.clear();
            cur_passthrough.copyTo(prev_passthrough);
            cur_passthrough.clear();

            Lit q_in = q[j];
            prev_c_out = getFalse();
            assert(prev_results.size() == m.size() + 1);
            for(int i = 1; i <= m.size(); i++){
                Lit m0 = prev_passthrough[i - 1];
                Lit m1 = prev_results[i];
                Lit c_in = prev_c_out;
                Lit c_out;
                Lit res = nBlock(m0, m1, q_in, c_in, c_out);
                prev_c_out = c_out;
                cur_passthrough.push(m0);
                cur_results.push(res);
                if(i == 1){
                    store_out.push(res);
                }
            }
            cur_results.push(prev_c_out);
        }
        //last row
        assert(cur_results[0] == store_out.last());
        for(int i = 1; i < cur_results.size(); i++){
            store_out.push(cur_results[i]);
        }
        assert(store_out.last() == prev_c_out);
        //store_out.push(prev_c_out);
        assert(store_out.size() == (m.size() + q.size()));
    }

    void Assert(Lit l){
        _addClause(l);
    }

    void AssertOr(Lit a){
        _addClause(a);
    }

    void AssertOr(Lit a, Lit b){
        _addClause(a, b);
    }

    void AssertOr(Lit a, Lit b, Lit c){
        _addClause(a, b, c);
    }

    void AssertOr(const std::list<Lit>& vals){

        tmp.clear();
        for(Lit l:vals)
            tmp.push(l);
        _addClause(tmp);

    }

    void AssertOr(const vec<Lit>& vals){
        tmp.clear();
        for(Lit l:vals){
            assert(l != lit_Undef);
            tmp.push(l);
        }
        _addClause(tmp);
    }

    template<typename... Args>
    void AssertOr(Lit a, Lit b, Args... args){
        store.clear();
        store.push(a);
        collect(store, b, args...);
        _addClause(store);
    }

    void AssertNand(Lit a){
        _addClause(~a);
    }

    void AssertNand(Lit a, Lit b){
        _addClause(~a, ~b);
    }

    void AssertNand(const std::list<Lit>& vals){

        tmp.clear();
        for(Lit l:vals)
            tmp.push(~l);
        _addClause(tmp);

    }

    void AssertNand(const vec<Lit>& vals){
        tmp.clear();
        for(Lit l:vals)
            tmp.push(~l);
        _addClause(tmp);
    }

    template<typename... Args>
    void AssertNand(Lit a, Lit b, Args... args){
        store.clear();
        store.push(a);
        collect(store, b, args...);

        tmp.clear();
        for(Lit l:store)
            tmp.push(~l);
        _addClause(tmp);
    }

    void AssertAnd(Lit a){
        _addClause(a);
    }

    void AssertAnd(Lit a, Lit b){
        _addClause(a);
        _addClause(b);
    }

    void AssertAnd(const std::list<Lit>& vals){
        for(Lit l:vals)
            _addClause(l);
    }

    void AssertAnd(const vec<Lit>& vals){
        for(Lit l:vals)
            _addClause(l);
    }

    template<typename... Args>
    void AssertAnd(Lit a, Lit b, Args... args){
        store.clear();
        store.push(a);
        collect(store, b, args...);
        for(Lit l:store)
            _addClause(l);
    }

    void AssertNor(Lit a){
        _addClause(~a);
    }

    void AssertNor(Lit a, Lit b){
        _addClause(~a);
        _addClause(~b);
    }

    void AssertNor(const std::list<Lit>& vals){
        for(Lit l:vals)
            _addClause(~l);
    }

    void AssertNor(const vec<Lit>& vals){
        for(Lit l:vals)
            _addClause(~l);
    }

    template<typename... Args>
    void AssertNor(Lit a, Lit b, Args... args){
        store.clear();
        store.push(a);
        collect(store, b, args...);
        for(Lit l:store)
            _addClause(~l);
    }

    void AssertXor(Lit a){
        _addClause(a);//is this correct?
    }

    void AssertXor(Lit a, Lit b){
        _addClause(a, b);
        _addClause(~a, ~b);
    }

    void AssertXor(const std::list<Lit>& vals){
        tmp.clear();
        for(Lit l:vals)
            tmp.push(l);
        AssertXor(tmp);
    }

    void AssertXor(const vec<Lit>& vals){
        if(vals.size() == 0)
            Assert(getFalse());
        else if(vals.size() == 1){
            Assert(vals[0]);
        }else if(vals.size() == 2){
            _addClause(vals[0], vals[1]);
            _addClause(~vals[0], ~vals[1]);
        }else{
            //this can probably be done better...
            Assert(Xor(vals));
        }
    }

    template<typename... Args>
    void AssertXor(Lit a, Lit b, Args... args){
        store.clear();
        store.push(a);
        collect(store, b, args...);
        AssertXor(store);
    }

    void AssertXnor(Lit a){
        _addClause(a);//is this correct?
    }

    void AssertXnor(Lit a, Lit b){
        _addClause(~a, b);
        _addClause(a, ~b);
    }

    void AssertXnor(const std::list<Lit>& vals){
        tmp.clear();
        for(Lit l:vals)
            tmp.push(l);
        AssertXnor(tmp);
    }

    void AssertXnor(const vec<Lit>& vals){
        if(vals.size() == 0)
            Assert(getFalse());
        else if(vals.size() == 1){
            Assert(vals[0]);
        }else if(vals.size() == 2){
            _addClause(~vals[0], vals[1]);
            _addClause(vals[0], ~vals[1]);
        }else{
            //this can probably be done better...
            Assert(Xnor(vals));
        }
    }

    template<typename... Args>
    void AssertXnor(Lit a, Lit b, Args... args){
        store.clear();
        store.push(a);
        collect(store, b, args...);
        AssertXnor(store);
    }

    void AssertImplies(Lit a, Lit b){
        _addClause(~a, b);
    }

    void AssertEqual(Lit a, Lit b){
        AssertXnor(a, b);
    }
/*
	void AssertEqual(const std::list<Lit> & vals){
		tmp.clear();
		for(Lit l:vals)
			tmp.push(l);
		//this is not a good definition.
		//equal should be true if all arguments are TRUE, or all arguments are FALSE.
		AssertXnor(tmp);
	}*/
/*
	void AssertEqual(const vec<Lit> & vals){
		if(vals.size()==0)
			Assert(getFalse());
		else if (vals.size()==1){
			Assert(vals[0]);
		}else if (vals.size()==2){
			_addClause(~vals[0], vals[1]);
			_addClause(vals[0], ~vals[1]);
		}else{
			//this is not a good definition.
		//equal should be true if all arguments are TRUE, or all arguments are FALSE.
			Assert(Xnor(vals));
		}
	}*/
/*
	template<typename... Args>
	void AssertEqual(Lit a, Lit b, Args... args)
	{
		store.clear();
		store.push(a);
		collect(store,b,args...);
		//this is not a good definition.
		//equal should be true if all arguments are TRUE, or all arguments are FALSE.
		AssertXnor(store);
	}*/

    int pow2roundup(uint32_t v){
        if(v == 0)
            return 0;
        uint32_t r = 0;
        while(v >>= 1) // unroll for more speed...
        {
            r++;
        }
        return r;
    }

public:
//from http://www.graphics.stanford.edu/~seander/bithacks.html
    bool _isPower2(int n){
        return (n > 0 && ((n & (n - 1)) == 0));
    }

    vec<Lit> assert_tmp;
    vec<Lit> assert_tmp2;


    Lit Equal(Lit a, Lit b){
        return Xnor(a, b);
    }

    Lit Equal(vec<Lit>& A, vec<Lit>& B){


        B.copyTo(assert_tmp);
        while(assert_tmp.size() < A.size())
            assert_tmp.push(getFalse());
        vec<Lit> xnors;
        //Lit all_eq=getTrue();
        for(int i = 0; i < A.size(); i++){
            xnors.push(Xnor(A[i], assert_tmp[i]));
        }

        /* for(int i = 0;i<A.size();i++){
             if(i==0){
                 all_eq =  Xnor(A[i], assert_tmp[i]);
             }else {
                 all_eq = And(all_eq, Xnor(A[i], assert_tmp[i]));
             }
         }*/
        // return all_eq;
        return And(xnors);
    }

    Lit LEQ(vec<Lit>& A, vec<Lit>& B){

        Lit p_lt_n = getFalse();
        Lit p_eq_n = getTrue();
        B.copyTo(assert_tmp);
        while(assert_tmp.size() < A.size())
            assert_tmp.push(getFalse());

        for(int j = A.size() - 1; j >= 0; j--){
            Lit pj_eq_n = Xnor(A[j], assert_tmp[j]);
            Lit pj_ltn = And(~A[j], assert_tmp[j]);
            Lit lt = And(p_eq_n, pj_ltn);
            p_lt_n = Or(lt, p_lt_n);

            p_eq_n = And(pj_eq_n, p_eq_n);

        }
        for(int j = A.size(); j < B.size(); j++){
            p_lt_n = Or(B[j], p_lt_n);
        }
        //The A must be less than or equal to B
        return Or(p_lt_n, p_eq_n);

    }

    Lit LT(vec<Lit>& A, vec<Lit>& B){

        Lit p_lt_n = getFalse();
        Lit p_eq_n = getTrue();
        B.copyTo(assert_tmp);
        while(assert_tmp.size() < A.size())
            assert_tmp.push(getFalse());

        for(int j = A.size() - 1; j >= 0; j--){
            Lit pj_eq_n = Xnor(A[j], assert_tmp[j]);
            Lit pj_ltn = And(~A[j], assert_tmp[j]);
            Lit lt = And(p_eq_n, pj_ltn);
            p_lt_n = Or(lt, p_lt_n);
            if(j > 0){
                p_eq_n = And(pj_eq_n, p_eq_n);
            }
        }
        for(int j = A.size(); j < B.size(); j++){
            p_lt_n = Or(B[j], p_lt_n);
        }
        //The A must be less than B
        return p_lt_n;
    }

    void AssertEqual(vec<Lit>& A, vec<Lit>& B){

        Lit p_lt_n = getFalse();
        Lit p_eq_n = getTrue();
        B.copyTo(assert_tmp);
        while(assert_tmp.size() < A.size())
            assert_tmp.push(getFalse());

        for(int i = 0; i < A.size(); i++){
            AssertEqual(A[i], B[i]);
        }
        for(int j = A.size(); j < B.size(); j++){
            Assert(~B[j]);//if B is larger than A, then all larger bits must be false.
        }
    }

    void AssertLEQ(const vec<Lit>& A, const vec<Lit>& B){

        Lit p_lt_n = getFalse();
        Lit p_eq_n = getTrue();
        B.copyTo(assert_tmp);
        while(assert_tmp.size() < A.size())
            assert_tmp.push(getFalse());

        for(int j = A.size() - 1; j >= 0; j--){
            Lit pj_eq_n = Xnor(A[j], assert_tmp[j]);
            Lit pj_ltn = And(~A[j], assert_tmp[j]);
            Lit lt = And(p_eq_n, pj_ltn);
            p_lt_n = Or(lt, p_lt_n);

            p_eq_n = And(pj_eq_n, p_eq_n);

        }
        //if B is bigger than A, and any of the larger B bits are true, then A is less than B
        for(int j = A.size(); j < B.size(); j++){
            p_lt_n = Or(B[j], p_lt_n);
        }
        //The A must be less than or equal to B
        _addClause(p_lt_n, p_eq_n);

    }

    void AssertLT(vec<Lit>& A, vec<Lit>& B){

        _addClause(LT(A, B));
    }

    Lit Mux(vec<Lit>& selector, vec<Lit>& data, int opt_mux = 0){

        vec<Lit> args;
        int n_args = pow2roundup(data.size() - 1) + 1;

        while(selector.size() < n_args){
            selector.push(mkLit(S.newVar()));
        }

        Lit out = mkLit(S.newVar());
        if((1 << n_args) > data.size()){
            //printf("Warning: non-power-of-2 number of operations incurs some small overhead\n");
        }
        if(opt_mux == 0 || opt_mux == 1){
            //build an explicit mux, without adding any new literals
            for(int i = 0; i < data.size(); i++){
                Lit in = data[i];
                args.clear();
                for(int j = 0; j < selector.size(); j++){
                    bool v = i & (1 << j);
                    if(v){
                        args.push(selector[j]);
                    }else{
                        args.push(~selector[j]);
                    }
                }
                for(int i = 0; i < args.size(); i++){
                    args[i] = ~args[i];
                }
                args.push(~in);
                args.push(out);
                _addClause(args);
                args.shrink(2);
                args.push(in);
                args.push(~out);
                _addClause(args);
            }
            if(opt_mux == 1){
                int max = 1 << selector.size();
                for(int i = data.size(); i < max; i++){

                    args.clear();
                    for(int j = 0; j < selector.size(); j++){
                        bool v = i & (1 << j);
                        if(v){
                            args.push(selector[j]);
                        }else{
                            args.push(~selector[j]);
                        }
                    }
                    for(int i = 0; i < args.size(); i++){
                        args[i] = ~args[i];
                    }
                    //explicitly rule out this assignment to args
                    _addClause(args);

                }
            }
        }else{
            //for opt_mux==2, consider using unary selectors, instead of binary.
            assert(false);
        }
        if(opt_mux == 0){
            //finally, if data.size() is not an exact power of two, assert that selector <= data.size()
            if(!_isPower2(data.size())){
                vec<Lit> max_val;
                for(int j = 0; j < selector.size(); j++){
                    bool v = data.size() & (1 << j);
                    if(v){
                        max_val.push(getTrue());
                    }else{
                        max_val.push(getFalse());
                    }
                }
                AssertLT(selector, max_val);
/*
			Lit p_lt_n=getFalse();
			Lit p_eq_n = getTrue();


			for(int j = selector.size()-1;j>=0;j--){
				Lit pj_eq_n = Xnor(selector[j],max_val[j]);
				Lit pj_ltn = And(~selector[j],max_val[j]);
				Lit lt = And(p_eq_n,pj_ltn);
				p_lt_n = Or(lt,p_lt_n);
				if(j>0) {
					p_eq_n = And(pj_eq_n, p_eq_n);
				}
			}
			//The p must be less than max_val
			_addClause(p_lt_n);*/
            }
        }
        return out;
    }

    //uses n^2 binary clauses to create a simple at-most-one constraint.
    //if you have more than 20 or so literals, strongly consider using a pseudo-Boolean constraint solver instead
    void AssertAMO(const vec<Lit>& lits){
        vec<Lit> set;
        S.cancelUntil(0);
        assert(S.decisionLevel() == 0);
        Lit constant_true_lit = lit_Undef;

        //first check if there are any literals which are constant true,
        //or if all literals are constant false.
        for(int i = 0; i < lits.size(); i++){
            Lit l = lits[i];
            if(S.value(l) == l_False && S.level(var(l)) == 0){
                //don't add this to the set
            }else if(S.value(l) == l_True && S.level(var(l)) == 0){
                if(constant_true_lit == lit_Undef){
                    constant_true_lit = l;
                }else{
                    //multiple literals are known to be true; this is a conflict
                    _addClause(~constant_true_lit, ~l);//this is a conflict
                    return;
                }
            }else{
                set.push(l);
            }
        }

        if(constant_true_lit == lit_Undef){
            if(set.size() < 2){
                //an empty set of literals, or a set with just one argument, always satisfies the AMO constraint
            }else{
                for(int i = 0; i < set.size(); i++){
                    for(int j = i + 1; j < set.size(); j++){
                        _addClause(~set[i], ~set[j]);
                    }
                }
            }
        }else{
            //all remaining elements of the set must be false, because constant_true_lit is true.
            for(int i = 0; i < set.size(); i++){
                if(set[i] != constant_true_lit){
                    _addClause(~constant_true_lit,
                               ~set[i]);//technically don't need to include ~constant_true_lit here, but it might make things cleaner to reason about elsewhere.
                    // (It will be eliminated by the solver anyhow, so this doesn't cost anything.)
                }
            }
        }

    }

    //uses n^2 binary clauses to create a simple exactly-one-constraint.
    //if you have more than 20 or so literals, strongly consider using a pseudo-Boolean constraint solver instead
    void AssertExactlyOne(const vec<Lit>& lits){
        AssertOr(lits);
        AssertAMO(lits);
    }

};
};

#endif /* CIRCUIT_H_ */
