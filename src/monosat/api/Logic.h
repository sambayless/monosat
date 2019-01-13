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

#ifndef LOGIC_H_
#define LOGIC_H_

#include "monosat/api/Circuit.h"
#include "monosat/core/Solver.h"
#include "monosat/simp/SimpSolver.h"
#include "monosat/core/SolverTypes.h"
#include "monosat/bv/BVTheorySolver.h"
#include "monosat/graph/GraphTheory.h"
#include "monosat/mtl/Vec.h"
#include <list>

//High-level C++ interface to Monosat
//Note: this is incomplete, use with caution!
namespace Monosat {

namespace Internal {
bool api_owns_solver = false;
//true if the solver was created (and should be deleted by) the api
bool api_owns_circuit = false;
//true if the circuit was created (and should be deleted by) the api
thread_local Circuit<SimpSolver>* circ_ptr = nullptr;
}

void clearSolver(){
    if(Internal::circ_ptr){
        if(Internal::api_owns_solver){
            SimpSolver* S = &Internal::circ_ptr->getSolver();
            delete (S);//but what if the solver was created elsewhere?
        }
        if(Internal::api_owns_circuit){
            delete (Internal::circ_ptr);
        }
        Internal::circ_ptr = nullptr;
        Internal::api_owns_solver = false;
        Internal::api_owns_circuit = false;
    }
}

SimpSolver* newSolver(){
    clearSolver();
    Internal::api_owns_solver = true;
    Internal::api_owns_circuit = true;
    SimpSolver* S = new SimpSolver();
    BVTheorySolver<int64_t>* bv = new BVTheorySolver<int64_t>(S);
    Internal::circ_ptr = new Circuit<SimpSolver>(*S);
}


void setSolver(SimpSolver* S, Circuit<SimpSolver>* circuit = nullptr){
    clearSolver();
    if(!S)
        return;
    Internal::api_owns_solver = false;
    if(!circuit){
        Internal::api_owns_circuit = true;
        circuit = new Circuit<SimpSolver>(*S);
    }
    Internal::circ_ptr = circuit;
}

namespace Internal {
Circuit<SimpSolver>& getCircuit(){
    if(!Internal::circ_ptr){
        newSolver();
    }
    return *Internal::circ_ptr;
}
}


SimpSolver* getSolver(){
    return &Internal::getCircuit().getSolver();
}

bool Solve(vec<Lit>& assumps){
    return getSolver()->solve(assumps, true, true);
}

bool Solve(Lit p){
    vec<Lit> s;
    s.push(p);
    return Solve(s);
}

bool Solve(Lit p, Lit q){
    vec<Lit> s;
    s.push(p);
    s.push(q);
    return Solve(s);
}

bool Solve(Lit p, Lit q, Lit r){
    vec<Lit> s;
    s.push(p);
    s.push(q);
    s.push(r);
    return Solve(s);
}

bool Solve(){
    static vec<Lit> ignore;
    return Solve(ignore);
}


Lit True(){
    return Internal::getCircuit().getTrue();
}

Lit False(){
    return Internal::getCircuit().getFalse();
}

Lit And(Lit a){
    return Internal::getCircuit().And(a);
}

Lit And(Lit a, Lit b){
    return Internal::getCircuit().And(a, b);
}

Lit And(const std::list<Lit>& vals){
    return Internal::getCircuit().And(vals);
}

Lit And(const vec<Lit>& vals){
    return Internal::getCircuit().And(vals);
}

template<typename... Args>
Lit And(Lit a, Lit b, Args... args){
    return Internal::getCircuit().And(a, b, args...);
}


Lit Or(Lit a){
    return Internal::getCircuit().Or(a);
}

Lit Or(Lit a, Lit b){
    return Internal::getCircuit().Or(a, b);
}

Lit Or(const std::list<Lit>& vals){
    return Internal::getCircuit().Or(vals);
}

Lit Or(const vec<Lit>& vals){
    return Internal::getCircuit().Or(vals);
}

template<typename... Args>
Lit Or(Lit a, Lit b, Args... args){
    return Internal::getCircuit().Or(a, b, args...);
}


Lit Nor(Lit a){
    return Internal::getCircuit().Nor(a);
}

Lit Nor(Lit a, Lit b){
    return Internal::getCircuit().Nor(a, b);
}

Lit Nor(const std::list<Lit>& vals){
    return Internal::getCircuit().Nor(vals);
}

Lit Nor(const vec<Lit>& vals){
    return Internal::getCircuit().Nor(vals);
}

template<typename... Args>
Lit Nor(Lit a, Lit b, Args... args){
    return Internal::getCircuit().Nor(a, b, args...);
}


Lit Nand(Lit a){
    return Internal::getCircuit().Nand(a);
}

Lit Nand(Lit a, Lit b){
    return Internal::getCircuit().Nand(a, b);
}

Lit Nand(const std::list<Lit>& vals){
    return Internal::getCircuit().Nand(vals);

}

Lit Nand(const vec<Lit>& vals){
    return Internal::getCircuit().Nand(vals);
}

template<typename... Args>
Lit Nand(Lit a, Lit b, Args... args){
    return Internal::getCircuit().Nand(a, b, args...);
}

Lit Xor(Lit a){
    return Internal::getCircuit().Xor(a);
}

Lit Xor(Lit a, Lit b){
    return Internal::getCircuit().Xor(a, b);
}

Lit Xor(const std::list<Lit>& vals){
    return Internal::getCircuit().Xor(vals);

}

Lit Xor(const vec<Lit>& vals){
    return Internal::getCircuit().Xor(vals);
}

template<typename... Args>
Lit Xor(Lit a, Lit b, Args... args){
    return Internal::getCircuit().Xor(a, b, args...);
}


Lit Xnor(Lit a){
    return Internal::getCircuit().Xnor(a);
}

Lit Xnor(Lit a, Lit b){
    return Internal::getCircuit().Xnor(a, b);
}

Lit Xnor(const std::list<Lit>& vals){
    return Internal::getCircuit().Xnor(vals);
}

Lit Xnor(const vec<Lit>& vals){
    return Internal::getCircuit().Xnor(vals);
}

template<typename... Args>
Lit Xnor(Lit a, Lit b, Args... args){
    return Internal::getCircuit().Xnor(a, b, args...);
}

Lit Implies(Lit a, Lit b){
    return Internal::getCircuit().Implies(a, b);
}

Lit Ite(Lit cond, Lit thn, Lit els){
    return Internal::getCircuit().Ite(cond, thn, els);
}


Lit halfAdder(Lit a, Lit b, Lit& carry_out){
    return Internal::getCircuit().halfAdder(a, b, carry_out);
}

void Add(vec<Lit>& a, vec<Lit>& b, vec<Lit>& store_out, Lit& carry_out){
    return Internal::getCircuit().add(a, b, store_out, carry_out);
}


void Assert(Lit l){
    return Internal::getCircuit().Assert(l);
}

void AssertOr(Lit a){
    return Internal::getCircuit().AssertOr(a);
}

void AssertOr(Lit a, Lit b){
    return Internal::getCircuit().AssertOr(a, b);
}

void AssertOr(const std::list<Lit>& vals){
    return Internal::getCircuit().AssertOr(vals);
}

void AssertOr(const vec<Lit>& vals){
    return Internal::getCircuit().AssertOr(vals);
}

template<typename... Args>
void AssertOr(Lit a, Lit b, Args... args){
    return Internal::getCircuit().AssertOr(a, b, args...);
}

void AssertNand(Lit a){
    return Internal::getCircuit().AssertNand(a);
}

void AssertNand(Lit a, Lit b){
    return Internal::getCircuit().AssertNand(a, b);
}

void AssertNand(const std::list<Lit>& vals){
    return Internal::getCircuit().AssertNand(vals);
}

void AssertNand(const vec<Lit>& vals){
    return Internal::getCircuit().AssertNand(vals);
}

template<typename... Args>
void AssertNand(Lit a, Lit b, Args... args){
    return Internal::getCircuit().AssertNand(a, b, args...);
}

void AssertAnd(Lit a){
    return Internal::getCircuit().AssertAnd(a);
}

void AssertAnd(Lit a, Lit b){
    return Internal::getCircuit().AssertAnd(a, b);
}

void AssertAnd(const std::list<Lit>& vals){
    return Internal::getCircuit().AssertAnd(vals);
}

void AssertAnd(const vec<Lit>& vals){
    return Internal::getCircuit().AssertAnd(vals);
}

template<typename... Args>
void AssertAnd(Lit a, Lit b, Args... args){
    return Internal::getCircuit().AssertAnd(a, b, args...);
}

void AssertNor(Lit a){
    return Internal::getCircuit().AssertNor(a);
}

void AssertNor(Lit a, Lit b){
    return Internal::getCircuit().AssertNor(a, b);
}

void AssertNor(const std::list<Lit>& vals){
    return Internal::getCircuit().AssertNor(vals);
}

void AssertNor(const vec<Lit>& vals){
    return Internal::getCircuit().AssertNor(vals);
}

template<typename... Args>
void AssertNor(Lit a, Lit b, Args... args){
    return Internal::getCircuit().AssertNor(a, b, args...);
}

void AssertXor(Lit a){
    return Internal::getCircuit().AssertXor(a);
}

void AssertXor(Lit a, Lit b){
    return Internal::getCircuit().AssertXor(a, b);
}

void AssertXor(const std::list<Lit>& vals){
    return Internal::getCircuit().AssertXor(vals);
}

void AssertXor(const vec<Lit>& vals){
    return Internal::getCircuit().AssertXor(vals);
}

template<typename... Args>
void AssertXor(Lit a, Lit b, Args... args){
    return Internal::getCircuit().AssertXor(a, b, args...);
}

void AssertXnor(Lit a){
    return Internal::getCircuit().AssertXnor(a);
}

void AssertXnor(Lit a, Lit b){
    return Internal::getCircuit().AssertXnor(a, b);
}

void AssertXnor(const std::list<Lit>& vals){
    return Internal::getCircuit().AssertXnor(vals);
}

void AssertXnor(const vec<Lit>& vals){
    return Internal::getCircuit().AssertXnor(vals);
}

template<typename... Args>
void AssertXnor(Lit a, Lit b, Args... args){
    return Internal::getCircuit().AssertXnor(a, b, args...);
}

void AssertImplies(Lit a, Lit b){
    return Internal::getCircuit().AssertImplies(a, b);
}

void AssertEqual(Lit a, Lit b){
    return Internal::getCircuit().AssertEqual(a, b);
}

void AssertEqual(const std::list<Lit>& vals){
    return Internal::getCircuit().AssertEqual(vals);
}

void AssertEqual(const vec<Lit>& vals){
    return Internal::getCircuit().AssertEqual(vals);
}

template<typename... Args>
void AssertEqual(Lit a, Lit b, Args... args){
    return Internal::getCircuit().AssertEqual(a, b, args...);
}

namespace Internal {
BVTheorySolver<int64_t>* getBVTheory(){
    SimpSolver& S = Internal::getCircuit().getSolver();
    BVTheorySolver<int64_t>* bv = (BVTheorySolver<int64_t>*) S.getBVTheory();
    if(!bv){
        bv = new BVTheorySolver<int64_t>(&S);
        S.setBVTheory(bv);//ensure that the solver has a bitvector theory
    }
    return bv;
}
}

GraphTheorySolver<int64_t>* newGraph(){
    SimpSolver& S = Internal::getCircuit().getSolver();
    GraphTheorySolver<int64_t>* G = new GraphTheorySolver<int64_t>(&S);

    G->setBVTheory(Internal::getBVTheory());//you only need this if you are using weighted edges
    return G;
}

//TODO: High level BitVector and Graph theory api...

};

#endif /* LOGIC_H_ */
