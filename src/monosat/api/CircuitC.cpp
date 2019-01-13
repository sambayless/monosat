/**************************************************************************************************
 The MIT License (MIT)

 Copyright (c) 2018, Sam Bayless

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

#include "Circuit.h"
#include "Monosat.h"
#include "CircuitC.h"
#include "MonosatInternal.h"

//bridge between CircuitC.h's c interface, and Circuit.h's C++ interface.
using namespace Monosat;
using namespace std;

//Convert an external integer representation of a lit,
//into an internal solver literal
Lit internalLit_(Monosat::SimpSolver* S, int l){
    return S->mapLit(toLit(l));
}

//Convert an internal solver lit into an external integer representation
int externalLit_(Monosat::SimpSolver* S, Lit l){
    return toInt(S->unmap(l));
}

Var internalVar_(Monosat::SimpSolver* S, int v){
    return S->mapVar(v);
}

//Convert an internal solver lit into an external integer representation
Var externalVar_(Monosat::SimpSolver* S, Var v){
    return S->unmap(v);
}


//Convert an external integer representation of a lit into dimacs_ format
inline int dimacs_(int externalLit_){
    Lit l = toLit(externalLit_);
    return sign(l) ? -(var(l) + 1) : (var(l) + 1);
}

inline int dimacs_(Monosat::SimpSolver* S, Lit internalLit_){
    Lit l = S->unmap(internalLit_);
    return sign(l) ? -(var(l) + 1) : (var(l) + 1);
}


vec<Lit> tmp_lits_a;
vec<Lit> tmp_lits_b;
vec<Lit> tmp_lits_c;

void toVec(Monosat::SimpSolver* S, int* lits, int n_lits, vec<Lit>& store){
    assert(lits);
    assert(n_lits >= 0);
    store.clear();
    for(int i = 0; i < n_lits; i++){
        int l = lits[i];
        assert(l >= 0);
        Lit lit = internalLit_(S, l);
        assert(lit != lit_Undef);
        assert(lit != lit_Error);
        store.push(lit);
    }
}

int newLit(Monosat::SimpSolver* S){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    return externalLit_(S, circuit.newLit());
}

int newLit_(Monosat::SimpSolver* S, bool decisionLit){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    return externalLit_(S, circuit.newLit(decisionLit));
}

int getTrue(Monosat::SimpSolver* S){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    return externalLit_(S, circuit.getTrue());
}

int getFalse(Monosat::SimpSolver* S){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    return externalLit_(S, circuit.getFalse());
}

int And_(Monosat::SimpSolver* S, int lit_a, int lit_b, int lit_out){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    return externalLit_(S, circuit.And_(internalLit_(S, lit_a), internalLit_(S, lit_b), internalLit_(S, lit_out)));
}

int Ands_(Monosat::SimpSolver* S, int* lits, int n_lits, int lit_out){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    toVec(S, lits, n_lits, tmp_lits_a);
    return externalLit_(S, circuit.And_(tmp_lits_a, internalLit_(S, lit_out)));
}

int Ands(Monosat::SimpSolver* S, int* lits, int n_lits){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    toVec(S, lits, n_lits, tmp_lits_a);
    return externalLit_(S, circuit.And(tmp_lits_a));
}

int And(Monosat::SimpSolver* S, int lit_a, int lit_b){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    return externalLit_(S, circuit.And(internalLit_(S, lit_a), internalLit_(S, lit_b)));
}

int Or_(Monosat::SimpSolver* S, int lit_a, int lit_b, int lit_out){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    return externalLit_(S, circuit.Or_(internalLit_(S, lit_a), internalLit_(S, lit_b), internalLit_(S, lit_out)));
}

int Ors_(Monosat::SimpSolver* S, int* lits, int n_lits, int lit_out){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    toVec(S, lits, n_lits, tmp_lits_a);
    return externalLit_(S, circuit.Or_(tmp_lits_a, internalLit_(S, lit_out)));
}

//If this gate is true, then all of vals must be true.
//But if this gate is false, vals may be true or false.
int ImpliesAnd_(Monosat::SimpSolver* S, int* lits, int n_lits, int lit_out){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    toVec(S, lits, n_lits, tmp_lits_a);
    return externalLit_(S, circuit.ImpliesAnd(tmp_lits_a, internalLit_(S, lit_out)));
}

//If this gate is true, then at least one of vals must be true.
//But if this gate is false, vals may be true or false.
int ImpliesOr(Monosat::SimpSolver* S, int* lits, int n_lits){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    toVec(S, lits, n_lits, tmp_lits_a);
    return externalLit_(S, circuit.ImpliesOr(tmp_lits_a));
}

int ImpliesOr_(Monosat::SimpSolver* S, int* lits, int n_lits, int lit_out){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    toVec(S, lits, n_lits, tmp_lits_a);
    return externalLit_(S, circuit.ImpliesOr(tmp_lits_a, internalLit_(S, lit_out)));
}

//This is an OR condition that holds only if implies is true
void AssertImpliesOr_(Monosat::SimpSolver* S, int implies, int* lits, int n_lits, int lit_out){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    toVec(S, lits, n_lits, tmp_lits_a);
    circuit.AssertImpliesOr_(internalLit_(S, implies), tmp_lits_a, internalLit_(S, lit_out));
}

void AssertImpliesAnd_(Monosat::SimpSolver* S, int implies, int* lits, int n_lits, int lit_out){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;

    toVec(S, lits, n_lits, tmp_lits_a);
    circuit.AssertImpliesAnd_(internalLit_(S, implies), tmp_lits_a, internalLit_(S, lit_out));
}

void AssertImpliesOr(Monosat::SimpSolver* S, int implies, int* lits, int n_lits){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    toVec(S, lits, n_lits, tmp_lits_a);
    circuit.AssertImpliesOr_(internalLit_(S, implies), tmp_lits_a, lit_Undef);
}

void AssertImpliesAnd(Monosat::SimpSolver* S, int implies, int* lits, int n_lits){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;

    toVec(S, lits, n_lits, tmp_lits_a);
    circuit.AssertImpliesAnd_(internalLit_(S, implies), tmp_lits_a, lit_Undef);
}

int Ors(Monosat::SimpSolver* S, int* lits, int n_lits){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    toVec(S, lits, n_lits, tmp_lits_a);
    return externalLit_(S, circuit.Or(tmp_lits_a));
}

int Or(Monosat::SimpSolver* S, int lit_a, int lit_b){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    return externalLit_(S, circuit.Or(internalLit_(S, lit_a), internalLit_(S, lit_b)));
}

int Nors(Monosat::SimpSolver* S, int* lits, int n_lits){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    toVec(S, lits, n_lits, tmp_lits_a);
    return externalLit_(S, circuit.Nor(tmp_lits_a));
}

int Nor(Monosat::SimpSolver* S, int lit_a, int lit_b){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    return externalLit_(S, circuit.Nor(internalLit_(S, lit_a), internalLit_(S, lit_b)));
}

int Nands(Monosat::SimpSolver* S, int* lits, int n_lits){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    toVec(S, lits, n_lits, tmp_lits_a);
    return externalLit_(S, circuit.Nand(tmp_lits_a));
}

int Nand(Monosat::SimpSolver* S, int lit_a, int lit_b){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    return externalLit_(S, circuit.Nand(internalLit_(S, lit_a), internalLit_(S, lit_b)));
}

int Xors(Monosat::SimpSolver* S, int* lits, int n_lits){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    toVec(S, lits, n_lits, tmp_lits_a);
    return externalLit_(S, circuit.Xor(tmp_lits_a));
}

int Xor(Monosat::SimpSolver* S, int lit_a, int lit_b){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    return externalLit_(S, circuit.Xor(internalLit_(S, lit_a), internalLit_(S, lit_b)));
}


int Xnors(Monosat::SimpSolver* S, int* lits, int n_lits){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    toVec(S, lits, n_lits, tmp_lits_a);
    return externalLit_(S, circuit.Xnor(tmp_lits_a));
}

int Xnor(Monosat::SimpSolver* S, int lit_a, int lit_b){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    return externalLit_(S, circuit.Xnor(internalLit_(S, lit_a), internalLit_(S, lit_b)));
}

int Implies(Monosat::SimpSolver* S, int lit_a, int lit_b){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    return externalLit_(S, circuit.Implies(internalLit_(S, lit_a), internalLit_(S, lit_b)));
}

int Implies_(Monosat::SimpSolver* S, int lit_a, int lit_b, int lit_out){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    return externalLit_(S, circuit.Implies_(internalLit_(S, lit_a), internalLit_(S, lit_b), internalLit_(S, lit_out)));
}

int Ite(Monosat::SimpSolver* S, int lit_cond, int lit_thn, int lit_els){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    return externalLit_(S, circuit.Ite(internalLit_(S, lit_cond), internalLit_(S, lit_thn), internalLit_(S, lit_els)));
}

int Ite_(Monosat::SimpSolver* S, int lit_cond, int lit_thn, int lit_els, int lit_out){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    return externalLit_(S, circuit.Ite_(internalLit_(S, lit_cond), internalLit_(S, lit_thn), internalLit_(S, lit_els),
                                        internalLit_(S, lit_out)));
}

int Add(Monosat::SimpSolver* S, int* lits_a, int* lits_b, int n_lits, int* lits_out){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    toVec(S, lits_a, n_lits, tmp_lits_a);
    toVec(S, lits_b, n_lits, tmp_lits_b);
    tmp_lits_c.clear();
    Lit carry = lit_Undef;
    circuit.Add(tmp_lits_a, tmp_lits_b, tmp_lits_c, carry);
    assert(lits_out);
    assert(tmp_lits_c.size() == n_lits);
    for(int i = 0; i < tmp_lits_c.size(); i++){
        lits_out[i] = externalLit_(S, tmp_lits_c[i]);
    }
    return externalLit_(S, carry);
}

int Add_(Monosat::SimpSolver* S, int* lits_a, int* lits_b, int n_lits, int* lits_out, int carry_lit){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    toVec(S, lits_a, n_lits, tmp_lits_a);
    toVec(S, lits_b, n_lits, tmp_lits_b);
    assert(lits_out);
    toVec(S, lits_out, n_lits, tmp_lits_c);
    Lit carry = lit_Undef;
    circuit.Add_(tmp_lits_a, tmp_lits_b, tmp_lits_c, carry);
    return externalLit_(S, carry);
}

int Subtract(Monosat::SimpSolver* S, int* lits_a, int* lits_b, int n_lits, int* lits_out){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    toVec(S, lits_a, n_lits, tmp_lits_a);
    toVec(S, lits_b, n_lits, tmp_lits_b);
    tmp_lits_c.clear();
    Lit carry = lit_Undef;
    circuit.Subtract(tmp_lits_a, tmp_lits_b, tmp_lits_c, carry);
    assert(lits_out);
    assert(tmp_lits_c.size() == n_lits);
    for(int i = 0; i < tmp_lits_c.size(); i++){
        lits_out[i] = externalLit_(S, tmp_lits_c[i]);
    }
    return externalLit_(S, carry);
}

int Subtract_(Monosat::SimpSolver* S, int* lits_a, int* lits_b, int n_lits, int* lits_out, int borrow_lit){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    toVec(S, lits_a, n_lits, tmp_lits_a);
    toVec(S, lits_b, n_lits, tmp_lits_b);
    assert(lits_out);

    toVec(S, lits_out, n_lits, tmp_lits_c);
    Lit carry = lit_Undef;
    circuit.Subtract_(tmp_lits_a, tmp_lits_b, tmp_lits_c, carry);
    return externalLit_(S, carry);
}

//perform two's complement negation (invert bits and add 1)
//int * lits_out must have enough room to hold n_lits
void Negate(Monosat::SimpSolver* S, int* lits, int n_lits, int* lits_out){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    toVec(S, lits, n_lits, tmp_lits_a);
    tmp_lits_b.clear();
    circuit.Negate(tmp_lits_a, tmp_lits_b);
    assert(lits_out);
    assert(tmp_lits_b.size() == n_lits);
    for(int i = 0; i < tmp_lits_b.size(); i++){
        lits_out[i] = externalLit_(S, tmp_lits_b[i]);
    }
}

//perform two's complement negation (invert bits and add 1)
void Negate_(Monosat::SimpSolver* S, int* lits, int n_lits, int* lits_out){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    toVec(S, lits, n_lits, tmp_lits_a);
    assert(lits_out);
    toVec(S, lits_out, n_lits, tmp_lits_b);
    circuit.Negate_(tmp_lits_a, tmp_lits_b);
}

void Assert(Monosat::SimpSolver* S, int lit){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    circuit.Assert(internalLit_(S, lit));
}

void AssertOrTertiary(Monosat::SimpSolver* S, int lit_a, int lit_b, int lit_c){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    circuit.AssertOr(internalLit_(S, lit_a), internalLit_(S, lit_b), internalLit_(S, lit_c));
}

void AssertOrs(Monosat::SimpSolver* S, int* lits, int n_lits){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    toVec(S, lits, n_lits, tmp_lits_a);
    circuit.AssertOr(tmp_lits_a);
}

void AssertOr(Monosat::SimpSolver* S, int lit_a, int lit_b){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    circuit.AssertOr(internalLit_(S, lit_a), internalLit_(S, lit_b));
}

void AssertNands(Monosat::SimpSolver* S, int* lits, int n_lits){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    toVec(S, lits, n_lits, tmp_lits_a);
    circuit.AssertNand(tmp_lits_a);
}

void AssertNand(Monosat::SimpSolver* S, int lit_a, int lit_b){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    circuit.AssertNand(internalLit_(S, lit_a), internalLit_(S, lit_b));
}

void AssertAnds(Monosat::SimpSolver* S, int* lits, int n_lits){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    toVec(S, lits, n_lits, tmp_lits_a);
    circuit.AssertAnd(tmp_lits_a);
}

void AssertAnd(Monosat::SimpSolver* S, int lit_a, int lit_b){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    circuit.AssertAnd(internalLit_(S, lit_a), internalLit_(S, lit_b));
}

void AssertNors(Monosat::SimpSolver* S, int* lits, int n_lits){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    toVec(S, lits, n_lits, tmp_lits_a);
    circuit.AssertNor(tmp_lits_a);
}

void AssertNor(Monosat::SimpSolver* S, int lit_a, int lit_b){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    circuit.AssertNor(internalLit_(S, lit_a), internalLit_(S, lit_b));
}

void AssertXor(Monosat::SimpSolver* S, int lit_a, int lit_b){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    circuit.AssertXor(internalLit_(S, lit_a), internalLit_(S, lit_b));
}

void AssertXors(Monosat::SimpSolver* S, int* lits, int n_lits){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    toVec(S, lits, n_lits, tmp_lits_a);
    circuit.AssertXor(tmp_lits_a);
}

void AssertXnors(Monosat::SimpSolver* S, int* lits, int n_lits){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    toVec(S, lits, n_lits, tmp_lits_a);
    circuit.AssertXnor(tmp_lits_a);
}

void AssertXnor(Monosat::SimpSolver* S, int lit_a, int lit_b){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    circuit.AssertXnor(internalLit_(S, lit_a), internalLit_(S, lit_b));
}

void AssertImplies(Monosat::SimpSolver* S, int lit_a, int lit_b){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    circuit.AssertImplies(internalLit_(S, lit_a), internalLit_(S, lit_b));
}

void AssertEqual(Monosat::SimpSolver* S, int lit_a, int lit_b){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    circuit.AssertEqual(internalLit_(S, lit_a), internalLit_(S, lit_b));
}

void AssertAllSame(Monosat::SimpSolver* S, int* lits, int n_lits){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    toVec(S, lits, n_lits, tmp_lits_a);
    //circuit.AssertEqual(tmp_lits_a);
    throw std::runtime_error("Not yet implemetned");
}


int Equals(Monosat::SimpSolver* S, int* lits_a, int* lits_b, int n_lits){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    toVec(S, lits_a, n_lits, tmp_lits_a);
    toVec(S, lits_b, n_lits, tmp_lits_b);
    throw std::runtime_error("Not yet implemetned");
    // return externalLit_(S,circuit.Equal(tmp_lits_a,tmp_lits_b));
}

int LEQ(Monosat::SimpSolver* S, int* lits_a, int* lits_b, int n_lits){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    toVec(S, lits_a, n_lits, tmp_lits_a);
    toVec(S, lits_b, n_lits, tmp_lits_b);
    return externalLit_(S, circuit.LEQ(tmp_lits_a, tmp_lits_b));
}

int LT(Monosat::SimpSolver* S, int* lits_a, int* lits_b, int n_lits){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    toVec(S, lits_a, n_lits, tmp_lits_a);
    toVec(S, lits_b, n_lits, tmp_lits_b);
    return externalLit_(S, circuit.LT(tmp_lits_a, tmp_lits_b));
}

void AssertEquals(Monosat::SimpSolver* S, int* lits_a, int* lits_b, int n_lits){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    toVec(S, lits_a, n_lits, tmp_lits_a);
    toVec(S, lits_b, n_lits, tmp_lits_b);
    circuit.AssertEqual(tmp_lits_a, tmp_lits_b);
}

void AssertLEQ(Monosat::SimpSolver* S, int* lits_a, int* lits_b, int n_lits){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    toVec(S, lits_a, n_lits, tmp_lits_a);
    toVec(S, lits_b, n_lits, tmp_lits_b);
    circuit.AssertLEQ(tmp_lits_a, tmp_lits_b);
}

void AssertLT(Monosat::SimpSolver* S, int* lits_a, int* lits_b, int n_lits){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    toVec(S, lits_a, n_lits, tmp_lits_a);
    toVec(S, lits_b, n_lits, tmp_lits_b);
    circuit.AssertLT(tmp_lits_a, tmp_lits_b);
}

//uses n^2 binary clauses to create a simple at-most-one constraint.
//if you have more than 20 or so literals, strongly consider using a pseudo-Boolean constraint solver instead
void AssertAMO(Monosat::SimpSolver* S, int* lits, int n_lits){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    toVec(S, lits, n_lits, tmp_lits_a);
    circuit.AssertAMO(tmp_lits_a);
}

//uses n^2 binary clauses to create a simple exactly-one-constraint.
//if you have more than 20 or so literals, strongly consider using a pseudo-Boolean constraint solver instead
void AssertExactlyOne(Monosat::SimpSolver* S, int* lits, int n_lits){
    MonosatData* d = (MonosatData*) S->_external_data;
    assert(d);
    Monosat::Circuit<Monosat::SimpSolver>& circuit = d->circuit;
    toVec(S, lits, n_lits, tmp_lits_a);
    circuit.AssertExactlyOne(tmp_lits_a);
}
