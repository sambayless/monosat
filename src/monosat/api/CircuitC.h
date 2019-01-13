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

#ifndef MONOSAT_CIRCUITC_H
#define MONOSAT_CIRCUITC_H

//Monosat circuit construction interface, in C
//
//Most functions have a few alternate forms:
//int And(SolverPtr S,int lit_a, int lit_b);
//int And_(SolverPtr S,int lit_a, int lit_b, int lit_out);
//AssertAnd(SolverPtr S,int lit_a, int lit_b);
//int Ands(SolverPtr S,int * lits, int n_lits);
//...

//This first of these creates an AND gate over lits_a, introducing a new literal to represent
//the output of thate gate and adding clauses to enforce the corresponding relationship.
//The new output literal is returned.
//The seccond function adds claues to enforce that lit_out is true iff lit_a and lit_b are true,
//but does not allocate a new output literal. For convenience, it also returns the output literal (lit_out).
//The third introduces no new output literal, and simply asserts that both lit_a and lit_b are true.
//The fourth creates an AND gate over all the literals in lits* (there must be exactly n_lits such literals),
//introducing a new literal to represent the output of that gate, and returns that literal.
//There are several combinations of the above, such as Ands_ or AssertAnds.


#include <stdint.h>
#include "Monosat.h"

#ifdef __cplusplus
extern "C"
{
#endif
int newLit(SolverPtr S);
int newLit_(SolverPtr S, bool decisionLit);
int getTrue(SolverPtr S);
int getFalse(SolverPtr S);

int And_(SolverPtr S, int lit_a, int lit_b, int lit_out);
int Ands_(SolverPtr S, int* lits, int n_lits, int lit_out);
void AssertImpliesAnd_(SolverPtr S, int implies, int* lits, int n_lits, int lit_out);

int Ands(SolverPtr S, int* lits, int n_lits);
int And(SolverPtr S, int lit_a, int lit_b);

int Or_(SolverPtr S, int lit_a, int lit_b, int lit_out);
int Ors_(SolverPtr S, int* lits, int n_lits, int lit_out);

//If this gate is true, then all of vals must be true.
//But if this gate is false, vals may be true or false.
int ImpliesAnd_(SolverPtr S, int* lits, int n_lits, int lit_out);
//If this gate is true, then at least one of vals must be true.
//But if this gate is false, vals may be true or false.
int ImpliesOr(SolverPtr S, int* lits, int n_lits);
int ImpliesOr_(SolverPtr S, int* lits, int n_lits, int lit_out);
//This is an OR condition that holds only if implies is true
void AssertImpliesOr_(SolverPtr S, int implies, int* lits, int n_lits, int lit_out);

void AssertImpliesOr(SolverPtr S, int implies, int* lits, int n_lits);
void AssertImpliesAnd(SolverPtr S, int implies, int* lits, int n_lits);
void AssertImpliesAnd_(SolverPtr S, int implies, int* lits, int n_lits, int lit_out);

int Ors(SolverPtr S, int* lits, int n_lits);
int Or(SolverPtr S, int lit_a, int lit_b);

int Nors(SolverPtr S, int* lits, int n_lits);
int Nor(SolverPtr S, int lit_a, int lit_b);

int Nands(SolverPtr S, int* lits, int n_lits);

int Nand(SolverPtr S, int lit_a, int lit_b);

int Xors(SolverPtr S, int* lits, int n_lits);
int Xor(SolverPtr S, int lit_a, int lit_b);


int Xnors(SolverPtr S, int* lits, int n_lits);
int Xnor(SolverPtr S, int lit_a, int lit_b);

int Implies(SolverPtr S, int lit_a, int lit_b);
int Implies_(SolverPtr S, int lit_a, int lit_b, int lit_out);

int Ite(SolverPtr S, int lit_cond, int lit_thn, int lit_els);
int Ite_(SolverPtr S, int lit_cond, int lit_thn, int lit_els, int lit_out);

int Add(SolverPtr S, int* lits_a, int* lits_b, int n_lits, int* lits_out);
int Add_(SolverPtr S, int* lits_a, int* lits_b, int n_lits, int* lits_out, int carry_lit);
int Subtract(SolverPtr S, int* lits_a, int* lits_b, int n_lits, int* lits_out);
int Subtract_(SolverPtr S, int* lits_a, int* lits_b, int n_lits, int* lits_out, int borrow_lit);

//perform two's complement negation (SolverPtr S,invert bits and add 1)
//int * lits_out must have enough room to hold n_lits
void Negate(SolverPtr S, int* lits, int n_lits, int* lits_out);
//perform two's complement negation (SolverPtr S,invert bits and add 1)
void Negate_(SolverPtr S, int* lits, int n_lits, int* lits_out);

void Assert(SolverPtr S, int lit);
void AssertOrTertiary(SolverPtr S, int lit_a, int lit_b, int lit_c);
void AssertOrs(SolverPtr S, int* lits, int n_lits);
void AssertOr(SolverPtr S, int lit_a, int lit_b);

void AssertNands(SolverPtr S, int* lits, int n_lits);
void AssertNand(SolverPtr S, int lit_a, int lit_b);

void AssertAnds(SolverPtr S, int* lits, int n_lits);
void AssertAnd(SolverPtr S, int lit_a, int lit_b);

void AssertNors(SolverPtr S, int* lits, int n_lits);
void AssertNor(SolverPtr S, int lit_a, int lit_b);

void AssertXor(SolverPtr S, int lit_a, int lit_b);
void AssertXors(SolverPtr S, int* lits, int n_lits);

void AssertXnors(SolverPtr S, int* lits, int n_lits);
void AssertXnor(SolverPtr S, int lit_a, int lit_b);

void AssertImplies(SolverPtr S, int lit_a, int lit_b);
void AssertEqual(SolverPtr S, int lit_a, int lit_b);

void AssertAllSame(SolverPtr S, int* lits, int n_lits);

int Equal(SolverPtr S, int a_lit, int b_lit);
int Equals(SolverPtr S, int* A_lits, int* B_lits, int n_lits);

int LEQ(SolverPtr S, int* A_lits, int* B_lits, int n_lits);
int LT(SolverPtr S, int* A_lits, int* B_lits, int n_lits);

void AssertEquals(SolverPtr S, int* A_lits, int* B_lits, int n_lits);
void AssertLEQ(SolverPtr S, int* A_lits, int* B_lits, int n_lits);
void AssertLT(SolverPtr S, int* A_lits, int* B_lits, int n_lits);

//uses n^2 binary clauses to create a simple at-most-one constraint.
//if you have more than 20 or so literals, strongly consider using a pseudo-Boolean constraint solver instead
void AssertAMO(SolverPtr S, int* lits, int n_lits);
//uses n^2 binary clauses to create a simple exactly-one-constraint.
//if you have more than 20 or so literals, strongly consider using a pseudo-Boolean constraint solver instead
void AssertExactlyOne(SolverPtr S, int* lits, int n_lits);


#ifdef __cplusplus
}
#endif

#endif //MONOSAT_CIRCUITC_H
