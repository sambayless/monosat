/*
 * Aiger.cpp
 *
 *  Created on: 2012-10-18
 *      Author: sam
 */


#include "Aiger.h"
#include <iostream>
extern "C" {
#include "aiger/aiger.h"
}
#include <map>
using namespace Minisat;

static int
lit2int (aiger * mgr, unsigned a, int startVar)
{
  int sign = aiger_sign (a) ? -1 : 1;
  int res = aiger_lit2var (a);

  if (res){
    res = (res+startVar)* sign;
  }else
    res = sign;	//assume true and false are 0 and 1
  assert(res!=0);
  return res;
}


Lit getLit(Solver & S, int dimacs){
	if (dimacs == 0) return lit_Undef;
	int v = abs(dimacs)-1;

	while (v >= S.nVars()) S.newVar();
    return( (dimacs > 0) ? mkLit(v) : ~mkLit(v) );
}

Lit getLit(Solver & S,const int startVar, aiger * mgr, unsigned lit, vec<Var> & in_latches){

	if(aiger_symbol * s = aiger_is_latch(mgr,aiger_strip(lit))){
		 int regnum = s - mgr->latches;
		  int res = aiger_lit2var (lit);

		  assert(aiger_lit2var( s->lit)==res);

			 //this is a previous state reg - use the variable already in the solver.
		  assert(regnum<in_latches.size());
		  Var v = in_latches[regnum];
		  return mkLit(v, aiger_sign (lit) );

	}else{
		bool sign = aiger_sign(lit);
		int dimacs = lit2int(mgr,lit,startVar);

		return getLit(S,dimacs);
	}


}

Lit  unroll(Solver& S,aiger * mgr,int k, vec<Var> & in_latches, vec<Var> & out_latches){
	if(k==0){
		return lit_Undef;
	}
	Lit safety = unroll(S,mgr,in_latches,out_latches);
	for(int i = 1;i<k;i++){
		static vec<Var> t;
		t.clear();
		out_latches.copyTo(t);
		safety = unroll(S,mgr,t,out_latches);
	}
	return safety;
}

Lit  unroll(Solver& S,aiger * mgr, vec<Var> & in_latches, vec<Var> & out_latches){
	vec<Lit> c;
	int startVar = S.nVars();

	assert(S.nVars()>=1);//assume first var is const true.
	assert(mgr->num_outputs==1);
	out_latches.clear();

  for (int i = 0; i < mgr->num_ands; i++)
    {
      aiger_and *a = mgr->ands + i;
      c.push( getLit(S,startVar, mgr, aiger_not (a->lhs),in_latches));
      c.push(getLit(S,startVar, mgr, a->rhs0,in_latches));
      S.addClause(c);
      c.clear();

      c.push( getLit(S,startVar, mgr, aiger_not (a->lhs),in_latches));
      c.push(getLit(S,startVar, mgr, a->rhs1,in_latches));
      S.addClause(c);
      c.clear();


      c.push(getLit(S,startVar, mgr, a->lhs,in_latches));
      c.push( getLit(S,startVar, mgr, aiger_not (a->rhs0),in_latches));
      c.push( getLit(S,startVar, mgr, aiger_not (a->rhs1),in_latches));
        S.addClause(c);
        c.clear();

    }

  //check if any latches have input latches as their 'next' function; if so, construct a variable for their function
  for(int i = 0;i<mgr->num_latches;i++){
	  getLit(S,startVar,mgr,mgr->latches[i].next, in_latches);
	  //****This loop IS critical, even though getLit's return value is ignored, because it forces the construction of a variable for the next function in the (rare) case that
	  //that variable hasn't yet been defined. This can happen if the next value of a latch is an input latch, and that input latch is not referenced anywhere else in the circuit.
	  //Constructing this variable is necessary before the final output latch variables are added below
  }

  //return the property literal
  //need to handle the case where the ouput _is_ a latch.
  Lit output =  getLit(S,startVar, mgr, mgr->outputs[0].lit,in_latches);

  for(int i = 0;i<mgr->num_latches;i++){


			Lit l = getLit(S,startVar,mgr,mgr->latches[i].next, in_latches);
			//we need the latches to be consecutive literals in the solver, and in the correct order, so create new lits for them at the end of the solver:

				Var outv = S.newVar();
				Lit outLatch = mkLit(outv,false);

				c.push(~l);
				c.push(outLatch);

				S.addClause(c);
				c.clear();

				c.push(l);
				c.push(~outLatch);
				S.addClause(c);
				c.clear();

				out_latches.push(outv);
  }


  return output;
}

bool isReset(vec<Lit> & assignment){
	for(int i = 0;i<assignment.size();i++){
		if(!sign(assignment[i])){
			return false;
		}
	}
	return true;
}

void prepare(Minisat::Solver& target,aiger * mgr, Minisat::vec<Minisat::Var> & latches){

	assert(target.nVars()==0);
	latches.clear();
	//add constant literal ~1 (inverted to match AIGER format)  (this will automatically get dropped from clauses by typical sat solvers anyways)
	target.addClause(mkLit( target.newVar(),true));
	//add enough latches to represent the inital registers
	for(int i = 0;i<mgr->num_latches;i++){
		Var v = target.newVar();
		latches.push(v);
	}

}

int zero(Solver& target,aiger * mgr, vec<Var> & latches){
	int constraint_size=0;

	//force these latches to zero - aiger always initializes latches to zero
	for(int i = 0;i<latches.size();i++){

		Lit l = mkLit(latches[i],true);
		target.addClause(l);
		constraint_size++;
	}
	return constraint_size;
}





