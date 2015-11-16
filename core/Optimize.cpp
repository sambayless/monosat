/*
 * Optimize.cpp
 *
 *  Created on: Nov 12, 2015
 *      Author: sam
 */

#ifndef OPTIMIZE_CPP_
#define OPTIMIZE_CPP_
#include "Optimize.h"
#include "core/Config.h"

#include <stdexcept>
#include <cstdarg>
namespace Monosat{

long optimize_linear(Monosat::SimpSolver * S, Monosat::BVTheorySolver<long> * bvTheory,const vec<Lit> & assumes,int bvID, int conflict_limit, bool & hit_cutoff, long & n_solves){
	vec<Lit> assume;
	for(Lit l:assumes)
		assume.push(l);
	vec<Lit> last_satisfying_assign;
	for(Var v = 0;v<S->nVars();v++){
		if(S->value(v)==l_True){
			last_satisfying_assign.push(mkLit(v));
		}else if(S->value(v)==l_False){
			last_satisfying_assign.push(mkLit(v,true));
		}else{
			//this variable was unassigned.
		}
	}

	long value = bvTheory->getUnderApprox(bvID);
	long last_decision_value=value;
	  if(opt_verb>=1){
		  printf("Min bv%d = %ld",bvID,value);
	  }
	  // int bvID,const Weight & to, Var outerVar = var_Undef, bool decidable=true
	  Lit last_decision_lit =  bvTheory->toSolver(bvTheory->newComparison(Comparison::leq,bvID,value,var_Undef,opt_decide_optimization_lits));
	  while(value>bvTheory->getUnderApprox(bvID,true)){
		  Lit decision_lit = bvTheory->toSolver(bvTheory->newComparison(Comparison::leq,bvID,value-1,var_Undef,opt_decide_optimization_lits));
		  assume.push(decision_lit);
		  n_solves++;

		  bool r;

		  if(opt_limit_conflicts>0 || conflict_limit>=0){
			  S->setConfBudget(std::min((int)opt_limit_conflicts,conflict_limit));
			  lbool res = S->solveLimited(assume,false,false);
			  if (res==l_Undef){
				  if(opt_verb>0){
					  printf("\nBudget exceeded during optimization, quiting early (model might not be optimal!)\n");
				  }
				  r=false;
			  }else{
				  r = res==l_True;
			  }
		  }else
			  r = S->solve(assume,false,false);


		  if (r){
			  last_satisfying_assign.clear();
				for(Var v = 0;v<S->nVars();v++){
					if(S->value(v)==l_True){
						last_satisfying_assign.push(mkLit(v));
					}else if(S->value(v)==l_False){
						last_satisfying_assign.push(mkLit(v,true));
					}else{
						//this variable was unassigned.
					}
				}
			  last_decision_lit=decision_lit;
			  last_decision_value=value-1;
			  if(S->value(decision_lit)!=l_True){
				  throw std::runtime_error("Error in optimization (comparison not enforced)");
			  }
			  for(Lit l:assume){
					if(S->value(l)!=l_True){
						throw std::runtime_error("Error in optimization (model is inconsistent with assumptions)");
					}
				}
			  long value2 = bvTheory->getUnderApprox(bvID);
			  if(value2>=value){
					throw std::runtime_error("Error in optimization (minimum values are inconsistent with model)");
			  }
			  value=value2;

			  assume.pop();
			  if(opt_verb>=1){
				  printf("\rMin bv%d = %ld",bvID,value);
			  }
		  }else{
			  assume.pop();

			  if(value<last_decision_value){
				  //if that last decrease in value was by more than 1
				  last_decision_lit =  bvTheory->toSolver(bvTheory->newComparison(Comparison::leq,bvID,value,var_Undef,opt_decide_optimization_lits));
				  last_decision_value=value;
			  }
			  assume.push(last_decision_lit);
			  r= S->solve(last_satisfying_assign,false,false);

			  if(!r){
				  throw std::runtime_error("Error in optimization (instance has become unsat)");
			  }
			  for(Lit l:assume){
					if(S->value(l)!=l_True){
						throw std::runtime_error("Error in optimization (model is inconsistent with assumptions)");
					}
				}
			  if(value!= bvTheory->getUnderApprox(bvID)){
				  throw std::runtime_error("Error in optimization (minimum values are inconsistent with model)");
			  }
			  break;
		  }
	  }
	  return value;
}

long optimize_binary(Monosat::SimpSolver * S, Monosat::BVTheorySolver<long> * bvTheory,const vec<Lit> & assumes,int bvID, int conflict_limit, bool & hit_cutoff, long & n_solves){
	vec<Lit> assume;
	for(Lit l:assumes)
		assume.push(l);
	vec<Lit> last_satisfying_assign;
	for(Var v = 0;v<S->nVars();v++){
		if(S->value(v)==l_True){
			last_satisfying_assign.push(mkLit(v));
		}else if(S->value(v)==l_False){
			last_satisfying_assign.push(mkLit(v,true));
		}else{
			//this variable was unassigned.
		}
	}

	  long min_val = bvTheory->getUnderApprox(bvID,true);
	  long max_val = bvTheory->getUnderApprox(bvID);
	  if(opt_verb>=1){
		  printf("Min bv%d = %ld",bvID,max_val);
	  }
	  //might be worth trying 0 and 1 and value-1 first...

	  long last_decision_value=max_val;
	  Lit last_decision_lit =  bvTheory->toSolver(bvTheory->newComparison(Comparison::leq,bvID,max_val,var_Undef,opt_decide_optimization_lits));
	  while(min_val < max_val){
		  long mid_point = min_val + (max_val - min_val) / 2;
		  assert(mid_point>=0);assert(mid_point>=min_val);assert(mid_point<max_val);

		  Lit decision_lit = bvTheory->toSolver(bvTheory->newComparison(Comparison::leq,bvID,mid_point,var_Undef,opt_decide_optimization_lits));
		  assume.push(decision_lit);
		  n_solves++;

		  bool r;

		  if(opt_limit_conflicts>0 || conflict_limit>=0){
			  S->setConfBudget(std::min((int)opt_limit_conflicts,conflict_limit));
			  lbool res = S->solveLimited(assume,false,false);
			  if (res==l_Undef){
				  if(opt_verb>0){
					  printf("\nBudget exceeded during optimization, quiting early (model might not be optimal!)\n");
				  }
				  r=false;
			  }else{
				  r = res==l_True;
			  }
		  }else
			  r = S->solve(assume,false,false);

		  assume.pop();
		  if (r){
			  last_satisfying_assign.clear();
				for(Var v = 0;v<S->nVars();v++){
					if(S->value(v)==l_True){
						last_satisfying_assign.push(mkLit(v));
					}else if(S->value(v)==l_False){
						last_satisfying_assign.push(mkLit(v,true));
					}else{
						//this variable was unassigned.
					}
				}
			  last_decision_lit=decision_lit;
			  last_decision_value=mid_point;
			  long new_value = bvTheory->getUnderApprox(bvID);
			  if(new_value>=max_val){
					throw std::runtime_error("Error in optimization (minimum values are inconsistent with model)");
			  }

			  assert(new_value<=mid_point);
			  assert(new_value<max_val);
			  max_val=new_value;
			  if(new_value<min_val){
				  //this can only happen if a budget was used and the solver quit early.
				  min_val=new_value;
				  assert(min_val>=bvTheory->getUnderApprox(bvID,true));
			  }
			  if(opt_verb>=1){
				  printf("\rMin bv%d = %ld",bvID,max_val);
			  }
		  }else{
			  min_val = mid_point+1;
			  //set the solver back to its last satisfying assignment
			  //this is technically not required, but it should be cheap, and will also reset the solvers decision phase heuristic
			  r= S->solve(last_satisfying_assign,false,false);
			  if(!r){
				  throw std::runtime_error("Error in optimization (instance has become unsat)");
			  }
		  }
	  }


	  if(max_val<last_decision_value){
		  //if that last decrease in value was by more than 1
		  last_decision_lit =  bvTheory->toSolver(bvTheory->newComparison(Comparison::leq,bvID,max_val,var_Undef,opt_decide_optimization_lits));
		  last_decision_value=max_val;
	  }
	  bool r;
	  assume.push(last_decision_lit);
	  r= S->solve(last_satisfying_assign,false,false);

	  if(!r){
		  throw std::runtime_error("Error in optimization (instance has become unsat)");
	  }
	  for(Lit l:assume){
			if(S->value(l)!=l_True){
				throw std::runtime_error("Error in optimization (model is inconsistent with assumptions)");
			}
		}
	  if(max_val!= bvTheory->getUnderApprox(bvID)){
		  throw std::runtime_error("Error in optimization (minimum values are inconsistent with model)");
	  }
	  return max_val;
}
};



#endif /* OPTIMIZE_CPP_ */
