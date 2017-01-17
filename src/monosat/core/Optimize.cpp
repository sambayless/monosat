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


#ifndef OPTIMIZE_CPP_
#define OPTIMIZE_CPP_
#include "monosat/core/Optimize.h"
#include <csignal>
#include <sys/resource.h>
#include <stdexcept>
#include <cstdarg>
#include <string>
#include <cstdint>
#include <limits>
namespace Monosat{

namespace Optimization{
	static long long time_limit=-1;
	static long long memory_limit=-1;

	static bool has_system_time_limit=false;
	static rlim_t system_time_limit;
	static bool has_system_mem_limit=false;
	static rlim_t system_mem_limit;

	static Solver * solver;

	static sighandler_t system_sigxcpu_handler = nullptr;

	//static initializer, following http://stackoverflow.com/a/1681655
	namespace {
	  struct initializer {
		initializer() {
			system_sigxcpu_handler=nullptr;
			time_limit=-1;
			memory_limit=-1;
			has_system_time_limit=false;
			has_system_mem_limit=false;
		}

		~initializer() {
			solver =nullptr;
		}
	  };
	  static initializer i;
	}
	 void disableResourceLimits(Solver * S);
	 static void SIGNAL_HANDLER_api(int signum) {
		if(solver){
			fprintf(stderr,"Monosat resource limit reached\n");
			Solver * s = solver;
			disableResourceLimits(s);
			s->interrupt();
		}
	}


	 void enableResourceLimits(Solver * S){
		if(!solver){
		assert(solver==nullptr);
		solver=S;

		struct rusage ru;
		getrusage(RUSAGE_SELF, &ru);
		__time_t cur_time = ru.ru_utime.tv_sec;

		rlimit rl;
		getrlimit(RLIMIT_CPU, &rl);
		time_limit = opt_limit_optimization_time;
		if(!has_system_time_limit){
			has_system_time_limit=true;
			system_time_limit=rl.rlim_cur;
		}

		if (time_limit < INT32_MAX && time_limit>=0) {
			assert(cur_time>=0);
			long long local_time_limit = 	 time_limit+cur_time;//make this a relative time limit
			if (rl.rlim_max == RLIM_INFINITY || (rlim_t) local_time_limit < rl.rlim_max) {
				rl.rlim_cur = local_time_limit;
				if (setrlimit(RLIMIT_CPU, &rl) == -1)
					fprintf(stderr,"WARNING! Could not set resource limit: CPU-time.\n");
			}
		}else{
			rl.rlim_cur = rl.rlim_max;
			if (setrlimit(RLIMIT_CPU, &rl) == -1)
				fprintf(stderr,"WARNING! Could not set resource limit: CPU-time.\n");
		}

		getrlimit(RLIMIT_AS, &rl);
		if(!has_system_mem_limit){
			has_system_mem_limit=true;
			system_mem_limit=rl.rlim_cur;
		}
		// Set limit on virtual memory:
		if (memory_limit < INT32_MAX && memory_limit>=0) {
			rlim_t new_mem_lim = (rlim_t) memory_limit * 1024 * 1024; //Is this safe?

			if (rl.rlim_max == RLIM_INFINITY || new_mem_lim < rl.rlim_max) {
				rl.rlim_cur = new_mem_lim;
				if (setrlimit(RLIMIT_AS, &rl) == -1)
					fprintf(stderr, "WARNING! Could not set resource limit: Virtual memory.\n");
			}else{
				rl.rlim_cur = rl.rlim_max;
				if (setrlimit(RLIMIT_AS, &rl) == -1)
					fprintf(stderr, "WARNING! Could not set resource limit: Virtual memory.\n");
			}
		}
		sighandler_t old_sigxcpu = signal(SIGXCPU, SIGNAL_HANDLER_api);
		if(old_sigxcpu != SIGNAL_HANDLER_api){
			system_sigxcpu_handler = old_sigxcpu;//store this value for later
		}
		}
	}

	 void disableResourceLimits(Solver * S){
		if(solver){

		solver=nullptr;
		rlimit rl;
		getrlimit(RLIMIT_CPU, &rl);
		if(has_system_time_limit){
			has_system_time_limit=false;
			if (rl.rlim_max == RLIM_INFINITY || (rlim_t)system_time_limit < rl.rlim_max) {
				rl.rlim_cur = system_time_limit;
				if (setrlimit(RLIMIT_CPU, &rl) == -1)
					fprintf(stderr,"WARNING! Could not set resource limit: CPU-time.\n");
			}else{
				rl.rlim_cur = rl.rlim_max;
				if (setrlimit(RLIMIT_CPU, &rl) == -1)
					fprintf(stderr,"WARNING! Could not set resource limit: CPU-time.\n");
			}
		}
		getrlimit(RLIMIT_AS, &rl);
		if(has_system_mem_limit){
			has_system_mem_limit=false;
			if (rl.rlim_max == RLIM_INFINITY || system_mem_limit < rl.rlim_max) {
				rl.rlim_cur = system_mem_limit;
				if (setrlimit(RLIMIT_AS, &rl) == -1)
					fprintf(stderr, "WARNING! Could not set resource limit: Virtual memory.\n");
			}else{
				rl.rlim_cur = rl.rlim_max;
				if (setrlimit(RLIMIT_AS, &rl) == -1)
					fprintf(stderr, "WARNING! Could not set resource limit: Virtual memory.\n");
			}
		}
		if (system_sigxcpu_handler){
			 signal(SIGXCPU, system_sigxcpu_handler);
			 system_sigxcpu_handler=nullptr;
		}
	}
	}
}

inline bool lt(int64_t lhs, int64_t rhs, bool invert){
	if(invert){
		return lhs>rhs;
	}else{
		return lhs<rhs;
	}
}
inline bool gt(int64_t lhs, int64_t rhs, bool invert){
	if(invert){
		return lhs<rhs;
	}else{
		return lhs>rhs;
	}
}
inline bool leq(int64_t lhs, int64_t rhs, bool invert){
	if(invert){
		return lhs>=rhs;
	}else{
		return lhs<=rhs;
	}
}
inline bool geq(int64_t lhs, int64_t rhs, bool invert){
	if(invert){
		return lhs<=rhs;
	}else{
		return lhs>=rhs;
	}
}

int64_t getApprox(Monosat::BVTheorySolver<int64_t> * bvTheory, int bvID,bool overApprox, bool level0=false){
	if(overApprox){
		return bvTheory->getOverApprox(bvID, level0);
	}else{
		return bvTheory->getUnderApprox(bvID, level0);
	}
}

int64_t optimize_linear_bv(Monosat::SimpSolver * S, Monosat::BVTheorySolver<int64_t> * bvTheory, bool invert,const vec<Lit> & assumes,int bvID, bool & hit_cutoff, int64_t & n_solves, bool & found_model){
	found_model=false;
	hit_cutoff=false;
	vec<Lit> assume;
	for(Lit l:assumes)
		assume.push(l);
	vec<Lit> last_satisfying_assign;
	for(Var v = 0;v<S->nVars();v++){
		if(!S->isEliminated(v)) {
			if (S->value(v) == l_True) {
				last_satisfying_assign.push(mkLit(v));
			} else if (S->value(v) == l_False) {
				last_satisfying_assign.push(mkLit(v, true));
			} else {
				//this variable was unassigned.
			}
		}
	}

	int64_t value = getApprox(bvTheory,bvID,!invert);
	int64_t last_decision_value=value;
	if(opt_verb>=1 || opt_verb_optimize>=1){
		printf("Best bv%d = %ld",bvID,value);
	}
	  // int bvID,const Weight & to, Var outerVar = var_Undef, bool decidable=true
	Lit last_decision_lit =  bvTheory->toSolver(bvTheory->newComparison(invert ? Comparison::geq : Comparison::leq,bvID,value,var_Undef,opt_decide_optimization_lits));
	while(gt(value,getApprox(bvTheory,bvID,invert,true),invert) && !hit_cutoff){
		Lit decision_lit = bvTheory->toSolver(bvTheory->newComparison(invert? Comparison::geq : Comparison::leq,bvID, invert?value+1 : value-1,var_Undef,opt_decide_optimization_lits));
		  assume.push(decision_lit);
		  n_solves++;



		  int conflict_limit = S->getConflictBudget();
		  if(conflict_limit<0)
			  conflict_limit=INT32_MAX;
		  int opt_lim = opt_limit_optimization_conflicts;
		  if(opt_lim<=0)
			  opt_lim=INT32_MAX;
		  int limit = std::min(opt_lim,conflict_limit);
		  if(limit>= INT32_MAX){
			limit=-1;//disable limit.
		  }
		  S->setConfBudget(limit);
		  bool r;
		  lbool res = S->solveLimited(assume,false,false);
		  found_model|=(res==l_True);
		  if (res==l_Undef){
			  hit_cutoff=true;
			  if(opt_verb>0){
				  printf("\nBudget exceeded during optimization, quiting early (model might not be optimal!)\n");
			  }
			  r=false;
		  }else{
			  r = res==l_True;
		  }


		  if (r){
			  last_satisfying_assign.clear();
			  for(Var v = 0;v<S->nVars();v++){
				  if(!S->isEliminated(v)) {
					  if (S->value(v) == l_True) {
						  last_satisfying_assign.push(mkLit(v));
					  } else if (S->value(v) == l_False) {
						  last_satisfying_assign.push(mkLit(v, true));
					  } else {
						  //this variable was unassigned.
					  }
				  }
			  }
			  last_decision_lit=decision_lit;
			  last_decision_value= invert ? value+1 : value-1;
			  if(S->value(decision_lit)!=l_True){
				  throw std::runtime_error("Error in optimization (comparison not enforced)");
			  }
			  for(Lit l:assume){
					if(S->value(l)!=l_True){
						throw std::runtime_error("Error in optimization (model is inconsistent with assumptions)");
					}
				}
			  int64_t value2 =getApprox(bvTheory, bvID, !invert); //bvTheory->getOverApprox(bvID);
			  if(geq(value2,value,invert)){
					throw std::runtime_error("Error in optimization (minimum values are inconsistent with model)");
			  }
			  value=value2;

			  assume.pop();
			  if(opt_verb>=1){
				  printf("\rMin bv%d = %ld",bvID,value);
			  }
		  }else{
			  assume.pop();

			if(lt(value,last_decision_value,invert)){
				  //if that last decrease in value was by more than 1
				last_decision_lit =  bvTheory->toSolver(bvTheory->newComparison(invert ? Comparison::geq : Comparison::leq,bvID,value,var_Undef,opt_decide_optimization_lits));
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
			if(lt(value, getApprox(bvTheory, bvID, !invert),invert)){//bvTheory->getOverApprox(bvID)){
				  throw std::runtime_error("Error in optimization (minimum values are inconsistent with model)");
			  }
			  //if(bvTheory->getOverApprox(bvID) < value){
            if(lt(getApprox(bvTheory, bvID, !invert),value,invert)){
				value = getApprox(bvTheory, bvID, !invert);//bvTheory->getOverApprox(bvID);
			  }
			  break;
		  }
	  }
	  return value;
}
int evalPB(SimpSolver & S,const Objective & o, bool over_approx, bool eval_at_level_0=false) {
	int sum_val = 0;
	assert(o.isPB());
	for (int i = 0; i < o.pb_lits.size(); i++) {
		Lit l = o.pb_lits[i];
		if (l == lit_Undef)
			continue;
		int weight = 1;
		if (i < o.pb_weights.size()) {
			weight = o.pb_weights[i];
		}
		lbool val = S.value(l);
		if (eval_at_level_0 && S.level(var(l)) > 0) {
			val = l_Undef;
		}
		if (val == l_True) {
			sum_val += weight;
		} else if (val == l_Undef) {
			if (over_approx) {
				sum_val += weight;
			}
		}
	}
	return sum_val;
}
	int64_t optimize_linear_pb(Monosat::SimpSolver * S, PB::PBConstraintSolver * pbSolver, bool invert, const vec<Lit> & assumes, const Objective & o , bool & hit_cutoff, int64_t & n_solves, bool & found_model){
		hit_cutoff=false;
		found_model=false;
		vec<Lit> discarded_pb_constraints;
		vec<Lit> tmp_clause;
		vec<Lit> assume;
		for(Lit l:assumes)
			assume.push(l);
		vec<Lit> last_satisfying_assign;
		for(Var v = 0;v<S->nVars();v++){
			if(!S->isEliminated(v)) {
				if (S->value(v) == l_True) {
					last_satisfying_assign.push(mkLit(v));
				} else if (S->value(v) == l_False) {
					last_satisfying_assign.push(mkLit(v, true));
				} else {
					//this variable was unassigned.
				}
			}
		}

		int value = evalPB(*S, o,!invert);
		int last_decision_value=value;

		// int bvID,const Weight & to, Var outerVar = var_Undef, bool decidable=true
		//Lit last_decision_lit =  lit_Undef;// pbSolver->addConditionalConstr(o.pb_lits, o.pb_weights,value, invert ? PB::Ineq::GEQ : PB::Ineq::LEQ);
		while(gt(value,evalPB(*S,o,invert,true),invert ) && !hit_cutoff){
			Lit decision_lit =  pbSolver->addConditionalConstr(o.pb_lits, o.pb_weights,invert ? value+1 : value-1, invert ? PB::Ineq::GEQ : PB::Ineq::LEQ);

			assume.push(decision_lit);
			n_solves++;



			int conflict_limit = S->getConflictBudget();
			if(conflict_limit<0)
				conflict_limit=INT32_MAX;
			int opt_lim = opt_limit_optimization_conflicts;
			if(opt_lim<=0)
				opt_lim=INT32_MAX;
			int limit = std::min(opt_lim,conflict_limit);
			if(limit>= INT32_MAX){
				limit=-1;//disable limit.
			}
			S->setConfBudget(limit);
			bool r;
			for(Lit l:discarded_pb_constraints){
				assume.push(~l);
			}
			lbool res = S->solveLimited(assume,false,false);
			found_model|=(res==l_True);
			if (res==l_Undef){
				hit_cutoff=true;
				if(opt_verb>0){
					printf("\nBudget exceeded during optimization, quiting early (model might not be optimal!)\n");
				}
				r=false;
			}else{
				r = res==l_True;
			}
			assume.shrink(discarded_pb_constraints.size()+1);

			if (r) {
				while(discarded_pb_constraints.size()){
					Lit l = discarded_pb_constraints.last();
					tmp_clause.clear();
					assert(S->value(l)==l_False);
					tmp_clause.push(~l);
					S->addClauseSafely(tmp_clause);
					discarded_pb_constraints.pop();
				}
				last_satisfying_assign.clear();
				for (Var v = 0; v < S->nVars(); v++) {
					if (!S->isEliminated(v) && v != var(decision_lit)) {
						if (S->value(v) == l_True) {
							last_satisfying_assign.push(mkLit(v));
						} else if (S->value(v) == l_False) {
							last_satisfying_assign.push(mkLit(v, true));
						} else {
							//this variable was unassigned.
						}
					}
				}

				last_decision_value = invert ? value + 1 : value - 1;
				//last_satisfying_assign.push(decision_lit);
				if (S->value(decision_lit) != l_True) {
					throw std::runtime_error("Error in optimization (comparison not enforced)");
				}
				for (Lit l:assume) {
					if (S->value(l) != l_True) {
						throw std::runtime_error("Error in optimization (model is inconsistent with assumptions)");
					}
				}
				int64_t value2 = evalPB(*S, o, !invert);//bvTheory->getOverApprox(bvID);
				//if(value2>=value){
				if (geq(value2, value, invert)) {
					throw std::runtime_error("Error in optimization (minimum values are inconsistent with model)");
				}
				value = value2;



			}else{


				if(lt(value,last_decision_value,invert) ){
					//if that last decrease in value was by more than 1
					//if(last_decision_lit!=lit_Undef)
					//    S->addClause(~last_decision_lit);
					//last_decision_lit = pbSolver->addConditionalConstr(o.pb_lits, o.pb_weights,value,invert ? PB::Ineq::GEQ : PB::Ineq::LEQ);
					//bvTheory->toSolver(bvTheory->newComparison(Comparison::leq,bvID,value,var_Undef,opt_decide_optimization_lits));
					last_decision_value=value;
					//last_satisfying_assign.push(last_decision_lit);
				}
				//assume.push(last_decision_lit);
				r= S->solve(last_satisfying_assign,false,false);

				if(!r){
					throw std::runtime_error("Error in optimization (instance has become unsat)");
				}
				for(Lit l:assume){
					if(S->value(l)!=l_True){
						throw std::runtime_error("Error in optimization (model is inconsistent with assumptions)");
					}
				}
				int over = evalPB(*S, o, !invert);
				if(lt(value, over,invert)){
					throw std::runtime_error("Error in optimization (minimum values are inconsistent with model)");
				}
				if(lt(over, value,invert)){
					value = over;
				}
				break;
			}
			discarded_pb_constraints.push(decision_lit);
		}
		while(discarded_pb_constraints.size()){
			Lit l = discarded_pb_constraints.last();
			tmp_clause.clear();
			tmp_clause.push(~l);
			S->addClauseSafely(tmp_clause);
			discarded_pb_constraints.pop();
		}
		return value;
	}


	int64_t optimize_binary_pb(Monosat::SimpSolver * S,  PB::PBConstraintSolver * pbSolver, bool invert, const vec<Lit> & assumes,const Objective & o,  bool & hit_cutoff, int64_t & n_solves, bool & found_model){
		hit_cutoff=false;
		found_model=false;
		vec<Lit> discarded_pb_constraints;
		vec<Lit> assume;
		vec<Lit> tmp_clause;
		for(Lit l:assumes)
			assume.push(l);
		vec<Lit> last_satisfying_assign;
		for(Var v = 0;v<S->nVars();v++){
			if(!S->isEliminated(v)) {
				if (S->value(v) == l_True) {
					last_satisfying_assign.push(mkLit(v));
				} else if (S->value(v) == l_False) {
					last_satisfying_assign.push(mkLit(v, true));
				} else {
					//this variable was unassigned.
				}
			}
		}

		int min_val = evalPB(*S, o, invert, true);// bvTheory->getUnderApprox(bvID,true);
		int max_val = evalPB(*S, o, !invert); //bvTheory->getOverApprox(bvID);

		int suggested_next_midpoint = -1;

		//try the minimum possible value first, just in case we get lucky.
		//it might also be a good idea to try min_val+1, and max_val-1.

		int underapprox_sat_val =  evalPB(*S, o, invert);
		if(lt(underapprox_sat_val,max_val,invert))
			suggested_next_midpoint =	underapprox_sat_val;

		bool first_round=true;
		int64_t last_decision_value=max_val;
		//Lit last_decision_lit =  lit_Undef;//  pbSolver->addConditionalConstr(o.pb_lits, o.pb_weights,max_val, PB::Ineq::LEQ); //bvTheory->toSolver(bvTheory->newComparison(Comparison::leq,bvID,max_val,var_Undef,opt_decide_optimization_lits));
		while(lt(min_val, max_val,invert) && !hit_cutoff){
			int64_t mid_point = min_val + (max_val - min_val) / 2;

			if (geq(mid_point, max_val,invert))
				mid_point = invert ? max_val+1 : max_val - 1;//can this ever happen?

			if(suggested_next_midpoint != -1 && geq(suggested_next_midpoint,min_val,invert) && lt(suggested_next_midpoint,mid_point,invert))
				mid_point =	suggested_next_midpoint;
			assert(mid_point>=0);assert(geq(mid_point,min_val,invert));assert(lt(mid_point,max_val,invert));
			Lit decision_lit =  pbSolver->addConditionalConstr(o.pb_lits, o.pb_weights,mid_point, invert? PB::Ineq::GEQ : PB::Ineq::LEQ);
			//bvTheory->toSolver(bvTheory->newComparison(Comparison::leq,bvID,mid_point,var_Undef,opt_decide_optimization_lits));
			assume.push(decision_lit);
			n_solves++;
			for(Lit l:discarded_pb_constraints){
				assume.push(~l);
			}
			{
				int conflict_limit = S->getConflictBudget();
				if(conflict_limit<0)
					conflict_limit=INT32_MAX;
				int opt_lim = opt_limit_optimization_conflicts;
				if(opt_lim<=0)
					opt_lim=INT32_MAX;
				int limit = std::min(opt_lim,conflict_limit);
				if(limit>= INT32_MAX){
					limit=-1;//disable limit.
				}
				S->setConfBudget(limit);
			}

			Optimization::enableResourceLimits(S);
			bool r;

			lbool res = S->solveLimited(assume,false,false);
			found_model|=(res==l_True);
			Optimization::disableResourceLimits(S);
			if (res==l_Undef){
				hit_cutoff=true;
				if(opt_verb>0){
					printf("\nBudget exceeded during optimization, quiting early (model might not be optimal!)\n");
				}
				r=false;
			}else{
				r = res==l_True;
			}

			assume.shrink(discarded_pb_constraints.size()+1);
			if (r){
				while(discarded_pb_constraints.size()){
					Lit l = discarded_pb_constraints.last();
					tmp_clause.clear();
					tmp_clause.push(~l);
					S->addClauseSafely(tmp_clause);//safe to add this clause at next restart, because it is guaranteed to match
					//the assignment in the solver currently, and hence to be compatible with last_satisfying_assign.
					discarded_pb_constraints.pop();
				}
				last_satisfying_assign.clear();
				for(Var v = 0;v<S->nVars();v++){
					if(!S->isEliminated(v)   && v != var(decision_lit) ) {
						if (S->value(v) == l_True) {
							last_satisfying_assign.push(mkLit(v));
						} else if (S->value(v) == l_False) {
							last_satisfying_assign.push(mkLit(v, true));
						} else {
							//this variable was unassigned.
						}
					}
				}
				//Lit old_dec_lit = last_decision_lit;

				//last_decision_lit=decision_lit;
				//last_satisfying_assign.push(decision_lit);
				last_decision_value=mid_point;
				int new_value = evalPB(*S, o, !invert);//bvTheory->getOverApprox(bvID);
				if(geq(new_value,max_val,invert)){
					throw std::runtime_error("Error in optimization (minimum values are inconsistent with model)");
				}
				int64_t underapprox_sat_val =  evalPB(*S, o, invert); //bvTheory->getUnderApprox(bvID);
				if(lt(underapprox_sat_val,max_val,invert))
					suggested_next_midpoint =	underapprox_sat_val;
				assert(leq(new_value,mid_point,invert));
				assert(lt(new_value,max_val,invert));
				max_val=new_value;
				if(lt(new_value,min_val,invert)){
					//this can only happen if a budget was used and the solver quit early.
					min_val=new_value;
					assert(geq(min_val, evalPB(*S, o, invert, true),invert));
				}
//            if (old_dec_lit != lit_Undef) {
//                discarded_pb_constraints.push(~old_dec_lit);
//                //S->addClause(~old_dec_lit); //why isn't this safe? Because it might imply different values for intermediate variable assignments in the PB constraint, which were recorded in last_satisfying_assign
//            }
			}else{

				if(invert){
					min_val = mid_point - 1;
				}else {
					min_val = mid_point + 1;
				}
				//set the solver back to its last satisfying assignment
				//this is technically not required, but it should be cheap, and will also reset the solvers decision phase heuristic
				r= S->solve(last_satisfying_assign,false,false);
				if(!r){
					throw std::runtime_error("Error in optimization (instance has become unsat)");
				}

			}
			discarded_pb_constraints.push(decision_lit);
		}


		if(lt(max_val,last_decision_value,invert) ){
			//if(last_decision_lit!=lit_Undef)
			//    S->addClause(~last_decision_lit);
			//if that last decrease in value was by more than 1
			//last_decision_lit = pbSolver->addConditionalConstr(o.pb_lits, o.pb_weights,max_val, invert ? PB::Ineq::GEQ : PB::Ineq::LEQ);
			//bvTheory->toSolver(bvTheory->newComparison(Comparison::leq,bvID,max_val,var_Undef,opt_decide_optimization_lits));
			last_decision_value=max_val;
			//last_satisfying_assign.push(last_decision_lit);
		}
		bool r;
		//assume.push(last_decision_lit);
		r= S->solve(last_satisfying_assign,false,false);
		if(!r){
			throw std::runtime_error("Error in optimization (instance has become unsat)");
		}
		for(Lit l:assume){
			if(S->value(l)!=l_True){
				throw std::runtime_error("Error in optimization (model is inconsistent with assumptions)");
			}
		}
		int over = evalPB(*S, o, !invert);
		if(lt(max_val,over,invert)){
			throw std::runtime_error("Error in optimization (minimum values are inconsistent with model)");
		}
		if(lt(over, max_val,invert)){
			max_val = over;
		}
		while(discarded_pb_constraints.size()){
			Lit l = discarded_pb_constraints.last();
			tmp_clause.clear();
			tmp_clause.push(~l);
			S->addClauseSafely(tmp_clause);
			discarded_pb_constraints.pop();
		}
		return max_val;
	}




int64_t optimize_binary_bv(Monosat::SimpSolver * S, Monosat::BVTheorySolver<int64_t> * bvTheory, bool maximize,const vec<Lit> & assumes,int bvID,  bool & hit_cutoff, int64_t & n_solves, bool & found_model){
	found_model=false;
	hit_cutoff=false;
	vec<Lit> assume;
	for(Lit l:assumes)
		assume.push(l);
	vec<Lit> last_satisfying_assign;
	for(Var v = 0;v<S->nVars();v++){
		if(!S->isEliminated(v)) {
			if (S->value(v) == l_True) {
				last_satisfying_assign.push(mkLit(v));
			} else if (S->value(v) == l_False) {
				last_satisfying_assign.push(mkLit(v, true));
			} else {
				//this variable was unassigned.
			}
		}
	}

	bool invert=maximize;


	int64_t min_val = getApprox(bvTheory, bvID, invert,true); //bvTheory->getUnderApprox(bvID,true);
	int64_t max_val = getApprox(bvTheory, bvID, !invert); //bvTheory->getOverApprox(bvID);
	  if(opt_verb>=1 || opt_verb_optimize>=1){
		  printf("Best bv%d = %ld",bvID,max_val);
	  }
	int64_t suggested_next_midpoint = -1;

		  //try the minimum possible value first, just in case we get lucky.
		  //it might also be a good idea to try min_val+1, and max_val-1.

	int64_t underapprox_sat_val = getApprox(bvTheory, bvID, invert); //bvTheory->getUnderApprox(bvID);
	int64_t overapprox_sat_val = getApprox(bvTheory, bvID, !invert); //bvTheory->getUnderApprox(bvID);
	if(!maximize) {
		if (lt(underapprox_sat_val, max_val, invert))
			suggested_next_midpoint = underapprox_sat_val;
	}else{
		if (lt(overapprox_sat_val, max_val, invert))
			suggested_next_midpoint = overapprox_sat_val;
	}
	  bool first_round=true;
	int64_t last_decision_value=max_val;
	Lit last_decision_lit =  bvTheory->toSolver(bvTheory->newComparison(invert ? Comparison::geq : Comparison::leq,bvID,max_val,var_Undef,opt_decide_optimization_lits));
	  while(lt(min_val , max_val,invert) && !hit_cutoff){
		int64_t mid_point;
		  if(!invert)
			  mid_point=min_val + (max_val - min_val) / 2;
		  else{
			  mid_point=max_val + (min_val - max_val) / 2;
			  if(mid_point==max_val){//extra check required because integer division rounds down
				  mid_point=min_val;
			  }
		  }
		if(suggested_next_midpoint !=-1 && geq(suggested_next_midpoint,min_val,invert) && lt(suggested_next_midpoint,mid_point,invert))
			  mid_point =	suggested_next_midpoint;
		assert(mid_point>=0);assert(geq(mid_point,min_val,invert));assert(lt(mid_point,max_val,invert));

		Lit decision_lit = bvTheory->toSolver(bvTheory->newComparison(invert? Comparison::geq : Comparison::leq,bvID,mid_point,var_Undef,opt_decide_optimization_lits));
		  assume.push(decision_lit);
		  n_solves++;

		  {
		  int conflict_limit = S->getConflictBudget();
		  if(conflict_limit<0)
			  conflict_limit=INT32_MAX;
		  int opt_lim = opt_limit_optimization_conflicts;
		  if(opt_lim<=0)
			  opt_lim=INT32_MAX;
		  int limit = std::min(opt_lim,conflict_limit);
		  if(limit>= INT32_MAX){
			limit=-1;//disable limit.
		  }
		  S->setConfBudget(limit);
		  }

		  Optimization::enableResourceLimits(S);
		  bool r;
		  lbool res = S->solveLimited(assume,false,false);
		  found_model|=(res==l_True);
		  Optimization::disableResourceLimits(S);
		  if (res==l_Undef){
			  hit_cutoff=true;
			  if(opt_verb>0|| opt_verb_optimize>=1){
				  printf("\nBudget exceeded during optimization, quiting early (model might not be optimal!)\n");
			  }
			  r=false;
		  }else{
			  r = res==l_True;
		  }

		  assume.pop();
		  if (r){
			  last_satisfying_assign.clear();
			  for(Var v = 0;v<S->nVars();v++){
				  if(!S->isEliminated(v)) {
					  if(S->value(v)==l_True){
						  last_satisfying_assign.push(mkLit(v));
					  }else if(S->value(v)==l_False){
						  last_satisfying_assign.push(mkLit(v,true));
					  }else {
						  //this variable was unassigned.
					  }
				  }
			  }
			  last_decision_lit=decision_lit;
			  last_decision_value=mid_point;
			int64_t new_value = getApprox(bvTheory, bvID, !invert); //bvTheory->getOverApprox(bvID);
			if(geq(new_value,max_val,invert)){
					throw std::runtime_error("Error in optimization (minimum values are inconsistent with model)");
			  }
			int64_t underapprox_sat_val =  getApprox(bvTheory, bvID, invert); //bvTheory->getUnderApprox(bvID);
			if(lt(underapprox_sat_val,max_val,invert))
				  suggested_next_midpoint =	underapprox_sat_val;
			assert(leq(new_value,mid_point,invert));
			assert(lt(new_value,max_val,invert));
			  max_val=new_value;
			if(leq(new_value,min_val,invert)){
				  //this can only happen if a budget was used and the solver quit early.
				  min_val=new_value;
                assert(gt(min_val,getApprox(bvTheory,bvID,invert,true),invert)); //assert(min_val>=bvTheory->getUnderApprox(bvID,true));
			  }
			  if(opt_verb>=1 || opt_verb_optimize>=1){
				  printf("\rBest bv%d = %ld",bvID,max_val);
			  }
		  }else{
			min_val = invert ? mid_point-1 : mid_point+1; //yes this is intentionally backward

			  //set the solver back to its last satisfying assignment
			  //this is technically not required, but it should be cheap, and will also reset the solvers decision phase heuristic
			  r= S->solve(last_satisfying_assign,false,false);
			  if(!r){
				  throw std::runtime_error("Error in optimization (instance has become unsat)");
			  }

		  }
	  }


	   if(lt(max_val,last_decision_value,invert)){
		  //if that last decrease in value was by more than 1
		  last_decision_lit =  bvTheory->toSolver(bvTheory->newComparison(invert ? Comparison::geq : Comparison::leq,bvID,max_val,var_Undef,opt_decide_optimization_lits));
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
	if(lt(max_val , getApprox(bvTheory,bvID,!invert),invert)){//bvTheory->getOverApprox(bvID)){
		  throw std::runtime_error("Error in optimization (minimum values are inconsistent with model)");
	  }
	//if(bvTheory->getOverApprox(bvID) < max_val){
    if(lt(getApprox(bvTheory,bvID,!invert),max_val,invert)){
		max_val = getApprox(bvTheory,bvID,!invert); //bvTheory->getOverApprox(bvID);
	  }
	  return max_val;
}

void copyModel(SimpSolver & S, vec<Lit> & dest){
    dest.clear();
    for(Var v = 0;v<S.nVars();v++){
        if(!S.isEliminated(v)  ) {
            if (S.value(v) == l_True) {
                dest.push(mkLit(v));
            } else if (S.value(v) == l_False) {
                dest.push(mkLit(v, true));
            } else {
                //this variable was unassigned.
            }
        }
    }
}

	lbool optimize_and_solve(SimpSolver & S,const vec<Lit> & assumes,const vec<Objective> & objectives,bool do_simp,  bool & found_optimal){
		vec<Lit> best_model;
		vec<Lit> assume;
		vec<int64_t> model_vals;
		for(Lit l:assumes)
			assume.push(l);
		static int solve_runs=0;
		found_optimal=true;
		solve_runs++;
		if(opt_verb>=1 || opt_verb_optimize>=1){
			if(solve_runs>1){
				printf("Solving(%d)...\n",solve_runs);
			}else{
				printf("Solving...\n");
			}
			fflush(stdout);
		}

		if (opt_verb > 1 && assume.size()) {
			printf("Assumptions: ");
			for (int i = 0; i < assume.size(); i++) {
				Lit l = assume[i];
				printf("%d, ", dimacs(l));
			}
			printf("\n");
		}
		if(!objectives.size()){

			return S.solveLimited(assume,opt_pre && do_simp, !opt_pre && do_simp);

		}else{
			bool any_pb=false;
			bool any_bv = false;
			for(Objective & o:objectives){
				any_pb|= (o.isPB() && o.pb_lits.size()>0 ); //don't need to create the pb theory solver if there is only 1 lit
				any_bv|= o.isBV();
			}

			if(any_bv && !S.getBVTheory()){
				throw std::runtime_error("No bitvector theory created (call initBVTheory())!");
			}
			if(any_pb && !S.getPB()){
				throw std::runtime_error("No pb solver created!");
			}

			bool r;
			if (opt_verb >= 2 || opt_verb_optimize>=2) {
				printf("Performing initial solve (before attempting optimization)...\n");
			}
			/*if(opt_begin_optimization_at_minimum ){
				bool r = S.propagateAssignment(assume);//apply initial propagation to get minimal bounds
				if(r) {
					for (int i = 0; i < objectives.size(); i++) {
						if (objectives[i].isBV()) {
							if (!objectives[i].maximize) {
								min_values[i] = bvTheory->getOverApprox(objectives[i].bvID, true);
								max_values[i] = bvTheory->getUnderApprox(objectives[i].bvID, true);
							} else {
								min_values[i] = bvTheory->getUnderApprox(objectives[i].bvID, true);
								max_values[i] = bvTheory->getOverApprox(objectives[i].bvID, true);
							}

						} else {
							if (!objectives[i].maximize) {
								min_values[i] = evalPB(S, objectives[i], true, true);
								max_values[i] = evalPB(S, objectives[i], false, true);
							} else {
								min_values[i] = evalPB(S, objectives[i], false, true);
								max_values[i] = evalPB(S, objectives[i], true, true);
							}
						}
					}
				}
			}else */
			r=true;
			bool ever_solved=false;
			if(opt_optimization_init_solve){
				//problem: this initial solve, in which the optimization function is free, can sometimes be much more expensive than the subsequent optimizaiton steps...
				lbool res = S.solveLimited(assume, opt_pre && do_simp, !opt_pre && do_simp);
				if (res == l_True) {
					r = true;
					ever_solved=true;
					copyModel(S,best_model);
				} else if (res == l_False) {
					r = false;
				} else {
					found_optimal = false;
					return l_Undef;
				}
			}else{
				r = S.propagateAssignment(assume);//apply initial propagation to get minimal bounds
			}
			if(r && objectives.size()){

				for(Lit l:assume){
					if(S.value(l)!=l_True){
						throw std::runtime_error("Error in optimization (model is inconsistent with assumptions)");
					}
				}
				if (opt_verb >= 2 || opt_verb_optimize>=2) {
					printf("Begining optimization...\n");
				}
				Monosat::BVTheorySolver<int64_t> * bvTheory = (Monosat::BVTheorySolver<int64_t> *) S.getBVTheory();
				Monosat::PB::PBConstraintSolver * pbSolver = S.getPB();
				vec<int64_t> min_values;
				vec<int64_t> max_values;
				min_values.growTo(objectives.size());
				max_values.growTo(objectives.size());

				for(int i = 0;i<objectives.size();i++){
					if(objectives[i].isBV()){
						if(!objectives[i].maximize){
							min_values[i]=bvTheory->getOverApprox(objectives[i].bvID,true);
							max_values[i]=bvTheory->getUnderApprox(objectives[i].bvID,true);
						}else{
							min_values[i]=bvTheory->getUnderApprox(objectives[i].bvID,true);
							max_values[i]=bvTheory->getOverApprox(objectives[i].bvID,true);
						}

					}else{
						if(!objectives[i].maximize){
							min_values[i]=evalPB(S,objectives[i],true,true);
							max_values[i]=evalPB(S, objectives[i],false,true);
						}else{
							min_values[i]=evalPB(S,objectives[i],false,true);
							max_values[i]=evalPB(S, objectives[i],true,true);
						}
					}
				}


				int64_t n_solves =  1;
				bool hit_cutoff=false;
				for (int i = 0;i<objectives.size() && !hit_cutoff;i++){
					if(objectives[i].isBV()) {
						int bvID = objectives[i].bvID;

						if (opt_verb >= 1 || opt_verb_optimize>=1) {
							printf("%s bv%d (%d of %d)\n", objectives[i].maximize ? "Maximizing": "Minimizing", bvID, i + 1, objectives.size());
						}
						int64_t val=0;
						if (!opt_binary_search_optimization_bv) {
							val = optimize_linear_bv(&S, bvTheory,objectives[i].maximize, assume, bvID, hit_cutoff, n_solves, found_model);

						} else {
							val = optimize_binary_bv(&S, bvTheory, objectives[i].maximize,assume, bvID, hit_cutoff, n_solves, found_model);
						}
						if(objectives[i].maximize){
							assert(val>=min_values[i]);assert(val<=max_values[i]);
							max_values[i] = val;
						}else{
							assert(val>=min_values[i]);assert(val<=max_values[i]);
							min_values[i] = val;
						}
						if (hit_cutoff) {
							found_optimal = false;
						}
						if (opt_limit_optimization_time_per_arg)
							hit_cutoff = false;//keep trying to minimize subsequent arguments
						if(ever_solved)
							copyModel(S,best_model);
						if(objectives[i].maximize){
							assume.push(bvTheory->toSolver(
									bvTheory->newComparison(Comparison::geq, bvID, max_values[i], var_Undef,
															opt_decide_optimization_lits)));
						}else {
							assume.push(bvTheory->toSolver(
									bvTheory->newComparison(Comparison::leq, bvID, min_values[i], var_Undef,
															opt_decide_optimization_lits)));
						}
						assert(min_values[i] >= 0);
						if (opt_verb >= 1 || opt_verb_optimize>=1) {
							printf("\rBest bv%d = %ld\n", bvID, objectives[i].maximize? max_values[i]: min_values[i]);
						}

					}else{

						if (opt_verb >= 1 || opt_verb_optimize>=1) {
							printf("%s pb (%d of %d)\n",objectives[i].maximize ? "Maximizing": "Minimizing", i + 1, objectives.size());
						}
						int64_t val=0;
						if (!opt_binary_search_optimization_pb) {
							val = optimize_linear_pb(&S, pbSolver,objectives[i].maximize, assume, objectives[i], hit_cutoff, n_solves, found_model);
						} else {
							val = optimize_binary_pb(&S, pbSolver,objectives[i].maximize, assume,  objectives[i], hit_cutoff, n_solves, found_model);
						}
						if(objectives[i].maximize){
							assert(val>=max_values[i]);
							max_values[i] = val;
						}else{
							assert(val<=min_values[i]);
							min_values[i] = val;
						}

						if (hit_cutoff) {
							found_optimal = false;
						}

						if (opt_limit_optimization_time_per_arg)
							hit_cutoff = false;//keep trying to minimize subsequent arguments
						if(objectives[i].maximize){
							assume.push(pbSolver->addConditionalConstr(objectives[i].pb_lits, objectives[i].pb_weights,
																	   max_values[i], PB::Ineq::GEQ));
						}else {
							assume.push(pbSolver->addConditionalConstr(objectives[i].pb_lits, objectives[i].pb_weights,
																	   min_values[i], PB::Ineq::LEQ));
						}
						//for bitvectors, we don't need to resolve.
						//however, for pb constraints, because we add clauses to the model to discard irelevant PB constraint
						//after finding a best solution, the solver must re-solve in order to keep the model consistent
						//(due to changing values for intermediate values in the PB constraints).
						if(ever_solved && !S.solve(assume)){
							throw std::runtime_error("Error in optimization (best values are inconsistent with model)");
						}
						if(ever_solved)
							copyModel(S,best_model);

						assert(min_values[i] >= 0);

					}
					//enforce that this bitvector stays at the best value that was found for it

					model_vals.clear();
					for(int j = 0;j<i;j++){
						if(objectives[j].isBV()) {
							int bvID = objectives[j].bvID;
							int64_t model_val = bvTheory->getOverApprox(bvID);
							model_vals.push(model_val);
						}else{
							int64_t model_val = evalPB(S,objectives[j],true);
							model_vals.push(model_val);
						}
					}

					for(int j = 0;j<i;j++){
						if(objectives[j].isBV()) {
							int bvID = objectives[j].bvID;
							int64_t best_value = objectives[i].maximize ? max_values[j] : min_values[j];
							int64_t model_val = model_vals[j];
							if (lt(best_value, model_val,objectives[i].maximize)) {
								throw std::runtime_error(
										"Error in optimization (best values are inconsistent with model for bv " +
										std::to_string(j) + " (bvid " + std::to_string(bvID) + " ): expected value <= " +
										std::to_string(best_value) + ", found " + std::to_string(model_val) + ")");
							} else if (lt(model_val, best_value,objectives[j].maximize)) {
								//if the best known value for any earlier bitvector, (which can happen if optimization is aborted early),
								//is found, enforce that this improved value must be kept in the future.
								Lit decision_lit = bvTheory->toSolver(
										bvTheory->newComparison(objectives[j].maximize ? Comparison::geq : Comparison::leq, bvID, model_val, var_Undef,
																opt_decide_optimization_lits));
								assume.push(decision_lit);
								if(objectives[j].maximize){
									max_values[j] = model_val;
								}else {
									min_values[j] = model_val;
								}
							}
						}else{
							int64_t best_value = objectives[i].maximize ? max_values[j] : min_values[j];
							int64_t model_val = model_vals[j];
							if (lt(best_value , model_val,objectives[i].maximize)) {
								throw std::runtime_error(
										"Error in optimization (best values are inconsistent with model for pb " +
										std::to_string(j) + " ): expected value <= " +
										std::to_string(best_value) + ", found " + std::to_string(model_val) + ")");
							} else if (lt(model_val , best_value,objectives[i].maximize)) {
								//if the best known value for any earlier bitvector, (which can happen if optimization is aborted early),
								//is found, enforce that this improved value must be kept in the future.
								Lit decision_lit = pbSolver->addConditionalConstr(objectives[j].pb_lits, objectives[j].pb_weights,model_val, objectives[i].maximize ? PB::Ineq::GEQ: PB::Ineq::LEQ);
								assume.push(decision_lit);
								if(objectives[j].maximize){
									max_values[j] = model_val;
								}else {
									min_values[j] = model_val;
								}
							}
						}
					}
				}
				assert(r);

				if(opt_verb>0 || opt_verb_optimize>=1){
					printf("Best values found (after %ld calls) : ",n_solves);
					for(int i = 0;i<min_values.size();i++){
						int64_t best_value = objectives[i].maximize ? max_values[i] : min_values[i];
						if(objectives[i].isBV()) {
							int bvID = objectives[i].bvID;
							printf("bv%d=%ld,", bvID, best_value);
						}else{

							printf("pb (");
							for(int j = 0;j<objectives[i].pb_lits.size();j++){
								Lit l = objectives[i].pb_lits[j];
								int weight = objectives[i].pb_weights[j];
								printf("%d * %dL +",weight,dimacs(l));
							}
							printf(") =%ld,", best_value);
						}
					}
					printf("\n");
				}
				if(opt_check_solution){
					if (!S.solve(best_model)){
						throw std::runtime_error("Error in optimization (best values are inconsistent with model)");
					}
					for(int i = 0;i<objectives.size();i++){
						if(objectives[i].isBV()) {
							int bvID = objectives[i].bvID;
							int64_t best_value = objectives[i].maximize ? max_values[i] : min_values[i];
							int64_t model_val = getApprox(bvTheory,bvID,!objectives[i].maximize);
							if (lt(best_value , model_val,objectives[i].maximize)) {
								throw std::runtime_error(
										"Error in optimization (best values are inconsistent with model)");
							}
						}else{
							int64_t best_value = objectives[i].maximize ? max_values[i] : min_values[i];
							int64_t model_val =evalPB(S, objectives[i],!objectives[i].maximize);
							if (lt(best_value , model_val,objectives[i].maximize)) {
								throw std::runtime_error(
										"Error in optimization (best values are inconsistent with model)");
							}
						}
					}
				}
			}

			return r? l_True:l_False;
		}


	}
};
#endif /* OPTIMIZE_CPP_ */
