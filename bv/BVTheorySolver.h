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

#ifndef BV_THEORY_SOLVER_H_
#define BV_THEORY_SOLVER_H_
#include <cstddef>
#include <gmpxx.h>
#include "utils/System.h"
#include "core/Theory.h"

#include "core/SolverTypes.h"
#include "mtl/Map.h"

#include "BVTheory.h"
#include "utils/System.h"
#include "core/TheorySolver.h"
#include <algorithm>
#include <vector>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <sstream>
#include <iostream>
#include <exception>

template<typename Weight>
class BVTheorySolver;
#include "core/Remap.h"
namespace Monosat {
using std::min;
using std::max;
enum class Comparison{
		lt,leq,gt,geq,none
};
inline Comparison operator ~(Comparison p) {
	switch (p){
		case Comparison::lt:
			return Comparison::gt;
		case Comparison::leq:
			return Comparison::geq;
		case Comparison::gt:
			return Comparison::lt;
		default:
			return Comparison::leq;
	}
}
inline Comparison operator -(Comparison p) {
	switch (p){
		case Comparison::lt:
			return Comparison::geq;
		case Comparison::leq:
			return Comparison::gt;
		case Comparison::gt:
			return Comparison::leq;
		default:
			return Comparison::lt;
	}
}
inline std::ostream& operator<<(std::ostream &out, const Comparison &p) {
	switch (p){
		case Comparison::lt:
			return out <<"<";
		case Comparison::leq:
			return out <<"<=";
		case Comparison::gt:
			return out <<">";
		default:
			return out <<">=";
	}
}

template<typename Weight>
class BVTheorySolver: public Theory {
public:
	enum class OperationType{
		none = 0,
		cause_is_bits = 1,
		refined_cause = 2,
		cause_is_comparison = 3,
		cause_is_bv_comparison = 4,
		cause_is_addition = 5,
		cause_is_addition_argument = 6,
		cause_is_condition = 7,
		cause_is_condition_argument = 8,
		cause_is_decision=9,
		cause_is_minmax = 10,
		cause_is_minmax_argument = 11,
		cause_is_popcount = 12,
		cause_is_theory=13,
		cause_is_invert = 14
	};
	struct Cause{

		OperationType type:5;//don't need this anymore.
		int index:27;
		Cause(const Cause & copy):type(copy.type),index(copy.index){

			}

		Cause():type(OperationType::none),index(-1){

		}
		bool hasCause(){
			//this can be improved, if we want.
			return type != OperationType::none;
		}

		void clear(){
			type= OperationType::none;
			index=-1;
		}
	};
	//An assignment is either an assignment to a variable, or is a tightening of a bitvector's bounds
	struct Assignment {
		bool isOperation:1;
		bool assign :1;
		int bvID:30;
		Var var;
		//ensure that only one of under or over can change per assignment, and only store under or over
		//replace var, bvID with operationID.
		Weight previous_under;
		Weight previous_over;
		Weight new_under;
		Weight new_over;

		Cause prev_under_cause;
		Cause prev_over_cause;
		Cause new_under_cause;
		Cause new_over_cause;
		Assignment( int bvID,Weight previous_under,Weight previous_over,Weight new_under,Weight new_over,
				Cause prev_under_cause,Cause prev_over_cause,Cause new_under_cause,Cause new_over_cause):isOperation(false),assign(0),bvID(bvID),var(var_Undef),previous_under(previous_under),previous_over(previous_over),new_under(new_under),new_over(new_over),
				prev_under_cause(prev_under_cause),prev_over_cause(prev_over_cause),new_under_cause(new_under_cause),new_over_cause(new_over_cause){

		}
		Assignment(bool isOperation, bool assign, int bvID, Var v):isOperation(isOperation),assign(assign),bvID(bvID),var(v),previous_under(-1),previous_over(-1),new_under(-1),new_over(-1),
				prev_under_cause(),prev_over_cause(),new_under_cause(),new_over_cause(){
		}

		bool isBoundAssignment()const{
			if(new_under>-1){
				assert(!isOperation);
			}
			return new_under>-1;
		}
	};

	class Operation{
			friend BVTheorySolver;
			int op_id=0;

		protected:
			BVTheorySolver & theory;

			Operation(BVTheorySolver & theory):theory(theory)
			{
				op_id = theory.operations.size();
				theory.operations.push(this);
				theory.comparison_needs_repropagation.growTo(theory.operations.size(),false);
			}
			virtual ~Operation(){

			}

		public:
			int getID()const{
				return op_id;
			}
			virtual int getBV()=0;
			virtual OperationType getType()const=0;

			virtual bool propagate(bool & changed, vec<Lit> & conflict)=0;
			virtual void enqueue(Lit l, bool alter_trail=true){
				throw std::runtime_error("No implementation");
			}
			virtual void backtrack(Assignment & e, bool rewind){

			}

			virtual void updateApprox(Var ignore_bv, Weight & under_new, Weight & over_new, Cause & under_cause_new, Cause & over_cause_new)=0;
			virtual bool checkApproxUpToDate(Weight & under,Weight&over){
				return true;
			}
			//Only required if the operation can assign literals
			virtual void buildReason(Lit p,CRef marker,vec<Lit> & conflict){
				throw std::runtime_error("No implementation");
			}
			virtual void analyzeReason(bool compareOver,Comparison op, Weight  to,  vec<Lit> & conflict)=0;
			virtual void move(int bvID){

			}
			virtual bool checkSolved()=0;
	};
//This is a convenience macro, to import a whole bunch of BVTheorySolver functions and methods into the local namespace so that they can be called without 'theory.'
//I'm relying on any unused lambda definitions here being completely optimized out by the compiler
#define importTheory(theory) vec<Weight> & under_approx =theory.under_approx;\
	vec<Weight> &over_approx=theory.over_approx;\
	auto value = [&](Lit l) {return theory.value(l);};\
	auto clip_over = [&](Weight & val, int bvID) {return theory.clip_over(val,bvID);};\
	auto clip_under = [&](Weight & val, int bvID) {return theory.clip_under(val,bvID);};\
	auto clip = [&](Weight & val, int bvID) {return theory.clip(val,bvID);};\
	auto assert_in_range = [&](Weight  val, int bvID){theory.assert_in_range(val, bvID);};\
	auto addAlteredBV = [&]( int newBV){theory.addAlteredBV(newBV);};\
	auto toSolver = [&](Lit l){return theory.toSolver(l);};\
	auto analyze = [&](vec<Lit> & c){theory.analyze(c);};\
	auto analyzeValueReason = [&](Comparison op, int bvID, Weight  to,  vec<Lit> & conflict){theory.analyzeValueReason(op,bvID,to,conflict);};\
	auto addAnalysis = [&](Comparison op, int bvID, Weight  to){return theory.addAnalysis(op,bvID,to);};\
	auto getUnderApprox = [&](int bvID, bool level0){return theory.getUnderApprox(bvID,level0);};\
	auto getOverApprox = [&](int bvID, bool level0){return theory.getOverApprox(bvID,level0);};\
	auto  applyOp = [&](Comparison op,int bvID,Weight to){return theory.applyOp(op,bvID,to);};\
	auto enqueue = [&](Lit l,CRef reason){return theory.enqueue(l,reason);};

	/*auto toSolver = [&](vec<Lit> & c){theory.toSolver(c);};\ */


	class BitOp:public Operation{
	public:
		using Operation::getID;
		using Operation::theory;


	public:
		int bvID=-1;
		int getBV()override{
			return bvID;
		}

		BitOp(BVTheorySolver & theory,int bvID):Operation(theory),bvID(bvID){

		}
		 void move(int bvID)override{
			 this->bvID=bvID;
		 }
		OperationType getType()const override{
			return OperationType::cause_is_bits;
		}
		void enqueue(Lit l, bool alter_trail){
			importTheory(theory);
			Var v = var(l);
			if(alter_trail){

				addAlteredBV(bvID);
			}
		}
		bool propagate(bool & changed,vec<Lit> & conflict) override{
			vec<Lit> & bv = theory.bitvectors[bvID];
			importTheory(theory);


			Weight & underApprox = under_approx[bvID];
			Weight & overApprox = over_approx[bvID];
			Weight under;
			Weight over;
			bool new_change;
			do{
				if(underApprox>overApprox){
					double startconftime = rtime(2);
					theory.stats_num_conflicts++;
					if(opt_verb>1){
						printf("bv approximation update conflict %ld\n", theory.stats_num_conflicts);
					}
					//propagationtime += startconftime - startproptime;
					//this is already a conflict;
					analyzeValueReason(Comparison::geq, bvID,underApprox,conflict);
					analyzeValueReason(Comparison::leq, bvID,overApprox,conflict);

					analyze(conflict);

					theory.stats_conflict_time+=rtime(2)-startconftime;
					return false;
				}

				under =0;
				over=(1L<<bv.size())-1;
				new_change=false;
				//for(int i = 0;i<bv.size();i++){
				for(int i = bv.size()-1;i>=0;i--){
					Weight bit = 1<<i;
					lbool val = value(bv[i]);
					Lit l = bv[i];
					if(val==l_True){
						under+=bit;
						assert_in_range(under,bvID);
						if(under>overApprox){
								double startconftime = rtime(2);
								//propagationtime += startconftime - startproptime;
								//this is a conflict
								for(int j = bv.size()-1;j>=i;j--){
									Weight bit = 1<<j;
									lbool val = value(bv[j]);
									if(val==l_True){
										conflict.push(toSolver(~bv[j]));
									}
								}
								theory.stats_num_conflicts++;theory.stats_bit_conflicts++;
								if(opt_verb>1){
									printf("bv bit conflict %ld\n", theory.stats_num_conflicts);
								}
								theory.buildComparisonReason(Comparison::leq,bvID,overApprox,conflict);

								theory.stats_conflict_time+=rtime(2)-startconftime;
								return false;

						}
					}else if (val==l_False){
						over-=bit;
						assert_in_range(under,bvID);
						if(over<underApprox){
								double startconftime = rtime(2);
								//propagationtime += startconftime - startproptime;
								//this is a conflict. Either this bit, or any previously assigned false bit, must be true, or the underapprox must be larger than it currently is.
								//is this really the best way to handle this conflict?
								for(int j = bv.size()-1;j>=i;j--){
									Weight bit = 1<<j;
									lbool val = value(bv[j]);
									if(val==l_False){
										conflict.push(toSolver(bv[j]));
									}
								}
								theory.stats_num_conflicts++;theory.stats_bit_conflicts++;
								if(opt_verb>1){
									printf("bv bit conflict %ld\n", theory.stats_num_conflicts);
								}
								theory.buildComparisonReason(Comparison::geq,bvID,underApprox,conflict);

								theory.stats_conflict_time+=rtime(2)-startconftime;
								return false;

						}
					}else{
						assert_in_range(over-bit,bvID);
						assert_in_range(under+bit,bvID);
						if(over-bit < underApprox){
							enqueue(l, theory.bvprop_marker);
							under+=bit;
						}else if (under+bit>overApprox){
							enqueue(~l, theory.bvprop_marker);
							over-=bit;
						}
					}
				}
				new_change = theory.updateApproximations(bvID);
				changed|=new_change;
				/*if(new_change){
					printf("iter %d, bv %d, under ",realprops , bvID); //: %d, over %d\n", bvID, underApprox,overApprox);
						std::cout<<underApprox << " over ";
						std::cout<<overApprox << "\n";
						fflush(stdout);
				}*/
			}while(new_change);//the bit assignment updates above can force a more precise over or under approximation, which can in turn lead to further bit assignments (I think this can happen?).
			return true;
		}


		bool checkApproxUpToDate(Weight & under_out,Weight&over_out)override{
			importTheory(theory);
			Weight under=0;
			Weight over=0;
			vec<Lit> & bv = theory.bitvectors[bvID];
			for(int i = 0;i<bv.size();i++){
				lbool val = value(bv[i]);
				if(val==l_True){
					Weight bit = 1<<i;
					under+=bit;
					over+=bit;
				}else if (val==l_False){
				}else{
					Weight bit = 1<<i;
					over+=bit;
				}
			}
			if(under>under_out)
				under_out=under;
			if(over<over_out)
				over_out=over;
			return true;
		}

		void updateApprox(Var ignore_bv, Weight & under_new, Weight & over_new, Cause & under_cause_new, Cause & over_cause_new) override{
			importTheory(theory);
			Weight under=0;
			Weight over=0;
			vec<Lit> & bv = theory.bitvectors[bvID];
			for(int i = 0;i<bv.size();i++){
				if(var(bv[i])==ignore_bv){
					Weight bit = 1<<i;
					over+=bit;
					continue;
				}
				lbool val = value(bv[i]);
				if(val==l_True){
					Weight bit = 1<<i;
					under+=bit;
					over+=bit;

				}else if (val==l_False){

				}else{
					Weight bit = 1<<i;
					over+=bit;
				}
			}
			//needed in case a decision was made, to preserve that decision's cause here...
			if(under >under_new){
				under_new=under;
				under_cause_new.clear();
				under_cause_new.type =getType();
				under_cause_new.index=getID();
			}
			if(over<over_new){
				over_new=over;
				over_cause_new.clear();
				over_cause_new.type =getType();
				over_cause_new.index=getID();
			}
		}

		void buildReason(Lit p,CRef marker, vec<Lit> & reason)override{
			importTheory(theory);
			assert (marker == theory.bvprop_marker);
			Var v = var(p);
			reason.push(toSolver(p));
			theory.rewind_trail_pos(theory.analysis_trail_pos-1);
			assert(value(p)==l_Undef);
			Weight underApprox = under_approx[bvID];
			Weight overApprox = over_approx[bvID];


			assert(under_approx>=0); assert(overApprox>=0);
			vec<Lit> & bv = theory.bitvectors[bvID];

			int bitpos=-1;
			for(int i = bv.size()-1;i>=0;i--){
				lbool val = value(bv[i]);
				Lit l = bv[i];
				if(var(bv[i])==v){
					bitpos=i;
					break;
				}
			}
			//Weight bit = 1<<bitpos;

			assert(bitpos>=0);
			Weight under = 0;
			Weight over=(1L<<bv.size())-1;

			for(int i = bv.size()-1;i>=0;i--){
				Weight bit = 1<<i;
				lbool val = value(bv[i]);
				Lit l = bv[i];
				if(val==l_True){
					under+=bit;
					//assert(under<=overApprox);
				}else if (val==l_False){
					over-=bit;
					//assert(over>=underApprox);
				}else if (bitpos==i){

					if(over-bit<underApprox){
						assert(bitpos==i);
						//this is a conflict. Either this bit, or any previously assigned false bit, must be true, or the underapprox must be larger than it currently is.
						//is this really the best way to handle this conflict?
						for(int j = bv.size()-1;j>=i;j--){
							Weight bit = 1<<j;
							lbool val = value(bv[j]);
							if(val==l_False){
								reason.push(toSolver(bv[j]));
							}
						}
						theory.buildComparisonReason(Comparison::geq,bvID,underApprox,reason);
						break;

					}
					if(under+bit>overApprox){
						assert(bitpos==i);
						//this is a conflict
						for(int j = bv.size()-1;j>=i;j--){
							Weight bit = 1<<j;
							lbool val = value(bv[j]);
							if(val==l_True){
								reason.push(toSolver(~bv[j]));
							}
						}
						theory.buildComparisonReason(Comparison::leq,bvID,overApprox,reason);
						break;
					}
				}
			}

		}


		void analyzeReason(bool compareOver,Comparison op, Weight  to,  vec<Lit> & conflict) override{
			importTheory(theory);
			vec<Lit> & bv = theory.bitvectors[bvID];
			if(compareOver){
				Weight over = over_approx[bvID];
				for(int i =0;i<bv.size();i++){
					Lit bl = bv[i];
					if(value(bl)==l_False ){
						Weight bit = 1<<i;
						if(theory.comp(op,over+bit,to)&& theory.level(var(bl))>0){
							//then we can skip this bit, because we would still have had a conflict even if it was assigned true.
							over+=bit;
						}else{
							assert(value(bl)==l_False);
							conflict.push(toSolver(bl));
						}
					}
				}
			}else{
				Weight under = under_approx[bvID];
				for(int i =0;i<bv.size();i++){
					Lit bl = bv[i];
					lbool val = value(bl);

					if(value(bl)==l_True){
						Weight bit = 1<<i;
						if(theory.comp(op,under-bit,to)  && theory.level(var(bl))>0){
							//then we can skip this bit, because we would still have had a conflict even if it was assigned false.
							under-=bit;
						}else{
							assert(value(~bl)==l_False);
							conflict.push(toSolver(~bl));
						}
					}
				}
			}

		}

		bool checkSolved()override{
			importTheory(theory);
			vec<Lit> & bv = theory.bitvectors[bvID];
			Weight over=0;
			Weight under=0;
			assert(bv.size());

			for(int i = 0;i<bv.size();i++){
				lbool val = value(bv[i]);
				if(val==l_True){
					Weight bit = 1<<i;
					under+=bit;
					over+=bit;
				}else if (val==l_False){

				}else{
					Weight bit = 1<<i;
					over+=bit;
				}
			}
			if(over!=over_approx[bvID])
				return false;
			if(under!=under_approx[bvID])
				return false;

			return true;
		}
	};
	class ComparisonOp:public Operation{
		public:
			using Operation::getID;
			using Operation::theory;




		public:
			int bvID:30;
			const int is_lt:1;
			const int is_strict:1;
			const Lit l;
			const Weight w;

			ComparisonOp(BVTheorySolver & theory,int bvID, Comparison op,Weight w,Lit l):Operation(theory),bvID(bvID),is_lt(op==Comparison::lt || op==Comparison::leq),is_strict(op==Comparison::lt || op==Comparison::gt),l(l),w(w){

			}
			int getBV()override{
				return bvID;
			}

			Comparison getOp(){
				if(is_lt && is_strict){
					return Comparison::lt;
				}else if(is_lt){
					return Comparison::leq;
				}else if(is_strict){
					return Comparison::gt;
				}else{
					return Comparison::geq;
				}
			}
			 void move(int bvID)override{
				 this->bvID=bvID;
			 }
			OperationType getType()const override{
				return OperationType::cause_is_comparison;
			}
			void enqueue(Lit l, bool alter_trail)override{
				importTheory(theory);
				Var v= var(l);
				if(alter_trail){
					//theory.trail.push( { true, !sign(l),getID(), v });
					theory.bv_needs_propagation[bvID]=true;


					addAlteredBV(bvID);
				}
			}
			void backtrack(Assignment & e, bool rewind)override{
				if(!rewind){
					if (theory.comparison_needs_repropagation[getID()]){
						if(!theory.alteredBV[bvID]){
							theory.alteredBV[bvID]=true;
							theory.altered_bvs.push(bvID);
						}
						theory.bv_needs_propagation[bvID]=true;
						theory.requiresPropagation=true;
						theory.S->needsPropagation(theory.getTheoryIndex());
					}
				}
			}
			bool propagate(bool & changed,vec<Lit> & conflict) override{
				return propagate(changed,conflict,true) && propagate(changed,conflict,false);
			}
			bool propagate(bool & changed,vec<Lit> & conflict, bool propagate_over_approx){
				importTheory(theory);
				Weight & underApprox = under_approx[bvID];
				Weight & overApprox = over_approx[bvID];
				Comparison op = getOp();
				const Weight & to = w;
				if(propagate_over_approx){
					if((op==Comparison::lt && overApprox<to) ||
							(op==Comparison::leq && overApprox<=to)){
						if(value(l)==l_True){
							//do nothing

						}else if (value(l)==l_False){
							double startconftime = rtime(2);

							theory.stats_num_conflicts++;theory.stats_compare_conflicts++;
							if(opt_verb>1){
								printf("bv comparison conflict %ld\n", theory.stats_num_conflicts);
							}
							assert(value(l)==l_False);
							assert(theory.dbg_value(l)==l_False);
							conflict.push(toSolver(l));
							theory.buildComparisonReason(op,bvID,to,conflict);

							theory.stats_conflict_time+=rtime(2)-startconftime;
							return false;
						}else {
							assert(value(l)==l_Undef);
							enqueue(l,theory.comparisonprop_marker);
						}
					}else if((op==Comparison::gt && overApprox<=to) ||
							(op==Comparison::geq && overApprox<to)){
						if(value(l)==l_True){
							double startconftime = rtime(2);

							theory.stats_num_conflicts++;theory.stats_compare_conflicts++;
							if(opt_verb>1){
								printf("bv comparison conflict %ld\n", theory.stats_num_conflicts);
							}
							conflict.push(toSolver(~l));
							theory.buildComparisonReason(-op,bvID,to,conflict);

							theory.stats_conflict_time+=rtime(2)-startconftime;
							return false;
						}else if (value(l)==l_False){

						}else {
							assert(value(l)==l_Undef);
							enqueue(~l, theory.comparisonprop_marker);
						}
					}
				}else{
					if((op==Comparison::lt && underApprox>=to) ||
							(op==Comparison::leq && underApprox>to)){
						if(value(l)==l_True){
							double startconftime = rtime(2);
							//propagationtime += startconftime - startproptime;
							theory.stats_num_conflicts++;theory.stats_compare_conflicts++;
							if(opt_verb>1){
								printf("bv comparison conflict %ld\n", theory.stats_num_conflicts);
							}
							conflict.push(toSolver(~l));
							theory.buildComparisonReason(-op,bvID,to,conflict);

							theory.stats_conflict_time+=rtime(2)-startconftime;
							return false;
						}else if (value(l)==l_False){
							//do nothing

						}else {
							assert(value(l)==l_Undef);
							enqueue(~l, theory.comparisonprop_marker);
						}
					}else if((op==Comparison::gt && underApprox>to) ||
							(op==Comparison::geq && underApprox>=to)){
						if(value(l)==l_True){

						}else if (value(l)==l_False){
							double startconftime = rtime(2);
							//propagationtime += startconftime - startproptime;
							theory.stats_num_conflicts++;theory.stats_compare_conflicts++;
							if(opt_verb>1){
								printf("bv comparison conflict %ld\n", theory.stats_num_conflicts);
							}
							conflict.push(toSolver(l));
							theory.buildComparisonReason(op,bvID,to,conflict);

							theory.stats_conflict_time+=rtime(2)-startconftime;
							return false;
						}else {
							assert(value(l)==l_Undef);
							enqueue(l,theory.comparisonprop_marker);
						}
					}
				}
				return true;
			}
			void updateApprox(Var ignore_bv, Weight & under_new, Weight & over_new, Cause & under_cause_new, Cause & over_cause_new) override{
				updateApprox(ignore_bv,under_new,over_new,under_cause_new,over_cause_new,true);
				updateApprox(ignore_bv,under_new,over_new,under_cause_new,over_cause_new,false);
			}

			void updateApprox(Var ignore_bv, Weight & under_new, Weight & over_new, Cause & under_cause_new, Cause & over_cause_new, bool update_over_approx){
				importTheory(theory);
				if(update_over_approx){
					Comparison op = getOp();
					bool setOver=false;
					switch(op){
						case Comparison::lt:
							if(value(l)==l_True && over_new>=w){
								if(w-1>=theory.getUnderApprox(bvID,true)){
									over_new=w-1;
									setOver=true;
								 }
							}
							break;
						case Comparison::leq:
							if(value(l)==l_True && over_new>w){
								if(w>=theory.getUnderApprox(bvID,true)){
									over_new=w;
									setOver=true;
								 }
							}
							break;
						case Comparison::gt:
							if (value(l)==l_False && over_new>w){
								if(w>=theory.getUnderApprox(bvID,true)){
									over_new=w;
									setOver=true;
								 }
							}
							break;
						case Comparison::geq:
						default:
							if (value(l)==l_False && over_new>=w){
								if(w-1>=theory.getUnderApprox(bvID,true)){
									over_new=w-1;
									setOver=true;
								 }
							}
							break;
					}

					if(setOver){
						over_cause_new.clear();
						over_cause_new.type =getType();
						over_cause_new.index=getID();
					}
				}else{
					Comparison op = getOp();
					bool setUnder=false;
					switch(op){
						case Comparison::lt:
							 if (value(l)==l_False && under_new<w){
								 if(w<=theory.getOverApprox(bvID,true)){
									 under_new=w;
									 setUnder=true;
								 }
							}
							break;
						case Comparison::leq:
							 if (value(l)==l_False && under_new<=w){
								 if(w<=theory.getOverApprox(bvID,true)){
									 under_new=w+1;
									setUnder=true;
								 }
							}
							break;
						case Comparison::gt:
							if(value(l)==l_True && under_new<=w){
								 if(w+1<=theory.getOverApprox(bvID,true)){
									 under_new=w+1;
									 setUnder=true;
								 }
							}
							break;
						case Comparison::geq:
						default:
							if(value(l)==l_True && under_new<w){
								 if(w<=theory.getOverApprox(bvID,true)){
									 under_new=w;
									 setUnder=true;
								 }
							}
							break;
					}

					if(setUnder){
						under_cause_new.clear();
						under_cause_new.type =getType();
						under_cause_new.index=getID();
					}
				}
			}

			bool checkApproxUpToDate(Weight & under,Weight&over)override{
				importTheory(theory);
				Comparison op = getOp();

				lbool val = value(l);

				switch(op){
					case Comparison::lt:
						if(value( l)==l_True && over>=w){
							over=w-1;
						}else if (value(l)==l_False && under<w){
							under=w;
						}
						break;
					case Comparison::leq:
						if(value( l)==l_True && over>w){
							over=w;
						}else if (value(l)==l_False && under<=w){
							under=w+1;
						}
						break;
					case Comparison::gt:
						if(value( l)==l_True && under<=w){
							under=w+1;
						}else if (value(l)==l_False && over>w){
							over=w;
						}
						break;
					case Comparison::geq:
					default:
						if(value( l)==l_True && under<w){
							under=w;
						}else if (value(l)==l_False && over>=w){
							over=w-1;
						}
						break;
				}
				return true;
			}

			void buildReason(Lit p,CRef marker, vec<Lit> & reason)override{
				importTheory(theory);
				reason.push(toSolver(p));
				Var v = var(p);
				assert(marker==theory.comparisonprop_marker);
				Comparison op = getOp();
				if(sign(p)){
					op=-op;
				}
				theory.buildComparisonReason(op,bvID,w,reason);

			}


			void analyzeReason(bool compareOver,Comparison op, Weight  to,  vec<Lit> & conflict) override{
				importTheory(theory);
				Comparison cop = getOp();//invert this because we are switch the direction of comparison
				assert(l!=lit_Undef);
				if(value(l)==l_True){
					assert(value(~l)==l_False);
					conflict.push(toSolver(~l));
				}else{
					assert(value(l)==l_False);
					conflict.push(toSolver(l));
				}
			}

			bool checkSolved()override{
				importTheory(theory);
				Comparison op = getOp();
				switch (op){
					case Comparison::lt:
						if(value(l)==l_True && under_approx[bvID]>= w){
							return false;
						}else if (value(l)==l_False && over_approx[bvID]<w){
							return false;
						}
						break;
					case Comparison::leq:
						if(value(l)==l_True && under_approx[bvID]> w){
							return false;
						}else if (value(l)==l_False && over_approx[bvID]<=w){
							return false;
						}
						break;
					case Comparison::gt:
						if(value(l)==l_True && over_approx[bvID]<= w){
							return false;
						}else if (value(l)==l_False && under_approx[bvID]>w){
							return false;
						}
						break;
					case Comparison::geq:
					default:
						if(value(l)==l_True && over_approx[bvID] < w){
							return false;
						}else if (value(l)==l_False && under_approx[bvID]>= w){
							return false;
						}
						break;
				}

				return true;
			}
		};
	class ComparisonBVOp:public Operation{
			public:
		using Operation::getID;
		using Operation::theory;




	public:
		int bvID:30;
		const int is_lt:1;
		const int is_strict:1;
		const Lit l;

		ComparisonBVOp * other;

		ComparisonBVOp(BVTheorySolver & theory,int bvID, Comparison op,Lit l):Operation(theory),bvID(bvID),is_lt(op==Comparison::lt || op==Comparison::leq),is_strict(op==Comparison::lt || op==Comparison::gt),l(l){

		}
		int getBV()override{
			return bvID;
		}
		void setOtherOp(ComparisonBVOp * other){
			this->other=other;
		}

		int getCompareID()const{
			return other->bvID;
		}

		Comparison getOp(){
			if(is_lt && is_strict){
				return Comparison::lt;
			}else if(is_lt){
				return Comparison::leq;
			}else if(is_strict){
				return Comparison::gt;
			}else{
				return Comparison::geq;
			}
		}
		 void move(int bvID)override{
			 this->bvID=bvID;

		 }
		OperationType getType()const override{
			return OperationType::cause_is_bv_comparison;
		}
		void enqueue(Lit l, bool alter_trail)override{
			importTheory(theory);
			if(alter_trail){
				Var v= var(l);

				//theory.trail.push( { true, !sign(l),getID(), v });
				theory.bv_needs_propagation[bvID]=true;


				addAlteredBV(bvID);
				//addAlteredBV(getCompareID());
			}
		}
		void backtrack(Assignment & e, bool rewind)override{
			if(!rewind){
				if (theory.comparison_needs_repropagation[getID()]){
					if(!theory.alteredBV[bvID]){
						theory.alteredBV[bvID]=true;
						theory.altered_bvs.push(bvID);
					}
					theory.bv_needs_propagation[bvID]=true;
					theory.requiresPropagation=true;
					theory.S->needsPropagation(theory.getTheoryIndex());
				}
			}
		}
		bool propagate(bool & changed,vec<Lit> & conflict) override{
			return propagate(changed,conflict,true) && propagate(changed,conflict,false);
		}
		bool propagate(bool & changed,vec<Lit> & conflict, bool propagate_over_approx){
			importTheory(theory);
			Weight & underApprox = under_approx[bvID];
			Weight & overApprox = over_approx[bvID];
			Comparison op = getOp();
			int compareID = getCompareID();
			assert(compareID>=0);

			CRef marker;
			if(compareID>bvID){
				marker = theory.comparisonbv_2ndarg_marker;
				assert(theory.getOperationID(var(l))==getID()-1);
			}else{
				marker = theory.comparisonbv_1starg_marker;
				assert(theory.getOperationID(var(l))==getID());
			}


			lbool val = value(l);
			Lit l = this->l;
			if(value(l)==l_False){
				l=~l;
				op=-op;
			}
			assert(value(l)!=l_False);
			if(propagate_over_approx){
				Weight & under_compare = under_approx[compareID];

				if((op==Comparison::lt && overApprox<under_compare) ||
						(op==Comparison::leq && overApprox<=under_compare)){
					if(value(l)==l_True){
						//do nothing

					}/*else if (value(l)==l_False){

						assert(value(l)==l_False);
						assert(dbg_value(l)==l_False);
						conflict.push(toSolver(l));
						buildValueReasonBV(op,bvID,compareID,conflict);

						return false;
					}*/else {

						assert(value(l)==l_Undef);
						enqueue(l,marker);
					}
				}else if((op==Comparison::gt && overApprox<=under_compare) ||
						(op==Comparison::geq && overApprox<under_compare)){
					if(value(l)==l_True){
						double startconftime = rtime(2);
						//propagationtime += startconftime - startproptime;
						theory.stats_num_conflicts++;theory.stats_bv_compare_conflicts++;
						if(opt_verb>1){
							printf("bv bv comparison conflict %ld\n", theory.stats_num_conflicts);
						}
						conflict.push(toSolver(~l));
						theory.buildComparisonReasonBV(-op,bvID,compareID,conflict);

						theory.stats_conflict_time+=rtime(2)-startconftime;
						return false;
					}/*else if (value(l)==l_False){

					}*/else {
						assert(value(l)==l_Undef);
						enqueue(~l, marker);
					}
				}
				if(value(l)==l_True &&((op==Comparison::lt && underApprox>=under_compare) ||
						(op==Comparison::leq && underApprox>under_compare))){
					addAlteredBV(compareID);
				}
			}else{
				Weight & over_compare = over_approx[compareID];
				if((op==Comparison::lt && underApprox>=over_compare) ||
						(op==Comparison::leq && underApprox>over_compare)){
					if(value(l)==l_True){
						double startconftime = rtime(2);
						//propagationtime += startconftime - startproptime;
						theory.stats_num_conflicts++;theory.stats_bv_compare_conflicts++;
						if(opt_verb>1){
							printf("bv bv comparison conflict %ld\n", theory.stats_num_conflicts);
						}
						conflict.push(toSolver(~l));
						theory.buildComparisonReasonBV(-op,bvID,compareID,conflict);

						theory.stats_conflict_time+=rtime(2)-startconftime;
						return false;
					}/*else if (value(l)==l_False){
						//do nothing

					}*/else {
						assert(value(l)==l_Undef);
						enqueue(~l, marker);
					}
				}else if((op==Comparison::gt && underApprox>over_compare) ||
						(op==Comparison::geq && underApprox>=over_compare)){
					if(value(l)==l_True){

					}/*else if (value(l)==l_False){

						conflict.push(toSolver(l));
						buildValueReasonBV(op,bvID,compareID,conflict);

						return false;
					}*/else {

						assert(value(l)==l_Undef);
						enqueue(l, marker);
					}
				}
				//we can also update the other bv's approximation, possibly:
				if(value(l)==l_True &&( (op==Comparison::gt && overApprox<=over_compare) ||
						(op==Comparison::geq && overApprox<over_compare))){
					//the other bv needs to be updated
					addAlteredBV(compareID);
				}
			}
			return true;
		}
		void updateApprox(Var ignore_bv, Weight & under_new, Weight & over_new, Cause & under_cause_new, Cause & over_cause_new) override{
			updateApprox(ignore_bv,under_new,over_new,under_cause_new,over_cause_new,true);
			updateApprox(ignore_bv,under_new,over_new,under_cause_new,over_cause_new,false);
		}

		void updateApprox(Var ignore_bv, Weight & under_new, Weight & over_new, Cause & under_cause_new, Cause & over_cause_new, bool update_over_approx) {
			importTheory(theory);
			Comparison op = getOp();

			int compareID= getCompareID();

			if(update_over_approx){
				Weight w = over_approx[compareID];
				bool setOver=false;
				switch(op){
					case Comparison::lt:
						if(value(l)==l_True && over_new>=w){
							if(w-1>=getUnderApprox(bvID,true)){
								over_new=w-1;
								setOver=true;
							}
						}
						break;
					case Comparison::leq:
						if(value(l)==l_True && over_new>w){
							if(w>=getUnderApprox(bvID,true)){
								over_new=w;
								setOver=true;
							}
						}
						break;
					case Comparison::gt:
						if (value(l)==l_False && over_new>w){
							if(w>=getUnderApprox(bvID,true)){
								over_new=w;
								setOver=true;
							}
						}
						break;
					case Comparison::geq:
					default:
						if (value(l)==l_False && over_new>=w){
							if(w-1>=getUnderApprox(bvID,true)){
								over_new=w-1;
								setOver=true;
							}
						}
						break;
				}

				if(setOver){
					over_cause_new.clear();
					over_cause_new.type =getType();
					over_cause_new.index=getID();
				}
			}else{
				Weight w = under_approx[compareID];
				bool setUnder=false;
				switch(op){
					case Comparison::lt:
						 if (value(l)==l_False && under_new<w){
							 if(w<=getOverApprox(bvID,true)){
								 under_new=w;
								 setUnder=true;
							 }
						}
						break;
					case Comparison::leq:
						 if (value(l)==l_False && under_new<=w){
							 if(w+1<=getOverApprox(bvID,true)){
								under_new=w+1;
								setUnder=true;
							 }
						}
						break;
					case Comparison::gt:
						if(value(l)==l_True && under_new<=w){
							 if(w+1<=getOverApprox(bvID,true)){
								 under_new=w+1;
								 setUnder=true;
							 }
						}
						break;
					case Comparison::geq:
					default:
						if(value(l)==l_True && under_new<w){
							 if(w<=getOverApprox(bvID,true)){
								 under_new=w;
								 setUnder=true;
							 }
						}
						break;
				}

				if(setUnder){
					under_cause_new.clear();
					under_cause_new.type = getType();
					under_cause_new.index=getID();
				}
			}
		}

		bool checkApproxUpToDate(Weight & under,Weight&over)override{
			importTheory(theory);
			Comparison op = getOp();

			lbool val = value(l);
			int compareID = getCompareID();
			Weight under_w = under_approx[compareID];
			Weight over_w = over_approx[compareID];

			switch(op){

				case Comparison::lt:
					if(value( l)==l_True && over>=over_w){
						over=over_w-1;
					}else if (value(l)==l_False && under<under_w){
						under=under_w;
					}
					break;
				case Comparison::leq:
					if(value( l)==l_True && over>over_w){
						over=over_w;
					}else if (value(l)==l_False && under<=under_w){
						under=under_w+1;
					}
					break;
				case Comparison::gt:
					if(value( l)==l_True && under<=under_w){
						under=under_w+1;
					}else if (value(l)==l_False && over>over_w){
						over=over_w;
					}
					break;
				case Comparison::geq:
				default:
					if(value( l)==l_True && under<under_w){
						under=under_w;
					}else if (value(l)==l_False && over>=over_w){
						over=over_w-1;
					}
					break;
			}
			return true;
		}

		void buildReason(Lit p,CRef marker, vec<Lit> & reason)override{
			importTheory(theory);
			 if (marker == theory.comparisonbv_1starg_marker) {
				reason.push(toSolver(p));
				Var v = var(p);
				int compareBV = getCompareID();
				Comparison op = getOp();
				if (compareBV<0){
					throw std::runtime_error("Critical error in BV solver");
				}

				if(sign(p)){
					op=-op;
				}
				theory.buildComparisonReasonBV(op,bvID,compareBV,reason);
			}else if (marker == theory.comparisonbv_2ndarg_marker) {
				other->buildReason(p,theory.comparisonbv_1starg_marker,reason);
			}else if (marker == theory.comparisonbv_prop_marker){
				reason.push(toSolver(p));
				Var v = var(p);

				Comparison op = getOp();
				if(sign(p)){
					op=-op;
				}
				int compareBV = getCompareID();
				theory.buildComparisonReasonBV(op,bvID,compareBV,reason);
			}

		}


		void analyzeReason(bool compare_over,Comparison op, Weight  to,  vec<Lit> & conflict) override{
			importTheory(theory);
			Comparison cop = getOp();//invert this because we are switch the direction of comparison
			assert(l!=lit_Undef);
			int compareID = getCompareID();
			if(value(l)==l_True){
				assert(value(~l)==l_False);
				conflict.push(toSolver(~l));

				if (compare_over){
					addAnalysis(Comparison::leq,compareID, over_approx[compareID]);
					//buildValueReason(Comparison::leq,c.compareID, over_approx[c.compareID],conflict,trail_pos-1);
				}else{
					addAnalysis(Comparison::geq,compareID, under_approx[compareID]);
					//buildValueReason(Comparison::geq,c.compareID, under_approx[c.compareID],conflict,trail_pos-1);
				}

			}else{
				assert(value(l)==l_False);
				conflict.push(toSolver(l));
				if (compare_over){
					addAnalysis(Comparison::leq,compareID, over_approx[compareID]);
				}else{
					addAnalysis(Comparison::geq,compareID, under_approx[compareID]);
				}
			}

		}

		bool checkSolved()override{
			importTheory(theory);
			Comparison op = getOp();
			int toID =getCompareID();
			vec<Lit> & bv_compare = theory.bitvectors[toID];
			Weight over_compare=0;
			Weight under_compare=0;
			for(int i = 0;i<bv_compare.size();i++){
				lbool val = value(bv_compare[i]);
				if(val==l_True){
					Weight bit = 1<<i;
					under_compare+=bit;
					over_compare+=bit;
				}else if (val==l_False){

				}else{
					Weight bit = 1<<i;
					over_compare+=bit;
				}
			}

			switch (op){
				case Comparison::lt:
					if(value(l)==l_True && under_approx[bvID]>= over_compare){
						return false;
					}else if (value(l)==l_False && over_approx[bvID]<under_compare){
						return false;
					}
					break;
				case Comparison::leq:
					if(value(l)==l_True && under_approx[bvID]> over_compare){
						return false;
					}else if (value(l)==l_False && over_approx[bvID]<=under_compare){
						return false;
					}
					break;
				case Comparison::gt:
					if(value(l)==l_True && over_approx[bvID]<= under_compare){
						return false;
					}else if (value(l)==l_False && under_approx[bvID]>over_compare){
						return false;
					}
					break;
				case Comparison::geq:
				default:
					if(value(l)==l_True && over_approx[bvID] < under_compare){
						return false;
					}else if (value(l)==l_False && under_approx[bvID]>= over_compare){
						return false;
					}
					break;
			}

			return true;
		}
	};
	class PopCountOp:public Operation{
			public:
		using Operation::getID;
		using Operation::theory;

	public:
		Weight under=0;
		Weight over=0;
		int bvID=-1;

		vec<Var>  args;
		PopCountOp(BVTheorySolver & theory,int bvID):Operation(theory),bvID(bvID){

		}
		int getBV()override{
			return bvID;
		}
		OperationType getType()const override{
			return OperationType::cause_is_popcount;
		}


		void enqueue(Lit l, bool alter_trail)override{
			importTheory(theory);
			Var v= var(l);
			assert(args.contains(v));
			if(sign(l)){
				over--;
			}else{
				under++;
			}
			assert(over>=0);
			assert(under<=args.size());
			assert(under>=0);
			assert(over<=args.size());
			if(alter_trail){
				//theory.trail.push( { true, !sign(l),getID(), v });


				if(over<over_approx[bvID] || under>under_approx[bvID]){
					theory.bv_needs_propagation[bvID]=true;
					addAlteredBV(bvID);
				}
			}
		}
		virtual void backtrack(Assignment & e, bool rewind){
			importTheory(theory);
			if(e.isOperation){
				assert(e.bvID==getID());
				Var v = e.var;
				assert(v!=var_Undef);
				assert(args.contains(v));
				assert(value(mkLit(v))!=l_Undef);
				if(value(mkLit(v))==l_True){
					under--;
				}else{
					assert(value(mkLit(v))==l_False);
					over++;
				}
				assert(over>=0);
				assert(under<=args.size());
				assert(under>=0);
				assert(over<=args.size());


			}
		}
		bool propagate(bool & changed,vec<Lit> & conflict) override{
			importTheory(theory);
			if(under_approx[bvID]>over){
				//this is a conflict
				buildReason(conflict);

				return false;
			}else if (over_approx[bvID]<under){
				//this is a conflict
				buildReason(conflict);

				return false;
			}else if (under_approx[bvID]==over_approx[bvID] && under<over){
				assert(under<under_approx[bvID] || over>over_approx[bvID]);
				if(under<under_approx[bvID] && over==over_approx[bvID]){
					for(Var v:args){
						if(value(mkLit(v))==l_Undef)
							enqueue(mkLit(v,false),theory.popcount_marker);
					}
					assert(under==over);
				}else if(over>over_approx[bvID] && under==under_approx[bvID]){
					for(Var v:args){
						if(value(mkLit(v))==l_Undef)
							enqueue(mkLit(v,true),theory.popcount_marker);
					}
					assert(under==over);
				}
			}

			return true;
		}
		void updateApprox(Var ignore_bv, Weight & under_new, Weight & over_new, Cause & under_cause_new, Cause & over_cause_new) override{
			importTheory(theory);
#ifndef NDEBUG
			Weight dbg_min=0;
			Weight dbg_max = 0;
			for(Var v:args){
				if(value(mkLit(v))==l_True){
					dbg_min++;
					dbg_max++;
				}else if (value(mkLit(v))!=l_False){
					dbg_max++;
				}
			}
			assert(dbg_min==under);
			assert(dbg_max==over);
#endif
			if (over<over_new){
				over_new=over;
				over_cause_new.clear();
				over_cause_new.type =getType();
				over_cause_new.index=getID();
			}
			if (under>under_new){
				under_new=under;
				under_cause_new.clear();
				under_cause_new.type =getType();
				under_cause_new.index=getID();
			}
		}
		bool checkApproxUpToDate(Weight & under_new,Weight&over_new)override{
			importTheory(theory);

			Weight dbg_min=0;
			Weight dbg_max = 0;
			for(Var v:args){
				if(value(mkLit(v))==l_True){
					dbg_min++;
					dbg_max++;
				}else if (value(mkLit(v))!=l_False){
					dbg_max++;
				}
			}
			assert(dbg_min==under);
			assert(dbg_max==over);

			if (dbg_max<over_new){
				over_new=dbg_max;

			}
			if (dbg_min>under_new){
				under_new=dbg_min;

			}
			return true;
		}
		void buildReason(vec<Lit> & conflict){
			importTheory(theory);
			if(under_approx[bvID]>over){
				//this is a conflict

				for(Var v:args){
					if(value(mkLit(v))==l_False)
						conflict.push(toSolver(mkLit(v,false)));

				}

				theory.buildComparisonReason(Comparison::geq,bvID,under,conflict);
			}else if (over_approx[bvID]<under){
				//this is a conflict

				for(Var v:args){
					if(value(mkLit(v))==l_True)
						conflict.push(toSolver(mkLit(v,true)));

				}

				theory.buildComparisonReason(Comparison::leq,bvID,over,conflict);
			}

		}
		void buildReason(Lit p,CRef marker, vec<Lit> & reason)override{
			importTheory(theory);
			assert(marker==theory.popcount_marker);
			reason.push(toSolver(p));
			if(sign(p)){
				//the reason that p was assigned to false, was that either one of the other true literals had to be false,
				//or bvID had to be larger

				for(Var v:args){
					if(value(mkLit(v))==l_True)
						reason.push(toSolver(mkLit(v,true)));

				}
				//assert(over==over_approx[bvID]);


				theory.buildComparisonReason(Comparison::leq,bvID,over_approx[bvID],reason);
			}else{
				//the reason that p was assigned to true, was that either one of the other false literals had to be true,
				//or bvID had to be smaller

				for(Var v:args){
					if(value(mkLit(v))==l_False)
						reason.push(toSolver(mkLit(v,false)));

				}

				theory.buildComparisonReason(Comparison::geq,bvID,under_approx[bvID],reason);
			}

		}
		void analyzeReason(bool compareOver,Comparison op, Weight  to,  vec<Lit> & conflict) override{
			importTheory(theory);

			Weight  overApprox = over_approx[bvID];
			Weight  underApprox = under_approx[bvID];

			if(compareOver){
				for(Var v:args){
					if(value(mkLit(v))==l_False)
						conflict.push(toSolver(mkLit(v,false)));
				}

			}else{

				for(Var v:args){
					if(value(mkLit(v))==l_True)
						conflict.push(toSolver(mkLit(v,true)));
				}
			}
		}

		bool checkSolved()override{
			importTheory(theory);
			Weight dbg_under=0;
			Weight dbg_over=args.size();
			for(Var v:args){
				if(value(mkLit(v))==l_False)
					dbg_over--;
				else if(value(mkLit(v))==l_True)
					dbg_under++;
			}
			if(dbg_over!=over){
				return false;
			}
			if(dbg_over!=over_approx[bvID])
				return false;

			if(dbg_under!=under){
				return false;
			}
			if(dbg_under!=under_approx[bvID])
				return false;
			return true;
		}
	};

	class ConditionalArg;
	class ConditionalID:public Operation{
			public:
		using Operation::getID;
		using Operation::theory;


	public:

		int bvID=-1;
		const Lit condition;
		ConditionalArg * thenOp;
		ConditionalArg * elseOp;

		int getThenBV()const{
			return thenOp->bvID;
		}
		int getElseBV()const{
			return elseOp->bvID;
		}
		//int thenCID=-1;
		//int elseCID=-1;

		ConditionalID(BVTheorySolver & theory,int bvID, Lit condition):Operation(theory),bvID(bvID),condition(condition){

		}
		int getBV()override{
			return bvID;
		}
		bool hasITE()const{
			return condition!=lit_Undef;
		}
		void setThen(ConditionalArg * thn){
			thenOp=thn;
		}
		void setElse(ConditionalArg * els){
			elseOp=els;
		}
		void move(int bvID) override{

			this->bvID=bvID;
		}

		OperationType getType()const override{
			return OperationType::cause_is_condition;
		}
		void enqueue(Lit l, bool alter_trail)override{
			importTheory(theory);
			Var v= var(l);
			if(alter_trail){

				int bvThenID=getThenBV();
				int bvElseID=getElseBV();
				int bvResultID=bvID;

				assert(bvID==bvResultID);
				assert(var(condition)==var(l));

				int bvCondID;
				if(!sign(l)){
					bvCondID=bvThenID;
				}else{
					bvCondID=bvElseID;
				}

				//theory.trail.push( { false, !sign(l),getID(), v });


				theory.bv_needs_propagation[bvID]=true;
				addAlteredBV(bvID);
				theory.bv_needs_propagation[bvCondID]=true;
				addAlteredBV(bvCondID);
			}

		}
		bool propagate(bool & changed_outer,vec<Lit> & conflict) override{
			importTheory(theory);
			Weight & underApprox = under_approx[bvID];
			Weight & overApprox = over_approx[bvID];
			//assert(bvThenID<bvID);
			//assert(bvElseID<bvID);
			int bvThenID=thenOp->bvID;
			int bvElseID=elseOp->bvID;
			bool changed=true;
			while(changed){
				changed=false;

				Weight under;
				Weight over;
				if (value(condition)==l_True){
					under = under_approx[bvThenID];
					over = over_approx[bvThenID];
				}else if (value(condition)==l_False){
					under = under_approx[bvElseID];
					over = over_approx[bvElseID];
				}else{
					under = min( under_approx[bvThenID], under_approx[bvElseID]);
					over = max( over_approx[bvThenID], over_approx[bvElseID]);
				}

				if(underApprox>over){
					//then we have a conflict
					double startconftime = rtime(2);
					//propagationtime += startconftime - startproptime;
					theory.stats_num_conflicts++;theory.stats_addition_conflicts++;
					if(opt_verb>1){
						printf("bv condition conflict %ld\n", theory.stats_num_conflicts);
					}
					if(value(condition)==l_True){
						conflict.push(toSolver(~condition));
					}else if(value(condition)==l_False){
						conflict.push(toSolver(condition));
					}
					buildReason(conflict);

					theory.stats_conflict_time+=rtime(2)-startconftime;
					return false;
				}else if (overApprox<under){
					double startconftime = rtime(2);
					//propagationtime += startconftime - startproptime;
					if(opt_verb>1){
						printf("bv condition conflict %ld\n", theory.stats_num_conflicts);
					}
					theory.stats_num_conflicts++;theory.stats_addition_conflicts++;
					if(value(condition)==l_True){
						conflict.push(toSolver(~condition));
					}else if(value(condition)==l_False){
						conflict.push(toSolver(condition));
					}
					buildReason(conflict);

					theory.stats_conflict_time+=rtime(2)-startconftime;
					return false;
				}

				if(value(condition)==l_True){
					if((underApprox> under_approx[bvThenID]) || (overApprox < over_approx[bvThenID])){
						//the other bv needs to be updated
						addAlteredBV(bvThenID);

					}
				}else if (value(condition)==l_False){
					if((underApprox> under_approx[bvElseID]) || (overApprox < over_approx[bvElseID])){
						//the other bv needs to be updated
						addAlteredBV(bvElseID);
					}
				}else{
					//we can propagate the value of condition, based on the values of its arguments:
					//if thenID is out of range of result ID, then condition must be false.
					//if elseID is out of range of result ID, then condition must be true
					//(if both are the case, then we would already have hit a conflict above(?))

					if(under_approx[bvThenID]>over_approx[bvID] || over_approx[bvThenID]<under_approx[bvID]){
						//condition must be false

						if(under_approx[bvElseID]>over_approx[bvID] || over_approx[bvElseID]<under_approx[bvID]){
							//condition must be both true and false, which is a conflict.
							buildReason(conflict);

							//buildConditionalReason(bvID,bvElseID,underApprox,conflict);
							return false;
						}else{
							enqueue(~condition,  theory.conditionelse_prop_marker);
							changed=true;
						}

					}else if(under_approx[bvElseID]>over_approx[bvID] || over_approx[bvElseID]<under_approx[bvID]){
						enqueue(condition,theory.conditionthen_prop_marker);
						changed=true;
						//condition must be true
					}
				}
			}


			return true;
		}
		void updateApprox(Var ignore_bv, Weight & under_new, Weight & over_new, Cause & under_cause_new, Cause & over_cause_new) override{
			importTheory(theory);
			int bvThenID=thenOp->bvID;
			int bvElseID=elseOp->bvID;

			Weight under;
			Weight over;
			if (value(condition)==l_True){
				under = under_approx[bvThenID];
				over = over_approx[bvThenID];
			}else if (value(condition)==l_False){
				under = under_approx[bvElseID];
				over = over_approx[bvElseID];
			}else{
				under = min( under_approx[bvThenID], under_approx[bvElseID]);
				over = max( over_approx[bvThenID], over_approx[bvElseID]);
			}
		/*	clip_over(under,bvID);
			clip_over(over,bvID);*/
			if(under >under_new){
				under_new=under;
				under_cause_new.clear();
				under_cause_new.type =getType();
				under_cause_new.index=getID();
			}
			if(over<over_new){
				over_new=over;
				over_cause_new.clear();
				over_cause_new.type =getType();
				over_cause_new.index=getID();
			}

		}
		void buildReason( vec<Lit> & conflict){
			importTheory(theory);
			theory.dbg_no_pending_analyses();
			assert(theory.eq_bitvectors[bvID]==bvID);

			theory.stats_build_condition_reason++;
			Weight  over_cur = over_approx[bvID];
			Weight  under_cur = under_approx[bvID];
			//assert(checkApproxUpToDate(bvID));


			//the reason that the addition is over is the reason that

			int bvThenID=thenOp->bvID;
			int bvElseID=elseOp->bvID;

			//assert(bvThenID<bvID);
			//assert(bvElseID<bvID);

			Weight under;
			Weight over;
			if (value(condition)==l_True){
				under = under_approx[bvThenID];
				over = over_approx[bvThenID];
				if(under_cur>over){
					analyzeValueReason(Comparison::gt, bvID,over,conflict);
					analyzeValueReason(Comparison::leq, bvThenID,over,conflict);
				}else{
					assert(over_cur<under);
					analyzeValueReason(Comparison::lt, bvID,under,conflict);
					analyzeValueReason(Comparison::geq, bvThenID,under,conflict);
				}
			}else if (value(condition)==l_False){
				under = under_approx[bvElseID];
				over = over_approx[bvElseID];
				if(under_cur>over){
					analyzeValueReason(Comparison::gt, bvID,over,conflict);
					analyzeValueReason(Comparison::leq, bvElseID,over,conflict);
				}else{
					assert(over_cur<under);
					analyzeValueReason(Comparison::lt, bvID,under,conflict);
					analyzeValueReason(Comparison::geq, bvElseID,under,conflict);
				}
			}else{
				assert((under_approx[bvThenID]>over_approx[bvID] || over_approx[bvThenID]<under_approx[bvID])&&
						(under_approx[bvElseID]>over_approx[bvID] || over_approx[bvElseID]<under_approx[bvID]));
				if(under_approx[bvThenID]>over_approx[bvID]){
					analyzeValueReason(Comparison::lt, bvID,under_approx[bvThenID],conflict);
					analyzeValueReason(Comparison::geq, bvThenID,under_approx[bvThenID],conflict);
				}else if (over_approx[bvThenID]<under_approx[bvID]){
					analyzeValueReason(Comparison::gt, bvID,over_approx[bvThenID],conflict);
					analyzeValueReason(Comparison::leq, bvThenID,over_approx[bvThenID],conflict);
				}else{
					assert(false);
				}
				if(under_approx[bvElseID]>over_approx[bvID]){
					analyzeValueReason(Comparison::lt, bvID,under_approx[bvElseID],conflict);
					analyzeValueReason(Comparison::geq, bvElseID,under_approx[bvElseID],conflict);
				}else if ( over_approx[bvElseID]<under_approx[bvID]){
					analyzeValueReason(Comparison::gt, bvID,over_approx[bvElseID],conflict);
					analyzeValueReason(Comparison::leq, bvElseID,over_approx[bvElseID],conflict);
				}else{
					assert(false);
				}
			}

			analyze(conflict);
		}
		void buildConditionalPropReason( vec<Lit> & conflict){
			importTheory(theory);
			theory.dbg_no_pending_analyses();
			assert(theory.eq_bitvectors[bvID]==bvID);

			theory.stats_build_condition_reason++;
			Weight  over_cur = over_approx[bvID];
			Weight  under_cur = under_approx[bvID];
			//assert(checkApproxUpToDate(bvID));


			//the reason that the addition is over is the reason that
			//bvID > addition_under, or the reason that addition_under>= its current value.
			int bvThenID = getThenBV();
			int bvElseID = getElseBV();

			//assert(bvThenID<bvID);
			//assert(bvElseID<bvID);
			assert(value(condition)!=l_Undef);
			Weight under;
			Weight over;
			if (value(condition)==l_False){
				//the reason the condition lit was set to false, was because the 'then' bv was out of the feasible range.
				under = under_approx[bvThenID];
				over = over_approx[bvThenID];
				if(under_cur>over){
					analyzeValueReason(Comparison::gt, bvID,over,conflict);
					analyzeValueReason(Comparison::leq, bvThenID,over,conflict);
				}else{
					assert(over_cur<under);
					analyzeValueReason(Comparison::lt, bvID,under,conflict);
					analyzeValueReason(Comparison::geq, bvThenID,under,conflict);
				}
			}else if (value(condition)==l_True){
				under = under_approx[bvElseID];
				over = over_approx[bvElseID];
				if(under_cur>over){
					analyzeValueReason(Comparison::gt, bvID,over,conflict);
					analyzeValueReason(Comparison::leq, bvElseID,over,conflict);
				}else{
					assert(over_cur<under);
					analyzeValueReason(Comparison::lt, bvID,under,conflict);
					analyzeValueReason(Comparison::geq, bvElseID,under,conflict);
				}
			}else{
				assert(false);
			}

			analyze(conflict);
		}

		void buildReason(Lit p, CRef marker, vec<Lit> & reason)override{
			importTheory(theory);
			if (marker == theory.conditionthen_prop_marker) {
				//the 'then' condition must hold if the else condition is out of range
				assert(value(p)==l_True);
				reason.push(toSolver(p));
				Var v = var(p);

				theory.writeBounds(bvID);

				int bvThenID = getThenBV();
				int bvElseID = getElseBV();

				assert(under_approx[bvElseID]>over_approx[bvID] || over_approx[bvElseID]<under_approx[bvID]);
				buildConditionalPropReason(reason);
			}else if (marker ==theory.conditionelse_prop_marker){
				assert(value(p)==l_True);
				reason.push(toSolver(p));
				Var v = var(p);

				theory.writeBounds(bvID);


				int bvThenID = getThenBV();
				int bvElseID = getElseBV();
				assert(under_approx[bvThenID]>over_approx[bvID] || over_approx[bvThenID]<under_approx[bvID]);

				buildConditionalPropReason(reason);
			}else if(marker == theory.conditionarg_prop_marker){
				assert(value(p)==l_True);
				reason.push(toSolver(p));
				Var v = var(p);


				int resultID = bvID;
				int bvThenID = getThenBV();
				int bvElseID = getElseBV();

				if(p ==condition){
					ConditionalArg * op = elseOp;//the reason is that the else condition was falsified.
					assert(under_approx[bvElseID]>over_approx[resultID] || over_approx[bvElseID]<under_approx[resultID]);
					op->buildReason(reason);
					//buildConditionalArgReason(bvElseID,opID,reason);
				}else{
					assert(p==~condition);
					ConditionalArg * op = thenOp;
					assert(under_approx[bvThenID]>over_approx[resultID] || over_approx[bvThenID]<under_approx[resultID]);
					//buildConditionalArgReason(bvThenID,opID,reason);
					op->buildReason(reason);
				}
			}

		}

		bool checkApproxUpToDate(Weight & under,Weight&over)override{
			importTheory(theory);
			int thenID=thenOp->bvID;
			int elseID=elseOp->bvID;


			if(value(condition)==l_Undef){
				if(under_approx[bvID]< min(under_approx[thenID],under_approx[elseID]) ||
						over_approx[bvID]> max(over_approx[thenID],over_approx[elseID])){
					return false;
				}
				Weight underc = min( under_approx[thenID], under_approx[elseID]);
				Weight overc = max( over_approx[thenID], over_approx[elseID]);
				if(underc >under){
					under=underc;
				}
				if(overc<over){
					over=overc;
				}
			}else if(value(condition)==l_True){
				if(under_approx[bvID]!= under_approx[thenID] || over_approx[bvID]!=over_approx[thenID]){
					return false;
				}
				if(under_approx[thenID] >under){
					under=under_approx[thenID];
				}
				if(over_approx[thenID]<over){
					over=over_approx[thenID];
				}
			}else if(value(condition)==l_False){
				if(under_approx[bvID]!= under_approx[elseID] || over_approx[bvID]!=over_approx[elseID]){
					return false;
				}
				if(under_approx[elseID] >under){
					under=under_approx[elseID];
				}
				if(over_approx[elseID]<over){
					over=over_approx[elseID];
				}
			}

			return true;
		}

		void analyzeReason(bool compareOver,Comparison op, Weight  to,  vec<Lit> & conflict){
			importTheory(theory);
			int bvThenID=thenOp->bvID;
			int bvElseID=elseOp->bvID;

			if(compareOver){

				//if the condition is assigned (at this point in the trail)
				if (value(condition)==l_True ){
					conflict.push(toSolver(~condition));
					assert(over_approx[bvThenID]<=over_approx[bvID]);//this is the cause
					addAnalysis(op,bvThenID,to);
				}else if (value(condition)==l_False ){
					conflict.push(toSolver(condition));
					assert(over_approx[bvElseID]<=over_approx[bvID]);//this is the cause
					addAnalysis(op,bvElseID,to);
				}else{//else if (over_approx[bvThenID]>=over_approx[bvElseID]){
					assert(over_approx[bvThenID]<=over_approx[bvID]);
					assert( over_approx[bvElseID]<=over_approx[bvID]);//this is the cause
					addAnalysis(op,bvThenID,to);
					addAnalysis(op,bvElseID,to);
				}
			}else{

				//if the condition is assigned (at this point in the trail)
				if (value(condition)==l_True ){
					conflict.push(toSolver(~condition));
					assert(under_approx[bvThenID]>=under_approx[bvID]);//this is the cause
					addAnalysis(op,bvThenID,to);
				}else if (value(condition)==l_False){
					conflict.push(toSolver(condition));
					assert(under_approx[bvElseID]>=under_approx[bvID]);//this is the cause
					addAnalysis(op,bvElseID,to);
				}else{//if ( under_approx[bvThenID]<=under_approx[bvElseID]){
					assert(under_approx[bvThenID]>=under_approx[bvID]);
					assert(under_approx[bvElseID]>=under_approx[bvID]);
					addAnalysis(op,bvThenID,to);
				//}else if ( under_approx[bvElseID]<under_approx[bvThenID]){

					addAnalysis(op,bvElseID,to);
				}
			}
		}

		bool checkSolved()override{
			importTheory(theory);
			int thenID=thenOp->bvID;
			int elseID=elseOp->bvID;

			if(value(condition)!=l_False){
				if(under_approx[bvID]!= under_approx[thenID] || over_approx[bvID]!=over_approx[thenID]){
					return false;
				}
			}
			if(value(condition)!=l_True){
				if(under_approx[bvID]!= under_approx[elseID] || over_approx[bvID]!=over_approx[elseID]){
					return false;
				}
			}
			return true;
		}

	};

	class ConditionalArg:public Operation{
	public:
		using Operation::getID;
		using Operation::theory;

	public:

		int bvID=-1;
		const Lit condition;
		ConditionalArg * otherOp;
		ConditionalID * resultOp;
		int getBV()override{
			return bvID;
		}
		ConditionalArg(BVTheorySolver & mtheory,int bvID,Lit condition,ConditionalID * result):Operation(mtheory),bvID(bvID),condition(condition),resultOp(result){

		}
		void setOtherArg(ConditionalArg * otherOp){
			this->otherOp=otherOp;
		}

		OperationType getType()const override{
			return OperationType::cause_is_condition_argument;
		}

		void move(int bvID) override{

			this->bvID=bvID;
		}

		bool propagate( bool & changed_outer,vec<Lit> & conflict) override{
			importTheory(theory);
			int other_argID=otherOp->bvID;
			int resultID=resultOp->bvID;
			Weight & underApprox = under_approx[bvID];
			Weight & overApprox = over_approx[bvID];


			assert(other_argID>=0);assert(resultID>=0);

			if(value(condition)==l_True){
				if(under_approx[bvID]>over_approx[resultID] || over_approx[bvID]<under_approx[resultID]){
					//this is a conflict
					conflict.push(toSolver(~condition));
					buildReason(conflict);

					return false;
				}else if(under_approx[bvID]>under_approx[resultID] || over_approx[bvID]<over_approx[resultID]){
					addAlteredBV(resultID);
				}
			}else if (value(condition)==l_Undef){
				//can potentially propagate the falsehood of the conditional, if this bv is out of range of the result bv
				if(under_approx[bvID]>over_approx[resultID] || over_approx[bvID]<under_approx[resultID]){
					enqueue(~condition,theory.conditionarg_prop_marker);

					addAlteredBV(other_argID);
					addAlteredBV(resultID);
				}else{
					Weight under_min = min(under_approx[bvID],under_approx[other_argID]);
					Weight over_max = max(over_approx[bvID],over_approx[other_argID]);

					if (under_min>under_approx[resultID] || over_max<over_approx[resultID]){
						//ensure that it is a change to this bv's bounds that caused the change here...
						//assert(under_approx[bvID]>under_approx[other_argID] || over_approx[bvID]<over_approx[resultID]);

						addAlteredBV(resultID);
					}

				}
			}
			return true;

		}
		void updateApprox(Var ignore_bv, Weight & under_new, Weight & over_new, Cause & under_cause_new, Cause & over_cause_new) override{
			importTheory(theory);
			int other_argID=otherOp->bvID;
			int resultID=resultOp->bvID;


			assert(other_argID>=0);assert(resultID>=0);
			if(value(condition)==l_True){
				Weight under = under_approx[resultID];
				Weight over = over_approx[resultID];
	/*			clip_under(under,bvID);
				clip_under(over,bvID);*/
				if(under >under_new){
					under_new=under;
					under_cause_new.clear();
					under_cause_new.type =getType();
					under_cause_new.index=getID();
				}
				if(over<over_new){
					over_new=over;
					over_cause_new.clear();
					over_cause_new.type =getType();
					over_cause_new.index=getID();
				}
			}
		}

		void buildReason(vec<Lit> & conflict){
			importTheory(theory);
			int other_argID=otherOp->bvID;
			int resultID=resultOp->bvID;

			theory.dbg_no_pending_analyses();
			assert(theory.eq_bitvectors[bvID]==bvID);


			theory.stats_build_condition_arg_reason++;
			Weight  over_cur = over_approx[bvID];
			Weight  under_cur = under_approx[bvID];
			//assert(checkApproxUpToDate(bvID));


			assert(value(condition)!=l_Undef);

			Weight under_result = under_approx[resultID];
			Weight over_result = over_approx[resultID];

			if(under_cur>over_result){
				//buildTrivialClause(conflict);
				analyzeValueReason(Comparison::leq, resultID,over_approx[resultID],conflict);
				analyzeValueReason(Comparison::gt, bvID,over_result,conflict);
			}else{
				assert(over_cur<under_result);

				analyzeValueReason(Comparison::geq, resultID,under_approx[resultID],conflict);
				analyzeValueReason(Comparison::lt, bvID,under_result,conflict);

			}
			analyze(conflict);
		}

		void analyzeReason(bool compareOver,Comparison op, Weight  to,  vec<Lit> & conflict){
			importTheory(theory);
			int other_argID=otherOp->bvID;
			int resultID=resultOp->bvID;

			if(compareOver){
				Weight over_sum = over_approx[resultID];
				assert(over_sum<=over_approx[bvID]);
				assert(value(condition)==l_True);
				conflict.push(toSolver(~condition));
				addAnalysis(op,resultID,to);
			}else{
				Weight under_sum = under_approx[resultID];
				assert(under_sum>=under_approx[bvID]);
				assert(value(condition)==l_True);

				conflict.push(toSolver(~condition));
				addAnalysis(op,resultID,to);
			}

		}

		bool checkApproxUpToDate(Weight & under,Weight&over)override{
			importTheory(theory);
			int other_argID=otherOp->bvID;
			int resultID=resultOp->bvID;

			assert(other_argID>=0);assert(resultID>=0);
			if(value(condition)==l_True){
				Weight underc = under_approx[resultID];
				Weight overc = over_approx[resultID];
	/*			clip_under(under,bvID);
				clip_under(over,bvID);*/
				if(underc >under_approx[bvID]){
					return false;
				}
				if(overc< over_approx[bvID]){
					return false;
				}
				if(under_approx[resultID] >under){
					under=under_approx[resultID];
				}
				if(over_approx[resultID]<over){
					over=over_approx[resultID];
				}
			}
			return true;
		}

		bool checkSolved()override{
			importTheory(theory);
			int other_argID=otherOp->bvID;
			int resultID=resultOp->bvID;

			if(value(condition)!=l_False){
				if(under_approx[resultID]!= under_approx[bvID] || over_approx[resultID]!=over_approx[bvID]){
					return false;
				}
			}
			if(value(condition)!=l_True){
				if(under_approx[resultID]!= under_approx[other_argID] || over_approx[resultID]!=over_approx[other_argID]){
					return false;
				}
			}
			return true;
		}
	};

	class AdditionArg;
	class Addition:public Operation{
	public:
		using Operation::getID;
		using Operation::theory;


	public:
		AdditionArg * arg1=nullptr;
		AdditionArg * arg2=nullptr;
		int bvID=-1;

		Addition(BVTheorySolver & theory,int bvID):Operation(theory),bvID(bvID){

		}
		int getBV()override{
			return bvID;
		}
		bool hasAddition()const{
			return arg1;
		}

		void setArg1(AdditionArg * arg1){
			this->arg1=arg1;
		}
		void setArg2(AdditionArg * arg2){
			this->arg2=arg2;
		}


		void move( int bvID) override{

			this->bvID=bvID;
		}



		OperationType getType()const override{
			return OperationType::cause_is_addition;
		}

		bool propagate(bool & changed_outer,vec<Lit> & conflict) override{
			importTheory(theory);
			int aID=arg1->bvID;
			int bID=arg2->bvID;
			Weight & underApprox = under_approx[bvID];
			Weight & overApprox = over_approx[bvID];

				//assert(aID<bvID);
				//assert(bID<bvID);
				Weight under = under_approx[aID] +  under_approx[bID];
				Weight over = over_approx[aID] +  over_approx[bID];
				clip_over(under,bvID);
				clip_over(over,bvID);
				if(underApprox>over){
					//then we have a conflict
					double startconftime = rtime(2);
					//propagationtime += startconftime - startproptime;
					theory.stats_num_conflicts++;
					theory.stats_addition_conflicts++;
					if(opt_verb>1){
						printf("bv addition conflict %ld\n", theory.stats_num_conflicts);
					}
					buildReason(conflict);

					theory.stats_conflict_time+=rtime(2)-startconftime;
					return false;
				}else if (overApprox<under){
					double startconftime = rtime(2);
					//propagationtime += startconftime - startproptime;
					if(opt_verb>1){
						printf("bv addition conflict %ld\n", theory.stats_num_conflicts);
					}
					theory.stats_num_conflicts++;
					theory.stats_addition_conflicts++;
					buildReason(conflict);

					theory.stats_conflict_time+=rtime(2)-startconftime;
					return false;
				}
				Weight under_arg_b = underApprox - over_approx[bID];
				Weight over_arg_b = overApprox - under_approx[bID];
				clip_under(under_arg_b,bvID);
				clip_under(over_arg_b,bvID);
				//this check may be especially important when either aID or bID is really a constant...
				if((under_arg_b> under_approx[aID]) || (over_arg_b < over_approx[aID])){
					//the other bv needs to be updated
					addAlteredBV(aID);
				}
				Weight under_arg_a =underApprox - over_approx[aID];
				Weight over_arg_a =  overApprox - under_approx[aID];
				clip_under(under_arg_a,bvID);
				clip_under(over_arg_a,bvID);
				if((under_arg_a > under_approx[bID]) || (over_arg_a< over_approx[bID])){
					//the other bv needs to be updated
					addAlteredBV(bID);
				}
			return true;

		}
		void updateApprox(Var ignore_bv, Weight & under_new, Weight & over_new, Cause & under_cause_new, Cause & over_cause_new) override{
			importTheory(theory);
			int aID=arg1->bvID;
			int bID=arg2->bvID;
			//assert(aID<bvID);
			//assert(bID<bvID);
			Weight under = under_approx[aID] +  under_approx[bID];
			Weight over = over_approx[aID] +  over_approx[bID];
			clip_over(under,bvID);
			clip_over(over,bvID);
			if(under >under_new){
				under_new=under;
				under_cause_new.clear();
				under_cause_new.type =getType();
				under_cause_new.index=getID();
			}
			if(over<over_new){
				over_new=over;
				over_cause_new.clear();
				over_cause_new.type =getType();
				over_cause_new.index=getID();
			}
		}

		void buildReason(vec<Lit> & conflict){
			importTheory(theory);
			int aID=arg1->bvID;
			int bID=arg2->bvID;
			theory.dbg_no_pending_analyses();
			assert(theory.eq_bitvectors[bvID]==bvID);
			//rewind_trail_pos(trail.size()-1);
			theory.stats_build_addition_reason++;
			Weight  over_cur = over_approx[bvID];
			Weight  under_cur = under_approx[bvID];
			//assert(theory.checkApproxUpToDate(bvID));


			//the reason that the addition is over is the reason that
			//bvID > addition_under, or the reason that addition_under>= its current value.

			//assert(aID<bvID);
			//assert(bID<bvID);

			Weight under_add = under_approx[aID] +  under_approx[bID];
			Weight over_add = over_approx[aID] +  over_approx[bID];

			int width = theory.bitvectors[bvID].size();
			Weight max_val = ((1L)<<width)-1;
			if(under_add>max_val){
				under_add=max_val;
			}
			if(over_add>max_val){
				over_add=max_val;
			}

			if(under_cur>over_add){

				analyzeValueReason(Comparison::gt, bvID,over_add,conflict);

				analyzeValueReason(Comparison::leq, aID,over_approx[aID],conflict);
				analyzeValueReason(Comparison::leq, bID,over_approx[bID],conflict);
			}else{
				assert(over_cur<under_add);
				analyzeValueReason(Comparison::lt, bvID,under_add,conflict);

				analyzeValueReason(Comparison::geq, aID,under_approx[aID],conflict);
				analyzeValueReason(Comparison::geq, bID,under_approx[bID],conflict);
			}
			analyze(conflict);
		}

		void analyzeReason(bool compareOver,Comparison op, Weight  to,  vec<Lit> & conflict){
			importTheory(theory);
			if(compareOver){

				int aID=arg1->bvID;
				int bID=arg2->bvID;
				//assert(aID<bvID);
				//assert(bID<bvID);
				Weight over_bid = over_approx[bID];
				Weight over_aid = over_approx[aID];
				//then the reason is that aID is <= weight-under(bID), or bID <= weight-under(aID)
				//addAnalysis(Comparison::leq,aID,over_aid,conflict);//to-over_bid
				//addAnalysis(Comparison::leq,bID,over_bid,conflict);//to-over_aid

				addAnalysis(Comparison::leq,aID,over_approx[bvID]-over_bid);
				addAnalysis(Comparison::leq,bID,over_approx[bvID]-over_aid);
			}else{


				int aID=arg1->bvID;
				int bID=arg2->bvID;
				Weight under_bid = under_approx[bID];
				Weight under_aid = under_approx[aID];
				//addAnalysis(Comparison::geq,aID,under_aid,conflict);
				//addAnalysis(Comparison::geq,bID,under_bid,conflict);
				addAnalysis(Comparison::geq,aID,under_approx[bvID]-under_bid);
				addAnalysis(Comparison::geq,bID,under_approx[bvID]-under_aid);
				//buildValueReason(op,aID,to-under_bid,conflict,trail_pos-1);
				//buildValueReason(op,bID,to-under_aid,conflict,trail_pos-1);
			}
		}

		bool checkApproxUpToDate(Weight & under,Weight&over)override{
			importTheory(theory);
			int aID=arg1->bvID;
			int bID=arg2->bvID;
			//assert(aID<bvID);
			//assert(bID<bvID);
			Weight underadd = under_approx[aID] +  under_approx[bID];
			Weight overadd = over_approx[aID] +  over_approx[bID];
			if(underadd >under){
				under=underadd;
			}
			if(overadd<over){
				over=overadd;
			}
			return true;
		}

		bool checkSolved()override{
			importTheory(theory);
			int aID=arg1->bvID;
			int bID=arg2->bvID;
			int width = theory.bitvectors[bvID].size();
			Weight max_val = (1L<<width)-1;
			//assert(aID<bvID);
			//assert(bID<bvID);
			Weight underadd = under_approx[aID] +  under_approx[bID];
			Weight overadd = over_approx[aID] +  over_approx[bID];
			if(underadd>max_val){
				underadd=max_val;
			}
			if(overadd>max_val){
				overadd=max_val;
			}
			if(underadd >under_approx[bvID]){
				return false;
			}
			if(overadd<over_approx[bvID]){
				return false;
			}
			return true;
		}
	};


	class AdditionArg:public Operation{
	public:
		using Operation::getID;
		using Operation::theory;

		//int other_argID=-1;
		//int sumID=-1;
	public:
		AdditionArg * otherOp;
		Addition * resultOp;
		int bvID=-1;

		AdditionArg(BVTheorySolver & theory,int bvID, Addition* result):Operation(theory),bvID(bvID),resultOp(result){

		}
		int getBV()override{
			return bvID;
		}
		void setOtherArg(AdditionArg * otherArg){
			this->otherOp=otherArg;
		}

		void move( int bvID) override{

			this->bvID=bvID;
		}
		OperationType getType()const override{
			return OperationType::cause_is_addition_argument;
		}

		bool propagate(bool & changed_outer,vec<Lit> & conflict) override{
			importTheory(theory);
			int other_argID=otherOp->bvID;
			int sumID=resultOp->bvID;
			Weight & underApprox = under_approx[bvID];
			Weight & overApprox = over_approx[bvID];

			assert(other_argID>=0);assert(sumID>=0);

			Weight under = under_approx[sumID] -  over_approx[other_argID];
			Weight over = over_approx[sumID] -  under_approx[other_argID];
			clip_under(under,bvID);
			clip_under(over,bvID);
			if(underApprox>over){
				//then we have a conflict
				double startconftime = rtime(2);
				//propagationtime += startconftime - startproptime;
				theory.stats_num_conflicts++;theory.stats_addition_conflicts++;

				if(opt_verb>1){
					printf("bv addition arg conflict %ld\n", theory.stats_num_conflicts);
				}
				buildReason(conflict);

				theory.stats_conflict_time+=rtime(2)-startconftime;
				return false;
			}else if (overApprox<under){
				double startconftime = rtime(2);
				//propagationtime += startconftime - startproptime;
				theory.stats_num_conflicts++;theory.stats_addition_conflicts++;
				if(opt_verb>1){
					printf("bv addition arg conflict %ld\n", theory.stats_num_conflicts);
				}
				buildReason(conflict);

				theory.stats_conflict_time+=rtime(2)-startconftime;
				return false;
			}
			Weight under_arg = underApprox + under_approx[other_argID];
			clip_over(under_arg,bvID);
			Weight over_arg =  overApprox + over_approx[other_argID];
			clip_over(over_arg,bvID);
			//this check may be especially important when either aID or bID is really a constant...
			if((under_arg > under_approx[sumID]) || (over_arg< over_approx[sumID])){
				//the other bv needs to be updated

				addAlteredBV(sumID);
			}
			Weight under_sum =  under_approx[sumID]  - overApprox;
			clip_under(under_sum,bvID);
			Weight over_sum =  over_approx[sumID] - underApprox;
			clip_under(over_sum,bvID);
			if( (under_sum> under_approx[other_argID]) || (over_sum< over_approx[other_argID])){
				//the other bv needs to be updated

				addAlteredBV(other_argID);
			}

			return true;

		}
		void updateApprox(Var ignore_bv, Weight & under_new, Weight & over_new, Cause & under_cause_new, Cause & over_cause_new) override{
			importTheory(theory);
			int other_argID=otherOp->bvID;
			int sumID=resultOp->bvID;

			assert(other_argID>=0);assert(sumID>=0);
			//assert((other_argID!=bvID &&   under_approx[sumID] >=  under_approx[other_argID] + under_old  ) || (other_argID==bvID &&   under_approx[sumID] >= under_old + under_old  ));
			//assert((other_argID!=bvID &&  over_approx[sumID] <=  over_approx[other_argID] + over_old ) || (other_argID==bvID &&  over_approx[sumID] <=  over_old + over_old ));
			Weight under = under_approx[sumID] -  over_approx[other_argID];
			Weight over = over_approx[sumID] -  under_approx[other_argID];
			clip_under(under,bvID);
			clip_under(over,bvID);
			if(under >under_new){
				under_new=under;
				under_cause_new.clear();
				under_cause_new.type =OperationType::cause_is_addition_argument;
				under_cause_new.index=getID();
			}
			if(over<over_new){
				over_new=over;
				over_cause_new.clear();
				over_cause_new.type =OperationType::cause_is_addition_argument;
				over_cause_new.index=getID();
			}


		}

		void buildReason(vec<Lit> & conflict){
			importTheory(theory);
			int other_argID=otherOp->bvID;
			int sumID=resultOp->bvID;
			theory.dbg_no_pending_analyses();
			assert(theory.eq_bitvectors[bvID]==bvID);
			//rewind_trail_pos(trail.size()-1);

			theory.stats_build_addition_arg_reason++;
			Weight  over_cur = over_approx[bvID];
			Weight  under_cur = under_approx[bvID];
			//assert(checkApproxUpToDate(bvID));


			//the reason that the addition is over is the reason that
			//bvID > addition_under, or the reason that addition_under>= its current value.



			Weight under_add = under_approx[sumID] -  over_approx[other_argID];
			Weight over_add = over_approx[sumID] -  under_approx[other_argID];

			int width = theory.bitvectors[bvID].size();
			Weight max_val = ((1L)<<width)-1;
			if(under_add>max_val){
				under_add=max_val;
			}
			if(over_add>max_val){
				over_add=max_val;
			}

			if(under_cur>over_add){
				//buildTrivialClause(conflict);
				analyzeValueReason(Comparison::leq, sumID,over_approx[sumID],conflict);

				analyzeValueReason(Comparison::geq, other_argID,under_approx[other_argID],conflict);
				analyzeValueReason(Comparison::gt, bvID,over_add,conflict);
			}else{
				assert(over_cur<under_add);
				//buildTrivialClause(conflict);
				analyzeValueReason(Comparison::geq, sumID,under_approx[sumID],conflict);
				analyzeValueReason(Comparison::leq, other_argID,over_approx[other_argID],conflict);
				analyzeValueReason(Comparison::lt, bvID,under_add,conflict);

			}
			analyze(conflict);
		}

		void analyzeReason(bool compareOver,Comparison op, Weight  to,  vec<Lit> & conflict){
			importTheory(theory);
			if(compareOver){

				int other_argID=otherOp->bvID;
				int sumID=resultOp->bvID;
				Weight over_sumID = over_approx[sumID];
				Weight under_argID = under_approx[other_argID];

				//Weight over = over_approx[sumID] -  under_approx[other_argID];

				addAnalysis(Comparison::geq,other_argID,over_sumID-over_approx[bvID]);
				addAnalysis(Comparison::leq,sumID,over_approx[bvID]+under_argID);
				//buildValueReason(~op,other_argID,over_sumID-to,conflict,trail_pos-1);
				//buildValueReason(op,sumID,to+under_argID,conflict,trail_pos-1);

			}else{

				int other_argID=otherOp->bvID;
				int sumID=resultOp->bvID;
				Weight under_sumID = under_approx[sumID];
				Weight over_argID = over_approx[other_argID];
				//Weight under = under_approx[sumID] -  over_approx[other_argID];
				addAnalysis(Comparison::leq,other_argID,under_sumID-under_approx[bvID]);
				addAnalysis(Comparison::geq,sumID,under_approx[bvID]+over_argID);
				//buildValueReason(~op,other_argID,under_sumID-to,conflict,trail_pos-1);
				//buildValueReason(op,sumID,to+over_argID,conflict,trail_pos-1);
			}
		}

		bool checkApproxUpToDate(Weight & under,Weight&over)override{
			importTheory(theory);
			int other_argID=otherOp->bvID;
			int sumID=resultOp->bvID;

			Weight under_add = under_approx[sumID] -  over_approx[other_argID];
			Weight over_add = over_approx[sumID] -  under_approx[other_argID];

			if(under_add >under){
				under=under_add;
			}
			if(over_add<over){
				over=over_add;
			}
			return true;
		}
		bool checkSolved()override{
			importTheory(theory);
			int other_argID=otherOp->bvID;
			int sumID=resultOp->bvID;
			int width = theory.bitvectors[sumID].size();
			Weight max_val = (1L<<width)-1;
			//assert(aID<bvID);
			//assert(bID<bvID);
			Weight underadd = under_approx[bvID] +  under_approx[other_argID];
			Weight overadd = over_approx[bvID] +  over_approx[other_argID];
			if(underadd>max_val){
				underadd=max_val;
			}
			if(overadd>max_val){
				overadd=max_val;
			}
			if(underadd >under_approx[sumID]){
				return false;
			}
			if(overadd<over_approx[sumID]){
				return false;
			}
			return true;
		}
	};
	class MinMaxArg;
	class MinMaxData:public Operation{
	public:
		using Operation::getID;
		using Operation::theory;


	public:
		const bool min;
		int bvID;
		vec<MinMaxArg*> args;


		MinMaxData(BVTheorySolver & theory, int bvID, bool min):Operation(theory),min(min),bvID(bvID){

		}
		int getBV()override{
			return bvID;
		}
		void move( int bvID) override{
			this->bvID=bvID;
		}
		OperationType getType()const override{
			return OperationType::cause_is_minmax;
		}

		void addArgument(MinMaxArg * arg){
			args.push(arg);
		}

		bool propagate( bool & changed_outer,vec<Lit> & conflict) override{
			importTheory(theory);
			Weight & underApprox = under_approx[bvID];
			Weight & overApprox = over_approx[bvID];

			bool isMin = min;

			//If any argument has an upper bound larger than bvID's, then
			//need to update that arg
			assert(args.size()>0);
			//there are two classes of conflicts here (for the max case).
			//1) the underapprox of bvID is higher than the highest overapprox of its arguments
			//then the conflict reason is that at least one of them must increase, or bvID must decrease.
			Weight minmax_arg_under = under_approx[bvID];
			Weight minmax_arg_over = over_approx[bvID];
			for(int i = 0;i<args.size();i++){
				int argID = args[i]->bvID;
				Weight arg_over = over_approx[argID];
				Weight arg_under = under_approx[argID];

				if(isMin){
					if(arg_under<minmax_arg_under){
						minmax_arg_under=arg_under;
					}
					if(arg_over<minmax_arg_over){
						minmax_arg_over=arg_over;
					}
				}else{
					if(arg_under>minmax_arg_under){
						minmax_arg_under=arg_under;
					}
					if(arg_over>minmax_arg_over){
						minmax_arg_over=arg_over;
					}
				}


				if(isMin ?(arg_over<overApprox || arg_under<underApprox) :( arg_over>overApprox || arg_under>underApprox)){
					addAlteredBV(argID);
				}
			}

			if(isMin){
				if(minmax_arg_over < underApprox){
					//this is a conflict
					buildReason(conflict);
					return false;
				}
				if(minmax_arg_under > overApprox){
					//this is a conflict
					buildReason(conflict);
					return false;
				}
			}else{

				if(minmax_arg_over < underApprox){
					//this is a conflict
					buildReason(conflict);
					return false;
				}
				if(minmax_arg_under > overApprox){
					//this is a conflict
					buildReason(conflict);
					return false;
				}

			}
			return true;

		}
		void updateApprox(Var ignore_bv, Weight & under_new, Weight & over_new, Cause & under_cause_new, Cause & over_cause_new) override{
			importTheory(theory);
			bool isMin = min;

			if (!isMin){
				using std::max;
				Weight highest_over=over_approx[args[0]->bvID];
				Weight highest_under=under_approx[args[0]->bvID];
				for(int i = 0;i<args.size();i++){
					int argID = args[i]->bvID;
					Weight under = under_approx[argID];
					Weight over = over_approx[argID];
					highest_under = max(under,highest_under);
					highest_over = max(over,highest_over);
				}
				if (highest_under>under_new){
					under_new=highest_under;
					under_cause_new.clear();
					under_cause_new.type =getType();
					under_cause_new.index=getID();
				}
				if (highest_over<over_new){
					over_new=highest_over;
					over_cause_new.clear();
					over_cause_new.type =getType();
					over_cause_new.index=getID();
				}
			}else{
				using std::min;
				Weight lowest_over=over_approx[args[0]->bvID];
				Weight lowest_under=under_approx[args[0]->bvID];
				for(int i = 0;i<args.size();i++){
					int argID = args[i]->bvID;
					Weight under = under_approx[argID];
					Weight over = over_approx[argID];
					lowest_under = min(under,lowest_under);
					lowest_over = min(over,lowest_over);
				}
				if (lowest_under>under_new){
					under_new=lowest_under;
					under_cause_new.clear();
					under_cause_new.type =getType();
					under_cause_new.index=getID();
				}
				if (lowest_over<over_new){
					over_new=lowest_over;
					over_cause_new.clear();
					over_cause_new.type =getType();
					over_cause_new.index=getID();
				}
			}



		}

		void buildReason( vec<Lit> & conflict){
			importTheory(theory);
			theory.dbg_no_pending_analyses();
			assert(theory.eq_bitvectors[bvID]==bvID);
			//rewind_trail_pos(trail.size()-1);
			theory.stats_build_addition_reason++;
			Weight  overApprox = over_approx[bvID];
			Weight  underApprox = under_approx[bvID];
			//assert(theory.checkApproxUpToDate(bvID));

			bool isMin = min;

			Weight minmax_arg_under = under_approx[bvID];
			Weight minmax_arg_over = over_approx[bvID];
			int minArg=0;
			int maxArg=0;
			for(int i = 0;i<args.size();i++){
				int argID = args[i]->bvID;
				Weight arg_over = over_approx[argID];
				Weight arg_under = under_approx[argID];

				if(isMin){
					if(arg_under<minmax_arg_under){
						minmax_arg_under=arg_under;
					}
					if(arg_over<minmax_arg_over){
						minmax_arg_over=arg_over;
						minArg=i;
					}
				}else{
					if(arg_under>minmax_arg_under){
						minmax_arg_under=arg_under;
						maxArg=i;
					}
					if(arg_over>minmax_arg_over){
						minmax_arg_over=arg_over;

					}
				}
			}
			bool foundConflict=false;
			if(isMin){
				if(minmax_arg_over < underApprox){
					foundConflict=true;
					//this is a conflict
					analyzeValueReason(Comparison::gt, bvID,minmax_arg_over,conflict);
					analyzeValueReason(Comparison::leq, minArg,minmax_arg_over,conflict);
				}else if(minmax_arg_under > overApprox){
					//this is a conflict
					foundConflict=true;
					analyzeValueReason(Comparison::lt, bvID,minmax_arg_under,conflict);
					for(int i = 0;i<args.size();i++){
						int argID = args[i]->bvID;
						analyzeValueReason(Comparison::geq, argID,minmax_arg_under,conflict);
					}
				}
			}else{
				if(minmax_arg_over < underApprox){
					//this is a conflict
					foundConflict=true;
					analyzeValueReason(Comparison::gt, bvID,minmax_arg_over,conflict);
					for(int i = 0;i<args.size();i++){
						int argID = args[i]->bvID;
						analyzeValueReason(Comparison::leq, argID,minmax_arg_over,conflict);
					}

				}else if(minmax_arg_under > overApprox){
					//this is a conflict
					foundConflict=true;
					analyzeValueReason(Comparison::lt, bvID,minmax_arg_under,conflict);
					analyzeValueReason(Comparison::geq, maxArg,minmax_arg_under,conflict);
				}
			}
			assert(foundConflict);
			analyze(conflict);
		}
		void analyzeReason(bool compare_over,Comparison op, Weight  to,  vec<Lit> & conflict){

			importTheory(theory);
			Weight  overApprox = over_approx[bvID];
			Weight  underApprox = under_approx[bvID];


			bool isMin =min;


			Weight minmax_arg_under = under_approx[args[0]->bvID];
			Weight minmax_arg_over = over_approx[args[0]->bvID];
			int minArg=args[0]->bvID;
			int maxArg=args[0]->bvID;
			for(int i = 0;i<args.size();i++){
				int argID = args[i]->bvID;
				Weight arg_over = over_approx[argID];
				Weight arg_under = under_approx[argID];

				if(isMin){
					if(arg_under<minmax_arg_under){
						minmax_arg_under=arg_under;
					}
					if(arg_over<minmax_arg_over){
						minmax_arg_over=arg_over;
						minArg=argID;
					}
				}else{
					if(arg_under>minmax_arg_under){
						minmax_arg_under=arg_under;
						maxArg=argID;
					}
					if(arg_over>minmax_arg_over){
						minmax_arg_over=arg_over;

					}
				}
			}

			if(isMin){
				if(compare_over){

					//analyzeValueReason(Comparison::lt, bvID,minmax_arg_under,conflict);

					//addAnalysis(Comparison::lt,bvID,minmax_arg_under,conflict);
					assert(applyOp(op,minArg,to));
					addAnalysis(op,minArg,to);
				}else{
					for(int i = 0;i<args.size();i++){
						int argID = args[i]->bvID;
						//analyzeValueReason(Comparison::geq, argID,minmax_arg_under,conflict);
						assert(applyOp(op,argID,to));
						addAnalysis(op, argID,to);
					}

					//addAnalysis(Comparison::gt, bvID,minmax_arg_over,conflict);

				}
			}else{
				if(compare_over){


					//addAnalysis(Comparison::lt, bvID,minmax_arg_under,conflict);
					assert(minmax_arg_over<=to);
					for(int i = 0;i<args.size();i++){
						int argID = args[i]->bvID;
						assert(applyOp(op,argID,to));
						addAnalysis(op, argID,to);
					}
				}else{

					assert(applyOp(op,maxArg,to));
					//addAnalysis(Comparison::gt, bvID,minmax_arg_over,conflict);
					addAnalysis(op, maxArg,to);

				}
			}

		}
		bool checkApproxUpToDate(Weight & under,Weight&over)override{
			importTheory(theory);
			bool isMin = min;

			if (!isMin){
				using std::max;
				Weight highest_over=over_approx[args[0]->bvID];
				Weight highest_under=under_approx[args[0]->bvID];
				for(int i = 0;i<args.size();i++){
					int argID = args[i]->bvID;
					Weight under2 = under_approx[argID];
					Weight over2 = over_approx[argID];
					highest_under = max(under2,highest_under);
					highest_over = max(over2,highest_over);
				}
				if (highest_under>under){
					under=highest_under;

				}
				if (highest_over<over){
					over=highest_over;

				}
			}else{
				using std::min;
				Weight lowest_over=over_approx[args[0]->bvID];
				Weight lowest_under=under_approx[args[0]->bvID];
				for(int i = 0;i<args.size();i++){
					int argID = args[i]->bvID;
					Weight under2 = under_approx[argID];
					Weight over2 = over_approx[argID];
					lowest_under = min(under2,lowest_under);
					lowest_over = min(over2,lowest_over);
				}
				if (lowest_under>under){
					under=lowest_under;

				}
				if (lowest_over<over){
					over=lowest_over;

				}
			}
			return true;
		}

		bool checkSolved()override{
			importTheory(theory);
			bool isMin = min;


			Weight min_arg=over_approx[args[0]->bvID];
			Weight max_arg=under_approx[args[0]->bvID];

			for(int i = 0;i<args.size();i++){
				int argID = args[i]->bvID;
				if (over_approx[argID]<min_arg){
					min_arg=over_approx[argID];
				}
				if (under_approx[argID]>max_arg){
					max_arg=under_approx[argID];
				}
			}
			if(isMin){
				if(min_arg!=over_approx[bvID]){
					return false;
				}
			}else{
				if (max_arg!=under_approx[bvID]){
					return false;
				}
			}
			return true;
		}
	};

	class MinMaxArg:public Operation{
			public:
		using Operation::getID;
		using Operation::theory;
	public:
		bool min;
		//int resultID;
		MinMaxData * resultOp=nullptr;

		int bvID=-1;

		void move( int bvID) override{
			this->bvID=bvID;
		}

		MinMaxArg(BVTheorySolver & theory, int bvID,MinMaxData * resultOp):Operation(theory),bvID(bvID),resultOp(resultOp){
			this->min = resultOp->min;
		}
		int getBV()override{
			return bvID;
		}
		OperationType getType()const override{
			return OperationType::cause_is_minmax_argument;
		}

		bool propagate( bool & changed_outer,vec<Lit> & conflict) override{
			importTheory(theory);
			Weight & underApprox = under_approx[bvID];
			Weight & overApprox = over_approx[bvID];

			bool isMin = min;
			int resultID= resultOp->bvID;

			Weight res_over = over_approx[resultID];
			Weight res_under = under_approx[resultID];
			if(isMin ? (overApprox<res_under) : (underApprox>res_over) ){
				//conflict
				double startconftime = rtime(2);
				//propagationtime += startconftime - startproptime;
				theory.stats_num_conflicts++;theory.stats_addition_conflicts++;
				if(opt_verb>1){
					printf("bv minmax conflict %ld\n", theory.stats_num_conflicts);
				}
				buildReason(conflict);

				theory.stats_conflict_time+=rtime(2)-startconftime;
				return false;
			}
			if(overApprox<res_over || underApprox>res_under){

				addAlteredBV(resultID);
			}

			return true;

		}
			void updateApprox(Var ignore_bv, Weight & under_new, Weight & over_new, Cause & under_cause_new, Cause & over_cause_new) override{
				importTheory(theory);
				bool isMin = min;
				int resultID= resultOp->bvID;
				if (!isMin){
					using std::max;
					Weight highest_over = over_approx[resultID];

					if (highest_over<over_new){
						over_new=highest_over;
						over_cause_new.clear();
						over_cause_new.type =getType();
						over_cause_new.index=getID();
					}
				}else{
					using std::min;
					Weight lowest_under = under_approx[resultID];

					if (lowest_under>under_new){
						under_new=lowest_under;
						under_cause_new.clear();
						under_cause_new.type =getType();
						under_cause_new.index=getID();
					}
				}
			}

			void buildReason(vec<Lit> & conflict){
				importTheory(theory);
				theory.dbg_no_pending_analyses();
				bool isMin = min;
				int resultID= resultOp->bvID;
				assert(theory.eq_bitvectors[bvID]==bvID);
				//rewind_trail_pos(trail.size()-1);
				theory.stats_build_addition_reason++;
				Weight  overApprox = over_approx[bvID];
				Weight  underApprox = under_approx[bvID];
				//assert(theory.checkApproxUpToDate(bvID));

				int argID = bvID;
				Weight arg_under = under_approx[argID];
				Weight arg_over = over_approx[argID];

				if(isMin){
					if(arg_over>=underApprox)
						throw std::runtime_error("Error in min reason");

					theory.buildComparisonReasonBV(Comparison::lt,argID,bvID,conflict);
				}else{
					if(arg_under<=overApprox)
						throw std::runtime_error("Error in max reason");

					theory.buildComparisonReasonBV(Comparison::gt,argID,bvID,conflict);
				}
			}
			void analyzeReason(bool compare_over,Comparison op, Weight  to,  vec<Lit> & conflict){
				importTheory(theory);
				bool isMin = min;
				int resultID= resultOp->bvID;
				if (!isMin){
					using std::max;
					Weight highest_over = over_approx[resultID];
					assert(compare_over);
					assert(applyOp(op,resultID,to));
					addAnalysis(op, resultID,to);
				}else{
					using std::min;
					Weight lowest_under = under_approx[resultID];
					assert(!compare_over);
					assert(applyOp(op,resultID,to));
					addAnalysis(op, resultID,to);
				}
			}

			bool checkApproxUpToDate(Weight & under,Weight&over)override{
				importTheory(theory);
				bool isMin = min;
				int resultID= resultOp->bvID;
				if (!isMin){
					using std::max;
					Weight highest_over = over_approx[resultID];
					if (highest_over<over){
						over=highest_over;
					}
				}else{
					using std::min;
					Weight lowest_under = under_approx[resultID];
					if (lowest_under>under){
						under=lowest_under;
					}
				}
				return true;
			}

			bool checkSolved()override{
				bool isMin = min;
				int resultID= resultOp->bvID;

				return true;
			}
		};
	class Invert:public Operation{
				public:
			using Operation::getID;
			using Operation::theory;
		public:

			 int iter = 0;
			int bvID=-1;
			Invert * argOp;
			void move( int bvID) override{
				this->bvID=bvID;
			}

			Invert(BVTheorySolver & theory, int bvID):Operation(theory),bvID(bvID){

			}
			void setArg(Invert * arg){
				this->argOp = arg;
			}
			int getBV()override{
				return bvID;
			}
			OperationType getType()const override{
				return OperationType::cause_is_invert;
			}

			bool propagate( bool & changed_outer,vec<Lit> & conflict) override{
				importTheory(theory);
				Weight & underApprox = under_approx[bvID];
				Weight & overApprox = over_approx[bvID];

				int width =  theory.bitvectors[bvID].size();
				Weight max_val = ((1L)<<width)-1;
				Weight & otherUnder = under_approx[ argOp->bvID];
				Weight & otherOver = over_approx[ argOp->bvID];


				Weight thisNewOver = max_val -otherUnder;
				Weight thisNewUnder = max_val -otherOver;


				Weight otherNewOver = max_val -underApprox;
				Weight otherNewUnder = max_val -overApprox;

				if(thisNewOver <underApprox ){
					//conflict
					double startconftime = rtime(2);
					//propagationtime += startconftime - startproptime;
					theory.stats_num_conflicts++;theory.stats_addition_conflicts++;
					if(opt_verb>1){
						printf("bv inversion conflict %ld\n", theory.stats_num_conflicts);
					}

					//the reason is that either this bv needs to be less than underApprox
					theory.buildComparisonReason(Comparison::geq,bvID,underApprox,conflict);
					//or argID needs to be greater
					theory.buildComparisonReason(Comparison::leq,argOp->getBV(),otherOver,conflict);


					theory.stats_conflict_time+=rtime(2)-startconftime;
					return false;
				}else if(thisNewUnder >overApprox ){
					//conflict
					double startconftime = rtime(2);
					//propagationtime += startconftime - startproptime;
					theory.stats_num_conflicts++;theory.stats_addition_conflicts++;
					if(opt_verb>1){
						printf("bv inversion conflict %ld\n", theory.stats_num_conflicts);
					}
					//the reason is that either this bv needs to be less than underApprox
					theory.buildComparisonReason(Comparison::leq,bvID,overApprox,conflict);
					//or argID needs to be greater
					theory.buildComparisonReason(Comparison::geq,argOp->getBV(),otherUnder,conflict);


					theory.stats_conflict_time+=rtime(2)-startconftime;
					return false;
				}else if (otherNewUnder>otherUnder){
					addAlteredBV(argOp->getBV());
				}else if (otherNewUnder<otherOver){
					addAlteredBV(argOp->getBV());
				}


				return true;

			}
				void updateApprox(Var ignore_bv, Weight & under_new, Weight & over_new, Cause & under_cause_new, Cause & over_cause_new) override{
					importTheory(theory);

					++iter;
					//printf("inv %d iter %d\n",bvID,iter);
					if(iter>=176){
						int a=1;
					}
					int width =  theory.bitvectors[bvID].size();
					Weight max_val = ((1L)<<width)-1;
					Weight & otherUnder = under_approx[ argOp->bvID];
					Weight & otherOver = over_approx[ argOp->bvID];

					Weight thisOver = max_val -otherUnder;
					Weight thisUnder = max_val -otherOver;
					assert(thisOver>=0);
					assert(thisUnder>=0);
					assert(thisOver<=max_val);
					assert(thisUnder<=max_val);

					if (thisOver<over_new){
						over_new=thisOver;
						over_cause_new.clear();
						over_cause_new.type =getType();
						over_cause_new.index=getID();
					}
					if (thisUnder>under_new){
						under_new=thisUnder;
						under_cause_new.clear();
						under_cause_new.type =getType();
						under_cause_new.index=getID();
					}
				}


				void analyzeReason(bool compare_over,Comparison op, Weight  to,  vec<Lit> & conflict){
					importTheory(theory);
					int width =  theory.bitvectors[bvID].size();
					Weight max_val = ((1L)<<width)-1;
					Weight val = max_val-to;
					theory.analyzeValueReason(~op,argOp->getBV(),val,conflict);
				}

				bool checkApproxUpToDate(Weight & under,Weight&over)override{
					importTheory(theory);
					int width =  theory.bitvectors[bvID].size();
					Weight max_val = ((1L)<<width)-1;
					Weight & otherUnder = under_approx[ argOp->bvID];
					Weight & otherOver = over_approx[ argOp->bvID];

					Weight thisOver = max_val -otherUnder;
					Weight thisUnder = max_val -otherOver;
					assert(thisOver>=0);
					assert(thisUnder>=0);
					assert(thisOver<=max_val);
					assert(thisUnder<=max_val);

					if (thisOver<over){
						over=thisOver;
					}
					if (thisUnder>over){
						over=thisUnder;
					}
					return true;
				}

				bool checkSolved()override{
					importTheory(theory);
					int width =  theory.bitvectors[bvID].size();
					Weight max_val = ((1L)<<width)-1;
					Weight & otherUnder = under_approx[ argOp->bvID];
					Weight & otherOver = over_approx[ argOp->bvID];

					Weight thisOver = max_val -otherUnder;
					Weight thisUnder = max_val -otherOver;

					Weight & underApprox = under_approx[bvID];
					Weight & overApprox = over_approx[bvID];


					if(thisOver<underApprox){
						return false;
					}else if (thisUnder>overApprox){
						return false;
					}
					return true;
				}
			};
	vec<vec<int>> operation_ids;

	vec<Operation*> operations;
	bool hasOperation(Lit l){
		return vars[var(l)].operationID>-1;
	}
	Operation & getOperation(Lit l){
		int opID = vars[var(l)].operationID;
		assert(opID>-1);
		return getOperation(opID);
	}
	Operation & getOperation(int opID){
		return *operations[opID];
	}

	void addOperation(int bvID,Operation* op){
		/*assert(op->getID()==operations.size());
		operations.push(op);*/
		operation_ids[bvID].push(op->getID());

	}






	class BitVector{
		BVTheorySolver * outer;
		int id;
	public:
		BitVector():outer(nullptr),id(0){

			}
		BitVector(BVTheorySolver & outer, int id):outer(&outer),id(id){

		}
		BitVector(const BitVector & copy):outer(copy.outer),id(copy.id){

			}
		BitVector(const BitVector && move):outer(move.outer),id(move.id){

			}
		BitVector& operator=(BitVector   other) {
			outer= other.outer;
			id=other.id;
			return *this;
		}

		Lit& operator[] (const int i){
			return outer->getBits(id)[i];
		}

		int getID()const{
			return id;
		}
		Weight & getUnder(bool level0=false){
			return outer->getUnderApprox(id,level0);
		}
		Weight & getOver(bool level0=false){
			return outer->getOverApprox(id,level0);
		}

		vec<Lit> & getBits()const{
			return outer->getBits(id);
		}

		bool isConst()const{
			return outer->isConst(id);
		}

		int width()const{
			return getBits().size();
		}
	};

	//typedef void (*CallBack)(int) ;

	double rnd_seed;

private:
	TheorySolver * S;

	bool comp(Comparison op, Weight o1, Weight o2){
		switch (op){
			case Comparison::lt:
				return o1<o2;
			case Comparison::leq:
				return  o1<=o2;
			case Comparison::gt:
				return  o1>o2;
			default:
				return o1>=o2;
		}
	}

	BVMap * bvRemap=nullptr;
public:

	bool first_propagation=true;
	long n_bits =0;
	long n_consts = 0;
	long n_starting_consts=0;
	long n_additions=0;
	long n_popcounts=0;
	vec<lbool> assigns;
	CRef comparisonprop_marker;
	CRef comparisonbv_prop_marker;
	CRef comparisonbv_1starg_marker;
	CRef comparisonbv_2ndarg_marker;
	CRef conditionthen_prop_marker;
	CRef conditionelse_prop_marker;
	CRef conditionarg_prop_marker;
	CRef bvprop_marker;
	CRef popcount_marker;
	Lit const_true=lit_Undef;
	vec<const char*> symbols;

	vec<Assignment> trail;
	vec<int> trail_lim;

	vec<Cause> under_causes;
	vec<Cause> over_causes;
	//Bitvectors are unsigned, and have least significant bit at index 0
	vec<vec<Lit> > bitvectors;





	vec<vec<int> > cause_set;//for each bv, this is the set of all bitvectors that have a greater id and might need to have their approx updated when this bv's approx changes.
	vec<vec<int> > compares; //for each bitvector, comparisons are all to unique values, and in ascending order of compareTo.
	vec<vec<int>> bvcompares;
	vec<int> eq_bitvectors;//if a bv has been proven to be equivalent to another, lower index bv, put the lowest such index here.
	vec<bool> bv_needs_propagation;
	vec<bool> comparison_needs_repropagation;
	vec<int> repropagate_comparisons;




	struct ToAnalyze{
		int bvID;
		int next_analysis;
		Weight value;
		void clear(){
			bvID=-1;
			next_analysis=-1;
		}
	};
	vec<ToAnalyze> analyses;
	vec<int> pending_over_analyses;
	vec<int> pending_under_analyses;
	int n_pending_analyses = 0;

/*	vec<vec<Addition>> additions;
	vec<vec<AdditionArg>> addition_arguments;*/

/*
	struct PopCountData{
		int resultID;
		vec<Lit> args;
		Weight under=0;
		Weight over=0;
	};
	vec<PopCountData> pop_counts;//improve this
	vec<vec<int>> pop_count_ids;
*/


/*
	vec<MinMaxData> minmaxs;
	vec<vec<int>> minmax_ops;//for now, allowing at most one minmax operation per bvID
	vec<vec<int>> minmax_args;*/

	vec<Weight> under_approx0;//under approx at level 0
	vec<Weight> over_approx0;//over approx at level 0
	vec<Weight> under_approx;
	vec<Weight> over_approx;
	vec<int> theoryIds;
	vec<int> altered_bvs;
	vec<bool> alteredBV;

	vec<int> backtrack_altered;

	vec<bool> bvconst;

	int analysis_trail_pos = -1;

	vec<bool> in_backtrack_queue;

	vec<int> backtrack_queue;
	vec<BVTheory*> theories;
	vec<BVTheory*> actual_theories;
public:


	//vec<int> marker_map;

	bool requiresPropagation = true;

	vec<char> seen;
	vec<int> to_visit;
	vec<Lit> tmp_clause;
	//Data about local theory variables, and how they connect to the sat solver's variables
	struct VarData {

		int occursPositive :1;
		int occursNegative :1;
		//int isOperation:1;
		//int bvID:29;
		int operationID:30;
		Var solverVar;
	};

	vec<VarData> vars;
	int theory_index = 0;
public:

	double reachtime = 0;
	double unreachtime = 0;
	double pathtime = 0;
	double propagationtime = 0;
	double stats_conflict_time = 0;
	long stats_propagations = 0;
	long stats_bv_propagations =0;
	long stats_bv_skipped_propagations=0;
	long stats_num_conflicts = 0;
	long stats_bit_conflicts = 0;
	long stats_addition_conflicts = 0;
	long stats_compare_conflicts = 0;
	long stats_bv_compare_conflicts = 0;
	long stats_decisions = 0;
	long stats_num_reasons = 0;
	long stats_build_value_reason=0;
	long stats_build_value_bv_reason=0;
	long stats_build_condition_reason=0;
	long stats_build_condition_arg_reason=0;
	long stats_build_addition_reason=0;
	long stats_build_addition_arg_reason =0;
	double stats_update_time=0;
	double reachupdatetime = 0;
	double unreachupdatetime = 0;
	double stats_initial_propagation_time = 0;
	double stats_decision_time = 0;
	double stats_reason_initial_time = 0;
	double stats_reason_time = 0;
	long num_learnt_paths = 0;
	long learnt_path_clause_length = 0;
	long num_learnt_cuts = 0;
	long learnt_cut_clause_length = 0;
	long stats_pure_skipped = 0;
	long stats_mc_calls = 0;
	long propagations =-1;
	long stats_propagations_skipped = 0;
	long statis_bv_updates = 0;

	BVTheorySolver(TheorySolver * S ) :
			S(S){
		rnd_seed = drand(S->getRandomSeed());
		S->addTheory(this);
		S->setBVTheory(this);
		comparisonprop_marker = S->newReasonMarker(this);
		comparisonbv_prop_marker = S->newReasonMarker(this);
		comparisonbv_1starg_marker = S->newReasonMarker(this);
		comparisonbv_2ndarg_marker=S->newReasonMarker(this);
		conditionthen_prop_marker = S->newReasonMarker(this);
		conditionelse_prop_marker = S->newReasonMarker(this);
		conditionarg_prop_marker= S->newReasonMarker(this);
		bvprop_marker = S->newReasonMarker(this);
		popcount_marker = S->newReasonMarker(this);

	}
	~BVTheorySolver(){

	}
	TheorySolver * getSolver(){
		return S;
	}

	bool hasTheory(int bvID){
		assert(bvID>=0);
		return theoryIds[bvID]>=0;
	}
	BVTheory * getTheory(int bvID){
		assert(hasTheory(bvID));
		return theories[theoryIds[bvID]];
	}


	void addTheory(BVTheory* theory){
		theories.growTo(theory->getTheoryIndexBV()+1);
		actual_theories.push(theory);
		theories[theory->getTheoryIndexBV()]=theory;
	}

	void printStats(int detailLevel) {
		printf("BV Theory %d stats:\n", this->getTheoryIndex());

		printf("%d bitvectors, %ld bits, %d comparisons (bvcomparisons %d), %ld additions\n", bitvectors.size(),n_bits,compares.size()+bvcompares.size(),bvcompares.size(), n_additions				 );
		printf("constant bitvectors (at start, end of deduction): %ld, %ld\n",n_starting_consts ,n_consts);


		printf("Propagations: %ld (%f s, avg: %f s, %ld skipped), bv updates: %ld (%f s), bv propagations %ld (%ld skipped)\n", stats_propagations, propagationtime,
				(propagationtime) / ((double) stats_propagations + 1), stats_propagations_skipped,statis_bv_updates,stats_update_time,stats_bv_propagations,stats_bv_skipped_propagations);
		printf("Decisions: %ld (%f s, avg: %f s)\n", stats_decisions, stats_decision_time,
				(stats_decision_time) / ((double) stats_decisions + 1));
		printf("Conflicts: %ld (bits: %ld, additions: %ld, comparisons: %ld, bv comparisons: %ld), %f seconds\n", stats_num_conflicts,stats_bit_conflicts,stats_addition_conflicts,stats_compare_conflicts,stats_bv_compare_conflicts, stats_conflict_time);
		printf("Reasons: %ld (%f s, avg: %f s)\n", stats_num_reasons, stats_reason_time,
				(stats_reason_time) / ((double) stats_num_reasons + 1));
		printf("Build: value reason %ld, bv value reason %ld, addition reason %ld\n", stats_build_value_reason,stats_build_value_bv_reason,stats_build_addition_reason);

		fflush(stdout);
	}
	
	void writeTheoryWitness(std::ostream& write_to) {
		for(int bvID = 0;bvID<bitvectors.size();bvID++){
				vec<Lit> & bv = bitvectors[bvID];
				updateApproximations(bvID);
				Weight under=under_approx[bvID];
				write_to<<"bv"<< this->unmapBV(bvID) << " = " << under << "\n";

		}
	}
	
	inline int getTheoryIndex() {
		return theory_index;
	}
	inline void setTheoryIndex(int id) {
		theory_index = id;
	}

	bool hasEquivalentBV(int bvID){
		return eq_bitvectors[bvID]!=bvID;
	}

	int getEquivalentBV(int bvID){
		return eq_bitvectors[bvID];
	}

	BitVector duplicateBV(int bvID){
		BitVector bv = newBitvector(-1, bitvectors[bvID].size() ,-1, bvID);
		makeEquivalent(bv.getID(),bvID);
		return bv;
	}

	void makeEquivalent(int bvID1, int bvID2){
		if (bvID1==bvID2)
			return;
		while(eq_bitvectors[bvID1]!=bvID1){
			bvID1=eq_bitvectors[bvID1];
		}
		while(eq_bitvectors[bvID2]!=bvID2){
			bvID2=eq_bitvectors[bvID2];
		}
		if (bvID1<bvID2)
			makeEquivalent(bvID2,bvID1);


		if(bitvectors[bvID1].size() && bitvectors[bvID2].size()){
			assert(bitvectors[bvID1].size()==bitvectors[bvID2].size());
			for(int i = 0;i<bitvectors[bvID1].size();i++){
				makeEqual(bitvectors[bvID1][i],bitvectors[bvID2][i]);
			}
		}
		//move the conditionals of bvID1 into bvID2.
		for (int opID:operation_ids[bvID1]){
			Operation & op = getOperation(opID);
			op.move(bvID2);
		}

		for(int i = compares[bvID1].size()-1;i>=0;i--){
			int cID = compares[bvID1][i];
			Operation & op = getOperation(cID);
			op.move(bvID2);
			compares[bvID2].push(cID);
		}

		for(int i = bvcompares[bvID1].size()-1;i>=0;i--){
			int cID = bvcompares[bvID1][i];
			Operation & op = getOperation(cID);
			op.move(bvID2);
			bvcompares[bvID2].push(cID);
		}

		//cause_set[bvID1].push(bvID2);
		eq_bitvectors[bvID1]=bvID2;
		cause_set[bvID2].push(bvID1);
		for (int bv:cause_set[bvID1]){
			if(bv!=bvID2){
				cause_set[bvID2].push(bv);
			}
		}
		cause_set[bvID1].clear();
		//bv_needs_propagation[bvID1]=true;
/*		if(!alteredBV[bvID1]){
			alteredBV[bvID1]=true;
			altered_bvs.push(bvID1);
		}*/

		bv_needs_propagation[bvID2]=true;
		if(!alteredBV[bvID2]){
			alteredBV[bvID2]=true;
			altered_bvs.push(bvID2);
		}
		requiresPropagation=true;
		//merge the bitvector's causes.
	}
	bool hasBV(int bvID)const{
		return bvID>=0 && bvID<under_approx.size() && under_approx[bvID]>-1;
	}
	BitVector getBV(int bvID){
		return BitVector(*this,bvID);
	}
	/*inline bool isBVVar(Var v)const{
		return vars[v].bvID>=0;
	}
*/
	inline bool isOperationVar(Var v) const{
		assert(v < vars.size());
		return vars[v].operationID>-1 ;// && vars[v].isOperation;
	}
	inline int getOperationID(Var v) const {
		assert(isOperationVar(v));
		return vars[v].operationID;
	}

	/*inline bool isComparisonVar(Var v) const{
		assert(v < vars.size());
		return vars[v].operationID>-1 && !vars[v].isOperation;
	}
	inline int getComparisonID(Var v) const {
		assert(isComparisonVar(v));
		return vars[v].operationID;
	}



	inline ComparisonID& getComparison(int comparisonID){
		return comparisons[comparisonID];
	}
	inline int getbvID(Var v) const{
		return vars[v].bvID;
	}*/

	void makeEqual(Lit l1, Lit l2) {
		Lit o1 = toSolver(l1);
		Lit o2 = toSolver(l2);
		tmp_clause.clear();
		tmp_clause.push(~o1);
		tmp_clause.push(o2);
		S->addClauseSafely(tmp_clause);
		tmp_clause.clear();
		tmp_clause.push(o1);
		tmp_clause.push(~o2);
		S->addClauseSafely(tmp_clause);
	}
	void makeEqualInSolver(Lit o1, Lit o2) {
		tmp_clause.clear();
		tmp_clause.push(~o1);
		tmp_clause.push(o2);
		S->addClauseSafely(tmp_clause);
		tmp_clause.clear();
		tmp_clause.push(o1);
		tmp_clause.push(~o2);
		S->addClauseSafely(tmp_clause);
	}
	void addClause(Lit l1) {
		Lit o1 = toSolver(l1);
		tmp_clause.clear();
		tmp_clause.push(o1);
		S->addClauseSafely(tmp_clause);
	}
	void addClause(Lit l1, Lit l2) {
		Lit o1 = toSolver(l1);
		Lit o2 = toSolver(l2);
		tmp_clause.clear();
		tmp_clause.push(o1);
		tmp_clause.push(o2);

		S->addClauseSafely(tmp_clause);
	}
	void addClause(Lit l1, Lit l2, Lit l3) {
		Lit o1 = toSolver(l1);
		Lit o2 = toSolver(l2);
		Lit o3 = toSolver(l3);
		tmp_clause.clear();
		tmp_clause.push(o1);
		tmp_clause.push(o2);
		tmp_clause.push(o3);
		S->addClauseSafely(tmp_clause);
	}
	void addClause(vec<Lit> & c) {
		tmp_clause.clear();
		c.copyTo(tmp_clause);
		toSolver(tmp_clause);
		S->addClauseSafely(tmp_clause);
	}
	void addClauseSafely(vec<Lit> & c) {
		tmp_clause.clear();
		c.copyTo(tmp_clause);
		toSolver(tmp_clause);
		
		S->addClauseSafely(tmp_clause);
	}
	

	Var newVar(Var solverVar=var_Undef, int operationID=-1, bool connectToTheory = true, bool decidable=true) {
		if(solverVar==var_Undef){
			solverVar = S->newVar();
		}
		Var v = vars.size();
		if (connectToTheory) {
			S->newTheoryVar(solverVar, getTheoryIndex(), v);

		}else{
			while (S->nVars() <= solverVar)
				S->newVar();
		}
		S->setDecisionVar(solverVar, decidable);

		vars.push();
		vars[v].occursPositive=false;
		vars[v].occursNegative=false;
		/*vars[v].isOperation=isOperation;
		vars[v].bvID=bvID;*/
		vars[v].operationID=operationID;
		vars[v].solverVar = solverVar;

		assigns.push(l_Undef);
		if (connectToTheory) {
			assert(toSolver(v) == solverVar);
		}
		return v;
	}
	inline int level(Var v) {
		return S->level(toSolver(v));
	}
	inline int decisionLevel() {
		return trail_lim.size(); //S->decisionLevel();
	}
	inline int nVars() const {
		return vars.size(); //S->nVars();
	}
	inline Var toSolver(Var v) {
		//return v;
		assert(v < vars.size());
		//assert(S->hasTheory(vars[v].solverVar));
		//assert(S->getTheoryVar(vars[v].solverVar)==v);
		return vars[v].solverVar;
	}
	
	inline Lit toSolver(Lit l) {
		//assert(S->hasTheory(vars[var(l)].solverVar));
		//assert(S->getTheoryVar(vars[var(l)].solverVar)==var(l));
		return mkLit(vars[var(l)].solverVar, sign(l));
	}
	
	void toSolver(vec<Lit> & c) {
		for (int i = 0; i < c.size(); i++) {
			c[i] = toSolver(c[i]);
		}
	}
	
	inline lbool value(Var v) {
		if (assigns[v] != l_Undef)
			assert(S->value(toSolver(v)) == assigns[v]);
		
		return assigns[v]; //S->value(toSolver(v));
	}
	inline lbool value(Lit l) {
		if (assigns[var(l)] != l_Undef) {
			assert(S->value(toSolver(l)) == (assigns[var(l)] ^ sign(l)));
		}
		return assigns[var(l)] ^ sign(l);
	}
	inline lbool dbg_value(Var v) {
		return S->value(toSolver(v));
	}
	inline lbool dbg_value(Lit l) {
		return S->value(toSolver(l));
	}

	void enqueueEager(Lit l,int bvID, Weight prev_under, Weight prev_over, Cause prev_under_cause, Cause prev_over_cause, CRef reason_marker){
		enqueue(l,reason_marker);
		//if newly created lits are enqueued, then they must provide a reason eagerly (so that complex propagation dependencies can be handled correctly by the solver, instead of having to be sorted out in the theories).
		/*static int iter = 0;
		++iter;
		vec<Lit>  reason;
		buildReason(l,reason,reason_marker);

		Lit sl = toSolver(l);

		S->addClauseSafely(reason);

		//the lit must have been propagated by this clause (or, alternatively, the solver might now have a conflict).
		if(S->value(sl)==l_True)
			enqueueTheory(l);*/
	}

	inline bool enqueue(Lit l,CRef reason) {
		assert(assigns[var(l)]==l_Undef);
		if(l.x==55){
			int a=1;
		}
		Lit sl = toSolver(l);
		if(sl.x ==83){
			int a=1;
		}
#ifndef NDEBUG

	/*printf("learnt ");
	for (int i = 0;i<nVars();i++){
			if(value(i)!=l_Undef){
				if(i!=var(l)){
					Lit l = mkLit(i,value(i)==l_True);
					assert(value(l)==l_False);
					printf(" %d",dimacs(toSolver(l)));
				}
			}
		}
	printf(" %d", dimacs(sl));
	printf(" 0\n");*/
#endif

		if (S->enqueue(sl, reason)) {
			return true;
		} else {
			return false;
		}
	}

	void rewind_trail_pos(int trail_pos){
		if(analysis_trail_pos==trail_pos)
			return;
		if (analysis_trail_pos>=trail.size()){
			analysis_trail_pos=trail.size()-1;
		}
		if(trail.size()==0){
			analysis_trail_pos=-1;
			return;
		}

		if(trail_pos>=trail.size()){
			trail_pos=trail.size()-1;
		}
		if(analysis_trail_pos<-1){
			analysis_trail_pos=-1;
		}
		if(analysis_trail_pos<trail_pos){
			for(;analysis_trail_pos<trail_pos;analysis_trail_pos++){
				assert(analysis_trail_pos>=-1);assert(analysis_trail_pos<trail.size()-1);
				Assignment & e = trail[analysis_trail_pos+1];

				if (e.isBoundAssignment()){
					int bvID = e.bvID;
					if(bvID==6){
						int a =1;
					}
					under_approx[bvID]=e.new_under;
					over_approx[bvID]=e.new_over;
					under_causes[bvID]=e.new_under_cause;
					over_causes[bvID]=e.new_over_cause;
					if(hasTheory(bvID))
						getTheory(bvID)->rewindBV(bvID);
				}else{
					Var x = e.var;
					Lit p = mkLit(x,!e.assign);
					assert(value(x)==l_Undef);
					assigns[x] = sign(p) ? l_False : l_True;
					if(e.isOperation){
						int opID = e.bvID;
						if(opID>-1){
							getOperation(opID).enqueue(p,false);
						}
					}
				}
			}

		}else if(analysis_trail_pos>trail_pos){
			for(;analysis_trail_pos>trail_pos;analysis_trail_pos--){
				assert(analysis_trail_pos>=0);assert(analysis_trail_pos<trail.size());
				Assignment & e = trail[analysis_trail_pos];

				if (e.isBoundAssignment()){
					int bvID = e.bvID;
					if(bvID==6){
						int a =1;
					}
					under_approx[bvID]=e.previous_under;
					over_approx[bvID]=e.previous_over;
					under_causes[bvID]=e.prev_under_cause;
					over_causes[bvID]=e.prev_over_cause;
					if(hasTheory(bvID))
						getTheory(bvID)->rewindBV(bvID);
				}else{
					Var x = e.var;
					if(e.isOperation){
						int opID = e.bvID;
						if(opID>-1){
							getOperation(opID).backtrack(e,true);
						}
					}
					Lit p = mkLit(x,!e.assign);
					assert(value(p)==l_True);
					assigns[x] = l_Undef;
				}
			}
		}
		assert(analysis_trail_pos>=-1);
		assert(analysis_trail_pos<=trail.size());
		assert(analysis_trail_pos==trail_pos);
	}

	int rewindUntil(int for_bvID, Comparison op, Weight compareTo){
		if(analysis_trail_pos<-1){
			analysis_trail_pos=-1;
		}
		if(analysis_trail_pos>=trail.size()){
			analysis_trail_pos=trail.size()-1;
		}
		//rewind until the previous under or over approx of bvID violates the comparison.
		for(;analysis_trail_pos>=0;analysis_trail_pos--){
			Assignment & e = trail[analysis_trail_pos];
			if (e.isBoundAssignment()){
				int bvID = e.bvID;
				if(bvID==6){
					int a =1;
				}
				if (bvID==for_bvID){
					if (op==Comparison::lt){
						assert(over_approx[for_bvID]<compareTo);
						if (e.previous_over >= compareTo){
							break;
						}
					}else if (op==Comparison::leq){
						assert(over_approx[for_bvID]<=compareTo);
						if (e.previous_over > compareTo){
							break;
						}
					}else if (op==Comparison::gt){
						assert(under_approx[for_bvID]>compareTo);
						if (e.previous_under <= compareTo){
							break;
						}
					}else if (op==Comparison::geq){
						assert(under_approx[for_bvID]>=compareTo);
						if (e.previous_under < compareTo){
							break;
						}
					}
				}

				under_approx[bvID]=e.previous_under;
				over_approx[bvID]=e.previous_over;
				under_causes[bvID]=e.prev_under_cause;
				over_causes[bvID]=e.prev_over_cause;
				if(hasTheory(bvID))
					getTheory(bvID)->rewindBV(bvID);
			}else{
				Var x = e.var;
				Lit p = mkLit(x,!e.assign);
				if(e.isOperation){
					int opID = e.bvID;
					if(opID>-1){
						getOperation(opID).backtrack(e,true);
					}
				}
				assert(value(p)==l_True);
				assigns[x] = l_Undef;
			}
		}

		assert(analysis_trail_pos>=-1);
		assert(analysis_trail_pos<=trail.size());
		return analysis_trail_pos;
	}

	int rewindUntil(Var until_assign){
		if(analysis_trail_pos<-1){
			analysis_trail_pos=-1;
		}
		if(analysis_trail_pos>=trail.size()){
			analysis_trail_pos=trail.size()-1;
		}
		//rewind until the previous under or over approx of bvID violates the comparison.
		for(;analysis_trail_pos>=0;analysis_trail_pos--){
			Assignment & e = trail[analysis_trail_pos];
			if (e.isBoundAssignment()){
				int bvID = e.bvID;
				under_approx[bvID]=e.previous_under;
				over_approx[bvID]=e.previous_over;
				under_causes[bvID]=e.prev_under_cause;
				over_causes[bvID]=e.prev_over_cause;
				if(hasTheory(bvID))
					getTheory(bvID)->rewindBV(bvID);
			}else{
				Var x = e.var;
				if(x==until_assign){
					break;
				}
				Lit p = mkLit(x,!e.assign);
				if(e.isOperation){
					int opID = e.bvID;
					if(opID>-1){
						getOperation(opID).backtrack(e,true);
					}
				}
				assigns[x] = l_Undef;
			}
		}

		assert(analysis_trail_pos>=-1);
		assert(analysis_trail_pos<=trail.size());
		return analysis_trail_pos;
	}


	void backtrackUntil(int level) {
		//it is NOT safe to remove altered bitvectors here, because if a comparison was added at a higher level, and then
		//a conflict was discovered _Before_ the comparison was processed, then the comparison may never be propagated at all if the altered_bvs are cleared here.
/*		for(int bvID:altered_bvs){
			alteredBV[bvID]=false;
		}
		altered_bvs.clear();*/

		rewind_trail_pos(trail.size());
		//need to remove and add edges in the two graphs accordingly.
		if (trail_lim.size() > level) {
			
			int stop = trail_lim[level];
			for (int i = trail.size() - 1; i >= trail_lim[level]; i--) {
				
				Assignment & e = trail[i];
				if(e.isBoundAssignment()){
					int bvID = e.bvID;

					under_approx[bvID]=e.previous_under;
					over_approx[bvID]=e.previous_over;
					under_causes[bvID]=e.prev_under_cause;
					over_causes[bvID]=e.prev_over_cause;

					if(hasTheory(bvID))
						getTheory(bvID)->backtrackBV(bvID);//only enqueue the bitvector in the subtheory _after_ it's approximation has been updated!

					for(int changedID: cause_set[bvID]){
						under_approx[changedID]=under_approx[bvID];
						over_approx[changedID]=over_approx[bvID];
						under_approx0[changedID]=under_approx0[bvID];
						over_approx0[changedID]=over_approx0[bvID];

						if(hasTheory(changedID))
							getTheory(changedID)->backtrackBV(changedID);//only enqueue the bitvector in the subtheory _after_ it's approximation has been updated!
					}
				}else{
					assert(assigns[e.var]!=l_Undef);

					int opID = e.bvID;
					if(opID>-1)
						getOperation(opID).backtrack(e,false);


					assigns[e.var] = l_Undef;
				}
				//changed = true;
				trail.pop();
			}
			//trail.shrink(trail.size() - stop);
			trail_lim.shrink(trail_lim.size() - level);
			assert(trail_lim.size() == level);
			assert(dbg_uptodate());
			if(level==0){//decisionLevel()==0 This check can fail if the levels of the theory and sat solver are out of sync
				for(int cID:repropagate_comparisons){
					assert(comparison_needs_repropagation[cID]);
					comparison_needs_repropagation[cID]=false;
					//need to handle the case where we backtracked before ever propagating this comparison's bitvector...
					Operation & c = getOperation(cID);
					int bvID = c.getBV();
					if(!alteredBV[bvID]){
						alteredBV[bvID]=true;
						altered_bvs.push(bvID);
					}
					bv_needs_propagation[bvID]=true;
					requiresPropagation=true;
					S->needsPropagation(getTheoryIndex());
				}
				repropagate_comparisons.clear();
			}
			analysis_trail_pos = trail.size()-1;

		}

	};

	bool decidableBV(Comparison op, int bvID, Weight to){
		switch(op){
			case Comparison::gt:
				return  !(under_approx[bvID]>to);
			case Comparison::geq:
				return ! (under_approx[bvID]>=to);
			case Comparison::lt:
				return !(over_approx[bvID]<to);
			case Comparison::leq:
			default:
				return !(over_approx[bvID]<=to);
		}
	}

	Lit decideBV(Comparison op, int bvID, Weight to){
		switch(op){
			case Comparison::gt:
				if (under_approx[bvID]>to)
					return lit_Undef;
				break;
			case Comparison::geq:
				if (under_approx[bvID]>=to)
					return lit_Undef;
				break;
			case Comparison::lt:
				if (over_approx[bvID]<to)
					return lit_Undef;
				break;
			case Comparison::leq:
				if (over_approx[bvID]<=to)
					return lit_Undef;
				break;
			default:
				return lit_Undef;
		}
		while(eq_bitvectors[bvID]!=bvID)
			bvID=eq_bitvectors[bvID];
		Lit l = lit_Undef;

		while (S->decisionLevel() > trail_lim.size()) {
			newDecisionLevel();
		}

		if (opt_decide_bv_intrinsic){
			Weight under_old = under_approx[bvID];
			Weight over_old = over_approx[bvID];

			Cause under_cause_old = under_causes[bvID];
			Cause over_cause_old = over_causes[bvID];

			switch(op){
				case Comparison::gt:
					under_approx[bvID]=to+1;
					under_causes[bvID].clear();
					under_causes[bvID].type=OperationType::cause_is_decision;
					under_causes[bvID].index=decisionLevel();
					break;
				case Comparison::geq:
					under_approx[bvID]=to;
					under_causes[bvID].clear();
					under_causes[bvID].type=OperationType::cause_is_decision;
					under_causes[bvID].index=decisionLevel();
					break;
				case Comparison::lt:
					over_approx[bvID]=to-1;
					over_causes[bvID].clear();
					over_causes[bvID].type=OperationType::cause_is_decision;
					under_causes[bvID].index=decisionLevel();
					break;
				case Comparison::leq:
					over_approx[bvID]=to;
					over_causes[bvID].clear();
					over_causes[bvID].type=OperationType::cause_is_decision;
					under_causes[bvID].index=decisionLevel();
					break;
				default:
					throw std::runtime_error("Bad bv decision!");
			}

			if ((under_old != under_approx[bvID]) || (over_old != over_approx[bvID])){
				assert(under_approx[bvID]>=under_old);
				assert( over_approx[bvID]<=over_old);
				assert(analysis_trail_pos==trail.size()-1);
				trail.push({bvID,under_old, over_old, under_approx[bvID], over_approx[bvID], under_cause_old, over_cause_old, under_causes[bvID],over_causes[bvID]});
				assert(trail.last().isBoundAssignment());
				assert(trail.last().bvID==bvID);
				assert(trail.last().new_over ==  over_approx[bvID]);
				assert(trail.last().new_under == under_approx[bvID]);
				assert(trail.last().previous_over == over_old);
				assert(trail.last().previous_under == under_old);
				analysis_trail_pos=trail.size()-1;

				bv_needs_propagation[bvID]=true;
				if(! alteredBV[bvID]){
					alteredBV[bvID]=true;
					altered_bvs.push(bvID);
				}


				if(hasTheory(bvID))
					getTheory(bvID)->enqueueBV(bvID);//only enqueue the bitvector in the subtheory _after_ it's approximation has been updated!

				requiresPropagation=true;

				S->needsPropagation(getTheoryIndex());

				Lit d=lit_Undef;
				if((d = getComparison(op, bvID, to))!=lit_Undef){
					//printf("theory decision %d at level %d\n", dimacs(toSolver(d)),decisionLevel());
					return toSolver(d);
				}else{
					Lit l = S->theoryDecisionLit(getTheoryIndex());//for now, using lit_Error to signify a decision with no associated literal... is there a better option for this?
					//printf("theory decision %d at level %d\n", dimacs(l),decisionLevel());
					return l;
				}
			}else{
				under_causes[bvID]=under_cause_old;
				over_causes[bvID]=over_cause_old;
			}
			return lit_Undef;
		}else if(opt_decide_bv_bitwise){
			Weight refined_under = refine_ubound(bvID, to);
			if(refined_under>to)//can this ever not be the case?
				to = refined_under;
			//find the highest order bit of edgeID that is unassigned, and assign it to match edgeWeight.
			vec<Lit> & bits =  getBits(bvID);
			Weight bit_under = 0;
			for(int i = bits.size()-1;i>=0;i--){
				Lit b = bits[i];
				if(value(b)==l_Undef){
					Weight bit = 1L<<i;

					bool positive=true;
					if (bit_under+ bit > to ||  bit_under+ bit > over_approx[bvID]){
						positive=false;
					}
					l= positive ? b:~b;
					break;
				}else if (value(b)==l_True){
					Weight bit = 1L<<i;
					bit_under+=bit;
				}else{

				}
			}
		}else{
			l = newComparison(op, bvID, to,var_Undef,opt_cmp_lits_decidable);
		}
		return toSolver(l);
	}

	Lit decideTheory() {

		return lit_Undef;
	}
	bool supportsDecisions() {
		return false;
	}



	void newDecisionLevel() {
		trail_lim.push(trail.size());
	}
	;

	void buildReason(Lit p, vec<Lit> & reason,CRef marker) {
		static int iter = 0;
		if(++iter==39){//17
			int a =1;
		}
		assert(value(p)!=l_False);

		assert(marker != CRef_Undef);
		int pos = CRef_Undef - marker;
		//int d = marker_map[pos];
		//double initial_start = rtime(1);
		double start = rtime(1);

		rewind_trail_pos(trail.size());

		//if the reason is being constructed eagerly, then p won't be assigned yet, and so wont be on the trail, so we skip this.
		rewindUntil(var(p));

		assert(value(p)!=l_False);
		//now that we have backtracked, we need to update the under/over approximations.
		//there is likely room to improve this, so that only relevant bitvectors are updated, but it might be
		//complicated to do so... Note that there ARE cases where bitvectors beyond the one directly responsible for p may need to be updated:
		//if p is caused by the addition of two other bitvectors, for example, then they both need to be updated.
	/*	for(int bvID = 0;bvID<bitvectors.size();bvID++){
			updateApproximations(bvID);
		}*/
		assert(hasOperation(p));
		Operation & op = getOperation(p);
		op.buildReason(p,marker,reason);
		//note: the reason has already been transformed into the solvers variable namespace at this point,
		//do _not_ call 'toSolver' again


		double finish = rtime(1);
		stats_reason_time += finish - start;
		stats_num_reasons++;
		//stats_reason_initial_time+=start-initial_start;
		if(reason.size()<2){
			int a=1;
		}

	}


	void preprocess() {
		if(const_true==lit_Undef)
			const_true=True();

	}
	void setLiteralOccurs(Lit l, bool occurs) {
		/*if (isEdgeVar(var(l))) {
			//don't do anything
		} else {
			//this is a graph property detector var
			if (!sign(l) && vars[var(l)].occursPositive != occurs)
				detectors[getDetector(var(l))]->setOccurs(l, occurs);
			else if (sign(l) && vars[var(l)].occursNegative != occurs)
				detectors[getDetector(var(l))]->setOccurs(l, occurs);
		}*/
		
	}

	//Can only be used to make the assignment to a BV _more_ precise than it previously was.
	//Returns false if the assignment leads to a conflict,
	bool assignBV(int bvID, Comparison op, Weight to, Operation & cause){
		while(eq_bitvectors[bvID]!=bvID)
			bvID=eq_bitvectors[bvID];

		Weight under_old = under_approx[bvID];
		Weight over_old = over_approx[bvID];

		Weight under_new=under_old;
		Weight over_new=over_old;



		assert_in_range(under_new,bvID);
		assert_in_range(over_new,bvID);
		Cause under_cause_old = under_causes[bvID];
		Cause over_cause_old = over_causes[bvID];
		Cause over_cause_new = over_cause_old;
		Cause under_cause_new = under_cause_old;



		if (op==Comparison::lt){
			op=Comparison::leq;
			assert(to>0);
			to-=1;
			assert(to>=0);
		}

		if (op==Comparison::gt){
			op=Comparison::geq;
			to+=1;
			assert(to>0);
		}

		bool changed=false;
		if(op==Comparison::geq){
			if(under_approx[bvID]>= to){
				//this assignment leads to no refinement
				return true;
			}else{
				changed=true;
				under_new = to;
				under_cause_new.type=cause.getType();
				under_cause_new.index = cause.getID();
			}
		}else if(op==Comparison::leq){
			if(over_approx[bvID]<=to ){
				//this assignment leads to no refinement
				return true;
			}else{
				changed=true;
				over_new = to;
				over_cause_new.type=cause.getType();
				over_cause_new.index = cause.getID();
			}
		}


		if(changed){
			assert_in_range(under_new,bvID);
			assert_in_range(over_new,bvID);
			under_approx[bvID]=under_new;
			over_approx[bvID]=over_new;
			under_causes[bvID]= under_cause_new;
			over_causes[bvID]= over_cause_new;

			if(decisionLevel()==0){
				under_approx0[bvID]=under_approx[bvID];
				over_approx0[bvID]=over_approx[bvID];
				if (under_approx[bvID]==over_approx[bvID]){
					assert(!bvconst[bvID]);
					bvconst[bvID]=true;//this bitvector is a known constant value in the formula.
					n_consts++;
				}
			}

			assert(under_new>=under_old);
			assert(over_new<=over_old);
			assert(analysis_trail_pos==trail.size()-1);
			trail.push({bvID,under_old, over_old, under_new, over_new, under_cause_old, over_cause_old, under_cause_new,over_cause_new});
			assert(trail.last().isBoundAssignment());
			assert(trail.last().bvID==bvID);
			assert(trail.last().new_over == over_new);
			assert(trail.last().new_under == under_new);
			assert(trail.last().previous_over == over_old);
			assert(trail.last().previous_under == under_old);
			analysis_trail_pos=trail.size()-1;

			writeBounds(bvID);

			bv_needs_propagation[bvID]=true;
			if(! alteredBV[bvID]){
				alteredBV[bvID]=true;
				altered_bvs.push(bvID);
			}
			requiresPropagation=true;
			S->needsPropagation(getTheoryIndex());

			if(under_new>over_new)
				return false;
		}
		return true;

	}

	void enqueueTheory(Lit l) {
		Var v = var(l);
		rewind_trail_pos(trail.size());
		int lev = level(v);
		while (lev > trail_lim.size()) {
			newDecisionLevel();
		}
		
		if (assigns[var(l)] != l_Undef) {
			return;			//this is already enqueued.
		}
		assert(assigns[var(l)]==l_Undef);
		assigns[var(l)] = sign(l) ? l_False : l_True;
		requiresPropagation = true;
		//printf("enqueue %d\n", dimacs(l));
#ifndef NDEBUG
		{
			for (int i = 0; i < trail.size(); i++) {
				assert(trail[i].var != v);
			}
		}
#endif
		assert(analysis_trail_pos==trail.size()-1);
		//if(isBVVar(var(l))){
		if(hasOperation(l)){
			Operation & op = getOperation(l);
			trail.push( { true, !sign(l),op.getID(), v});
			op.enqueue(l,true);

		}else{
			trail.push( { false, !sign(l),-1, v});
		}
		analysis_trail_pos=trail.size()-1;

		//}
	}

private:
	Weight lowest(int bvID);
	Weight highest(int bvID);
	Weight refine_lbound(int bvID, Weight bound, Var ignore_bit=var_Undef);
	Weight refine_ubound(int bvID, Weight bound, Var ignore_bit=var_Undef);
	Weight refine_ubound_check(int bvID, Weight bound, Var ignore_bit);
	Weight refine_lbound_check(int bvID, Weight bound, Var ignore_bit);
	void dbg_evaluate(int bvID, int pos,vec<Weight> & vals,Weight val);
public:
	bool updateApproximations(int bvID, int ignoreCID=-1, Var ignore_bv=var_Undef){
		if(isConst(bvID))
			return false;
		if(bvID==1288){
			int a=1;
		}
		double update_start_time= rtime(3);
		statis_bv_updates++;
		static int iter = 0;
		++iter;

#ifndef NDEBUG
/*		for(int i = 0;i<vars.size();i++){
			if(value(i)==l_True){
				std::cout << "1";
			}else if (value(i)==l_False){
				std::cout << "0";
			}else{
				std::cout << "x";
			}
		}
		std::cout<<"\n";*/
#endif



		Weight under_old = under_approx[bvID];
		Weight over_old = over_approx[bvID];
		vec<Lit> & bv = bitvectors[bvID];
		Weight under_new=under_old;
		Weight over_new=over_old;

		bool any_changed=false;

		if(eq_bitvectors[bvID]!=bvID){
			//this bitvector is equivalent to some other (lower index) bv, so just copy its values.
			int eqBV = eq_bitvectors[bvID];
			assert(eqBV<bvID);
			Cause under_cause_old = under_causes[bvID];
			Cause over_cause_old = over_causes[bvID];
			under_new = under_approx[eqBV];
			over_new = over_approx[eqBV];
			under_approx[bvID]=under_approx[eqBV];
			over_approx[bvID]=over_approx[eqBV];
			under_approx0[bvID]=under_approx0[eqBV];
			over_approx0[bvID]=over_approx0[eqBV];

			return false;
		}
		assert_in_range(under_new,bvID);
		assert_in_range(over_new,bvID);
		Cause under_cause_old = under_causes[bvID];
		Cause over_cause_old = over_causes[bvID];
		Cause over_cause_new;
		Cause under_cause_new;

		for(int opID:operation_ids[bvID]){
			if(opID==ignoreCID)
				continue;
			getOperation(opID).updateApprox(ignore_bv, under_new,  over_new,  under_cause_new,  over_cause_new);
		}

		assert_in_range(under_new,bvID);
		assert_in_range(over_new,bvID);
		//special handling for comparisons against constants, which are kept in sorted order, allowing some short-circuiting
		for(int i = compares[bvID].size()-1;i>=0;i--){
			int cID = compares[bvID][i];
			if(cID==ignoreCID){
				continue;
			}
			ComparisonOp & c = (ComparisonOp&) getOperation(cID);
			c.updateApprox(ignore_bv, under_new,  over_new,  under_cause_new,  over_cause_new,true);
		}

		for(int cID:compares[bvID]){
			if(cID==ignoreCID){
				continue;
			}
			ComparisonOp & c = (ComparisonOp&) getOperation(cID);
			c.updateApprox(ignore_bv, under_new,  over_new,  under_cause_new,  over_cause_new,false);
		}

		assert_in_range(under_new,bvID);
		assert_in_range(over_new,bvID);

		if(over_new>over_approx0[bvID]){
			over_new=over_approx0[bvID];
			over_cause_new.clear();
		}
		if(under_new<under_approx0[bvID]){
			under_new=under_approx0[bvID];
			under_cause_new.clear();
		}
		//If no improvement was made in these bounds, then just stick with the old versions (and their causes).
		if(over_new>=over_old){
			over_new=over_old;
			over_cause_new=over_cause_old;
		}else{
			assert(over_cause_new.hasCause());
		}
		if(under_new<=under_old){
			under_new=under_old;
			under_cause_new=under_cause_old;
		}

		assert_in_range(under_new,bvID);
		assert_in_range(over_new,bvID);
		Weight refined_over = refine_lbound(bvID, over_new);
		if(refined_over>-1 && refined_over< over_new){
			//std::cout<< "Refined overapprox for bv " << bvID << " from " << over_new << " to " << refined_over << "\n";

			if(over_new<over_old){
				//need to record this previous change as a separate entry in the trail, for conflict analysis later...
				any_changed=true;
				assert(over_new<=over_old);
				assert(analysis_trail_pos==trail.size()-1);
				trail.push({bvID,under_old, over_old, under_old, over_new, under_cause_old, over_cause_old, under_cause_old,over_cause_new});
				assert(trail.last().isBoundAssignment());
				assert(trail.last().bvID==bvID);
				assert(trail.last().new_over == over_new);
				assert(trail.last().new_under == under_old);//intentionally not changing under here.
				assert(trail.last().previous_over == over_old);
				assert(trail.last().previous_under == under_old);
				analysis_trail_pos=trail.size()-1;


				//ONLY update over_old to over new here.
				over_old=over_new;
				over_cause_old=over_cause_new;
			}

			over_new=refined_over;
			over_cause_new.clear();
			over_cause_new.type =OperationType::refined_cause;
		}
		Weight refined_under = refine_ubound(bvID, under_new);
		if(refined_under>-1  && refined_under> under_new){
			//std::cout<< "Refined underapprox for bv " << bvID << " from " << under_new<< " to " << refined_under << "\n";

			if(under_new>under_old){
				//need to record this previous change as a separate entry in the trail, for conflict analysis later...
				assert(under_new>=under_old);
				assert(over_new<=over_old);
				assert(analysis_trail_pos==trail.size()-1);
				trail.push({bvID,under_old, over_old, under_new, over_new, under_cause_old, over_cause_old, under_cause_new,over_cause_new});
				assert(trail.last().isBoundAssignment());
				assert(trail.last().bvID==bvID);
				assert(trail.last().new_over == over_new);
				assert(trail.last().new_under == under_new);
				assert(trail.last().previous_over == over_old);
				assert(trail.last().previous_under == under_old);
				analysis_trail_pos=trail.size()-1;
				any_changed=true;
				under_old=under_new;
				over_old=over_new;
				under_cause_old=under_cause_new;
				over_cause_old=over_cause_new;
			}

			under_new=refined_under;
			under_cause_new.clear();
			under_cause_new.type =OperationType::refined_cause;
		}
		assert_in_range(under_new,bvID);
		assert_in_range(over_new,bvID);

		//int width = bitvectors[bvID].size();
		//Weight max_val = ((1L)<<width)-1;


		assert_in_range(under_new,bvID);
		assert_in_range(over_new,bvID);
		under_approx[bvID]=under_new;
		over_approx[bvID]=over_new;
		under_causes[bvID]= under_cause_new;
		over_causes[bvID]= over_cause_new;

		if(decisionLevel()==0){
			under_approx0[bvID]=under_approx[bvID];
			over_approx0[bvID]=over_approx[bvID];
			if (under_approx[bvID]==over_approx[bvID]){
				assert(!bvconst[bvID]);
				bvconst[bvID]=true;//this bitvector is a known constant value in the formula.
				n_consts++;
			}
		}

		bool changed = (under_old != under_approx[bvID]) || (over_old != over_approx[bvID]);
		if(changed){
			any_changed=true;
			assert(under_new>=under_old);
			assert(over_new<=over_old);
			assert(analysis_trail_pos==trail.size()-1);
			trail.push({bvID,under_old, over_old, under_new, over_new, under_cause_old, over_cause_old, under_cause_new,over_cause_new});
			assert(trail.last().isBoundAssignment());
			assert(trail.last().bvID==bvID);
			assert(trail.last().new_over == over_new);
			assert(trail.last().new_under == under_new);
			assert(trail.last().previous_over == over_old);
			assert(trail.last().previous_under == under_old);
			analysis_trail_pos=trail.size()-1;

		}else{
			//ensure that the cause isn't altered if the approx was not changed.
			under_causes[bvID] = under_cause_old;
			over_causes[bvID] = over_cause_old;
		}
	/*	if(any_changed && getSymbol(bvID)){
			std::cout<< "q bv " << getSymbol(bvID) << " " << under_approx[bvID] << " <= bv <=" <<  over_approx[bvID] << "\n" ;
		}*/
		writeBounds(bvID);
		stats_update_time+= rtime(3) -update_start_time;
		return 	any_changed;//return whether either weight has changed.
	}

	void writeBounds(int bvID){
		if(opt_write_learnt_clauses && opt_write_bv_bounds){
			static int bound_num=0;
			if(++opt_n_learnts==44231){
				int a=1;
			}

			int bv =  unmapBV( bvID);
			fprintf(opt_write_learnt_clauses,"learnt bound ");
			fprintf(opt_write_learnt_clauses," %d <= ", unmapBV( bvID));

			std::stringstream ss;
			ss<<over_approx[bvID] << " ";
			//fprintf(opt_write_learnt_clauses,"learnt fact bv %d %s %d\n", bvID, ss.str().c_str(),dimacs(toSolver(l)) );
			fprintf(opt_write_learnt_clauses,"%s",ss.str().c_str());

			if(over_approx[bvID]!=over_approx0[bvID]){

				for (int i = 0;i<nVars();i++){
					if(value(i)!=l_Undef){
						Lit l = mkLit(i,value(i)==l_True);
						assert(value(l)==l_False);
						fprintf(opt_write_learnt_clauses," %d",(dimacs(S->unmap(toSolver(l)))));
					}
				}
			}
			fprintf(opt_write_learnt_clauses," 0\n");
			++bound_num;
			if(bound_num==10){
				int a=1;
			}

			if(++opt_n_learnts==44231){
				int a=1;
			}
			fprintf(opt_write_learnt_clauses,"learnt bound ");
			fprintf(opt_write_learnt_clauses," %d >= ", unmapBV( bvID));
			//std::cout<<under_approx[bvID] << " ";
			std::stringstream ss2;
			ss2<<under_approx[bvID] << " ";
			//fprintf(opt_write_learnt_clauses,"learnt fact bv %d %s %d\n", unmapBV( bvID), ss.str().c_str(),(dimacs(S.unmaptoSolver(l)) ));
			fprintf(opt_write_learnt_clauses,"%s",ss2.str().c_str());
			if(under_approx[bvID]!=under_approx0[bvID]){
				for (int i = 0;i<nVars();i++){
					if(value(i)!=l_Undef){
						Lit l = mkLit(i,value(i)==l_True);
						assert(value(l)==l_False);
						fprintf(opt_write_learnt_clauses," %d",dimacs(S->unmap(toSolver(l))));
					}
				}
			}
			fprintf(opt_write_learnt_clauses," 0\n");
			if(++bound_num==10){
				int  a=1;
			}
			//fprintf(opt_write_learnt_clauses,"%d\n", bound_num);
			fflush(opt_write_learnt_clauses);
		}
	}

	Weight & getUnderApprox(int bvID, bool level0=false){
		if(level0){
			return under_approx0[bvID];
		}else{
			return under_approx[bvID];
		}
	}
	Weight & getOverApprox(int bvID, bool level0=false){
		if(level0){
			return over_approx0[bvID];
		}else{
			return over_approx[bvID];
		}
	}
	int getWidth(int bvID){
		return getBits(bvID).size();
	}
	vec<Lit> & getBits(int bvID){
		//can this be avoided?
		while(eq_bitvectors[bvID]!=bvID)
			bvID=eq_bitvectors[bvID];
		return bitvectors[bvID];
	}

	void buildTrivialClause(vec<Lit> & conflict){
		for (int i = 0;i<nVars();i++){
			if(value(i)!=l_Undef){
				Lit l = mkLit(i,value(i)==l_True);
				assert(value(l)==l_False);
				conflict.push(toSolver(l));
			}
		}
	}

	bool checkApproxUpToDate(int bvID, Weight * under_store=nullptr,Weight*over_store=nullptr){
#ifndef NDEBUG
		if(bvID==2966){
			int a=1;
		}
		Weight under =under_approx0[bvID];
		Weight over=over_approx0[bvID];
		vec<Lit> & bv = bitvectors[bvID];
		if(eq_bitvectors[bvID]!=bvID && eq_bitvectors[bvID]>-1 ){
			int eqID = eq_bitvectors[bvID];
			under=under_approx[eqID];
			over=over_approx[eqID];
			assert(over_approx[bvID]==over);
			assert(under_approx[bvID]==under);
			return true;
		}else{


			for (int opID:operation_ids[bvID]){
				getOperation(opID).checkApproxUpToDate(under,over);
			}
			for(int i = compares[bvID].size()-1;i>=0;i--){
				int cID = compares[bvID][i];

				ComparisonOp & c = (ComparisonOp&) getOperation(cID);
				c.checkApproxUpToDate(under,over);
			}
			Weight refined_over = refine_lbound(bvID, over);
			if(refined_over>-1 && refined_over< over){

				over=refined_over;
			}
			Weight refined_under = refine_ubound(bvID,under);
			if(refined_under>-1  && refined_under> under){

				under=refined_under;
			}

			int width = bitvectors[bvID].size();
			Weight max_val = (1L<<width)-1;
			if(under>over_approx0[bvID]){
				under=over_approx0[bvID];
			}
			if(over>over_approx0[bvID]){
				over=over_approx0[bvID];
			}
			if(under<under_approx0[bvID]){
				under=under_approx0[bvID];
			}
			if(over<under_approx0[bvID]){
				over=under_approx0[bvID];
			}
		}
		if (under_causes[bvID].type !=OperationType::cause_is_decision && under_causes[bvID].type !=OperationType::cause_is_theory)
			assert(under==under_approx[bvID]);
		if(over_causes[bvID].type !=OperationType::cause_is_decision  && over_causes[bvID].type !=OperationType::cause_is_theory)
			assert(over==over_approx[bvID]);

		assert(under_approx0[bvID]<=under_approx[bvID]);
		assert(over_approx0[bvID]>=over_approx[bvID]);

		if(eq_bitvectors[bvID]!=bvID && eq_bitvectors[bvID]>-1 ){
			int eqID = eq_bitvectors[bvID];
			Weight under_expect=under_approx[eqID];
			Weight over_expect=over_approx[eqID];
			assert(under_expect==under);
			assert(over_expect==over);
			Weight under_approx0_expect = under_approx0[eqID];
			Weight over_approx0_expect = over_approx0[eqID];

			assert(under_approx0_expect==under_approx0[bvID]);
			assert(over_approx0_expect==over_approx0[bvID]);
		}

		if(under_store){
			(*under_store)=under;
		}
		if(over_store){
			(*over_store)=over;
		}
#endif
		return true;
	}

	bool dbg_synced(){
#ifndef NDEBUG
		for(Var v = 0;v<nVars();v++){
			assert(value(v)==dbg_value(v));
		}
#endif
		return true;
	}

	bool checkAllApproxUpToDate(){
#ifndef NDEBUG
		for (int i = 0;i< bitvectors.size();i++){
			assert(checkApproxUpToDate(i));
		}
#endif
		return true;
	}

/*	bool dbg_uptodate(){
#ifndef NDEBUG
		assert(checkAllApproxUpToDate());
		for(int cID = 0;cID<comparisons.size();cID++){
			ComparisonID & c = comparisons[cID];
			int bvID = c.bvID;
			Comparison op =c.op();
			if(!c.bvCompare()){
				Weight w = c.w;
				if(op==Comparison::lt && getOverApprox(bvID)<w){

				}

			}else{

			}

		}
#endif
		return true;
	}*/
	void addAlteredBV( int newBV){
		if (altered_bvs.size()==0){
			altered_bvs.push(newBV);
			alteredBV[newBV]=true;
		}else{
			int bvID = altered_bvs.last();
			if(!alteredBV[newBV]){
				alteredBV[newBV]=true;
				altered_bvs.last()=newBV;
				assert(altered_bvs.last()==newBV);
				altered_bvs.push(bvID);
				assert(altered_bvs.last()==bvID);
			}
		}
	}

	bool propagateTheory(vec<Lit> & conflict){
		return propagateTheory(conflict,false);
	}
	bool propagateTheory(vec<Lit> & conflict, bool force_propagation) {
		static int realprops = 0;
		stats_propagations++;
		if(stats_propagations==55){
			int a=1;
		}
		if (!force_propagation && !requiresPropagation ) {
			stats_propagations_skipped++;
			assert(dbg_uptodate());
			return true;
		}

		propagations++;

		if (!force_propagation && (propagations % opt_bv_prop_skip != 0)){
			stats_propagations_skipped++;
			return true;
		}

		rewind_trail_pos(trail.size());
		if(++realprops==14){
			int a =1;
		}
		//printf("bv prop %d\n",stats_propagations);
		if(stats_propagations==22){
			int a =1;
		}
		bool any_change = false;
		double startproptime = rtime(2);
		//static vec<int> detectors_to_check;

		while (S->decisionLevel() > trail_lim.size()) {
			newDecisionLevel();
		}
		
		conflict.clear();
		

		while(altered_bvs.size()){
			int bvID = altered_bvs.last();

			if(eq_bitvectors[bvID]!=bvID)
			{
				altered_bvs.pop(); //this isn't safe to do, because recently enforced comparisons need to be propagated, even if the under/over approx didn't change.
				alteredBV[bvID]=false;
				continue;
			}

		//for(int bvID = 0;bvID<bitvectors.size();bvID++){
			assert(alteredBV[bvID]);
			Weight  underApprox_prev = under_approx[bvID];
			Weight  overApprox_prev = over_approx[bvID];
			Cause prev_under_cause = under_causes[bvID];
			Cause prev_over_cause = over_causes[bvID];

			bool changed = updateApproximations(bvID);//can split this into changedUpper and changedLower...
			changed |=bv_needs_propagation[bvID];
			if(!changed){
				assert( under_approx[bvID]<= over_approx[bvID]);
				stats_bv_skipped_propagations++;
				altered_bvs.pop(); //this isn't safe to do, because recently enforced comparisons need to be propagated, even if the under/over approx didn't change.
				alteredBV[bvID]=false;
				continue;
			}



			stats_bv_propagations++;
			vec<Lit> & bv = bitvectors[bvID];

			/*printf("iter %d, bv %d, under ",realprops , bvID); //: %d, over %d\n", bvID, underApprox,overApprox);
			std::cout<<underApprox << " over ";
			std::cout<<overApprox << "\n";
			fflush(stdout);*/


			//also need to iterate through the additions that this bv is an argument of...
			for(int opID:operation_ids[bvID]){
				if(opID==1288){
					int a=1;
				}
				if(!getOperation(opID).propagate(changed,conflict)){
					return false;
				}
			}

			vec<int> & compare = compares[bvID];
			//update over approx lits.
			//It might be worth doing a binary search here...
			for(int i = 0;i<compare.size();i++){
				int cID = compare[i];

				assert( getOperation(cID).getType()==OperationType::cause_is_comparison);
				ComparisonOp & c = (ComparisonOp &) getOperation(cID);

				if(!c.propagate(changed,conflict,true))
					return false;
			}

			for(int i=compare.size()-1;i>=0;i--){
				int cID = compare[i];
				assert( getOperation(cID).getType()==OperationType::cause_is_comparison);
				ComparisonOp & c = (ComparisonOp &) getOperation(cID);

				if(!c.propagate(changed,conflict,false))
					return false;
			}

			//comparisons to bitvectors.
			/*vec<int> & bvcompare = bvcompares[bvID];

			for(int i = 0;i<bvcompare.size();i++){
				int cID = bvcompare[i];
				assert( getOperation(cID).getType()==OperationType::cause_is_bv_comparison);
				ComparisonBVOp & c = (ComparisonBVOp &) getOperation(cID);

				if(!c.propagate(changed,conflict,false))
					return false;

			}

			for(int i=bvcompare.size()-1;i>=0;i--){
				int cID = bvcompare[i];
				assert( getOperation(cID).getType()==OperationType::cause_is_bv_comparison);
				ComparisonBVOp & c = (ComparisonBVOp &) getOperation(cID);

				if(!c.propagate(changed,conflict,false))
					return false;

			}*/
			bv_needs_propagation[bvID]=false;
			if(changed){
				if(hasTheory(bvID))
					getTheory(bvID)->enqueueBV(bvID);//only enqueue the bitvector in the subtheory _after_ it's approximation has been updated!
			}

			if(changed){
				//then any additions this is an argument of need to be updated.
				for(int changedID: cause_set[bvID]){
					under_approx[changedID]=under_approx[bvID];
					over_approx[changedID]=over_approx[bvID];
					under_approx0[changedID]=under_approx0[bvID];
					over_approx0[changedID]=over_approx0[bvID];
					if(hasTheory(changedID))
						getTheory(changedID)->enqueueBV(changedID);//only enqueue the bitvector in the subtheory _after_ it's approximation has been updated!
				}
			}

			assert(altered_bvs.last()==bvID);
			assert(alteredBV[bvID]==true);
			altered_bvs.pop();
			alteredBV[bvID]=false;
		}
		
		requiresPropagation = false;
		propagationtime += rtime(2) - startproptime;;
		assert(dbg_uptodate());
		if(first_propagation){
			first_propagation=false;
			n_starting_consts=n_consts;
		}
		return true;
	};

	Weight ceildiv(Weight a, Weight b);



	void buildComparisonReasonBV(Comparison op, int bvID,int comparebvID, vec<Lit> & conflict){
		dbg_no_pending_analyses();
		assert(eq_bitvectors[bvID]==bvID);


		//rewind_trail_pos(trail_pos);
		stats_build_value_bv_reason++;
		Weight  over_cur = over_approx[bvID];
		Weight  under_cur = under_approx[bvID];


		Weight  over_comp = over_approx[comparebvID];
		Weight  under_comp = under_approx[comparebvID];

		if (op==Comparison::lt){
			assert(over_cur<under_comp);

		}
		else if (op==Comparison::leq){
			assert(over_cur<=under_comp);

		}
		else if (op==Comparison::geq){
			assert(under_cur>=over_comp);

		}
		else if (op==Comparison::gt){
			assert(under_cur>over_comp);
		}

		if(isConst(bvID) && isConst(comparebvID)){
			//both bitvectors are constants, so no reason is required
			return;
		}else if (isConst(comparebvID)){
			assert(under_comp==over_comp);
			//if the other bitvector is a constant,
			//then the reason for the assignment is because this bitvector is (<,<=,>,>=) than the underapproximation of that constant.
			//newComparison(op,bvID,under_comp);
			//int cID = getComparisonID(op,bvID,under_comp);
			analyzeValueReason(op,bvID,under_comp,conflict);
		}else if (isConst(bvID)){
			assert(under_cur==over_cur);
			//if the other bitvector is a constant,
			//then the reason for the assignment is because this bitvector is (<,<=,>,>=) than the underapproximation of that constant.
			op=~op;
			//newComparison(op,comparebvID,over_cur);
			//int cID = getComparisonID(op,comparebvID,over_cur);
			analyzeValueReason(op,comparebvID,over_cur,conflict);
		}else{
			//neither bitvector is constant. Pick a value between them, and learn relative to that.
			Weight midval;
			Comparison cOp;

			if (op==Comparison::lt){
				assert(over_cur>=0);
				assert(under_comp>=0);
				assert(over_cur<under_comp);

				midval = ceildiv( under_comp- over_cur,2)+over_cur; //integer ceiling division
				assert(midval<=under_comp);
				assert(over_cur<midval);
				cOp=Comparison::geq;
				assert(midval!=over_cur);//because we took the integer ceiling division
			}
			else if (op==Comparison::leq){
				assert(over_cur<=under_comp);
				Weight m = under_comp- over_cur;
				midval = ceildiv( under_comp- over_cur,2)+over_cur; //integer ceiling division
				assert(midval<=under_comp);
				assert(over_cur<=midval);
				cOp=Comparison::geq;
			}
			else if (op==Comparison::geq){
				assert(under_cur>=over_comp);
				midval = (under_cur- over_comp)/2 + over_comp;//integer floor division
				assert(midval>=over_comp);
				assert(under_cur>=midval);
				cOp=Comparison::leq;
			}
			else{
				assert((op==Comparison::gt));
				assert(under_cur>over_comp);
				midval = (under_cur- over_comp)/2 + over_comp;//integer floor division
				assert(midval>=over_comp);
				assert(under_cur>midval);
				cOp=Comparison::leq;
			}

			//newComparison(op,bvID,midval);
			//int cID = getComparisonID(op,bvID,midval);
			analyzeValueReason(op,bvID,midval,conflict);

			//newComparison(op,comparebvID,midval);
			//cID = getComparisonID(op,comparebvID,midval);
			analyzeValueReason(cOp,comparebvID,midval,conflict);
		}
		analyze(conflict);

	}

	void analyzeValueReason(Comparison op, int bvID, Weight  to,  vec<Lit> & conflict){
		if(isConst(bvID)){
			// a constant bitvector needs no reason
			return;
		}
		while(eq_bitvectors[bvID]!=bvID)
			bvID=eq_bitvectors[bvID];
		writeBounds(bvID);
		stats_build_value_reason++;
		static int iter = 0;
		++iter;
		//printf("reason %d: %d\n",iter,bvID);
		if(iter==430){
			int a=1;
		}

		//rewind_trail_pos(trail_pos);
		//trail_pos  = rewindUntil(bvID,op,to);


		bool compare_over;

		vec<Lit> & bv = bitvectors[bvID];
		Weight  over_cur = over_approx[bvID];
		Weight  under_cur = under_approx[bvID];
		//assert(checkApproxUpToDate(bvID));
		if (op==Comparison::lt){
			assert(over_cur<to);
			compare_over=true;
		/*	assert(to>0);
			op=Comparison::leq;
			to-=1;*/

		}
		if (op==Comparison::leq){
			assert(over_cur<=to);
			compare_over=true;
			//to = over_cur;
		}
		if (op==Comparison::gt){
			assert(under_cur>to);
			compare_over=false;
	/*		op=Comparison::geq;
			to+=1;*/
		}
		if (op==Comparison::geq){
			assert(under_cur>=to);
			compare_over=false;
			//to = under_cur;
		}
		if(!compare_over && under_approx[bvID]<= getUnderApprox(bvID, true) ){
			//no reason necessary; this is the lowest possible value.
			return;
		}else if(compare_over && over_approx[bvID]>=getOverApprox(bvID, true)  ){
			//no reason necessary; this is the lowest possible value.
			return;
		}
		if(opt_write_learnt_clauses && opt_write_bv_analysis){
			static int n_analysis=0;
			if(++n_analysis==430){
				int a =1;
			}
			std::stringstream ss;
			ss<< op << " " << to;
			fprintf(opt_write_learnt_clauses,"learnt analysis bv %d %s ",unmapBV( bvID),ss.str().c_str());
		}

		if(compare_over)
			assert(over_causes[bvID].hasCause());
		if(!compare_over)
			assert(under_causes[bvID].hasCause());

		if (compare_over &&  over_causes[bvID].type==OperationType::cause_is_decision){
			//the reason that the bvID's over approx is <= its current value
			//is because it was a decision.
			//Create a literal on the fly to explain this...
			int lev = over_causes[bvID].index;
			assert(lev>=0);
			assert(lev<=decisionLevel());
			Lit reason = newComparison(Comparison::leq,bvID, over_cur,var_Undef,opt_cmp_lits_decidable);
			assert(value(reason)!=l_False);
			conflict.push(toSolver(~reason));

			//S->prependToTrail(toSolver(reason),lev);//this is a decision that was made, without a corresponding literal in the solver at the time it was made.
			// need to  ensure that this lit can be properly analyzed, so prepend it to the trail at this decision level.

		}else if (!compare_over &&  under_causes[bvID].type==OperationType::cause_is_decision){
			int lev = under_causes[bvID].index;
			assert(lev>=0);
			assert(lev<=decisionLevel());
			Lit reason = newComparison(Comparison::geq,bvID, under_cur,var_Undef,opt_cmp_lits_decidable);
			assert(value(reason)!=l_False);
			conflict.push(toSolver(~reason));

			//S->prependToTrail(toSolver(reason),lev);//this is a decision that was made, without a corresponding literal in the solver at the time it was made.
			//need to ensure that this lit can be properly analyzed, so prepend it to the trail at this decision level.

		}else if (compare_over &&  over_causes[bvID].type==OperationType::refined_cause){
			//then the reason the underapprox is too large is because of the assignment to the bits
			//can this analysis be improved upon?
			for(int i =0;i<bv.size();i++){
				Lit bl = bv[i];
				if(value(bl)==l_False){

						assert(value(bl)==l_False);
						conflict.push(toSolver(bl));

				}else if(value(bl)==l_True){

					assert(value(~bl)==l_False);
					conflict.push(toSolver(~bl));

				}
			}
		/*	assert(trail[analysis_trail_pos].isBoundAssignment());
			assert(trail[analysis_trail_pos].bvID==bvID);
			Weight prev_weight = trail[analysis_trail_pos].previous_over;*/
			//search back in the trail
			int trail_pos = analysis_trail_pos;
			rewindUntil(bvID,Comparison::leq, over_approx[bvID]);
			assert(trail[analysis_trail_pos].isBoundAssignment());
			assert(trail[analysis_trail_pos].bvID==bvID);
			Weight prev_weight = trail[analysis_trail_pos].previous_over;
			rewind_trail_pos(trail_pos);
			addAnalysis(Comparison::leq,bvID,prev_weight);


			//buildValueReason(Comparison::leq,bvID,over_approx[bvID],conflict,trail_pos-1);

		}else if (!compare_over  && under_causes[bvID].type==OperationType::refined_cause){
			//then the reason the underapprox is too large is because of the assignment to the bits

			for(int i =0;i<bv.size();i++){
				Lit bl = bv[i];
				lbool val = value(bl);
				lbool dbgval = dbg_value(bl);
				if(value(bl)==l_False){

						assert(value(bl)==l_False);
						conflict.push(toSolver(bl));

				}else if(value(bl)==l_True){

					assert(value(~bl)==l_False);
					conflict.push(toSolver(~bl));

				}
			}

			int trail_pos = analysis_trail_pos;
			rewindUntil(bvID,Comparison::geq, under_approx[bvID]);
			assert(trail[analysis_trail_pos].isBoundAssignment());
			assert(trail[analysis_trail_pos].bvID==bvID);
			Weight prev_weight = trail[analysis_trail_pos].previous_under;
			rewind_trail_pos(trail_pos);
			addAnalysis(Comparison::geq,bvID,prev_weight);

			//buildValueReason(Comparison::geq,bvID,under_approx[bvID],conflict,trail_pos-1);
		}else if (compare_over){
			getOperation(over_causes[bvID].index).analyzeReason(compare_over,op,to,conflict);
		}else if (!compare_over ){
			getOperation(under_causes[bvID].index).analyzeReason(compare_over,op,to,conflict);
		}

		if(opt_write_learnt_clauses && opt_write_bv_analysis){
			fprintf(opt_write_learnt_clauses," lits " );
			for(int i = 0;i<conflict.size();i++){
				fprintf(opt_write_learnt_clauses,"%d ",dimacs(S->unmap(toSolver(conflict[i]))));
			}
			fprintf(opt_write_learnt_clauses, "\n");
		}
	}
	inline void assert_in_range(Weight val, int bvID){
#ifndef NDEBUG
		int width = bitvectors[bvID].size();
		Weight max_val = ((1L)<<width)-1;
		assert(val>=0);
		assert(val<=max_val);
#endif
	}
	inline void clip(Weight & val, int bvID){
		clip_under(val,bvID);
		clip_over(val,bvID);
		assert_in_range(val,bvID);
	}
	inline void clip_under(Weight & val, int bvID){
		if(val<0)
			val=0;
	}
	inline void clip_over(Weight & val, int bvID){
		int width = bitvectors[bvID].size();
		Weight max_val = ((1L)<<width)-1;
		if(val>max_val)
			val=max_val;
	}

	bool addAnalysis(Comparison op, int bvID, Weight  to){
		while(eq_bitvectors[bvID]!=bvID)
			bvID=eq_bitvectors[bvID];
		if(op==Comparison::geq && to<=0){
			return false;
		}else if (op==Comparison::gt && to<0){
			return false;
		}
		int width = bitvectors[bvID].size();
		Weight max_val = ((1L)<<width)-1;
		if (op==Comparison::leq && to>= max_val){
			return false;
		}else if (op==Comparison::lt && to>max_val){
			return false;
		}
		clip(to,bvID);

		int aID = analyses.size();
		bool compare_over;
		if (op==Comparison::lt){
			compare_over=true;
			op=Comparison::leq;
			to--;
		}
		if (op==Comparison::leq){
			compare_over=true;
		}
		if (op==Comparison::geq){
			compare_over=false;
		}
		if (op==Comparison::gt){
			compare_over=false;
			op=Comparison::geq;
			to++;
		}
		//add this to bv comparison to the stack of analyses to perform.
		if(compare_over){
			if(to>= over_approx0[bvID])
				return false;//no analysis required.


			//need to check whether this analysis has already been requested for this bitvector, and if not, insert it into the analysis chain in the right position.
			//should this use a binary search?
			int cID = pending_over_analyses[bvID];

			int prevID = -1;
			while(cID>-1 && analyses[cID].value < to){
				ToAnalyze & c = analyses[cID];
				prevID = cID;
				if(cID<0 || cID==c.next_analysis){
					throw std::logic_error("Cycle in BV Theory Solver!");//this is a cycle
				}
				cID = c.next_analysis;
			}
			if(cID>=0){
				if( analyses[cID].value==to){
					//this analysis was already requested; don't do anything
					if( opt_write_learnt_clauses && opt_write_bv_analysis){
						std::stringstream ss;
						ss<< op << " " << to;
						fprintf(opt_write_learnt_clauses, " bv %d %s ",unmapBV( bvID),ss.str().c_str());
					}
					return true;
				}
			}
			if(prevID>=0){
				assert(analyses[prevID].value<to);
			}
			n_pending_analyses++;
			//insert this analysis between prevID (if it exists) and cID
			analyses.push({bvID,cID,to});//to
			if(prevID>-1){
				analyses[prevID].next_analysis=aID;
			}else{
				pending_over_analyses[bvID]=aID;
			}
		}else{
			if(to<= under_approx0[bvID])
				return false;//no analysis required.
			//need to check whether this analysis has already been requested for this bitvector, and if not, insert it into the analysis chain in the right position.
			//should this use a binary search?

			int cID = pending_under_analyses[bvID];

			int prevID = -1;
			while(cID>-1 && analyses[cID].value > to){
				ToAnalyze & c = analyses[cID];
				prevID = cID;
				if(cID<0 || cID==c.next_analysis){
					throw std::logic_error("Cycle in BV Theory Solver!");//this is a cycle
				}
				cID = c.next_analysis;
			}
			if(cID>=0){
				if( analyses[cID].value==to){
					//this analysis was already requested; don't do anything
					if( opt_write_learnt_clauses && opt_write_bv_analysis){
						std::stringstream ss;
						ss<< op << " " << to;
						fprintf(opt_write_learnt_clauses, " bv %d %s ",unmapBV( bvID),ss.str().c_str());
					}
					return true;
				}
			}
			if(prevID>=0){
				assert(analyses[prevID].value>to);
			}
			//insert this analysis between prevID (if it exists) and cID
			n_pending_analyses++;
			analyses.push({bvID,cID,to});//op
			if(prevID>-1){
				analyses[prevID].next_analysis=aID;
			}else{
				pending_under_analyses[bvID]=aID;
			}
		}

		if( opt_write_learnt_clauses && opt_write_bv_analysis){
			std::stringstream ss;
			ss<< op << " " << to;
			fprintf(opt_write_learnt_clauses, " bv %d %s ",unmapBV( bvID),ss.str().c_str());
		}
		return true;
	}

	void buildComparisonReason(Comparison op, int bvID, Weight  to,  vec<Lit> & conflict){
		dbg_no_pending_analyses();
		//rewind_trail_pos(trail.size()-1);
		int trail_pos = rewindUntil(bvID,op,to);
		//analyses.push({bvID,op,conflict});
		analyzeValueReason(op,bvID,to,conflict);
		analyze(conflict);
	}

	void analyze(vec<Lit> & conflict){
		static long iter = 0;
		int prev_pos = analysis_trail_pos;
		while(n_pending_analyses>0){

			assert(analysis_trail_pos>=0);
			if(analysis_trail_pos<0 || analysis_trail_pos>=trail.size()){
				throw std::runtime_error("Error in BV Theory Solver!");
			}
			Assignment & e = trail[analysis_trail_pos];

			if (e.isBoundAssignment()){
				int bvID = e.bvID;


				while(pending_over_analyses[bvID]>-1 && e.previous_over > analyses[pending_over_analyses[bvID]].value){
					int aID = pending_over_analyses[bvID];
					ToAnalyze & a =  analyses[aID];
					int nID = a.next_analysis;
					Weight w = a.value;
					assert(nID!=aID);
					pending_over_analyses[bvID]=nID;
					n_pending_analyses--;
					analyses[aID].clear();
					analyzeValueReason(Comparison::leq,bvID,w,conflict);
				}
				while(pending_under_analyses[bvID]>-1 && e.previous_under <  analyses[pending_under_analyses[bvID]].value){
					int aID = pending_under_analyses[bvID];
					ToAnalyze & a =  analyses[aID];
					int nID = a.next_analysis;
					Weight w = a.value;
					assert(nID!=aID);
					pending_under_analyses[bvID]=nID;
					n_pending_analyses--;
					analyses[aID].clear();
					analyzeValueReason(Comparison::geq,bvID,w,conflict);
				}
			}
			rewind_trail_pos(analysis_trail_pos-1);
			assert(analysis_trail_pos<prev_pos);
			prev_pos=analysis_trail_pos;
		}
		analyses.clear();
		dbg_no_pending_analyses();
		rewind_trail_pos(trail.size()-1);
		//now walk back through the trail to find the
	}

	inline void dbg_no_pending_analyses(){
#ifndef NDEBUG
		for (int i = 0;i<pending_under_analyses.size();i++){
			assert(pending_under_analyses[i]==-1);
			assert(pending_over_analyses[i]==-1);
		}
#endif

	}

	bool solveTheory(vec<Lit> & conflict) {
		requiresPropagation = true;		//Just to be on the safe side... but this shouldn't really be required.
		bool ret = propagateTheory(conflict);
		//Under normal conditions, this should _always_ hold (as propagateTheory should have been called and checked by the parent solver before getting to this point).
		assert(ret);
		return ret;
	}

	bool check_propagated()override{
		return !requiresPropagation;
	}
	bool check_solved() {

		for(int bvID = 0;bvID<bitvectors.size();bvID++){
			int eqBV = bvID;
			while(eq_bitvectors[eqBV]>-1 && eq_bitvectors[eqBV]!=eqBV ){
				eqBV=eq_bitvectors[eqBV];
			}

			if(over_approx[eqBV]!= under_approx[eqBV]){
				return false;
			}


		}
		for(Operation * op:operations){
			if(!op->checkSolved())
				return false;
		}

		return true;
	}
	
	bool dbg_solved() {

		return true;
	}




	void setBitvectorTheory(int bvID, int theoryID){
		theoryIds[bvID]=theoryID;
	}


	BitVector newMaxBV(int resultID,vec<int> & args){

		MinMaxData * op = new MinMaxData(*this,resultID,false);
		addOperation(resultID,op);
		for (int argID:args){
			MinMaxArg * arg = new MinMaxArg(*this,argID,op);
			addOperation(argID,arg);
			op->addArgument(arg);
		}
		return getBV(resultID);
	}


	BitVector newMinBV(int resultID, vec<int> & args){

		MinMaxData * op = new MinMaxData(*this,resultID,true);
		addOperation(resultID,op);
		for (int argID:args){
			MinMaxArg * arg = new MinMaxArg(*this,argID,op);
			addOperation(argID,arg);
			op->addArgument(arg);
		}
		return getBV(resultID);
	}

	//ITE with BV arguments
	BitVector newConditionalBV(Lit outer_condition_lit, int bvThenID, int bvElseID,int resultID){
		//need to tie this to to the conditional...

	/*	if(conditionals[resultID].hasCondition()){
			throw std::invalid_argument("Bitvectors can have at most one defined conditional (ITE)");
		}*/
		int bitwidth = getBV(resultID).width();
		if(resultID<=bvThenID || resultID<=bvElseID){
			throw std::invalid_argument("Condition result must have a strictly greater id than its arguments");
		}
		if(bitwidth !=  getBV(bvThenID).width()){
			throw std::invalid_argument("Bit widths must match for bitvectors");
		}
		if(bitwidth !=  getBV(bvElseID).width()){
			throw std::invalid_argument("Bit widths must match for bitvectors");
		}
		Var outerVar = var(outer_condition_lit);
		Lit condition = mkLit(newVar(outerVar, operations.size(),true));
		ConditionalID * op = new ConditionalID(*this, resultID,condition);
		addOperation(resultID,op);


		ConditionalArg * thn = new ConditionalArg(*this, bvThenID,condition,op);
		addOperation(bvThenID,thn);
		ConditionalArg * els = new ConditionalArg(*this, bvElseID,~condition,op);
		addOperation(bvElseID,els);

		op->setThen(thn);
		op->setElse(els);
		thn->setOtherArg(els);
		els->setOtherArg(thn);

		bv_needs_propagation[resultID]=true;
		if(!alteredBV[resultID]){
			alteredBV[resultID]=true;
			altered_bvs.push(resultID);
		}
		bv_needs_propagation[bvThenID]=true;
		if(!alteredBV[bvThenID]){
			alteredBV[bvThenID]=true;
			altered_bvs.push(bvThenID);
		}
		bv_needs_propagation[bvElseID]=true;
		if(!alteredBV[bvElseID]){
			alteredBV[bvElseID]=true;
			altered_bvs.push(bvElseID);
		}
		requiresPropagation=true;

		return getBV(resultID);

	}


	BitVector newSubtractionBV(int resultID, int aID, int bID){
		//n(self.getID(), args[1].getID(), args[0].getID())
		return newAdditionBV(aID,bID,resultID);
	}

	BitVector newAdditionBV(int resultID, int aID, int bID){
		if(!hasBV(aID)){
			throw std::runtime_error("Undefined bitvector ID " + std::to_string(aID));
		}
		if(!hasBV(bID)){
			throw std::runtime_error("Undefined bitvector ID " + std::to_string(bID));
		}
		if(!hasBV(resultID)){
			throw std::runtime_error("Undefined bitvector ID " + std::to_string(resultID));
		}
		while(eq_bitvectors[resultID]!=resultID)
			resultID=eq_bitvectors[resultID];
		while(eq_bitvectors[aID]!=aID)
			aID=eq_bitvectors[aID];
		while(eq_bitvectors[bID]!=bID)
			bID=eq_bitvectors[bID];

		int bitwidth = getBV(resultID).width();
	/*	if(resultID<=aID || resultID<=bID){
			throw std::invalid_argument("Addition result must have a strictly greater id than its arguments");
		}*/
		if(bitwidth !=  getBV(aID).width()){
			throw std::invalid_argument("Bit widths must match for bitvectors");
		}
		if(bitwidth !=  getBV(bID).width()){
			throw std::invalid_argument("Bit widths must match for bitvectors");
		}

		Addition * add = new Addition(*this, resultID);
		addOperation(resultID,add);

		AdditionArg * arg1 = new AdditionArg(*this, aID,add);
		addOperation(aID,arg1);
		AdditionArg * arg2 = new AdditionArg(*this, bID,add);
		addOperation(bID,arg2);

		add->setArg1(arg1);
		add->setArg2(arg2);

		arg1->setOtherArg(arg2);
		arg2->setOtherArg(arg1);

		bv_needs_propagation[resultID]=true;
		if(!alteredBV[resultID]){
			alteredBV[resultID]=true;
			altered_bvs.push(resultID);
		}
		bv_needs_propagation[aID]=true;
		if(!alteredBV[aID]){
			alteredBV[aID]=true;
			altered_bvs.push(aID);
		}
		bv_needs_propagation[bID]=true;
		if(!alteredBV[bID]){
			alteredBV[bID]=true;
			altered_bvs.push(bID);
		}
		requiresPropagation=true;
		return getBV(resultID);
	}
	BitVector newInvertBV(int resultID, int argID){

			if(!hasBV(resultID)){
				throw std::runtime_error("Undefined bitvector ID " + std::to_string(resultID));
			}
			while(eq_bitvectors[resultID]!=resultID)
				resultID=eq_bitvectors[resultID];


			int bitwidth = getBV(resultID).width();


			Invert * op1 = new Invert(*this, resultID);
			addOperation(resultID,op1);
			Invert * op2 = new Invert(*this,  argID);
			addOperation(argID,op2);

			op1->setArg(op2);
			op2->setArg(op1);

			bv_needs_propagation[argID]=true;
			bv_needs_propagation[resultID]=true;
			if(!alteredBV[resultID]){
				alteredBV[resultID]=true;
				altered_bvs.push(resultID);
			}
			if(!alteredBV[argID]){
				alteredBV[argID]=true;
				altered_bvs.push(argID);
			}
			requiresPropagation=true;

			for(int i = 0;i<bitvectors[resultID].size();i++){
				Lit l1 = bitvectors[resultID][i];
				Lit l2 = bitvectors[argID][i];
				//the bits of result, arg must have opposite assignments.
				addClause(l1,l2);
				addClause(~l1,~l2);
			}

			return getBV(resultID);
		}
	BitVector newPopCountBV(int resultID, vec<Var> & args){

		if(!hasBV(resultID)){
			throw std::runtime_error("Undefined bitvector ID " + std::to_string(resultID));
		}
		while(eq_bitvectors[resultID]!=resultID)
			resultID=eq_bitvectors[resultID];


		int bitwidth = getBV(resultID).width();

		if(bitwidth <  ceil(log2(args.size())) ){
			throw std::invalid_argument("Bit width is too small to hold population count");
		}
		PopCountOp * op = new PopCountOp(*this, resultID);
		addOperation(resultID,op);

		for (Var v :args){
			Var a = newVar(v, op->getID(),true);
			op->args.push(a);
			op->over++;
		}

		bv_needs_propagation[resultID]=true;
		if(!alteredBV[resultID]){
			alteredBV[resultID]=true;
			altered_bvs.push(resultID);
		}
		requiresPropagation=true;
		return getBV(resultID);
	}
	void setSymbol(int bvID, const char* symbol){
		symbols[bvID]=symbol;
	}

	const char * getSymbol(int bvID){
		return symbols[bvID];
	}
	int nBitvectors()const{
		return bitvectors.size();
	}
	BitVector newBitvector(int bvID, vec<Var> & vars){
		if(bvID<0){
			bvID = nBitvectors();
		}
		if(bvID==55349){
			int a=1;
		}
		if(vars.size()> (sizeof(Weight) *8)-1){
				throw std::runtime_error("Bit widths larger than " + std::to_string(sizeof(Weight)*8-1) + " are not currently supported (was " + std::to_string(vars.size()) + ")");
			}
		n_bits+=vars.size();
		//bv_callbacks.growTo(id+1,nullptr);
		bitvectors.growTo(bvID+1);
		theoryIds.growTo(bvID+1,-1);
		symbols.growTo(bvID+1,nullptr);
		under_approx.growTo(bvID+1,-1);
		over_approx.growTo(bvID+1,-1);
		under_approx0.growTo(bvID+1,-1);
		over_approx0.growTo(bvID+1,-1);
		alteredBV.growTo(bvID+1);
		bvcompares.growTo(bvID+1);
		compares.growTo(bvID+1);
		bvconst.growTo(bvID+1);
		operation_ids.growTo(bvID+1);

		under_causes.growTo(bvID+1);
		over_causes.growTo(bvID+1);
		cause_set.growTo(bvID+1);
		eq_bitvectors.growTo(bvID+1,-1);
		eq_bitvectors[bvID]=bvID;

		pending_under_analyses.growTo(bvID+1,-1);
		pending_over_analyses.growTo(bvID+1,-1);
		bv_needs_propagation.growTo(bvID+1);
		bv_needs_propagation[bvID]=true;

		if(under_approx[bvID]>-1){
			throw std::invalid_argument("Redefined bitvector ID " + std::to_string(bvID) );
		}
		under_approx[bvID]=0;
		if(vars.size()>0)
			over_approx[bvID]=(1L<<vars.size())-1;
		else
			over_approx[bvID]=0;
		under_approx0[bvID]=under_approx[bvID];
		over_approx0[bvID]=over_approx[bvID];

		BitOp * op = new BitOp(*this, bvID);
		addOperation(bvID,op);
		for(int i = 0;i<vars.size();i++){
			bitvectors[bvID].push(mkLit(newVar(vars[i],op->getID())));
		}




		alteredBV[bvID]=true;
		altered_bvs.push(bvID);
		requiresPropagation=true;
		S->needsPropagation(getTheoryIndex());
		return BitVector(*this,bvID);
	}

	void makeConst(int bvID, Weight c){
		while(eq_bitvectors[bvID]!=bvID)
			bvID=eq_bitvectors[bvID];
		BitVector bv = getBV(bvID);
		for (int i = bv.width()-1;i>=0;i--){
			Weight v = ((Weight)1)<<i;
			Lit l = bv.getBits()[i];
			if(c>=v){
				c-=v;

				addClause(l);
			}else{
				addClause(~l);
			}
		}
	}
	//These are untested!

	//Construct a new bitvector from the concatenation of two bitvectors
	BitVector concat(BitVector a, BitVector b){
		int w = a.width()+b.width();
		BitVector c = newBitvector(-1,w);
		return concat(a,b,c);
	}
	BitVector concat(BitVector a, BitVector b, BitVector result){
		int w = a.width()+b.width();
		assert(result.width()==w);
		//use the SAT solver to make these equal at the bit level.
		//(discarding potentially more precise bounds information for now)
		for(int i = 0;i<a.width();i++){
			this->makeEqual(result[i], a[i]);
		}
		for(int i = 0;i<b.width();i++){
			this->makeEqual(result[i+a.width()], b[i]);
		}
		return result;
	}
	//A BitVector with width upper - lower + 1.
	BitVector slice(BitVector a, int lower, int upper){
		assert(upper>lower);
		assert(lower>=0);
		assert(upper<a.width());
		int w = upper-lower+1;
		assert(w>=0);
		BitVector c = newBitvector(-1,w);
		return slice(a,lower,upper, c);
	}
	BitVector slice(BitVector a, int lower, int upper, BitVector result){
		assert(upper>lower);
		assert(lower>=0);
		assert(upper<a.width());
		upper+=1;
		int w = upper-lower;
		assert(w>=0);
		assert(result.width()==w);
		for(int i = lower;i<upper;i++){
			this->makeEqual(result[i-lower], a[i]);
		}
		return result;
	}

	BitVector bitwiseNot(BitVector a){
		BitVector c = newBitvector(-1,a.width());
		return bitwiseNot(a,c);
	}

	BitVector bitwiseNot(BitVector a, BitVector out){
		assert(a.width()==out.width());
		for(int i = 0;i<a.width();i++){
			this->makeEqual(out[i],~a[i]);
		}
		return out;
	}

	BitVector bitwiseAnd(BitVector a, BitVector b){
		assert(a.width()==b.width());
		BitVector c = newBitvector(-1,a.width());
		return bitwiseAnd(a,b,c);
	}

	BitVector bitwiseAnd(BitVector a, BitVector b, BitVector out){
		assert(a.width()==b.width());
		assert(a.width()==out.width());
		for(int i = 0;i<a.width();i++){
			this->addClause(out[i],~a[i],~b[i]);
			this->addClause(~out[i],a[i]);
			this->addClause(~out[i],b[i]);
		}
		return out;
	}

	BitVector bitwiseNand(BitVector a, BitVector b){
		assert(a.width()==b.width());
		BitVector c = newBitvector(-1,a.width());
		return bitwiseNand(a,b,c);
	}
	BitVector bitwiseNand(BitVector a, BitVector b, BitVector out){
		assert(a.width()==b.width());
		assert(a.width()==out.width());
		for(int i = 0;i<a.width();i++){
			this->addClause(~out[i],~a[i],~b[i]);
			this->addClause(out[i],a[i]);
			this->addClause(out[i],b[i]);
		}
		return out;
	}

	BitVector bitwiseOr(BitVector a, BitVector b){
		assert(a.width()==b.width());
		BitVector c = newBitvector(-1,a.width());
		return bitwiseOr(a,b,c);
	}

	BitVector bitwiseOr(BitVector a, BitVector b, BitVector out){
		assert(a.width()==b.width());
		assert(a.width()==out.width());

		for(int i = 0;i<a.width();i++){
			this->addClause(~out[i],a[i],b[i]);
			this->addClause(out[i],~a[i]);
			this->addClause(out[i],~b[i]);
		}
		return out;
	}

	BitVector bitwiseNor(BitVector a, BitVector b){
		assert(a.width()==b.width());
		BitVector c = newBitvector(-1,a.width());
		return bitwiseNor(a,b,c);
	}

	BitVector bitwiseNor(BitVector a, BitVector b, BitVector out){
		assert(a.width()==b.width());
		assert(a.width()==out.width());
		for(int i = 0;i<a.width();i++){
			this->addClause(out[i],a[i],b[i]);
			this->addClause(~out[i],~a[i]);
			this->addClause(~out[i],~b[i]);
		}
		return out;
	}

	BitVector bitwiseXor(BitVector a, BitVector b){
		assert(a.width()==b.width());
		BitVector c = newBitvector(-1,a.width());
		return bitwiseXor(a,b,c);
	}

	BitVector bitwiseXor(BitVector a, BitVector b, BitVector out){
		assert(a.width()==b.width());
		assert(a.width()==out.width());

		for(int i = 0;i<a.width();i++){
			this->addClause(~a[i],b[i],out[i]);
			this->addClause(a[i],~b[i],out[i]);
			this->addClause(a[i],b[i],~out[i]);
			this->addClause(~a[i],~b[i],~out[i]);
		}
		return out;
	}

	BitVector bitwiseXnor(BitVector a, BitVector b){
		assert(a.width()==b.width());
		BitVector c = newBitvector(-1,a.width());
		return bitwiseXnor(a,b,c);
	}

	BitVector bitwiseXnor(BitVector a, BitVector b, BitVector out){
		assert(a.width()==b.width());

		for(int i = 0;i<a.width();i++){
			this->addClause(~a[i],b[i],~out[i]);
			this->addClause(a[i],~b[i],~out[i]);
			this->addClause(a[i],b[i],out[i]);
			this->addClause(~a[i],~b[i],out[i]);
		}
		return out;
	}

	Lit True(){
		if (const_true==lit_Undef){
			backtrackUntil(0);
			const_true = mkLit(newVar(var(S->True()),-1,true,false),sign(S->True()));
			enqueueTheory(const_true);
		}
		return const_true;
	}
	BitVector newBitvector_Anon(int bvID,int bitwidth){
		if(bvID<0){
			bvID = nBitvectors();
		}
		if(bvID==792){
			int a=1;
		}

		//bv_callbacks.growTo(id+1,nullptr);

		bitvectors.growTo(bvID+1);
		theoryIds.growTo(bvID+1,-1);
		symbols.growTo(bvID+1,nullptr);
		under_approx.growTo(bvID+1,-1);
		over_approx.growTo(bvID+1,-1);
		under_approx0.growTo(bvID+1,-1);
		over_approx0.growTo(bvID+1,-1);
		bvconst.growTo(bvID+1);
		alteredBV.growTo(bvID+1);
		bvcompares.growTo(bvID+1);
		compares.growTo(bvID+1);
		operation_ids.growTo(bvID+1);
		under_causes.growTo(bvID+1);
		over_causes.growTo(bvID+1);
		cause_set.growTo(bvID+1);
		eq_bitvectors.growTo(bvID+1,-1);
		eq_bitvectors[bvID]=bvID;

		pending_under_analyses.growTo(bvID+1,-1);
		pending_over_analyses.growTo(bvID+1,-1);
		bv_needs_propagation.growTo(bvID+1);
		bv_needs_propagation[bvID]=true;
		//bv_callbacks.growTo(bvID+1);
		if(under_approx[bvID]>-1){
			throw std::invalid_argument("Redefined bitvector ID " + std::to_string(bvID) );
		}
		under_approx[bvID]=0;
		over_approx[bvID]=(1L<<bitwidth)-1;
		under_approx0[bvID]=under_approx[bvID];
		over_approx0[bvID]=over_approx[bvID];



		for(int i = 0;i<bitwidth;i++){
			bitvectors[bvID].push(mkLit(newVar(var_Undef,bvID,true,false)));
		}
		BitOp * op = new BitOp(*this, bvID);
		addOperation(bvID,op);



		bv_needs_propagation.growTo(bvID+1);
		bv_needs_propagation[bvID]=true;
		alteredBV[bvID]=true;
		altered_bvs.push(bvID);
		requiresPropagation=true;
		S->needsPropagation(getTheoryIndex());
		return BitVector(*this,bvID);
	}

	BitVector newBitvector(int bvID, int bitwidth,Weight constval=-1, int equivalentBV=-1){
		if(bvID<0){
			bvID = nBitvectors();
		}
		if(bvID==55349){
			int a=1;
		}
		if(bitwidth> (sizeof(Weight) *8)-1){
			throw std::runtime_error("Bit widths larger than " + std::to_string((sizeof(Weight) *8)-1) + " are not currently supported  (was " + std::to_string(bitwidth) + ")");
		}
		if (constval<0)
			n_bits+=bitwidth;
		//bv_callbacks.growTo(id+1,nullptr);

		bitvectors.growTo(bvID+1);
		theoryIds.growTo(bvID+1,-1);
		symbols.growTo(bvID+1,nullptr);
		under_approx.growTo(bvID+1,-1);
		over_approx.growTo(bvID+1,-1);
		under_approx0.growTo(bvID+1,-1);
		over_approx0.growTo(bvID+1,-1);
		bvconst.growTo(bvID+1);
		alteredBV.growTo(bvID+1);
		bvcompares.growTo(bvID+1);
		compares.growTo(bvID+1);
		operation_ids.growTo(bvID+1);
		under_causes.growTo(bvID+1);
		over_causes.growTo(bvID+1);
		cause_set.growTo(bvID+1);
		eq_bitvectors.growTo(bvID+1,-1);
		eq_bitvectors[bvID]=bvID;

		pending_under_analyses.growTo(bvID+1,-1);
		pending_over_analyses.growTo(bvID+1,-1);
		bv_needs_propagation.growTo(bvID+1);
		bv_needs_propagation[bvID]=true;
		//bv_callbacks.growTo(bvID+1);
		if(under_approx[bvID]>-1){
			throw std::invalid_argument("Redefined bitvector ID " + std::to_string(bvID) );
		}
		under_approx[bvID]=0;
		over_approx[bvID]=(1L<<bitwidth)-1;
		under_approx0[bvID]=under_approx[bvID];
		over_approx0[bvID]=over_approx[bvID];

		if(equivalentBV<0){
			if (constval>=0){
				bitvectors[bvID].growTo(bitwidth);
				Weight val = constval;
				//for now, just bitblast this constant value.
				for (int i = bitwidth-1;i>=0;i--){
					Weight v = 1L<<i;
					if (val>=v){
						val-=v;
						bitvectors[bvID][i] = True();
					}else
						bitvectors[bvID][i] = ~True();
				}
				assert(val==0);
			}else{
				for(int i = 0;i<bitwidth;i++){
					bitvectors[bvID].push(mkLit(newVar(var_Undef,bvID)));
				}
			}
			BitOp * op = new BitOp(*this, bvID);
			addOperation(bvID,op);
		}else{
			bitvectors[equivalentBV].copyTo(bitvectors[bvID]);//is this really the right thing to do?
		}

		bv_needs_propagation.growTo(bvID+1);
		bv_needs_propagation[bvID]=true;
		alteredBV[bvID]=true;
		altered_bvs.push(bvID);
		requiresPropagation=true;
		S->needsPropagation(getTheoryIndex());
		return BitVector(*this,bvID);
	}

	bool isConst(int bvID)const{
		return bvconst[bvID];
	}

	bool hasBV(int bvID){
		return bvID>=0 && bvID<under_approx.size() && under_approx[bvID]>-1;
	}
private:
	Lit getComparison(Comparison op, int bvID,const Weight & w){
		//could do a binary search here:
		int cID = getComparisonID(op,bvID,w);

		if(cID<0){
			cID = getComparisonID(-op,bvID,w);//this can be improved upon...
			if(cID<0)
				return lit_Undef;
			else{
				ComparisonOp & cOp = (ComparisonOp &)getOperation(cID);
				return ~cOp.l;
			}
		}else{
			ComparisonOp & cOp = (ComparisonOp &)getOperation(cID);
			return cOp.l;
		}
	}
	int getComparisonID(Comparison op, int bvID,const Weight & w){
		//could do a binary search here:
#ifndef NDEBUG
		int expect = -1;
		for(int i=0;i<compares[bvID].size();i++){
			int cID = compares[bvID][i];
			ComparisonOp & cOp = (ComparisonOp &)getOperation(cID);
			if (cOp.getOp() == op && cOp.w == w){
				expect= cID;
				break;
			}
		}
#endif
		vec<int> & compare = compares[bvID];
		if(compare.size()){
			dbg_compares_sorted(bvID);
			int pos = binary_search_Weight(compare,w);
			if(pos>=0 && pos < compare.size() && ((ComparisonOp &)getOperation(compare[pos])).w==w){
				//we found a comparison to the same bitvector. Now lets see if there is a comparison with the same operator to that bitvector
				while(pos< compare.size() && ((ComparisonOp &)getOperation(compare[pos])).w == w){
					if( ((ComparisonOp &)getOperation(compare[pos])).getOp()==op){
						assert(compare[pos]==expect);
						return compare[pos];
					}
					pos++;
				}
			}
		}
		assert(expect==-1);
		return -1;
	}
	bool dbg_compares_sorted(int bvID){

#ifndef NDEBUG
		vec<int> & compare = compares[bvID];
		for(int i = 1;i<compare.size();i++){
			int cID0 = compare[i-1];
			int cID1 = compare[i];
			assert(cID0 != cID1);
			assert(((ComparisonOp &) getOperation(cID0)).w <= ((ComparisonOp &) getOperation(cID1)).w);
		}
#endif

		return true;
	}
	bool dbg_bvcompares_sorted(int bvID){

#ifndef NDEBUG
		vec<int> & bvcompare = bvcompares[bvID];
		for(int i = 1;i<bvcompare.size();i++){
			int cID0 = bvcompare[i-1];
			int cID1 = bvcompare[i];
			assert(cID0 != cID1);

			assert(((ComparisonBVOp &) getOperation(cID0)).getCompareID() <= ((ComparisonBVOp &) getOperation(cID1)).getCompareID());
		}
#endif

		return true;
	}

	//Returns a CID with the same weight if one exists,  returns the index this CID should be insert after otherwise
	int binary_search_Weight(vec<int> & compares, Weight w)
	{
		int low = 0;
		int high = compares.size() - 1;
		while (low <= high)
		{
			int midpoint = low + (high - low)/2;
			int cID = compares[midpoint];
			ComparisonOp & cOp = ((ComparisonOp &) getOperation(cID));
			if (cOp.w == w )
			{
				//ensure we get the first such matching CID
				while(midpoint>0 && ((ComparisonOp &) getOperation(compares[midpoint-1])).w == w){
					midpoint--;
				}
				return midpoint;
			}
			else if (w <cOp.w)
				high = midpoint - 1;
			else
				low = midpoint + 1;
		}
		//this can probably be done more cleanly...
		if (compares.size()==0)
			return -1;
		else if (low<compares.size() &&  ((ComparisonOp &) getOperation(compares[low])).w > w){
			return low-1;
		}else if (low>compares.size()-1){
			return compares.size()-1;
		}
		return low;
	}

	//Returns a CID comparing to the same bitvector if one exists,  returns the index this CID should be insert after otherwise
	int binary_search_CID(vec<int> & compares,  int compareID)
	{
		int low = 0;
		int high = compares.size() - 1;
		while (low <= high)
		{
			int midpoint = low + (high - low)/2;
			int cID = compares[midpoint];
			ComparisonBVOp & cOp = ((ComparisonBVOp &) getOperation(cID));
			if (cOp.getCompareID() == compareID )
			{
				//ensure we get the first such matching CID
				while(midpoint>0 &&((ComparisonBVOp &) getOperation(compares[midpoint-1])).getCompareID() == compareID){
					midpoint--;
				}
				return midpoint;
			}
			else if (compareID <cOp.getCompareID())
				high = midpoint - 1;
			else
				low = midpoint + 1;
		}
		if (compares.size()==0)
			return -1;
		else if (low<compares.size() && ((ComparisonBVOp &) getOperation(compares[low])).getCompareID()  > compareID){
			return low-1;
		}else if (low>compares.size()-1){
			return compares.size()-1;
		}
		return low;
	}

	Lit getComparisonBV(Comparison op, int bvID,int compareID){

		if(bvID>compareID){
			Lit c =getComparisonBV(-op,compareID,bvID);
			if(c!=lit_Undef){
				return ~c;
			}else{
				return lit_Undef;
			}
		}

		vec<int> & bvcompare = bvcompares[bvID];
		assert(dbg_bvcompares_sorted(bvID));
		//The bv's are kept in sorted order.
		//Do a binary search to find this bvcompare, if it exists.
#ifndef NDEBUG
		Lit expect = lit_Undef;
		for(int i=0;i<bvcompare.size();i++){
			int cID = bvcompare[i];
			ComparisonBVOp & cOp = ((ComparisonBVOp &) getOperation(cID));
			if ( cOp.getCompareID() == compareID && cOp.getOp()==op){
				expect =cOp.l;
				break;
			}
		}
#endif
		if(bvcompare.size()){
			int pos = binary_search_CID(bvcompare,compareID);

			if(pos>=0 && pos < bvcompare.size() && ((ComparisonBVOp &) getOperation(bvcompare[pos])).getCompareID()==compareID){
				//we found a comparison to the same bitvector. Now lets see if there is a comparison with the same operator to that bitvector
				while(pos< bvcompare.size() && ((ComparisonBVOp &) getOperation(bvcompare[pos])).getCompareID() == compareID){
					ComparisonBVOp & cOp = ((ComparisonBVOp &) getOperation(bvcompare[pos]));
					if( cOp.getOp()==op){
						assert(cOp.l==expect);
						return cOp.l;
					}
					pos++;
				}
			}
		}
		assert(expect==lit_Undef);
		return lit_Undef;
	}

public:

	Lit newComparison(Comparison op, int bvID,const Weight & to, Var outerVar = var_Undef, bool decidable=true) {
		Lit l;
		if(!hasBV(bvID)){
			throw std::runtime_error("Undefined bitvector ID " + std::to_string(bvID) );
		}
		while(eq_bitvectors[bvID]!=bvID)
			bvID=eq_bitvectors[bvID];

		if( outerVar == var_Undef){
			//canonicalize the comparison operator to <=
			if(op==Comparison::gt){
				return ~newComparison(Comparison::leq, bvID,to, outerVar, decidable);
			}else if (op == Comparison::geq){
				return ~newComparison(Comparison::leq, bvID,to-1, outerVar, decidable);
			}else if (op == Comparison::lt){
				return newComparison(Comparison::leq, bvID,to-1, outerVar, decidable);
			}
		}
		int comparisonID = operations.size();


		if((l = getComparison(op, bvID, to))!=lit_Undef){
			if(outerVar != var_Undef){
				makeEqualInSolver(mkLit(outerVar),toSolver(l));
			}
			return l;
		}else{
			l = mkLit(newVar(outerVar, comparisonID,true,decidable));

		}
		static int iter = 0;
		if(++iter==29){
			int a =1;
		}

		if(opt_write_learnt_clauses){
			std::stringstream ss;
			ss << op << " " << to ;
			fprintf(opt_write_learnt_clauses,"learnt fact bv %d %s %d\n", unmapBV( bvID), ss.str().c_str(),dimacs(S->unmap(toSolver(l)) ));
			fflush(opt_write_learnt_clauses);
			//std::cout << "learnt fact " << "bv " << bvID <<  <<" " << dimacs(toSolver(l)) << "\n";
		}


		ComparisonOp * compOp = new ComparisonOp(*this, bvID,op,to,l);

		operations.push(compOp);

		dbg_compares_sorted(bvID);
		vec<int> & compare = compares[bvID];
		//insert this value in order.
		int insertPos = binary_search_Weight(compare,to)+1;
		assert(insertPos>=0);assert(insertPos<=compare.size());
		compare.push(comparisonID);
		if(insertPos<compare.size()-1){

			int curVal = compare[insertPos];
			for(int i = insertPos+1;i<compare.size();i++){
				int newVal = compare[i];
				compare[i]=curVal;
				curVal=newVal;
			}
			compare[insertPos]=comparisonID;
		}

		dbg_compares_sorted(bvID);


		comparison_needs_repropagation.growTo(operations.size());
		if(!comparison_needs_repropagation[comparisonID]){
			//we will need to force the associated bitvector to re-update after backtracking.
			comparison_needs_repropagation[comparisonID]=true;
			repropagate_comparisons.push(comparisonID);
		}

		bv_needs_propagation[bvID]=true;
		if(!alteredBV[bvID]){
			alteredBV[bvID]=true;
			altered_bvs.push(bvID);
		}
		//set the value of this immediately, if needed
		requiresPropagation=true;
		S->needsPropagation(getTheoryIndex());

		Weight & underApprox = under_approx[bvID];
		Weight & overApprox = over_approx[bvID];
		Cause  prev_under_cause = under_causes[bvID];
		Cause  prev_over_cause = over_causes[bvID];
		{

			switch (op){
				case Comparison::lt:
					if (overApprox<to){
						if(value(l)==l_True){
							//do nothing
						}else if (value(l)==l_False){
							assert(false);//this should not happen!
						}else {
							assert(value(l)==l_Undef);
							if (over_causes[bvID].type==OperationType::cause_is_decision && overApprox==to-1 ){

								int lev = over_causes[bvID].index;
										assert(lev>=0);
										assert(lev<=decisionLevel());
								if(S->value(toSolver(l))==l_Undef)
									S->instantiateLazyDecision(toSolver(l),lev);
							}else{
								enqueue(l,comparisonprop_marker);
							}

							//enqueueEager(l,bvID,underApprox,overApprox,prev_under_cause,prev_over_cause, comparisonprop_marker);
						}
					}

					if (underApprox>=to){
						if(value(l)==l_True){
							assert(false);
						}else if (value(l)==l_False){
							//do nothing
						}else {
							assert(value(l)==l_Undef);
							//enqueueEager(~l,bvID,underApprox,overApprox,prev_under_cause,prev_over_cause, comparisonprop_marker);

							if (under_causes[bvID].type==OperationType::cause_is_decision && underApprox==to ){
								int lev = under_causes[bvID].index;
								assert(lev>=0);
								assert(lev<=decisionLevel());
								if(S->value(toSolver(l))==l_Undef)
									S->instantiateLazyDecision(toSolver(~l),lev);
							}else
								enqueue(~l,comparisonprop_marker);
						}
					}
					break;
				case Comparison::leq:
					if (overApprox<=to){
						if(value(l)==l_True){
							//do nothing
						}else if (value(l)==l_False){
							assert(false);//this should not happen!
						}else {
							assert(value(l)==l_Undef);
							//enqueueEager(l,bvID,underApprox,overApprox,prev_under_cause,prev_over_cause, comparisonprop_marker);
							//if this was a decision
							if (over_causes[bvID].type==OperationType::cause_is_decision && overApprox==to ){
								int lev = over_causes[bvID].index;
								assert(lev>=0);
								assert(lev<=decisionLevel());
								if(S->value(toSolver(l))==l_Undef)
									S->instantiateLazyDecision(toSolver(l),lev);
							}else
								enqueue(l,comparisonprop_marker);
						}
					}

					if (underApprox>to){
						if(value(l)==l_True){
							assert(false);
						}else if (value(l)==l_False){
							//do nothing
						}else {
							assert(value(l)==l_Undef);
							//enqueueEager(~l,bvID,underApprox,overApprox,prev_under_cause,prev_over_cause, comparisonprop_marker);
							if (under_causes[bvID].type==OperationType::cause_is_decision && underApprox==to+1 ){
								int lev = under_causes[bvID].index;
								assert(lev>=0);
								assert(lev<=decisionLevel());
								if(S->value(toSolver(l))==l_Undef)
									S->instantiateLazyDecision(toSolver(~l),lev);
							}else
								enqueue(~l,comparisonprop_marker);

						}
					}
					break;
				case Comparison::gt:
					if (overApprox<=to){
						if(value(l)==l_True){
							assert(false);
						}else if (value(l)==l_False){

						}else {
							assert(value(l)==l_Undef);
							//enqueueEager(~l,bvID,underApprox,overApprox,prev_under_cause,prev_over_cause, comparisonprop_marker);
							if (over_causes[bvID].type==OperationType::cause_is_decision && overApprox==to ){
								int lev = over_causes[bvID].index;
								assert(lev>=0);
								assert(lev<=decisionLevel());
								if(S->value(toSolver(l))==l_Undef)
									S->instantiateLazyDecision(toSolver(~l),lev);
							}else
								enqueue(~l,comparisonprop_marker);
						}
					}

					if (underApprox>to){
						if(value(l)==l_True){

						}else if (value(l)==l_False){
							assert(false);
						}else {
							assert(value(l)==l_Undef);
							//enqueueEager(l,bvID,underApprox,overApprox,prev_under_cause,prev_over_cause, comparisonprop_marker);
							if (under_causes[bvID].type==OperationType::cause_is_decision && underApprox==to+1 ){
								int lev = under_causes[bvID].index;
								assert(lev>=0);
								assert(lev<=decisionLevel());
								if(S->value(toSolver(l))==l_Undef)
									S->instantiateLazyDecision(toSolver(l),lev);
							}else
								enqueue(l,comparisonprop_marker);

						}
					}
					break;
				case Comparison::geq:
				default:
					if (overApprox<to){
						if(value(l)==l_True){
							assert(false);
						}else if (value(l)==l_False){

						}else {
							assert(value(l)==l_Undef);
							//enqueueEager(~l,bvID,underApprox,overApprox,prev_under_cause,prev_over_cause, comparisonprop_marker);

							if (over_causes[bvID].type==OperationType::cause_is_decision && overApprox==to-1 ){
								int lev = over_causes[bvID].index;
								assert(lev>=0);
								assert(lev<=decisionLevel());
								if(S->value(toSolver(l))==l_Undef)
									S->instantiateLazyDecision(toSolver(~l),lev);
							}else
								enqueue(~l,comparisonprop_marker);
						}
					}

					if (underApprox>=to){
						if(value(l)==l_True){

						}else if (value(l)==l_False){
							assert(false);
						}else {
							assert(value(l)==l_Undef);
							//enqueueEager(l,bvID,underApprox,overApprox,prev_under_cause,prev_over_cause, comparisonprop_marker);
							if (under_causes[bvID].type==OperationType::cause_is_decision && underApprox==to ){
								int lev = under_causes[bvID].index;
								assert(lev>=0);
								assert(lev<=decisionLevel());
								if(S->value(toSolver(l))==l_Undef)
									S->instantiateLazyDecision(toSolver(l),lev);
							}else
								enqueue(l,comparisonprop_marker);
						}
					}
					break;
			}
		}
		return l;
	}



	Lit newComparisonBV(Comparison op, int bvID,int toID, Var outerVar = var_Undef) {

		if(!hasBV(bvID)){
			throw std::runtime_error("Undefined bitvector ID " + std::to_string(toID));
		}
		if(!hasBV(toID)){
			throw std::runtime_error("Undefined bitvector ID " + std::to_string(toID));
		}
		while(eq_bitvectors[bvID]!=bvID)
			bvID=eq_bitvectors[bvID];
		while(eq_bitvectors[toID]!=toID)
			toID=eq_bitvectors[toID];
		if(bvID<toID){

			Lit l = newComparisonBV(~op,toID,bvID,outerVar);//is this correct?
		/*	if(outerVar !=var_Undef){
				makeEqualInSolver(mkLit(outerVar),toSolver(~l));
			}*/
			return l;
		}


		if (bvID==toID){
			if(outerVar==var_Undef){
				//then this is a constant
				if (op==Comparison::leq || op==Comparison::geq){
					return True();
				}else{
					return ~True();
				}
			}else{
				Lit l = mkLit(newVar(outerVar, bvID,-1,false));
				if (op==Comparison::leq || op==Comparison::geq){
					addClause(l);//const true
				}else{
					addClause(~l);//const false
				}
				return l;
			}
		}

		Lit l;
		int comparisonID = operations.size();
		if((l = getComparisonBV(op,bvID, toID))!=lit_Undef){
			if(outerVar != var_Undef){
				makeEqualInSolver(mkLit(outerVar),toSolver(l));
			}
			return l;
		}else{
			l = mkLit(newVar(outerVar, comparisonID));
		}
		ComparisonBVOp * compOp1 = new ComparisonBVOp(*this, bvID,op,l);
		addOperation(bvID,compOp1);
		//operations.push(compOp1);

		dbg_bvcompares_sorted(bvID);


		//insert this value in order.
		//could do a binary search here...


		//insert this value in order.
		{
			vec<int> & compare= bvcompares[bvID];
			int insertPos = binary_search_CID(compare,toID)+1;
			assert(insertPos>=0);assert(insertPos<=compare.size());
			compare.push(comparisonID);
			if(insertPos<compare.size()-1){

				int curVal = compare[insertPos];
				for(int i = insertPos+1;i<compare.size();i++){
					int newVal = compare[i];
					compare[i]=curVal;
					curVal=newVal;
				}
				compare[insertPos]=comparisonID;
			}
		}


/*
		for(int i=0;i<bvcompares[bvID].size()-1;i++){
			int cid = bvcompares[bvID][i];
			if(comparisons[cid].compareID>= toID){
				for(int j = bvcompares[bvID].size()-1; j>i ;j--){
					bvcompares[bvID][j]=bvcompares[bvID][j-1];
				}
				bvcompares[bvID][i]=comparisonID;
				break;
			}
		}
*/

/*		if(!cause_set[bvID].contains(toID))
			cause_set[bvID].push(toID);*/
		//Also need to attach an equivalent (but reversed) comparator to the other bitvector

		comparisonID = operations.size();

		ComparisonBVOp * compOp2 = new ComparisonBVOp(*this, toID,~op,l);
		addOperation(toID,compOp2);
		//operations.push(compOp2);

		compOp1->setOtherOp(compOp2);
		compOp2->setOtherOp(compOp1);
		dbg_bvcompares_sorted(bvID);
		dbg_bvcompares_sorted(toID);
		//insert this value in order.

		{
			vec<int> & compare= bvcompares[toID];
			int insertPos = binary_search_CID(compare,bvID)+1;
			assert(insertPos>=0);assert(insertPos<=compare.size());
			compare.push(comparisonID);
			if(insertPos<compare.size()-1){

				int curVal = compare[insertPos];
				for(int i = insertPos+1;i<compare.size();i++){
					int newVal = compare[i];
					compare[i]=curVal;
					curVal=newVal;
				}
				compare[insertPos]=comparisonID;
			}
		}
		dbg_bvcompares_sorted(toID);
/*
		for(int i=0;i<bvcompares[toID].size()-1;i++){
			int cid = bvcompares[toID][i];
			if(comparisons[cid].compareID>= bvID){
				for(int j = bvcompares[toID].size()-1; j>i ;j--){
					bvcompares[toID][j]=bvcompares[toID][j-1];
				}
				bvcompares[toID][i]=comparisonID;
				break;
			}
		}
*/

		//Add (some) obvious implied relationships to the SAT solver:
		//if we have two relationships,a ==( x>y), b==( x<=y), then we also have (~a or ~b) and (a or b) in the sat solver.
		//If we have relationships a ==x>y, b== x>=y, then we have (a -> b).
		//if we have relationships a == x < y, b== y < z, c== x < z, then we have ((a and b) -> c)
		//more?
		comparison_needs_repropagation.growTo(operations.size());
		if(!comparison_needs_repropagation[comparisonID-1]){
			//we will need to force the associated bitvector to re-update after backtracking.
			comparison_needs_repropagation[comparisonID-1]=true;
			repropagate_comparisons.push(comparisonID-1);
		}
		if(!comparison_needs_repropagation[comparisonID]){
			// we will need to force the associated bitvector to re-update after backtracking.
			comparison_needs_repropagation[comparisonID]=true;
			repropagate_comparisons.push(comparisonID);
		}


		bv_needs_propagation[bvID]=true;
		if(!alteredBV[bvID]){
			alteredBV[bvID]=true;
			altered_bvs.push(bvID);
		}

		bv_needs_propagation[toID]=true;
		if(!alteredBV[toID]){
			alteredBV[toID]=true;
			altered_bvs.push(toID);
		}

		//set the value of this immediately, if needed
		requiresPropagation=true;
		S->needsPropagation(getTheoryIndex());
		Weight & underApprox = under_approx[bvID];
		Weight & overApprox = over_approx[bvID];

		Weight & underCompare = under_approx[toID];
		Weight & overCompare = over_approx[toID];
		Cause  prev_under_cause = under_causes[bvID];
		Cause  prev_over_cause = over_causes[bvID];
		switch (op){
			case Comparison::lt:
				if (overApprox<underCompare){
					if(value(l)==l_True){
						//do nothing
					}else if (value(l)==l_False){
						assert(false);//this should not happen!
					}else {
						assert(value(l)==l_Undef);

						enqueueEager(l,bvID,underApprox,overApprox,prev_under_cause,prev_over_cause, comparisonbv_prop_marker);
					}
				}

				if (underApprox>=overCompare){
					if(value(l)==l_True){
						assert(false);
					}else if (value(l)==l_False){
						//do nothing
					}else {
						assert(value(l)==l_Undef);
						enqueueEager(~l,bvID,underApprox,overApprox,prev_under_cause,prev_over_cause, comparisonbv_prop_marker);
					}
				}
				break;
			case Comparison::leq:
				if (overApprox<=underCompare){
					if(value(l)==l_True){
						//do nothing
					}else if (value(l)==l_False){
						assert(false);//this should not happen!
					}else {
						assert(value(l)==l_Undef);
						enqueueEager(l,bvID,underApprox,overApprox,prev_under_cause,prev_over_cause, comparisonbv_prop_marker);
					}
				}

				if (underApprox>overCompare){
					if(value(l)==l_True){
						assert(false);
					}else if (value(l)==l_False){
						//do nothing
					}else {
						assert(value(l)==l_Undef);
						enqueueEager(~l,bvID,underApprox,overApprox,prev_under_cause,prev_over_cause, comparisonbv_prop_marker);
					}
				}
				break;
			case Comparison::gt:
				if (overApprox<=underCompare){
					if(value(l)==l_True){
						assert(false);
					}else if (value(l)==l_False){

					}else {
						assert(value(l)==l_Undef);
						enqueueEager(~l,bvID,underApprox,overApprox,prev_under_cause,prev_over_cause, comparisonbv_prop_marker);
					}
				}

				if (underApprox>overCompare){
					if(value(l)==l_True){

					}else if (value(l)==l_False){
						assert(false);
					}else {
						assert(value(l)==l_Undef);
						enqueueEager(l,bvID,underApprox,overApprox,prev_under_cause,prev_over_cause, comparisonbv_prop_marker);
					}
				}
				break;
			case Comparison::geq:
			default:
				if (overApprox<underCompare){
					if(value(l)==l_True){
						assert(false);
					}else if (value(l)==l_False){

					}else {
						assert(value(l)==l_Undef);
						enqueueEager(~l,bvID,underApprox,overApprox,prev_under_cause,prev_over_cause, comparisonbv_prop_marker);
					}
				}

				if (underApprox>=overCompare){
					if(value(l)==l_True){

					}else if (value(l)==l_False){
						assert(false);
					}else {
						assert(value(l)==l_Undef);
						enqueueEager(l,bvID,underApprox,overApprox,prev_under_cause,prev_over_cause, comparisonbv_prop_marker);
					}
				}
				break;
		}

		return l;
	}

	void printSolution() {

	}

	bool dbg_uptodate(){
#ifndef NDEBUG
		//dbg_synced();
		for(int bvID = 0;bvID<bitvectors.size();bvID++){
			Weight under;
			Weight over;
			assert(checkApproxUpToDate(bvID,&under,&over));

			//Weight & under = under_approx[bvID];
			//Weight & over = over_approx[bvID];
			for(int cID:compares[bvID]){
				if(comparison_needs_repropagation[cID])
					continue;//this is a newly added comparison which might not have been propagated yet
			ComparisonOp & c =(ComparisonOp & ) getOperation(cID);
			Comparison op = c.getOp();

			lbool val = value(c.l);

			switch(op){
				case Comparison::lt:
					if(over<c.w){
						assert(value(c.l)==l_True);
					}else if (under>=c.w){
						assert(value(c.l)==l_False);
					}

					break;
				case Comparison::leq:
					if(over<=c.w){
						assert(value(c.l)==l_True);
					}else if (under>c.w){
						assert(value(c.l)==l_False);
					}
					break;
				case Comparison::gt:
					if(over<=c.w){
						assert(value(c.l)==l_False);
					}else if (under>c.w){
						assert(value(c.l)==l_True);
					}
					break;
				case Comparison::geq:
				default:
					if(over<c.w){
						assert(value(c.l)==l_False);
					}else if (under>=c.w){
						assert(value(c.l)==l_True);
					}
				break;
			}

		}
		for(int cID:bvcompares[bvID]){
			ComparisonBVOp & c =(ComparisonBVOp & ) getOperation(cID);
			Comparison op = c.getOp();

			lbool val = value(c.l);

			Weight under_w = under_approx[c.getCompareID()];
			Weight over_w = over_approx[c.getCompareID()];

			switch(op){

				case Comparison::lt:
					if( over<under_w){
						assert(value(c.l)==l_True);
					}else if (under>=over_w){
						assert(value(c.l)==l_False);
					}
					break;
				case Comparison::leq:
					if( over<=under_w){
						assert(value(c.l)==l_True);
					}else if (under>over_w){
						assert(value(c.l)==l_False);
					}
					break;
				case Comparison::gt:
					if( over<=under_w){
						assert(value(c.l)==l_False);
					}else if (under>over_w){
						assert(value(c.l)==l_True);
					}
					break;
				case Comparison::geq:
				default:
					if( over<under_w){
						assert(value(c.l)==l_False);
					}else if (under>=over_w){
						assert(value(c.l)==l_True);
					}
					break;
			}

		}
		}
#endif
		return true;
	}

public:
	//If bitvectors were renumbered during parsing, this obtains the original numbering
	int unmapBV(int bvID){
		if(bvRemap){
			return bvRemap->unmapBV(bvID);
		}else{
			return bvID;
		}
	}
	void setBVMap(BVMap * map){
		bvRemap=map;
	}

	bool applyOp(Comparison op,int bvID,Weight to){
		switch(op){
			case Comparison::geq:
				return under_approx[bvID]>=to;
			case Comparison::gt:
				return under_approx[bvID]>to;
			case Comparison::leq:
				return over_approx[bvID]<=to;
			case Comparison::lt:
				return over_approx[bvID]<to;
			default:
				assert(false);
				return false;
		}

	}

};


template<typename Weight>
inline Weight BVTheorySolver<Weight>::ceildiv(Weight a, Weight b){
	Weight rem = a%b;
	Weight div = a/b;
	return  (rem !=0) ? (div+1) : div;
}
template<>
inline mpq_class BVTheorySolver<mpq_class>::ceildiv(mpq_class a, mpq_class b){
	//what is the correct way to do this, if any?
	assert(false);//this is not properly implemented yet...
	return a/b;

}
template<>
inline mpf_class BVTheorySolver<mpf_class>::ceildiv(mpf_class a, mpf_class b){
	return ceil(a/b);
}
template<>
inline double BVTheorySolver<double>::ceildiv(double a, double b){
	return std::ceil(a/b);
}

template<>
inline float BVTheorySolver<float>::ceildiv(float a, float b){
	return std::ceil(a/b);
}
template<>
inline mpq_class BVTheorySolver<mpq_class>::lowest(int bvID){
	throw std::runtime_error("Unimplemented");
}
template<>
inline mpq_class BVTheorySolver<mpq_class>::highest(int bvID){
	throw std::runtime_error("Unimplemented");
}
template<>
inline mpq_class BVTheorySolver<mpq_class>::refine_ubound(int bvID, mpq_class bound, Var ignore_bit){
	throw std::runtime_error("Unimplemented");
}
template<>
inline double BVTheorySolver<double>::lowest(int bvID){
	throw std::runtime_error("Unimplemented");
}
template<>
inline double BVTheorySolver<double>::highest(int bvID){
	throw std::runtime_error("Unimplemented");
}
template<>
inline double BVTheorySolver<double>::refine_ubound(int bvID, double bound, Var ignore_bit){
	throw std::runtime_error("Unimplemented");
}
template<>
inline mpq_class BVTheorySolver<mpq_class>::refine_lbound(int bvID, mpq_class bound, Var ignore_bit){
	throw std::runtime_error("Unimplemented");
}
template<>
inline double BVTheorySolver<double>::refine_lbound(int bvID, double bound, Var ignore_bit){
	throw std::runtime_error("Unimplemented");
}
template<typename Weight>
Weight BVTheorySolver<Weight>::lowest(int bvID){
	    Weight val= 0;
	    vec<Lit> & bv = bitvectors[bvID];
	    for (int i = 0;i<bv.size();i++){
	        if (value(bv[i])==l_True)
	            val += ((Weight )1)<<i;
	    }
	    return val;
	}
template<typename Weight>
	Weight BVTheorySolver<Weight>::highest(int bvID){
	    Weight val= 0;
	    vec<Lit> & bv = bitvectors[bvID];
	    for (int i = 0;i<bv.size();i++){
	        if (value(bv[i])!=l_False)
	            val += ((Weight )1)<<i;
	    }
	    return val;
	}

template<>
inline mpq_class BVTheorySolver<mpq_class>::refine_ubound_check(int bvID, mpq_class bound, Var ignore_bit){
	throw std::runtime_error("Unimplemented");
}
template<>
inline double BVTheorySolver<double>::refine_ubound_check(int bvID, double bound, Var ignore_bit){
	throw std::runtime_error("Unimplemented");
}


template<>
inline mpq_class BVTheorySolver<mpq_class>::refine_lbound_check(int bvID, mpq_class bound, Var ignore_bit){
	throw std::runtime_error("Unimplemented");
}
template<>
inline double BVTheorySolver<double>::refine_lbound_check(int bvID, double bound, Var ignore_bit){
	throw std::runtime_error("Unimplemented");
}

template<>
inline void BVTheorySolver<double>::dbg_evaluate(int bvID, int pos,vec<double> & vals,double val){

}
template<>
inline void BVTheorySolver<mpq_class>::dbg_evaluate(int bvID, int pos,vec<mpq_class> & vals,mpq_class val){

}

template<typename Weight>
void BVTheorySolver<Weight>::dbg_evaluate(int bvID, int pos,vec<Weight> & vals,Weight val){
    if(pos==-1){
        vals.push(val);
        return;
    }
    vec<Lit> & bv = bitvectors[bvID];
    if (value(bv[pos]) ==l_True or value(bv[pos]) ==l_Undef){
    	Weight bit = (1<<pos);
		dbg_evaluate(bvID,pos-1,vals, val+bit);
    }
    if (value(bv[pos]) ==l_False or value(bv[pos]) ==l_Undef)
    	dbg_evaluate(bvID,pos-1,vals, val);
}
template<typename Weight>
Weight BVTheorySolver<Weight>::refine_ubound_check(int bvID, Weight bound, Var ignore_bit){
#ifndef NDEBUG
	//test all values of mbits, find the lowest one >= i
	static vec<Weight> vals;
	vals.clear();

    dbg_evaluate(bvID,bitvectors[bvID].size()-1,vals,0);
    Weight lowest_bound = -1;
    for(Weight x:vals){
    	if (x>= bound && (lowest_bound<0 || x<lowest_bound)){
    		lowest_bound=x;
    	}
    }
    return lowest_bound;

#endif
    return 0;
}

template<typename Weight>
Weight BVTheorySolver<Weight>::refine_lbound_check(int bvID, Weight bound, Var ignore_bit){
#ifndef NDEBUG
	//test all values of mbits, find the lowest one >= i
	static vec<Weight> vals;
	vals.clear();

    dbg_evaluate(bvID,bitvectors[bvID].size()-1,vals,0);
    Weight highest_bound = -1;
    for(Weight x:vals){
    	if (x<= bound && (highest_bound<0 || x>highest_bound)){
    		highest_bound=x;
    	}
    }
    return highest_bound;

#endif
    return 0;
}

	template<typename Weight>
	Weight BVTheorySolver<Weight>::refine_ubound(int bvID, Weight bound, Var ignore_bit){
#ifndef NDEBUG
		Weight expected = refine_ubound_check(bvID,bound,ignore_bit);
#endif
		//Weight under_old = under_approx[bvID];
		//Weight over_old = over_approx[bvID];
		vec<Lit> & bv = bitvectors[bvID];
		bool done = false;
		//under_causes[bvID].clear();
		//over_causes[bvID].clear();
		int last_set_x =-1;
		int j = bv.size()-1;


		Weight proposed_bound=0;
		for(int i = bv.size()-1; i>=0 && ! done;i--){

			Weight ibit = ((Weight )1)<<i;
			if ((value(bv[i])!=l_True || var(bv[i])==ignore_bit) && (bound & (ibit))){
				bool found=false;
				while(j>=i){
					Weight jbit = ((Weight )1)<<j;
					if ((value(bv[j])==l_Undef || var(bv[j])==ignore_bit) && !(proposed_bound & jbit) ){
						//bv[j] must also not have been used yet as a flip value for this case (else we should try the (!(bound & (jbit))) case)...
						last_set_x=j;
					}else if  (!(bound & (jbit)) && ((value(bv[j])==l_True || var(bv[j])==ignore_bit || (proposed_bound & jbit)))){
						found=true;
						last_set_x=-1;
						//if j is > than i, then we are done.
						if (j>i){
							done=true;
							break;
						}
					}
					j--;
				}
				if (last_set_x  >-1){
					found=true;
					//#we were forced to set an xbit to 1 to match this ibit.
					Weight xbit = ((Weight )1)<<last_set_x;
					assert(!(proposed_bound & xbit));
					proposed_bound|=xbit;
					//mp[last_set_x]='s'
					if (last_set_x>i){
						//#mp is > i, so we are done
						done=true;
						break;
					}else{
						assert(last_set_x==i);
						j=bv.size()-1;//can we do better than this?
					}
					last_set_x=-1;
				}

				if (! found){
					assert(expected==-1);
					return -1;
				}
			}
		}
		Weight refined_bound=lowest(bvID);
	    j =bv.size()-1;
		while (j>=0){
			Weight jbit = ((Weight )1)<<j;
			if (proposed_bound & jbit){
				assert(!(refined_bound & jbit));
				refined_bound |=jbit;
				if (!(bound & (jbit))){
					break;
				}

			}
			j--;
		}


	    if(refined_bound< bound){
	    	assert(expected==-1);
	        return -1;
	    }
	        //#return ubound(ibits,mbits2)
#ifndef NDEBUG
		if(refined_bound!=expected){
			assert(false);
			exit(5);
		}
#endif
	    return refined_bound;
	}
template<typename Weight>
Weight BVTheorySolver<Weight>::refine_lbound(int bvID, Weight obound, Var ignore_bit){
#ifndef NDEBUG
	Weight expected = refine_lbound_check(bvID,obound,ignore_bit);
#endif
/*	static int iter = 0;
	if(++iter==22739){
		int a=1;
	}*/
	//Weight under_old = under_approx[bvID];
	//Weight over_old = over_approx[bvID];
	vec<Lit> & bv = bitvectors[bvID];
	bool done = false;
	//under_causes[bvID].clear();
	//over_causes[bvID].clear();

	int j = bv.size()-1;
	int last_set_x =-1;
	Weight bound = ~obound;

	Weight proposed_bound=0;
	for(int i = bv.size()-1; i>=0 && ! done;i--){

		Weight ibit = ((Weight )1)<<i;
		if ((value(bv[i])!=l_False || var(bv[i])==ignore_bit) && (bound & (ibit))){
			bool found=false;
			while(j>=i){
				Weight jbit = ((Weight )1)<<j;
				if ((value(bv[j])==l_Undef || var(bv[j])==ignore_bit) && !(proposed_bound & jbit) ){
					//bv[j] must also not have been used yet as a flip value for this case (else we should try the (!(bound & (jbit))) case)...
					last_set_x=j;
				}else if  (!(bound & (jbit)) && ((value(bv[j])==l_False || var(bv[j])==ignore_bit || (proposed_bound & jbit)))){
					found=true;
					last_set_x=-1;
					//if j is > than i, then we are done.
					if (j>i){
						done=true;
						break;
					}
				}
				j--;
			}
			if (last_set_x  >-1){
				found=true;
				//#we were forced to set an xbit to 1 to match this ibit.
				Weight xbit = ((Weight )1)<<last_set_x;
				assert(!(proposed_bound & xbit));
				proposed_bound|=xbit;
				//mp[last_set_x]='s'
				if (last_set_x>i){
					//#mp is > i, so we are done
					done=true;
					break;
				}else{
					assert(last_set_x==i);
					j=bv.size()-1;//can we do better than this?
				}
				last_set_x =-1;
			}

			if (! found){
				assert(expected==-1);
				return -1;
			}
		}
	}
	Weight refined_bound= 0;//(bvID);
	for (int i = 0;i<bv.size();i++){
		if (value(bv[i])==l_False)
			refined_bound += ((Weight )1)<<i;
	}
	j =bv.size()-1;
	while (j>=0){
		Weight jbit = ((Weight )1)<<j;
		if (proposed_bound & jbit){
			assert(!(refined_bound & jbit));
			refined_bound |=jbit;
			if (!(bound & (jbit))){
				break;
			}

		}
		j--;
	}

	refined_bound = (1L<<bv.size()) + ~refined_bound;
	if(refined_bound> obound){
		assert(expected==-1);
		return -1;
	}

		//#return ubound(ibits,mbits2)

#ifndef NDEBUG
		if(refined_bound!=expected){
			assert(false);
			exit(5);
		}
#endif
	return refined_bound;
}
template<typename Weight>
using BitVector = typename BVTheorySolver<Weight>::BitVector;

}
;

#endif
