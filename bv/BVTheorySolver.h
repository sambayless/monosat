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

	struct Cause{
		int cause_is_bits:1;
		int refined_cause:1;
		int cause_is_comparison:1;
		int cause_is_addition:1;
		int cause_is_addition_argument:1;
		int cause_is_condition:1;
		int cause_is_condition_argument:1;
		int cause_is_decision:1;
		int index:26;

		Cause(const Cause & copy):cause_is_bits(copy.cause_is_bits),refined_cause(copy.refined_cause),cause_is_comparison(copy.cause_is_comparison),cause_is_addition(copy.cause_is_addition),cause_is_addition_argument(copy.cause_is_addition_argument),cause_is_condition(copy.cause_is_condition),cause_is_condition_argument(copy.cause_is_condition_argument),cause_is_decision(copy.cause_is_decision),index(copy.index){

		}

		Cause():cause_is_bits(0),refined_cause(0),cause_is_comparison(0),cause_is_addition(0),cause_is_addition_argument(0),cause_is_condition(0),cause_is_condition_argument(0),cause_is_decision(0),index(-1){

		}
		/*Cause(bool bits, bool addition, int comparison=-1):cause_is_bits(bits),refined_cause(0),cause_is_addition(addition),cause_is_addition_arg(0),comparison_cause(comparison){

		}*/
		bool hasCause(){
			//this can be improved, if we want.
			return cause_is_bits|| cause_is_addition || cause_is_addition_argument|| cause_is_condition || cause_is_condition_argument|| cause_is_comparison || cause_is_decision|| refined_cause;
		}
		void clear(){
			cause_is_decision=0;
			cause_is_bits=0;
			refined_cause=0;
			cause_is_addition=0;
			cause_is_addition_argument=0;
			cause_is_condition=0;
			cause_is_condition_argument=0;
			cause_is_comparison=0;
			index=-1;
		}
	};
	//An assignment is either an assignment to a variable, or is a tightening of a bitvector's bounds
	struct Assignment {

		bool isOperation:1;
		bool isComparator :1;
		bool assign :1;
		int bvID:29;
		Var var;
		Weight previous_under;
		Weight previous_over;
		Weight new_under;
		Weight new_over;

		Cause prev_under_cause;
		Cause prev_over_cause;
		Cause new_under_cause;
		Cause new_over_cause;
		Assignment( int bvID,Weight previous_under,Weight previous_over,Weight new_under,Weight new_over,
				Cause prev_under_cause,Cause prev_over_cause,Cause new_under_cause,Cause new_over_cause):isOperation(false),isComparator(false),assign(0),bvID(bvID),var(var_Undef),previous_under(previous_under),previous_over(previous_over),new_under(new_under),new_over(new_over),
				prev_under_cause(prev_under_cause),prev_over_cause(prev_over_cause),new_under_cause(new_under_cause),new_over_cause(new_over_cause){

		}
		Assignment(bool isComparator, bool assign, int bvID, Var v):isOperation(false),isComparator(isComparator),assign(assign),bvID(bvID),var(v),previous_under(-1),previous_over(-1),new_under(-1),new_over(-1),
				prev_under_cause(),prev_over_cause(),new_under_cause(),new_over_cause(){
		}
		Assignment(bool isComparator, bool assign, int bvID, Var v, bool isOperation):isOperation(isOperation),isComparator(isComparator),assign(assign),bvID(bvID),var(v),previous_under(-1),previous_over(-1),new_under(-1),new_over(-1),
				prev_under_cause(),prev_over_cause(),new_under_cause(),new_over_cause(){
		}
		bool isBoundAssignment()const{
			if(new_under>-1){
				assert(!isOperation);
				assert(!isComparator);
			}
			return new_under>-1;
		}
	};


	struct ComparisonID{
		Weight w;
		Lit l;
		int compareID;
		int bvID:30;
		int is_lt:1;
		int is_strict:1;

		ComparisonID(Weight w, int compareID, Lit l,int bvID, Comparison op ):w(w),l(l),compareID(compareID),bvID(bvID),is_lt(op==Comparison::lt || op==Comparison::leq),is_strict(op==Comparison::lt || op==Comparison::gt){

		}

		bool bvCompare(){
			return compareID>-1;
		}

		Comparison op(){
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
	};


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

		int getID(){
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

public:

	bool first_propagation=true;
	long n_bits =0;
	long n_consts = 0;
	long n_starting_consts=0;
	long n_additions=0;
	vec<lbool> assigns;
	CRef comparisonprop_marker;
	CRef conditionthen_prop_marker;
	CRef conditionelse_prop_marker;
	CRef conditionarg_prop_marker;
	CRef bvprop_marker;
	Lit const_true=lit_Undef;
	vec<const char*> symbols;

	vec<Assignment> trail;
	vec<int> trail_lim;

	vec<Cause> under_causes;
	vec<Cause> over_causes;
	//Bitvectors are unsigned, and have least significant bit at index 0
	vec<vec<Lit> > bitvectors;


	struct ConditionalID{
		Lit l=lit_Undef;
		int cID=-1;
		int bvThenID=-1;
		int bvElseID=-1;
		int bvResultID=-1;
		bool hasITE()const{
			return l!=lit_Undef;
		}
	};

	struct ConditionalArg{
		Lit l=lit_Undef;
		int cID=-1;
		int other_argID=-1;
		int resultID=-1;
	};

	vec<ConditionalID> all_conditionals;
	vec<vec<ConditionalID>> conditionals;//ITE bvs
	vec<vec<ConditionalArg>> conditional_arguments;

	vec<ComparisonID> comparisons;
	vec<vec<int> > cause_set;//for each bv, this is the set of all bitvectors that have a greater id and might need to have their approx updated when this bv's approx changes.
	vec<vec<int> > compares; //for each bitvector, comparisons are all to unique values, and in ascending order of compareTo.
	vec<vec<int>> bvcompares;
	vec<int> eq_bitvectors;//if a bv has been proven to be equivalent to another, lower index bv, put the lowest such index here.
	vec<bool> bv_needs_propagation;
	vec<bool> comparison_needs_repropagation;
	vec<int> repropagate_comparisons;

	struct Addition{
		int aID=-1;
		int bID=-1;

		bool hasAddition()const{
			return aID>=0;
		}

	};
	struct AdditionArg{
		int other_argID=-1;
		int sumID=-1;
	};

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

	vec<vec<Addition>> additions;//each bitvector is the result of at most one addition.
	vec<vec<AdditionArg>> addition_arguments;
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
		int isOperation:1;
		int bvID:29;
		int operationID;
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
		rnd_seed = opt_random_seed;
		S->addTheory(this);
		comparisonprop_marker = S->newReasonMarker(this);
		conditionthen_prop_marker = S->newReasonMarker(this);
		conditionelse_prop_marker = S->newReasonMarker(this);
		conditionarg_prop_marker= S->newReasonMarker(this);
		bvprop_marker = S->newReasonMarker(this);

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
				write_to<<"bv"<<bvID << " = " << under << "\n";

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
		for(int i = 0;i<conditionals[bvID1].size();i++){
			conditionals[bvID2].push(conditionals[bvID1][i]);
			//fix the other bitvector's references.
			int arg1 = conditionals[bvID2][i].bvThenID;
			int arg2 = conditionals[bvID2][i].bvElseID;

			for(int j = 0;j<conditional_arguments[arg1].size();j++){
				if(conditional_arguments[arg1][j].resultID==bvID1){
					conditional_arguments[arg1][j].resultID=bvID2;
				}
			}
			for(int j = 0;j<conditional_arguments[arg2].size();j++){
				if(conditional_arguments[arg2][j].resultID==bvID1){
					conditional_arguments[arg2][j].resultID=bvID2;
				}
			}
		}

		//move the arguments of bvID1 into bvID2.
		for(int i = 0;i<additions[bvID1].size();i++){
			additions[bvID2].push(additions[bvID1][i]);
			//fix the other bitvector's references.
			int arg1 = additions[bvID2][i].aID;
			int arg2 = additions[bvID2][i].bID;

			for(int j = 0;j<addition_arguments[arg1].size();j++){
				if(addition_arguments[arg1][j].sumID==bvID1){
					addition_arguments[arg1][j].sumID=bvID2;
				}
			}
			for(int j = 0;j<addition_arguments[arg2].size();j++){
				if(addition_arguments[arg2][j].sumID==bvID1){
					addition_arguments[arg2][j].sumID=bvID2;
				}
			}
		}
		//additions[bvID1].clear(); //don't clear this, because I want to be able to double check it, later...
		for(int i = 0;i<addition_arguments[bvID1].size();i++){
			addition_arguments[bvID2].push(addition_arguments[bvID1][i]);
			int sumID = addition_arguments[bvID1][i].sumID;
			int other_argID = addition_arguments[bvID1][i].other_argID;

			for(int j = 0;j<addition_arguments[other_argID].size();j++){
				if(addition_arguments[other_argID][j].other_argID==bvID1){
					addition_arguments[other_argID][j].other_argID=bvID2;
				}
			}
			for(int j = 0;j<additions[sumID].size();j++){
				if(additions[sumID][j].aID==bvID1){
					additions[sumID][j].aID=bvID2;
				}
				if(additions[sumID][j].bID==bvID1){
					additions[sumID][j].bID=bvID2;
				}
			}
		}

		for(int i = compares[bvID1].size()-1;i>=0;i--){
			int cID = compares[bvID1][i];
			comparisons[cID].bvID=bvID2;
			compares[bvID2].push(cID);
		}

		for(int i = bvcompares[bvID1].size()-1;i>=0;i--){
			int cID = bvcompares[bvID1][i];
			comparisons[cID].bvID=bvID2;
			int compareID = comparisons[cID].compareID;
			//need to update that bv's comparison also
			for(int j = 0;j<bvcompares[compareID].size();j++){
				int cID2 = bvcompares[compareID][j];
				if(comparisons[cID2].compareID==bvID1){
					comparisons[cID2].compareID=bvID2;
				}
			}
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
	inline bool isBVVar(Var v)const{
		return vars[v].bvID>=0;
	}

	inline bool isOperationVar(Var v) const{
		assert(v < vars.size());
		return vars[v].operationID>-1  && vars[v].isOperation;
	}
	inline int getOperationID(Var v) const {
		assert(isOperationVar(v));
		return vars[v].operationID;
	}

	inline bool isComparisonVar(Var v) const{
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
	}

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
	

	Var newVar(Var solverVar=var_Undef, int bvID=-1, int comparisonID=-1, bool connectToTheory = true, bool decidable=true,bool isOperation=false) {
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
		vars[v].isOperation=isOperation;
		vars[v].bvID=bvID;
		vars[v].operationID=comparisonID;
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
				}else{
					Var x = e.var;
					Lit p = mkLit(x,!e.assign);
					assert(value(x)==l_Undef);
					assigns[x] = sign(p) ? l_False : l_True;
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
				}else{
					Var x = e.var;
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
			}else{
				Var x = e.var;
				Lit p = mkLit(x,!e.assign);
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
			}else{
				Var x = e.var;
				if(x==until_assign){
					break;
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
					if (e.isComparator) {
						int cID = e.bvID;
						assert(cID>=0);
						if (comparison_needs_repropagation[cID]){
							ComparisonID & c = comparisons[cID];
							int bvID = c.bvID;
							if(!alteredBV[bvID]){
								alteredBV[bvID]=true;
								altered_bvs.push(bvID);
							}
							requiresPropagation=true;
							S->needsPropagation(getTheoryIndex());
						}
					}else if (e.isOperation){
						int cID = e.bvID;
						assert(cID>=0);
						//does anything need to happen here?
					}

					assigns[e.var] = l_Undef;
				}
				//changed = true;
			}
			trail.shrink(trail.size() - stop);
			trail_lim.shrink(trail_lim.size() - level);
			assert(trail_lim.size() == level);
			assert(dbg_uptodate());
			if(level==0){//decisionLevel()==0 This check can fail if the levels of the theory and sat solver are out of sync
				for(int cID:repropagate_comparisons){
					assert(comparison_needs_repropagation[cID]);
					comparison_needs_repropagation[cID]=false;
					//need to handle the case where we backtracked before ever propagating this comparison's bitvector...
					ComparisonID & c = comparisons[cID];
					int bvID = c.bvID;
					if(!alteredBV[bvID]){
						alteredBV[bvID]=true;
						altered_bvs.push(bvID);
					}
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
					under_causes[bvID].cause_is_decision=true;
					under_causes[bvID].index=decisionLevel();
					break;
				case Comparison::geq:
					under_approx[bvID]=to;
					under_causes[bvID].clear();
					under_causes[bvID].cause_is_decision=true;
					under_causes[bvID].index=decisionLevel();
					break;
				case Comparison::lt:
					over_approx[bvID]=to-1;
					over_causes[bvID].clear();
					over_causes[bvID].cause_is_decision=true;
					under_causes[bvID].index=decisionLevel();
					break;
				case Comparison::leq:
					over_approx[bvID]=to;
					over_causes[bvID].clear();
					over_causes[bvID].cause_is_decision=true;
					under_causes[bvID].index=decisionLevel();
				break;
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
	



	void newDecisionLevel() {
		trail_lim.push(trail.size());
	}
	;

	void buildReason(Lit p, vec<Lit> & reason,CRef marker) {
		static int iter = 0;
		if(++iter==164){//17
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



		if (marker == comparisonprop_marker) {
			reason.push(p);
			Var v = var(p);
			int bvID = getbvID(v);

			assert(isComparisonVar(v));
			int comparisonID = getComparisonID(v);

			Comparison op = comparisons[comparisonID].op();
			if (comparisons[comparisonID].compareID<0){
				if(sign(p)){
					op=-op;
				}
				buildComparisonReason(op,bvID,comparisons[comparisonID].w,reason);
			}else{
				int compareBV = comparisons[comparisonID].compareID;
				if(sign(p)){
					op=-op;
				}
				buildComparisonReasonBV(op,bvID,compareBV,reason);
			}


		} else if (marker == bvprop_marker) {
			Var v = var(p);

			assert(!isComparisonVar(v));
			int bvID = getbvID(v);

			reason.push(p);
			rewind_trail_pos(analysis_trail_pos-1);
			assert(value(p)==l_Undef);
			Weight underApprox = under_approx[bvID];
			Weight overApprox = over_approx[bvID];


			assert(under_approx>=0); assert(overApprox>=0);
			vec<Lit> & bv = bitvectors[bvID];

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
								reason.push(bv[j]);
							}
						}
						buildComparisonReason(Comparison::geq,bvID,underApprox,reason);
						break;

					}
					if(under+bit>overApprox){
						assert(bitpos==i);
						//this is a conflict
						for(int j = bv.size()-1;j>=i;j--){
							Weight bit = 1<<j;
							lbool val = value(bv[j]);
							if(val==l_True){
								reason.push(~bv[j]);
							}
						}
						buildComparisonReason(Comparison::leq,bvID,overApprox,reason);
						break;
					}
				}
			}

		}else if (marker == conditionthen_prop_marker) {
			//the 'then' condition must hold if the else condition is out of range
			assert(value(p)==l_True);
			reason.push(p);
			Var v = var(p);
			int bvID = getbvID(v);
			writeBounds(bvID);
			assert(isOperationVar(v));
			int opID = getOperationID(v);

			ConditionalID & cID = all_conditionals[opID];
			assert(cID.l==p);
			int bvThenID = cID.bvThenID;
			int bvElseID = cID.bvElseID;

			assert(under_approx[bvElseID]>over_approx[bvID] || over_approx[bvElseID]<under_approx[bvID]);
			buildConditionalPropReason(bvID,opID,reason);
		}else if (marker ==conditionelse_prop_marker){
			assert(value(p)==l_True);
			reason.push(p);
			Var v = var(p);
			int bvID = getbvID(v);
			writeBounds(bvID);
			assert(isOperationVar(v));
			int opID = getOperationID(v);

			ConditionalID & cID = all_conditionals[opID];
			assert(cID.l==~p);
			int bvThenID = cID.bvThenID;
			int bvElseID = cID.bvElseID;
			assert(under_approx[bvThenID]>over_approx[bvID] || over_approx[bvThenID]<under_approx[bvID]);

			buildConditionalPropReason(bvID,opID,reason);

		}else if (marker==conditionarg_prop_marker){
			assert(value(p)==l_True);
			reason.push(p);
			Var v = var(p);
			int bvID = getbvID(v);

			assert(isOperationVar(v));
			int opID = getOperationID(v);

			ConditionalID & cID = all_conditionals[opID];
			int resultID = cID.bvResultID;
			int bvThenID = cID.bvThenID;
			int bvElseID = cID.bvElseID;

			assert(bvID==resultID);
			if(p ==cID.l){
				assert(under_approx[bvElseID]>over_approx[resultID] || over_approx[bvElseID]<under_approx[resultID]);
				buildConditionalArgReason(bvElseID,opID,reason);
			}else{
				assert(p==~cID.l);
				assert(under_approx[bvThenID]>over_approx[resultID] || over_approx[bvThenID]<under_approx[resultID]);
				buildConditionalArgReason(bvThenID,opID,reason);

			}


		}else {
			assert(false);
		}

		toSolver(reason);
		double finish = rtime(1);
		stats_reason_time += finish - start;
		stats_num_reasons++;
		//stats_reason_initial_time+=start-initial_start;
		if(reason.size()<2){
			int a=1;
		}
	}


	void preprocess() {
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
		if(isBVVar(var(l))){

			if(isComparisonVar(var(l))) {

				int bvID = getbvID(var(l));
				int comparisonID = getComparisonID(var(l)); //v-min_edge_var;
				//status.comparisonAltered(bvID, comparisonID);
				trail.push( { true, !sign(l),comparisonID, v });
				bv_needs_propagation[bvID]=true;
				analysis_trail_pos=trail.size()-1;
				//trail.push( { true, !sign(l),edgeID, v });
				if(!alteredBV[bvID]){
					alteredBV[bvID]=true;
					if (altered_bvs.size()==0)
						altered_bvs.push(bvID);
					else{
						//preserve last element in alteredbvs, in case this enqueue happened during propagation
						int lastID = altered_bvs.last();
						altered_bvs.last()=bvID;
						assert(altered_bvs.last()==bvID);
						altered_bvs.push(lastID);
						assert(altered_bvs.last()==lastID);
					}

				}
			}else if (isOperationVar(var(l))){
				//for now, this means a conditional ite
				int bvID = getbvID(var(l));
				int cID = getOperationID(var(l));
				auto & c = all_conditionals[cID];

				Lit condition=c.l;

				int bvThenID=c.bvThenID;
				int bvElseID=c.bvElseID;
				int bvResultID=c.bvResultID;

				assert(bvID==bvResultID);
				assert(var(condition)==var(l));

				int bvCondID;
				if(!sign(l)){
					bvCondID=bvThenID;
				}else{
					bvCondID=bvElseID;
				}

				trail.push( { false, !sign(l),cID, v,true });
				analysis_trail_pos=trail.size()-1;

				bv_needs_propagation[bvID]=true;

				if(!alteredBV[bvID]){
					alteredBV[bvID]=true;
					if (altered_bvs.size()==0)
						altered_bvs.push(bvID);
					else{
						//preserve last element in alteredbvs, in case this enqueue happened during propagation
						int lastID = altered_bvs.last();
						altered_bvs.last()=bvID;
						assert(altered_bvs.last()==bvID);
						altered_bvs.push(lastID);
						assert(altered_bvs.last()==lastID);
					}

				}

				bv_needs_propagation[bvCondID]=true;
				if(!alteredBV[bvCondID]){
					alteredBV[bvCondID]=true;
					if (altered_bvs.size()==0)
						altered_bvs.push(bvCondID);
					else{
						//preserve last element in alteredbvs, in case this enqueue happened during propagation
						int lastID = altered_bvs.last();
						altered_bvs.last()=bvCondID;
						assert(altered_bvs.last()==bvCondID);
						altered_bvs.push(lastID);
						assert(altered_bvs.last()==lastID);
					}

				}


			}else{

				int bvID = getbvID(v);
				trail.push( { false, !sign(l),bvID, v});

				analysis_trail_pos=trail.size()-1;
				if(!alteredBV[bvID]){
					alteredBV[bvID]=true;
					if (altered_bvs.size()==0)
						altered_bvs.push(bvID);
					else{
						//preserve last element in alteredbvs, in case this enqueue happened during propagation
						int lastID = altered_bvs.last();
						altered_bvs.last()=bvID;
						assert(altered_bvs.last()==bvID);
						altered_bvs.push(lastID);
						assert(altered_bvs.last()==lastID);
					}

				}

			}
		}
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
		Weight under_new=0;
		Weight over_new=0;

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

			/*bool changed = (under_old != under_approx[bvID]) || (over_old != over_approx[bvID]);
			if(changed){
				any_changed=true;
				assert(under_new>=under_old);
				assert(over_new<=over_old);
				assert(analysis_trail_pos==trail.size()-1);
				trail.push({bvID,under_old, over_old, under_new, over_new, under_cause_old, over_cause_old, under_causes[bvID],over_causes[bvID]});
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
			}*/
			return false;
		}
		assert_in_range(under_new,bvID);
		assert_in_range(over_new,bvID);
		Cause under_cause_old = under_causes[bvID];
		Cause over_cause_old = over_causes[bvID];
		Cause over_cause_new;
		Cause under_cause_new;

		//under_approx[bvID]=0;
		//over_approx[bvID]=0;
		for(int i = 0;i<bv.size();i++){
			if(var(bv[i])==ignore_bv){
				Weight bit = 1<<i;
				over_new+=bit;
				continue;
			}
			lbool val = value(bv[i]);
			if(val==l_True){
				Weight bit = 1<<i;
				under_new+=bit;
				over_new+=bit;
				under_cause_new.cause_is_bits=true;
			}else if (val==l_False){
				over_cause_new.cause_is_bits=true;
			}else{
				Weight bit = 1<<i;
				over_new+=bit;

			}
		}
		//needed in case a decision was made, to preserve that decision's cause here...
		if( under_old>= under_new){
			under_new=under_old;
			under_cause_new = under_cause_old;
		}

		if(over_old<=over_new){
			over_new=over_old;
			over_cause_new=over_cause_old;
		}

		for(int i = 0;i<conditionals[bvID].size();i++){
			Lit condition = conditionals[bvID][i].l;
			int bvThenID = conditionals[bvID][i].bvThenID;
			int bvElseID = conditionals[bvID][i].bvElseID;
			//assert(bvThenID<bvID);
			//assert(bvElseID<bvID);
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
				under_cause_new.cause_is_condition=true;
				under_cause_new.index=i;
			}
			if(over<over_new){
				over_new=over;
				over_cause_new.clear();
				over_cause_new.cause_is_condition=true;
				over_cause_new.index=i;
			}
		}

		//also need to iterate through the additions that this bv is an argument of...
		for(int i = 0;i<conditional_arguments[bvID].size();i++){

			Lit condition=conditional_arguments[bvID][i].l;
			int other_argID=conditional_arguments[bvID][i].other_argID;
			int resultID=conditional_arguments[bvID][i].resultID;

			assert(other_argID>=0);assert(resultID>=0);
			if(value(condition)==l_True){
				Weight under = under_approx[resultID];
				Weight over = over_approx[resultID];
	/*			clip_under(under,bvID);
				clip_under(over,bvID);*/
				if(under >under_new){
					under_new=under;
					under_cause_new.clear();
					under_cause_new.cause_is_condition_argument=true;
					under_cause_new.index=i;
				}
				if(over<over_new){
					over_new=over;
					over_cause_new.clear();
					over_cause_new.cause_is_condition_argument=true;
					over_cause_new.index=i;
				}
			}
		}


		for(int i = 0;i<additions[bvID].size();i++){
			int aID = additions[bvID][i].aID;
			int bID = additions[bvID][i].bID;
			//assert(aID<bvID);
			//assert(bID<bvID);
			Weight under = under_approx[aID] +  under_approx[bID];
			Weight over = over_approx[aID] +  over_approx[bID];
			clip_over(under,bvID);
			clip_over(over,bvID);
			if(under >under_new){
				under_new=under;
				under_cause_new.clear();
				under_cause_new.cause_is_addition=true;
				under_cause_new.index=i;
			}
			if(over<over_new){
				over_new=over;
				over_cause_new.clear();
				over_cause_new.cause_is_addition=true;
				over_cause_new.index=i;
			}
		}

		//also need to iterate through the additions that this bv is an argument of...
		for(int i = 0;i<addition_arguments[bvID].size();i++){
			int other_argID = addition_arguments[bvID][i].other_argID;
			int sumID = addition_arguments[bvID][i].sumID;
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
				under_cause_new.cause_is_addition_argument=true;
				under_cause_new.index=i;
			}
			if(over<over_new){
				over_new=over;
				over_cause_new.clear();
				over_cause_new.cause_is_addition_argument=true;
				over_cause_new.index=i;
			}
		}
		assert_in_range(under_new,bvID);
		assert_in_range(over_new,bvID);
		for(int i = compares[bvID].size()-1;i>=0;i--){
			int cID = compares[bvID][i];
			if(cID==ignoreCID){
				continue;
			}
			ComparisonID & c = comparisons[cID];
			Comparison op = c.op();

			bool setOver=false;
			Weight w = c.w;
			assert(!c.bvCompare());
			switch(op){
				case Comparison::lt:
					if(value( c.l)==l_True && over_new>=w){
						if(w-1>=getUnderApprox(bvID,true)){
							over_new=w-1;
							setOver=true;
						 }
					}
					break;
				case Comparison::leq:
					if(value( c.l)==l_True && over_new>w){
						if(w>=getUnderApprox(bvID,true)){
							over_new=w;
							setOver=true;
						 }
					}
					break;
				case Comparison::gt:
					if (value(c.l)==l_False && over_new>w){
						if(w>=getUnderApprox(bvID,true)){
							over_new=w;
							setOver=true;
						 }
					}
					break;
				case Comparison::geq:
				default:
					if (value(c.l)==l_False && over_new>=w){
						if(w-1>=getUnderApprox(bvID,true)){
							over_new=w-1;
							setOver=true;
						 }
					}
					break;
			}

			if(setOver){
				over_cause_new.clear();
				over_cause_new.cause_is_comparison=true;
				over_cause_new.index=cID;
			}
		}


		for(int cID:compares[bvID]){
			if(cID==ignoreCID){
				continue;
			}
			ComparisonID & c = comparisons[cID];
			Comparison op = c.op();
			bool setUnder=false;
			assert(!c.bvCompare());
			Weight w = c.w;

			switch(op){
				case Comparison::lt:
					 if (value(c.l)==l_False && under_new<w){
						 if(w<=getOverApprox(bvID,true)){
							 under_new=w;
							 setUnder=true;
						 }
					}
					break;
				case Comparison::leq:
					 if (value(c.l)==l_False && under_new<=w){
						 if(w<=getOverApprox(bvID,true)){
							 under_new=w+1;
							setUnder=true;
						 }
					}
					break;
				case Comparison::gt:
					if(value( c.l)==l_True && under_new<=w){
						 if(w+1<=getOverApprox(bvID,true)){
							 under_new=w+1;
							 setUnder=true;
						 }
					}
					break;
				case Comparison::geq:
				default:
					if(value( c.l)==l_True && under_new<w){
						 if(w<=getOverApprox(bvID,true)){
							 under_new=w;
							 setUnder=true;
						 }
					}
					break;
			}

			if(setUnder){
				under_cause_new.clear();
				under_cause_new.cause_is_comparison=true;
				under_cause_new.index=cID;
			}

		}

		assert_in_range(under_new,bvID);
		assert_in_range(over_new,bvID);
		for(int i = bvcompares[bvID].size()-1;i>=0;i--){
			int cID = bvcompares[bvID][i];
			if(cID==ignoreCID){
				continue;
			}
			ComparisonID & c = comparisons[cID];
			Comparison op = c.op();

			bool setOver=false;

			assert(c.bvCompare());
			Weight w = over_approx[c.compareID];
			switch(op){
				case Comparison::lt:
					if(value( c.l)==l_True && over_new>=w){
						if(w-1>=getUnderApprox(bvID,true)){
							over_new=w-1;
							setOver=true;
						}
					}
					break;
				case Comparison::leq:
					if(value( c.l)==l_True && over_new>w){
						if(w>=getUnderApprox(bvID,true)){
							over_new=w;
							setOver=true;
						}
					}
					break;
				case Comparison::gt:
					if (value(c.l)==l_False && over_new>w){
						if(w>=getUnderApprox(bvID,true)){
							over_new=w;
							setOver=true;
						}
					}
					break;
				case Comparison::geq:
				default:
					if (value(c.l)==l_False && over_new>=w){
						if(w-1>=getUnderApprox(bvID,true)){
							over_new=w-1;
							setOver=true;
						}
					}
					break;
			}

			if(setOver){
				over_cause_new.clear();
				over_cause_new.cause_is_comparison=true;
				over_cause_new.index=cID;
			}
		}

		for(int cID:bvcompares[bvID]){
			if(cID==ignoreCID){
				continue;
			}
			ComparisonID & c = comparisons[cID];
			Comparison op = c.op();
			bool setUnder=false;
			assert(c.bvID>=0);
			Weight w = under_approx[c.compareID];

			switch(op){
				case Comparison::lt:
					 if (value(c.l)==l_False && under_new<w){
						 if(w<=getOverApprox(bvID,true)){
							 under_new=w;
							 setUnder=true;
						 }
					}
					break;
				case Comparison::leq:
					 if (value(c.l)==l_False && under_new<=w){
						 if(w+1<=getOverApprox(bvID,true)){
							under_new=w+1;
							setUnder=true;
						 }
					}
					break;
				case Comparison::gt:
					if(value( c.l)==l_True && under_new<=w){
						 if(w+1<=getOverApprox(bvID,true)){
							 under_new=w+1;
							 setUnder=true;
						 }
					}
					break;
				case Comparison::geq:
				default:
					if(value( c.l)==l_True && under_new<w){
						 if(w<=getOverApprox(bvID,true)){
							 under_new=w;
							 setUnder=true;
						 }
					}
					break;
			}

			if(setUnder){
				under_cause_new.clear();
				under_cause_new.cause_is_comparison=true;
				under_cause_new.index=cID;
			}

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
			over_cause_new.refined_cause=true;
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
			under_cause_new.refined_cause=true;
		}
		assert_in_range(under_new,bvID);
		assert_in_range(over_new,bvID);

		//int width = bitvectors[bvID].size();
		//Weight max_val = ((1L)<<width)-1;
		if(under_new>over_approx0[bvID]){
			under_new=over_approx0[bvID];
		}
		if(over_new>over_approx0[bvID]){
			over_new=over_approx0[bvID];
		}
		if(under_new<under_approx0[bvID]){
			under_new=under_approx0[bvID];
		}
		if(over_new<under_approx0[bvID]){
			over_new=under_approx0[bvID];
		}
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


			fprintf(opt_write_learnt_clauses,"learnt bound ");
			fprintf(opt_write_learnt_clauses," %d <= ", bvID);

			std::stringstream ss;
			ss<<over_approx[bvID] << " ";
			//fprintf(opt_write_learnt_clauses,"learnt fact bv %d %s %d\n", bvID, ss.str().c_str(),dimacs(toSolver(l)) );
			fprintf(opt_write_learnt_clauses,"%s",ss.str().c_str());

			if(over_approx[bvID]!=over_approx0[bvID]){

				for (int i = 0;i<nVars();i++){
					if(value(i)!=l_Undef){
						Lit l = mkLit(i,value(i)==l_True);
						assert(value(l)==l_False);
						fprintf(opt_write_learnt_clauses," %d",dimacs(toSolver(l)));
					}
				}
			}
			fprintf(opt_write_learnt_clauses," 0\n");
			++bound_num;
			if(++opt_n_learnts==44231){
				int a=1;
			}
			fprintf(opt_write_learnt_clauses,"learnt bound ");
			fprintf(opt_write_learnt_clauses," %d >= ", bvID);
			//std::cout<<under_approx[bvID] << " ";
			std::stringstream ss2;
			ss2<<under_approx[bvID] << " ";
			//fprintf(opt_write_learnt_clauses,"learnt fact bv %d %s %d\n", bvID, ss.str().c_str(),dimacs(toSolver(l)) );
			fprintf(opt_write_learnt_clauses,"%s",ss2.str().c_str());
			if(under_approx[bvID]!=under_approx0[bvID]){
				for (int i = 0;i<nVars();i++){
					if(value(i)!=l_Undef){
						Lit l = mkLit(i,value(i)==l_True);
						assert(value(l)==l_False);
						fprintf(opt_write_learnt_clauses," %d",dimacs(toSolver(l)));
					}
				}
			}
			fprintf(opt_write_learnt_clauses," 0\n");
			if(++bound_num>=80){
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
				conflict.push(l);
			}
		}
	}

	bool checkApproxUpToDate(int bvID, Weight * under_store=nullptr,Weight*over_store=nullptr){
#ifndef NDEBUG
		if(bvID==2966){
			int a=1;
		}
		Weight under =0;
		Weight over=0;
		vec<Lit> & bv = bitvectors[bvID];
		if(eq_bitvectors[bvID]!=bvID && eq_bitvectors[bvID]>-1 ){
			int eqID = eq_bitvectors[bvID];
			under=under_approx[eqID];
			over=over_approx[eqID];
			assert(over_approx[bvID]==over);
			assert(under_approx[bvID]==under);
			return true;
		}else{

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

			for(int i = 0;i< conditionals[bvID].size();i++){
				int thenID = conditionals[bvID][i].bvThenID;
				int elseID = conditionals[bvID][i].bvElseID;
				Lit condition = conditionals[bvID][i].l;

				if(value(condition)==l_Undef){
					if(under_approx[bvID]< min( under_approx[thenID],under_approx[elseID]) ||
							over_approx[bvID]> max( over_approx[thenID],over_approx[elseID])){
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
			}

			for(int i = 0;i< conditional_arguments[bvID].size();i++){
				Lit condition=conditional_arguments[bvID][i].l;
				int other_argID=conditional_arguments[bvID][i].other_argID;
				int resultID=conditional_arguments[bvID][i].resultID;
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
			}

			for(int i = 0;i<additions[bvID].size();i++){
				int aID = additions[bvID][i].aID;
				int bID = additions[bvID][i].bID;
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
			}

			//also need to iterate through the additions that this bv is an argument of...
			for(int i = 0;i<addition_arguments[bvID].size();i++){
				int other_argID = addition_arguments[bvID][i].other_argID;
				int sumID = addition_arguments[bvID][i].sumID;

				Weight under_add = under_approx[sumID] -  over_approx[other_argID];
				Weight over_add = over_approx[sumID] -  under_approx[other_argID];

				if(under_add >under){
					under=under_add;
				}
				if(over_add<over){
					over=over_add;
				}
			}

			for(int cID:compares[bvID]){
				ComparisonID & c = comparisons[cID];
				Comparison op = c.op();

				lbool val = value(c.l);
				if(!c.bvCompare()){
					switch(op){
						case Comparison::lt:
							if(value( c.l)==l_True && over>=c.w){
								over=c.w-1;
							}else if (value(c.l)==l_False && under<c.w){
								under=c.w;
							}
							break;
						case Comparison::leq:
							if(value( c.l)==l_True && over>c.w){
								over=c.w;
							}else if (value(c.l)==l_False && under<=c.w){
								under=c.w+1;
							}
							break;
						case Comparison::gt:
							if(value( c.l)==l_True && under<=c.w){
								under=c.w+1;
							}else if (value(c.l)==l_False && over>c.w){
								over=c.w;
							}
							break;
						case Comparison::geq:
						default:
							if(value( c.l)==l_True && under<c.w){
								under=c.w;
							}else if (value(c.l)==l_False && over>=c.w){
								over=c.w-1;
							}
							break;
					}
				}else{


				}
			}
			for(int cID:bvcompares[bvID]){
				ComparisonID & c = comparisons[cID];
				Comparison op = c.op();

				lbool val = value(c.l);
				if(c.bvCompare()){
					Weight under_w = under_approx[c.compareID];
					Weight over_w = over_approx[c.compareID];

					switch(op){

						case Comparison::lt:
							if(value( c.l)==l_True && over>=over_w){
								over=over_w-1;
							}else if (value(c.l)==l_False && under<under_w){
								under=under_w;
							}
							break;
						case Comparison::leq:
							if(value( c.l)==l_True && over>over_w){
								over=over_w;
							}else if (value(c.l)==l_False && under<=under_w){
								under=under_w+1;
							}
							break;
						case Comparison::gt:
							if(value( c.l)==l_True && under<=under_w){
								under=under_w+1;
							}else if (value(c.l)==l_False && over>over_w){
								over=over_w;
							}
							break;
						case Comparison::geq:
						default:
							if(value( c.l)==l_True && under<under_w){
								under=under_w;
							}else if (value(c.l)==l_False && over>=over_w){
								over=over_w-1;
							}
							break;
					}
				}else{


				}
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
		if (!under_causes[bvID].cause_is_decision)
			assert(under==under_approx[bvID]);
		if(!over_causes[bvID].cause_is_decision)
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
	bool propagateTheory(vec<Lit> & conflict){
		return propagateTheory(conflict,false);
	}
	bool propagateTheory(vec<Lit> & conflict, bool force_propagation) {
		static int realprops = 0;
		stats_propagations++;

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
		if(stats_propagations==17){
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
			if(bvID==12 || bvID == 8 || bvID==11){
				int a=1;
			}
		//for(int bvID = 0;bvID<bitvectors.size();bvID++){
			assert(alteredBV[bvID]);
			Weight  underApprox_prev = under_approx[bvID];
			Weight  overApprox_prev = over_approx[bvID];
			Cause prev_under_cause = over_causes[bvID];
			Cause prev_over_cause = over_causes[bvID];
			if(stats_bv_propagations==17){
				int a=1;
			}
			bool changed = updateApproximations(bvID);
			changed |=bv_needs_propagation[bvID];
			if(!changed){
				assert( under_approx[bvID]<= over_approx[bvID]);
				stats_bv_skipped_propagations++;
				altered_bvs.pop(); //this isn't safe to do, because recently enforced comparisons need to be propagated, even if the under/over approx didn't change.
				alteredBV[bvID]=false;
				continue;
			}
			bv_needs_propagation[bvID]=false;

			Weight & underApprox = under_approx[bvID];
			Weight & overApprox = over_approx[bvID];
			stats_bv_propagations++;
			if(stats_bv_propagations==18){
				int a=1;
			}

			/*printf("iter %d, bv %d, under ",realprops , bvID); //: %d, over %d\n", bvID, underApprox,overApprox);
			std::cout<<underApprox << " over ";
			std::cout<<overApprox << "\n";
			fflush(stdout);*/

			vec<Lit> & bv = bitvectors[bvID];
			Weight under;
			Weight over;
			bool new_change;
			do{
				if(underApprox>overApprox){
					double startconftime = rtime(2);
					stats_num_conflicts++;
					if(opt_verb>1){
						printf("bv approximation update conflict %ld\n", stats_num_conflicts);
					}
					propagationtime += startconftime - startproptime;
					//this is already a conflict;
					analyzeValueReason(Comparison::geq, bvID,underApprox,conflict);
					analyzeValueReason(Comparison::leq, bvID,overApprox,conflict);

					analyze(conflict);
					toSolver(conflict);
					stats_conflict_time+=rtime(2)-startconftime;
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
								propagationtime += startconftime - startproptime;
								//this is a conflict
								for(int j = bv.size()-1;j>=i;j--){
									Weight bit = 1<<j;
									lbool val = value(bv[j]);
									if(val==l_True){
										conflict.push(~bv[j]);
									}
								}
								stats_num_conflicts++;stats_bit_conflicts++;
								if(opt_verb>1){
									printf("bv bit conflict %ld\n", stats_num_conflicts);
								}
								buildComparisonReason(Comparison::leq,bvID,overApprox,conflict);
								toSolver(conflict);
								stats_conflict_time+=rtime(2)-startconftime;
								return false;

						}
					}else if (val==l_False){
						over-=bit;
						assert_in_range(under,bvID);
						if(over<underApprox){
								double startconftime = rtime(2);
								propagationtime += startconftime - startproptime;
								//this is a conflict. Either this bit, or any previously assigned false bit, must be true, or the underapprox must be larger than it currently is.
								//is this really the best way to handle this conflict?
								for(int j = bv.size()-1;j>=i;j--){
									Weight bit = 1<<j;
									lbool val = value(bv[j]);
									if(val==l_False){
										conflict.push(bv[j]);
									}
								}
								stats_num_conflicts++;stats_bit_conflicts++;
								if(opt_verb>1){
									printf("bv bit conflict %ld\n", stats_num_conflicts);
								}
								buildComparisonReason(Comparison::geq,bvID,underApprox,conflict);
								toSolver(conflict);
								stats_conflict_time+=rtime(2)-startconftime;
								return false;

						}
					}else{
						assert_in_range(over-bit,bvID);
						assert_in_range(under+bit,bvID);
						if(over-bit < underApprox){
							enqueue(l, bvprop_marker);
							under+=bit;
						}else if (under+bit>overApprox){
							enqueue(~l, bvprop_marker);
							over-=bit;
						}
					}
				}
				new_change = updateApproximations(bvID);
				changed|=new_change;
				/*if(new_change){
					printf("iter %d, bv %d, under ",realprops , bvID); //: %d, over %d\n", bvID, underApprox,overApprox);
						std::cout<<underApprox << " over ";
						std::cout<<overApprox << "\n";
						fflush(stdout);
				}*/
			}while(new_change);//the bit assignment updates above can force a more precise over or under approximation, which can in turn lead to further bit assignments (I think this can happen?).

			for(int i = 0;i<conditionals[bvID].size();i++){
				int cID =  conditionals[bvID][i].cID;
				Lit condition = conditionals[bvID][i].l;
				int bvThenID = conditionals[bvID][i].bvThenID;
				int bvElseID = conditionals[bvID][i].bvElseID;
				//assert(bvThenID<bvID);
				//assert(bvElseID<bvID);

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
						propagationtime += startconftime - startproptime;
						stats_num_conflicts++;stats_addition_conflicts++;
						if(opt_verb>1){
							printf("bv condition conflict %ld\n", stats_num_conflicts);
						}
						if(value(condition)==l_True){
							conflict.push(~condition);
						}else if(value(condition)==l_False){
							conflict.push(condition);
						}
						buildConditionalReason(bvID,cID,conflict);
						toSolver(conflict);
						stats_conflict_time+=rtime(2)-startconftime;
						return false;
					}else if (overApprox<under){
						double startconftime = rtime(2);
						propagationtime += startconftime - startproptime;
						if(opt_verb>1){
							printf("bv condition conflict %ld\n", stats_num_conflicts);
						}
						stats_num_conflicts++;stats_addition_conflicts++;
						if(value(condition)==l_True){
							conflict.push(~condition);
						}else if(value(condition)==l_False){
							conflict.push(condition);
						}
						buildConditionalReason(bvID,cID,conflict);
						toSolver(conflict);
						stats_conflict_time+=rtime(2)-startconftime;
						return false;
					}

					if(value(condition)==l_True){
						if((underApprox> under_approx[bvThenID]) || (overApprox < over_approx[bvThenID])){
							//the other bv needs to be updated
							if(!alteredBV[bvThenID]){
								alteredBV[bvThenID]=true;
								assert(altered_bvs.last()==bvID);
								altered_bvs.last()=bvThenID;
								assert(altered_bvs.last()==bvThenID);
								altered_bvs.push(bvID);
								assert(altered_bvs.last()==bvID);
							}
						}
					}else if (value(condition)==l_False){
						if((underApprox> under_approx[bvElseID]) || (overApprox < over_approx[bvElseID])){
							//the other bv needs to be updated
							if(!alteredBV[bvElseID]){
								alteredBV[bvElseID]=true;
								assert(altered_bvs.last()==bvID);
								altered_bvs.last()=bvElseID;
								assert(altered_bvs.last()==bvElseID);
								altered_bvs.push(bvID);
								assert(altered_bvs.last()==bvID);
							}
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
								buildConditionalReason(bvID,cID,conflict);
								toSolver(conflict);
								//buildConditionalReason(bvID,bvElseID,underApprox,conflict);
								return false;
							}else{
								enqueue(~condition,  conditionelse_prop_marker);
								changed=true;
							}

						}else if(under_approx[bvElseID]>over_approx[bvID] || over_approx[bvElseID]<under_approx[bvID]){
							enqueue(condition,conditionthen_prop_marker);
							changed=true;
							//condition must be true
						}
					}
				}

			}

			for(int i = 0;i<conditional_arguments[bvID].size();i++){
				int cID =  conditional_arguments[bvID][i].cID;
				Lit condition = conditional_arguments[bvID][i].l;
				int other_argID = conditional_arguments[bvID][i].other_argID;
				int resultID = conditional_arguments[bvID][i].resultID;
				assert(other_argID>=0);assert(resultID>=0);

				if(value(condition)==l_True){
					if(under_approx[bvID]>over_approx[resultID] || over_approx[bvID]<under_approx[resultID]){
						//this is a conflict
						conflict.push(~condition);
						buildConditionalArgReason(bvID,cID,conflict);
						toSolver(conflict);
						return false;
					}else if(under_approx[bvID]>under_approx[resultID] || over_approx[bvID]<over_approx[resultID]){
						if(!alteredBV[resultID]){
							alteredBV[resultID]=true;
							assert(altered_bvs.last()==bvID);
							altered_bvs.last()=resultID;
							assert(altered_bvs.last()==resultID);
							altered_bvs.push(bvID);
							assert(altered_bvs.last()==bvID);
						}
					}



				}else if (value(condition)==l_Undef){
					//can potentially propagate the falsehood of the conditional, if this bv is out of range of the result bv
					if(under_approx[bvID]>over_approx[resultID] || over_approx[bvID]<under_approx[resultID]){
						enqueue(~condition,conditionarg_prop_marker);
						if(!alteredBV[other_argID]){
							//exit(8);
							alteredBV[other_argID]=true;
							assert(altered_bvs.last()==bvID);
							altered_bvs.last()=other_argID;
							assert(altered_bvs.last()==other_argID);
							altered_bvs.push(bvID);
							assert(altered_bvs.last()==bvID);
						}
						if(!alteredBV[resultID]){
							alteredBV[resultID]=true;
							assert(altered_bvs.last()==bvID);
							altered_bvs.last()=resultID;
							assert(altered_bvs.last()==resultID);
							altered_bvs.push(bvID);
							assert(altered_bvs.last()==bvID);
						}
					}else{
						Weight under_min = min(under_approx[bvID],under_approx[other_argID]);
						Weight over_max = max(over_approx[bvID],over_approx[other_argID]);

						if (under_min>under_approx[resultID] || over_max<over_approx[resultID]){
							//ensure that it is a change to this bv's bounds that caused the change here...
							//assert(under_approx[bvID]>under_approx[other_argID] || over_approx[bvID]<over_approx[resultID]);

							if(!alteredBV[resultID]){
								alteredBV[resultID]=true;
								assert(altered_bvs.last()==bvID);
								altered_bvs.last()=resultID;
								assert(altered_bvs.last()==resultID);
								altered_bvs.push(bvID);
								assert(altered_bvs.last()==bvID);
							}
						}

					}
				}
			}



			for(int i = 0;i<additions[bvID].size();i++){
				int aID = additions[bvID][i].aID;
				int bID = additions[bvID][i].bID;
				//assert(aID<bvID);
				//assert(bID<bvID);
				Weight under = under_approx[aID] +  under_approx[bID];
				Weight over = over_approx[aID] +  over_approx[bID];
				clip_over(under,bvID);
				clip_over(over,bvID);
				if(underApprox>over){
					//then we have a conflict
					double startconftime = rtime(2);
					propagationtime += startconftime - startproptime;
					stats_num_conflicts++;stats_addition_conflicts++;
					if(opt_verb>1){
						printf("bv addition conflict %ld\n", stats_num_conflicts);
					}
					buildAdditionReason(bvID,i,conflict);
					toSolver(conflict);
					stats_conflict_time+=rtime(2)-startconftime;
					return false;
				}else if (overApprox<under){
					double startconftime = rtime(2);
					propagationtime += startconftime - startproptime;
					if(opt_verb>1){
						printf("bv addition conflict %ld\n", stats_num_conflicts);
					}
					stats_num_conflicts++;stats_addition_conflicts++;
					buildAdditionReason(bvID,i,conflict);
					toSolver(conflict);
					stats_conflict_time+=rtime(2)-startconftime;
					return false;
				}
				Weight under_arg_b = underApprox - over_approx[bID];
				Weight over_arg_b = overApprox - under_approx[bID];
				clip_under(under_arg_b,bvID);
				clip_under(over_arg_b,bvID);
				//this check may be especially important when either aID or bID is really a constant...
				if((under_arg_b> under_approx[aID]) || (over_arg_b < over_approx[aID])){
					//the other bv needs to be updated
					if(!alteredBV[aID]){
						alteredBV[aID]=true;
						assert(altered_bvs.last()==bvID);
						altered_bvs.last()=aID;
						assert(altered_bvs.last()==aID);
						altered_bvs.push(bvID);
						assert(altered_bvs.last()==bvID);
					}
				}
				Weight under_arg_a =underApprox - over_approx[aID];
				Weight over_arg_a =  overApprox - under_approx[aID];
				clip_under(under_arg_a,bvID);
				clip_under(over_arg_a,bvID);
				if((under_arg_a > under_approx[bID]) || (over_arg_a< over_approx[bID])){
					//the other bv needs to be updated
					if(!alteredBV[bID]){
						alteredBV[bID]=true;
						assert(altered_bvs.last()==bvID);
						altered_bvs.last()=bID;
						assert(altered_bvs.last()==bID);
						altered_bvs.push(bvID);
						assert(altered_bvs.last()==bvID);
					}
				}
			}

			//also need to iterate through the additions that this bv is an argument of...
			for(int i = 0;i<addition_arguments[bvID].size();i++){
				int other_argID = addition_arguments[bvID][i].other_argID;
				int sumID = addition_arguments[bvID][i].sumID;
				assert(other_argID>=0);assert(sumID>=0);

				Weight under = under_approx[sumID] -  over_approx[other_argID];
				Weight over = over_approx[sumID] -  under_approx[other_argID];
				clip_under(under,bvID);
				clip_under(over,bvID);
				if(underApprox>over){
					//then we have a conflict
					double startconftime = rtime(2);
					propagationtime += startconftime - startproptime;
					stats_num_conflicts++;stats_addition_conflicts++;
					if(stats_num_conflicts==416){
						int a=1;
					}
					if(opt_verb>1){
						printf("bv addition arg conflict %ld\n", stats_num_conflicts);
					}
					buildAdditionArgReason(bvID,i,conflict);
					toSolver(conflict);
					stats_conflict_time+=rtime(2)-startconftime;
					return false;
				}else if (overApprox<under){
					double startconftime = rtime(2);
					propagationtime += startconftime - startproptime;
					stats_num_conflicts++;stats_addition_conflicts++;
					if(opt_verb>1){
						printf("bv addition arg conflict %ld\n", stats_num_conflicts);
					}
					buildAdditionArgReason(bvID,i,conflict);
					toSolver(conflict);
					stats_conflict_time+=rtime(2)-startconftime;
					return false;
				}
				Weight under_arg = underApprox + under_approx[other_argID];
				clip_over(under_arg,bvID);
				Weight over_arg =  overApprox + over_approx[other_argID];
				clip_over(over_arg,bvID);
				//this check may be especially important when either aID or bID is really a constant...
				if((under_arg > under_approx[sumID]) || (over_arg< over_approx[sumID])){
					//the other bv needs to be updated

					if(!alteredBV[sumID]){
						//exit(8);
						alteredBV[sumID]=true;
						assert(altered_bvs.last()==bvID);
						altered_bvs.last()=sumID;
						assert(altered_bvs.last()==sumID);
						altered_bvs.push(bvID);
						assert(altered_bvs.last()==bvID);
					}
				}
				Weight under_sum =  under_approx[sumID]  - overApprox;
				clip_under(under_sum,bvID);
				Weight over_sum =  over_approx[sumID] - underApprox;
				clip_under(over_sum,bvID);
				if( (under_sum> under_approx[other_argID]) || (over_sum< over_approx[other_argID])){
					//the other bv needs to be updated

					if(!alteredBV[other_argID]){
						//exit(8);
						alteredBV[other_argID]=true;
						assert(altered_bvs.last()==bvID);
						altered_bvs.last()=other_argID;
						assert(altered_bvs.last()==other_argID);
						altered_bvs.push(bvID);
						assert(altered_bvs.last()==bvID);
					}
				}
			}


			vec<int> & compare = compares[bvID];
			//update over approx lits
			for(int i = 0;i<compare.size();i++){
				int cID = compare[i];
				ComparisonID & c = comparisons[cID];
				Comparison op = c.op();
				Weight & to = c.w;
				Lit l =  c.l;

				if((op==Comparison::lt && overApprox<to) ||
						(op==Comparison::leq && overApprox<=to)){
					if(value(l)==l_True){
						//do nothing

					}else if (value(l)==l_False){
						double startconftime = rtime(2);
						propagationtime += startconftime - startproptime;
						stats_num_conflicts++;stats_compare_conflicts++;
						if(opt_verb>1){
							printf("bv comparison conflict %ld\n", stats_num_conflicts);
						}
						assert(value(l)==l_False);
						assert(dbg_value(l)==l_False);
						conflict.push(l);
						buildComparisonReason(op,bvID,to,conflict);
						toSolver(conflict);
						stats_conflict_time+=rtime(2)-startconftime;
						return false;
					}else {
						assert(value(l)==l_Undef);
						enqueue(l,comparisonprop_marker);
					}
				}else if((op==Comparison::gt && overApprox<=to) ||
						(op==Comparison::geq && overApprox<to)){
					if(value(l)==l_True){
						double startconftime = rtime(2);
						propagationtime += startconftime - startproptime;
						stats_num_conflicts++;stats_compare_conflicts++;
						if(opt_verb>1){
							printf("bv comparison conflict %ld\n", stats_num_conflicts);
						}
						conflict.push(~l);
						buildComparisonReason(-op,bvID,to,conflict);
						toSolver(conflict);
						stats_conflict_time+=rtime(2)-startconftime;
						return false;
					}else if (value(l)==l_False){

					}else {
						assert(value(l)==l_Undef);
						enqueue(~l, comparisonprop_marker);
					}
				}
			}


			for(int i=compare.size()-1;i>=0;i--){
				int cID = compare[i];
				ComparisonID & c = comparisons[cID];
				assert(!c.bvCompare());
				Comparison op = c.op();
				Weight & to = c.w;
				Lit l =  c.l;

				if((op==Comparison::lt && underApprox>=to) ||
						(op==Comparison::leq && underApprox>to)){
					if(value(l)==l_True){
						double startconftime = rtime(2);
						propagationtime += startconftime - startproptime;
						stats_num_conflicts++;stats_compare_conflicts++;
						if(opt_verb>1){
							printf("bv comparison conflict %ld\n", stats_num_conflicts);
						}
						conflict.push(~l);
						buildComparisonReason(-op,bvID,to,conflict);
						toSolver(conflict);
						stats_conflict_time+=rtime(2)-startconftime;
						return false;
					}else if (value(l)==l_False){
						//do nothing

					}else {
						assert(value(l)==l_Undef);
						enqueue(~l, comparisonprop_marker);
					}
				}else if((op==Comparison::gt && underApprox>to) ||
						(op==Comparison::geq && underApprox>=to)){
					if(value(l)==l_True){

					}else if (value(l)==l_False){
						double startconftime = rtime(2);
						propagationtime += startconftime - startproptime;
						stats_num_conflicts++;stats_compare_conflicts++;
						if(opt_verb>1){
							printf("bv comparison conflict %ld\n", stats_num_conflicts);
						}
						conflict.push(l);
						buildComparisonReason(op,bvID,to,conflict);
						toSolver(conflict);
						stats_conflict_time+=rtime(2)-startconftime;
						return false;
					}else {
						assert(value(l)==l_Undef);
						enqueue(l,comparisonprop_marker);
					}
				}
			}

			//comparisons to bitvectors.
			vec<int> & bvcompare = bvcompares[bvID];

			for(int i = 0;i<bvcompare.size();i++){
				int cID = bvcompare[i];
				ComparisonID & c = comparisons[cID];
				assert(c.bvCompare());
				Comparison op = c.op();
				int compareID = c.compareID;
				assert(compareID>=0);
				Lit l =  c.l;
				lbool val = value(l);
				if(value(l)==l_False){
					l=~l;
					op=-op;
				}
				assert(value(l)!=l_False);
				Weight & under_compare = under_approx[compareID];

				if((op==Comparison::lt && overApprox<under_compare) ||
						(op==Comparison::leq && overApprox<=under_compare)){
					if(value(l)==l_True){
						//do nothing

					}/*else if (value(l)==l_False){

						assert(value(l)==l_False);
						assert(dbg_value(l)==l_False);
						conflict.push(l);
						buildValueReasonBV(op,bvID,compareID,conflict);
						toSolver(conflict);
						return false;
					}*/else {
						assert(value(l)==l_Undef);
						enqueue(l,comparisonprop_marker);
					}
				}else if((op==Comparison::gt && overApprox<=under_compare) ||
						(op==Comparison::geq && overApprox<under_compare)){
					if(value(l)==l_True){
						double startconftime = rtime(2);
						propagationtime += startconftime - startproptime;
						stats_num_conflicts++;stats_bv_compare_conflicts++;
						if(opt_verb>1){
							printf("bv bv comparison conflict %ld\n", stats_num_conflicts);
						}
						conflict.push(~l);
						buildComparisonReasonBV(-op,bvID,compareID,conflict);
						toSolver(conflict);
						stats_conflict_time+=rtime(2)-startconftime;
						return false;
					}/*else if (value(l)==l_False){

					}*/else {
						assert(value(l)==l_Undef);
						enqueue(~l, comparisonprop_marker);
					}
				}
				if(value(l)==l_True &&((op==Comparison::lt && underApprox>=under_compare) ||
						(op==Comparison::leq && underApprox>under_compare))){
					//the other bv needs to be updated
					if(!alteredBV[compareID]){
						alteredBV[compareID]=true;
						assert(altered_bvs.last()==bvID);
						altered_bvs.last()=compareID;
						assert(altered_bvs.last()==compareID);
						altered_bvs.push(bvID);
						assert(altered_bvs.last()==bvID);
					}
				}
			}


			for(int i=bvcompare.size()-1;i>=0;i--){
				int cID = bvcompare[i];
				ComparisonID & c = comparisons[cID];
				Comparison op = c.op();
				int compareID = c.compareID;
				assert(compareID>=0);
				Lit l =  c.l;
				lbool val = value(l);
				if(value(l)==l_False){
					l=~l;
					op=-op;
				}
				assert(value(l)!=l_False);
				//updateApproximations(compareID);
				Weight & over_compare = over_approx[compareID];
				if((op==Comparison::lt && underApprox>=over_compare) ||
						(op==Comparison::leq && underApprox>over_compare)){
					if(value(l)==l_True){
						double startconftime = rtime(2);
						propagationtime += startconftime - startproptime;
						stats_num_conflicts++;stats_bv_compare_conflicts++;
						if(opt_verb>1){
							printf("bv bv comparison conflict %ld\n", stats_num_conflicts);
						}
						conflict.push(~l);
						buildComparisonReasonBV(-op,bvID,compareID,conflict);
						toSolver(conflict);
						stats_conflict_time+=rtime(2)-startconftime;
						return false;
					}/*else if (value(l)==l_False){
						//do nothing

					}*/else {
						assert(value(l)==l_Undef);
						enqueue(~l, comparisonprop_marker);
					}
				}else if((op==Comparison::gt && underApprox>over_compare) ||
						(op==Comparison::geq && underApprox>=over_compare)){
					if(value(l)==l_True){

					}/*else if (value(l)==l_False){

						conflict.push(l);
						buildValueReasonBV(op,bvID,compareID,conflict);
						toSolver(conflict);
						return false;
					}*/else {
						assert(value(l)==l_Undef);
						enqueue(l, comparisonprop_marker);
					}
				}
				//we can also update the other bv's approximation, possibly:
				if(value(l)==l_True &&( (op==Comparison::gt && overApprox<=over_compare) ||
						(op==Comparison::geq && overApprox<over_compare))){
					//the other bv needs to be updated
					if(!alteredBV[compareID]){
						alteredBV[compareID]=true;
						assert(altered_bvs.last()==bvID);
						altered_bvs.last()=compareID;
						assert(altered_bvs.last()==compareID);
						altered_bvs.push(bvID);
						assert(altered_bvs.last()==bvID);
					}
				}

			}
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

	void buildConditionalPropReason(int bvID,int conditionalID, vec<Lit> & conflict){
		dbg_no_pending_analyses();
		assert(eq_bitvectors[bvID]==bvID);

		stats_build_condition_reason++;
		Weight  over_cur = over_approx[bvID];
		Weight  under_cur = under_approx[bvID];
		//assert(checkApproxUpToDate(bvID));


		//the reason that the addition is over is the reason that
		//bvID > addition_under, or the reason that addition_under>= its current value.
		int bvThenID = all_conditionals[conditionalID].bvThenID;
		int bvElseID = all_conditionals[conditionalID].bvElseID;
		Lit condition =  all_conditionals[conditionalID].l;
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

	void buildConditionalReason(int bvID,int conditionalID, vec<Lit> & conflict){
		dbg_no_pending_analyses();
		assert(eq_bitvectors[bvID]==bvID);

		stats_build_condition_reason++;
		Weight  over_cur = over_approx[bvID];
		Weight  under_cur = under_approx[bvID];
		//assert(checkApproxUpToDate(bvID));


		//the reason that the addition is over is the reason that

		int bvThenID = all_conditionals[conditionalID].bvThenID;
		int bvElseID = all_conditionals[conditionalID].bvElseID;
		Lit condition =  all_conditionals[conditionalID].l;
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
	void buildConditionalArgReason(int bvID, int conditionalID, vec<Lit> & conflict){
			dbg_no_pending_analyses();
			assert(eq_bitvectors[bvID]==bvID);


			stats_build_condition_arg_reason++;
			Weight  over_cur = over_approx[bvID];
			Weight  under_cur = under_approx[bvID];
			//assert(checkApproxUpToDate(bvID));

			Lit condition =   all_conditionals[conditionalID].l;
			assert(value(condition)!=l_Undef);
			//the reason that the addition is over is the reason that
			//bvID > addition_under, or the reason that addition_under>= its current value.
			assert(all_conditionals[conditionalID].bvThenID == bvID || all_conditionals[conditionalID].bvElseID == bvID);
			int other_argID = all_conditionals[conditionalID].bvThenID == bvID ?  all_conditionals[conditionalID].bvElseID : all_conditionals[conditionalID].bvThenID;
			int resultID = all_conditionals[conditionalID].bvResultID;

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


	void buildAdditionArgReason(int bvID, int argindex, vec<Lit> & conflict){
		dbg_no_pending_analyses();
		assert(eq_bitvectors[bvID]==bvID);
		//rewind_trail_pos(trail.size()-1);

		stats_build_addition_arg_reason++;
		Weight  over_cur = over_approx[bvID];
		Weight  under_cur = under_approx[bvID];
		//assert(checkApproxUpToDate(bvID));


		//the reason that the addition is over is the reason that
		//bvID > addition_under, or the reason that addition_under>= its current value.
		int other_argID = addition_arguments[bvID][argindex].other_argID;
		int sumID = addition_arguments[bvID][argindex].sumID;


		Weight under_add = under_approx[sumID] -  over_approx[other_argID];
		Weight over_add = over_approx[sumID] -  under_approx[other_argID];

		int width = bitvectors[bvID].size();
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

	void buildAdditionReason(int bvID,int addition_index, vec<Lit> & conflict){
		dbg_no_pending_analyses();
		assert(eq_bitvectors[bvID]==bvID);
		//rewind_trail_pos(trail.size()-1);
		stats_build_addition_reason++;
		Weight  over_cur = over_approx[bvID];
		Weight  under_cur = under_approx[bvID];
		assert(checkApproxUpToDate(bvID));


		//the reason that the addition is over is the reason that
		//bvID > addition_under, or the reason that addition_under>= its current value.
		int aID = additions[bvID][addition_index].aID;
		int bID = additions[bvID][addition_index].bID;

		//assert(aID<bvID);
		//assert(bID<bvID);

		Weight under_add = under_approx[aID] +  under_approx[bID];
		Weight over_add = over_approx[aID] +  over_approx[bID];

		int width = bitvectors[bvID].size();
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

	void buildComparisonReasonBV(Comparison op, int bvID,int comparebvID, vec<Lit> & conflict,int trail_pos=-1){
		dbg_no_pending_analyses();
		assert(eq_bitvectors[bvID]==bvID);
		if (trail_pos<0){
			trail_pos=trail.size();
		}
		rewind_trail_pos(trail_pos);
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
		if(iter==363){
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
			//to = over_cur;
		}
		if (op==Comparison::leq){
			assert(over_cur<=to);
			compare_over=true;
			//to = over_cur;
		}
		if (op==Comparison::geq){
			assert(under_cur>=to);
			compare_over=false;
			//to = under_cur;
		}
		if (op==Comparison::gt){
			assert(under_cur>to);
			compare_over=false;
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
			if(++n_analysis==3237){
				int a =1;
			}
			std::stringstream ss;
			ss<< op << " " << to;
			fprintf(opt_write_learnt_clauses,"learnt analysis bv %d %s ",bvID,ss.str().c_str());
		}

		if(compare_over)
			assert(over_causes[bvID].hasCause());
		if(!compare_over)
			assert(under_causes[bvID].hasCause());

		if (compare_over &&  over_causes[bvID].cause_is_decision){
			//the reason that the bvID's over approx is <= its current value
			//is because it was a decision.
			//Create a literal on the fly to explain this...
			int lev = over_causes[bvID].index;
			assert(lev>=0);
			assert(lev<=decisionLevel());
			Lit reason = newComparison(Comparison::leq,bvID, over_cur,var_Undef,opt_cmp_lits_decidable);
			assert(value(reason)!=l_False);
			conflict.push(~reason);

			//S->prependToTrail(toSolver(reason),lev);//this is a decision that was made, without a corresponding literal in the solver at the time it was made.
			// need to  ensure that this lit can be properly analyzed, so prepend it to the trail at this decision level.

		}else if (!compare_over &&  under_causes[bvID].cause_is_decision){
			int lev = under_causes[bvID].index;
			assert(lev>=0);
			assert(lev<=decisionLevel());
			Lit reason = newComparison(Comparison::geq,bvID, under_cur,var_Undef,opt_cmp_lits_decidable);
			assert(value(reason)!=l_False);
			conflict.push(~reason);

			//S->prependToTrail(toSolver(reason),lev);//this is a decision that was made, without a corresponding literal in the solver at the time it was made.
			//need to ensure that this lit can be properly analyzed, so prepend it to the trail at this decision level.

		}else if (compare_over &&  over_causes[bvID].refined_cause){
			//then the reason the underapprox is too large is because of the assignment to the bits
			//can this analysis be improved upon?
			for(int i =0;i<bv.size();i++){
				Lit bl = bv[i];
				if(value(bl)==l_False){

						assert(value(bl)==l_False);
						conflict.push(bl);

				}else if(value(bl)==l_True){

					assert(value(~bl)==l_False);
					conflict.push(~bl);

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
			addAnalysis(Comparison::leq,bvID,prev_weight,conflict);


			//buildValueReason(Comparison::leq,bvID,over_approx[bvID],conflict,trail_pos-1);

		}else if (!compare_over  && under_causes[bvID].refined_cause){
			//then the reason the underapprox is too large is because of the assignment to the bits

			for(int i =0;i<bv.size();i++){
				Lit bl = bv[i];
				lbool val = value(bl);
				lbool dbgval = dbg_value(bl);
				if(value(bl)==l_False){

						assert(value(bl)==l_False);
						conflict.push(bl);

				}else if(value(bl)==l_True){

					assert(value(~bl)==l_False);
					conflict.push(~bl);

				}
			}

			int trail_pos = analysis_trail_pos;
			rewindUntil(bvID,Comparison::geq, under_approx[bvID]);
			assert(trail[analysis_trail_pos].isBoundAssignment());
			assert(trail[analysis_trail_pos].bvID==bvID);
			Weight prev_weight = trail[analysis_trail_pos].previous_under;
			rewind_trail_pos(trail_pos);
			addAnalysis(Comparison::geq,bvID,prev_weight,conflict);

			//buildValueReason(Comparison::geq,bvID,under_approx[bvID],conflict,trail_pos-1);

		}else if (compare_over &&  over_causes[bvID].cause_is_bits){
			//then the reason the underapprox is too large is because of the assignment to the bits
			Weight over = over_approx[bvID];
			for(int i =0;i<bv.size();i++){
				Lit bl = bv[i];
				if(value(bl)==l_False ){
					Weight bit = 1<<i;
					if(comp(op,over+bit,to)&& level(var(bl))>0){
						//then we can skip this bit, because we would still have had a conflict even if it was assigned true.
						over+=bit;
					}else{
						assert(value(bl)==l_False);
						conflict.push(bl);
					}
				}
			}

		}else if (!compare_over  && under_causes[bvID].cause_is_bits){
			//then the reason the underapprox is too large is because of the assignment to the bits
			Weight under = under_approx[bvID];
			for(int i =0;i<bv.size();i++){
				Lit bl = bv[i];
				lbool val = value(bl);
				lbool dbgval = dbg_value(bl);
				if(value(bl)==l_True){
					Weight bit = 1<<i;
					if(comp(op,under-bit,to)  && level(var(bl))>0){
						//then we can skip this bit, because we would still have had a conflict even if it was assigned false.
						under-=bit;
					}else{
						assert(value(~bl)==l_False);
						conflict.push(~bl);
					}
				}
			}


		}else if(compare_over && over_causes[bvID].cause_is_addition){
			int index = over_causes[bvID].index;
			assert(index>=0);
			assert(index< additions[bvID].size());
			int aID = additions[bvID][index].aID;
			int bID = additions[bvID][index].bID;
			//assert(aID<bvID);
			//assert(bID<bvID);
			Weight over_bid = over_approx[bID];
			Weight over_aid = over_approx[aID];
			//then the reason is that aID is <= weight-under(bID), or bID <= weight-under(aID)
			//addAnalysis(Comparison::leq,aID,over_aid,conflict);//to-over_bid
			//addAnalysis(Comparison::leq,bID,over_bid,conflict);//to-over_aid

			addAnalysis(Comparison::leq,aID,over_approx[bvID]-over_bid,conflict);
			addAnalysis(Comparison::leq,bID,over_approx[bvID]-over_aid,conflict);

		}else if(!compare_over && under_causes[bvID].cause_is_addition){
			int index = under_causes[bvID].index;
			assert(index>=0);
			assert(index< additions[bvID].size());
			int aID = additions[bvID][index].aID;
			int bID = additions[bvID][index].bID;
			Weight under_bid = under_approx[bID];
			Weight under_aid = under_approx[aID];
			//addAnalysis(Comparison::geq,aID,under_aid,conflict);
			//addAnalysis(Comparison::geq,bID,under_bid,conflict);
			addAnalysis(Comparison::geq,aID,under_approx[bvID]-under_bid,conflict);
			addAnalysis(Comparison::geq,bID,under_approx[bvID]-under_aid,conflict);
			//buildValueReason(op,aID,to-under_bid,conflict,trail_pos-1);
			//buildValueReason(op,bID,to-under_aid,conflict,trail_pos-1);


		}else if (compare_over && over_causes[bvID].cause_is_addition_argument){
			int argindex = over_causes[bvID].index;
			int other_argID = addition_arguments[bvID][argindex].other_argID;
			int sumID = addition_arguments[bvID][argindex].sumID;
			Weight over_sumID = over_approx[sumID];
			Weight under_argID = under_approx[other_argID];

			//Weight over = over_approx[sumID] -  under_approx[other_argID];

			addAnalysis(Comparison::geq,other_argID,over_sumID-over_approx[bvID],conflict);
			addAnalysis(Comparison::leq,sumID,over_approx[bvID]+under_argID,conflict);
			//buildValueReason(~op,other_argID,over_sumID-to,conflict,trail_pos-1);
			//buildValueReason(op,sumID,to+under_argID,conflict,trail_pos-1);

		}else if (!compare_over && under_causes[bvID].cause_is_addition_argument){
			int argindex = under_causes[bvID].index;
			int other_argID = addition_arguments[bvID][argindex].other_argID;
			int sumID = addition_arguments[bvID][argindex].sumID;
			Weight under_sumID = under_approx[sumID];
			Weight over_argID = over_approx[other_argID];
			//Weight under = under_approx[sumID] -  over_approx[other_argID];
			addAnalysis(Comparison::leq,other_argID,under_sumID-under_approx[bvID],conflict);
			addAnalysis(Comparison::geq,sumID,under_approx[bvID]+over_argID,conflict);
			//buildValueReason(~op,other_argID,under_sumID-to,conflict,trail_pos-1);
			//buildValueReason(op,sumID,to+over_argID,conflict,trail_pos-1);

		}else if(compare_over && over_causes[bvID].cause_is_condition){


			int index = over_causes[bvID].index;
			assert(index>=0);
			assert(index< conditionals[bvID].size());

			assert( conditionals[bvID][index].bvResultID==bvID);
			int bvThenID = conditionals[bvID][index].bvThenID ;
			int bvElseID = conditionals[bvID][index].bvElseID;
			Lit condition = conditionals[bvID][index].l;
			//if the condition is assigned (at this point in the trail)
			if (value(condition)==l_True ){
				conflict.push(~condition);
				assert(over_approx[bvThenID]<=over_approx[bvID]);//this is the cause
				addAnalysis(op,bvThenID,to,conflict);
			}else if (value(condition)==l_False ){
				conflict.push(condition);
				assert(over_approx[bvElseID]<=over_approx[bvID]);//this is the cause
				addAnalysis(op,bvElseID,to,conflict);
			}else{//else if (over_approx[bvThenID]>=over_approx[bvElseID]){
				assert(over_approx[bvThenID]<=over_approx[bvID]);
				assert( over_approx[bvElseID]<=over_approx[bvID]);//this is the cause
				addAnalysis(op,bvThenID,to,conflict);
				addAnalysis(op,bvElseID,to,conflict);
			}
		}else if(!compare_over && under_causes[bvID].cause_is_condition){
			int index = under_causes[bvID].index;
			assert(index>=0);
			assert(index< conditionals[bvID].size());

			assert( conditionals[bvID][index].bvResultID==bvID);
			int bvThenID = conditionals[bvID][index].bvThenID ;
			int bvElseID = conditionals[bvID][index].bvElseID;
			Lit condition = conditionals[bvID][index].l;
			//if the condition is assigned (at this point in the trail)
			if (value(condition)==l_True ){
				conflict.push(~condition);
				assert(under_approx[bvThenID]>=under_approx[bvID]);//this is the cause
				addAnalysis(op,bvThenID,to,conflict);
			}else if (value(condition)==l_False){
				conflict.push(condition);
				assert(under_approx[bvElseID]>=under_approx[bvID]);//this is the cause
				addAnalysis(op,bvElseID,to,conflict);
			}else{//if ( under_approx[bvThenID]<=under_approx[bvElseID]){
				assert(under_approx[bvThenID]>=under_approx[bvID]);
				assert(under_approx[bvElseID]>=under_approx[bvID]);
				addAnalysis(op,bvThenID,to,conflict);
			//}else if ( under_approx[bvElseID]<under_approx[bvThenID]){

				addAnalysis(op,bvElseID,to,conflict);
			}
		}else if (compare_over && over_causes[bvID].cause_is_condition_argument){
			int argindex = over_causes[bvID].index;
			assert(argindex< conditional_arguments[bvID].size());
			Lit condition =  conditional_arguments[bvID][argindex].l;
			int other_argID = conditional_arguments[bvID][argindex].other_argID;
			int resultID = conditional_arguments[bvID][argindex].resultID;

			Weight over_sum = over_approx[resultID];
			assert(over_sum<=over_approx[bvID]);
			assert(value(condition)==l_True);

			conflict.push(~condition);
			addAnalysis(op,resultID,to,conflict);

		}else if (!compare_over && under_causes[bvID].cause_is_condition_argument){
			int argindex = under_causes[bvID].index;
			assert(argindex< conditional_arguments[bvID].size());
			Lit condition =  conditional_arguments[bvID][argindex].l;
			int other_argID = conditional_arguments[bvID][argindex].other_argID;
			int resultID = conditional_arguments[bvID][argindex].resultID;

			Weight under_sum = under_approx[resultID];
			assert(under_sum>=under_approx[bvID]);
			assert(value(condition)==l_True);

			conflict.push(~condition);
			addAnalysis(op,resultID,to,conflict);

		}

		int cID=-1;
		if (compare_over && over_causes[bvID].cause_is_comparison){
			cID= over_causes[bvID].index;
		}else if (!compare_over && under_causes[bvID].cause_is_comparison){
			cID= under_causes[bvID].index;
		}
		if(cID>-1){
			ComparisonID & c = comparisons[cID];
			Comparison cop = c.op();//invert this because we are switch the direction of comparison
			assert(c.l!=lit_Undef);
			if(value(c.l)==l_True){
				assert(value(~c.l)==l_False);
				conflict.push(~c.l);
				if(c.bvCompare()){
					if (compare_over){
						addAnalysis(Comparison::leq,c.compareID, over_approx[c.compareID],conflict);
						//buildValueReason(Comparison::leq,c.compareID, over_approx[c.compareID],conflict,trail_pos-1);
					}else{
						addAnalysis(Comparison::geq,c.compareID, under_approx[c.compareID],conflict);
						//buildValueReason(Comparison::geq,c.compareID, under_approx[c.compareID],conflict,trail_pos-1);
					}
				}
			}else{
				assert(value(c.l)==l_False);
				conflict.push(c.l);
				if(c.bvCompare()){
					if (compare_over){
						addAnalysis(Comparison::leq,c.compareID, over_approx[c.compareID],conflict);
					}else{
						addAnalysis(Comparison::geq,c.compareID, under_approx[c.compareID],conflict);
					}
				}
			}
		}
		if(opt_write_learnt_clauses && opt_write_bv_analysis){
			fprintf(opt_write_learnt_clauses," lits " );
			for(int i = 0;i<conflict.size();i++){
				fprintf(opt_write_learnt_clauses,"%d ",dimacs(toSolver(conflict[i])));
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

	bool addAnalysis(Comparison op, int bvID, Weight  to,  vec<Lit> & conflict){

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
						fprintf(opt_write_learnt_clauses, " bv %d %s ",bvID,ss.str().c_str());
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
						fprintf(opt_write_learnt_clauses, " bv %d %s ",bvID,ss.str().c_str());
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
			fprintf(opt_write_learnt_clauses, " bv %d %s ",bvID,ss.str().c_str());
		}
		return true;
	}

	void buildComparisonReason(Comparison op, int bvID, Weight  to,  vec<Lit> & conflict){
		dbg_no_pending_analyses();
		rewind_trail_pos(trail.size()-1);
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

	bool check_solved() {
		for(int bvID = 0;bvID<bitvectors.size();bvID++){
			int eqBV = bvID;
			while(eq_bitvectors[eqBV]>-1 && eq_bitvectors[eqBV]!=eqBV ){
				eqBV=eq_bitvectors[eqBV];
			}

			if(over_approx[eqBV]!= under_approx[eqBV]){
				return false;
			}


			vec<Lit> & bv = bitvectors[bvID];
			Weight over=0;
			Weight under=0;
			if(bv.size()){
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
			}else{
				over=over_approx[eqBV];
				under=under_approx[eqBV];
			}

			for(int i = 0;i< conditionals[bvID].size();i++){
				int thenID = conditionals[bvID][i].bvThenID;
				int elseID = conditionals[bvID][i].bvElseID;
				Lit condition = conditionals[bvID][i].l;
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
			}

			for(int i = 0;i< additions[bvID].size();i++){
				int aID = additions[bvID][i].aID;
				int bID = additions[bvID][i].bID;
				int width = bitvectors[bvID].size();
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
				if(underadd >under){
					return false;
				}
				if(overadd<over){
					return false;
				}
			}

			for(int cID:compares[bvID]){
				ComparisonID & c = comparisons[cID];
				assert(c.compareID<0);
				Comparison op = c.op();
				switch (op){
					case Comparison::lt:
						if(value(c.l)==l_True && under>= c.w){
							return false;
						}else if (value(c.l)==l_False && over<c.w){
							return false;
						}
						break;
					case Comparison::leq:
						if(value(c.l)==l_True && under> c.w){
							return false;
						}else if (value(c.l)==l_False && over<=c.w){
							return false;
						}
						break;
					case Comparison::gt:
						if(value(c.l)==l_True && over<= c.w){
							return false;
						}else if (value(c.l)==l_False && under>c.w){
							return false;
						}
						break;
					case Comparison::geq:
					default:
						if(value(c.l)==l_True && over < c.w){
							return false;
						}else if (value(c.l)==l_False && under>= c.w){
							return false;
						}
						break;
				}

			}

			for(int cID:bvcompares[bvID]){
				ComparisonID & c = comparisons[cID];
				assert(c.compareID>=0);
				Comparison op = c.op();
				int toID = c.compareID;
				vec<Lit> & bv_compare = bitvectors[toID];
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
						if(value(c.l)==l_True && under>= over_compare){
							return false;
						}else if (value(c.l)==l_False && over<under_compare){
							return false;
						}
						break;
					case Comparison::leq:
						if(value(c.l)==l_True && under> over_compare){
							return false;
						}else if (value(c.l)==l_False && over<=under_compare){
							return false;
						}
						break;
					case Comparison::gt:
						if(value(c.l)==l_True && over<= under_compare){
							return false;
						}else if (value(c.l)==l_False && under>over_compare){
							return false;
						}
						break;
					case Comparison::geq:
					default:
						if(value(c.l)==l_True && over < under_compare){
							return false;
						}else if (value(c.l)==l_False && under>= over_compare){
							return false;
						}
						break;
				}
			}
		}
		return true;
	}
	
	bool dbg_solved() {

		return true;
	}




	void setBitvectorTheory(int bvID, int theoryID){
		theoryIds[bvID]=theoryID;
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
		int conditionID = all_conditionals.size();//need this for correctness if this bv is merged with another, later...
		Var outerVar = var(outer_condition_lit);
		Lit condition = mkLit(newVar(outerVar, resultID,conditionID,true,true,true));

		int cID = all_conditionals.size();
		all_conditionals.push();
		all_conditionals.last().cID=cID;
		all_conditionals.last().l = condition;
		all_conditionals.last().bvResultID = resultID;
		all_conditionals.last().bvThenID = bvThenID;
		all_conditionals.last().bvElseID = bvElseID;

		conditionals[resultID].push();
		conditionals[resultID].last().cID=cID;
		conditionals[resultID].last().l = condition;
		conditionals[resultID].last().bvResultID = resultID;
		conditionals[resultID].last().bvThenID = bvThenID;
		conditionals[resultID].last().bvElseID = bvElseID;

		conditional_arguments[bvThenID].push();
		conditional_arguments[bvThenID].last().cID = cID;
		conditional_arguments[bvThenID].last().l = condition;
		conditional_arguments[bvThenID].last().other_argID=bvElseID;
		conditional_arguments[bvThenID].last().resultID=resultID;

		conditional_arguments[bvElseID].push();
		conditional_arguments[bvElseID].last().cID = cID;
		conditional_arguments[bvElseID].last().l = ~condition;
		conditional_arguments[bvElseID].last().other_argID=bvThenID;
		conditional_arguments[bvElseID].last().resultID=resultID;

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

		additions[resultID].push();
		additions[resultID].last().aID=aID;
		additions[resultID].last().bID=bID;

		addition_arguments[aID].push();
		addition_arguments[aID].last().other_argID=bID;
		addition_arguments[aID].last().sumID=resultID;

		addition_arguments[bID].push();
		addition_arguments[bID].last().other_argID=aID;
		addition_arguments[bID].last().sumID=resultID;

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
		additions.growTo(bvID+1);
		addition_arguments.growTo(bvID+1);
		conditionals.growTo(bvID+1);
		conditional_arguments.growTo(bvID+1);
		under_causes.growTo(bvID+1);
		over_causes.growTo(bvID+1);
		cause_set.growTo(bvID+1);
		eq_bitvectors.growTo(bvID+1,-1);
		eq_bitvectors[bvID]=bvID;
		pending_under_analyses.growTo(bvID+1,-1);
		pending_over_analyses.growTo(bvID+1,-1);

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


		for(int i = 0;i<vars.size();i++){
			bitvectors[bvID].push(mkLit(newVar(vars[i],bvID)));
		}

		bv_needs_propagation.growTo(bvID+1);
		bv_needs_propagation[bvID]=true;
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
			const_true = mkLit(newVar(var(S->True()),-1,-1,true,false),sign(S->True()));
			enqueueTheory(const_true);
		}
		return const_true;
	}

	BitVector newBitvector(int bvID, int bitwidth,Weight constval=-1, int equivalentBV=-1){
		if(bvID<0){
			bvID = nBitvectors();
		}
		if(bvID==55349){
			int a=1;
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
		additions.growTo(bvID+1);
		addition_arguments.growTo(bvID+1);
		conditionals.growTo(bvID+1);
		conditional_arguments.growTo(bvID+1);
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
			else
				return ~comparisons[cID].l;
		}else
			return comparisons[cID].l;
	}
	int getComparisonID(Comparison op, int bvID,const Weight & w){
		//could do a binary search here:
#ifndef NDEBUG
		int expect = -1;
		for(int i=0;i<compares[bvID].size();i++){
			int cID = compares[bvID][i];
			if (comparisons[cID].op() == op && comparisons[cID].w == w){
				expect= cID;
				break;
			}
		}
#endif
		vec<int> & compare = compares[bvID];
		if(compare.size()){
			dbg_compares_sorted(bvID);
			int pos = binary_search_Weight(compare,w);
			if(pos>=0 && pos < compare.size() && comparisons[compare[pos]].w==w){
				//we found a comparison to the same bitvector. Now lets see if there is a comparison with the same operator to that bitvector
				while(pos< compare.size() && comparisons[compare[pos]].w == w){
					if( comparisons[compare[pos]].op()==op){
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
			assert(comparisons[cID0].w <= comparisons[cID1].w);
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
			assert(comparisons[cID0].compareID <= comparisons[cID1].compareID);
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
			if (comparisons[cID].w == w )
			{
				//ensure we get the first such matching CID
				while(midpoint>0 && comparisons[compares[midpoint-1]].w == w){
					midpoint--;
				}
				return midpoint;
			}
			else if (w <comparisons[cID].w)
				high = midpoint - 1;
			else
				low = midpoint + 1;
		}
		//this can probably be done more cleanly...
		if (compares.size()==0)
			return -1;
		else if (low<compares.size() && comparisons[compares[low]].w > w){
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
			if (comparisons[cID].compareID == compareID )
			{
				//ensure we get the first such matching CID
				while(midpoint>0 && comparisons[compares[midpoint-1]].compareID == compareID){
					midpoint--;
				}
				return midpoint;
			}
			else if (compareID <comparisons[cID].compareID)
				high = midpoint - 1;
			else
				low = midpoint + 1;
		}
		if (compares.size()==0)
			return -1;
		else if (low<compares.size() && comparisons[compares[low]].compareID > compareID){
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
			if (comparisons[cID].compareID == compareID && comparisons[cID].op()==op){
				expect = comparisons[cID].l;
				break;
			}
		}
#endif
		if(bvcompare.size()){
			int pos = binary_search_CID(bvcompare,compareID);
			if(pos>=0 && pos < bvcompare.size() && comparisons[bvcompare[pos]].compareID==compareID){
				//we found a comparison to the same bitvector. Now lets see if there is a comparison with the same operator to that bitvector
				while(pos< bvcompare.size() && comparisons[bvcompare[pos]].compareID == compareID){
					if( comparisons[bvcompare[pos]].op()==op){
						assert(comparisons[bvcompare[pos]].l==expect);
						return comparisons[bvcompare[pos]].l;
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

		int comparisonID = comparisons.size();
		if(comparisonID==8){
			int a=1;
		}
		if((l = getComparison(op, bvID, to))!=lit_Undef){
			if(outerVar != var_Undef){
				makeEqualInSolver(mkLit(outerVar),toSolver(l));
			}
			return l;
		}else{
			l = mkLit(newVar(outerVar, bvID,comparisonID,true,decidable));

		}
		static int iter = 0;
		if(++iter==29){
			int a =1;
		}

		if(opt_write_learnt_clauses){
			std::stringstream ss;
			ss << op << " " << to ;
			fprintf(opt_write_learnt_clauses,"learnt fact bv %d %s %d\n", bvID, ss.str().c_str(),dimacs(toSolver(l)) );
			fflush(opt_write_learnt_clauses);
			//std::cout << "learnt fact " << "bv " << bvID <<  <<" " << dimacs(toSolver(l)) << "\n";
		}



		//updateApproximations(bvID);
		comparisons.push(ComparisonID(to,-1,l,bvID,op));


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


		comparison_needs_repropagation.growTo(comparisons.size());
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
							if (over_causes[bvID].cause_is_decision && overApprox==to-1 ){

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

							if (under_causes[bvID].cause_is_decision && underApprox==to ){
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
							if (over_causes[bvID].cause_is_decision && overApprox==to ){
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
							if (under_causes[bvID].cause_is_decision && underApprox==to+1 ){
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
							if (over_causes[bvID].cause_is_decision && overApprox==to ){
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
							if (under_causes[bvID].cause_is_decision && underApprox==to+1 ){
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

							if (over_causes[bvID].cause_is_decision && overApprox==to-1 ){
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
							if (under_causes[bvID].cause_is_decision && underApprox==to ){
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
		int comparisonID = comparisons.size();
		if((l = getComparisonBV(op,bvID, toID))!=lit_Undef){
			if(outerVar != var_Undef){
				makeEqualInSolver(mkLit(outerVar),toSolver(l));
			}
			return l;
		}else{
			l = mkLit(newVar(outerVar, bvID,comparisonID));
		}

		//updateApproximations(bvID);
		//updateApproximations(toID);

	/*	if(!cause_set[toID].contains(bvID))
			cause_set[toID].push(bvID);*/
		dbg_bvcompares_sorted(bvID);
		comparisons.push(ComparisonID(-1,toID,l,bvID,op));

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
		dbg_bvcompares_sorted(bvID);

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
		dbg_bvcompares_sorted(toID);
		comparisonID = comparisons.size();
		comparisons.push(ComparisonID(-1,bvID,l,toID,~op));

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
		comparison_needs_repropagation.growTo(comparisons.size());
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

						enqueueEager(l,bvID,underApprox,overApprox,prev_under_cause,prev_over_cause, comparisonprop_marker);
					}
				}

				if (underApprox>=overCompare){
					if(value(l)==l_True){
						assert(false);
					}else if (value(l)==l_False){
						//do nothing
					}else {
						assert(value(l)==l_Undef);
						enqueueEager(~l,bvID,underApprox,overApprox,prev_under_cause,prev_over_cause, comparisonprop_marker);
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
						enqueueEager(l,bvID,underApprox,overApprox,prev_under_cause,prev_over_cause, comparisonprop_marker);
					}
				}

				if (underApprox>overCompare){
					if(value(l)==l_True){
						assert(false);
					}else if (value(l)==l_False){
						//do nothing
					}else {
						assert(value(l)==l_Undef);
						enqueueEager(~l,bvID,underApprox,overApprox,prev_under_cause,prev_over_cause, comparisonprop_marker);
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
						enqueueEager(~l,bvID,underApprox,overApprox,prev_under_cause,prev_over_cause, comparisonprop_marker);
					}
				}

				if (underApprox>overCompare){
					if(value(l)==l_True){

					}else if (value(l)==l_False){
						assert(false);
					}else {
						assert(value(l)==l_Undef);
						enqueueEager(l,bvID,underApprox,overApprox,prev_under_cause,prev_over_cause, comparisonprop_marker);
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
						enqueueEager(~l,bvID,underApprox,overApprox,prev_under_cause,prev_over_cause, comparisonprop_marker);
					}
				}

				if (underApprox>=overCompare){
					if(value(l)==l_True){

					}else if (value(l)==l_False){
						assert(false);
					}else {
						assert(value(l)==l_Undef);
						enqueueEager(l,bvID,underApprox,overApprox,prev_under_cause,prev_over_cause, comparisonprop_marker);
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
			ComparisonID & c = comparisons[cID];
			Comparison op = c.op();

			lbool val = value(c.l);
			if(!c.bvCompare()){
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
			}else{


			}
		}
		for(int cID:bvcompares[bvID]){
			ComparisonID & c = comparisons[cID];
			Comparison op = c.op();

			lbool val = value(c.l);
			if(c.bvCompare()){
				Weight under_w = under_approx[c.compareID];
				Weight over_w = over_approx[c.compareID];

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
			}else{


			}
		}
		}
#endif
		return true;
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
