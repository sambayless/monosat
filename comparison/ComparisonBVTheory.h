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

#ifndef COMPARISON_THEORY_H_
#define COMPARISON_THEORY_H_

#include "utils/System.h"
#include "core/Theory.h"

#include "core/SolverTypes.h"
#include "mtl/Map.h"


#include "utils/System.h"
#include "core/Solver.h"

#include <vector>

#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <sstream>

using namespace dgl;
namespace Monosat {


template<typename Weight>
class ComparisonBVTheorySolver: public Theory {
public:

	struct Assignment {
		bool isEdge :1;
		bool assign :1;
		int edgeID:30;
		Var var;
		Assignment(bool isEdge, bool assign, int edgeID, Var v):isEdge(isEdge),assign(assign),edgeID(edgeID),var(v){

		}
	};

	struct Comparison{
		Weight w;
		Lit l;

		int bvID:31;
		int is_lt:1;
	};

	double rnd_seed;

private:
	Solver * S;
	int local_q = 0;
public:
	int id;

	vec<lbool> assigns;
	CRef comparisonprop_marker;
	CRef bvprop_marker;

	vec<Assignment> trail;
	vec<int> trail_lim;

	//Bitvectors are unsigned, and have least significant bit at index 0
	vec<vec<Lit> > bitvectors;
	vec<Comparison> comparisons;
	vec<vec<int> > comparisons_lt; //for each bitvector, comparisons are all to unique values, and in ascending order of compareTo.
	vec<vec<int> > comparisons_gt; //for each bitvector, comparisons are all to unique values, and in ascending order of compareTo.

	vec<Weight> under_approx;
	vec<Weight> over_approx;

	vec<int> altered_bvs;
	vec<bool> alteredBV;
public:


	vec<int> marker_map;

	bool requiresPropagation = true;

	vec<char> seen;
	vec<int> to_visit;
	vec<Lit> tmp_clause;
	//Data about local theory variables, and how they connect to the sat solver's variables
	struct VarData {

		int occursPositive :1;
		int occursNegative :1;

		int bvID:30;
		int comparisonID;
		Var solverVar;
	};

	vec<VarData> vars;
	int theory_index = 0;
public:

	double reachtime = 0;
	double unreachtime = 0;
	double pathtime = 0;
	double propagationtime = 0;
	long stats_propagations = 0;
	long stats_num_conflicts = 0;
	long stats_decisions = 0;
	long stats_num_reasons = 0;

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
	long stats_propagations_skipped = 0;


	ComparisonBVTheorySolver(Solver * S_, int _id = -1) :
			S(S_), id(_id){
		rnd_seed = opt_random_seed;


		S->addTheory(this);
		comparisonprop_marker = S->newReasonMarker(this);
		 bvprop_marker = S->newReasonMarker(this);

	}
	~ComparisonBVTheorySolver(){

	}
	Solver * getSolver(){
		return S;
	}
	void printStats(int detailLevel) {

	}
	
	void writeTheoryWitness(std::ostream& write_to) {

	}
	
	inline int getTheoryIndex() {
		return theory_index;
	}
	inline void setTheoryIndex(int id) {
		theory_index = id;
	}


	inline bool isComparisonVar(Var v) const{
		assert(v < vars.size());
		return vars[v].comparisonID>-1;
	}
	inline int getComparisonID(Var v) const {
		assert(isComparisonVar(v));
		return vars[v].comparisonID;
	}
	inline Comparison& getComparison(int comparisonID){
		return comparisons[comparisonID];
	}
	inline int getbvID(Var v) const{

		return vars[v].bvID;
	}

	inline int getDetector(Var v) const {
		assert(!isComparisonVar(v));
		return vars[v].detector_edge;
	}
	

	
	void makeEqual(Lit l1, Lit l2) {
		Lit o1 = toSolver(l1);
		Lit o2 = toSolver(l2);
		S->addClause(~o1, o2);
		S->addClause(o1, ~o2);
	}
	void makeEqualInSolver(Lit l1, Lit l2) {
		S->addClause(~l1, l2);
		S->addClause(l1, ~l2);
	}
	void addClause(Lit l1) {
		Lit o1 = toSolver(l1);
		S->addClause(o1);
	}
	void addClause(Lit l1, Lit l2) {
		Lit o1 = toSolver(l1);
		Lit o2 = toSolver(l2);
		S->addClause(o1, o2);
	}
	void addClause(Lit l1, Lit l2, Lit l3) {
		Lit o1 = toSolver(l1);
		Lit o2 = toSolver(l2);
		Lit o3 = toSolver(l3);
		S->addClause(o1, o2, o3);
	}
	void addClause(vec<Lit> & c) {
		tmp_clause.clear();
		c.copyTo(tmp_clause);
		toSolver(tmp_clause);
		S->addClause(tmp_clause);
	}
	void addClauseSafely(vec<Lit> & c) {
		tmp_clause.clear();
		c.copyTo(tmp_clause);
		toSolver(tmp_clause);
		
		S->addClauseSafely(tmp_clause);
	}
	

	Var newVar(Var solverVar=var_Undef, int bvID=-1, int comparisonID=-1, bool connectToTheory = true) {
		if(solverVar==var_Undef){
			solverVar = S->newVar();
		}
		while (S->nVars() <= solverVar)
			S->newVar();
		Var v = vars.size();

		vars.push();
		vars[v].occursPositive=false;
		vars[v].occursNegative=false;
		vars[v].bvID=bvID;
		vars[v].comparisonID=comparisonID;
		vars[v].solverVar = solverVar;

		assigns.push(l_Undef);
		if (connectToTheory) {
			S->setTheoryVar(solverVar, getTheoryIndex(), v);
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
	inline bool enqueue(Lit l, CRef reason) {
		assert(assigns[var(l)]==l_Undef);
		
		Lit sl = toSolver(l);
		if (S->enqueue(sl, reason)) {
			enqueueTheory(l);
			return true;
		} else {
			return false;
		}
	}




	void backtrackUntil(int level) {
		static int it = 0;
		
		bool changed = false;
		//need to remove and add edges in the two graphs accordingly.
		if (trail_lim.size() > level) {
			
			int stop = trail_lim[level];
			for (int i = trail.size() - 1; i >= trail_lim[level]; i--) {
				
				Assignment & e = trail[i];
				assert(assigns[e.var]!=l_Undef);
				if (e.isEdge) {
					assert(dbg_value(e.var)==l_Undef);


					if (e.assign) {
						//g_unders[bvID]->disableTransition(edgeID, input,output);

					} else {
						//g_overs[bvID]->enableTransition(edgeID,input,output);

					}
				} else {
					//This is a reachability literal				  
					//detectors[getDetector(e.var)]->unassign(mkLit(e.var, !e.assign));
				}
				assigns[e.var] = l_Undef;
				//changed = true;
			}
			trail.shrink(trail.size() - stop);
			trail_lim.shrink(trail_lim.size() - level);
			assert(trail_lim.size() == level);
			
			if (changed) {
				//requiresPropagation = true;

			}
			

		}

		
	};

	Lit decideTheory() {

		return lit_Undef;
	}
	
	void backtrackUntil(Lit p) {
		//need to remove and add edges in the two graphs accordingly.
		assert(value(p)==l_True);
		int i = trail.size() - 1;
		for (; i >= 0; i--) {
			Assignment e = trail[i];
			if (var(p) == e.var) {
				assert(sign(p) != e.assign);
				break;
			}
			if (e.isEdge) {

				assert(assigns[e.var]!=l_Undef);
				assigns[e.var] = l_Undef;
				if (e.assign) {
				//	g_unders[bvID]->disableTransition(edgeID,input,output);
				} else {
				//	g_overs[bvID]->enableTransition(edgeID,input,output);
				}
			} else {

				assigns[e.var] = l_Undef;
				//detectors[getDetector(e.var)]->unassign(mkLit(e.var, !e.assign));
			}
		}
		
		trail.shrink(trail.size() - (i + 1));
		//if(i>0){
		requiresPropagation = true;

	}
	;

	void newDecisionLevel() {
		trail_lim.push(trail.size());
	}
	;

	void buildReason(Lit p, vec<Lit> & reason) {
		CRef marker = S->reason(var(toSolver(p)));
		assert(marker != CRef_Undef);
		int pos = CRef_Undef - marker;
		int d = marker_map[pos];
		//double initial_start = rtime(1);
		double start = rtime(1);
		backtrackUntil(p);
		
		if (marker == comparisonprop_marker) {
			reason.push(p);
			Var v = var(p);
			int bvID = getbvID(v);
			if(isComparisonVar(v)){
				int comparisonID = getComparisonID(v);
				if(comparisons[comparisonID].is_lt){
					if(!sign(p)){
						buildValueLTReason(bvID,comparisonID,reason);
					}else{
						buildValueGEQReason(bvID,comparisonID,reason);
					}
				}else{
					if(!sign(p)){
						buildValueGTReason(bvID,comparisonID,reason);
					}else{
						buildValueLEQReason(bvID,comparisonID,reason);
					}
				}
			}else{

			}

		} else if (marker == bvprop_marker) {
			reason.push(p);
			Var v = var(p);

		}  else {
			assert(false);
		}

		toSolver(reason);
		double finish = rtime(1);
		stats_reason_time += finish - start;
		stats_num_reasons++;
		//stats_reason_initial_time+=start-initial_start;
		
	}


	void preprocess() {

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
		
		int lev = level(v);
		
		assert(decisionLevel() <= lev);
		
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
		

		
		if (!isComparisonVar(var(l))) {
			

			trail.push( { false, !sign(l),-1, v });
			int bvID = getbvID(v);
			if(!alteredBV[bvID]){
				alteredBV[bvID]=true;
				altered_bvs.push(bvID);
			}

		} else {

			int bvID = getbvID(var(l));
			int edgeID = getComparisonID(var(l)); //v-min_edge_var;


			//trail.push( { true, !sign(l),edgeID, v });

		}
	}
	;

	void updateApproximations(int bvID){
		vec<Lit> & bv = bitvectors[bvID];
		under_approx[bvID]=0;
		over_approx[bvID]=0;
		for(int i = 0;i<bv.size();i++){
			lbool val = value(bv[i]);
			if(val==l_True){
				Weight bit = 1<<i;
				under_approx[bvID]+=bit;
				over_approx[bvID]+=bit;
			}else if (val==l_False){

			}else{
				Weight bit = 1<<i;
				over_approx[bvID]+=bit;
			}
		}
	}

	bool checkApproxUpToDate(int bvID){
#ifndef NDEBUG
		vec<Lit> & bv = bitvectors[bvID];
		Weight under =0;
		Weight over=0;
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
		assert(under==under_approx[bvID]);
		assert(over==over_approx[bvID]);
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

	bool propagateTheory(vec<Lit> & conflict) {

		stats_propagations++;
		assert(dbg_synced());
		if (!requiresPropagation) {
			stats_propagations_skipped++;

			return true;
		}
		
		bool any_change = false;
		double startproptime = rtime(1);
		//static vec<int> detectors_to_check;
		
		conflict.clear();

		
		while(altered_bvs.size()){
			int bvID = altered_bvs.last();
			assert(alteredBV[bvID]);

			updateApproximations(bvID);

			Weight & underApprox = under_approx[bvID];
			Weight & overApprox = over_approx[bvID];

			vec<int> & compares_lt = comparisons_lt[bvID];
			//update over approx lits
			for(int i = 0;i<compares_lt.size();i++){
				int comparisonID = compares_lt[i];
				Weight & lt = getComparison(comparisonID).w;
				Lit l =  getComparison(comparisonID).l;
				if (overApprox<lt){
					if(value(l)==l_True){
						//do nothing
					}else if (value(l)==l_False){
						assert(value(l)==l_False);
						assert(dbg_value(l)==l_False);
						conflict.push(l);
						buildValueLTReason(bvID,comparisonID,conflict);
						toSolver(conflict);
						return false;
					}else {
						assert(value(l)==l_Undef);
						enqueue(l, comparisonprop_marker);
					}

				}else{
					break;
				}
			}

			//update under approx lits
			for(int i = compares_lt.size()-1;i>=0;i--){
				int comparisonID = compares_lt[i];
				Weight & lt = getComparison(comparisonID).w;
				Lit l =  getComparison(comparisonID).l;
				if (underApprox>=lt){
					if(value(l)==l_True){

						conflict.push(~l);
						buildValueGEQReason(bvID,comparisonID,conflict);
						toSolver(conflict);
						return false;
					}else if (value(l)==l_False){
						//do nothing
					}else {
						assert(value(l)==l_Undef);
						enqueue(~l, comparisonprop_marker);
					}

				}else{
					break;
				}
			}

			vec<int> & compares_gt = comparisons_gt[bvID];
			//update over approx lits
			for(int i = 0;i<compares_gt.size();i++){
				int comparisonID = compares_gt[i];
				Weight & gt = getComparison(comparisonID).w;
				Lit l =  getComparison(comparisonID).l;
				if (overApprox<=gt){
					if(value(l)==l_True){
						conflict.push(~l);
						buildValueLEQReason(bvID,comparisonID,conflict);
						toSolver(conflict);
						return false;
					}else if (value(l)==l_False){

					}else {
						assert(value(l)==l_Undef);
						enqueue(~l, comparisonprop_marker);
					}

				}else{
					break;
				}
			}

			//update under approx lits
			for(int i = compares_gt.size()-1;i>=0;i--){
				int comparisonID = compares_gt[i];
				Weight & gt = getComparison(comparisonID).w;
				Lit l =  getComparison(comparisonID).l;
				if (underApprox>gt){
					if(value(l)==l_True){

					}else if (value(l)==l_False){
						conflict.push(l);
						buildValueGTReason(bvID,comparisonID,conflict);
						toSolver(conflict);
						return false;
					}else {
						assert(value(l)==l_Undef);
						enqueue(l, comparisonprop_marker);
					}

				}else{
					break;
				}
			}
			alteredBV[bvID]=false;
			altered_bvs.pop();
		}
		


		
		requiresPropagation = false;

		double elapsed = rtime(1) - startproptime;
		propagationtime += elapsed;

		return true;
	};

	void buildValueLTReason(int bvID, int comparisonID, vec<Lit> & conflict){

		Weight & lt = getComparison(comparisonID).w;
		Lit l =  getComparison(comparisonID).l;

		vec<Lit> & bv = bitvectors[bvID];
		Weight  over = over_approx[bvID];
		assert(checkApproxUpToDate(bvID));
		assert(over<lt);
		//the reason that the value is less than the weight 'lt' is that the _overapprox_ of the weight is less than lt.
		for(int i =0;i<bv.size();i++){
			Lit bl = bv[i];
			if(value(bl)==l_False && level(var(bl))>0){
				Weight bit = 1<<i;
				if(over+bit<lt){
					//then we can skip this bit, because we would still have had a conflict even if it was assigned true.
					over+=bit;
				}else{
					conflict.push(bl);
				}
			}
		}

	}

	void buildValueGEQReason(int bvID, int comparisonID, vec<Lit> & conflict){
		Weight & lt = getComparison(comparisonID).w;
		Lit l =  getComparison(comparisonID).l;

		vec<Lit> & bv = bitvectors[bvID];
		Weight  under = under_approx[bvID];
		assert(checkApproxUpToDate(bvID));
		assert(under>=lt);
		//the reason that the value is less than the weight 'lt' is that the _overapprox_ of the weight is less than lt.

		for(int i =0;i<bv.size();i++){
			Lit bl = bv[i];
			if(value(bl)==l_True && level(var(bl))>0){
				Weight bit = 1<<i;
				if(under-bit>=lt){
					//then we can skip this bit, because we would still have had a conflict even if it was assigned false.
					under-=bit;
				}else{
					conflict.push(~bl);
				}
			}
		}
	}

	void buildValueLEQReason(int bvID, int comparisonID, vec<Lit> & conflict){

			Weight & lt = getComparison(comparisonID).w;
			Lit l =  getComparison(comparisonID).l;

			vec<Lit> & bv = bitvectors[bvID];
			Weight  over = over_approx[bvID];
			assert(checkApproxUpToDate(bvID));
			assert(over<=lt);
			//the reason that the value is less than the weight 'lt' is that the _overapprox_ of the weight is less than lt.
			for(int i =0;i<bv.size();i++){
				Lit bl = bv[i];
				if(value(bl)==l_False && level(var(bl))>0){
					Weight bit = 1<<i;
					if(over+bit<=lt){
						//then we can skip this bit, because we would still have had a conflict even if it was assigned true.
						over+=bit;
					}else{
						conflict.push(bl);
					}
				}
			}

		}

		void buildValueGTReason(int bvID, int comparisonID, vec<Lit> & conflict){
			Weight & lt = getComparison(comparisonID).w;
			Lit l =  getComparison(comparisonID).l;

			vec<Lit> & bv = bitvectors[bvID];
			Weight  under = under_approx[bvID];
			assert(checkApproxUpToDate(bvID));
			assert(under>lt);
			//the reason that the value is less than the weight 'lt' is that the _overapprox_ of the weight is less than lt.

			for(int i =0;i<bv.size();i++){
				Lit bl = bv[i];
				if(value(bl)==l_True && level(var(bl))>0){
					Weight bit = 1<<i;
					if(under-bit>lt){
						//then we can skip this bit, because we would still have had a conflict even if it was assigned false.
						under-=bit;
					}else{
						conflict.push(~bl);
					}
				}
			}
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
			vec<Lit> & bv = bitvectors[bvID];
			Weight over=0;
			Weight under=0;
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
			for(int cID:comparisons_lt[bvID]){
				Comparison & c = comparisons[cID];
				if(value(c.l)==l_True && under>= c.w){
					return false;
				}else if (value(c.l)==l_False && over<c.w){
					return false;
				}
			}

			for(int cID:comparisons_gt[bvID]){
				Comparison & c = comparisons[cID];
				if(value(c.l)==l_True && over<= c.w){
					return false;
				}else if (value(c.l)==l_False && under>c.w){
					return false;
				}
			}

		}
		return true;
	}
	
	bool dbg_solved() {

		return true;
	}


	CRef newReasonMarker(int detectorID) {
		CRef reasonMarker = S->newReasonMarker(this);
		int mnum = CRef_Undef - reasonMarker;
		marker_map.growTo(mnum + 1);
		marker_map[mnum] = detectorID;
		return reasonMarker;
	}
	void newBitvector(int bvID, vec<Var> & vars){
		bitvectors.growTo(bvID+1);
		under_approx.growTo(bvID+1,-1);
		over_approx.growTo(bvID+1,-1);
		comparisons_lt.growTo(bvID+1);
		comparisons_gt.growTo(bvID+1);
		alteredBV.growTo(bvID+1);
		if(under_approx[bvID]>-1){
			assert(false);
			fprintf(stderr,"Redefined bitvector, Aborting!");
			exit(1);
		}
		under_approx[bvID]=0;
		over_approx[bvID]=0;

		for(int i = 0;i<vars.size();i++){
			bitvectors[bvID].push(mkLit(newVar(vars[i],bvID)));
		}

		alteredBV[bvID]=true;
		altered_bvs.push(bvID);
	}
	void newBitvector(int bvID, int bitwidth){
		bitvectors.growTo(bvID+1);
		under_approx.growTo(bvID+1,-1);
		over_approx.growTo(bvID+1,-1);
		comparisons_lt.growTo(bvID+1);
		comparisons_gt.growTo(bvID+1);
		alteredBV.growTo(bvID+1);
		if(under_approx[bvID]>-1){
			assert(false);
			fprintf(stderr,"Redefined bitvector, Aborting!");
			exit(1);
		}
		under_approx[bvID]=0;
		over_approx[bvID]=0;

		for(int i = 0;i<bitwidth;i++){
			bitvectors[bvID].push(mkLit(newVar(var_Undef,bvID)));
		}

		alteredBV[bvID]=true;
		altered_bvs.push(bvID);
	}


	bool hasBV(int bvID){
		return bvID>=0 && bvID<under_approx.size() && under_approx[bvID]>-1;
	}

	Lit getComparisonLT(int bvID, Weight & lt){
		//could do a binary search here:
		for(int i=0;i<comparisons_lt[bvID].size()-1;i++){
			int cID = comparisons_lt[bvID][i];
			if (comparisons[cID].w == lt){
				return comparisons[cID].l;
			}
		}

		return lit_Undef;
	}
	Lit getComparisonGT(int bvID, Weight & gt){
		//could do a binary search here:
		for(int i=0;i<comparisons_gt[bvID].size()-1;i++){
			int cID = comparisons_gt[bvID][i];
			if (comparisons[cID].w == gt){
				return comparisons[cID].l;
			}
		}

		return lit_Undef;
	}

	Lit newComparisonLT(int bvID, Weight & lt, Var outerVar = var_Undef) {
		Lit l;
		if(!hasBV(bvID)){
			exit(1);
		}
		int comparisonID = comparisons.size();
		if((l = getComparisonLT(bvID, lt))!=lit_Undef){
			if(outerVar != var_Undef){
				makeEqualInSolver(mkLit(outerVar),toSolver(l));
			}
			return l;
		}else{
			l = mkLit(newVar(outerVar, bvID,comparisonID));
		}


		comparisons.push({lt,l,bvID,true});
		comparisons_lt[bvID].push(comparisonID);
		//insert this value in order.
		//could do a binary search here...
		for(int i=0;i<comparisons_lt[bvID].size()-1;i++){
			int cid = comparisons_lt[bvID][i];
			if(comparisons[cid].w>= lt){
				for(int j = comparisons_lt[bvID].size()-1; j>i ;j--){
					comparisons_lt[bvID][j]=comparisons_lt[bvID][j-1];
				}
				comparisons_lt[bvID][i]=comparisonID;
				break;
			}
		}
		if(!alteredBV[bvID]){
			alteredBV[bvID]=true;
			altered_bvs.push(bvID);
		}
		return l;
	}
	Lit newComparisonGT(int bvID, Weight & gt, Var outerVar = var_Undef) {
			Lit l;
			int comparisonID = comparisons.size();
			if((l = getComparisonGT(bvID, gt))!=lit_Undef){
				if(outerVar != var_Undef){
					makeEqualInSolver(mkLit(outerVar),toSolver(l));
				}
				return l;
			}else{
				l = mkLit(newVar(outerVar, bvID,comparisonID));
			}



			comparisons.push({gt,l,bvID,false});
			comparisons_gt[bvID].push(comparisonID);
			//insert this value in order.
			//could do a binary search here...
			for(int i=0;i<comparisons_gt[bvID].size()-1;i++){
				int cid = comparisons_gt[bvID][i];
				if(comparisons[cid].w>= gt){
					for(int j = comparisons_gt[bvID].size()-1; j>i ;j--){
						comparisons_gt[bvID][j]=comparisons_gt[bvID][j-1];
					}
					comparisons_gt[bvID][i]=comparisonID;
					break;
				}
			}
			if(!alteredBV[bvID]){
				alteredBV[bvID]=true;
				altered_bvs.push(bvID);
			}
			return l;
		}
	void printSolution() {

	}

};

}
;

#endif
