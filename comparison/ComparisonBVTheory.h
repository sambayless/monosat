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
		Weight lt;
		Lit l;

		int valueID;
	};

	double rnd_seed;

private:
	Solver * S;
	int local_q = 0;
public:
	int id;

	vec<lbool> assigns;
	CRef underprop_marker;
	CRef overprop_marker;

	vec<Assignment> trail;
	vec<int> trail_lim;

	//Bitvectors are unsigned, and have least significant bit at index 0
	vec<vec<Lit> > bitvectors;
	vec<Comparison> comparisons;
	vec<vec<int> > comparisonsBV; //for each bitvector, comparisons are all to unique values, and in ascending order of compareTo.

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

		int valueID:30;
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
		 underprop_marker = S->newReasonMarker(this);
		 overprop_marker = S->newReasonMarker(this);

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
		return vars[v].isEdge;
	}
	inline int getComparisonID(Var v) const {
		assert(isComparisonVar(v));
		return vars[v].detector_edge;
	}
	inline Comparison& getComparison(int comparisonID)const{
		return comparisons[comparisonID];
	}
	inline int getValueID(Var v) const{

		return vars[v].valueID;
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
	

	Var newVar(Var solverVar=var_Undef,  int comparisonID=-1,int valueID=-1, bool connectToTheory = true) {
		if(solverVar==var_Undef){
			solverVar = S->newVar();
		}
		while (S->nVars() <= solverVar)
			S->newVar();
		Var v = vars.size();

		vars.push();
		vars[v].occursPositive=false;
		vars[v].occursNegative=false;
		vars[v].valueID=valueID;
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
					int valueID = getvalueID(e.var);
					int edgeID = getEdgeID(e.var); //e.var-min_edge_var;
					int input = getInput(e.var);
					int output = getOutput(e.var);
					assert(edgeID==e.edgeID);

					if (e.assign) {
						//g_unders[valueID]->disableTransition(edgeID, input,output);

					} else {
						//g_overs[valueID]->enableTransition(edgeID,input,output);

					}
				} else {
					//This is a reachability literal				  
					//detectors[getDetector(e.var)]->unassign(mkLit(e.var, !e.assign));
				}
				assigns[e.var] = l_Undef;
				changed = true;
			}
			trail.shrink(trail.size() - stop);
			trail_lim.shrink(trail_lim.size() - level);
			assert(trail_lim.size() == level);
			
			if (changed) {
				requiresPropagation = true;

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
				int valueID = getvalueID(e.var);
				int edgeID = getEdgeID(e.var); //e.var-min_edge_var;
				int input = getInput(e.var);
				int output = getOutput(e.var);
				assert(assigns[e.var]!=l_Undef);
				assigns[e.var] = l_Undef;
				if (e.assign) {
				//	g_unders[valueID]->disableTransition(edgeID,input,output);
				} else {
				//	g_overs[valueID]->enableTransition(edgeID,input,output);
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
		
		//assert(d < detectors.size());
		//detectors[d]->buildReason(p, reason, marker);
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
			int valueID = getValueID(v);
			if(!alteredBV[valueID]){
				alteredBV[valueID]=true;
				altered_bvs.push(valueID);
			}

		} else {

			int valueID = getValueID(var(l));
			int edgeID = getComparisonID(var(l)); //v-min_edge_var;


			//trail.push( { true, !sign(l),edgeID, v });

		}
	}
	;

	void updateApproximations(int valueID){
		vec<Lit> & bv = bitvectors[valueID];
		under_approx[valueID]=0;
		over_approx[valueID]=0;
		for(int i = 0;i<bv.size();i++){
			lbool val = value(bv[i]);
			if(val==l_True){
				Weight bit = 1<<i;
				under_approx[valueID]+=bit;
				over_approx[valueID]+=bit;
			}else if (val==l_False){

			}else{
				Weight bit = 1<<i;
				over_approx[valueID]+=bit;
			}
		}
	}

	bool checkApproxUpToDate(int valueID){
#ifndef NDEBUG
		vec<Lit> & bv = bitvectors[valueID];
		Weight under =0;
		Weight over=0;
		for(int i = 0;i<bv.size();i++){
			lbool val = value(bv[i]);
			if(val==l_True){
				Weight bit = 1<<i;
				under+=bit;
				under+=bit;
			}else if (val==l_False){

			}else{
				Weight bit = 1<<i;
				under+=bit;
			}
		}
		assert(under==under_approx[valueID]);
		assert(over==over_approx[valueID]);
#endif
		return true;
	}

	bool propagateTheory(vec<Lit> & conflict) {

		stats_propagations++;

		if (!requiresPropagation) {
			stats_propagations_skipped++;

			return true;
		}
		
		bool any_change = false;
		double startproptime = rtime(1);
		//static vec<int> detectors_to_check;
		
		conflict.clear();

		
		while(altered_bvs.size()){
			int valueID = altered_bvs.last();
			assert(alteredBV[valueID]);

			updateApproximations(valueID);

			Weight & underApprox = under_approx[valueID];
			Weight & overApprox = over_approx[valueID];

			vec<int> & compares = comparisonsBV[valueID];
			//update over approx lits
			for(int i = 0;i<compares.size();i++){
				int comparisonID = compares[i];
				Weight & lt = getComparison(comparisonID).lt;
				Lit l =  getComparison(comparisonID).l;
				if (overApprox<lt){
					if(value(l)==l_True){
						//do nothing
					}else if (value(l)==l_False){
						assert(value(l)==l_False);
						conflict.push(l);
						buildValueGEQReason(valueID,comparisonID,conflict);
					}else {
						assert(value(l)==l_Undef);
						enqueue(l, overprop_marker);
					}

				}else{
					break;
				}
			}

			//update under approx lits
			for(int i = compares.size()-1;i>=0;i--){
				int comparisonID = compares[i];
				Weight & lt = getComparison(comparisonID).lt;
				Lit l =  getComparison(comparisonID).l;
				if (underApprox>=lt){
					if(value(l)==l_True){
						assert(value(l)==l_False);
						conflict.push(~l);
						buildValueLTReason(valueID,comparisonID,conflict);
					}else if (value(l)==l_False){
						//do nothing
					}else {
						assert(value(l)==l_Undef);
						enqueue(~l, underprop_marker);
					}

				}else{
					break;
				}
			}
			alteredBV[valueID]=false;
			altered_bvs.pop();
		}
		


		
		requiresPropagation = false;

		double elapsed = rtime(1) - startproptime;
		propagationtime += elapsed;

		return true;
	};

	void buildValueLTReason(int valueID, int comparisonID, vec<Lit> & conflict){

		Weight & lt = getComparison(comparisonID).lt;
		Lit l =  getComparison(comparisonID).l;

		vec<Lit> & bv = bitvectors[valueID];
		Weight  over = over_approx[valueID];
		assert(checkApproxUpToDate(valueID));

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

	void buildValueGEQReason(int valueID, int comparisonID, vec<Lit> & conflict){
		Weight & lt = getComparison(comparisonID).lt;
		Lit l =  getComparison(comparisonID).l;

		vec<Lit> & bv = bitvectors[valueID];
		Weight  under = under_approx[valueID];
		assert(checkApproxUpToDate(valueID));

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

	bool solveTheory(vec<Lit> & conflict) {
		requiresPropagation = true;		//Just to be on the safe side... but this shouldn't really be required.
		bool ret = propagateTheory(conflict);
		//Under normal conditions, this should _always_ hold (as propagateTheory should have been called and checked by the parent solver before getting to this point).
		assert(ret);
		return ret;
	}

	bool check_solved() {

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
	void newValue(int valueID, int bitwidth){
		bitvectors.growTo(valueID+1);
		under_approx.growTo(valueID+1,-1);
		over_approx.growTo(valueID+1,-1);
		comparisonsBV.growTo(valueID+1);

		if(under_approx[valueID]>-1){
			assert(false);
			fprintf(stderr,"Redefined bitvector, Aborting!");
			exit(1);
		}
		under_approx[valueID]=0;
		over_approx[valueID]=0;

		for(int i = 0;i<bitwidth;i++){
			bitvectors[i].push(mkLit(newVar(var_Undef,valueID)));
		}


	}

	Lit getComparisonLit(int valueID, Weight & lt){
		//could do a binary search here:
		for(int i=0;i<comparisonsBV[valueID].size()-1;i++){
			int cID = comparisonsBV[valueID][i];
			if (comparisons[cID].lt == lt){
				return comparisons[cID].l;
			}
		}

		return lit_Undef;
	}

	Lit newComparison(int valueID, Weight & lt, Var outerVar = var_Undef) {
		Lit l;
		if((l = getComparisonLit(valueID, lt))!=lit_Undef){
			if(outerVar != var_Undef){
				makeEqualInSolver(mkLit(outerVar),toSolver(l));
			}
			return l;
		}

		int comparisonID = comparisons.size();
		comparisons.push({lt,l,valueID});
		comparisonsBV[valueID].push(comparisonID);
		//insert this value in order.
		//could do a binary search here...
		for(int i=0;i<comparisonsBV[valueID].size()-1;i++){
			int cid = comparisonsBV[valueID][i];
			if(comparisons[cid].lt>= lt){
				for(int j = comparisonsBV[valueID].size()-1; j>i ;j--){
					comparisonsBV[j]=comparisonsBV[j-1];
				}
				comparisonsBV[i]=comparisonID;
				break;
			}
		}

		return l;
	}

	void printSolution() {

	}

};

}
;

#endif
