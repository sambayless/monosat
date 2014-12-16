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
#include "mtl/Vec.h"
#include "FSMReachDetector.h"

#include "FSMTheory.h"
#include "core/Config.h"

using namespace Monosat;


FSMReachDetector::FSMReachDetector(int detectorID, FSMTheorySolver * outer, DynamicFSM &g_under,
		DynamicFSM &g_over, int source, vec<vec<int>> & str, double seed) :
		FSMDetector(detectorID), outer(outer), g_under(g_under), g_over(g_over), source(source),strings(str), rnd_seed(seed){

	underReachStatus = new FSMReachDetector::ReachStatus(*this, true);
	overReachStatus = new FSMReachDetector::ReachStatus(*this, false);

	underapprox_detector = new NFAReach<FSMReachDetector::ReachStatus>(g_under,source,str,*underReachStatus);
	overapprox_detector = new NFAReach<FSMReachDetector::ReachStatus>(g_over,source,str,*overReachStatus);

	underprop_marker = outer->newReasonMarker(getID());
	overprop_marker = outer->newReasonMarker(getID());
}

void FSMReachDetector::addReachLit(int state, int strID, Var outer_reach_var){

	reach_lits.growTo(strings.size());
	reach_lits[strID].growTo(g_under.nodes(),lit_Undef);



	if (reach_lits[strID][state] != lit_Undef) {
		Lit r = reach_lits[strID][state];
		//force equality between the new lit and the old reach lit, in the SAT solver
		outer->makeEqualInSolver(outer->toSolver(r), mkLit(outer_reach_var));
		return;
	}

	g_under.invalidate();
	g_over.invalidate();

	Var reach_var = outer->newVar(outer_reach_var, getID());

	if (first_reach_var == var_Undef) {
		first_reach_var = reach_var;
	} else {
		assert(reach_var >= first_reach_var);
	}
	is_changed.growTo(indexOf(reach_var));
	Lit reachLit = mkLit(reach_var, false);
	assert(reach_lits[strID][state] == lit_Undef);
	//if(reach_lits[to]==lit_Undef){
	reach_lits[strID][state] = reachLit;
	while (reach_lit_map.size() <= reach_var - first_reach_var) {
		reach_lit_map.push({-1,-1});
	}
	reach_lit_map[reach_var - first_reach_var] = {strID,state};

}

bool FSMReachDetector::propagate(vec<Lit> & conflict) {
	static int iter = 0;
	if (++iter == 87) {
		int a = 1;
	}

	bool skipped_positive = false;
	if (underapprox_detector && (!opt_detect_pure_theory_lits || unassigned_positives > 0)) {
		double startdreachtime = rtime(2);
		stats_under_updates++;
		underapprox_detector->update();
		double reachUpdateElapsed = rtime(2) - startdreachtime;
		//outer->reachupdatetime+=reachUpdateElapsed;
		stats_under_update_time += rtime(2) - startdreachtime;
	} else {
		skipped_positive = true;
		//outer->stats_pure_skipped++;
		stats_skipped_under_updates++;
	}
	bool skipped_negative = false;
	if (overapprox_detector && (!opt_detect_pure_theory_lits || unassigned_negatives > 0)) {
		double startunreachtime = rtime(2);
		stats_over_updates++;
		overapprox_detector->update();
		double unreachUpdateElapsed = rtime(2) - startunreachtime;
		stats_over_update_time += rtime(2) - startunreachtime;
	} else {
		skipped_negative = true;
		stats_skipped_over_updates++;
	}

	if (opt_rnd_shuffle) {
		randomShuffle(rnd_seed, changed);
	}

	while (changed.size()) {
		int sz = changed.size();
		Var v = changed.last().v;
		int u = changed.last().u;
		int str = changed.last().str;
		assert(is_changed[indexOf(v)]);
		Lit l;

		if (underapprox_detector && !skipped_positive && underapprox_detector->reachable(u)) {
			l = mkLit(v, false);
		} else if (overapprox_detector && !skipped_negative && !overapprox_detector->reachable(u)) {
			l = mkLit(v, true);
		} else {
			assert(sz == changed.size());
			assert(changed.last().u == u);
			is_changed[indexOf(v)] = false;
			changed.pop();
			//this can happen if the changed node's reachability status was reported before a backtrack in the solver.
			continue;
		}

		bool reach = !sign(l);
		if (outer->value(l) == l_True) {
			//do nothing
		} else if (outer->value(l) == l_Undef) {

			//trail.push(Assignment(false,reach,detectorID,0,var(l)));
			if (reach)
				outer->enqueue(l, underprop_marker);
			else
				outer->enqueue(l, overprop_marker);
		} else if (outer->value(l) == l_False) {
			conflict.push(l);

			if (reach) {

				//conflict
				//The reason is a path in g from to s in d
				buildReachReason(u,str, conflict);
				//add it to s
				//return it as a conflict

			} else {
				//The reason is a cut separating s from t
				buildNonReachReason(u,str, conflict);

			}
			return false;
		}

		assert(sz == changed.size());//This can be really tricky - if you are not careful, an a reach detector's update phase was skipped at the beginning of propagate, then if the reach detector is called during propagate it can push a change onto the list, which can cause the wrong item to be removed here.
		assert(changed.last().u == u);
		is_changed[indexOf(v)] = false;
		changed.pop();
	}

	return true;
}

void FSMReachDetector::ReachStatus::reaches(int string,int state,int edgeID,int label){
	Lit l = detector.reach_lits[string][state];
	if (l != lit_Undef && !detector.is_changed[state]) {
		lbool assign = detector.outer->value(l);
		detector.is_changed[detector.indexOf(var(l))] = true;
		detector.changed.push( { var(l), state,string});
	}
}

void FSMReachDetector::buildReason(Lit p, vec<Lit> & reason, CRef marker) {
	if (marker == underprop_marker) {
		reason.push(p);
		Var v = var(p);
		int u = getState(v);
		int str = getString(v);
		buildReachReason(u,str, reason);
	} else if (marker == overprop_marker) {
		reason.push(p);
		Var v = var(p);
		int t = getState(v);
		int str = getString(v);
		buildNonReachReason(t,str, reason);
	}  else {
		assert(false);
	}
}

void FSMReachDetector::buildReachReason(int node,int str, vec<Lit> & conflict){

}
void FSMReachDetector::buildNonReachReason(int node,int str, vec<Lit> & conflict){

}
void FSMReachDetector::printSolution(std::ostream& out){

}
bool FSMReachDetector::checkSatisfied(){
	NFAReach<> check(g_under,source,strings);
	check.update();
	for(int str = 0;str<reach_lits.size();str++){
		vec<int> & string = strings[str];

		for(int to = 0;to<reach_lits[str].size();to++){
			if(reach_lits[str][to]!=lit_Undef){
				Lit l = reach_lits[str][to];
				if(outer->value(l)==l_Undef){
					return false;
				}

				else if (outer->value(l)==l_False && check.reachable(to)){
					return false;
				}else if (outer->value(l)==l_True && !check.reachable(to)){
					return false;
				}

			}
		}

	}

	return true;
}


