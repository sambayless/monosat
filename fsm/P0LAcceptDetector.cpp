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
#include "P0LAcceptDetector.h"

#include "LSystemTheory.h"
#include "core/Config.h"

using namespace Monosat;


P0LAcceptDetector::P0LAcceptDetector(int detectorID, LSystemSolver * outer, LSystem &g_under,
		LSystem &g_over,vec<vec<int>> & str, double seed) :
		FSMDetector(detectorID), outer(outer), g_under(g_under), g_over(g_over), strings(str), rnd_seed(seed){

	underapprox_detector = new NP0LAccept(g_under,str);
	overapprox_detector = new NP0LAccept(g_over,str);

	underprop_marker = outer->newReasonMarker(getID());
	overprop_marker = outer->newReasonMarker(getID());
}

void P0LAcceptDetector::addProducesLit(int strID, Var outer_reach_var){

	accept_lits.growTo(strings.size(),lit_Undef);

	if (accept_lits[strID] != lit_Undef) {
		Lit r = accept_lits[strID];
		//force equality between the new lit and the old reach lit, in the SAT solver
		outer->makeEqualInSolver(outer->toSolver(r), mkLit(outer_reach_var));
		return;
	}

	g_under.invalidate();
	g_over.invalidate();

	Var accept_var = outer->newVar(outer_reach_var, getID());

	if (first_var == var_Undef) {
		first_var = accept_var;
	} else {
		assert(accept_var >= first_var);
	}
	int index = accept_var - first_var;

	//is_changed.growTo(index+1);
	Lit acceptLit = mkLit(accept_var, false);
	all_lits.push({acceptLit,strID});
	assert(accept_lits[strID]== lit_Undef);
	//if(reach_lits[to]==lit_Undef){
	accept_lits[strID] = acceptLit;
	while (accept_lit_map.size() <= accept_var - first_var) {
		accept_lit_map.push({-1});
	}
	accept_lit_map[accept_var - first_var] = {strID};

}




bool P0LAcceptDetector::propagate(vec<Lit> & conflict) {
	static int iter = 0;
	if (++iter == 87) {
		int a = 1;
	}
	changed.clear();
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


	for(auto & t:all_lits){

		Lit l = t.l;
		int string = t.str;


		if(underapprox_detector->acceptsString(string)){
			if (outer->value(l) == l_True) {
				//do nothing
			} else if (outer->value(l) == l_Undef) {
				outer->enqueue(l, underprop_marker);
			}else{
				conflict.push(l);
				buildAcceptReason(string, conflict);
				return false;
			}
		}else if (!overapprox_detector->acceptsString(string)){
			l=~l;
			if (outer->value(l) == l_True) {
				//do nothing
			} else if (outer->value(l) == l_Undef) {
				outer->enqueue(l, overprop_marker);
			}else{
				conflict.push(l);
				buildNonAcceptReason(string, conflict);
				return false;
			}
		}else{

		}
	}

	return true;
}


void P0LAcceptDetector::buildReason(Lit p, vec<Lit> & reason, CRef marker) {
	if (marker == underprop_marker) {
		reason.push(p);
		Var v = var(p);

		int str = getString(v);
		buildAcceptReason(str, reason);
	} else if (marker == overprop_marker) {
		reason.push(p);
		Var v = var(p);

		int str = getString(v);
		buildNonAcceptReason(str, reason);
	}  else {
		assert(false);
	}
}

void P0LAcceptDetector::buildAcceptReason(int str, vec<Lit> & conflict){

}
void P0LAcceptDetector::buildNonAcceptReason(int str, vec<Lit> & conflict){
	static vec<int> blocking;
	blocking.clear();
	overapprox_detector->blockingEdges(str,blocking);

	for(int ruleID:blocking){
		Var v = outer->getRuleVar(ruleID);
		if (v!=var_Undef){
			assert(outer->value(v)==l_False);
			if (outer->level(v)>0){
				//learn v
				conflict.push(mkLit(v));//rely on the sat solver to remove duplicates, here...
			}
		}
	}

}
void P0LAcceptDetector::printSolution(std::ostream& out){

}
bool P0LAcceptDetector::checkSatisfied(){
	NP0LAccept check(g_under,strings);


	for(int str = 0;str<accept_lits.size();str++){
		vec<int> & string = strings[str];


			if(accept_lits[str]!=lit_Undef){
				Lit l = accept_lits[str];
				if(outer->value(l)==l_Undef){
					return false;
				}

				else if (outer->value(l)==l_False && check.acceptsString(str)){
					return false;
				}else if (outer->value(l)==l_True && !check.acceptsString(str)){
					return false;
				}

			}


	}

	return true;
}


