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
#include "FSMGeneratorAcceptorDetector.h"

#include "FSMTheory.h"
#include "core/Config.h"

using namespace Monosat;

FSMGeneratorAcceptorDetector::FSMGeneratorAcceptorDetector(int detectorID, FSMTheorySolver * outer, DynamicFSM &g_under,
		DynamicFSM &g_over,DynamicFSM & acceptor_under,DynamicFSM & acceptor_over,int gen_source, int accept_source,  double seed) :
		FSMDetector(detectorID), outer(outer), g_under(g_under), g_over(g_over), acceptor_under(acceptor_under), acceptor_over(acceptor_over),gen_source(gen_source), accept_source(accept_source), rnd_seed(seed){

	underReachStatus = new FSMGeneratorAcceptorDetector::AcceptStatus(*this, true);
	overReachStatus = new FSMGeneratorAcceptorDetector::AcceptStatus(*this, false);

	underapprox_detector = new NFALinearGeneratorAcceptor<FSMGeneratorAcceptorDetector::AcceptStatus>(g_under,acceptor_under,gen_source,accept_source,*underReachStatus);
	overapprox_detector = new NFALinearGeneratorAcceptor<FSMGeneratorAcceptorDetector::AcceptStatus>(g_over,acceptor_over,gen_source,accept_source,*overReachStatus);

	underprop_marker = outer->newReasonMarker(getID());
	overprop_marker = outer->newReasonMarker(getID());
}

void FSMGeneratorAcceptorDetector::addAcceptLit(int generatorFinalState, int acceptorFinalState, Var outer_reach_var){



	if(first_destination==-1)
		first_destination= generatorFinalState;

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
	all_lits.push(acceptLit);


	while (accept_lit_map.size() <= accept_var - first_var) {
		accept_lit_map.push({-1,-1});
	}
	accept_lit_map[accept_var - first_var] = {acceptLit,generatorFinalState,acceptorFinalState};
	all_accept_lits.push( {acceptLit,generatorFinalState,acceptorFinalState});
}

void FSMGeneratorAcceptorDetector::AcceptStatus::accepts(int string,int state,int edgeID,int label, bool accepts){
/*	Lit l = detector.accept_lits[string][state];

	if (l != lit_Undef){// && !detector.is_changed[detector.indexOf(var(l))]) {
		if(!accepts){
			l=~l;
		}
		if (polarity == accepts){
			lbool assign = detector.outer->value(l);
			//detector.is_changed[detector.indexOf(var(l))] = true;
			detector.changed.push( { l, state,string});
		}
	}*/
}



bool FSMGeneratorAcceptorDetector::propagate(vec<Lit> & conflict) {
	static int iter = 0;
	if (++iter == 87) {
		int a = 1;
	}
	for(auto & t:all_accept_lits){

			Lit l =t.l;

			int gen_to = t.gen_to;
			int accept_to = t.accept_to;
			//assert(is_changed[indexOf(var(l))]);

			if(outer->value(l)!=l_True && underapprox_detector->accepts(gen_to,accept_to)){
				if (outer->value(l) == l_True) {
					//do nothing
				} else if (outer->value(l) == l_Undef) {
					outer->enqueue(l, underprop_marker);
				}else{
					conflict.push(l);
					buildAcceptReason(gen_to,accept_to, conflict);
					return false;
				}
			}else if (outer->value(l)!=l_False && !overapprox_detector->accepts(gen_to,accept_to)){
				l=~l;
				if (outer->value(l) == l_True) {
					//do nothing
				} else if (outer->value(l) == l_Undef) {
					outer->enqueue(l, overprop_marker);
				}else{
					conflict.push(l);
					buildNonAcceptReason(gen_to,accept_to, conflict);
					return false;
				}
			}else{

			}
		}

	return true;
}


void FSMGeneratorAcceptorDetector::buildReason(Lit p, vec<Lit> & reason, CRef marker) {
	if (marker == underprop_marker) {
		reason.push(p);
		Var v = var(p);
		int gen_final = getGeneratorFinal(v);
		int accept_final = getAcceptorFinal(v);
		buildAcceptReason(gen_final,accept_final, reason);
	} else if (marker == overprop_marker) {
		reason.push(p);
		Var v = var(p);
		int gen_final = getGeneratorFinal(v);
		int accept_final = getAcceptorFinal(v);
		buildNonAcceptReason(gen_final,accept_final, reason);
	}  else {
		assert(false);
	}
}

void FSMGeneratorAcceptorDetector::buildAcceptReason(int genFinal, int acceptFinal, vec<Lit> & conflict){
	static int iter = 0;
	++iter;
//find a path - ideally, the one that traverses the fewest unique transitions - from source to node, learn that one of the transitions on that path must be disabled.

	static vec<NFATransition> path;
	path.clear();
	underapprox_detector->getGeneratorPath(genFinal,acceptFinal,path);

	assert(underapprox_detector->accepts(genFinal,acceptFinal));
	for(auto & t:path){
		int edgeID = t.edgeID;
		int input = t.input;
		assert(input==0);
		int output = t.output;
		Var v = outer->getTransitionVar(g_over.getID(),edgeID,0,output);
		assert(outer->value(v)==l_True);
		conflict.push(mkLit(v,true));
	}
	//note: if there are repeated edges in this conflict, they will be cheaply removed by the sat solver anyhow, so that is not a major problem.

}
void FSMGeneratorAcceptorDetector::buildNonAcceptReason(int genFinal, int acceptFinal, vec<Lit> & conflict){

	static int iter = 0;

}
void FSMGeneratorAcceptorDetector::printSolution(std::ostream& out){

	for(auto & t:all_accept_lits){
		Lit l =t.l;

		int gen_to = t.gen_to;
		int accept_to = t.accept_to;

		if(outer->value(l)==l_True){
			static vec<NFATransition> path;
			path.clear();
			assert(underapprox_detector->accepts(gen_to,accept_to));
			underapprox_detector->getGeneratorPath(gen_to,accept_to,path);

			out<<"Generated string: ";
			for(auto & t:path){
				int edgeID = t.edgeID;
				int output = t.output;
				Var v = outer->getTransitionVar(g_over.getID(),edgeID,0,output);
				assert(outer->value(v)==l_True);
				out<< output;
			}
			out<<"\n";
		}else{
			static vec<NFATransition> path;
			path.clear();
			if(g_under.generates(gen_source,gen_to, path)){
				out<<"Generated string: ";
				for(auto & t:path){
					int edgeID = t.edgeID;
					int output = t.output;
					Var v = outer->getTransitionVar(g_over.getID(),edgeID,0,output);
					assert(outer->value(v)==l_True);
					out<< output;
				}
				out<<"\n";
			}

		}
	}
}
bool FSMGeneratorAcceptorDetector::checkSatisfied(){
	printSolution(std::cout);
	NFALinearGeneratorAcceptor<> check(g_under,acceptor_under,gen_source,accept_source);

	//g_under.draw(gen_source,first_destination );
	for(auto & t:all_accept_lits){
		Lit l =t.l;

		int gen_to = t.gen_to;
		int accept_to = t.accept_to;


		if(outer->value(l)==l_Undef){
			return false;
		}else if (outer->value(l)==l_False && check.accepts(gen_to,accept_to)){
			return false;
		}else if (outer->value(l)==l_True && !check.accepts(gen_to,accept_to)){
			return false;
		}


	}


	return true;
}


