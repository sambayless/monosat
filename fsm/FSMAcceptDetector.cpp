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
#include "FSMAcceptDetector.h"

#include "FSMTheory.h"
#include "core/Config.h"

using namespace Monosat;


FSMAcceptDetector::FSMAcceptDetector(int detectorID, FSMTheorySolver * outer, DynamicFSM &g_under,
		DynamicFSM &g_over, int source, vec<vec<int>> & str, double seed) :
		FSMDetector(detectorID), outer(outer), g_under(g_under), g_over(g_over), source(source),strings(str), rnd_seed(seed){

	underReachStatus = new FSMAcceptDetector::AcceptStatus(*this, true);
	overReachStatus = new FSMAcceptDetector::AcceptStatus(*this, false);

	underapprox_detector = new NFAAccept<FSMAcceptDetector::AcceptStatus>(g_under,source,str,*underReachStatus);
	overapprox_detector = new NFAAccept<FSMAcceptDetector::AcceptStatus>(g_over,source,str,*overReachStatus);

	underprop_marker = outer->newReasonMarker(getID());
	overprop_marker = outer->newReasonMarker(getID());
}

void FSMAcceptDetector::addAcceptLit(int state, int strID, Var outer_reach_var){

	accept_lits.growTo(strings.size());
	accept_lits[strID].growTo(g_under.nodes(),lit_Undef);
	if(first_destination==-1)
		first_destination= state;

	if (accept_lits[strID][state] != lit_Undef) {
		Lit r = accept_lits[strID][state];
		//force equality between the new lit and the old reach lit, in the SAT solver
		outer->makeEqualInSolver(outer->toSolver(r), mkLit(outer_reach_var));
		return;
	}

	g_under.invalidate();
	g_over.invalidate();

	Var accept_var = outer->newVar(g_over.getID(),outer_reach_var, getID());

	if (first_var == var_Undef) {
		first_var = accept_var;
	} else {
		assert(accept_var >= first_var);
	}
	int index = accept_var - first_var;

	//is_changed.growTo(index+1);
	Lit acceptLit = mkLit(accept_var, false);
	all_lits.push(acceptLit);
	assert(accept_lits[strID][state] == lit_Undef);
	//if(reach_lits[to]==lit_Undef){
	accept_lits[strID][state] = acceptLit;
	while (accept_lit_map.size() <= accept_var - first_var) {
		accept_lit_map.push({-1,-1});
	}
	accept_lit_map[accept_var - first_var] = {strID,state};

}

void FSMAcceptDetector::AcceptStatus::accepts(int string,int state,int edgeID,int label, bool accepts){
	Lit l = detector.accept_lits[string][state];

	if (l != lit_Undef){// && !detector.is_changed[detector.indexOf(var(l))]) {
		if(!accepts){
			l=~l;
		}
		if (polarity == accepts){
			lbool assign = detector.outer->value(l);
			//detector.is_changed[detector.indexOf(var(l))] = true;
			detector.changed.push( { l, state,string});
		}
	}
}



bool FSMAcceptDetector::propagate(vec<Lit> & conflict) {
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


	while (changed.size()) {
			int sz = changed.size();
			Lit l = changed.last().l;

			int u = changed.last().u;
			int str = changed.last().str;
			//assert(is_changed[indexOf(var(l))]);
			if(sign(l)){
				assert(!overapprox_detector->acceptsString(str,u));
			}else{
				assert(underapprox_detector->acceptsString(str,u));
			}

			if (underapprox_detector && !skipped_positive && !sign(l)) {

			} else if (overapprox_detector && !skipped_negative && sign(l)) {

			} else {
				assert(sz == changed.size());
				assert(changed.last().u == u);
				//is_changed[indexOf(var(l))] = false;
				changed.pop();
				//this can happen if the changed node's reachability status was reported before a backtrack in the solver.
				continue;
			}

			bool reach = !sign(l);
			if (outer->value(l) == l_True) {
				//do nothing
			} else if (outer->value(l) == l_Undef) {

				if (reach)
					outer->enqueue(l, underprop_marker);
				else
					outer->enqueue(l, overprop_marker);
			} else if (outer->value(l) == l_False) {
				conflict.push(l);

				if (reach) {
					buildAcceptReason(u,str, conflict);
				} else {
					//The reason is a cut separating s from t
					buildNonAcceptReason(u,str, conflict);
				}

				return false;
			}

			assert(sz == changed.size());//This can be really tricky - if you are not careful, an a reach detector's update phase was skipped at the beginning of propagate, then if the reach detector is called during propagate it can push a change onto the list, which can cause the wrong item to be removed here.
			assert(changed.last().u == u);
			//is_changed[indexOf(var(l))] = false;
			changed.pop();
		}
	assert(changed.size()==0);
	return true;
}


void FSMAcceptDetector::buildReason(Lit p, vec<Lit> & reason, CRef marker) {
	if (marker == underprop_marker) {
		reason.push(p);
		Var v = var(p);
		int u = getState(v);
		int str = getString(v);
		buildAcceptReason(u,str, reason);
	} else if (marker == overprop_marker) {
		reason.push(p);
		Var v = var(p);
		int t = getState(v);
		int str = getString(v);
		buildNonAcceptReason(t,str, reason);
	}  else {
		assert(false);
	}
}

void FSMAcceptDetector::buildAcceptReason(int node,int str, vec<Lit> & conflict){
	static int iter = 0;
	++iter;
//find a path - ideally, the one that traverses the fewest unique transitions - from source to node, learn that one of the transitions on that path must be disabled.
/*	g_under.draw(source);
	vec<int> & string = strings[str];
	printf("Accepts: \"");
	for(int s:string){
		printf("%d ",s);
	}
	printf("\"\n");*/
	static vec<NFATransition> path;
	path.clear();
	bool hasPath =underapprox_detector->getPath(str,node,path);
	assert(hasPath);
	assert(underapprox_detector->acceptsString(str,node));
	for(auto & t:path){
		int edgeID = t.edgeID;
		int input = t.input;
		Var v = outer->getTransitionVar(g_over.getID(),edgeID,input,0);
		assert(outer->value(v)==l_True);
		conflict.push(mkLit(v,true));
	}
	//note: if there are repeated edges in this conflict, they will be cheaply removed by the sat solver anyhow, so that is not a major problem.

}
void FSMAcceptDetector::buildNonAcceptReason(int node,int str, vec<Lit> & conflict){

	static int iter = 0;
//optionally, remove all transitions from the graph that would not be traversed by this string operating on the level 0 overapprox graph.

	//This doesn't work:
	//then, find a cut in the remaining graph through the disabled transitions, separating source from node.

	//Instead:
	//The cut has to be made through disabled transitions in the unrolled graph, with only appropriate transitions for the string (+emoves) traversible in each time frame.
	//graph must be unrolled to length of string.

	//instead of actually unrolling the graph, I am going to traverse it backwards, 'unrolling it' implicitly.
	static vec<int> to_visit;
	static vec<int> next_visit;
	vec<int> & string = strings[str];
/*
	g_over.draw(source);

	printf("%d: Doesn't accept: \"",iter);
	for(int s:string){
		printf("%d ",s);
	}
	printf("\"\n");
*/



	if(++iter==10189){
		int a=1;
	}

	assert(!overapprox_detector->acceptsString(str,node));
	//int strpos = string.size()-1;
	to_visit.clear();
	next_visit.clear();

	static vec<bool> cur_seen;
	static vec<bool> next_seen;
	cur_seen.clear();
	cur_seen.growTo(g_under.states());

	next_seen.clear();
	next_seen.growTo(g_under.states());

	int str_pos = string.size();
	if(str_pos==0){
		//special handling for empty string. Only e-move edges are considered.
		assert(node!=source);//otherwise, this would have been accepted.

		if(!g_over.emovesEnabled()){
			//no way to get to the node without consuming strings.
			return;
		}
		cur_seen[node]=true;
		to_visit.push(node);
		for(int j = 0;j<to_visit.size();j++){
			int u = to_visit[j];

			assert(str_pos>=0);

			for (int i = 0;i<g_under.nIncoming(u);i++){
				int edgeID = g_under.incoming(u,i).id;
				int from = g_under.incoming(u,i).node;

				if(g_over.emovesEnabled()){
					if (g_over.transitionEnabled(edgeID,0,0)){
						if (!cur_seen[from]){
							cur_seen[from]=true;
							to_visit.push(from);//emove transition, if enabled
						}
					}else{
						Var v = outer->getTransitionVar(g_over.getID(),edgeID,0,0);
						if (v!=var_Undef){
							assert(outer->value(v)==l_False);
							if (outer->level(v)>0){
								//learn v
								conflict.push(mkLit(v));
							}
						}
					}
				}


			}
		}

		for(int s:to_visit){
			assert(cur_seen[s]);
			cur_seen[s]=false;
		}
		to_visit.clear();

		return;
	}


	for(int s:next_visit){
		assert(next_seen[s]);
		next_seen[s]=false;
	}
	next_visit.clear();

	next_visit.push(node);
	next_seen[node]=true;


	while(next_visit.size()){
		str_pos --;
		assert(str_pos>=0);
		next_visit.swap(to_visit);
		next_seen.swap(cur_seen);

		for(int s:next_visit){
			assert(next_seen[s]);
			next_seen[s]=false;
		}
		next_visit.clear();

		int l = string[str_pos];

		for(int j = 0;j<to_visit.size();j++){
			int u = to_visit[j];

			assert(str_pos>=0);


			for (int i = 0;i<g_under.nIncoming(u);i++){
				int edgeID = g_under.incoming(u,i).id;
				int from = g_under.incoming(u,i).node;

				if(g_over.emovesEnabled()){
					if (g_over.transitionEnabled(edgeID,0,0)){
						if (!cur_seen[from]){
							cur_seen[from]=true;
							to_visit.push(from);//emove transition, if enabled
						}
					}else{
						Var v = outer->getTransitionVar(g_over.getID(),edgeID,0,0);
						if (v!=var_Undef){
							assert(outer->value(v)==l_False);
							if (outer->level(v)>0){
								//learn v
								conflict.push(mkLit(v));
							}
						}
					}
				}

				if (g_over.transitionEnabled(edgeID,l,0)){
					if (!next_seen[from] && str_pos>0){
						next_seen[from]=true;
						next_visit.push(from);
					}
				}else{
					Var v = outer->getTransitionVar(g_over.getID(),edgeID,l,0);
					if (v!=var_Undef){
						assert(outer->value(v)==l_False);
						if (outer->level(v)>0){
							//learn v
							conflict.push(mkLit(v));//rely on the sat solver to remove duplicates, here...
						}
					}
				}
			}
		}

	}
/*	printf("conflict: ");
	for (int i = 1;i<conflict.size();i++){
		Lit l = conflict[i];
		int edgeID = outer->getEdgeID(var(l));
		int f = g_under.getEdge(edgeID).from;
		int t = g_under.getEdge(edgeID).to;
		printf("(%d->%d %d),",f,t,outer->getLabel(var(l)));
	}*/
}
void FSMAcceptDetector::printSolution(std::ostream& out){
	g_under.draw(source);
}
bool FSMAcceptDetector::checkSatisfied(){
	NFAAccept<> check(g_under,source,strings);

	g_under.draw(source,first_destination );
	for(int str = 0;str<accept_lits.size();str++){
		vec<int> & string = strings[str];
		check.run(str);
		for(int to = 0;to<accept_lits[str].size();to++){
			if(accept_lits[str][to]!=lit_Undef){
				Lit l = accept_lits[str][to];
				if(outer->value(l)==l_Undef){
					return false;
				}

				else if (outer->value(l)==l_False && check.accepting(to)){
					return false;
				}else if (outer->value(l)==l_True && !check.accepting(to)){
					return false;
				}

			}
		}

	}

	return true;
}


