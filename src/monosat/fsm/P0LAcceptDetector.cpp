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

#include "LSystemTheory.h"

using namespace Monosat;


P0LAcceptDetector::P0LAcceptDetector(int detectorID, LSystemSolver * outer, LSystem &g_under,
		LSystem &g_over,vec<vec<int>> & str, double seed) :
		FSMDetector(detectorID), outer(outer), g_under(g_under), g_over(g_over), strings(str), rnd_seed(seed){

	underapprox_detector = new NP0LAccept(g_under,str);
	overapprox_detector = new NP0LAccept(g_over,str);

	underprop_marker = outer->newReasonMarker(getID());
	overprop_marker = outer->newReasonMarker(getID());
	buildAcceptors();
}

void P0LAcceptDetector::addProducesLit(int atom,int strID, Var outer_reach_var){
	while(accept_lits.size()<strings.size()){
		accept_lits.push();
		accept_lits.last().growTo(g_under.nCharacters()+1,lit_Undef);
	}


	if (accept_lits[strID][atom+1] != lit_Undef) {
		Lit r = accept_lits[strID][atom+1] ;
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
	all_lits.push({acceptLit,strID,atom});
	assert(accept_lits[strID][atom+1] == lit_Undef);
	//if(reach_lits[to]==lit_Undef){
	accept_lits[strID][atom+1]  = acceptLit;
	while (accept_lit_map.size() <= accept_var - first_var) {
		accept_lit_map.push({-1,-1});
	}
	accept_lit_map[accept_var - first_var] = {strID,atom};

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
		int atom = t.atom;

		if(outer->value(l)!=l_True && underapprox_detector->acceptsString(string)){
			if (outer->value(l) == l_True) {
				//do nothing
			} else if (outer->value(l) == l_Undef) {
				outer->enqueue(l, underprop_marker);
			}else{
				conflict.push(l);
				buildAcceptReason(atom,string, conflict);
				return false;
			}
		}else if (outer->value(l)!=l_False && !overapprox_detector->acceptsString(string)){
			l=~l;
			if (outer->value(l) == l_True) {
				//do nothing
			} else if (outer->value(l) == l_Undef) {
				outer->enqueue(l, overprop_marker);
			}else{
				conflict.push(l);
				buildNonAcceptReason(atom,string, conflict);
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
		int atom = getAtom(v);
		int str = getString(v);
		buildAcceptReason(atom,str, reason);
	} else if (marker == overprop_marker) {
		reason.push(p);
		Var v = var(p);
		int atom = getAtom(v);
		int str = getString(v);
		buildNonAcceptReason(atom,str, reason);
	}  else {
		assert(false);
	}
}

void P0LAcceptDetector::buildAcceptReason(int atom,int str, vec<Lit> & conflict){
	//the reason the string was accepted is simply the set of all unique ruleIDs that were traversed accepting the string.
	static vec<int> store_edges;
	store_edges.clear();
	bool r = underapprox_detector->getUsedRules(str,store_edges);
	assert(r);
	for(int ruleID:store_edges){
		Var v = outer->getRuleVar(ruleID);
		if (v!=var_Undef){
			assert(outer->value(v)==l_True);
			if (outer->level(v)>0){
				//learn v
				conflict.push(~mkLit(v));
			}
		}
	}

}

bool P0LAcceptDetector::find_path(int source, int dest,vec<int> & string,int depth,vec<Bitset> & suffixTable, vec<int> & path,vec<int> * blocking_edges){
	int atom =0;
	int str_pos=0;
	int emove_count=0;

	struct PathElement{
		int emove:1;//every move is either an emove or a string move.
		int used_rule:1;
		int on_path:1;
		//int edgeID:29;
		int state:29;
		int next_edge;
		int next_transition;
		PathElement(bool emove, bool used_rule, bool on_path, int state, int next_edge,int next_transition):emove(emove),used_rule(used_rule),on_path(on_path),state(state),next_edge(next_edge),next_transition(next_transition){

		}
	};

	vec<PathElement> stack;
	stack.push(PathElement(true,false,false,source,0,0));

	while(true){
		if(path.size()>=string.size()){
			//don't consider a string this long.
			//backtrack until the path is backtrack size.
			int backtrack = path.size()-1;
			assert(backtrack<path.size());
			while(backtrack<path.size()){
				PathElement & p= stack.last();
				if(p.emove){
					emove_count--;
				}else{
					str_pos--;
				}
				if (p.on_path){
					path.pop();
				}
	/*			if(p.used_rule){
					int rID = used_rule_set.last();
					used_edges[rID]=false;
					used_rule_set.pop();
				}*/
				stack.pop();
			}
			continue;
		}

		int cur_index = stack.size()-1;
		PathElement & p= stack[cur_index];

		int s = p.state;
		int next_edge = p.next_edge;
		int next_transition = p.next_transition;
		if(next_edge==0 && next_transition ==0){
			if(p.emove){
				emove_count++;
				assert(emove_count<acceptor.states());
			}else{
				str_pos++;


				if(str_pos==string.size() && (s==dest || dest<0) ){
					if(path.size()==1 && path[0]==atom){
						return true;
					}else{
						int backtrack = check_accepts(-1,depth+1,blocking_edges);
						if(backtrack== path.size()){
							//path accepted; return true.
							return true;
						}else{
							//backtrack until the path is backtrack size.
							assert(backtrack<path.size());
							while(backtrack<path.size()){
								PathElement & p2= stack.last();
								if(p2.emove){
									emove_count--;
								}else{
									str_pos--;
								}
								if (p2.on_path){
									path.pop();
								}
					/*			if(p2.used_rule){
									int rID = used_rule_set.last();
									used_edges[rID]=false;
									used_rule_set.pop();
								}*/
								stack.pop();
							}
							continue;
						}
					}
				}
			}
		}
		if(next_edge>=acceptor.nIncident(s)){

			if(p.emove){
				emove_count--;
			}else{
				str_pos--;
			}
			if (p.on_path){
				path.pop();
			}
/*			if(p.used_rule){
				int rID = used_rule_set.last();
				used_edges[rID]=false;
				used_rule_set.pop();
			}*/
			stack.pop();
			if(stack.size()==0){
				return false;
			}
		}else{

			int l = string[str_pos];
			//now check if the label is active
			int edgeID= acceptor.incident(s,next_edge).id;
			int to = acceptor.incident(s,next_edge).node;

			if(next_transition<acceptor.outAlphabet()){
				if(emove_count+1 < acceptor.states() && suffixTable[str_pos][to]){
					int o =next_transition;
						if(acceptor.transitionEnabled(edgeID,0,o)){

							if(o>0)
								path.push(o);
							bool used = false;
			/*				int rID = getRule(edgeID,0);
							if(rID>=0 && ! used_edges[rID]){
								used=true;
								used_edges[rID]=true;
								used_rule_set.push(rID);
							}*/
							stack.push(PathElement(true,used,o>0,to,0,0));
						}



				}
			}else{
				if(str_pos< string.size() && suffixTable[str_pos+1][to]){
					int o =next_transition-acceptor.outAlphabet();
					//for(int o = 0;o<acceptor.outAlphabet();o++){
					if (acceptor.transitionEnabled(edgeID,l,o)){
						if(o>0)
							path.push(o);
						bool used = false;
					/*	int rID = getRule(edgeID,l);
						if(rID>=0 && ! used_edges[rID]){
							used=true;
							used_edges[rID]=true;
							used_rule_set.push(rID);
						}*/
						stack.push(PathElement(false,used,o>0,to,0,0));
					}
					//}
				}
			}
			if(next_transition+1>=acceptor.outAlphabet()*2){
				stack[cur_index].next_transition=0;
				stack[cur_index].next_edge++;
			}else{
				stack[cur_index].next_transition++;
			}
		}
	}

}
bool P0LAcceptDetector::path_rec(int atom,int s, int dest,vec<int> & string,int str_pos,int emove_count,int depth,vec<Bitset> & suffixTable, vec<int> & path,vec<int> * blocking_edges){
	if(str_pos==string.size() && (s==dest || dest<0) ){
		//string accepted by lsystem.
		if(path.size()==0){
			return false;
		}else if(path.size()==1 && path[0]==atom){
			return true;
		}else
			return accepts_rec(atom,-1,depth+1,blocking_edges);

	}
	if (emove_count>=acceptor.states()){
		return false;//this is not a great way to solve the problem of avoiding infinite e-move cycles...
	}
	if(!suffixTable[str_pos][s])
		return false;
	int l = string[str_pos];
	for(int j = 0;j<acceptor.nIncident(s);j++){
		//now check if the label is active
		int edgeID= acceptor.incident(s,j).id;
		int to = acceptor.incident(s,j).node;
		if(suffixTable[str_pos][to]){
			for(int o = 0;o<acceptor.outAlphabet();o++){
				if(acceptor.transitionEnabled(edgeID,0,o)){
					//assert(suffixTable[str_pos][to]);
					if(o>0)
						path.push(o);

					if(path_rec(atom,to,dest,string,str_pos,emove_count+1,depth,suffixTable,path,blocking_edges)){//str_pos is NOT incremented!
						return true;
					}else if (o>0){
						path.pop();
					}
				}
			}
		}

		if(str_pos< string.size() && suffixTable[str_pos+1][to]){
			for(int o = 0;o<acceptor.outAlphabet();o++){
				if (acceptor.transitionEnabled(edgeID,l,o)){
					if(o>0)
						path.push(o);

					if(path_rec(atom,to,dest,string,str_pos+1,0,depth,suffixTable,path,blocking_edges)){//str_pos is incremented
						return true;
					}else if (o>0){
						path.pop();
					}
				}
			}
		}
	}
	return false;
}
int P0LAcceptDetector::check_accepts(int str,int depth,vec<int> * blocking_edges){
		if(depth>9)
			return 0;
		int atom = 0;
		static int iter = 0;
		++iter;

		assert(stringset.size()>depth || depth==0);
		assert(fsmstrings.size()>str);

		vec<int> & string = depth==0 ? fsmstrings[str] :stringset[depth];


		vec<Bitset> & suffixTable = suffixTables[depth];

		int backtrack=0;
		if((backtrack = acceptor.accepts_prefix(0,0,string)) < string.size()){
			if(blocking_edges){
				analyzeNFT(atom,0,0,string,*blocking_edges,suffixTable);
			}
			return backtrack;
		}


		//build suffix table of states that can reach the final state from the nth suffix of the string
		acceptor.buildSuffixTable(0,0,string,suffixTable);//<string.size());//{
		//	return suffix;
		//}



		//now find all paths through the nfa using a dfs, but filtering the search using the suffix table so that all paths explored are valid paths.

		int str_pos = 0;
		stringset[depth+1].clear();
		if( find_path(0,0,string,depth,suffixTable,stringset[depth+1],blocking_edges)){
			return string.size();
		}else{
			if(blocking_edges){
				analyzeNFT(atom,0,0,string,*blocking_edges,suffixTable);
			}
			return string.size()-1;
		}

	}

bool P0LAcceptDetector::accepts_rec(int atom,int str,int depth,vec<int> * blocking_edges){
	if(depth>9)
		return false;
	assert(stringset.size()>depth || depth==0);
	assert(strings.size()>str);
	assert(depth>=0);
	vec<int> & string = stringset[depth];

	vec<Bitset> & suffixTable = suffixTables[depth];

/*	if(!acceptor.accepts(0,0,string)){
		if(blocking_edges){
			analyzeNFA(0,0,string,*blocking_edges);
		}
		return false;
	}*/

	//build suffix table of states that can reach the final state from the nth suffix of the string
	acceptor.buildSuffixTable(0,0,string,suffixTable);



	//now find all paths through the nfa using a dfs, but filtering the search using the suffix table so that all paths explored are valid paths.


	int str_pos = 0;
	stringset[depth+1].clear();
	if(! path_rec(atom,0,0,string,0,0,depth,suffixTable,stringset[depth+1],blocking_edges)){
		if(blocking_edges){
			analyzeNFT(atom,0,0,string,*blocking_edges,suffixTable);
		}
		return false;
	}
	return true;

}

void P0LAcceptDetector::analyzeNFT(int atom,int source, int final,vec<int> & string,vec<int> & blocking,vec<Bitset> & suffixTable){
		static vec<int> to_visit;
		static vec<int> next_visit;

		//explain why this transducer only produces the set of strings it produces, and not others, given this string as input.
		//build FULL, level-0 overapprox suffix table to filter the dfs exploration below.
		bool r = acceptor_over.buildSuffixTable(source,final,string,suffixTable);
		if(!r){
			return;//no blocking edges to be found
		}
		//int strpos = string.size()-1;
		to_visit.clear();
		next_visit.clear();

		static vec<bool> cur_seen;
		static vec<bool> next_seen;
		cur_seen.clear();
		cur_seen.growTo(acceptor.states());

		next_seen.clear();
		next_seen.growTo(acceptor.states());
		int node = source;
		int str_pos = 0;
		if(string.size()==0){
			//special handling for empty string. Only e-move edges are considered.
			assert(node!=final);//otherwise, this would have been accepted.

			if(!acceptor.emovesEnabled()){
				//no way to get to the node without consuming strings.
				return;
			}
			cur_seen[node]=true;
			to_visit.push(node);
			for(int j = 0;j<to_visit.size();j++){
				int u = to_visit[j];

				assert(str_pos>=0);

				for (int i = 0;i<acceptor.nIncident(u);i++){
					int edgeID = acceptor.incident(u,i).id;
					int from = acceptor.incident(u,i).node;

					if(acceptor.emovesEnabled()){
						if (acceptor.transitionEnabled(edgeID,0,-1)){
							if (!cur_seen[from]){
								cur_seen[from]=true;
								to_visit.push(from);//emove transition, if enabled
							}
						}else{
							int ruleID = getRule(edgeID,0);
							if(ruleID>=0 && !edge_blocking[ruleID]){
								blocking.push(ruleID);
								edge_blocking[ruleID]=true;
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

			assert(str_pos>=0);
			next_visit.swap(to_visit);
			next_seen.swap(cur_seen);

			for(int s:next_visit){
				assert(next_seen[s]);
				next_seen[s]=false;
			}
			next_visit.clear();

			int l = string[str_pos];//add 1 to each character, for the nfa

			for(int j = 0;j<to_visit.size();j++){
				int u = to_visit[j];

				assert(str_pos>=0);


				for (int i = 0;i<acceptor.nIncident(u);i++){
					int edgeID = acceptor.incident(u,i).id;
					int to = acceptor.incident(u,i).node;

					if(acceptor.emovesEnabled()){
						if(suffixTable[str_pos][to]){
							if (acceptor.transitionEnabled(edgeID,0,-1)){
								if (!cur_seen[to]){
									cur_seen[to]=true;
									to_visit.push(to);//emove transition, if enabled
								}
							}else{
								int ruleID = getRule(edgeID,0);
								if(ruleID>=0 &&!edge_blocking[ruleID]){
									blocking.push(ruleID);
									edge_blocking[ruleID]=true;
								}
							}
						}
					}
					if(str_pos+1 < string.size()){
						if(suffixTable[str_pos+1][to]){
							if (acceptor.transitionEnabled(edgeID,l,-1)){
								if (!next_seen[to]){
									next_seen[to]=true;
									next_visit.push(to);
								}
							}else{
								int ruleID = getRule(edgeID,l);
								if(ruleID>=0 && !edge_blocking[ruleID]){
									blocking.push(ruleID);
									edge_blocking[ruleID]=true;
								}
							}
						}
					}
				}
			}
			str_pos ++;
		}
	}
void P0LAcceptDetector::buildAcceptors(){
	assert(outer->decisionLevel()==0);

	for(vec<int> & str:strings){
		fsmstrings.push();
		for(int c:str){
			fsmstrings.last().push(c+1);//need to add one here...
		}
	}

		LSystem & f = g_over;
		assert(f.strictlyProducing);
		//bool buildSuffixTable(int startState, int finalState, vec<int> & string, vec<Bitset> & table){
		for(int i = 0;i<=f.nCharacters();i++){
			acceptor.addInCharacter();
			acceptor.addOutCharacter();

		}
		int start = acceptor.addState();
		for(int c = 0;c<f.nCharacters();c++){
			int  cState = acceptor.addState();
			//Add an arc that takes epsilon and outputs this character.
			acceptor.addTransition(start,cState,-1,0,c+1);
			//In the future: combine rules into a prefix tree first.

			for(int rID : f.getRules(c)){
				vec<int> & rule = f.getRule(rID);
				int from = cState;
				int firstRule = -1;
				for(int i = 0;i<rule.size();i++){
					int o = rule[i];
					int next;
					if(i<rule.size()-1){
						next = acceptor.addState();
					}else{
						next=start;
					}

					//Add an arc that takes this rule character, and outputs epsilon
					int edgeID = acceptor.addTransition(from,next,-1,o+1,0);
					if(firstRule<0){
						firstRule=edgeID;
						ruleMap.growTo(rID+1);
						ruleMap[rID].edgeID = edgeID;
						ruleMap[rID].inChar=o+1;
						ruleMap[rID].outChar=0;

						rules.growTo(edgeID+1);
						rules[edgeID].growTo(acceptor.inAlphabet(),-1);//this vector size assumes that we aren't interested in output labels!
						rules[edgeID][o+1]=rID;
					}
					from = next;
				}
			}
		}
	acceptor.copyTo(acceptor_over);
	updateAcceptor();
}
void P0LAcceptDetector::updateAcceptor(){
	LSystem & g = g_over;
	if (last_modification > 0 && g.modifications == last_modification) {
			stats_skipped_updates++;
			return;
		}
		static int iteration = 0;
		int local_it = ++iteration;
		stats_full_updates++;

		if (last_deletion == g.deletions) {
			stats_num_skipable_deletions++;
		}
		bool level0 = outer->decisionLevel()==0;
		if (last_modification <= 0 || g.changed() || last_history_clear != g.historyclears) {
			for(int i = 0;i<g.nRules();i++){
				if(g.hasRule(i)){
					setRuleEnabled(i,g.ruleEnabled(i),level0);
				}
			}
		}else{
			for (int i = history_qhead; i < g.history.size(); i++) {
				int edgeid = g.history[i].id;

				if (g.history[i].addition && g.ruleEnabled(edgeid) && !ruleEnabled(edgeid)) {

					setRuleEnabled(edgeid,true,level0);

				} else if (!g.history[i].addition && !g.ruleEnabled(edgeid) && ruleEnabled(edgeid)) {

					setRuleEnabled(edgeid,false,level0);

				}
			}
		}



		last_modification = g.modifications;
		last_deletion = g.deletions;
		last_addition = g.additions;

		history_qhead = g.history.size();
		last_history_clear = g.historyclears;

}

void P0LAcceptDetector::buildNonAcceptReason(int atom,int str, vec<Lit> & conflict){
	static vec<int> store_edges;
	store_edges.clear();

	updateAcceptor();
	stringset.growTo(strings[str].size()+2);
	suffixTables.growTo(strings[str].size()+2);

	stringset[0].clear();
	for(int c:strings[str]){
		stringset[0].push(c+1);
	}

	store_edges.clear();
	edge_blocking.clear();
	edge_blocking.growTo(g_over.nRules());
	bool r = check_accepts(str,0,&store_edges)<strings[str].size();
	assert(!r);

	for(int ruleID:store_edges){
		Var v = outer->getRuleVar(ruleID);
		if (v!=var_Undef){
			assert(outer->value(v)==l_False);
			if (outer->level(v)>0){
				//learn v
				conflict.push(mkLit(v));
			}
		}
	}

}
void P0LAcceptDetector::printSolution(std::ostream& out){

	for(int c = 0;c<g_under.nCharacters();c++){
		int n_rules = 0;
		for(int ruleID:g_under.getRules(c)){
			if(g_under.ruleEnabled(ruleID)){
				n_rules++;
				out<<c<<" -> ";
				for(int d:g_under.getRule(ruleID)){
					out<<d<<",";
				}
				out<<"\n";
			}
		}
	}
	int depth = 5;


}
bool P0LAcceptDetector::checkSatisfied(){
	NP0LAccept check(g_under,strings);
	printSolution(std::cout);

	for(auto & t:all_lits){
		int str = t.str;
		int atom = t.atom;
		Lit l = t.l;
		vec<int> & string = strings[str];

		if(outer->value(l)==l_Undef){
			return false;
		}

		else if (outer->value(l)==l_False && check.acceptsString(str)){
			return false;
		}else if (outer->value(l)==l_True && !check.acceptsString(str)){
			return false;
		}




	}

	return true;
}


