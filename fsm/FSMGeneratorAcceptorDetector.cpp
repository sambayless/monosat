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

	if(!opt_fsm_as_graph){
		underReachStatus = new FSMGeneratorAcceptorDetector::AcceptStatus(*this, true);
		overReachStatus = new FSMGeneratorAcceptorDetector::AcceptStatus(*this, false);

		underapprox_detector = new NFALinearGeneratorAcceptor<FSMGeneratorAcceptorDetector::AcceptStatus>(g_under,acceptor_under,gen_source,accept_source,*underReachStatus);
		overapprox_detector = new NFALinearGeneratorAcceptor<FSMGeneratorAcceptorDetector::AcceptStatus>(g_over,acceptor_over,gen_source,accept_source,*overReachStatus);

		underprop_marker = outer->newReasonMarker(getID());
		overprop_marker = outer->newReasonMarker(getID());
		forcededge_marker= outer->newReasonMarker(getID());

		cur_seen.growTo(acceptor_over.states());
		gen_cur_seen.growTo(g_over.states());

		next_seen.growTo(acceptor_over.states());
		gen_next_seen.growTo(g_over.states());
		seen_chars.growTo(g_over.outAlphabet()+1);
	}else{
		graph = new GraphTheorySolver<long>(outer->getSolver());
		outer->getSolver()->addTheory(graph);

		//fully unroll the fsm you get by feeding the generator into the acceptor
		//(can improve on this by pruning some nodes if we know the final states already).

		constructAllPaths();

	}


}

void FSMGeneratorAcceptorDetector::constructAllPaths(){
	assert(graph);
	//this assumes that the generator is LINEAR!
	assert(g_over.isLinear());
	int gen_pos=0;
	cur_seen.growTo(acceptor_over.states());
	gen_cur_seen.growTo(g_over.states());

	next_seen.growTo(acceptor_over.states());
	gen_next_seen.growTo(g_over.states());
	seen_chars.growTo(g_over.outAlphabet()+1);

	nodes.growTo(g_over.states());
	for(int i = 0;i<g_over.states();i++){
		for(int j = 0;j<acceptor_over.states();j++){
			nodes[i].push(graph->newNode());
		}
	}


	assert(outer->decisionLevel()==0);


	for(int s:cur){
		assert(cur_seen);
		cur_seen[s]=false;
	}
	cur.clear();
	assert(next.size()==0);


	cur_seen[accept_source]=true;
	cur.push(accept_source);

	for(int s:gen_cur){
		assert(gen_cur_seen);
		gen_cur_seen[s]=false;
	}
	gen_cur.clear();
	assert(next.size()==0);
	gen_cur_seen[gen_source]=true;
	gen_cur.push(gen_source);
	vec<Transition> chars;
	chars.clear();
	DynamicFSM & g = acceptor_over;
	DynamicFSM & gen = g_over;

	//initial emove pass:
/*	if(g.emovesEnabled()){
		for(int i = 0;i<cur.size();i++){
			int s = cur[i];
			for(int j = 0;j<g.nIncident(s);j++){
				//now check if the label is active
				int edgeID= g.incident(s,j).id;
				int to = g.incident(s,j).node;
				int from = i;
				if(!cur_seen[to] && g.transitionEnabled(edgeID,0,0)){
					cur_seen[to]=true;
					//suffixTable[gen_pos].set(to);
					cur.push(to);
				}
				if (g.transitionEnabled(edgeID,0,0)){
					Var transitionVar = outer->getTransition(g.getID() , edgeID, 0,0).outerVar;
					Var outerVar = outer->toSolver(outer->newAuxVar());
					Lit edgeLit = mkLit(outerVar);
					graph->newEdge(nodes[gen_source][from],nodes[gen_source][to],outerVar);
					outer->makeEqualInSolver(mkLit(transitionVar),edgeLit);
				}

			}
		}
	}*/
	//initial emove pass:
	if(gen.emovesEnabled()){
		for(int i = 0;i<gen_cur.size();i++){
			int s = gen_cur[i];
			for(int j = 0;j<gen.nIncident(s);j++){
				//now check if the label is active
				int edgeID= gen.incident(s,j).id;
				int to = gen.incident(s,j).node;
				int from = i;
				if(!gen_cur_seen[to] && gen.transitionEnabled(edgeID,0,0)){
					gen_cur_seen[to]=true;
					gen_cur.push(to);

				}
				if (gen.transitionEnabled(edgeID,0,0)){
					Var transitionVar = outer->getTransition(gen.getID() , edgeID, 0,0).outerVar;
					Var outerVar =  outer->toSolver(outer->newAuxVar());
					Lit edgeLit = mkLit(outerVar);
					graph->newEdge(nodes[from][accept_source],nodes[to][accept_source],outerVar);
					outer->makeEqualInSolver(mkLit(transitionVar),edgeLit);
				}
			}
		}
	}


	int gen_prev_state = gen_source;
	//use the linear generator to produce a (set) of strings. Because the generator is linear, it is only ever in one state, which greatly simplifies the reasoning here...
	while(gen_pos<g_over.states()){
		int gen_state=0;
		stepGeneratorForward(chars,seen_chars,gen_state);//get set of next strings

		if(chars.size()==0){
			for(int i = 0;i<cur.size();i++){
				int s = cur[i];
				for(int j = 0;j<g.nIncident(s);j++){
					//now check if the label is active
					int edgeID= g.incident(s,j).id;
					int to = g.incident(s,j).node;
					int from = s;
					if(!cur_seen[to] && g.transitionEnabled(edgeID,0,0)){
						cur_seen[to]=true;
						cur.push(to);
						//suffixTable[gen_pos].set(to);


					}
					if (g.transitionEnabled(edgeID,0,0)){
						Var transitionVar = outer->getTransition(g.getID() , edgeID, 0,0).outerVar;
						Var outerVar = outer->toSolver(outer->newAuxVar());
						Lit edgeLit = mkLit(outerVar);
						graph->newEdge(nodes[gen_prev_state][from],nodes[gen_prev_state][to],outerVar);
						outer->makeEqualInSolver(mkLit(transitionVar),edgeLit);
					}
				}
			}
		}else{

				for(int i = 0;i<cur.size();i++){
					int s = cur[i];
					for(int j = 0;j<g.nIncident(s);j++){
						//now check if the label is active
						int edgeID= g.incident(s,j).id;
						int to = g.incident(s,j).node;
						int from = s;
						if(!cur_seen[to] && g.transitionEnabled(edgeID,0,0)){
							cur_seen[to]=true;
							cur.push(to);
						}
						if ( g.transitionEnabled(edgeID,0,0)){
							Var transitionVar = outer->getTransition(g.getID() , edgeID, 0,0).outerVar;
							Var outerVar = outer->toSolver(outer->newAuxVar());
							Lit edgeLit = mkLit(outerVar);
							graph->newEdge(nodes[gen_prev_state][from],nodes[gen_prev_state][to],outerVar);
							outer->makeEqualInSolver(mkLit(transitionVar),edgeLit);
						}
						for(auto transition:chars)
						{
							//int genEdge = transition.edgeID;
							int l = transition.out;
							assert(l>0);
							if (!next_seen[to] && g.transitionEnabled(edgeID,l,0)){
								next_seen[to]=true;
								next.push(to);
							}

							if (g.transitionEnabled(edgeID,l,0)){
								Lit genLit = mkLit( outer->getTransition(gen.getID() , transition.edgeID, transition.in,transition.out).outerVar);
								Lit acceptLit = mkLit( outer->getTransition(g.getID() ,edgeID,l,0).outerVar);
								Var outerVar = outer->toSolver(outer->newAuxVar());
								Lit edgeLit = mkLit(outerVar);
								Var inner_edge = var(graph->newEdge(nodes[gen_prev_state][from],nodes[gen_state][to],outerVar));
								//EdgeLit == (genVar && acceptVar)
								outer->getSolver()->addClause(genLit,~edgeLit);
								outer->getSolver()->addClause(acceptLit,~edgeLit);
								outer->getSolver()->addClause(edgeLit,~genLit,~acceptLit);
							}

						}
					}
				}

		}

		next.swap(cur);
		next_seen.swap(cur_seen);

		for(int s:next){
			assert(next_seen[s]);
			next_seen[s]=false;
		}
		next.clear();
		if(chars.size()==0){
			//must eventually happen because the generator is linear.
			break;
		}

		for(auto c :chars){
			assert(seen_chars[c.out]);
			seen_chars[c.out]=false;
		}
		chars.clear();
		gen_prev_state = gen_state;
		gen_pos++;
	}
	//g_over.draw();
	//acceptor_over.draw();
	//graph->drawFull();
}
void FSMGeneratorAcceptorDetector::stepGeneratorForward(vec<Transition> & store, vec<bool> & store_seen, int & cur_gen_state){
		DynamicFSM & g = g_over;
		cur_gen_state=0;

		for(int i = 0;i<gen_cur.size();i++){
			int s = gen_cur[i];

			for(int j = 0;j<g.nIncident(s);j++){
				//now check if the label is active
				int edgeID= g.incident(s,j).id;
				int to = g.incident(s,j).node;

				if(g.transitionEnabled(edgeID,0,0)){
					if(!gen_cur_seen[to] ){
						gen_cur_seen[to]=true;
						gen_cur.push(to);
					}
					cur_gen_state = to;
				}
				for(int l = 1;l<g.outAlphabet();l++){
					if (g.transitionEnabled(edgeID,0,l)){
						if(!gen_next_seen[to]){
							gen_next_seen[to]=true;
							gen_next.push(to);
						}
						cur_gen_state = to;
						if(!store_seen[l]){
							store_seen[l]=true;
							store.push({edgeID,0,l});
						}

					}
				}
			}
		}

		gen_next.swap(gen_cur);
		gen_next_seen.swap(gen_cur_seen);

		for(int s:gen_next){
			assert(gen_next_seen[s]);
			gen_next_seen[s]=false;
		}
		gen_next.clear();



	}
void FSMGeneratorAcceptorDetector::preprocess(){
	if(graph){
		graph->implementConstraints();
	}
}
void FSMGeneratorAcceptorDetector::addAcceptLit(int generatorFinalState, int acceptorFinalState, Var outer_reach_var){

	if(opt_fsm_as_graph){
		assert(graph);
		graph->reaches(nodes[gen_source][accept_source],nodes[generatorFinalState][acceptorFinalState], outer_reach_var);
		return;
	}

	lit_backward_map.growTo(g_over.states() * acceptor_over.states(),var_Undef);

	if(lit_backward_map[generatorFinalState + acceptorFinalState*g_over.states()] !=var_Undef){
		Var v = lit_backward_map[generatorFinalState*g_over.states() + acceptorFinalState];
		outer->makeEqualInSolver(mkLit(outer->toSolver(v)),mkLit( outer_reach_var));
		return;
	}

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

	lit_backward_map[generatorFinalState + acceptorFinalState*g_over.states()]= var(acceptLit);

	while (accept_lit_map.size() <= accept_var - first_var) {
        accept_lit_map.push({{-1},-1});
	}
	accept_lit_map[accept_var - first_var] = {acceptLit,generatorFinalState,acceptorFinalState};
	all_accept_lits.push( {acceptLit,generatorFinalState,acceptorFinalState});

	updatePrefixTable(generatorFinalState,acceptorFinalState);
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
	static vec<ForcedTransition> forced_edges;

	for(auto & t:all_accept_lits){
			forced_edges.clear();
			Lit l =t.l;
			Lit lit =t.l;
			int gen_to = t.gen_to;
			int accept_to = t.accept_to;
			//assert(is_changed[indexOf(var(l))]);
			//g_over.draw(gen_source,gen_to);
			if(!opt_fsm_negate_underapprox && outer->value(l)!=l_True && underapprox_detector->accepts(gen_to,accept_to,false,&forced_edges)){
				if (outer->value(l) == l_True) {
					//do nothing
				} else if (outer->value(l) == l_Undef) {
					outer->enqueue(l, underprop_marker);
				}else{
					conflict.push(l);
					buildAcceptReason(gen_to,accept_to, conflict);
					return false;
				}
			}else if(opt_fsm_negate_underapprox && outer->value(l)!=l_True && !overapprox_detector->accepts(gen_to,accept_to,true,&forced_edges)){
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

				if (outer->value(~l) == l_True) {
					//do nothing
				} else if (outer->value(~l) == l_Undef) {
					outer->enqueue(~l, overprop_marker);
				}else{
					conflict.push(~l);
					buildNonAcceptReason(gen_to,accept_to, conflict);
					return false;
				}
			}else{

			}
			if(outer->value(lit)==l_False){
				for(auto & t:forced_edges){
					//int edgeID = t.edgeID;
					//int input = t.input;
					//int output = t.output;
					int generator_state = t.generator_state;

					int label = t.character;

					for(int i = 0;i<g_over.nIncident(generator_state);i++){
						int edgeID = g_over.incident(generator_state,i).id;
						if(g_over.transitionEnabled(edgeID,0,label)){
							Var v = outer->getTransitionVar(g_over.getID(),edgeID,0,label);
							Lit f = mkLit(v,true);
							//Var test = var(lit);
							assert(outer->value(lit)==l_False);
							//int lev = outer->level(var(lit));
							setForcedVar(v,~lit);
							if(outer->value(f)==l_Undef){
								outer->enqueue(f, forcededge_marker);
							}else if(outer->value(f)==l_False){
								conflict.push(f);
								buildForcedEdgeReason(gen_to,accept_to,edgeID,label, conflict);
							}
						}
					}



				}
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
	}else if (marker ==forcededge_marker){
		reason.push(p);
		Var v = var(p);

		Lit forL = getForcedVar(v);
		//lbool val = outer->value(forL);
		Var forV = var(forL);
		//lbool val2 = outer->value(forV);
		assert(outer->value(forL)==l_True);
		assert(outer->value(forV)==l_False);
		reason.push(~forL);
		int gen_final = getGeneratorFinal(forV);
		int accept_final = getAcceptorFinal(forV);
		int forcedEdge = outer->getEdgeID(v);
		int forcedLabel = outer->getOutput(v);
		buildForcedEdgeReason(gen_final,accept_final,forcedEdge,  forcedLabel,  reason);
	}else{
		assert(false);
	}
}

void FSMGeneratorAcceptorDetector::buildAcceptReason(int genFinal, int acceptFinal, vec<Lit> & conflict){
	static int iter = 0;
	++iter;
//find a path - ideally, the one that traverses the fewest unique transitions - from source to node, learn that one of the transitions on that path must be disabled.
	if(!opt_fsm_negate_underapprox){
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
	}else{
		static vec<ForcedTransition> forced_edges;
		forced_edges.clear();
		assert(!overapprox_detector->accepts(genFinal,acceptFinal,true,&forced_edges));
		assert( !overapprox_detector->accepts(genFinal,acceptFinal,true));
		//run an NFA to find all transitions that accepting prefixes use.
		//printf("conflict %d\n",iter);
		//g_over.draw(gen_source,genFinal);
		static vec<NFATransition> path;
		path.clear();
		overapprox_detector->getGeneratorPath(genFinal,acceptFinal,path,false,true);

		for(auto & t:path){
			//printf("%d(c%d), ", t.edgeID, t.output);
			int edgeID = t.edgeID;
			int input = t.input;
			assert(input==0);
			int output = t.output;
			Var v = outer->getTransitionVar(g_over.getID(),edgeID,0,output);
			assert(outer->value(v)!=l_False);
			stats_forced_edges++;
			if(outer->value(v)==l_True){
				conflict.push(mkLit(v,true));
			}
		}
		//printf("\n");

	}
}
bool FSMGeneratorAcceptorDetector::stepGenerator(int final,int forcedEdge,int forcedLabel, vec<int> & store, vec<bool> & store_seen, int & cur_gen_state, vec<NFATransition> * path){
	DynamicFSM & g = g_over;


	for(int i = 0;i<gen_cur.size();i++){
		int s = gen_cur[i];
		cur_gen_state = s;
		for(int j = 0;j<g.nIncident(s);j++){
			//now check if the label is active
			int edgeID= g.incident(s,j).id;
			int to = g.incident(s,j).node;

			if(g.transitionEnabled(edgeID,0,0)){
				if(!gen_cur_seen[to] ){
					gen_cur_seen[to]=true;
					gen_cur.push(to);
				}
				if(path){
					path->push({edgeID,0,0});
				}
			}

			for(int l = 1;l<g.outAlphabet();l++){
				if (g.transitionEnabled(edgeID,0,l) || (edgeID==forcedEdge && l==forcedLabel)){
					if(!gen_next_seen[to]){
						gen_next_seen[to]=true;
						gen_next.push(to);
					}
					if(path && !( (edgeID==forcedEdge && l==forcedLabel))){
						path->push({edgeID,0,l});
					}
					if(!store_seen[l]){
						store_seen[l]=true;
						store.push(l);
					}

				}
			}

		}
	}

	gen_next.swap(gen_cur);
	gen_next_seen.swap(gen_cur_seen);

	for(int s:gen_next){
		assert(gen_next_seen[s]);
		gen_next_seen[s]=false;
	}
	gen_next.clear();

	return gen_cur_seen[final];

}
bool FSMGeneratorAcceptorDetector::isAttractor(int acceptorState){
	DynamicFSM & accept = acceptor_over;
		if(acceptorState<0){
			return true;
		}else{
			//should really fix this to work correctly for epsilon transitions between multiple acceptor states...
			for(int i = 0;i<accept.nIncident(acceptorState);i++){
				int edgeID = accept.incident(acceptorState,i).id;
				int to =  accept.incident(acceptorState,i).node;
				if(to==acceptorState){
					for(int c = 1;c<accept.inAlphabet();c++){
						if(!accept.transitionEnabled(edgeID,c,-1)){
							return false;
						}
					}
					return true;
				}
			}
		}
		return false;
	}
bool FSMGeneratorAcceptorDetector::find_gen_path(int gen_final, int accept_final,int forcedEdge,int forcedLabel, vec<NFATransition> & path,bool invertAcceptance , bool all_paths){
		bool accepting_state_is_attractor= !invertAcceptance && isAttractor(accept_final) ;
		path.clear();
		for(int s:cur){
			assert(cur_seen);
			cur_seen[s]=false;
		}
		cur.clear();
		assert(next.size()==0);
		cur_seen[accept_source]=true;
		cur.push(accept_source);


		for(int s:gen_cur){
			assert(gen_cur_seen);
			gen_cur_seen[s]=false;
		}
		gen_cur.clear();
		assert(next.size()==0);
		gen_cur_seen[gen_source]=true;
		gen_cur.push(gen_source);
		chars.clear();
		DynamicFSM & g =  acceptor_over;
		DynamicFSM & gen =  g_over;
		bool any_non_acceptors=accept_source!=accept_final;
		bool any_non_acceptors_next=false;
		//initial emove pass:
		if(g.emovesEnabled()){
			for(int i = 0;i<cur.size();i++){
				int s = cur[i];
				for(int j = 0;j<g.nIncident(s);j++){
					//now check if the label is active
					int edgeID= g.incident(s,j).id;
					int to = g.incident(s,j).node;
					if(!cur_seen[to] && g.transitionEnabled(edgeID,0,0)){
						cur_seen[to]=true;
						cur.push(to);
						if(to!=accept_final)
							any_non_acceptors=true;
					}

				}
			}
		}
		//initial emove pass:
		if(gen.emovesEnabled()){
			for(int i = 0;i<gen_cur.size();i++){
				int s = gen_cur[i];
				for(int j = 0;j<gen.nIncident(s);j++){
					//now check if the label is active
					int edgeID= gen.incident(s,j).id;
					int to = gen.incident(s,j).node;
					if(gen.transitionEnabled(edgeID,0,0)){
						if(!gen_cur_seen[to]){
							gen_cur_seen[to]=true;
							gen_cur.push(to);
						}
						path.push({edgeID,0,0});
					}
				}
			}
		}

		bool prev_accepting=accepting_state_is_attractor ? true:gen_cur_seen[gen_final];
		bool accepted=false;
		//use the linear generator to produce a (set) of strings. Because the generator is linear, it is only ever in one state, which greatly simplifies the reasoning here...
		while(!accepted){
			int prev_path_size = path.size();
			int  cur_gen_state=0;
			bool accepting = stepGenerator(gen_final, forcedEdge, forcedLabel, chars,seen_chars, cur_gen_state,& path);//get set of next strings
			if(accepting_state_is_attractor){
				accepting =true;
			}

			if(chars.size()==0){
				for(int i = 0;i<cur.size();i++){
					int s = cur[i];
					for(int j = 0;j<g.nIncident(s);j++){
						//now check if the label is active
						int edgeID= g.incident(s,j).id;
						int to = g.incident(s,j).node;
						if(!cur_seen[to] && g.transitionEnabled(edgeID,0,0)){
							cur_seen[to]=true;
							cur.push(to);
							if(to!=accept_final)
								any_non_acceptors=true;

						}
					}
				}
			}else{
				for(int l:chars)
				{
					assert(l>0);
					for(int i = 0;i<cur.size();i++){
						int s = cur[i];
						for(int j = 0;j<g.nIncident(s);j++){
							//now check if the label is active
							int edgeID= g.incident(s,j).id;
							int to = g.incident(s,j).node;
							if(!cur_seen[to] && g.transitionEnabled(edgeID,0,0)){
								cur_seen[to]=true;
								cur.push(to);
								if(to!=accept_final)
									any_non_acceptors=true;
								//status.reaches(str,to,edgeID,0);
							}

							if (!next_seen[to] && g.transitionEnabled(edgeID,l,0)){
								//status.reaches(str,to,edgeID,l);
								next_seen[to]=true;
								next.push(to);
								if(to!=accept_final)
									any_non_acceptors_next=true;
							}
						}
					}
				}
			}
			if(!invertAcceptance && ! all_paths){
				if(prev_accepting && cur_seen[accept_final]){
					accepted=true;
					path.shrink(path.size()-prev_path_size);
				}
			}else if (invertAcceptance && ! all_paths){
				if(prev_accepting && any_non_acceptors){
					accepted=true;
					path.shrink(path.size()-prev_path_size);
				}
			}else if (all_paths){
				if(prev_accepting && !any_non_acceptors){
					accepted=true;
					path.shrink(path.size()-prev_path_size);
				}
			}

			next.swap(cur);
			next_seen.swap(cur_seen);

			for(int s:next){
				assert(next_seen[s]);
				next_seen[s]=false;
			}
			next.clear();
			if(chars.size()==0){
				//must eventually happen because the generator is linear.
				break;
			}

			for(int l :chars){
				assert(seen_chars[l]);
				seen_chars[l]=false;
			}
			chars.clear();
			prev_accepting = accepting;
			any_non_acceptors= any_non_acceptors_next;
			any_non_acceptors_next=false;
		}

		return accepted;
	}


void FSMGeneratorAcceptorDetector::buildForcedEdgeReason(int genFinal, int acceptFinal,int forcedEdge, int forcedLabel, vec<Lit> & conflict){
	//the reason that an edge was forced by an _unreachability_ condition, was that if it was negated then that state would be reachable.
	//so the 'reason' is the same as the reachability reason for that state IF the forced edge was flipped.

	if(!opt_fsm_negate_underapprox){
		assert(false);
		throw std::logic_error("Bad fsm option");
	}else{

		static vec<NFATransition> path;
		path.clear();
		find_gen_path(genFinal,acceptFinal,forcedEdge,forcedLabel,path,false,true);

		for(auto & t:path){
			//printf("%d(c%d), ", t.edgeID, t.output);
			int edgeID = t.edgeID;
			int input = t.input;
			assert(input==0);
			int output = t.output;
			Var v = outer->getTransitionVar(g_over.getID(),edgeID,0,output);
			assert(outer->value(v)!=l_False);
			stats_forced_edges++;
			if(outer->value(v)==l_True){
				conflict.push(mkLit(v,true));
			}
		}
		//printf("\n");

	}
}

void FSMGeneratorAcceptorDetector::buildNonAcceptReason(int genFinal, int acceptFinal, vec<Lit> & conflict){

	static vec<NFATransition> path;
	path.clear();

	assert( !overapprox_detector->accepts(genFinal,acceptFinal,false));
	//run an NFA to find all transitions that accepting prefixes use.
	//printf("conflict %d\n",iter);
	//g_over.draw(gen_source,genFinal);
	//acceptor_over.draw(accept_source,acceptFinal);
	overapprox_detector->getGeneratorPath(genFinal,acceptFinal,path,true,true);
	static vec<bool> seen_states;
	/*seen_states.clear();
	seen_states.growTo(g_over.states());
	if(true||g_over.mustBeDeterministic()){

		//this reasoning is only correct if the SAT solver itself is already enforcing that atleast one outgoing edge from each state must be enabled.
		//otherwise, it is incorrect.

		for(auto & t:path){
			//printf("%d(c%d), ", t.edgeID, t.output);
			int edgeID = t.edgeID;
			if(edgeID==10){
				int a=1;
			}
			int input = t.input;
			assert(input==0);
			int output = t.output;
			Var v = outer->getTransitionVar(g_over.getID(),edgeID,input,output);
			assert(outer->value(v)!=l_False);
			if(outer->value(v)==l_True){
				if(outer->level(v)==0){
					//this edge cannot be disabled
				}else{
					conflict.push(mkLit(v,true));
				}
			}
		}

	}else{
		for(auto & t:path){
			//printf("%d(c%d), ", t.edgeID, t.output);
			int edgeID = t.edgeID;
			int input = t.input;
			assert(input==0);
			int output = t.output;
			Var v = outer->getTransitionVar(g_over.getID(),edgeID,0,output);
			assert(outer->value(v)!=l_False);
			//for each transition(from,to) that leads to a non-accepting state, learn that
			//at least one currently disabled edge (from,anywhere) must be enabled.

			if(outer->value(v)==l_True && outer->level(v)==0){
				//this edge cannot be disabled
			}else{

				int from = g_over.getEdge(edgeID).from;
				if(seen_states[from]){
					continue;
				}
				seen_states[from]=true;

				for (int j = 0;j<g_over.nIncident(from);j++){
					int edgeID = g_over.incident(from,j).id;
					for(int i = 0;i<g_over.inAlphabet();i++){
						for (int o = 0;o<g_over.outAlphabet();o++){

							if(!g_over.transitionEnabled(edgeID,i,o)  ){
								Var v2 =  outer->getTransition(g_over.getID(),edgeID,i,o).v;
								if(v2!=var_Undef){
									conflict.push(mkLit(v2,false));
								}
							}
						}
					}
				}
			}
		}
	}
*/


	buildSuffixCut(genFinal,acceptFinal,conflict,false,false);

}

void FSMGeneratorAcceptorDetector::updatePrefixTable(int gen_final, int accept_final){
	if(last_prefix_update<0 && outer->decisionLevel()==0){
		//(int gen_final,int accept_final,vec<Bitset> & suffixTable, bool accepting_state_is_attractor, bool invertAcceptance){
		Var v = getDetectorVar(gen_final,accept_final);
		assert(v!=var_Undef);
		int index = v - first_var;
		prefixTables.growTo(index+1);
		overapprox_detector->buildPrefixTable(gen_final,accept_final,prefixTables[index],false,false);
		last_prefix_update = g_over.modifications;
	}
}

bool FSMGeneratorAcceptorDetector::buildSuffixCut(int gen_final,int accept_final,vec<Lit> & cut, bool accepting_state_is_attractor, bool invertAcceptance){
		//run the nfa backwards from the end of the generator, and collect the set of reachable fsa states at each step in the (linear) generator.
		DynamicFSM & gen = g_over;
		DynamicFSM & accept = acceptor_over;
		cur_seen.growTo(accept.states());
		gen_cur_seen.growTo(gen.states());

		next_seen.growTo(accept.states());
		gen_next_seen.growTo(gen.states());
		seen_chars.growTo(gen.outAlphabet()+1);
		vec<Bitset> & prefixTable = getPrefixTable(gen_final,accept_final);

		for(int s:cur){
			assert(cur_seen);
			cur_seen[s]=false;
		}
		cur.clear();
		assert(next.size()==0);
		int gen_pos = gen.states()-1;
		if(!invertAcceptance){
			cur_seen[accept_final]=true;
			cur.push(accept_final);
		}else{
			for(int i = 0;i<accept.states();i++){
				if(i!=accept_final){
					cur_seen[i]=true;
					cur.push(i);
				}
			}
		}

		for(int s:gen_cur){
			assert(gen_cur_seen);
			gen_cur_seen[s]=false;
		}
		gen_cur.clear();
		assert(next.size()==0);
		gen_cur_seen[gen_final]=true;
		gen_cur.push(gen_final);


		chars.clear();
		DynamicFSM & g = accept;
		bool any_non_acceptors=!cur_seen[accept_final];
		bool any_non_acceptors_next=false;
		//initial emove pass:
		if(g.emovesEnabled()){
			for(int i = 0;i<cur.size();i++){
				int s = cur[i];
				for(int j = 0;j<g.nIncoming(s);j++){
					//now check if the label is active
					int edgeID= g.incoming(s,j).id;
					int to = g.incoming(s,j).node;
					if(!cur_seen[to] && g.transitionEnabled(edgeID,0,0)){
						cur_seen[to]=true;
						cur.push(to);
						if(to!=accept_final)
							any_non_acceptors=true;
					}else if (!g.transitionEnabled(edgeID,0,0)){
						Var v = outer->getTransitionVar(g.getID(),edgeID,0,0);
						if(v!=var_Undef &&  outer->value(v)==l_False && outer->level(v)>0){
							cut.push(mkLit(v));
						}
					}

				}
			}
		}
		//initial emove pass:
		if(gen.emovesEnabled()){
			for(int i = 0;i<gen_cur.size();i++){
				int s = gen_cur[i];
				for(int j = 0;j<gen.nIncoming(s);j++){
					//now check if the label is active
					int edgeID= gen.incoming(s,j).id;
					int to = gen.incoming(s,j).node;
					if(!gen_cur_seen[to] && gen.transitionEnabled(edgeID,0,0)){
						gen_cur_seen[to]=true;
						gen_cur.push(to);
					}else if (!gen.transitionEnabled(edgeID,0,0)){// && prefixTable[gen_final][to]
						Var v = outer->getTransitionVar(gen.getID(),edgeID,0,0);
						if(v!=var_Undef &&  outer->value(v)==l_False && outer->level(v)>0){
							cut.push(mkLit(v));
						}
					}
				}
			}
		}

		bool prev_accepting=accepting_state_is_attractor ? true:gen_cur_seen[gen_source];
		bool accepted=false;


		//use the linear generator to produce a (set) of strings. Because the generator is linear, it is only ever in one state, which greatly simplifies the reasoning here...
		while(!accepted){
			bool accepting = stepGeneratorBackward(gen_final, prefixTable,cut, chars,seen_chars);//get set of next strings
			if(accepting_state_is_attractor){
				accepting =true;
			}
			if(chars.size()==0){
				for(int i = 0;i<cur.size();i++){
					int s = cur[i];
					for(int j = 0;j<g.nIncoming(s);j++){
						//now check if the label is active
						int edgeID= g.incoming(s,j).id;
						int to = g.incoming(s,j).node;
						if(!cur_seen[to] && g.transitionEnabled(edgeID,0,0)){
							cur_seen[to]=true;
							cur.push(to);

							if(to!=accept_final)
								any_non_acceptors=true;

						}else if (!g.transitionEnabled(edgeID,0,0)){
							Var v = outer->getTransitionVar(g.getID(),edgeID,0,0);
							if(v!=var_Undef &&  outer->value(v)==l_False && outer->level(v)>0){
								cut.push(mkLit(v));
							}
						}
					}
				}
			}else{
				for(int l:chars)
				{
					assert(l>0);
					for(int i = 0;i<cur.size();i++){
						int s = cur[i];
						for(int j = 0;j<g.nIncoming(s);j++){
							//now check if the label is active
							int edgeID= g.incoming(s,j).id;
							int to = g.incoming(s,j).node;
							if(!cur_seen[to] && g.transitionEnabled(edgeID,0,0)){
								cur_seen[to]=true;
								cur.push(to);

								if(to!=accept_final)
									any_non_acceptors=true;
								//status.reaches(str,to,edgeID,0);
							}else if (!g.transitionEnabled(edgeID,0,0)){
								Var v = outer->getTransitionVar(g.getID(),edgeID,0,0);
								if(v!=var_Undef &&  outer->value(v)==l_False && outer->level(v)>0){
									cut.push(mkLit(v));
								}
							}

							if (!next_seen[to] && g.transitionEnabled(edgeID,l,0)){

								//status.reaches(str,to,edgeID,l);
								next_seen[to]=true;
								next.push(to);
								if(to!=accept_final)
									any_non_acceptors_next=true;
							}else if (!g.transitionEnabled(edgeID,l,0)){
								Var v = outer->getTransitionVar(g.getID(),edgeID,l,0);
								if(v!=var_Undef &&  outer->value(v)==l_False && outer->level(v)>0){
									cut.push(mkLit(v));
								}
							}
						}
					}
				}
			}
			if(prev_accepting && cur_seen[accept_source]){
					accepted=true;
				}
			/*if(!invertAcceptance){

			}else{
				if(prev_accepting && any_non_acceptors){
					accepted=true;
				}
			}*/
			next.swap(cur);
			next_seen.swap(cur_seen);

			for(int s:next){
				assert(next_seen[s]);
				next_seen[s]=false;
			}
			next.clear();
			if(chars.size()==0){
				//must eventually happen because the generator is linear.
				break;
			}

			for(int l :chars){
				assert(seen_chars[l]);
				seen_chars[l]=false;
			}
			chars.clear();
			prev_accepting = accepting;
			any_non_acceptors= any_non_acceptors_next;
			any_non_acceptors_next=false;
			gen_pos--;
		}

		return accepted;

	}

bool FSMGeneratorAcceptorDetector::stepGeneratorBackward(int final, vec<Bitset> & prefixTable, vec<Lit> & cut, vec<int> & store, vec<bool> & store_seen, vec<NFATransition> * path){
		DynamicFSM & g = g_over;

		for(int i = 0;i<gen_cur.size();i++){
			int s = gen_cur[i];
			for(int j = 0;j<g.nIncoming(s);j++){
				//now check if the label is active
				int edgeID= g.incoming(s,j).id;
				int to = g.incoming(s,j).node;

				if(g.transitionEnabled(edgeID,0,0)){
					if(!gen_cur_seen[to] ){
						gen_cur_seen[to]=true;
						gen_cur.push(to);
					}
					if(path){
						path->push({edgeID,0,0});
					}
				}else {//if (prefixTable[s][to])
					Var v = outer->getTransitionVar(g.getID(),edgeID,0,0);
					if(v!=var_Undef &&  outer->value(v)==l_False && outer->level(v)>0){
						cut.push(mkLit(v));
					}
				}

				int edge_assigned_true=-1;
				for(int l = 1;l<g_under.outAlphabet();l++){
					if (g_under.transitionEnabled(edgeID,0,l)){
						edge_assigned_true=l;
						break;
					}
				}
				if(edge_assigned_true>=0){
					//if a character has been decided true, then since exactly one character is learnt in each position, we can just learn !character, rather than learning the negation of the set of disabled characters
					Var v = outer->getTransitionVar(g.getID(),edgeID,0,edge_assigned_true);
					if(v!=var_Undef &&  outer->value(v)==l_True && outer->level(v)>0){
						cut.push(mkLit(v,true));
					}
				}

				for(int l = 1;l<g.outAlphabet();l++){
					//g.mustBeDeterministic() &&
					 if (g.transitionEnabled(edgeID,0,l)){
						if(!gen_next_seen[to]){
							gen_next_seen[to]=true;
							gen_next.push(to);
						}
						if(path){
							path->push({edgeID,0,l});
						}
						if(!store_seen[l]){
							store_seen[l]=true;
							store.push(l);
						}

					}else if (edge_assigned_true==-1) {//if (prefixTable[s][to])
						Var v = outer->getTransitionVar(g.getID(),edgeID,0,l);
						if(v!=var_Undef &&  outer->value(v)==l_False && outer->level(v)>0){
							cut.push(mkLit(v));
						}
					}
				}

			}
		}

		gen_next.swap(gen_cur);
		gen_next_seen.swap(gen_cur_seen);

		for(int s:gen_next){
			assert(gen_next_seen[s]);
			gen_next_seen[s]=false;
		}
		gen_next.clear();

		return gen_cur_seen[gen_source];

	}


Lit FSMGeneratorAcceptorDetector::decide(int level) {


	double startdecidetime = rtime(2);
	for(auto & t:all_accept_lits){

		Lit l =t.l;
		if(outer->value(l)==l_False){
			int gen_to = t.gen_to;
			int accept_to = t.accept_to;
			static vec<NFATransition> path;
			path.clear();

			if(overapprox_detector->getGeneratorPath(gen_to,accept_to,path,true,false)){
				for(auto & t:path){
					Var v = outer->getTransitionVar(g_over.getID(),t.edgeID,t.input,t.output);
					if(outer->value(v)==l_Undef){
						stats_decisions++;
						stats_decide_time += rtime(2) - startdecidetime;
						return mkLit(v,false);
					}
				}
			}
		}
	}
	stats_decide_time += rtime(2) - startdecidetime;
	return lit_Undef;
}

void FSMGeneratorAcceptorDetector::printSolution(std::ostream& out){
	if(graph){
		//graph->drawCurrent();
		//graph->printSolution();
		return;
	}

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


