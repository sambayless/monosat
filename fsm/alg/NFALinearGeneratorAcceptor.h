/*
 * NFAGenerateAccept.h
 *
 *  Created on: Jan 3, 2015
 *      Author: sam
 */

#ifndef NFAGENERATEACCEPT_H_
#define NFAGENERATEACCEPT_H_

#include "../DynamicFSM.h"
#include "mtl/Bitset.h"
#include "mtl/Vec.h"
#include "NFATypes.h"
using namespace Monosat;
struct ForcedTransition{
	int generator_state;
	int character;
};

template<class Status=FSMNullStatus>
class NFALinearGeneratorAcceptor{
	DynamicFSM & gen;
	DynamicFSM & accept;
	Status & status;
	int gen_last_modification=-1;
	int gen_last_addition=-1;
	int gen_last_deletion=-1;
	int gen_history_qhead=0;
	int gen_last_history_clear=0;

	int accept_last_modification=-1;
	int accept_last_addition=-1;
	int accept_last_deletion=-1;
	int accept_history_qhead=0;
	int accept_last_history_clear=0;


	int stats_full_updates=0;
	int stats_fast_updates=0;
	int stats_fast_failed_updates=0;
	int stats_skip_deletes=0;
	int stats_skipped_updates=0;
	int stats_num_skipable_deletions=0;
	double mod_percentage=0;

	double stats_full_update_time=0;
	double stats_fast_update_time=0;


	vec<int> next;
	vec<int> cur;

	vec<bool> next_seen;
	vec<bool> cur_seen;

	vec<int> gen_cur;
	vec<int> gen_next;
	vec<bool> gen_next_seen;
	vec<bool> gen_cur_seen;
	int gen_source;
	int accept_source;

	vec<Bitset> suffixTable;
	vec<int> chars;
	vec<bool> seen_chars;

	struct Check{
		int gen_final;
		int accept_final;
	};

	vec<Check> toCheck;

public:



	NFALinearGeneratorAcceptor(DynamicFSM &  gen,DynamicFSM &  accept,int gen_source, int accept_source,Status & status=fsmNullStatus):gen(gen),accept(accept),status(status),gen_source(gen_source),accept_source(accept_source){
		//Generator must be _linear_ (but may be non-deterministic.)
		cur_seen.growTo(accept.states());
		gen_cur_seen.growTo(gen.states());

		next_seen.growTo(accept.states());
		gen_next_seen.growTo(gen.states());
		seen_chars.growTo(gen.outAlphabet()+1);
		vec<int> stack;
		stack.push(gen_source);
		while(stack.size()){
			int n = stack.last();
			stack.pop();
			int to = -1;
			for(int i = 0;i<gen.nIncident(n);i++){
				int t = gen.incident(n,i).node;
				int edgeID = gen.incident(n,i).id;
				if(gen.transitionEnabled(edgeID,-1,-1)){
					if(to!=-1){
						fprintf(stderr,"Error! FSM generator is not linear!");
						exit(1);
						break;
					}
					to = t;
				}
			}
			if(to!=-1){
				stack.push(to);
			}
		}
	}

private:

	bool stepGeneratorBackward(int final, vec<int> & store, vec<bool> & store_seen, int & cur_gen_state, vec<NFATransition> * path=nullptr){
			DynamicFSM & g = gen;


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
					}
					for(int l = 1;l<g.outAlphabet();l++){
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

	bool stepGenerator(int final, vec<int> & store, vec<bool> & store_seen, int & cur_gen_state, vec<NFATransition> * path=nullptr){
		DynamicFSM & g = gen;


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

	void addCheck(int generatorFinalState, int acceptorFinalState){
		toCheck.push({generatorFinalState,acceptorFinalState});
	}
	bool isAttractor(int acceptorState){
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
public:
	bool buildSuffixTable(int gen_final,int accept_final,vec<Bitset> & suffixTable, bool accepting_state_is_attractor, bool invertAcceptance){
		//run the nfa backwards from the end of the generator, and collect the set of reachable fsa states at each step in the (linear) generator.

		suffixTable.growTo(gen.states());
		for(int i = 0;i<suffixTable.size();i++){
			suffixTable[i].clear();
			suffixTable[i].growTo(accept.states());
		}

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
			suffixTable[gen_pos].set(accept_final);
		}else{
			for(int i = 0;i<accept.states();i++){
				if(i!=accept_final){
					cur_seen[i]=true;
					cur.push(i);
					suffixTable[gen_pos].set(i);
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
		bool any_non_acceptors=!cur_seen[accept_final];//should this be 'accept_source'?
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
						suffixTable[gen_pos].set(to);
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
				for(int j = 0;j<gen.nIncoming(s);j++){
					//now check if the label is active
					int edgeID= gen.incoming(s,j).id;
					int to = gen.incoming(s,j).node;
					if(!gen_cur_seen[to] && gen.transitionEnabled(edgeID,0,0)){
						gen_cur_seen[to]=true;
						gen_cur.push(to);

					}
				}
			}
		}

		bool prev_accepting=accepting_state_is_attractor ? true:gen_cur_seen[gen_source];
		bool accepted=false;


		//use the linear generator to produce a (set) of strings. Because the generator is linear, it is only ever in one state, which greatly simplifies the reasoning here...
		while(!accepted){
			int ignore;
			bool accepting = stepGeneratorBackward(gen_final, chars,seen_chars,ignore);//get set of next strings
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
							suffixTable[gen_pos].set(to);
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
						for(int j = 0;j<g.nIncoming(s);j++){
							//now check if the label is active
							int edgeID= g.incoming(s,j).id;
							int to = g.incoming(s,j).node;
							if(!cur_seen[to] && g.transitionEnabled(edgeID,0,0)){
								cur_seen[to]=true;
								cur.push(to);
								suffixTable[gen_pos].set(to);
								if(to!=accept_final)
									any_non_acceptors=true;
								//status.reaches(str,to,edgeID,0);
							}

							if (!next_seen[to] && g.transitionEnabled(edgeID,l,0)){
								if(gen_pos>0){
									suffixTable[gen_pos-1].set(to);
								}
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
	bool buildPrefixTable(int gen_final,int accept_final,vec<Bitset> & suffixTable, bool accepting_state_is_attractor, bool invertAcceptance){
			//run the nfa backwards from the end of the generator, and collect the set of reachable fsa states at each step in the (linear) generator.

			suffixTable.growTo(gen.states());
			for(int i = 0;i<suffixTable.size();i++){
				suffixTable[i].clear();
				suffixTable[i].growTo(accept.states());
			}

			for(int s:cur){
				assert(cur_seen);
				cur_seen[s]=false;
			}
			cur.clear();
			assert(next.size()==0);
			int gen_pos = gen.states()-1;
			if(!invertAcceptance){
				cur_seen[accept_source]=true;
				cur.push(accept_source);
				suffixTable[gen_pos].set(accept_source);
			}else{
				for(int i = 0;i<accept.states();i++){
					if(i!=accept_source){
						cur_seen[i]=true;
						cur.push(i);
						suffixTable[gen_pos].set(i);
					}
				}
			}

			for(int s:gen_cur){
				assert(gen_cur_seen);
				gen_cur_seen[s]=false;
			}
			gen_cur.clear();
			assert(next.size()==0);
			gen_cur_seen[gen_source]=true;
			gen_cur.push(gen_source);
			chars.clear();
			DynamicFSM & g = accept;
			bool any_non_acceptors=!cur_seen[accept_final];//check this!
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
							suffixTable[gen_pos].set(to);
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
						if(!gen_cur_seen[to] && gen.transitionEnabled(edgeID,0,0)){
							gen_cur_seen[to]=true;
							gen_cur.push(to);

						}
					}
				}
			}

			bool prev_accepting=accepting_state_is_attractor ? true:gen_cur_seen[gen_source];
			bool accepted=false;


			//use the linear generator to produce a (set) of strings. Because the generator is linear, it is only ever in one state, which greatly simplifies the reasoning here...
			while(!accepted){
				int ignore;
				bool accepting = stepGenerator(gen_final, chars,seen_chars,ignore);//get set of next strings
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
								suffixTable[gen_pos].set(to);
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
									suffixTable[gen_pos].set(to);
									if(to!=accept_final)
										any_non_acceptors=true;
									//status.reaches(str,to,edgeID,0);
								}

								if (!next_seen[to] && g.transitionEnabled(edgeID,l,0)){
									if(gen_pos>0){
										suffixTable[gen_pos-1].set(to);
									}
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
				gen_pos++;
			}

			return accepted;

		}
private:
	bool find_accepts(int gen_final, int accept_final, bool invertAcceptance, vec<ForcedTransition> * forced_edges=nullptr){
		bool accepting_state_is_attractor= !invertAcceptance && isAttractor(accept_final) ;
		bool hasSuffix = false;

		if(forced_edges){

			if(!buildSuffixTable(gen_final, accept_final, suffixTable,accepting_state_is_attractor,invertAcceptance)){
				return false;
			}
			hasSuffix=true;
		}
		int gen_pos = 0;
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
		DynamicFSM & g = accept;
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
					if(hasSuffix && !suffixTable[gen_pos][to])
						continue;
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
					if(!gen_cur_seen[to] && gen.transitionEnabled(edgeID,0,0)){
						gen_cur_seen[to]=true;
						gen_cur.push(to);

					}
				}
			}
		}

		bool prev_accepting=accepting_state_is_attractor ? true:gen_cur_seen[gen_final];
		bool accepted=false;


		//use the linear generator to produce a (set) of strings. Because the generator is linear, it is only ever in one state, which greatly simplifies the reasoning here...
		while(!accepted){
			int  cur_gen_state=0;
			bool accepting = stepGenerator(gen_final, chars,seen_chars,cur_gen_state);//get set of next strings
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
						if(hasSuffix && !suffixTable[gen_pos][to])
							continue;
						if( g.transitionEnabled(edgeID,0,0)){
							if(!cur_seen[to]){
							cur_seen[to]=true;
							cur.push(to);
							if(to!=accept_final)
								any_non_acceptors=true;
							}

						}
					}
				}
			}else{
				for(int l:chars)
				{
					bool character_cannot_lead_to_accepting_state=true;
					assert(l>0);
					for(int i = 0;i<cur.size();i++){
						int s = cur[i];
						for(int j = 0;j<g.nIncident(s);j++){
							//now check if the label is active
							int edgeID= g.incident(s,j).id;
							int to = g.incident(s,j).node;

							if(g.transitionEnabled(edgeID,0,0)){
								if(!cur_seen[to]){
									if(!hasSuffix || suffixTable[gen_pos][to]){
										cur_seen[to]=true;
										cur.push(to);
										if(to!=accept_final)
											any_non_acceptors=true;
									}
								}

								//status.reaches(str,to,edgeID,0);
							}

							if (g.transitionEnabled(edgeID,l,0)){
								//status.reaches(str,to,edgeID,l);
								if(!hasSuffix || (gen_pos+1<suffixTable.size() && suffixTable[gen_pos+1][to])){
									if(!next_seen[to]){
										next_seen[to]=true;
										next.push(to);
									}
									character_cannot_lead_to_accepting_state=false;
									if(to!=accept_final)
										any_non_acceptors_next=true;
								}


							}
						}
					}
					if(forced_edges && character_cannot_lead_to_accepting_state){
													//this is an edge that _must_ be disabled, because it leads to a state in the nfa that cannot reach the acceptor.
													//forced_edges->push(NFATransition{edgeID,l,0});
						forced_edges->push({cur_gen_state,l});

					}

				}
			}
			if(!invertAcceptance){
				if(prev_accepting && cur_seen[accept_final]){
					accepted=true;
				}
			}else{
				if(prev_accepting && any_non_acceptors){
					accepted=true;
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
			gen_pos++;
		}

		return accepted;
	}

	bool find_gen_path(int gen_final, int accept_final, vec<NFATransition> & path,bool invertAcceptance = false, bool all_paths=false){
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
		DynamicFSM & g = accept;
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
			bool accepting = stepGenerator(gen_final, chars,seen_chars, cur_gen_state,& path);//get set of next strings
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



public:

	void update(){

		if (gen_last_modification > 0 && gen.modifications == gen_last_modification && accept_last_modification > 0 && accept.modifications == accept_last_modification) {
			stats_skipped_updates++;
			return;
		}
		static int iteration = 0;
		int local_it = ++iteration;
		stats_full_updates++;

		if (gen_last_modification <= 0 || gen.changed() || gen_last_history_clear != gen.historyclears ||
				accept_last_modification <= 0 || accept.changed() || accept_last_history_clear != accept.historyclears) {
			gen_next_seen.clear();
			gen_next_seen.growTo(gen.states());
			gen_cur_seen.clear();
			gen_cur_seen.growTo(gen.states());

			next_seen.clear();
			next_seen.growTo(accept.states());
			cur_seen.clear();
			cur_seen.growTo(accept.states());
		}

		accept_last_modification = accept.modifications;
		accept_last_deletion = accept.deletions;
		accept_last_addition = accept.additions;
		accept_history_qhead = accept.history.size();
		accept_last_history_clear = accept.historyclears;

		gen_last_modification = gen.modifications;
		gen_last_deletion = gen.deletions;
		gen_last_addition = gen.additions;
		gen_history_qhead = gen.history.size();
		gen_last_history_clear = gen.historyclears;
	}




public:


	bool accepts(int genFinal, int acceptFinal, bool invertAcceptor=false,vec<ForcedTransition> * forced_edges=nullptr){
		return find_accepts(genFinal, acceptFinal,invertAcceptor,forced_edges);

	}
	bool getGeneratorPath(int genFinal, int acceptFinal, vec<NFATransition> & path, bool invert_acceptor=false, bool all_paths = false){
		return find_gen_path(genFinal, acceptFinal,path,invert_acceptor,all_paths);
	}



};



#endif /* NFAGENERATEACCEPT_H_ */
