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
#ifndef FSM_GENERATORACCEPTDETECTOR_H_
#define FSM_GENERATORACCEPTDETECTOR_H_
#include "utils/System.h"

#include "dgl/DynamicGraph.h"

#include "DynamicFSM.h"

#include "core/SolverTypes.h"
#include "mtl/Map.h"


#include "utils/System.h"
#include "FSMDetector.h"
#include "alg/NFALinearGeneratorAcceptor.h"
#include "graph/GraphTheory.h"
using namespace dgl;
namespace Monosat {

class FSMTheorySolver;

class FSMGeneratorAcceptorDetector: public FSMDetector {
public:
	FSMTheorySolver* outer;
	DynamicFSM &g_under;
	DynamicFSM &g_over;
	DynamicFSM &acceptor_under;
	DynamicFSM &acceptor_over;

	int first_destination=-1;
	int gen_source;
	int accept_source;
	double rnd_seed;

	struct AcceptStatus {
		FSMGeneratorAcceptorDetector & detector;
		bool polarity;
		void accepts(int string, int state,int edgeID,int label, bool accepts);

		AcceptStatus(FSMGeneratorAcceptorDetector & _outer, bool _polarity) :
				detector(_outer), polarity(_polarity) {
		}
	};

	AcceptStatus *underReachStatus = nullptr;
	AcceptStatus *overReachStatus = nullptr;

	NFALinearGeneratorAcceptor<AcceptStatus> * underapprox_detector = nullptr;
	NFALinearGeneratorAcceptor<AcceptStatus> * overapprox_detector = nullptr;



	CRef underprop_marker;
	CRef overprop_marker;
	CRef forcededge_marker;
	struct Change {
		Lit l;
		int u;
		int str;
	};
	//vec<bool> is_changed;
	vec<Change> changed;

	Var first_var=var_Undef;
	struct AcceptLit{
		Lit l;
		int gen_to;
		int accept_to;;
	};
	vec<AcceptLit> accept_lit_map;
	vec<AcceptLit> all_accept_lits;
	vec<Lit> all_lits;

	//stats
	
	int stats_full_updates = 0;
	int stats_fast_updates = 0;
	int stats_fast_failed_updates = 0;
	int stats_skip_deletes = 0;
	int stats_skipped_updates = 0;
	int stats_num_skipable_deletions = 0;
	int stats_learnt_components = 0;
	int stats_learnt_components_sz = 0;
	long stats_forced_edges=0;
	double mod_percentage = 0.2;
	int stats_pure_skipped = 0;
	int stats_shrink_removed = 0;
	double stats_full_update_time = 0;
	double stats_fast_update_time = 0;

	Map<Var,Lit> forcedVars;
	vec<Var> lit_backward_map;

	GraphTheorySolver<long> * graph=nullptr;//for if we reduce the nfa to a graph
	vec<vec<int> > nodes;

	void printStats() {
		//printf("Reach detector\n");
		if(graph){
			graph->printStats(1);
		}else{
			FSMDetector::printStats();
			if (opt_detect_pure_theory_lits)
				printf("\tPropagations skipped by pure literal detection: %d\n", stats_pure_skipped);
			if (opt_shrink_theory_conflicts) {
				printf("\t%d lits removed by shrinking conflicts\n", stats_shrink_removed);
			}
			if (opt_learn_unreachable_component) {
				printf("\t%d components learned, average component size: %f\n", stats_learnt_components,
						stats_learnt_components_sz / (float) stats_learnt_components);
			}
			printf("Forced edge assignments: %ld\n", stats_forced_edges);
		}
	}

	inline void setForcedVar(Var edgeVar, Lit forcedBy){
		if(forcedVars.has(edgeVar)){
			forcedVars.remove(edgeVar);
		}

	//	printf("Forced edge %d by %d, %d\n", edgeVar, var(forcedBy), sign(forcedBy)?-1:1);
		forcedVars.insert(edgeVar,forcedBy);
	}

	inline Lit getForcedVar(Var edgeVar){
		return forcedVars[edgeVar];
	}

	inline int indexOf(Var v)const{
		int index = v - first_var;
		assert(index < accept_lit_map.size());
		return index;
	}

	int getGeneratorFinal(Var reachVar) {
		assert(reachVar >= first_var);
		int index = indexOf(reachVar);
		assert(accept_lit_map[index].gen_to >= 0);
		return accept_lit_map[index].gen_to;
	}
	int getAcceptorFinal(Var reachVar) {
		assert(reachVar >= first_var);
		int index = indexOf(reachVar);

		assert(accept_lit_map[index].accept_to >= 0);
		return accept_lit_map[index].accept_to;
	}

	Lit decide(int level);
	bool propagate(vec<Lit> & conflict);
	void buildAcceptReason(int genFinal, int acceptFinal, vec<Lit> & conflict);
	void buildNonAcceptReason(int genFinal, int acceptFinal, vec<Lit> & conflict);
	void buildForcedEdgeReason(int genFinal, int acceptFinal,int forcedEdge, int forcedLabel,  vec<Lit> & conflict);
	 void preprocess() ;
	void buildReason(Lit p, vec<Lit> & reason, CRef marker);
	bool checkSatisfied();
	void printSolution(std::ostream& write_to);

	void addAcceptLit(int state, int strID, Var reach_var);

	Var getDetectorVar(int gen_final, int accept_final){
		return lit_backward_map[gen_final +accept_final*g_over.states()];
	}

	FSMGeneratorAcceptorDetector(int _detectorID, FSMTheorySolver * _outer, DynamicFSM &g_under, DynamicFSM &g_over,DynamicFSM & acceptor_under,DynamicFSM & acceptor_over,
			int gen_source,int acceptor_source, double seed = 1);
	virtual ~FSMGeneratorAcceptorDetector() {
		
	}
	
	const char* getName() {
		return "NFA Generator Acceptor Detector";
	}
private:

	struct Transition{
		int edgeID;
		int in;
		int out;
	};

	void constructAllPaths();
	void stepGeneratorForward(vec<Transition> & store, vec<bool> & store_seen, int & cur_gen_state);
	vec<int> next;
	vec<int> cur;

	vec<bool> next_seen;
	vec<bool> cur_seen;

	vec<int> gen_cur;
	vec<int> gen_next;
	vec<bool> gen_next_seen;
	vec<bool> gen_cur_seen;
	vec<int> chars;
	vec<bool> seen_chars;
	vec<vec<Bitset>>  prefixTables;

	int last_prefix_update=-1;
	bool isAttractor(int acceptorState);
	bool find_gen_path(int gen_final, int accept_final,int forcedEdge,int forcedLabel, vec<NFATransition> & path,bool invertAcceptance = false, bool all_paths=false);
	bool stepGenerator(int final,int forcedEdge,int forcedLabel, vec<int> & store, vec<bool> & store_seen, int & cur_gen_state, vec<NFATransition> * path=nullptr);
	bool buildSuffixCut(int gen_final,int accept_final,vec<Lit> & cut, bool accepting_state_is_attractor, bool invertAcceptance);
	bool stepGeneratorBackward(int final,vec<Bitset> & prefixTable, vec<Lit> & cut,  vec<int> & store, vec<bool> & store_seen, vec<NFATransition> * path=nullptr);

	void updatePrefixTable(int gen_final, int accept_final);

	vec<Bitset> & getPrefixTable(int gen_final, int accept_final){
		updatePrefixTable(gen_final,accept_final);
		Var v = getDetectorVar(gen_final,accept_final);
		assert(v!=var_Undef);
		int index = v - first_var;
		return prefixTables[index];
	}

	bool isAttractor(DynamicFSM & accept, int acceptorState){
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

};
}
;
#endif /* REACHDETECTOR_H_ */
