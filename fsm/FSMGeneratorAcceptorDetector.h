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

	NFALinearGeneratorAcceptor<AcceptStatus> * underapprox_detector;
	NFALinearGeneratorAcceptor<AcceptStatus> * overapprox_detector;



	CRef underprop_marker;
	CRef overprop_marker;

	struct Change {
		Lit l;
		int u;
		int str;
	};
	//vec<bool> is_changed;
	vec<Change> changed;

	vec<vec<Lit>> accept_lits;
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
	double mod_percentage = 0.2;
	int stats_pure_skipped = 0;
	int stats_shrink_removed = 0;
	double stats_full_update_time = 0;
	double stats_fast_update_time = 0;




	void printStats() {
		//printf("Reach detector\n");
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


	bool propagate(vec<Lit> & conflict);
	void buildAcceptReason(int genFinal, int acceptFinal, vec<Lit> & conflict);
	void buildNonAcceptReason(int genFinal, int acceptFinal, vec<Lit> & conflict);

	void buildReason(Lit p, vec<Lit> & reason, CRef marker);
	bool checkSatisfied();
	void printSolution(std::ostream& write_to);

	void addAcceptLit(int state, int strID, Var reach_var);



	FSMGeneratorAcceptorDetector(int _detectorID, FSMTheorySolver * _outer, DynamicFSM &g_under, DynamicFSM &g_over,DynamicFSM & acceptor_under,DynamicFSM & acceptor_over,
			int gen_source,int acceptor_source, double seed = 1);
	virtual ~FSMGeneratorAcceptorDetector() {
		
	}
	
	const char* getName() {
		return "NFA Generator Acceptor Detector";
	}
	
	

};
}
;
#endif /* REACHDETECTOR_H_ */
