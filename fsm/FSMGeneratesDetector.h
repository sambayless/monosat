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
#ifndef FSM_GENERATEDETECTOR_H_
#define FSM_GENERATEDETECTOR_H_
#include "utils/System.h"

#include "dgl/DynamicGraph.h"

#include "DynamicFSM.h"

#include "core/SolverTypes.h"
#include "mtl/Map.h"


#include "utils/System.h"
#include "FSMDetector.h"
#include "alg/NFAGenerate.h"

using namespace dgl;
namespace Monosat {

class FSMTheorySolver;

class FSMGeneratesDetector: public FSMDetector {
public:
	FSMTheorySolver* outer;
	DynamicFSM &g_under;
	DynamicFSM &g_over;


	int source;
	vec<vec<int>> & strings;
	double rnd_seed;

	struct GenerateStatus {
		FSMGeneratesDetector & detector;
		bool polarity;
		void generates(int string, bool generates);

		GenerateStatus(FSMGeneratesDetector & _outer, bool _polarity) :
				detector(_outer), polarity(_polarity) {
		}
	};

	GenerateStatus *underReachStatus = nullptr;
	GenerateStatus *overReachStatus = nullptr;

	NFAGenerate<GenerateStatus> * underapprox_detector;
	NFAGenerate<GenerateStatus> * overapprox_detector;



	CRef underprop_marker;
	CRef overprop_marker;

	struct Change {
		Lit l;

		int str;
	};
	//vec<bool> is_changed;
	vec<Change> changed;

	vec<Lit> generate_lits;
	Var first_var=var_Undef;
	struct GenerateLit{
		int str;
	};
	vec<GenerateLit> generate_lit_map;
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
	
	void unassign(Lit l) {
		FSMDetector::unassign(l);
		int index = indexOf(var(l));
		if (index >= 0 && index < generate_lit_map.size() && generate_lit_map[index].str != -1) {

			int str =  generate_lit_map[index].str;
			//if (!is_changed[index]) {
            changed.push( { {var(l)}, str });
			//	is_changed[index] = true;
			//}
		}
	}
	
	inline int indexOf(Var v)const{
		int index = v - first_var;
		assert(index < generate_lit_map.size());
		return index;
	}


	int getString(Var reachVar) {
		assert(reachVar >= first_var);
		int index = indexOf(reachVar);


		return generate_lit_map[index].str;
	}

	bool propagate(vec<Lit> & conflict);
	void buildGeneratesReason(int str, vec<Lit> & conflict);
	void buildNonGeneratesReason(int str, vec<Lit> & conflict);

	void buildReason(Lit p, vec<Lit> & reason, CRef marker);
	bool checkSatisfied();
	void printSolution(std::ostream& write_to);

	void addGeneratesLit( int strID, Var reach_var);



	FSMGeneratesDetector(int _detectorID, FSMTheorySolver * _outer, DynamicFSM &g_under, DynamicFSM &g_over,
			int _source, vec<vec<int>> &  strs, double seed = 1);
	virtual ~FSMGeneratesDetector() {
		
	}
	
	const char* getName() {
		return "NFA Generates Detector";
	}
	
private:

	struct UsedTransition{
		int edge=-1;
		int label=-1;
	};

	vec<UsedTransition> used_transitions;
	bool unique_path_conflict(int s,int string,int str_pos,int emove_count, vec<NFATransition> & path,vec<Lit> & conflict);
	

};
}
;
#endif /* REACHDETECTOR_H_ */
