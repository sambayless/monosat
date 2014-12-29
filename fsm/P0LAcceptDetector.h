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
#ifndef POL_ACCEPTDETECTOR_H_
#define POL_ACCEPTDETECTOR_H_
#include "utils/System.h"

#include "dgl/DynamicGraph.h"

#include "LSystem.h"

#include "core/SolverTypes.h"
#include "mtl/Map.h"


#include "utils/System.h"
#include "FSMDetector.h"
#include "alg/NP0LAcceptor.h"

using namespace dgl;
namespace Monosat {

class LSystemSolver;

class P0LAcceptDetector: public FSMDetector {
public:
	LSystemSolver* outer;
	LSystem &g_under;
	LSystem &g_over;



	int first_destination=-1;
	vec<vec<int>> & strings;
	double rnd_seed;


	NP0LAccept * underapprox_detector;
	NP0LAccept * overapprox_detector;



	CRef underprop_marker;
	CRef overprop_marker;

	struct Change {
		Lit l;

		int str;
	};
	//vec<bool> is_changed;
	vec<Change> changed;

	vec<Lit> accept_lits;
	Var first_var=var_Undef;
	struct AcceptLit{
		int str;

	};
	vec<AcceptLit> accept_lit_map;
	struct AcceptLitPair{
		Lit l;
		int str;
	};
	vec<AcceptLitPair> all_lits;
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
		if (index >= 0 && index < accept_lit_map.size() && accept_lit_map[index].str != -1) {

			int str =  accept_lit_map[index].str;
			//if (!is_changed[index]) {
				changed.push( { var(l), str });
			//	is_changed[index] = true;
			//}
		}
	}
	
	inline int indexOf(Var v)const{
		int index = v - first_var;
		assert(index < accept_lit_map.size());
		return index;
	}

	int getString(Var reachVar) {
		assert(reachVar >= first_var);
		int index = indexOf(reachVar);

		return accept_lit_map[index].str;
	}

	bool propagate(vec<Lit> & conflict);
	void buildAcceptReason(int str, vec<Lit> & conflict);
	void buildNonAcceptReason(int str, vec<Lit> & conflict);

	void buildReason(Lit p, vec<Lit> & reason, CRef marker);
	bool checkSatisfied();
	void printSolution(std::ostream& write_to);

	void addProducesLit( int strID, Var reach_var);



	P0LAcceptDetector(int _detectorID, LSystemSolver * _outer, LSystem &g_under, LSystem &g_over,
			 vec<vec<int>> &  strs, double seed = 1);
	virtual ~P0LAcceptDetector() {
		
	}
	
	const char* getName() {
		return "P0-LSystem Accepts Detector";
	}
	
	

};
}
;
#endif /* REACHDETECTOR_H_ */
