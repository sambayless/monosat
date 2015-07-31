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

#ifndef DETECTOR_H_
#define DETECTOR_H_
#include "core/SolverTypes.h"
#include "GraphTheoryTypes.h"
#include "mtl/Vec.h"
#include "core/Config.h"
#include <cstdio>
#include <iostream>
namespace Monosat {

//Graph properties (such as shortest paths, minimum spanning tree weights) are computed by individual 'detectors',
//atached to the central graph theory. This allows them to share their edge atoms and a lot of other infrastructure.
class Detector {
public:
	
	int detectorID;

	int unassigned_positives = 0;
	int unassigned_negatives = 0;

	//Stats
	double stats_under_update_time = 0;
	double stats_over_update_time = 0;
	long stats_under_updates = 0;
	long stats_over_updates = 0;
	long stats_under_conflicts = 0;
	long stats_over_conflicts = 0;
	double stats_under_conflict_time = 0;
	double stats_over_conflict_time = 0;
	long stats_skipped_under_updates = 0;
	long stats_skipped_over_updates = 0;
	long stats_decisions = 0;
	double stats_decide_time = 0;
	long n_stats_priority_decisions =0;
	long n_stats_vsids_decisions=0;
	long stats_under_clause_length = 0;
	long stats_over_clause_length = 0;

	int getID() {
		return detectorID;
	}
	
	virtual const char* getName() {
		return "<unknown>";
	}
	
	virtual void printStats() {
		if (opt_verb > 0) {
			printf("Detector %d (%s):\n", getID(), getName());
			//printf("Updates: %d (under), %d over\n", stats_under_updates, stats_over_updates);
			printf("\tUnder-approx updates: %ld (%ld skipped) (%f s total, %f s avg)\n", stats_under_updates,
					stats_skipped_under_updates, (double) stats_under_update_time,
					(double) stats_under_update_time / (double) (stats_under_updates + 1));
			printf("\tOver-approx updates: %ld (%ld skipped)  (%f s total, %f s avg)\n", stats_over_updates,
					stats_skipped_over_updates, (double) stats_over_update_time,
					(double) stats_over_update_time / (double) (stats_over_updates + 1));
			printf("\tTheory Decisions: %ld (%f s total, %f s avg, %ld priority)\n", stats_decisions, (double) stats_decide_time,
					(double) stats_decide_time / (double) (stats_decisions + 1),n_stats_priority_decisions);
			printf(
					"\tConflicts (under,over): %ld (clause literals: %ld), %ld, (clause literals: %ld), (under time %f s, over time %f s)\n",
					stats_under_conflicts, stats_under_clause_length, stats_over_conflicts, stats_over_clause_length,
					stats_under_conflict_time, stats_over_conflict_time);
			
		}
	}
	
	virtual void printSolution(std::ostream & write_to = std::cout) {
	}
	
	virtual bool propagate(vec<Lit> & conflict)=0;
	virtual bool propagate(vec<Lit> & conflict	, bool backtrackOnly, Lit & conflictLit){
		return propagate(conflict);
	}

	virtual void buildReason(Lit p, vec<Lit> & reason, CRef marker)=0;
	virtual bool checkSatisfied()=0;
	virtual void preprocess() {
		
	}
	virtual void backtrack(int level) {
		//do nothing
	}
	virtual void undecideEdgeWeight(int edgeID){

	}
	virtual void undecide(Lit l){

	}
	virtual void buildModel(){

	}

	virtual Lit decide()=0;
	virtual bool supportsEdgeDecisions(){
		return false;
	}
	virtual void suggestDecision(Lit l){
		//do nothing
	}
	virtual void setOccurs(Lit l, bool occurs) {
		if (!occurs) {
			if (sign(l))
				unassigned_negatives--;
			else
				unassigned_positives--;
		} else {
			if (sign(l))
				unassigned_negatives++;
			else
				unassigned_positives++;
		}
		assert(unassigned_positives >= 0);
		assert(unassigned_negatives >= 0);
	}
	virtual void assign(Lit l) {
		if (sign(l))
			unassigned_negatives--;
		else
			unassigned_positives--;
		assert(unassigned_positives >= 0);
		assert(unassigned_negatives >= 0);
		
	}
	virtual void unassign(Lit l) {
		if (sign(l))
			unassigned_negatives++;
		else
			unassigned_positives++;
	}
	//virtual vec<int> & getLitMap();
	Detector(int detectorID) :
			detectorID(detectorID), unassigned_positives(0), unassigned_negatives(0) {
	}
	;
	virtual ~Detector() {
		
	}
	
	/*	virtual void addLit(Lit l){
	 unassigned_negatives++;
	 unassigned_positives++;
	 }*/
	virtual void addVar(Var v) {
		unassigned_negatives++;
		unassigned_positives++;
	}
};

enum class DetectorComparison{
	lt,leq,gt,geq,eq,ne
};

template<typename Weight>
class EdgeDecider{
public:
	virtual bool decideEdgeWeight(int edgeID, Weight & store, DetectorComparison & op)=0;
	virtual ~EdgeDecider() {

	}
};

/*class LevelDetector: public Detector {
	
	vec<int> level_trail;

public:
	LevelDetector(int detectorID) :
			Detector(detectorID) {
	}
	;

	virtual ~LevelDetector() {
	}
	
	virtual void backtrack(int level) {
		while (level_trail.size() && level < level_trail.last()) {
			level_trail.pop();
			localBacktrack();
		}
	}
	virtual void localBacktrack() {
		//do nothing
	}
	//local decision level within the detector
	int decisionLevel() {
		return level_trail.size();
	}
	
	void newDecisionLevel(int outer_level) {
		
		level_trail.push(outer_level);
	}
	
};*/

}
;

#endif /* DETECTOR_H_ */
