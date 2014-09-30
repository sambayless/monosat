/*
 * Detector.h
 *
 *  Created on: 2014-01-08
 *      Author: sam
 */

#ifndef DETECTOR_H_
#define DETECTOR_H_
#include "core/SolverTypes.h"
#include "GraphTheoryTypes.h"
#include "mtl/Vec.h"
#include "core/Config.h"
#include <cstdio>
namespace Monosat{

class Detector{
public:


	int detectorID;

	int unassigned_positives=0;
	int unassigned_negatives=0;

	//Stats
	double stats_under_update_time=0;
	double stats_over_update_time=0;
	long stats_under_updates=0;
	long stats_over_updates=0;
	long stats_under_conflicts=0;
	long stats_over_conflicts=0;
	long stats_under_conflict_time=0;
	long stats_over_conflict_time=0;
	long stats_skipped_under_updates=0;
	long stats_skipped_over_updates=0;

	int getID(){
		return detectorID;
	}

	virtual const char* getName(){
		return "<unknown>";
	}

	virtual void printStats(){
		if(opt_verb>0){
			printf("Detector %d (%s):\n", getID(), getName());
			//printf("Updates: %d (under), %d over\n", stats_under_updates, stats_over_updates);
			printf("\tUnder-approx updates: %d (%d skipped) (%f s total, %f s avg)\n", stats_under_updates,stats_skipped_under_updates,(double)stats_under_update_time, (double)stats_under_update_time/(double)(stats_under_updates+1) );
			printf("\tOver-approx updates: %d (%d skipped)  (%f s total, %f s avg)\n", stats_over_updates,stats_skipped_over_updates,(double)stats_over_update_time, (double)stats_over_update_time/(double)(stats_over_updates+1) );


		}
	}

	virtual void printSolution(){

	}

	virtual bool propagate(vec<Lit> & conflict)=0;
	/*virtual void buildReachReason(int node,vec<Lit> & conflict)=0;
	virtual void buildNonReachReason(int node,vec<Lit> & conflict)=0;
	virtual void buildForcedEdgeReason(int reach_node, int forced_edge_id,vec<Lit> & conflict)=0;*/
	virtual void buildReason(Lit p, vec<Lit> & reason, CRef marker)=0;
	virtual bool checkSatisfied()=0;
	virtual void preprocess(){

	}
	virtual Lit decide()=0;
	virtual void setOccurs(Lit l, bool occurs){
		if(!occurs){
			 if(sign(l))
			   unassigned_negatives--;
			 else
			   unassigned_positives--;
		}else{
			 if(sign(l))
			   unassigned_negatives++;
			 else
			   unassigned_positives++;
		}
		assert(unassigned_positives>=0);
		assert(unassigned_negatives>=0);
	}
	virtual void assign(Lit l){
		 if(sign(l))
		   unassigned_negatives--;
		 else
		   unassigned_positives--;
		 assert(unassigned_positives>=0);
		 assert(unassigned_negatives>=0);

	}
	virtual	 void unassign(Lit l){
		 if(sign(l))
		   unassigned_negatives++;
		 else
		   unassigned_positives++;
	}
	//virtual vec<int> & getLitMap();
	Detector(int _detectorID):detectorID(_detectorID),unassigned_positives(0),unassigned_negatives(0){};
	virtual ~Detector(){

	}

/*	virtual void addLit(Lit l){
		unassigned_negatives++;
		unassigned_positives++;
	}*/
	virtual void addVar(Var v){
		unassigned_negatives++;
		unassigned_positives++;
	}
};
};


#endif /* DETECTOR_H_ */
