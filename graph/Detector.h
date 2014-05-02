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
namespace Minisat{

class Detector{
public:


	int detectorID;
	int unassigned_positives;
	int unassigned_negatives;
	int getID(){
		return detectorID;
	}
	virtual void printStats(){

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
protected:
	virtual void addLit(Lit l){
		unassigned_negatives++;
		unassigned_positives++;
	}
};
};


#endif /* DETECTOR_H_ */
