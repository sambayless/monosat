/*
 * Detector.h
 *
 *  Created on: 2014-01-08
 *      Author: sam
 */

#ifndef GEOMETRY_DETECTOR_H_
#define GEOMETRY_DETECTOR_H_
#include "core/SolverTypes.h"

#include "mtl/Vec.h"
using namespace Minisat;

class GeometryDetector{
public:
	int detectorID;

	int unassigned_positives=0;
	int unassigned_negatives=0;

	int getID(){
		return detectorID;
	}
	virtual void addPointContainmentLit(Lit l,vec<double> & point)=0;
	virtual bool propagate(vec<Lit> & conflict)=0;
	virtual void buildReason(Lit p, vec<Lit> & reason, CRef marker)=0;
	virtual bool checkSatisfied()=0;
	virtual void backtrackUntil(int level){

	}
	virtual void printStats(){

	}
	virtual void backtrackUntil(Lit p){

	}

	virtual void preprocess(){

	}
	virtual Lit decide()=0;
	//virtual vec<int> & getLitMap();
	GeometryDetector(int _detectorID):detectorID(_detectorID){};
	virtual ~GeometryDetector(){

	}
	virtual void addVar(Var v){
		unassigned_negatives++;
		unassigned_positives++;
	}
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
};



#endif /* GEOMETRY_DETECTOR_H_ */
