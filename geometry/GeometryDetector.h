/*
 * Detector.h
 *
 *  Created on: 2014-01-08
 *      Author: sam
 */

#ifndef GEOMETRY_DETECTOR_H_
#define GEOMETRY_DETECTOR_H_
#include "core/SolverTypes.h"
#include "GeometryTypes.h"
#include "mtl/Vec.h"
using namespace Minisat;

class GeometryDetector{
public:
	int detectorID;
	int getID(){
		return detectorID;
	}
	virtual void addPointContainmentLit(Lit l,vec<double> & point)=0;
	virtual bool propagate(vec<Lit> &trail ,vec<Lit> & conflict)=0;
	virtual void buildReason(Lit p, vec<Lit> & reason, CRef marker)=0;
	virtual bool checkSatisfied()=0;
	virtual void backtrackUntil(int level){

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
};



#endif /* GEOMETRY_DETECTOR_H_ */
