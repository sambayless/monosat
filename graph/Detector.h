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
struct Change{
		Lit l;
		int u;
	};
class Detector{
public:
	CRef reach_marker;
	CRef non_reach_marker;
	CRef forced_reach_marker;
	int detectorID;
	int getID(){
		return detectorID;
	}
	virtual bool propagate(vec<Assignment> &trail ,vec<Lit> & conflict)=0;
	virtual void buildReachReason(int node,vec<Lit> & conflict)=0;
	virtual void buildNonReachReason(int node,vec<Lit> & conflict)=0;
	virtual void buildForcedEdgeReason(int reach_node, int forced_edge_id,vec<Lit> & conflict)=0;
	virtual void buildReason(Lit p, vec<Lit> & reason, CRef marker)=0;

	virtual	vec<Change> & getChanged()=0;
	virtual Lit decide()=0;
	Detector(int _detectorID):detectorID(_detectorID){};
	virtual ~Detector(){

	}
};
};


#endif /* DETECTOR_H_ */
