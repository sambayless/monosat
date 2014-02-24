/*
 * GeometryTheory.cpp
 *
 *  Created on: 2014-02-23
 *      Author: sam
 */

#include "GeometryTheory.h"
#include "ConvexHullDetector.h"
using namespace Minisat;
void GeometryTheorySolver::createConvexHull(int id, const int D){

	GeometryDetector * d;
	int detectorID = detectors.size();

	//should probably swap this out for some fancy template recursion...
	switch (D){
		case 2:
			d = new ConvexHullDetector<2,double>(detectorID,this,drand(rnd_seed));break;
		case 3:
			//d = new ConvexHullDetector<3,double>(detectorID,this,drand(rnd_seed));break;
		default:
			assert(false);
			fprintf(stderr,"Unsupported number of dimensions (%d), exiting\n",D);fflush(stderr);exit(1);
	};
	detectors.push(d);
	polygonList.growTo(id+1);
	polygonList[id]=d;
}

void GeometryTheorySolver::addHullPoint(int hullID,Lit l, vec<double> & point){
	assert(polygonList.size()>hullID);
	assert(polygonList[hullID]!=NULL);
	GeometryDetector*d = polygonList[hullID];
	if(point.size()==2){
		 ConvexHullDetector<2,double>* hull = (ConvexHullDetector<2,double>*) d;

		 hull->addPoint(Point2D(point),l);
	}else{
		fprintf(stderr,"Unsupported number of dimensions (%d), exiting\n",point.size());fflush(stderr);exit(1);
	}

}
