
#ifndef CONVEX_DETECTOR_H_
#define CONVEX_DETECTOR_H_
#include "core/SolverTypes.h"
#include "PointSet.h"
#include "GeometryDetector.h"
#include "GeometryTheory.h"
#include "ConvexHull.h"
#include "MonotoneConvexHull.h"
#include "Polygon.h"
using namespace Minisat;

template<unsigned int D, class T=double>
class ConvexHullDetector:public GeometryDetector{
public:
		GeometryTheorySolver * outer;
		//int within;
		PointSet<D,T>  over;
		PointSet<D,T>  under;
		ConvexHull<D,T>* over_hull;
		ConvexHull<D,T>* under_hull;

		double rnd_seed;
		CRef point_contained_marker;
		CRef point_not_contained_marker;
		CRef area_geq_marker;
		CRef area_not_geq_marker;
		vec<Lit> pointContainedLits;
		struct AreaLit{
			T areaGreaterEqThan;
			Lit l;
		};
		vec<AreaLit> areaDetectors;

		bool propagate(vec<Lit> & trail,vec<Lit> & conflict);
		void buildAreaGEQReason(T area, vec<Lit> & conflict);
		void buildAreaLTReason(T area,vec<Lit> & conflict);
		void buildPointContainedReason(vec<Lit> & conflict);
		void buildPointNotContainedReason(vec<Lit> & conflict);

		void buildReason(Lit p, vec<Lit> & reason, CRef marker);
		bool checkSatisfied();
		Lit decide();
		void addAreaDetectorLit(double areaGreaterEqThan, Var v);

		ConvexHullDetector(int _detectorID, GeometryTheorySolver * _outer,  double seed=1);
		virtual ~ConvexHullDetector(){

		}

};


template<unsigned int D, class T>
ConvexHullDetector<D,T>::ConvexHullDetector(int _detectorID, GeometryTheorySolver * _outer,double seed):
GeometryDetector(_detectorID),outer(_outer),rnd_seed(seed){

	point_contained_marker=outer->newReasonMarker(getID());
	point_not_contained_marker=outer->newReasonMarker(getID());

	over_hull = new MonotoneConvexHull<D,T>(over);
	under_hull = new MonotoneConvexHull<D,T>(under);
}

template<unsigned int D, class T>
void ConvexHullDetector<D,T>::addAreaDetectorLit(double areaGreaterEqThan, Var v){
	Lit l = mkLit(v,false);

	areaDetectors.push({areaGreaterEqThan,l});

}


template<unsigned int D, class T>
void ConvexHullDetector<D,T>::buildReason(Lit p, vec<Lit> & reason, CRef marker){

		if(marker==point_contained_marker){
			reason.push(p);

			//buildPointContainedReason(reason);

		}else if(marker==point_not_contained_marker){
			reason.push(p);

			//buildPointNotContainedReason(reason);

		}else{
			assert(false);
		}
}

template<unsigned int D, class T>
bool ConvexHullDetector<D,T>::propagate(vec<Lit> & trail,vec<Lit> & conflict){



		over_hull->update();

		under_hull->update();

		if(areaDetectors.size()){
			Polygon<D,T> & p = over_hull->getHull();
			double over_area = p.getArea();
			double under_area = under_hull->getHull().getArea();
			assert(under_area<=over_area);
			for(int i = 0;i<areaDetectors.size();i++){

				Lit l = areaDetectors[i].l;
				T areaGEQ = areaDetectors[i].areaGreaterEqThan;
				if(under_area>=areaGEQ){
					//l is true
					if(outer->S->value(l)==l_True){
						//do nothing
					}else if(outer->S->value(l)==l_Undef){

						outer->S->uncheckedEnqueue(l,area_geq_marker) ;
					}else if (outer->S->value(l)==l_False){
						conflict.push(l);
						//buildAreaGEQReason(under_area,conflict);
						return false;
					}
				}else if (over_area<areaGEQ){
					l=~l;
					//l is true
					if(outer->S->value(l)==l_True){
						//do nothing
					}else if(outer->S->value(l)==l_Undef){
						outer->S->uncheckedEnqueue(l,area_not_geq_marker) ;
					}else if (outer->S->value(l)==l_False){
						conflict.push(l);
						//buildAreaLTReason(over_area,conflict);
						return false;
					}
				}


			}
		}
			return true;
		}
template<unsigned int D, class T>
bool ConvexHullDetector<D,T>::checkSatisfied(){

	return true;
}
template<unsigned int D, class T>
Lit ConvexHullDetector<D,T>::decide(){

	return lit_Undef;
}



#endif
