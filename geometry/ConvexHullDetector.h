
#ifndef CONVEX_DETECTOR_H_
#define CONVEX_DETECTOR_H_
#include "core/SolverTypes.h"
#include "PointSet.h"
#include "GeometryTheory.h"
#include "GeometryDetector.h"

#include "ConvexHull.h"
#include "MonotoneConvexHull.h"
#include "QuickConvexHull.h"
#include "Polygon.h"
using namespace Minisat;

template<unsigned int D, class T>
class ConvexHullDetector:public GeometryDetector{
public:
		GeometryTheorySolver<D,T> * outer;
		//int within;
		PointSet<D,T> & over;
		PointSet<D,T> & under;
		ConvexHull<D,T>* over_hull;
		ConvexHull<D,T>* under_hull;

		double rnd_seed;
		CRef point_contained_marker;
		CRef point_not_contained_marker;
		CRef area_geq_marker;
		CRef area_not_geq_marker;

		Var lowest_point_var;
		vec<int> point_lit_map;
		vec<Lit> lit_point_map;

		struct PointContainedLit{
			Point<D,T> p;
			Lit l;
		};
		vec<PointContainedLit> pointContainedLits;
		struct AreaLit{
			T areaGreaterEqThan;
			Lit l;
		};
		vec<AreaLit> areaDetectors;

		int qhead;


		bool propagate(vec<Lit> & conflict);
		void buildAreaGEQReason(T area, vec<Lit> & conflict);
		void buildAreaLTReason(T area,vec<Lit> & conflict);
		void buildPointContainedReason(vec<Lit> & conflict);
		void buildPointNotContainedReason(vec<Lit> & conflict);


		void buildReason(Lit p, vec<Lit> & reason, CRef marker);
		bool checkSatisfied();
		Lit decide();
		void addAreaDetectorLit(double areaGreaterEqThan, Var v);
		void addPointContainmentLit(Point<D,T> p,Var outerVar);

		ConvexHullDetector(int detectorID,PointSet<D,T> & under, PointSet<D,T> & over, GeometryTheorySolver<D,T> * outer,  double seed=1);
		virtual ~ConvexHullDetector(){

		}

};


template<unsigned int D, class T>
ConvexHullDetector<D,T>::ConvexHullDetector(int detectorID,PointSet<D,T> & under, PointSet<D,T> & over, GeometryTheorySolver<D,T> * outer,double seed):
GeometryDetector(detectorID),under(under),over(over),outer(outer),rnd_seed(seed){

	point_contained_marker=outer->newReasonMarker(getID());
	point_not_contained_marker=outer->newReasonMarker(getID());
	if(hullAlg== ConvexHullAlg::ALG_QUICKHULL){
		if(D==2){
			over_hull = new QuickConvexHull<D,T>(over);
			under_hull = new QuickConvexHull<D,T>(under);
		}else if (D==3){
			over_hull = new QuickConvexHull<D,T>(over); //new MonotoneConvexHull<D,T>(over);
			under_hull = new QuickConvexHull<D,T>(under);
		}
	}else if (hullAlg== ConvexHullAlg::ALG_MONOTONE_HULL){
		if(D==2){
			over_hull = new MonotoneConvexHull<D,T>(over);
			under_hull = new MonotoneConvexHull<D,T>(under);
		}
	}
	lowest_point_var=var_Undef;
	qhead=0;
}

template<unsigned int D, class T>
void ConvexHullDetector<D,T>::addPointContainmentLit(Point<D,T> p,Var outerVar){

		under.invalidate();
		over.invalidate();

		Var containVar = outer->newVar(outerVar,getID());
		pointContainedLits.push({p,mkLit(containVar,false)});
	}

template<unsigned int D, class T>
void ConvexHullDetector<D,T>::addAreaDetectorLit(double areaGreaterEqThan, Var outerVar){
	Var v = outer->newVar(outerVar,getID());
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
bool ConvexHullDetector<D,T>::propagate(vec<Lit> & conflict){



		over_hull->update();

		under_hull->update();

		if(areaDetectors.size()){

			double over_area = over_hull->getHull().getArea();
			double under_area = under_hull->getHull().getArea();
			assert(under_area<=over_area);
			for(int i = 0;i<areaDetectors.size();i++){

				Lit l = areaDetectors[i].l;
				T areaGEQ = areaDetectors[i].areaGreaterEqThan;
				if(under_area>=areaGEQ){
					//l is true
					if(outer->value(l)==l_True){
						//do nothing
					}else if(outer->value(l)==l_Undef){

						outer->enqueue(l,area_geq_marker) ;
					}else if (outer->value(l)==l_False){
						conflict.push(l);
						//buildAreaGEQReason(under_area,conflict);
						return false;
					}
				}else if (over_area<areaGEQ){
					l=~l;
					//l is true
					if(outer->value(l)==l_True){
						//do nothing
					}else if(outer->value(l)==l_Undef){
						outer->enqueue(l,area_not_geq_marker) ;
					}else if (outer->value(l)==l_False){
						conflict.push(l);
						//buildAreaLTReason(over_area,conflict);
						return false;
					}
				}


			}
		}
		//If we are making many queries, it is probably worth it to pre-process the polygon and then make the queries.
		for(int i =0;i<pointContainedLits.size();i++){
			Point<D,T> point = pointContainedLits[i].p;
			Lit l = pointContainedLits[i].l;
			ConvexPolygon<D,T> & p_over = over_hull->getHull();
			ConvexPolygon<D,T> & p_under = under_hull->getHull();

			if(p_under.contains(point)){
				//l is true
				if(outer->value(l)==l_True){
					//do nothing
				}else if(outer->value(l)==l_Undef){

					outer->enqueue(l,point_contained_marker) ;
				}else if (outer->value(l)==l_False){
					conflict.push(l);
					//buildAreaGEQReason(under_area,conflict);
					return false;
				}
			}else if (!p_over.contains(point)){
				l=~l;
				//l is true
				if(outer->value(l)==l_True){
					//do nothing
				}else if(outer->value(l)==l_Undef){
					outer->enqueue(l,point_not_contained_marker) ;
				}else if (outer->value(l)==l_False){
					conflict.push(l);
					//buildAreaLTReason(over_area,conflict);
					return false;
				}
			}
		}

		return true;
	}
template<unsigned int D, class T>
bool ConvexHullDetector<D,T>::checkSatisfied(){


	/*for (auto & p:over.getEnabledPoints()){

	}*/
	MonotoneConvexHull<D,T> cv(over);
	T area = cv.getArea();
	for(auto & a: areaDetectors){
		T area_cmp = a.areaGreaterEqThan;
		Lit l = a.l;
		if(outer->value(l)==l_True){
			if(area_cmp>=area){
				return false;
			}
		}else if(outer->value(l)==l_False){
			if(area_cmp<area){
				return false;
			}
		}
	}

	return true;
}
template<unsigned int D, class T>
Lit ConvexHullDetector<D,T>::decide(){

	return lit_Undef;
}

template<unsigned int D, class T>
void ConvexHullDetector<D,T>::buildAreaGEQReason(T area, vec<Lit> & conflict){

}
template<unsigned int D, class T>
void ConvexHullDetector<D,T>::buildAreaLTReason(T area,vec<Lit> & conflict){

}
template<unsigned int D, class T>
void ConvexHullDetector<D,T>::buildPointContainedReason(vec<Lit> & conflict){

}
template<unsigned int D, class T>
void ConvexHullDetector<D,T>::buildPointNotContainedReason(vec<Lit> & conflict){

}

#endif
