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
#ifndef CONVEX_DETECTOR_H_
#define CONVEX_DETECTOR_H_
#include "monosat/core/SolverTypes.h"

#include "GeometryDetector.h"
#include "GeometryTheory.h"
#include "Heightmap.h"

#include "Polygon.h"
using namespace Monosat;

template<unsigned int D, class T = double>
class HeightmapDetector: public GeometryDetector {
public:
	GeometryTheorySolver<D, T> * outer;
	//int within;
	Heightmap<D, T>* over_hull;
	Heightmap<D, T>* under_hull;

	double rnd_seed;
	CRef point_contained_marker;
	CRef point_not_contained_marker;
	CRef area_geq_marker;
	CRef area_not_geq_marker;

	Var lowest_point_var;
	vec<int> point_lit_map;
	vec<Lit> lit_point_map;

	struct PointContainedLit {
		Point<D, T> p;
		Lit l;
	};
	vec<PointContainedLit> pointContainedLits;
	struct AreaLit {
		T areaGreaterEqThan;
		Lit l;
	};
	vec<AreaLit> areaDetectors;

	int qhead;

	bool propagate(vec<Lit> & trail, vec<Lit> & conflict);
	void buildAreaGEQReason(T area, vec<Lit> & conflict);
	void buildAreaLTReason(T area, vec<Lit> & conflict);
	void buildPointContainedReason(vec<Lit> & conflict);
	void buildPointNotContainedReason(vec<Lit> & conflict);
	/*
	 void backtrackUntil(int level){
	 assert(qhead<outer->S->trail.size());
	 for (int c = qhead; c >= outer->S->trail_lim[level]; c--){
	 Lit l = outer->S->trail[c];
	 Var      x  = var(outer->S->trail[c]);
	 int pointIndex = x-lowest_point_var;
	 if(pointIndex>=0 && pointIndex<point_lit_map.size() && point_lit_map[pointIndex]!=-1){
	 //then this is an assignment to a hull point;
	 int pointID = point_lit_map[pointIndex];
	 if(sign(l)){
	 //this point is excluded from the hull
	 assert(!under.isEnabled(pointID));
	 over.setPointEnabled(pointID,true);
	 }else{
	 assert(over.isEnabled(pointID));
	 under.setPointEnabled(pointID,false);
	 }
	 }
	 }
	 qhead = outer->S->trail_lim[level];
	 }
	 virtual void backtrackUntil(Lit p){
	 assert(qhead<outer->S->trail.size());
	 for ( ; qhead >= 0; qhead--){
	 Lit l = outer->S->trail[qhead];
	 if(l==p){
	 break;
	 }
	 Var      x  = var(l);
	 int pointIndex = x-lowest_point_var;
	 if(pointIndex>=0 && pointIndex<point_lit_map.size() && point_lit_map[pointIndex]!=-1){
	 //then this is an assignment to a hull point;
	 int pointID = point_lit_map[pointIndex];
	 if(sign(l)){
	 //this point is excluded from the hull
	 assert(!under.isEnabled(pointID));
	 over.setPointEnabled(pointID,true);
	 }else{
	 assert(over.isEnabled(pointID));
	 under.setPointEnabled(pointID,false);
	 }
	 }
	 }

	 }*/
	void buildReason(Lit p, vec<Lit> & reason, CRef marker);
	bool checkSatisfied();
	Lit decide();
	void addAreaDetectorLit(double areaGreaterEqThan, Var v);
	void addPoint(const Point<D, T> & point, Lit l) {
		int id = over.addPoint(point);
		int underID = under.addPoint(point);
		assert(underID == id);
		if (lowest_point_var == var_Undef) {
			lowest_point_var = var(l);
		}
		assert(var(l) >= lowest_point_var);
		int index = var(l) - lowest_point_var;
		point_lit_map.growTo(index + 1, -1);
		point_lit_map[index] = id;
		lit_point_map.growTo(id + 1, lit_Undef);
		lit_point_map[id] = l;
	}
	
	void addPointContainmentLit(Lit l, vec<double> & point) {
		assert(point.size() == D);
		Point<D, T> p(point);
		pointContainedLits.push( { p, l });
	}
	
	HeightmapDetector(int _detectorID, GeometryTheorySolver * _outer, double seed = 1);
	virtual ~HeightmapDetector() {
		
	}
	
};

template<unsigned int D, class T>
HeightmapDetector<D, T>::HeightmapDetector(int _detectorID, GeometryTheorySolver * _outer, double seed) :
		GeometryDetector(_detectorID), outer(_outer), rnd_seed(seed) {
	
	point_contained_marker = outer->newReasonMarker(getID());
	point_not_contained_marker = outer->newReasonMarker(getID());
	if (hullAlg == ALG_QUICKHULL) {
		if (D == 2) {
			over_hull = new QuickConvexHull<D, T>(over);
			under_hull = new QuickConvexHull<D, T>(under);
		} else if (D == 3) {
			over_hull = new QuickConvexHull<D, T>(over); //new MonotoneConvexHull<D,T>(over);
			under_hull = new QuickConvexHull<D, T>(under);
		}
	} else if (hullAlg == ALG_MONOTONE_HULL) {
		if (D == 2) {
			over_hull = new MonotoneConvexHull<D, T>(over);
			under_hull = new MonotoneConvexHull<D, T>(under);
		}
	}
	lowest_point_var = var_Undef;
	qhead = 0;
}

template<unsigned int D, class T>
void HeightmapDetector<D, T>::addAreaDetectorLit(double areaGreaterEqThan, Var v) {
	Lit l = mkLit(v, false);
	
	areaDetectors.push( { areaGreaterEqThan, l });
	
}

template<unsigned int D, class T>
void HeightmapDetector<D, T>::buildReason(Lit p, vec<Lit> & reason, CRef marker) {
	
	if (marker == point_contained_marker) {
		reason.push(p);
		
		//buildPointContainedReason(reason);
		
	} else if (marker == point_not_contained_marker) {
		reason.push(p);
		
		//buildPointNotContainedReason(reason);
		
	} else {
		assert(false);
	}
}

template<unsigned int D, class T>
bool HeightmapDetector<D, T>::propagate(vec<Lit> & trail, vec<Lit> & conflict) {
	bool any_changed = false;
	for (; qhead < trail.size(); qhead++) {
		Lit l = trail[qhead];
		int pointIndex = var(l) - lowest_point_var;
		if (pointIndex >= 0 && pointIndex < point_lit_map.size() && point_lit_map[pointIndex] != -1) {
			any_changed = true;
			//then this is an assignment to a hull point;
			int pointID = point_lit_map[pointIndex];
			if (sign(l)) {
				//this point is excluded from the hull
				assert(!under.isEnabled(pointID));
				over.setPointEnabled(pointID, false);
			} else {
				assert(over.isEnabled(pointID));
				under.setPointEnabled(pointID, true);
			}
		}
	}
	
	if (any_changed) {
		
		over_hull->update();
		
		under_hull->update();
		
		if (areaDetectors.size()) {
			
			double over_area = over_hull->getHull().getArea();
			double under_area = under_hull->getHull().getArea();
			assert(under_area <= over_area);
			for (int i = 0; i < areaDetectors.size(); i++) {
				
				Lit l = areaDetectors[i].l;
				T areaGEQ = areaDetectors[i].areaGreaterEqThan;
				if (under_area >= areaGEQ) {
					//l is true
					if (outer->S->value(l) == l_True) {
						//do nothing
					} else if (outer->S->value(l) == l_Undef) {
						
						outer->S->uncheckedEnqueue(l, area_geq_marker);
					} else if (outer->S->value(l) == l_False) {
						conflict.push(l);
						//buildAreaGEQReason(under_area,conflict);
						return false;
					}
				} else if (over_area < areaGEQ) {
					l = ~l;
					//l is true
					if (outer->S->value(l) == l_True) {
						//do nothing
					} else if (outer->S->value(l) == l_Undef) {
						outer->S->uncheckedEnqueue(l, area_not_geq_marker);
					} else if (outer->S->value(l) == l_False) {
						conflict.push(l);
						//buildAreaLTReason(over_area,conflict);
						return false;
					}
				}
				
			}
		}
		//If we are making many queries, it is probably worth it to pre-process the polygon and then make the queries.
		for (int i = 0; i < pointContainedLits.size(); i++) {
			Point<D, T> point = pointContainedLits[i].p;
			Lit l = pointContainedLits[i].l;
			ConvexPolygon<D, T> & p_over = over_hull->getHull();
			ConvexPolygon<D, T> & p_under = under_hull->getHull();
			
			if (p_under.contains(point)) {
				//l is true
				if (outer->S->value(l) == l_True) {
					//do nothing
				} else if (outer->S->value(l) == l_Undef) {
					
					outer->S->uncheckedEnqueue(l, point_contained_marker);
				} else if (outer->S->value(l) == l_False) {
					conflict.push(l);
					//buildAreaGEQReason(under_area,conflict);
					return false;
				}
			} else if (!p_over.contains(point)) {
				l = ~l;
				//l is true
				if (outer->S->value(l) == l_True) {
					//do nothing
				} else if (outer->S->value(l) == l_Undef) {
					outer->S->uncheckedEnqueue(l, point_not_contained_marker);
				} else if (outer->S->value(l) == l_False) {
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
bool HeightmapDetector<D, T>::checkSatisfied() {
	
	return true;
}
template<unsigned int D, class T>
Lit HeightmapDetector<D, T>::decide() {
	
	return lit_Undef;
}

#endif
