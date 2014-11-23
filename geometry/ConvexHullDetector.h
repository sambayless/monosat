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
#include "GeometryTypes.h"
#include "core/SolverTypes.h"
#include "PointSet.h"
#include "GeometryDetector.h"
#include "core/Config.h"
#include "ConvexHull.h"
#include "MonotoneConvexHull.h"
#include "QuickConvexHull.h"
#include "Polygon.h"
#include "LineSegment.h"
#include "Line.h"
#include <vector>
using namespace Monosat;

template<unsigned int D, class T>
class ConvexHullDetector: public GeometryDetector {
	bool hasSet = false;
public:
	GeometryTheorySolver<D, T> * outer;
	//int within;
	PointSet<D, T> & under;
	PointSet<D, T> & over;
	ConvexHull<D, T>* under_hull = nullptr;
	ConvexHull<D, T>* over_hull = nullptr;

	double rnd_seed;
	CRef point_contained_marker;
	CRef point_not_contained_marker;
	CRef convex_intersection_marker;
	CRef convex_not_intersection_marker;
	CRef line_intersection_marker;
	CRef line_not_intersection_marker;
	CRef area_geq_marker;
	CRef area_not_geq_marker;
	CRef point_on_hull_marker;
	CRef point_not_on_hull_marker;

	Var lowest_point_var;
	vec<int> point_lit_map;
	vec<Lit> lit_point_map;
	vec<bool> under_hull_members;
	vec<bool> over_hull_members;

	NConvexPolygon<D, T> tmp_polygon;

	struct PointContainedLit {
		Point<D, T> p;
		Lit l;
		bool inclusive;
		//indices of 3 containing vertices in the under approximation that contain this point (or -1 if no such triangle exists).
		NConvexPolygon<D, T> under_containing_triangle;
		NConvexPolygon<D, T> over_containing_triangle;
		
	};
	std::vector<PointContainedLit> pointContainedLits;

	struct LineIntersectionLit {
		LineSegment<D, T> line;
		Lit l;
		bool inclusive;
		NConvexPolygon<D, T> under_intersecting_polygon;
		NConvexPolygon<D, T> over_intersecting_polygon;
		
	};
	std::vector<LineIntersectionLit> lineIntersectionLits;

	struct ConvexIntersectionLit {
		NConvexPolygon<D, T> polygon;
		Lit l;
		bool inclusive;
	};
	std::vector<ConvexIntersectionLit> convexIntersectionLits;

	struct PointOnHullLit {
		Var pointVar;
		int pointIndex;
		Point<D, T> p;
		Lit l;
	};
	std::vector<PointOnHullLit> pointOnHullLits;

	struct AreaLit {
		T areaGreaterEqThan;
		Lit l;
	};
	std::vector<AreaLit> areaDetectors;

	static long stats_propagations;
	static long stats_bound_checks;
	static long stats_bounds_skips_under;
	static long stats_bounds_skips_over;
	static long stats_under_clause_length;
	static long stats_over_clause_length;
	static long stats_under_clauses;
	static long stats_over_clauses;
	static long stats_line_intersection_skips_under;
	static long stats_line_intersection_skips_over;

	void printStats() {
		printf("Convex hull %d: ", getID());
		cout << under_hull->getHull() << "\n";
		
		printf("propagations: %ld\n", stats_propagations);
		printf("bound checks %ld (skipped under,over: %ld,%ld)\n", stats_bound_checks, stats_bounds_skips_under,
				stats_bounds_skips_over);
		printf("containment: bounds skipped %ld, triangle skipped %ld, checks %ld, %ld depth\n",
				ConvexPolygon<D, T>::stats_bounds_avoided, ConvexPolygon<D, T>::stats_triangle_avoided,
				ConvexPolygon<D, T>::stats_split_full_checks, ConvexPolygon<D, T>::stats_split_checks_depths);
		printf("intersections: bounds skipped %ld, checks skipped under: %ld, skipped over: %ld\n",
				ConvexPolygon<D, T>::stats_bounds_intersections_avoided, stats_line_intersection_skips_under,
				stats_line_intersection_skips_over);
		
		printf("Skipped updates: %ld/%ld under, %ld/%ld over\n",
				((MonotoneConvexHull<D, T>*) under_hull)->stats_skipped_updates,
				((MonotoneConvexHull<D, T>*) under_hull)->stats_updates,
				((MonotoneConvexHull<D, T>*) over_hull)->stats_skipped_updates,
				((MonotoneConvexHull<D, T>*) under_hull)->stats_updates);
	}
	void printSolution() {
		printf("Convex hull %d: ", getID());
		cout << under_hull->getHull() << "\n";
	}
	bool propagate(vec<Lit> & conflict);
	void buildAreaGEQReason(T area, vec<Lit> & conflict);
	void buildAreaLTReason(T area, vec<Lit> & conflict);
	void buildPointContainedReason(const Point<D, T> & p, vec<Lit> & conflict, bool inclusive);
	void buildPointNotContainedReason(const Point<D, T> & p, vec<Lit> & conflict, bool inclusive);
	void buildPointOnHullOrDisabledReason(Var pointVar, const Point<D, T> & p, vec<Lit> & conflict);
	void buildPointNotOnHullOrDisabledReason(Var pointVar, const Point<D, T> & p, vec<Lit> & conflict);
	void buildConvexIntersectsReason(ConvexPolygon<D, T> & polygon, vec<Lit> & conflict, bool inclusive);
	void buildConvexNotIntersectsReason(ConvexPolygon<D, T> & polygon, vec<Lit> & conflict, bool inclusive);

	void buildLineIntersectsReason(LineSegment<D, T> & polygon, vec<Lit> & conflict, bool inclusive);
	void buildLineNotIntersectsReason(LineSegment<D, T> & polygon, vec<Lit> & conflict, bool inclusive);

	void buildReason(Lit p, vec<Lit> & reason, CRef marker);
	bool checkSatisfied();
	Lit decide();
	void addAreaDetectorLit(T areaGreaterEqThan, Var v);
	void addPointContainmentLit(Point<D, T> p, Var outerVarr, bool inclusive);
	void addPointOnHullLit(int pointsetIndex, Var outerVar);
	//void addLineIntersection(LineSegment<D,T> line, Var outerVar);
	void addConvexIntersectionLit(NConvexPolygon<D, T> &polygon, Var outerVarr, bool inclusive);
	void addLineIntersectionLit(LineSegment<D, T> line, Var outerVar, bool inclusive);

	ConvexHullDetector(int detectorID, PointSet<D, T> & under, PointSet<D, T> & over,
			GeometryTheorySolver<D, T> * outer, double seed = 1);
	~ConvexHullDetector() {
		if (under_hull)
			delete (under_hull);
		if (over_hull)
			delete (over_hull);
	}
	
	ConvexHull<D, T>* getConvexHull(bool overApprox) {
		if (overApprox) {
			return over_hull;
		} else {
			return under_hull;
		}
	}
	ConvexHull<D, T>* getPointset(bool overApprox) {
		if (overApprox) {
			return over;
		} else {
			return under;
		}
	}
	
private:
	void findFarPoints(Line<D, T> & testline, std::vector<Point<D, T> > & test_set, std::vector<Point<D, T> > & min_set,
			ConvexPolygon<D, T> &hull, bool includeCollinearPoints) {
		test_set.clear();
		T hullside = 0;
		for (int i = 0; hullside == 0 && i < hull.size(); i++) {
			hullside = testline.whichSide(hull[i]);
		}
		assert(hullside != 0);
		
		for (int i = 0; i < over.size(); i++) {
			if (!over.pointEnabled(i)) {
				int side = testline.whichSide(over[i]);
				if (side != hullside && (includeCollinearPoints || side != 0)) {
					test_set.push_back(over[i]);
					if (hasSet && test_set.size() >= min_set.size()) {
						return; //shortcut
					}
				}
			}
		}
		
		if (test_set.size() < min_set.size() || !hasSet) {
			hasSet = true;
			min_set = test_set;
			//test_set.copyTo(min_set);
		}
	}
	void findPointsSameSide(Line<D, T> & testline, bool includeCollinearPoints, Point<D, T> & point,
			std::vector<Point<D, T> > & points_out, ConvexPolygon<D, T> &hull, bool inclusive) {
		points_out.clear();
		T hullside = 0;
		hullside = testline.whichSide(point);
		assert(hullside != 0);
		for (int i = 0; i < over.size(); i++) {
			if (!over.pointEnabled(i)) {
				int side = testline.whichSide(over[i]);
				if (side == hullside && (includeCollinearPoints || side != 0)) {
					points_out.push_back(over[i]);
				}
			}
		}
	}
	
	bool findSeparatingAxis2d(ConvexPolygon<2, T> & hull1, ConvexPolygon<2, T> & hull2, PointSet<2, T> & pointset1,
			std::vector<std::pair<Point<2, T>, T>> &projection_out,
			std::vector<std::pair<Point<2, T>, T>> &projection_out2, bool inclusive);
	void buildConvexIntersectsReason2d(ConvexPolygon<2, T> & line, vec<Lit> & conflict, bool inclusive);
	void buildConvexNotIntersectsReason2d(ConvexPolygon<2, T> & line, vec<Lit> & conflict, bool inclusive);
	void buildPointContainedReason2d(const Point<2, T> & p, vec<Lit> & conflict, bool inclusive);
	void buildPointNotContainedReason2d(const Point<2, T> & p, vec<Lit> & conflict, bool inclusive);
	void buildPointOnHullOrDisabledReason2d(Var pointVar, const Point<2, T> & p, vec<Lit> & conflict);
	void buildPointNotOnHullOrDisabledReason2d(Var pointVar, const Point<2, T> & p, vec<Lit> & conflict);
	inline bool checkContainingTriangle(NConvexPolygon<D, T> & containing, const Point<2, T> & point,
			ConvexPolygon<D, T> & hull, PointSet<D, T> & pointset, bool inclusive) {
		if (containing.size() == 0)
			return false;
		for (auto &p : containing) {
			int index = outer->getPointsetIndex(p.getID());
			if (!pointset.pointEnabled(index)) {
				containing.clear();
				return false;
			}
		}
		assert(containing.contains(point, inclusive));
		assert(hull.contains(point, inclusive));
		if (&pointset == &under) {
			stats_bounds_skips_under++;
		} else if (&pointset == &over) {
			stats_bounds_skips_over++;
		}
		return true;
		
	}
	inline bool checkLineIntersection(NConvexPolygon<D, T> & containing, LineSegment<D, T> & line,
			ConvexPolygon<D, T> & hull, PointSet<D, T> & pointset, bool inclusive) {
		if (containing.size() == 0)
			return false;
		for (auto &p : containing) {
			int index = outer->getPointsetIndex(p.getID());
			if (!pointset.pointEnabled(index)) {
				containing.clear();
				return false;
			}
		}
		assert(containing.intersects(line, inclusive));
		assert(hull.intersects(line, inclusive));
		if (&pointset == &under) {
			stats_line_intersection_skips_under++;
		} else if (&pointset == &over) {
			stats_line_intersection_skips_over++;
		}
		return true;
		
	}
};
template<unsigned int D, class T>
long ConvexHullDetector<D, T>::stats_bounds_skips_under = 0;
template<unsigned int D, class T>
long ConvexHullDetector<D, T>::stats_bounds_skips_over = 0;

template<unsigned int D, class T>
long ConvexHullDetector<D, T>::stats_line_intersection_skips_under = 0;
template<unsigned int D, class T>
long ConvexHullDetector<D, T>::stats_line_intersection_skips_over = 0;

template<unsigned int D, class T>
long ConvexHullDetector<D, T>::stats_propagations = 0;
template<unsigned int D, class T>
long ConvexHullDetector<D, T>::stats_bound_checks = 0;

template<unsigned int D, class T>
long ConvexHullDetector<D, T>::stats_under_clause_length = 0;
template<unsigned int D, class T>
long ConvexHullDetector<D, T>::stats_over_clause_length = 0;
template<unsigned int D, class T>
long ConvexHullDetector<D, T>::stats_under_clauses = 0;
template<unsigned int D, class T>
long ConvexHullDetector<D, T>::stats_over_clauses = 0;

template<unsigned int D, class T>
ConvexHullDetector<D, T>::ConvexHullDetector(int detectorID, PointSet<D, T> & under, PointSet<D, T> & over,
		GeometryTheorySolver<D, T> * outer, double seed) :
		GeometryDetector(detectorID), outer(outer), under(under), over(over), rnd_seed(seed) {
	
	point_contained_marker = outer->newReasonMarker(getID());
	point_not_contained_marker = outer->newReasonMarker(getID());
	point_on_hull_marker = outer->newReasonMarker(getID());
	point_not_on_hull_marker = outer->newReasonMarker(getID());
	area_geq_marker = outer->newReasonMarker(getID());
	area_not_geq_marker = outer->newReasonMarker(getID());
	convex_intersection_marker = outer->newReasonMarker(getID());
	convex_not_intersection_marker = outer->newReasonMarker(getID());
	line_intersection_marker = outer->newReasonMarker(getID());
	line_not_intersection_marker = outer->newReasonMarker(getID());
	
	if (hullAlg == ConvexHullAlg::ALG_QUICKHULL) {
		over_hull = new QuickConvexHull<D, T>(over);
		under_hull = new QuickConvexHull<D, T>(under);
	} else if (hullAlg == ConvexHullAlg::ALG_MONOTONE_HULL) {
		
		over_hull = new MonotoneConvexHull<D, T>(over);
		under_hull = new MonotoneConvexHull<D, T>(under);
		
	}
	lowest_point_var = var_Undef;
	
}

template<unsigned int D, class T>
void ConvexHullDetector<D, T>::addPointContainmentLit(Point<D, T> p, Var outerVar, bool inclusive) {
	
	under.invalidate();
	over.invalidate();
	
	Var containVar = outer->newVar(outerVar, getID());
	pointContainedLits.push_back( { p, mkLit(containVar, false), inclusive });
}

template<unsigned int D, class T>
void ConvexHullDetector<D, T>::addLineIntersectionLit(LineSegment<D, T> line, Var outerVar, bool inclusive) {
	
	under.invalidate();
	over.invalidate();
	
	Var containVar = outer->newVar(outerVar, getID());
	lineIntersectionLits.push_back( { line, mkLit(containVar, false), inclusive });
}

template<unsigned int D, class T>
void ConvexHullDetector<D, T>::addConvexIntersectionLit(NConvexPolygon<D, T> & polygon, Var outerVar, bool inclusive) {
	
	under.invalidate();
	over.invalidate();
	
	Var containVar = outer->newVar(outerVar, getID());
	convexIntersectionLits.push_back( { polygon, mkLit(containVar, false), inclusive });
}

template<unsigned int D, class T>
void ConvexHullDetector<D, T>::addAreaDetectorLit(T areaGreaterEqThan, Var outerVar) {
	Var v = outer->newVar(outerVar, getID());
	Lit l = mkLit(v, false);
	
	areaDetectors.push_back( { areaGreaterEqThan, l });
	
}
template<unsigned int D, class T>
void ConvexHullDetector<D, T>::addPointOnHullLit(int pointsetIndex, Var outerVar) {
	Var v = outer->newVar(outerVar, getID());
	Lit l = mkLit(v, false);
	Point<D, T> & point = under[pointsetIndex];
	Var pointsetVar = outer->getPointVar(point.getID());
	pointOnHullLits.push_back( { pointsetVar, pointsetIndex, point, v });
}

template<unsigned int D, class T>
void ConvexHullDetector<D, T>::buildReason(Lit p, vec<Lit> & reason, CRef marker) {
	reason.push(p);
	if (marker == area_geq_marker) {
		
		for (auto & a : areaDetectors) {
			if (var(a.l) == var(p)) {
				buildAreaGEQReason(a.areaGreaterEqThan, reason);
				return;
			}
		}
		
	} else if (marker == area_not_geq_marker) {
		for (auto & a : areaDetectors) {
			if (var(a.l) == var(p)) {
				buildAreaLTReason(a.areaGreaterEqThan, reason);
				return;
			}
		}
		
	} else if (marker == point_contained_marker) {
		for (auto & a : pointContainedLits) {
			if (var(a.l) == var(p)) {
				buildPointContainedReason(a.p, reason, a.inclusive);
				return;
			}
		}
		
	} else if (marker == point_not_contained_marker) {
		for (auto & a : pointContainedLits) {
			if (var(a.l) == var(p)) {
				buildPointNotContainedReason(a.p, reason, a.inclusive);
				return;
			}
		}
		
	} else if (marker == line_intersection_marker) {
		for (auto & a : lineIntersectionLits) {
			if (var(a.l) == var(p)) {
				buildLineIntersectsReason(a.line, reason, a.inclusive);
				return;
			}
		}
		
	} else if (marker == line_not_intersection_marker) {
		for (auto & a : lineIntersectionLits) {
			if (var(a.l) == var(p)) {
				buildLineNotIntersectsReason(a.line, reason, a.inclusive);
				return;
			}
		}
	} else if (marker == convex_intersection_marker) {
		for (auto & a : convexIntersectionLits) {
			if (var(a.l) == var(p)) {
				buildConvexIntersectsReason(a.polygon, reason, a.inclusive);
				return;
			}
		}
		
	} else if (marker == convex_not_intersection_marker) {
		for (auto & a : convexIntersectionLits) {
			if (var(a.l) == var(p)) {
				buildConvexNotIntersectsReason(a.polygon, reason, a.inclusive);
				return;
			}
		}
	} else if (marker == point_on_hull_marker) {
		for (auto & a : pointOnHullLits) {
			if (var(a.l) == var(p)) {
				buildPointOnHullOrDisabledReason(a.pointVar, a.p, reason);
				return;
			}
		}
	} else if (marker == point_not_on_hull_marker) {
		for (auto & a : pointOnHullLits) {
			if (var(a.l) == var(p)) {
				buildPointNotOnHullOrDisabledReason(a.pointVar, a.p, reason);
				return;
			}
		}
	} else {
		assert(false);
	}
	assert(false);
}

template<unsigned int D, class T>
bool ConvexHullDetector<D, T>::propagate(vec<Lit> & conflict) {
	
	static int iter = 0;
	if (++iter == 5764) {
		int a = 1;
	}
	stats_propagations++;
	over_hull->update();
	under_hull->update();
	ConvexPolygon<D, T> & p_over = over_hull->getHull();
	ConvexPolygon<D, T> & p_under = under_hull->getHull();
	/*printf("under: ");
	 for(auto & p:p_under){
	 printf("(%f,%f),",p.x,p.y);
	 }
	 printf("\nOver ");
	 for(auto & p:p_over){
	 printf("(%f,%f),",p.x,p.y);
	 }
	 printf("\n");*/

#ifndef NDEBUG
	cout << "Round " << iter << ", " << getID() << ":\n";
	cout << "enabled_" << getID() << " = [";
	for (int i = 0; i < under.size(); i++) {
		if (under.pointEnabled(i)) {
			cout << under[i] << ", ";
		}
	}
	cout << "]\ndisabled_" << getID() << " = [";
	for (int i = 0; i < over.size(); i++) {
		if (!over.pointEnabled(i)) {
			cout << over[i] << ", ";
		}
	}
	cout << "]\n";
	cout << "under_" << getID() << " = " << under_hull->getHull() << "\n";
	cout << "over_" << getID() << " = " << over_hull->getHull() << "\n";
#endif
	
	if (areaDetectors.size()) {
		
		T over_area = p_over.getArea();
		T under_area = p_under.getArea();
		//double tover =  over_area.get_d();
		//double tunder = under_area.get_d();
		assert(under_area <= over_area);
		for (int i = 0; i < areaDetectors.size(); i++) {
			
			Lit l = areaDetectors[i].l;
			T areaGEQ = areaDetectors[i].areaGreaterEqThan;
			if (under_area >= areaGEQ) {
				//l is true
				if (outer->value(l) == l_True) {
					//do nothing
				} else if (outer->value(l) == l_Undef) {
					
					outer->enqueue(l, area_geq_marker);
				} else if (outer->value(l) == l_False) {
					conflict.push(l);
					buildAreaGEQReason(under_area, conflict);
					return false;
				}
			} else if (over_area < areaGEQ) {
				l = ~l;
				//l is true
				if (outer->value(l) == l_True) {
					//do nothing
				} else if (outer->value(l) == l_Undef) {
					outer->enqueue(l, area_not_geq_marker);
				} else if (outer->value(l) == l_False) {
					conflict.push(l);
					buildAreaLTReason(over_area, conflict);
					return false;
				}
			}
		}
	}
	
	for (int i = 0; i < pointContainedLits.size(); i++) {
		Point<D, T> & point = pointContainedLits[i].p;
		stats_bound_checks++;
		Lit l = pointContainedLits[i].l;
		
		//Before doing a full, expensive check for containment, check if the triangle of points that contained this point last time (if there were such points) are all still enabled (in which case, containment is guaranteed).
		if (checkContainingTriangle(pointContainedLits[i].under_containing_triangle, point, p_under, under,
				pointContainedLits[i].inclusive)
				|| p_under.contains(point, &pointContainedLits[i].under_containing_triangle,
						pointContainedLits[i].inclusive)) {
			//l is true
			if (outer->value(l) == l_True) {
				//do nothing
			} else if (outer->value(l) == l_Undef) {
				
				outer->enqueue(l, point_contained_marker);
			} else if (outer->value(l) == l_False) {
				conflict.push(l);
				buildPointContainedReason(point, conflict, pointContainedLits[i].inclusive);
				return false;
			}
		} else if (!checkContainingTriangle(pointContainedLits[i].over_containing_triangle, point, p_over, over,
				pointContainedLits[i].inclusive)
				&& !p_over.contains(point, &pointContainedLits[i].over_containing_triangle,
						pointContainedLits[i].inclusive)) {
			//}else if (!p_over.contains(point)){
			l = ~l;
			//l is true
			if (outer->value(l) == l_True) {
				//do nothing
			} else if (outer->value(l) == l_Undef) {
				outer->enqueue(l, point_not_contained_marker);
			} else if (outer->value(l) == l_False) {
				conflict.push(l);
				buildPointNotContainedReason(point, conflict, pointContainedLits[i].inclusive);
				return false;
			}
		}
	}
	for (int i = 0; i < lineIntersectionLits.size(); i++) {
		
		LineSegment<D, T> & line = lineIntersectionLits[i].line;
		Lit l = lineIntersectionLits[i].l;
		
		if (checkLineIntersection(lineIntersectionLits[i].under_intersecting_polygon, line, p_under, under,
				lineIntersectionLits[i].inclusive)
				|| p_under.intersects(line, &lineIntersectionLits[i].under_intersecting_polygon, nullptr,
						lineIntersectionLits[i].inclusive)) {
			//l is true
			if (outer->value(l) == l_True) {
				//do nothing
			} else if (outer->value(l) == l_Undef) {
				
				outer->enqueue(l, line_intersection_marker);
			} else if (outer->value(l) == l_False) {
				conflict.push(l);
				buildLineIntersectsReason(line, conflict, lineIntersectionLits[i].inclusive);
				return false;
			}
		} else if (!checkLineIntersection(lineIntersectionLits[i].over_intersecting_polygon, line, p_over, over,
				lineIntersectionLits[i].inclusive)
				&& !p_over.intersects(line, &lineIntersectionLits[i].over_intersecting_polygon, nullptr,
						lineIntersectionLits[i].inclusive)) {
			//}else if (!p_over.intersects(line)){
			l = ~l;
			//l is true
			if (outer->value(l) == l_True) {
				//do nothing
			} else if (outer->value(l) == l_Undef) {
				outer->enqueue(l, line_not_intersection_marker);
			} else if (outer->value(l) == l_False) {
				conflict.push(l);
				buildLineNotIntersectsReason(line, conflict, lineIntersectionLits[i].inclusive);
				return false;
			}
		}
	}
	
	for (int i = 0; i < convexIntersectionLits.size(); i++) {
		ConvexPolygon<D, T> & polygon = convexIntersectionLits[i].polygon;
		Lit l = convexIntersectionLits[i].l;
		
		if (p_under.intersects(polygon, convexIntersectionLits[i].inclusive)) {
			//l is true
			if (outer->value(l) == l_True) {
				//do nothing
			} else if (outer->value(l) == l_Undef) {
				
				outer->enqueue(l, convex_intersection_marker);
			} else if (outer->value(l) == l_False) {
				conflict.push(l);
				buildConvexIntersectsReason(polygon, conflict, convexIntersectionLits[i].inclusive);
				return false;
			}
		} else if (!p_over.intersects(polygon, convexIntersectionLits[i].inclusive)) {
			l = ~l;
			//l is true
			if (outer->value(l) == l_True) {
				//do nothing
			} else if (outer->value(l) == l_Undef) {
				outer->enqueue(l, convex_not_intersection_marker);
			} else if (outer->value(l) == l_False) {
				conflict.push(l);
				buildConvexNotIntersectsReason(polygon, conflict, convexIntersectionLits[i].inclusive);
				return false;
			}
		}
	}
	if (pointOnHullLits.size()) {
#ifndef NDEBUG
		for (bool b : under_hull_members)
			assert(!b);
		for (bool b : over_hull_members)
			assert(!b);
#endif
		under_hull_members.growTo(under.size());
		over_hull_members.growTo(over.size());
		for (auto & p : p_under) {
			int pointsetIndex = outer->getPointsetIndex(p.getID());
			under_hull_members[pointsetIndex] = true;
		}
		for (auto & p : p_over) {
			int pointsetIndex = outer->getPointsetIndex(p.getID());
			over_hull_members[pointsetIndex] = true;
		}
		for (int i = 0; i < pointOnHullLits.size(); i++) {
			Point<D, T> & point = pointOnHullLits[i].p;
			Lit l = pointOnHullLits[i].l;
			
			//check if a point is a member of the hull
			int pointsetIndex = outer->getPointsetIndex(point.getID());
			bool under_member = under_hull_members[pointsetIndex];
			bool over_member = over_hull_members[pointsetIndex];
			
			//this is backward, I think...
			if (!over.pointEnabled(pointsetIndex) || under_member) {
				//l is true
				if (outer->value(l) == l_True) {
					//do nothing
				} else if (outer->value(l) == l_Undef) {
					
					outer->enqueue(l, point_contained_marker);
				} else if (outer->value(l) == l_False) {
					conflict.push(l);
					buildPointOnHullOrDisabledReason(pointOnHullLits[i].pointVar, point, conflict);
					return false;
				}
			} else if (!over_member) {
				assert(p_over.contains(point, false));
				l = ~l;
				//l is true
				if (outer->value(l) == l_True) {
					//do nothing
				} else if (outer->value(l) == l_Undef) {
					outer->enqueue(l, point_not_contained_marker);
				} else if (outer->value(l) == l_False) {
					conflict.push(l);
					buildPointNotOnHullOrDisabledReason(pointOnHullLits[i].pointVar, point, conflict);
					return false;
				}
			}
		}
		for (auto & p : p_under) {
			under_hull_members[p.getID()] = false;
		}
		for (auto & p : p_over) {
			over_hull_members[p.getID()] = false;
		}
	}
	return true;
}
template<unsigned int D, class T>
bool ConvexHullDetector<D, T>::checkSatisfied() {
	
	std::vector<Point<D, T>> enabled_points;
	MonotoneConvexHull<D, T> cv_under(under);
	MonotoneConvexHull<D, T> cv_over(over);
	ConvexPolygon<D, T> & h_under = under_hull->getHull();
	ConvexPolygon<D, T> & h_over = over_hull->getHull();
	/*	printf("Hull %d: ", getID());
	 for (auto & p:h_under){
	 cout << "(" << p.x <<","<<p.y <<")";
	 }
	 printf("\n");*/

	//T expectOver = h_under.getArea();
	//T expectUnder =h_over.getArea();
	T area_under = cv_under.getHull().getArea();
	T area_over = cv_over.getHull().getArea();
	//assert(equal_epsilon(area, expectOver));
	//assert(equal_epsilon(area, expectUnder));
	for (auto & a : areaDetectors) {
		T & area_cmp = a.areaGreaterEqThan;
		Lit l = a.l;
		if (outer->value(l) == l_True) {
			if (area_over < area_cmp) {
				
				cout << getID() << "Failed on area under" << area_over << "," << area_cmp << "\n";
				return false;
			}
		} else if (outer->value(l) == l_False) {
			if (area_under >= area_cmp) {
				cout << getID() << "Failed on area over" << area_over << "," << area_cmp << "\n";
				return false;
			}
		}
	}
	for (auto & a : pointContainedLits) {
		
		Lit l = a.l;
		if (!h_over.contains(a.p, a.inclusive)) {
			if (outer->value(l) != l_False) {
				cout << getID() << "Failed on point inclusion " << a.p << ", inclusive " << a.inclusive << "\n";
				return false;
			}
		} else if (h_under.contains(a.p, a.inclusive)) {
			if (outer->value(l) != l_True) {
				cout << getID() << "Failed on point exclusion " << a.p << ", inclusive " << a.inclusive << "\n";
				return false;
			}
		}
	}
	for (auto & a : lineIntersectionLits) {
		
		Lit l = a.l;
		if (!h_over.intersects(a.line, a.inclusive)) {
			if (outer->value(l) != l_False) {
				cout << getID() << "Failed on line inclusion " << a.line << ", inclusive " << a.inclusive << "\n";
				return false;
			}
		} else if (h_under.intersects(a.line, a.inclusive)) {
			if (outer->value(l) != l_True) {
				cout << getID() << "Failed on line exclusion " << a.line << ", inclusive " << a.inclusive << "\n";
				return false;
			}
		}
		
	}
	for (auto & a : convexIntersectionLits) {
		
		Lit l = a.l;
		if (!h_over.intersects(a.polygon, a.inclusive)) {
			if (outer->value(l) != l_False) {
				cout << getID() << "Failed on polygon inclusion " << a.polygon << ", inclusive " << a.inclusive << "\n";
				return false;
			}
		} else if (h_under.intersects(a.polygon, a.inclusive)) {
			if (outer->value(l) != l_True) {
				cout << getID() << "Failed on polygon exclusion " << a.polygon << ", inclusive " << a.inclusive << "\n";
				return false;
			}
		}
		/*	if(outer->value(l)==l_True){
		 if(!h_under.intersects(a.polygon,a.inclusive)){
		 cout << getID() << "Failed on polygon inclusion " << a.polygon << ", inclusive " << a.inclusive <<"\n";
		 return false;
		 }
		 }else if(outer->value(l)==l_False){
		 if(h_over.intersects(a.polygon,a.inclusive)){
		 cout << getID() << "Failed on polygon exclusion " << a.polygon << ", inclusive " << a.inclusive <<"\n";
		 return false;
		 }
		 }*/
	}
	for (auto & a : pointOnHullLits) {
		Lit l = a.l;
		//If the point is disabled, then it is counted as being in the hull.
		if (!over.pointEnabled(outer->getPointsetIndex(a.p.getID()))) {
			if (outer->value(l) == l_False) {
				cout << getID() << "Failed on point-in-hull exclusion " << a.p << "\n";
				return false;
			} else {
				continue;
			}
		}
		if (outer->value(l) == l_True) {
			bool hasPoint = false;
			for (auto & p : h_over) {
				if (p.getID() == a.p.getID()) {
					hasPoint = true;
					break;
				}
			}
			if (!hasPoint) {
				cout << getID() << "Failed on point-in-hull exclusion " << a.p << "\n";
				return false;
			}
			
		} else if (outer->value(l) == l_False) {
			/*bool hasPoint = false;
			 for(auto & p:h_under){
			 if(p.getID()==a.p.getID()){
			 hasPoint=true;
			 break;
			 }
			 }
			 if(hasPoint)
			 return false;*/
			if (!h_under.contains(a.p, false)) {
				cout << getID() << "Failed on point-in-hull inclusion " << a.p << "\n";
				return false;
			}
		}
	}
	return true;
}
template<unsigned int D, class T>
Lit ConvexHullDetector<D, T>::decide() {
	
	return lit_Undef;
}

template<unsigned int D, class T>
void ConvexHullDetector<D, T>::buildAreaGEQReason(T area, vec<Lit> & conflict) {
	//the reason that the area is greater or equal to the current value is the set of points in the convex hull (all of which are enabled).
	assert(under_hull->getHull().getArea() >= area);
	under_hull->update();
	for (auto & p : under_hull->getHull()) {
		int pID = p.getID();
		Lit l = mkLit(outer->getPointVar(pID), true);
		assert(outer->value(l)==l_False);
		conflict.push(l);
	}
}
template<unsigned int D, class T>
void ConvexHullDetector<D, T>::buildAreaLTReason(T area, vec<Lit> & conflict) {
	//the reason that the area is less than some value is that some point that is OUTSIDE the convex hull is not enabled.
	for (int i = 0; i < over.size(); i++) {
		if (!over.pointEnabled(i) && !over_hull->getHull().contains(over[i], true)) {
			conflict.push(mkLit(outer->getPointVar(over[i].getID()), false));
		}
	}
}
template<unsigned int D, class T>
void ConvexHullDetector<D, T>::buildPointContainedReason(const Point<D, T> & s, vec<Lit> & conflict, bool inclusive) {
	if (D == 2) {
		buildPointContainedReason2d((const Point<2, T>&) s, conflict, inclusive);
	}
	
}

template<unsigned int D, class T>
void ConvexHullDetector<D, T>::buildPointNotContainedReason(const Point<D, T> & s, vec<Lit> & conflict,
		bool inclusive) {
	if (D == 2) {
		buildPointNotContainedReason2d((const Point<2, T>&) s, conflict, inclusive);
	}
}
template<unsigned int D, class T>
void ConvexHullDetector<D, T>::buildLineIntersectsReason(LineSegment<D, T> & line, vec<Lit> & conflict,
		bool inclusive) {
	/*tmp_polygon.clear();
	 tmp_polygon.addVertex(line.a);
	 tmp_polygon.addVertex(line.b);*/
	buildConvexIntersectsReason(line, conflict, inclusive);
}
template<unsigned int D, class T>
void ConvexHullDetector<D, T>::buildLineNotIntersectsReason(LineSegment<D, T> & line, vec<Lit> & conflict,
		bool inclusive) {
	/*tmp_polygon.clear();
	 tmp_polygon.addVertex(line.a);
	 tmp_polygon.addVertex(line.b);*/
	buildConvexNotIntersectsReason(line, conflict, inclusive);
}

template<unsigned int D, class T>
void ConvexHullDetector<D, T>::buildConvexIntersectsReason(ConvexPolygon<D, T> & polygon, vec<Lit> & conflict,
		bool inclusive) {
	if (D == 2) {
		buildConvexIntersectsReason2d((ConvexPolygon<2, T>&) polygon, conflict, inclusive);
	}
}
template<unsigned int D, class T>
void ConvexHullDetector<D, T>::buildConvexNotIntersectsReason(ConvexPolygon<D, T> & polygon, vec<Lit> & conflict,
		bool inclusive) {
	if (D == 2) {
		buildConvexNotIntersectsReason2d((ConvexPolygon<2, T>&) polygon, conflict, inclusive);
	}
}

template<unsigned int D, class T>
void ConvexHullDetector<D, T>::buildPointOnHullOrDisabledReason(Var pointVar, const Point<D, T> & p,
		vec<Lit> & conflict) {
	if (D == 2) {
		buildPointOnHullOrDisabledReason2d(pointVar, (const Point<2, T>&) p, conflict);
	}
}

template<unsigned int D, class T>
void ConvexHullDetector<D, T>::buildPointNotOnHullOrDisabledReason(Var pointVar, const Point<D, T> & p,
		vec<Lit> & conflict) {
	if (D == 2) {
		buildPointNotOnHullOrDisabledReason2d(pointVar, (const Point<2, T>&) p, conflict);
	}
}

template<unsigned int D, class T>
void ConvexHullDetector<D, T>::buildPointContainedReason2d(const Point<2, T> & s, vec<Lit> & conflict, bool inclusive) {
	//A point is contained because there exists three enabled points that form a containing triangle
	//(might have to handle the edge case of a 2 point hull depending on how we interpret such a hull)
	//its not immediately clear how we would pick which 3 points to include, given that there may be a choice
	ConvexPolygon<2, T> &hull = under_hull->getHull();
	assert(hull.contains(s, inclusive));
	
	NConvexPolygon<2, T> triangle;
	hull.contains(s, &triangle, inclusive);
	
	assert(triangle.contains(s, inclusive));
	/*	printf("Learn must disable one of: ");
	 for(auto & p:triangle){
	 printf("(%f, %f),", p.x,p.y);
	 }
	 printf("\n");*/
	for (auto & p : triangle) {
		int id = p.getID();
		Var v = outer->getPointVar(id);
		assert(outer->value(v)==l_True);
		conflict.push(mkLit(v, true));
	}
	stats_under_clauses++;
	stats_under_clause_length += conflict.size();
}

template<unsigned int D, class T>
void ConvexHullDetector<D, T>::buildPointNotContainedReason2d(const Point<2, T> & s, vec<Lit> & conflict,
		bool inclusive) {
	//find ANY line/plane/hyperplane that passes through s, and does NOT intersect the hull.
	//let the side of the plane that contains the hull be the 'near' side
	//any set of points that are strictly on the near side of the plane (not including points on the plane itself)
	//will also have a convex hull that is strictly on the near side of the plane.
	//it follows that at least one point that is either on the plane, or on the 'far' side of the plane, must be included in the point set
	//in order for p to be contained in the convex hull.
	
	//we can pick any plane; for now, we will search for a plane that minimizes the number of points on the far side (or on the plane exactly)
	
	//first, find the extreme points of the hull relative to the point.
	//we could find this by adding the point to the pointset and finding the convex hull, and then finding the neighbours of the point...
	
	//for 2D, we do this by forming the line (s,p), where p is any disabled point. Then we check if that line intersects the hull; if it doesn't
	
	Line<2, T> testline;
	
	testline.a = s;
	ConvexPolygon<2, T> &hull = over_hull->getHull();
	assert(!hull.contains(s, inclusive));
	/*	for(auto & p:hull){
	 printf("(%f, %f), ", p.x, p.y);
	 }
	 printf("\n");*/
	//edge cases:
	if (hull.size() <= 1) {
		//then we just report that _some_ disabled vertex must be enabled. For hull.size()==1, we can do better than this using the Separating axis theorem.
		for (int i = 0; i < over.size(); i++) {
			Point<2, T> & p = over[i];
			if (!over.pointEnabled(i)) {
				conflict.push(mkLit(outer->getPointVar(p.getID()), false));
			}
		}
		stats_over_clauses++;
		stats_over_clause_length += conflict.size();
		return;
	}/*else if(hull.size()==2){
	 //then we just report that some disabled vertex on the same side of the line as this vertex must be disabled.
	 LineSegment<2,T> testline(hull[0],hull[1]);

	 //now check to see if test line intersects the hull (it can still intersect on the otherside
	 if(!testline.intersects(hull)){
	 //then this is a point outside the hull that we can try testing.
	 findFarPoints(testline,test_set,min_set,hull,true);
	 }
	 return;
	 }*/
	
	Line<2, T> top_line;
	Line<2, T> bottom_line;
	top_line.a = s;
	bottom_line.a = s;
	
	bool found = false;
	for (auto & p : hull) {
		if (p != s) {
			found = true;
			top_line.b = p;
			bottom_line.b = p;
		}
	}
	if (!found) {
		assert(!inclusive);	//because if the collision is inclusive, then a hull consisting only of points s is considered to contain s.
				
		//the hull consists entirely of points that are exactly the point we are checking against.
		//In that case, just learn that some point that is not equal to s must be added
		for (int i = 0; i < over.size(); i++) {
			Point<2, T> & p = over[i];
			
			if (!over.pointEnabled(i)) {
				if (p == s)
					continue;
				assert(!hull.contains(p, inclusive));
				Point<2, T> & p = over[i];
				assert(p.getID() >= 0);
				assert(!over.pointEnabled(outer->getPointsetIndex(p.getID())));
				Lit l = mkLit(outer->getPointVar(p.getID()), false);
				assert(outer->value(l)==l_False);
				conflict.push(l);
			} else {
				assert(hull.contains(p, true));
			}
		}
		return;
	}
	
	/*	top_line.b = hull[0];
	 bottom_line.b = hull.back();*/

	std::vector<Point<2, T> > min_set;
	std::vector<Point<2, T> > test_set;
	hasSet = false;
	
	for (int i = 0; i < over.size(); i++) {
		Point<2, T> & p = over[i];
		
		if (!over.pointEnabled(i)) {
			if (p == s)
				continue;
			if (!hull.contains(p, inclusive)) {
				
				testline.b = p;
				
				//now check to see if test line intersects the hull (it can still intersect on the otherside
				if (!testline.intersects(hull, inclusive)) {
					//then this is a point outside the hull that we can try testing.
					findFarPoints(testline, test_set, min_set, hull, true);
				}
			}
		} else {
			assert(hull.contains(p, true));
		}
	}
	//find the extreme points on the hull itself
	for (auto & p : hull) {
		if (top_line.whichSide(p) > 0) {
			top_line.b = p;
		}
		
		if (bottom_line.whichSide(p) < 0) {
			bottom_line.b = p;
		}
	}
	
	findFarPoints(top_line, test_set, min_set, hull, true);
	
	findFarPoints(bottom_line, test_set, min_set, hull, true);
	/*	printf("Learn must add: ");
	 for(auto & p:min_set){
	 printf("(%f, %f), ", p.x,p.y);
	 }
	 printf("\n");*/
	for (int i = 0; i < min_set.size(); i++) {
		Point<2, T> & p = min_set[i];
		assert(p.getID() >= 0);
		assert(!over.pointEnabled(outer->getPointsetIndex(p.getID())));
		Lit l = mkLit(outer->getPointVar(p.getID()), false);
		assert(outer->value(l)==l_False);
		conflict.push(l);
		
	}
	stats_over_clauses++;
	stats_over_clause_length += conflict.size();
}
template<unsigned int D, class T>
void ConvexHullDetector<D, T>::buildPointOnHullOrDisabledReason2d(Var pointVar, const Point<2, T> & s,
		vec<Lit> & conflict) {
	//If the point IS on the hull (or disabled), it is either because it is currently disabled...
	
	assert(outer->getPointVar(s.getID()) == pointVar);
	Lit enabledLit = mkLit(pointVar, false);
	
	if (outer->value(enabledLit) == l_False) {
		conflict.push(enabledLit);
		return;
	}
	
	//or because there is some line we can form that passes through this point, such that no enabled point is on the far side of the line, and the entire hull (except for this point)
	//is on the near side.
	//Note that as we are forcing the convex hull to contain all collinear points on the boundary,
	//this code is handles the case where enabled points are exactly on this line in the opposite way as  'buildPointNotContainedReason' does
	
	Line<2, T> testline;
	
	testline.a = s;
	ConvexPolygon<2, T> &hull = over_hull->getHull();
	
	//edge cases:
	if (hull.size() < 3) {
		//then we just report that _some_ disabled vertex must be enabled
		for (int i = 0; i < over.size(); i++) {
			Point<2, T> & p = over[i];
			if (!over.pointEnabled(i)) {
				conflict.push(mkLit(outer->getPointVar(p.getID()), false));
			}
		}
		return;
	}
	
	Line<2, T> top_line;
	Line<2, T> bottom_line;
	top_line.a = s;
	bottom_line.a = s;
	top_line.b = hull[0];
	bottom_line.b = hull.back();
	
	std::vector<Point<2, T> > min_set;
	std::vector<Point<2, T> > test_set;
	hasSet = false;
	
	for (int i = 0; i < over.size(); i++) {
		Point<2, T> & p = over[i];
		if (!over.pointEnabled(i)) {
			if (!hull.contains(p, true)) {
				testline.b = p;
				
				//now check to see if test line intersects the hull (it can still intersect on the otherside
				if (!testline.intersects(hull, true)) {
					//then this is a point outside the hull that we can try testing.
					//Note that as we are forcing the convex hull to contain all collinear points on the boundary,
					//this code is handles the case where enabled points are exactly on this line in the opposite way as  'buildPointNotContainedReason' does
					
					findFarPoints(testline, test_set, min_set, hull, false);
				}
			}
		} else {
			assert(hull.contains(p, true));
		}
	}
	//find the extreme points on the hull itself
	for (auto & p : hull) {
		if (top_line.whichSide(p) >= 0) {
			top_line.b = p;
		}
		
		if (bottom_line.whichSide(p) <= 0) {
			bottom_line.b = p;
		}
	}
	findFarPoints(top_line, test_set, min_set, hull, true);
	findFarPoints(bottom_line, test_set, min_set, hull, true);
	
	for (int i = 0; i < min_set.size(); i++) {
		Point<2, T> & p = min_set[i];
		assert(!over.pointEnabled(p.getID()));
		Lit l = mkLit(outer->getPointVar(p.getID()), false);
		assert(outer->value(l)==l_False);
		conflict.push(l);
		
	}
}

template<unsigned int D, class T>
void ConvexHullDetector<D, T>::buildPointNotOnHullOrDisabledReason2d(Var pointVar, const Point<2, T> & s,
		vec<Lit> & conflict) {
	//If the point is not (on the hull or disabled), it is because it is both currently enabled...
	assert(outer->getPointVar(s.getID()) == pointVar);
	Lit enabledLit = mkLit(pointVar, false);
	assert(outer->value(enabledLit)==l_True);
	conflict.push(~enabledLit);
	//AND, it is because there exist three other enabled points that form a triangle that contain this point:
	ConvexPolygon<2, T> &hull = under_hull->getHull();
#ifndef NDEBUG
	for (auto & p : hull)
		assert(p.getID() != s.getID());
#endif
	NConvexPolygon<2, T> triangle;
	hull.contains(s, &triangle, true);
	//findContainingTriangle2d(hull,s,triangle,true);
	assert(triangle.contains(s, true));
	for (auto & p : triangle) {
		int id = p.getID();
		Var v = outer->getPointVar(id);
		assert(outer->value(v)==l_True);
		conflict.push(mkLit(v, true));
	}
	
}

template<unsigned int D, class T>
void ConvexHullDetector<D, T>::buildConvexNotIntersectsReason2d(ConvexPolygon<2, T> & polygon, vec<Lit> & conflict,
		bool inclusive) {
	ConvexPolygon<2, T> & h1 = over_hull->getHull();
	ConvexPolygon<2, T> & h2 = polygon;
	/*printf("h1: ");
	 for (auto & p:h1){
	 cout<<p << " ";
	 }
	 printf("\nh2: ");
	 for (auto & p:h2){
	 cout<<p << " ";
	 }
	 printf("\n");*/
	/*

	 static int iter = 0;
	 if(++iter==18){
	 int a=1;
	 }
	 */

	assert(!h1.intersects(h2, inclusive));
	assert(!h2.intersects(h1, inclusive));
	
	if (h1.size() == 0) {
		//then the reason that the shapes don't collide is that at least one point in each of them must be enabled
		//can we improve on this?
		for (int i = 0; i < h1.size(); i++) {
			if (!over.pointEnabled(i)) {
				Lit l = mkLit(outer->getPointVar(over[i].getID()));
				assert(outer->value(l)==l_False);
				conflict.push(l);
			}
		}
		return;
	} else if (h2.size() == 0) {
		return;
	} else if (h2.size() == 1) {
		buildPointNotContainedReason(h2[0], conflict, inclusive);
	} else if (h1.size() == 1) {
		assert(h1[0] != h2[0]);		//identical points are considered to collide
	}
	
	//the reason that the two convex hulls do NOT collide is that there exists a separating axis between them.
	//Find that axis, and then find all the positions of the disabled points of both polygons along that axis.
	//let h1 be left of h2.
	//Then either some disabled point in h1 that is to the right of the rightmost point in h1 must be enabled,
	//or some disabled point in h2 that is to the left of the leftmost point in h2 must be enabled.
	
	//Note: It may be possible to improve on this analysis!
	static std::vector<std::pair<Point<2, T>, T>> projection;
	projection.clear();
	static std::vector<std::pair<Point<2, T>, T>> projection2;
	projection2.clear();
	bool found = findSeparatingAxis2d(h1, h2, over, projection, projection2, inclusive);
	assert(found);
	if (!found || (projection.size() == 0 || projection2.size() == 0)) {
		cout << "Error! Failed to find separating axis between hulls! Aborting.\n";
		cout << h1 << "\n";
		cout << h2 << "\n";
		exit(4);
	}
	T leftmost1 = numeric<T>::infinity();
	T rightmost1 = -numeric<T>::infinity();
	
	for (auto & p : projection) {
		int pID = p.first.getID();
		int pointset = outer->getPointset(pID);
		assert(pointset == over.getID());
		int pointsetIndex = outer->getPointsetIndex(pID);
		if (over.pointEnabled(pointsetIndex)) {
			if (p.second < leftmost1)
				leftmost1 = p.second;
			if (p.second > rightmost1)
				rightmost1 = p.second;
		}
	}
	
	T leftmost2 = numeric<T>::infinity();
	T rightmost2 = -numeric<T>::infinity();
	
	for (auto & p : projection2) {
		if (p.second < leftmost2)
			leftmost2 = p.second;
		if (p.second > rightmost2)
			rightmost2 = p.second;
	}
	bool h1_is_left;
	bool h1_is_right;
	if (inclusive) {
		assert(rightmost1 < leftmost2 || rightmost2 < leftmost1);
		h1_is_left = rightmost1 < leftmost2;
		h1_is_right = rightmost1 > leftmost2;
	} else {
		assert(rightmost1 <= leftmost2 || rightmost2 <= leftmost1);
		h1_is_left = rightmost1 <= leftmost2;
		h1_is_right = rightmost1 >= leftmost2;
	}
	
	for (auto & p : projection) {
		int pID = p.first.getID();
		int pointset = outer->getPointset(pID);
		int pointsetIndex = outer->getPointsetIndex(pID);
		assert(pointset == over.getID());
		if (!over.pointEnabled(pointsetIndex)) {
			if (!over.pointEnabled(pointsetIndex)) {
				if (h1_is_left && p.second > rightmost1) {
					//then if we enable this point, the two hulls will move closer to each other.
					//we can probably improve on this, by only considering points that are also >= the rightmost _disabled_ point of pointset2...
					Lit l = mkLit(outer->getPointVar(pID));
					assert(outer->value(l)==l_False);
					conflict.push(l);
				}
				if (h1_is_right && p.second < leftmost1) {
					Lit l = mkLit(outer->getPointVar(pID));
					assert(outer->value(l)==l_False);
					conflict.push(l);
				}
			}
		}
	}
	
}

template<unsigned int D, class T>
void ConvexHullDetector<D, T>::buildConvexIntersectsReason2d(ConvexPolygon<2, T> & polygon, vec<Lit> & conflict,
		bool inclusive) {
	/*	static int iter = 0;
	 if(++iter==218){
	 int a=1;
	 }
	 printf("iter %d\n",iter);*/
	//If the two polygons intersect, there are two possible cases (these are not mutually exclusive)
	//1) 1 polygon has at least one point contained in the other polygon. Learn that either that point must be disabled, or at least
	//		one point from its containing triangle must be disabled (as in the point containment theory).
	//2) There exists an edge (or, more generally, a line segment between any two vertices, as opposed to just an edge) from each poygon,
	//		such that the two edges intersect. Learn that one of the end points must be disabled.
	ConvexPolygon<2, T> & h1 = under_hull->getHull();
	ConvexPolygon<2, T> & h2 = polygon;
	//h1.intersects(h2,inclusive);
	/*printf("h1: ");
	 for (auto & p:h1){
	 cout<<p << " ";
	 }
	 printf("\nh2: ");
	 for (auto & p:h2){
	 cout<<p << " ";
	 }
	 printf("\n");*/
	//so, first, check each edge segment to see if there is an intersecting line.
	//I'm doing this first, on the hypothesis that these line intersection constraints may form better learned clauses.
	//(since cases 1 and 2 overlap, we have to make a choice about which to favour...)
	assert(h1.intersects(h2, inclusive));
	assert(h2.intersects(h1, inclusive));
	static NConvexPolygon<2, T> intersect1;
	
	intersect1.clear();
	
	bool check = h1.intersects(h2, &intersect1, nullptr, inclusive);
	assert(check);
	if (!check || intersect1.size() == 0) {
		exit(3);
	}
	for (auto & p : intersect1) {
		conflict.push(~mkLit(outer->getPointVar(p.getID())));
	}

	//it is probably possible to improve on this quadratic time search...
	/*for(int i = 0; i<h1.size();i++){
	 Point<2,T> & prev = h1[i-1];
	 Point<2,T> & p = h1[i];
	 LineSegment<2,T> edge1(prev,p);
	 for(int j = 0;  j<h2.size();j++){
	 Point<2,T> & prev2 = h2[j-1];
	 Point<2,T> & p2 = h2[j];
	 LineSegment<2,T> edge2(prev2,p2);
	 if(edge1.intersects(edge2,inclusive)){
	 //learn that one of the endpoints of the intersecting lines must be disabled
	 conflict.push(~mkLit(outer->getPointVar(prev.getID())));
	 conflict.push(~mkLit(outer->getPointVar(p.getID())));
	 return;
	 }
	 }
	 }

	 //if no intersecting line segment was found, then it follows that one of the polygons is wholly contained in the other.
	 //so treat this as a contained point problem - pick one of the points in the contained polygon (arbitrarily), and a containing triangle
	 //from the other polygon, and learn that one of these 4 points must be disabled.


	 for(int i = 0; i<h1.size();i++){
	 if(h2.contains(h1[i],inclusive)){
	 static NConvexPolygon<2,T> triangle;
	 h2.contains(h1[i],&triangle,inclusive);
	 //findContainingTriangle2d(h2,h1[i],triangle,inclusive);
	 assert(triangle.contains(h1[i],inclusive));
	 conflict.push(~mkLit(outer->getPointVar(h1[i].getID())));
	 return;
	 }
	 }

	 for(int i = 0; i<h2.size();i++){
	 if(h1.contains(h2[i],inclusive)){
	 static NConvexPolygon<2,T> triangle;
	 h1.contains(h2[i],&triangle,inclusive);
	 //findContainingTriangle2d(h1,h2[i],triangle,inclusive);
	 assert(triangle.contains(h2[i],inclusive));
	 for(auto & p: triangle){
	 int id = p.getID();
	 Var v = outer->getPointVar(id);
	 assert(outer->value(v)==l_True);
	 conflict.push(mkLit(v,true));
	 }
	 assert(conflict.size());
	 return;
	 }
	 }
	 */
	//if collisions are not inclusive, then there is one more possibility - that the polygons collide, but only on edges whose endpoints are shared with each other.
	//so, find a line that intersects the polygon, whose endpoints are exactly on the polygon. the collision
	assert(conflict.size());
}

template<unsigned int D, class T>
bool ConvexHullDetector<D, T>::findSeparatingAxis2d(ConvexPolygon<2, T> & hull1, ConvexPolygon<2, T> & hull2,
		PointSet<2, T> & pointset1, std::vector<std::pair<Point<2, T>, T>> &projection_out1_t,
		std::vector<std::pair<Point<2, T>, T>> &projection_out2_t, bool inclusive) {
	
	ConvexPolygon<2, T> & h1 = (hull1.size() <= hull2.size()) ? hull1 : hull2;
	ConvexPolygon<2, T> & h2 = (hull1.size() <= hull2.size()) ? hull2 : hull1;
	std::vector<std::pair<Point<2, T>, T>> &projection_out1 =
			(hull1.size() <= hull2.size()) ? projection_out1_t : projection_out2_t;
	std::vector<std::pair<Point<2, T>, T>> &projection_out2 =
			(hull1.size() <= hull2.size()) ? projection_out2_t : projection_out1_t;
	
	projection_out1.clear();
	projection_out2.clear();
	if (h1.size() == 0 || h2.size() == 0) {
		return false;
	} else if (h1.size() == 1 && h2.size() == 1) {
		//what is the correct response for exclusive collisions, here?
		Point<2, T> un_normalized_normal = h2[0] - h1[0];
		//now place all the remaining (disabled) points on the axis as well
		for (int i = 0; i < pointset1.size(); i++) {
			auto & p = pointset1[i];
			T projection = un_normalized_normal.dot(p);
			projection_out1.push_back( { p, projection });
		}
		for (int i = 0; i < h2.size(); i++) {
			auto & p = h2[i];
			T projection = un_normalized_normal.dot(p);
			projection_out2.push_back( { p, projection });
		}
		return true;
	}
	
	//Separating Axis Theorem for collision detection between two convex polygons
	//loop through each edge in _each_ polygon and project both polygons onto that edge's normal.
	//If any of the projections are non-intersection, then these don't collide; else, they do collide
	if (h1.size() > 1) {
		for (int i = 0; i < h1.size(); i++) {
			auto & p = h1[i];
			auto & prev = (h1)[i - 1];
			Point<2, T> edge = p - prev;
			Point<2, T> un_normalized_normal(-edge.y, edge.x);
			projection_out1.clear();
			projection_out2.clear();
			//now project both polygons onto to this normal and see if they overlap, by finding the minimum and maximum distances
			//Note that since we are NOT normalizing the normal vector, the projection is distorted along that vector
			//(this still allows us to check overlaps, but means that the minimum distance found between the two shapes may be incorrect)
			T left = numeric<T>::infinity();
			T right = -numeric<T>::infinity();
			for (auto & p : h1) {
				T projection = un_normalized_normal.dot(p);
				projection_out1.push_back( { p, projection });
				if (projection < left) {
					left = projection;
				}
				if (projection > right) {
					right = projection;
				}
			}
			bool seenLeft = false;
			bool seenRight = false;
			
			for (auto & p : h2) {
				T projection = un_normalized_normal.dot(p);
				projection_out2.push_back( { p, projection });
				if (inclusive) {
					if (projection >= left && projection <= right) {
						seenRight = true;
						seenLeft = true;
						break;
					} else if (projection < left) {
						seenLeft = true;
						if (seenRight) {
							break;
						}
					} else if (projection > right) {
						seenRight = true;
						if (seenLeft) {
							break;
						}
					} else if (seenLeft && projection > left) {
						seenRight = true;
						break;
					} else if (seenRight && projection < right) {
						seenRight = true;
						break;
					}
				} else {
					if (projection > left) {
						seenRight = true;
						if (seenLeft) {
							break;
						}
					}
					if (projection < right) {
						seenLeft = true;
						if (seenRight) {
							break;
						}
					}
					
				}
			}
			
			if (!(seenLeft && seenRight)) {
				//now place all the remaining (disabled) points on the axis as well
				for (int i = 0; i < pointset1.size(); i++) {
					if (!pointset1.pointEnabled(i)) {
						auto & p = pointset1[i];
						T projection = un_normalized_normal.dot(p);
						if (hull1.size() <= hull2.size()) {
							projection_out1.push_back( { p, projection });
						} else {
							projection_out2.push_back( { p, projection });
						}
					}
				}
				
				return true;
			}
		}
	}
	if (h2.size() > 1) {
		//now test the axis produced by the other polygon
		for (int i = 0; i < h2.size(); i++) {
			auto & p = (h2)[i];
			auto & prev = (h2)[i - 1];
			Point<2, T> edge = p - prev;
			Point<2, T> un_normalized_normal(-edge.y, edge.x);
			projection_out1.clear();
			projection_out2.clear();
			
			T left = numeric<T>::infinity();//we can't use this, because although these are defined for gmp mpq class, in at least the current release, they are meaningless.
			T right = -numeric<T>::infinity();
			for (auto & p : h2) {
				T projection = un_normalized_normal.dot(p);
				projection_out2.push_back( { p, projection });
				if (projection < left) {
					left = projection;
				}
				if (projection > right) {
					right = projection;
				}
			}
			
			bool seenLeft = false;
			bool seenRight = false;
			bool overlaps = false;
			for (auto & p : h1) {
				T projection = un_normalized_normal.dot(p);
				projection_out1.push_back( { p, projection });
				if (inclusive) {
					if (projection >= left && projection <= right) {
						seenRight = true;
						seenLeft = true;
						break;
					} else if (projection < left) {
						seenLeft = true;
						if (seenRight) {
							break;
						}
					} else if (projection > right) {
						seenRight = true;
						if (seenLeft) {
							break;
						}
					} else if (seenLeft && projection > left) {
						seenRight = true;
						break;
					} else if (seenRight && projection < right) {
						seenRight = true;
						break;
					}
				} else {
					if (projection > left) {
						seenRight = true;
						if (seenLeft) {
							break;
						}
					}
					if (projection < right) {
						seenLeft = true;
						if (seenRight) {
							break;
						}
					}
					
				}
			}
			if (!(seenLeft && seenRight)) {
				for (int i = 0; i < pointset1.size(); i++) {
					if (!pointset1.pointEnabled(i)) {
						auto & p = pointset1[i];
						T projection = un_normalized_normal.dot(p);
						if (hull1.size() <= hull2.size()) {
							projection_out1.push_back( { p, projection });
						} else {
							projection_out2.push_back( { p, projection });
						}
					}
				}
				return true;
			}
		}
	}
	projection_out1.clear();
	projection_out2.clear();
	return false;
}

#endif
