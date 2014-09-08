/*
 * SymbolicPolygon.h
 *
 *  Created on: Aug 2, 2014
 *      Author: sam
 */

#ifndef SYMBOLIC_POLYGON
#define SYMBOLIC_POLYGON

#include <assert.h>
#include <vector>
#include "Polygon.h"
#include "GeometryDetector.h"
#include "GeometryTypes.h"
#include "core/SolverTypes.h"
#include "PointSet.h"
#include "GeometryDetector.h"
#include "core/Config.h"
#include "Polygon.h"
#include "LineSegment.h"
#include "Line.h"
#include "HalfPlane.h"
#include "SymbolicScalar.h"
#include "Ray.h"

template<unsigned int D,class T>
class SymbolicPolygon:public GeometryDetector{
protected:
	TheorySolver * solver;
	double rnd_seed;
private:

	CRef point_contained_marker=CRef_Undef;
	CRef point_not_contained_marker=CRef_Undef;
	CRef convex_intersection_marker=CRef_Undef;
	CRef convex_not_intersection_marker=CRef_Undef;
	CRef line_intersection_marker=CRef_Undef;
	CRef line_not_intersection_marker=CRef_Undef;
	CRef area_geq_marker=CRef_Undef;
	CRef area_not_geq_marker=CRef_Undef;
	CRef point_on_hull_marker=CRef_Undef;
	CRef point_not_on_hull_marker=CRef_Undef;

	Var lowest_point_var=var_Undef;
	vec<int> point_lit_map;
	vec<Lit> lit_point_map;
	vec<bool> under_hull_members;
	vec<bool> over_hull_members;

	bool requires_propagation=false;

	NConvexPolygon<D,T> tmp_polygon;
	HalfPlane<2,T> separating_halfplane;

	struct PointContainedLit{
		Point<D,T> p;
		Lit l;
		bool inclusive;
		//indices of 3 containing vertices in the under approximation that contain this point (or -1 if no such triangle exists).
		NConvexPolygon<D,T> under_containing_triangle;
		NConvexPolygon<D,T> over_containing_triangle;

	};
	std::vector<PointContainedLit> pointContainedLits;


	struct LineIntersectionLit{
		LineSegment<D,T> line;
		Lit l;
		bool inclusive;
		NConvexPolygon<D,T> under_intersecting_polygon;
		NConvexPolygon<D,T> over_intersecting_polygon;

	};
	std::vector<LineIntersectionLit> lineIntersectionLits;

	struct ConvexIntersectionLit{
		NConvexPolygon<D,T> polygon;
		Lit l;
		bool inclusive;
	};
	std::vector<ConvexIntersectionLit> convexIntersectionLits;

	struct PointOnHullLit{
		Var pointVar;
		int pointIndex;
		Point<D,T> p;
		Lit l;
	};
	std::vector<PointOnHullLit> pointOnHullLits;

	struct AreaLit{
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

	void init_lits(){
		if(!requires_propagation){
			init_lits();
			requires_propagation=true;
			point_contained_marker=solver->newReasonMarker(getID());
			point_not_contained_marker=solver->newReasonMarker(getID());
			point_on_hull_marker=solver->newReasonMarker(getID());
			point_not_on_hull_marker=solver->newReasonMarker(getID());
			area_geq_marker=solver->newReasonMarker(getID());
			area_not_geq_marker=solver->newReasonMarker(getID());
			convex_intersection_marker=solver->newReasonMarker(getID());
			convex_not_intersection_marker=solver->newReasonMarker(getID());
			line_intersection_marker=solver->newReasonMarker(getID());
			line_not_intersection_marker=solver->newReasonMarker(getID());
		}
	}

public:

	SymbolicPolygon(int detectorID,TheorySolver * solver,  double seed=1):GeometryDetector(detectorID),solver(solver),rnd_seed(seed){

	}

	virtual ~SymbolicPolygon(){

	}

	virtual bool isConstant(){
		return false;
	}

	virtual Polygon<D,T>& getOverApprox()=0;
	virtual Polygon<D,T>& getUnderApprox()=0;
	void addAreaGEQLit(T areaGreaterEqThan, Var solverVar){
		Var v = solver->newVar(solverVar,getID());
		Lit l = mkLit(v,false);
		areaDetectors.push_back({areaGreaterEqThan,l});
		init_lits();
	}

	void addPointContainmentLit(Point<D,T> p,Var solverVar, bool inclusive){
		Var containVar = solver->newVar(solverVar,getID());
		pointContainedLits.push_back({p,mkLit(containVar,false),inclusive});
		init_lits();
	}

	void addLineIntersectionLit(LineSegment<D,T> line, Var solverVar, bool inclusive){
		Var containVar = solver->newVar(solverVar,getID());
		lineIntersectionLits.push_back({line,mkLit(containVar,false),inclusive});
		init_lits();
	}

	void addConvexIntersectionLit(NConvexPolygon<D,T> & polygon, Var solverVar, bool inclusive){
		Var containVar = solver->newVar(solverVar,getID());
		convexIntersectionLits.push_back({polygon,mkLit(containVar,false),inclusive});
		init_lits();
	}

	void addAreaDetectorLit(T areaGreaterEqThan, Var solverVar){
		Var v = solver->newVar(solverVar,getID());
		Lit l = mkLit(v,false);
		areaDetectors.push_back({areaGreaterEqThan,l});
		init_lits();
	}

	inline bool checkContainingTriangle(NConvexPolygon<D, T> & containing,const Point<2,T> & point, Polygon<D,T> & hull,  bool inclusive){
		/*if(containing.size()==0)
			return false;
		for (auto &p: containing){
			int index = solver->getPointsetIndex(p.getID());
			if(!pointset.pointEnabled(index)){
				containing.clear();
				return false;
			}
		}
		assert(containing.contains(point,inclusive));
		assert(hull.contains(point,inclusive));

		return true;
*/
	return false;
	}
	inline bool checkLineIntersection(NConvexPolygon<D, T> & containing, LineSegment<D,T> & line, Polygon<D,T> & hull,  bool inclusive){
	/*	if(containing.size()==0)
			return false;
		for (auto &p: containing){
			int index = solver->getPointsetIndex(p.getID());
			if(!pointset.pointEnabled(index)){
				containing.clear();
				return false;
			}
		}
		assert(containing.intersects(line,inclusive));
		assert(hull.intersects(line,inclusive));

		return true;*/
		return false;

	}
	bool propagate(vec<Lit> & conflict){
		if(!requires_propagation)
			return true;
		static int iter = 0;
		if(++iter==5764){
			int a=1;
		}

		stats_propagations++;

		Polygon<D,T> & p_over = getOverApprox();
		Polygon<D,T> & p_under = getUnderApprox();


		if(areaDetectors.size()){

			T over_area = p_over.getArea();
			T under_area = p_under.getArea();
			//double tover =  over_area.get_d();
			//double tunder = under_area.get_d();
			assert(under_area<=over_area);
			for(int i = 0;i<areaDetectors.size();i++){

				Lit l = areaDetectors[i].l;
				T areaGEQ = areaDetectors[i].areaGreaterEqThan;
				if(under_area>=areaGEQ){
					//l is true
					if(solver->value(l)==l_True){
						//do nothing
					}else if(solver->value(l)==l_Undef){

						solver->enqueue(l,area_geq_marker) ;
					}else if (solver->value(l)==l_False){
						conflict.push(l);
						buildAreaGTReason(under_area,conflict);
						return false;
					}
				}else if (over_area<areaGEQ){
					l=~l;
					//l is true
					if(solver->value(l)==l_True){
						//do nothing
					}else if(solver->value(l)==l_Undef){
						solver->enqueue(l,area_not_geq_marker) ;
					}else if (solver->value(l)==l_False){
						conflict.push(l);
						buildAreaLEQReason(over_area,conflict);
						return false;
					}
				}
			}
		}



		for(int i =0;i<pointContainedLits.size();i++){
			Point<D,T> & point = pointContainedLits[i].p;
			stats_bound_checks++;
			Lit l = pointContainedLits[i].l;

			//Before doing a full, expensive check for containment, check if the triangle of points that contained this point last time (if there were such points) are all still enabled (in which case, containment is guaranteed).
			if(checkContainingTriangle(pointContainedLits[i].under_containing_triangle,point,p_under,pointContainedLits[i].inclusive) || p_under.contains(point,&pointContainedLits[i].under_containing_triangle,pointContainedLits[i].inclusive)){
				//l is true
				if(solver->value(l)==l_True){
					//do nothing
				}else if(solver->value(l)==l_Undef){

					solver->enqueue(l,point_contained_marker) ;
				}else if (solver->value(l)==l_False){
					conflict.push(l);
					buildPointContainedReason(point,conflict,pointContainedLits[i].inclusive);
					return false;
				}
			}else if (!checkContainingTriangle(pointContainedLits[i].over_containing_triangle,point,p_over, pointContainedLits[i].inclusive) && !p_over.contains(point,&pointContainedLits[i].over_containing_triangle,pointContainedLits[i].inclusive)){
			//}else if (!p_over.contains(point)){
				l=~l;
				//l is true
				if(solver->value(l)==l_True){
					//do nothing
				}else if(solver->value(l)==l_Undef){
					solver->enqueue(l,point_not_contained_marker);
				}else if (solver->value(l)==l_False){
					conflict.push(l);
					buildPointNotContainedReason(point,conflict,pointContainedLits[i].inclusive);
					return false;
				}
			}
		}

		for(int i =0;i<convexIntersectionLits.size();i++){
			ConvexPolygon<D,T> & polygon = convexIntersectionLits[i].polygon;
			Lit l = convexIntersectionLits[i].l;

			if(p_under.intersects(polygon,convexIntersectionLits[i].inclusive)){
				//l is true
				if(solver->value(l)==l_True){
					//do nothing
				}else if(solver->value(l)==l_Undef){

					solver->enqueue(l,convex_intersection_marker) ;
				}else if (solver->value(l)==l_False){
					conflict.push(l);
					buildIntersectsReason(polygon,conflict,convexIntersectionLits[i].inclusive);
					return false;
				}
			}else if (!p_over.intersects(polygon,convexIntersectionLits[i].inclusive)){
				l=~l;
				//l is true
				if(solver->value(l)==l_True){
					//do nothing
				}else if(solver->value(l)==l_Undef){
					solver->enqueue(l,convex_not_intersection_marker) ;
				}else if (solver->value(l)==l_False){
					conflict.push(l);
					buildNotIntersectsReason(polygon,conflict,convexIntersectionLits[i].inclusive);
					return false;
				}
			}
		}

		return true;
	}

	void buildReason(Lit p, vec<Lit> & reason, CRef marker){

	}

	virtual void buildAreaGTReason(T area, vec<Lit> & conflict)=0;
	virtual void buildAreaLEQReason(T area,vec<Lit> & conflict)=0;
	virtual void buildPointContainedReason(const Point<D,T> & s, vec<Lit> & conflict, bool inclusive)=0;
	virtual void buildPointNotContainedReason(const Point<D,T> & s, vec<Lit> & conflict, bool inclusive)=0;
	virtual void buildIntersectsReason(Shape<D,T> & p,vec<Lit> & conflict, bool inclusive)=0;
	virtual void buildNotIntersectsReason(Shape<D,T> & p,vec<Lit> & conflict, bool inclusive)=0;

	Shape<D,T> & findSeparatingRegion(Polygon<D,T> & p1, Polygon<D,T> & p2, bool inclusive, bool favourP2=false){
		assert(!p1.intersects(p2,false));
		if(D==2){
			Shape<D,T>&   s = (Shape<D,T>& ) findSeparatingRegion2D((Polygon<2,T>& )p1,(Polygon<2,T>& )p2,inclusive,favourP2);
			assert(s.intersects(p1,true));
			assert(!s.intersects(p2,inclusive));
			return s;
		}
		assert(false);
	}
private:

	//Find a separating region that completely contains p1, and doesn't contain p2.
	Shape<2,T> & findSeparatingRegion2d(Polygon<2,T> & p1, Polygon<2,T> & p2, bool inclusive, bool favourP2=false){
		//first, attempt to find a separating plane.
		if(findSeparatingHalfPlane2d(p1,p2,separating_halfplane,inclusive,favourP2)){
			return separating_halfplane;
		}
		//try to find a 'low complexity' separating region that divides these in half.

		//strategy: decompose both polygons into convex polygons. Find separating planes between each pair.
		//Then find all the intersection points between these planes. The separating region will be selected from among those points.

	/*	static std::vector<NPolygon<2,T>> decomposition1;
		p1.convexDecomposition(decomposition1);

		static std::vector<NPolygon<2,T>> decomposition2;
		p2.convexDecomposition(decomposition2);

		//now, for all pairs (p1,p2) find separating planes

		static std::vector<HalfPlane<2,T>> planes;
		//continue to add separating half planes until p1 and p2 are fully separated.
		//how do you add separating planes?
		for (auto & convex1:decomposition1){
			assert(convex1.isConvex());
			for (auto & convex2:decomposition2){
				assert(convex2.isConvex());
				planes.push_back();
				findSeparatingHalfPlane2d(p1,p2,planes.back(),inclusive,favourP2);

			}
		}*/

		//if nothing else works, return either p1 or p2. These are always valid separating regions.
		if(favourP2){
			return p2.inverse() ;
		}else{
			return p1;
		}

	}
	//Only guaranteed to succeed if p1 and p2 are both convex; MAY succeed otherwise.
	//returns a half plane that is between the two polygons; an effort will be made to place it either in the middle, or closest to p2
	bool findSeparatingHalfPlane2D(Polygon<2,T> & p1, Polygon<2,T> & p2, HalfPlane<2,T> & store, bool inclusive, bool favourP2=false){

			if(p1.size()==0 || p2.size()==0){
				return false;
			}else if (p1.size()==1 && p2.size()==1){
				 Point<2,T> un_normalized_normal = p2[0]-p1[0];
				 store.normal = un_normalized_normal;
				 if(favourP2){
					 //place the half plane as close as possible to p2.
					 store.distance = store.normal.dot(p2[0]);
				 }else{
					 //place the half plane in between p1 and p2
					 store.distance = (store.normal.dot(p2[0])+store.normal.dot(p1[0]))/2;
				 }
				return true;
			}
			for (int q = 0;q<2;q++){
				Polygon<2, T> & h1 = q==0 ? p1 : p2; //(p1.size()<=p2.size())?p1:p2;
				Polygon<2, T> & h2 = q==1 ?  p1 : p2; //(p1.size()<=p2.size())?p2:p1;
					 //Separating Axis Theorem for collision detection between two convex polygons
					 //loop through each edge in _each_ polygon and project both polygons onto that edge's normal.
					 //If any of the projections are non-intersection, then these don't collide; else, they do collide
					 if(h1.size()>1){
						for(int i = 0;i<h1.size();i++){
						 auto & p = h1[i];
						 auto & prev = h1[i-1];
						 Point<2,T> edge = p-prev;
						 Point<2,T> un_normalized_normal(-edge.y, edge.x);

						 //now project both polygons onto to this normal and see if they overlap, by finding the minimum and maximum distances
						 //Note that since we are NOT normalizing the normal vector, the projection is distorted along that vector
						 //(this still allows us to check overlaps, but means that the minimum distance found between the two shapes may be incorrect)
						 T left = numeric<T>::infinity();
						 T right = -numeric<T>::infinity();
						 for (auto & p:h1){
							 T projection = un_normalized_normal.dot(p);
							 if (projection < left) {
								  left = projection;
							 }
							if (projection > right) {
								  right = projection;
							 }
						 }

						 bool seenLeft = false;
						 bool seenRight=false;
						 bool overlaps = false;

						 T left_h2 = numeric<T>::infinity();
						 T right_h2 = -numeric<T>::infinity();

						 for (auto & p:h2){
								 T projection = un_normalized_normal.dot(p);
								 if(inclusive){
									 if (projection >= left && projection <= right ) {
										 seenRight=true;
										 seenLeft=true;
										 break;
									 }else if (projection < left ){
										 seenLeft=true;
										 if(seenRight){
											 break;
										 }
									 }else if (projection>right){
										 seenRight=true;
										 if (seenLeft){
											 break;
										 }
									 }else if (seenLeft && projection > left ){
										 seenRight=true;
										 break;
									 }else if (seenRight && projection < right ){
										 seenRight=true;
										 break;
									 }
							 }else{
								 if(projection>left){
									 seenRight=true;
									 if (seenLeft){
										 break;
									 }
								 }
								 if (projection<right){
									 seenLeft=true;
									 if(seenRight){
										 break;
									 }
								 }
							 }

								 if (projection < left_h2) {
									 left_h2 = projection;
								 }
								if (projection > right_h2) {
									right_h2 = projection;
								 }

						 }
						 if(!(seenLeft&&seenRight)){
							 store.normal = un_normalized_normal;
							 if(favourP2){
								 if(q==0){
									 if(seenLeft){
										 store.distance = left_h2;
									 }else{
										 store.distance = right_h2;
									 }
								 }else{
									 if(seenLeft){
										 store.distance = left;
									 }else{
										 store.distance = right;
									 }
								 }
							 }else{
								 if(seenLeft){
									 store.distance = (left+left_h2)/2;
								 }else{
									 store.distance = (right+right_h2)/2;
								 }
							 }

							 return true;
						 }
					 }
				}
			}

			 return false;

	}
};

template<unsigned int D, class T>
long SymbolicPolygon<D,T>::stats_bounds_skips_under=0;
template<unsigned int D, class T>
long SymbolicPolygon<D,T>::stats_bounds_skips_over=0;

template<unsigned int D, class T>
long SymbolicPolygon<D,T>::stats_line_intersection_skips_under=0;
template<unsigned int D, class T>
long SymbolicPolygon<D,T>::stats_line_intersection_skips_over=0;

template<unsigned int D, class T>
long SymbolicPolygon<D,T>::stats_propagations=0;
template<unsigned int D, class T>
long SymbolicPolygon<D,T>::stats_bound_checks=0;

template<unsigned int D, class T>
long SymbolicPolygon<D,T>::stats_under_clause_length=0;
template<unsigned int D, class T>
long SymbolicPolygon<D,T>::stats_over_clause_length=0;
template<unsigned int D, class T>
long SymbolicPolygon<D,T>::stats_under_clauses=0;
template<unsigned int D, class T>
long SymbolicPolygon<D,T>::stats_over_clauses=0;



template<unsigned int D,class T>
class ConstantPolygon:public SymbolicPolygon<D,T>{
	Polygon<D,T>* polygon;

public:

	ConstantPolygon(Polygon<D,T>* polygon,int detectorID,TheorySolver * solver,  double seed=1): SymbolicPolygon<D,T>(detectorID,solver,seed), polygon(polygon){

	}

	Polygon<D,T>& getOverApprox(){
		return *polygon;
	}
	Polygon<D,T>& getUnderApprox(){
		return *polygon;
	}
	void buildAreaGTReason(T area, vec<Lit> & conflict){
		assert(polygon->getArea()>=area);
		//this polygon is constant, so learn nothing.
	}

	bool isConstant(){
		return true;
	}
	void buildAreaLEQReason(T area,vec<Lit> & conflict){
		assert(polygon->getArea()<area);
		//this polygon is constant, so learn nothing.
	}

	void buildPointContainedReason(const Point<D,T> & s, vec<Lit> & conflict, bool inclusive){
		assert(polygon->contains(s,inclusive));
		//this polygon is constant, so learn nothing.

	}


	void buildPointNotContainedReason(const Point<D,T> & s, vec<Lit> & conflict, bool inclusive){
		assert(!polygon->contains(s,inclusive));
		//this polygon is constant, so learn nothing.
	}

	void buildIntersectsReason(Shape<D,T> & p,vec<Lit> & conflict, bool inclusive){
		assert(p.intersects(*polygon,inclusive));
		//this polygon is constant, so learn nothing.
	}

	void buildNotIntersectsReason(Shape<D,T> & p,vec<Lit> & conflict, bool inclusive){
		assert(!p.intersects(*polygon,inclusive));
		//this polygon is constant, so learn nothing.
	}


};

/*
enum class BinaryPolygonOperation{
	op_union,op_intersect,op_difference,op_minkowski_sum
} ;*/

template<unsigned int D,class T>
class BinaryOperationPolygon:public SymbolicPolygon<D,T>{

	SymbolicPolygon<D,T> & a;
	SymbolicPolygon<D,T> & b;

	PolygonOperationType operation;
	NPolygon<D,T> over_approx;
	NPolygon<D,T> under_approx;
	NPolygon<D,T> polygon_tmp;
	NPolygon<D,T> polygon_tmp2;
	HalfPlane<D,T> tmp_plane;
	std::vector<HalfPlane<D,T>> tmp_planes;
	Polygon<D,T> & applyOperation(Polygon<D,T> & p_a, Polygon<D,T> & p_b, NPolygon<D,T> * store){
		switch (operation){
			case PolygonOperationType::op_union:{
				return *p_a.binary_union(p_b,store);
			}break;
			case PolygonOperationType::op_intersect:{
				return *p_a.binary_intersect(p_b,store);
			}break;
			case PolygonOperationType::op_difference:{
				return *p_a.binary_difference(p_b,store);
			}break;
			case PolygonOperationType::op_minkowski_sum:{
				return *p_a.binary_minkowski_sum(p_b,store);
			}break;
			default:
				assert(false);
				break;
		}
	}

public:

	BinaryOperationPolygon(SymbolicPolygon<D,T> & a, SymbolicPolygon<D,T> & b,PolygonOperationType operation,int detectorID,TheorySolver * solver,  double seed=1): SymbolicPolygon<D,T>(detectorID,solver,seed),  a(a),b(b),operation(operation){

	}
	bool isConstant(){
		return a.isConstant() && b.isConstant();
	}
	Polygon<D,T> & getOverApprox(){
		if(operation==PolygonOperationType::op_difference){
			return *a.getOverApprox().binary_difference(b.getUnderApprox(),&over_approx);
		}else
			return applyOperation(a.getOverApprox(),b.getOverApprox(),&over_approx);
	}
	Polygon<D,T> & getUnderApprox(){
		if(operation==PolygonOperationType::op_difference){
			return *a.getUnderApprox().binary_difference(b.getOverApprox(),&under_approx);
		}else
			return applyOperation(a.getUnderApprox(),b.getUnderApprox(),&under_approx);
	}

	void buildAreaGTReason(T area, vec<Lit> & conflict){
		switch (operation){
			case PolygonOperationType::op_union:{
				//either need to get larger
				a.buildAreaGTReason(0,conflict);
				b.buildAreaGTReason(0,conflict);
			}break;
			case PolygonOperationType::op_intersect:{
				//both need to get larger
				a.buildAreaGTReason(0,conflict);
				b.buildAreaGTReason(0,conflict);
			}break;
			case PolygonOperationType::op_difference:{
				//T intersect_area = a.getOverApprox().binary_intersect(b.getUnderApprox()).getArea();
				//a has to get larger, or b has to get smaller.
				a.buildAreaGTReason(0,conflict);
				b.buildAreaLEQReason(0,conflict);
			}break;
				case PolygonOperationType::op_minkowski_sum:{
					//either need to get larger
					a.buildAreaGTReason(0,conflict);
					b.buildAreaGTReason(0,conflict);
			}break;
			default:
				assert(false);
				break;
		}
	}
	void buildAreaLEQReason(T area,vec<Lit> & conflict){
		switch (operation){
			case PolygonOperationType::op_union:{
				//either need to get larger
				a.buildAreaLEQReason(0,conflict);
				b.buildAreaLEQReason(0,conflict);
			}break;
			case PolygonOperationType::op_intersect:{
				//both need to get larger
				a.buildAreaLEQReason(0,conflict);
				b.buildAreaLEQReason(0,conflict);
			}break;
			case PolygonOperationType::op_difference:{
				//T intersect_area = a.getOverApprox().binary_intersect(b.getUnderApprox()).getArea();
				//a has to get larger, or b has to get smaller.
				a.buildAreaLEQReason(0,conflict);
				b.buildAreaGTReason(0,conflict);
			}break;
				case PolygonOperationType::op_minkowski_sum:{
					//either need to get larger
					a.buildAreaLEQReason(0,conflict);
					b.buildAreaLEQReason(0,conflict);
			}break;
			default:
				assert(false);
				break;
		}
	}

	void buildIntersectsReason(Shape<D,T> & hp,vec<Lit> & conflict, bool inclusive){
		assert(getUnderApprox().intersects(hp,inclusive));
		switch (operation){
			case PolygonOperationType::op_union:{
				//either a or b (or both) intersects hp. pick one, and explain it.
				if(a.getUnderApprox().intersects(hp,inclusive)){
					a.buildIntersectsReason(hp, conflict,inclusive);
				}else{
					assert(b.getUnderApprox().intersects(hp));
					b.buildIntersectsReason(hp, conflict,inclusive);
				}
			}break;
			case PolygonOperationType::op_intersect:{
				//both a and b must both intersect hp, so in order to stop them from intersecting, one of them must change.
				assert(a.getUnderApprox().intersects(hp,inclusive));
				assert(b.getUnderApprox().intersects(hp,inclusive));
				a.buildIntersectsReason(hp, conflict,inclusive);
				b.buildIntersectsReason(hp, conflict,inclusive);
			}break;
			case PolygonOperationType::op_difference:{
				//a must intersect hp, and b must NOT intersect hp, so in order to stop them from intersecting, one of them must change.
				assert(a.getUnderApprox().intersects(hp,inclusive) && !b.getOverApprox().intersects(hp,inclusive) );
				a.buildIntersectsReason(hp, conflict,inclusive);
				b.buildNotIntersectsReason(hp, conflict,inclusive);
			}break;
				case PolygonOperationType::op_minkowski_sum:{
				if(hp.getType()==HALF_PLANE){
					buildMinkowskiHalfPlaneIntersectReason((HalfPlane<D,T>&)hp,conflict,inclusive);
				}else if(hp.getType()==LINE){
					assert(false);
					fprintf(stderr,"unimplemented collision learning\n");
					exit(3);
				}else{
					buildMinkowskiPolygonIntersectReason((Polygon<D,T>&)hp,conflict,inclusive);
				}/*else{
					assert(false);
					fprintf(stderr,"unimplemented collision learning\n");
					exit(3);
				}*/
			}break;
			default:
				assert(false);
				break;
		}
	}

	void buildNotIntersectsReason(Shape<D,T> & hp,vec<Lit> & conflict, bool inclusive){
		assert(!getOverApprox().intersects(hp,inclusive));
		switch (operation){
			case PolygonOperationType::op_union:{
				//neither a nor b intersects hp.
				assert(!a.getOverApprox().intersects(hp,inclusive));
				assert(!b.getOverApprox().intersects(hp,inclusive));
				a.buildNotIntersectsReason(hp, conflict,inclusive);
				b.buildNotIntersectsReason(hp, conflict,inclusive);
			}break;
			case PolygonOperationType::op_intersect:{
				//either a or b doesn't intersect hp. pick one that doesn't intersect and explain it
				assert(!a.getOverApprox().intersects(hp,inclusive)||b.getOverApprox().intersects(hp,inclusive));

				if(!a.getOverApprox().intersects(hp,inclusive)){
					a.buildNotIntersectsReason(hp, conflict,inclusive);
				}else{
					assert(b.getOverApprox().intersects(hp,inclusive));
					b.buildNotIntersectsReason(hp, conflict,inclusive);
				}
			}break;
			case PolygonOperationType::op_difference:{
				//if A-B doesn't intersect HP, then either
				//A doesnt intersect const(HP-B),
				//or B DOES intersect (HP intersect A)

				if(!a.getOverApprox().intersects(hp,inclusive)){
					//then regardless of b, a failing to intersect hp is a sufficient explanation
					a.buildNotIntersectsReason(hp, conflict,inclusive);
				}else{
					//then it must be the case that B intersects HP. (this analysis can be made more precise though...)
					b.buildIntersectsReason(hp, conflict,inclusive);
				}
			}break;
			case PolygonOperationType::op_minkowski_sum:{
				if(hp.getType()==HALF_PLANE){
					buildMinkowskiHalfPlaneNotIntersectReason((HalfPlane<D,T>&)hp,conflict,inclusive);
				}else if(hp.getType()==LINE){
					assert(false);
					fprintf(stderr,"unimplemented collision learning\n");
					exit(3);
				}else{
					buildMinkowskiPolygonNotIntersectReason((Polygon<D,T>&)hp,conflict,inclusive);
				}/*else{
					assert(false);
					fprintf(stderr,"unimplemented collision learning\n");
					exit(3);
				}*/
			}break;
			default:
				assert(false);
				break;
		}
	}

	void buildPointContainedReason(const Point<D,T> & s, vec<Lit> & conflict, bool inclusive){
		assert(getUnderApprox().contains(s));
		switch (operation){
			case PolygonOperationType::op_union:{
				//either a or b (or both) contains s.
				if(a.getUnderApprox().contains(s,inclusive)){
					a.buildPointContainedReason(s, conflict,inclusive);
				}else{
					assert(b.getUnderApprox()->contains(s,inclusive));
					b.buildPointContainedReason(s, conflict,inclusive);
				}
			}break;
			case PolygonOperationType::op_intersect:{
				//both a and b must both contain s
				assert(a.getUnderApprox().contains(s,inclusive));
				assert(b.getUnderApprox().contains(s,inclusive));
				a.buildPointContainedReason(s, conflict,inclusive);
				b.buildPointContainedReason(s, conflict,inclusive);
			}break;
			case PolygonOperationType::op_difference:{
				//a must contain s, and b must NOT contain s
				assert(a.getUnderApprox().contains(s,inclusive) && !b.getOverApprox().contains(s,inclusive) );
				a.buildPointContainedReason(s, conflict,inclusive);
				b.buildPointNotContainedReason(s, conflict,inclusive);
			}break;
			case PolygonOperationType::op_minkowski_sum:{
				//the reason that A + B contains point s is that
				//(A-s) intersects B
				//or, equivalently, (B-s) intersects A

				Polygon<D,T> &translatedB = *b.getUnderApprox().translate(s,&polygon_tmp);
				assert(translatedB.intersects(*a.getUnderApprox(),inclusive));
				a.buildIntersectsReason(translatedB,conflict,inclusive);
				Polygon<D,T> &translatedA = *a.getUnderApprox().translate(s,&polygon_tmp);
				assert(translatedA.intersects(*b.getUnderApprox(),inclusive));
				b.buildIntersectsReason(translatedA,conflict,inclusive);

			}break;
			default:
				assert(false);
				break;
		}
	}

	void buildPointNotContainedReason(const Point<D,T> & s, vec<Lit> & conflict, bool inclusive){
		assert(!getOverApprox().contains(s));
		switch (operation){
			case PolygonOperationType::op_union:{
				//neither a nor b intersects hp.
				assert(!a.getOverApprox().contains(s,inclusive));
				assert(!b.getOverApprox().contains(s,inclusive));
				a.buildPointNotContainedReason(s, conflict,inclusive);
				b.buildPointNotContainedReason(s, conflict,inclusive);
			}break;
			case PolygonOperationType::op_intersect:{
				//either a or b doesn't intersect hp. pick one that doesn't intersect and explain it
				assert(!a.getOverApprox().contains(s,inclusive)||b.getOverApprox().contains(s,inclusive));

				if(!a.getOverApprox().contains(s,inclusive)){
					a.buildPointNotContainedReason(s, conflict,inclusive);
				}else{
					assert(b.getOverApprox().contains(s,inclusive));
					b.buildPointNotContainedReason(s, conflict,inclusive);
				}
			}break;
			case PolygonOperationType::op_difference:{
				//if A-B doesn't intersect HP, then either
				//A doesnt intersect const(HP-B),
				//or B DOES intersect (HP intersect A)

				if(!a.getOverApprox().contains(s,inclusive)){
					//then regardless of b, a failing to intersect hp is a sufficient explanation
					a.buildPointNotContainedReason(s, conflict,inclusive);
				}else{
					//then it must be the case that B intersects HP. (this analysis can be made more precise though...)
					b.buildPointContainedReason(s, conflict,inclusive);
				}
			}break;
			case PolygonOperationType::op_minkowski_sum:{
				//the reason that A + B does not contain point s is that
				//(A-s) does not intersect B
				Polygon<D,T> &translatedB = *b.getOverApprox().translate(s,&polygon_tmp);
				assert(!translatedB.intersects(*a.getOverApprox(),inclusive));
				a.buildNotIntersectsReason(translatedB,conflict,inclusive);
				Polygon<D,T> &translatedA = *a.getOverApprox().translate(s,&polygon_tmp);
				assert(!translatedA.intersects(*b.getOverApprox(),inclusive));
				b.buildNotIntersectsReason(translatedA,conflict,inclusive);
			}break;
			default:
				assert(false);
				break;
		}
	}

	void buildMinkowskiPolygonIntersectReason(Polygon<D,T> & hp, vec<Lit> & conflict, bool inclusive){
		assert(getUnderApprox().intersects(hp));
		assert(operation== PolygonOperationType::op_minkowski_sum);

		//Two polygons a, b intersect if the minkowski sum of a and the inverse of b contains the (0,0):
		// a + (-b) contains 0
		//The minkowski sum of a and b is associative.
		//so ((a+b) + -c) contains 0 <==> (a + (b+ -c)) contains 0
		Polygon<D,T> & inverse = *hp.inverse(&polygon_tmp2);
		Polygon<D,T> & addition = *a.getOverApprox().binary_minkowski_sum(inverse,  &polygon_tmp);
		assert(addition.intersects(b.getOverApprox()));
		b.buildIntersectsReason(addition,conflict,inclusive);

		addition = *b.getOverApprox().binary_minkowski_sum(inverse,  &polygon_tmp);
		assert(addition.intersects(a.getOverApprox()));
		a.buildIntersectsReason(addition,conflict,inclusive);
	}
	void buildMinkowskiPolygonNotIntersectReason(Polygon<D,T> & h2, vec<Lit> & conflict, bool inclusive){
		assert(!getOverApprox().intersects(h2,inclusive));
		assert(operation==PolygonOperationType::op_minkowski_sum);
		Polygon<D,T> & inverse = *h2.inverse(&polygon_tmp2);
		Polygon<D,T> & addition = *a.getOverApprox().binary_minkowski_sum(inverse,  &polygon_tmp);
		assert(!addition.intersects(b.getOverApprox()));
		b.buildNotIntersectsReason(addition,conflict,inclusive);

		addition = *b.getOverApprox().binary_minkowski_sum(inverse,  &polygon_tmp);
		assert(!addition.intersects(a.getOverApprox()));
		a.buildNotIntersectsReason(addition,conflict,inclusive);
	/*	//first, attempt to find a separating axis. This is not guaranteed to exist unless both polygons are convex
		Polygon<D,T> & h1 = getOverApprox();
		 //Separating Axis Theorem for collision detection between two convex polygons
		 //loop through each edge in _each_ polygon and project both polygons onto that edge's normal.
		 //If any of the projections are non-intersection, then these don't collide; else, they do collide
		if(h2.size()==0){
			//no intersection is possible
			return;
		}else if (h2.size()==1){
			if(!inclusive){
				return;
			}else{
				buildPointNotContainedReason(h2[0],conflict,inclusive);
				return;
			}
		}else if(h1.size()==0){
			//not sure how to handle best this.
			buildAreaLEQReason(0,conflict);
			return;
		}

		 for(int i = 0;i<h1.size();i++){
			 auto & p = h1[i];
			 auto & prev = (h1)[i-1];
			 Point<2,T> edge = p-prev;
			 Point<2,T> un_normalized_normal(-edge.y, edge.x);

			 //now project both polygons onto to this normal and see if they overlap, by finding the minimum and maximum distances
			 //Note that since we are NOT normalizing the normal vector, the projection is distorted along that vector
			 //(this still allows us to check overlaps, but means that the minimum distance found between the two shapes may be incorrect)
			 T left = numeric<T>::infinity();
			 T right = -numeric<T>::infinity();
			 for (auto & p:h1){
				 T projection = un_normalized_normal.dot(p);
				 //projection_out1.push_back({p,projection});
				 if (projection < left) {
					  left = projection;
				 }
				 if (projection > right) {
					  right = projection;
				 }
			 }
			 bool seenLeft = false;
			 bool seenRight=false;
			 T h2_min = numeric<T>::infinity();
			 T h2_max = -numeric<T>::infinity();

			 for (auto & p:h2){
				 T projection = un_normalized_normal.dot(p);
				 //projection_out2.push_back({p,projection});
				 if(projection< h2_min){
					 h2_min=projection;
				 }
				 if(projection> h2_max){
					 h2_max=projection;
				 }
				 if(inclusive){
				 if (projection >= left && projection <= right ) {
					 seenRight=true;
					 seenLeft=true;
					 break;
				 }else if (projection < left ){
					 seenLeft=true;
					 if(seenRight){
						 break;
					 }
				 }else if (projection>right){
					 seenRight=true;
					 if (seenLeft){
						 break;
					 }
				 }else if (seenLeft && projection > left ){
					 seenRight=true;
					 break;
				 }else if (seenRight && projection < right ){
					 seenRight=true;
					 break;
				 }
			 }else{
				 if(projection>left){
					 seenRight=true;
					 if (seenLeft){
						 break;
					 }
				 }
				 if (projection<right){
					 seenLeft=true;
					 if(seenRight){
						 break;
					 }
				 }
			 }
			 }

			 if(!(seenLeft&&seenRight)){
				 //this is a separating axis
				 HalfPlane<D,T> & hp = tmp_plane;
				 hp.normal = un_normalized_normal;
				 assert(seenLeft||seenRight);
				 if(seen_left){
					 hp.distance = h2_max;
				 }else{
					 hp.distance = h2_min;
				 }

				 this->buildIntersectsReason(hp,conflict,inclusive);
				 return;
			 }

		 }

		 //now test the axis produced by the other polygon
		 for(int i = 0;i<h2.size();i++){
			 auto & p = (h2)[i];
			 auto & prev = (h2)[i-1];
			 Point<2,T> edge = p-prev;
			 Point<2,T> un_normalized_normal(-edge.y, edge.x);


			 T left =numeric<T>::infinity();
			 T right = -numeric<T>::infinity();
			 for (auto & p:h2){
				 T projection = un_normalized_normal.dot(p);

				 if (projection < left) {
					  left = projection;
				 }
				 if (projection > right) {
					  right = projection;
				 }
			 }
			 T h1_min = numeric<T>::infinity();
			 T h1_max = -numeric<T>::infinity();
			 bool seenLeft = false;
			 bool seenRight=false;
			 bool overlaps = false;
			 for (auto & p:h1){
				 T projection = un_normalized_normal.dot(p);
				 if(inclusive){
				 if (projection >= left && projection <= right ) {
					 seenRight=true;
					 seenLeft=true;
					 break;
				 }else if (projection < left ){
					 seenLeft=true;
					 if(seenRight){
						 break;
					 }
				 }else if (projection>right){
					 seenRight=true;
					 if (seenLeft){
						 break;
					 }
				 }else if (seenLeft && projection > left ){
					 seenRight=true;
					 break;
				 }else if (seenRight && projection < right ){
					 seenRight=true;
					 break;
				 }
			 }else{
				 if(projection>left){
					 seenRight=true;
					 if (seenLeft){
						 break;
					 }
				 }
				 if (projection<right){
					 seenLeft=true;
					 if(seenRight){
						 break;
					 }
				 }

			 }
			 }
			 if(!(seenLeft&&seenRight)){
				 //this is a separating axis
				 HalfPlane<D,T> & hp = tmp_plane;
				 hp.normal = un_normalized_normal;
				 assert(seenLeft||seenRight);
				 if(seen_left){
					 hp.distance = right;
				 }else{
					 hp.distance = left;
				 }
				 this->buildIntersectsReason(hp,conflict,inclusive);
				 return;
			 }
		 }

		assert(!this->getOverApprox().isConvex() || !h2.isConvex());
*/
		//or is the following correct:
		//a + b not intersects h2 <--> a not intersects h2+(inverted b), and b not intersects h2+(inverted a)


		//if either or both shape is not convex, then I will handle that by finding a set of half planes that contain h2 and not h1
		//then for each halfplane_i in that set, form  (this intersect ((union half planes) - halfplane_i)), build the reason for that partition not to intersect b.
		/*if(tmp_planes.size()<h2.size()-1){
			tmp_planes.resize(h2.size()-1);
		}
		//build half planes around each edge
		for(int i = 0;i<h2.size();i++){
			auto & p0 = h2[i-1];
			auto & p1 = h2[i];
			Point<2,T> edge = p1-p0;
			Point<2,T> un_normalized_normal(-edge.y, edge.x);
			tmp_planes[i].normal.x = -edge.y;
			tmp_planes[i].normal.x = edge.x;
			tmp_planes[i].distance=p0.dot(un_normalized_normal);
			assert(!tmp_planes[i].intersetcs(h1,inclusive));
			assert(tmp_planes[i].intersetcs(h2,true));
		}
*/


	}

	void buildMinkowskiHalfPlaneIntersectReason(HalfPlane<D,T> & hp, vec<Lit> & conflict, bool inclusive){
		assert(getUnderApprox().intersects(hp,inclusive));
		assert(operation==PolygonOperationType::op_minkowski_sum);

		//one option is to convert this to point containment:
		//there exists a point in HP that intersects the minkowski sum of a and b.
		//analyze that point containment.

		//instead, for now, stick with half plane intersection
		//compute two new, translated half planes with the same normal: one for a, one for b. at least one of them must intersect.

		//is this correct??

		T projection = hp.normal.dot(a.getUnderApprox()[0]);
		T distance = hp.distance-projection;
		HalfPlane<D,T> translated(hp.normal,distance);
		//project an arbitrary point from A into the normal
		if(translated.intersects(b.getUnderApprox(),inclusive)){
			b.buildIntersectsReason(hp, conflict,inclusive);
		}else{
			T projection = hp.normal.dot(b.getUnderApprox()[0]);
			translated.distance = hp.distance-projection;
			assert(translated.intersects(a.getUnderApprox(),inclusive));
			//remove this check later...
			if(!translated.intersects(a.getUnderApprox(),inclusive)){
				fprintf(stderr,"Intersection error!\n");
				exit(3);
			}
			a.buildIntersectsReason(hp, conflict,inclusive);
		}
	}


	void buildMinkowskiHalfPlaneNotIntersectReason(HalfPlane<D,T> & hp, vec<Lit> & conflict, bool inclusive){
		assert(!getOverApprox().intersects(hp,inclusive));
		assert(operation==PolygonOperationType::op_minkowski_sum);

		//project all the points of a, b onto the hp axis. Find the farthest ones from hp
		//create two new half planes, with the same normal but translated distances.
		Polygon<D,T> & a_under = a.getUnderApprox();
		T max_projection =  -numeric<T>::inf;
		for(auto & p:a_under){
			T projection = hp.normal.dot(p);
			if(projection>max_projection){
				max_projection=projection;
			}
		}
		T distance = hp.distance-max_projection;
		HalfPlane<D,T> translated(hp.normal,distance);
		assert(!translated.intersects(b.getUnderApprox(),inclusive));
		b.buildIntersectsReason(translated, conflict,inclusive);

		Polygon<D,T> & b_under = b.getUnderApprox();
		max_projection =  -numeric<T>::inf;
		for(auto & p:b_under){
			T projection = hp.normal.dot(p);
			if(projection>max_projection){
				max_projection=projection;
			}
		}
		distance = hp.distance-max_projection;
		translated.distance =distance;
		assert(!translated.intersects(a.getUnderApprox(),inclusive));
		a.buildNotIntersectsReason(translated, conflict,inclusive);
	}
};

template<unsigned int D,class T>
class ConditionalPolygon:public SymbolicPolygon<D,T>{
	SymbolicPolygon<D,T> & thn;
	SymbolicPolygon<D,T> & els;
	Lit condition;

	NPolygon<D,T> over_approx;
	NPolygon<D,T> under_approx;
public:

	ConditionalPolygon(Lit condition, SymbolicPolygon<D,T> & thn, SymbolicPolygon<D,T> & els,int detectorID,TheorySolver * solver,  double seed=1): SymbolicPolygon<D,T>(detectorID,solver,seed), thn(thn),els(els){
		condition= mkLit( this->solver->newVar(var(condition),this->getID()),sign(condition));
	}

	bool isConstant(){
		if (this->solver->isConstant(condition)){
			if(this->solver->value(condition)==l_True){
				return thn.isConstant();
			}else{
				return els.isConstant();
			}
		}
		return false;
	}

	Polygon<D,T>& getOverApprox(){
		lbool val = this->solver->value(condition);
		if(val==l_True){
			return thn.getOverApprox();
		}else if (val==l_False){
			return els.getOverApprox();
		}
		//return the union of the two polygons
		return *thn.getOverApprox().binary_union(els.getOverApprox(),&over_approx);
	}
	Polygon<D,T>& getUnderApprox(){
		lbool val = this->solver->value(condition);
		if(val==l_True){
			return thn.getUnderApprox();
		}else if (val==l_False){
			return els.getUnderApprox();
		}
		//return the intersection of the two polygons
		return *thn.getUnderApprox().binary_intersect(els.getUnderApprox(),&under_approx);
	}

	void buildAreaGTReason(T area, vec<Lit> & conflict){

	}
	void buildAreaLEQReason(T area,vec<Lit> & conflict){

	}
	void buildPointContainedReason(const Point<D,T> & s, vec<Lit> & conflict, bool inclusive){
		assert(getUnderApprox().contains(s,inclusive));
		lbool val = this->solver->value(condition);
		if(val==l_True){
			thn.buildPointContainedReason(s,conflict,inclusive);
			if(opt_minimize_collision_clauses && !els.getUnderApprox().contains(s,inclusive)){
				els.buildPointContainedReason(s,conflict,inclusive);
			}else{
				conflict.push(~condition);
			}
		}else if (val==l_False){
			els.buildPointContainedReason(s,conflict,inclusive);
			if(opt_minimize_collision_clauses  && !thn.getUnderApprox().contains(s,inclusive)){
				thn.buildPointContainedReason(s,conflict,inclusive);
			}else{
				conflict.push(condition);
			}
		}else{
			assert(val==l_Undef);
			//either a or b (or both) intersects hp. pick one, and explain it.
			if(thn.getUnderApprox().contains(s)){
				thn.buildPointContainedReason(s, conflict,inclusive);
			}else{
				assert(els.getUnderApprox().contains(s));
				els.buildPointContainedReason(s, conflict,inclusive);
			}
		}
	}

	void buildPointNotContainedReason(const Point<D,T> & hp, vec<Lit> & conflict, bool inclusive){
		assert(!getOverApprox().contains(hp,inclusive));
		lbool val = this->solver->value(condition);
		if(val==l_True){
			thn.buildPointNotContainedReason(hp,conflict,inclusive);
			if(opt_minimize_collision_clauses && !els.getOverApprox().contains(hp,inclusive)){
				els.buildPointNotContainedReason(hp,conflict,inclusive);
			}else{
				conflict.push(~condition);
			}
		}else if (val==l_False){
			els.buildPointNotContainedReason(hp,conflict,inclusive);
			if(opt_minimize_collision_clauses && !thn.getOverApprox().contains(hp,inclusive)){
				thn.buildPointNotContainedReason(hp,conflict,inclusive);
			}else{
				conflict.push(condition);
			}
		}else{
			assert(val==l_Undef);
			//either thn or els doesn't intersect hp. pick one that doesn't intersect and explain it
			assert(!thn.getOverApprox().contains(hp)||els.getOverApprox().contains(hp));

			if(!thn.getOverApprox().contains(hp)){
				thn.buildPointNotContainedReason(hp, conflict,inclusive);
			}else{
				assert(els.getOverApprox().contains(hp));
				els.buildPointNotContainedReason(hp, conflict,inclusive);
			}
		}
	}

	void buildIntersectsReason(Shape<D,T> & hp,vec<Lit> & conflict, bool inclusive){
		assert(getUnderApprox().intersects(hp,inclusive));
		lbool val = this->solver->value(condition);
		if(val==l_True){
			thn.buildIntersectsReason(hp,conflict,inclusive);
			if(opt_minimize_collision_clauses && !els.getUnderApprox().intersects(hp,inclusive)){
				els.buildIntersectsReason(hp,conflict,inclusive);
			}else{
				conflict.push(~condition);
			}
		}else if (val==l_False){
			els.buildIntersectsReason(hp,conflict,inclusive);
			if(opt_minimize_collision_clauses && !thn.getUnderApprox().intersects(hp,inclusive)){
				thn.buildIntersectsReason(hp,conflict,inclusive);
			}else{
				conflict.push(condition);
			}
		}else{
			assert(val==l_Undef);
			//either a or b (or both) intersects hp. pick one, and explain it.
			if(thn.getUnderApprox().intersects(hp)){
				thn.buildIntersectsReason(hp, conflict,inclusive);
			}else{
				assert(els.getUnderApprox().intersects(hp));
				els.buildIntersectsReason(hp, conflict,inclusive);
			}
		}
	}

	void buildNotIntersectsReason(Shape<D,T> & hp,vec<Lit> & conflict, bool inclusive){
		assert(!getOverApprox().intersects(hp,inclusive));
		lbool val = this->solver->value(condition);
		if(val==l_True){
			thn.buildNotIntersectsReason(hp,conflict,inclusive);
			if(opt_minimize_collision_clauses && !els.getOverApprox().intersects(hp,inclusive)){
				els.buildNotIntersectsReason(hp,conflict,inclusive);
			}else{
				conflict.push(~condition);
			}
		}else if (val==l_False){
			els.buildNotIntersectsReason(hp,conflict,inclusive);
			if(opt_minimize_collision_clauses && !thn.getOverApprox().intersects(hp,inclusive)){
				thn.buildNotIntersectsReason(hp,conflict,inclusive);
			}else{
				conflict.push(condition);
			}
		}else{
			assert(val==l_Undef);
			//either thn or els doesn't intersect hp. pick one that doesn't intersect and explain it
			assert(!thn.getOverApprox().intersects(hp)||els.getOverApprox().intersects(hp));

			if(!thn.getOverApprox().intersects(hp)){
				thn.buildNotIntersectsReason(hp, conflict,inclusive);
			}else{
				assert(els.getOverApprox().intersects(hp));
				els.buildNotIntersectsReason(hp, conflict,inclusive);
			}
		}
	}

};

template<unsigned int D,class T>
class TranslatedPolygon:public SymbolicPolygon<D,T>{

	SymbolicPolygon<D,T> & p;

	Point<D,T> axis;

	NPolygon<D,T> over_approx;
	NPolygon<D,T> under_approx;
	SymbolicScalar<T> & translation;
	//The idea here (suggested by Jacob Bayless) applies only to convex polygons; you can extrude the polygon infinitely far in the axis of translation in one direction, then separately in the other direction.
	//now you can translate those two around independently, and point containment in either of them is monotonic; you can then assert (outside the theory solver)
	//that those two polygons have the same positions, and take the intersection of them, to recover the original translated polygon

	//we cannot directly apply this to concave polygons - but we don't need to.
	//So long as all concave polygons are the product of a union/difference of convex polygons, we can apply the translation to those convex polygons instead.
	bool isConstant(){
		return translation.isConstant()&& p.isConstant();
	}

};

template<unsigned int D,class T>
class ExtrudedPolygon:public SymbolicPolygon<D,T>{

	SymbolicPolygon<D,T> & polygon;

	Point<D,T> axis;//this is constant

	NPolygon<D,T> over_approx;
	NPolygon<D,T> under_approx;

	SymbolicScalar<T> & extrusion_amount;

	Polygon<D,T> & extrude(Polygon<D,T> & polygon, T & amount,NPolygon<D,T> & store){
		//Given an input polygon, the output is a polygon with twice as many vertices.

		return store;
	}

public:
	//Extrudes a polygon a fixed or variable distance along a fixed axis.
	ExtrudedPolygon(SymbolicPolygon<D,T> & polygon,SymbolicScalar<T> & extrusion_amount,Point<D,T> axis,int detectorID,TheorySolver * solver,  double seed=1): SymbolicPolygon<D,T>(detectorID,solver,seed), polygon(polygon),extrusion_amount(extrusion_amount),axis(axis){

	}

	Polygon<D,T>& getOverApprox(){
		return extrude( polygon.getOverApprox(), extrusion_amount.getOverApprox(),over_approx);
	}
	Polygon<D,T>& getUnderApprox(){
		return extrude( polygon.getUnderApprox(), extrusion_amount.getUnderApprox(),under_approx);
	}
	void buildAreaGTReason(T area, vec<Lit> & conflict){
		assert(polygon->getArea()>=area);

	}
	bool isConstant(){
		return extrusion_amount.isConstant()&& polygon.isConstant();
	}
	void buildAreaLEQReason(T area,vec<Lit> & conflict){
		assert(polygon->getArea()<area);

	}

	void buildPointContainedReason(const Point<D,T> & s, vec<Lit> & conflict, bool inclusive){
		assert(polygon->contains(s,inclusive));

	}
	void buildPointNotContainedReason(const Point<D,T> & s, vec<Lit> & conflict, bool inclusive){

	}
	void buildIntersectsReason(Shape<D,T> & p,vec<Lit> & conflict, bool inclusive){
		//It is not enough to just cast rays from each point of polygon - it is possible for the shape to intersect the extruded shape
		//without intersecting a casted ray.

		//I haven't tested this yet, but my hope is that if I raycast from all the vertices in _both_ polygons, then all intersections will
		//be caught - though this only works for tests against polygons.

		//alternatively, we should be able to find the closest point of intersection on the shape, and cast a ray
		//from that point to somewhere in the polygon, and figure out the extrusion depth that way. That seems like it should always
		//work, but it much more involved.

		T & under_extent = this->extrusion_amount.getUnderApprox();

		if(under_extent==0){
			//this is the original polygon
			polygon.buildIntersectsReason(p,conflict,inclusive);
			return;
		}else{

		}

		if(p.getType()==POLYGON || p.getType()==CONVEX_POLYGON){// || p.getType()==LINE_SEGMENT){
			Polygon<D,T> & p2 = (Polygon<D,T> &)p;

			static LineSegment<D,T> check;
			static Point<D,T> intersection;
			for(auto & p:polygon->getUnderApprox()){
				//cast a ray from this point in the direction of the extrusion, and see if it intersects
				check.a = p;
				check.b = check.a + this->axis*under_extent;
			/*	if(under_extent>0){
					check.b = ray.a+this->axis;
				}else{
					assert(under_extent!=0);
					ray.b = ray.a-this->axis;
				}*/

				if(check.intersects(p,intersection,true)){
					//then this is a sufficient condition to explain the collision.
					//so long as this vertex is contained in the polygon, and the extrusion amount doesn't lower, then this collision will remain.
					//there are many ways to consider optimizing how we split this constraint up to optimize clause learning, below is just the easiest

					T length = check.dot(intersection);
					assert(length<=under_extent);
					if(!inclusive){
						assert(length<under_extent);
					}
					//the reason for the collision is that
					//either the extrusion amount must be less than length
					this->extrusion_amount.buildValueGreaterThanReason(length,conflict,inclusive);
					//OR this vertex must no longer be contained in the polygon
					this->polygon.buildPointContainedReason(p,conflict,true);

					return;
				}
			}

			//if no ray cast from this polygon intersected the other polygon, it is still posible that a ray cast from the other polygon could intersect this one.

			for(auto & p:p2){
				//cast a ray from this point in the direction of the extrusion, and see if it intersects
				check.a = p;
				check.b = check.a - this->axis*under_extent;
			/*	if(under_extent>0){
					check.b = ray.a+this->axis;
				}else{
					assert(under_extent!=0);
					ray.b = ray.a-this->axis;
				}*/

				if(check.intersects(this->polygon,intersection,true)){
					//then this is a sufficient condition to explain the collision.
					//so long as this vertex is contained in the polygon, and the extrusion amount doesn't lower, then this collision will remain.
					//there are many ways to consider optimizing how we split this constraint up to optimize clause learning, below is just the easiest

					T length = check.dot(intersection);
					assert(length<=under_extent);
					if(!inclusive){
						assert(length<under_extent);
					}
					//the reason for the collision is that
					//either the extrusion amount must be less than length
					this->extrusion_amount.buildValueGreaterThanReason(length,conflict,inclusive);
					//OR the intersection point must no longer be contained in the polygon
					this->polygon.buildPointContainedReason(intersection,conflict,true);

					return;
				}
			}

			assert(false);
		}
		assert(false);
	}

	//Analyze why this extrusion amount applied to two CONSTANT, non symbolic shapes, is not sufficient to make them intersect.
	void analyzeExtrusionNotIntersectsReason(Polygon<D,T> & from, Shape<D,T> & to, SymbolicScalar<T> & extrusion_distance, Point<D,T> & axis, vec<Lit> & conflict, bool inclusive){
		if(to.getType()==POLYGON || to.getType()==CONVEX_POLYGON){// || p.getType()==LINE_SEGMENT){
			Polygon<D,T> & p2 = (Polygon<D,T> &)to;
			T min_required_length = numeric<T>::inf;
			bool any_intersect=false;

			static Ray<D,T> check;
			static Point<D,T> intersection;

			for(auto & p:from.getUnderApprox()){
				//cast a ray from this point in the direction of the extrusion, and see if it intersects
				check.a = p;
				check.b = check.a + axis;

				if(check.intersects(p,intersection,true)){
					//then this is a sufficient condition to explain the collision.
					//so long as this vertex is contained in the polygon, and the extrusion amount doesn't lower, then this collision will remain.
					//there are many ways to consider optimizing how we split this constraint up to optimize clause learning, below is just the easiest
					any_intersect=true;
					T length = check.dot(intersection);
					if(length<min_required_length){
						min_required_length=length;
					}
				}
			}

			//if no ray cast from this polygon intersected the other polygon, it is still posible that a ray cast from the other polygon could intersect this one.

			for(auto & p:p2){
				//cast a ray from this point in the direction of the extrusion, and see if it intersects
				check.a = p;
				check.b = check.a - axis;

				if(check.intersects(from,intersection,true)){
					//then this is a sufficient condition to explain the collision.
					//so long as this vertex is contained in the polygon, and the extrusion amount doesn't lower, then this collision will remain.
					//there are many ways to consider optimizing how we split this constraint up to optimize clause learning, below is just the easiest
					any_intersect=true;
					T length = check.dot(intersection);
					if(length<min_required_length){
						min_required_length=length;
					}
				}
			}
			if(any_intersect){
				assert(min_required_length< numeric<T>::inf);
				//the extrusion amount must be _atleast_ this length to cause a collision between these two constant polygons.
				assert(extrusion_distance.getOverApprox()<min_required_length);

				extrusion_distance.buildValueLessThanReason(min_required_length,conflict,inclusive);

			}else{
				assert(min_required_length== numeric<T>::inf);
				//no amount of extrusion along this axis will cause these two to collide
			}
		}
	}

	void buildNotIntersectsReason(Shape<D,T> & p,vec<Lit> & conflict, bool inclusive){
		T & over_extent = this->extrusion_amount.getOverApprox();
		assert(!polygon.getOverApprox().intersects(p,inclusive));
		//one part of the reason for the non-intersection is always that the underlying polygon doesn't intersect ...

		//actually, the correct way to do this, is to convert this into constant shape p, extrueded, doesn't intersect the symbolic polygon.
		//find a separating polygon bewteen these, and then deal with them independently.


		if(polygon.isConstant() ){

			analyzeExtrusionNotIntersectsReason(polygon.getOverApprox(),p,extrusion_amount,this->axis,conflict,inclusive);
		}else if (extrusion_amount.isConstant() || over_extent==0){
			//the over_extent==0 check is only correct is negative extrusions are not allowed.
			polygon.buildNotIntersectsReason(p,conflict,inclusive);
		}else{
			//need to separate these two and test them individually.

			//convert this into an analysis of why constant shape p, extruded along -axis, doesn't intersect the symbolic polygon.
			//find a constant separating polygon between these, and then deal with them independently.



		}
	}
};

#endif


