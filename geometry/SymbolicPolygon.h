/*
 * SymbolicPolygon.h
 *
 *  Created on: Aug 2, 2014
 *      Author: sam
 */

#ifndef SYMBOLIC_POLYGON
#define SYMBOLIC_POLYGON
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
#include <vector>
template<unsigned int D,class T>
class SymbolicPolygon:public GeometryDetector{
protected:
	GeometryTheorySolver<D,T> * solver;
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

	SymbolicPolygon(int detectorID,GeometryTheorySolver<D,T> * solver,  double seed=1):GeometryDetector(detectorID),solver(solver),rnd_seed(seed){

	}

	virtual ~SymbolicPolygon(){

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

	ConstantPolygon(Polygon<D,T>* polygon,int detectorID,GeometryTheorySolver<D,T> * solver,  double seed=1): SymbolicPolygon<D,T>(detectorID,solver,seed), polygon(polygon){

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

	BinaryOperationPolygon(SymbolicPolygon<D,T> & a, SymbolicPolygon<D,T> & b,PolygonOperationType operation,int detectorID,GeometryTheorySolver<D,T> * solver,  double seed=1): SymbolicPolygon<D,T>(detectorID,solver,seed),  a(a),b(b),operation(operation){

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
		Polygon<D,T> & inverse = *hp.inverse(&polygon_tmp2);
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

	ConditionalPolygon(Lit condition, SymbolicPolygon<D,T> & thn, SymbolicPolygon<D,T> & els,int detectorID,GeometryTheorySolver<D,T> * solver,  double seed=1): SymbolicPolygon<D,T>(detectorID,solver,seed), thn(thn),els(els){
		condition= mkLit( this->solver->newVar(var(condition),this->getID()),sign(condition));
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



#endif


