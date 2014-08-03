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

	Var lowest_point_var;
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

	virtual Polygon<D,T>* getOverApprox()=0;
	virtual Polygon<D,T>* getUnderApprox()=0;
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

		Polygon<D,T> & p_over = *getOverApprox();
		Polygon<D,T> & p_under = *getUnderApprox();


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
						buildAreaGEQReason(under_area,conflict);
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
						buildAreaLTReason(over_area,conflict);
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
		for(int i =0;i<lineIntersectionLits.size();i++){

			LineSegment<D,T> & line = lineIntersectionLits[i].line;
			Lit l = lineIntersectionLits[i].l;

			if(checkLineIntersection(lineIntersectionLits[i].under_intersecting_polygon,line,p_under,lineIntersectionLits[i].inclusive) || p_under.intersects(line,&lineIntersectionLits[i].under_intersecting_polygon,nullptr,lineIntersectionLits[i].inclusive)){
				//l is true
				if(solver->value(l)==l_True){
					//do nothing
				}else if(solver->value(l)==l_Undef){

					solver->enqueue(l,line_intersection_marker) ;
				}else if (solver->value(l)==l_False){
					conflict.push(l);
					buildLineIntersectsReason(line,conflict,lineIntersectionLits[i].inclusive);
					return false;
				}
			}else if (!checkLineIntersection(lineIntersectionLits[i].over_intersecting_polygon,line,p_over, lineIntersectionLits[i].inclusive) && !p_over.intersects(line,&lineIntersectionLits[i].over_intersecting_polygon,nullptr,lineIntersectionLits[i].inclusive)){
			//}else if (!p_over.intersects(line)){
				l=~l;
				//l is true
				if(solver->value(l)==l_True){
					//do nothing
				}else if(solver->value(l)==l_Undef){
					solver->enqueue(l,line_not_intersection_marker) ;
				}else if (solver->value(l)==l_False){
					conflict.push(l);
					buildLineNotIntersectsReason(line,conflict,lineIntersectionLits[i].inclusive);
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
					buildConvexIntersectsReason(polygon,conflict,convexIntersectionLits[i].inclusive);
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
					buildConvexNotIntersectsReason(polygon,conflict,convexIntersectionLits[i].inclusive);
					return false;
				}
			}
		}

		return true;
	}

	void buildReason(Lit p, vec<Lit> & reason, CRef marker){

	}


	void buildAreaGEQReason(T area, vec<Lit> & conflict){

	}

	void buildAreaLTReason(T area,vec<Lit> & conflict){

	}

	void buildPointContainedReason(const Point<D,T> & s, vec<Lit> & conflict, bool inclusive){


	}


	void buildPointNotContainedReason(const Point<D,T> & s, vec<Lit> & conflict, bool inclusive){

	}

	void buildLineIntersectsReason(LineSegment<D,T> & line,vec<Lit> & conflict, bool inclusive){

	}

	void buildLineNotIntersectsReason(LineSegment<D,T> & line,vec<Lit> & conflict, bool inclusive){

	}


	void buildConvexIntersectsReason(ConvexPolygon<D,T> & polygon,vec<Lit> & conflict, bool inclusive){

	}

	void buildConvexNotIntersectsReason(ConvexPolygon<D,T> & polygon,vec<Lit> & conflict, bool inclusive){

	}

	void buildPointOnHullOrDisabledReason(Var pointVar,const Point<D,T> & p, vec<Lit> & conflict){

	}


	void buildPointNotOnHullOrDisabledReason(Var pointVar,const Point<D,T> & p, vec<Lit> & conflict){

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

	ConstantPolygon(Polygon<D,T>* polygon,int detectorID,GeometryTheorySolver<D,T> * solver,  double seed=1): SymbolicPolygon<D,T>(detectorID,solver,seed), polygon(polygon){

	}

	Polygon<D,T>* getOverApprox(){
		return polygon;
	}
	Polygon<D,T>* getUnderApprox(){
		return polygon;
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

	Polygon<D,T> * applyOperation(Polygon<D,T> * p_a, Polygon<D,T> * p_b, NPolygon<D,T> * store){
		switch (operation){
			case PolygonOperationType::op_union:{
				return p_a->binary_union(p_b,store);
			}break;
			case PolygonOperationType::op_intersect:{
				return p_a->binary_intersect(p_b,store);
			}break;
			case PolygonOperationType::op_difference:{
				return p_a->binary_difference(p_b,store);
			}break;
			case PolygonOperationType::op_minkowski_sum:{
				return p_a->binary_minkowski_sum(p_b,store);
			}break;
			default:
				assert(false);
				break;
		}
		return nullptr;
	}

public:

	BinaryOperationPolygon(SymbolicPolygon<D,T> & a, SymbolicPolygon<D,T> & b,PolygonOperationType operation,int detectorID,GeometryTheorySolver<D,T> * solver,  double seed=1): SymbolicPolygon<D,T>(detectorID,solver,seed),  a(a),b(b),operation(operation){

	}

	Polygon<D,T>* getOverApprox(){
		return applyOperation(a.getOverApprox(),b.getOverApprox(),&over_approx);
	}
	Polygon<D,T>* getUnderApprox(){
		return applyOperation(a.getUnderApprox(),b.getUnderApprox(),&under_approx);
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

	Polygon<D,T>* getOverApprox(){
		lbool val = this->solver->value(condition);
		if(val==l_True){
			return thn.getOverApprox();
		}else if (val==l_False){
			return els.getOverApprox();
		}
		//return the union of the two polygons
		return thn.getOverApprox()->binary_union(els.getOverApprox(),&over_approx);
	}
	Polygon<D,T>* getUnderApprox(){
		lbool val = this->solver->value(condition);
		if(val==l_True){
			return thn.getUnderApprox();
		}else if (val==l_False){
			return els.getUnderApprox();
		}
		//return the intersection of the two polygons
		return thn.getUnderApprox()->binary_intersect(els.getUnderApprox(),&under_approx);
	}
};



/*
template<unsigned int D,class T>
class SymbolicPolygon:public Shape<D,T>{
	Polygon<D,T>* over_approx=nullptr;
	Polygon<D,T>* under_approx=nullptr;

	PolygonDefinition<D,T> def;


};
*/

#endif


