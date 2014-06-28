
#ifndef CONVEX_DETECTOR_H_
#define CONVEX_DETECTOR_H_
#include "core/SolverTypes.h"
#include "PointSet.h"
#include "GeometryDetector.h"
#include "core/Config.h"
#include "ConvexHull.h"
#include "MonotoneConvexHull.h"
#include "QuickConvexHull.h"
#include "Polygon.h"
#include "Line.h"
using namespace Minisat;



template<unsigned int D, class T>
class ConvexHullDetector:public GeometryDetector{
		bool hasSet=false;
public:
		GeometryTheorySolver<D,T> * outer;
		//int within;
		PointSet<D,T> & under;
		PointSet<D,T> & over;
		ConvexHull<D,T>* under_hull;
		ConvexHull<D,T>* over_hull;

		double rnd_seed;
		CRef point_contained_marker;
		CRef point_not_contained_marker;
		CRef area_geq_marker;
		CRef area_not_geq_marker;
		CRef point_on_hull_marker;
		CRef point_not_on_hull_marker;

		Var lowest_point_var;
		vec<int> point_lit_map;
		vec<Lit> lit_point_map;
		vec<bool> under_hull_members;
		vec<bool> over_hull_members;
		struct PointContainedLit{
			Point<D,T> p;
			Lit l;
		};
		vec<PointContainedLit> pointContainedLits;
		struct PointOnHullLit{
			Var pointVar;
			Point<D,T> p;
			Lit l;
		};
		vec<PointOnHullLit> pointOnHullLits;

		struct AreaLit{
			T areaGreaterEqThan;
			Lit l;
		};
		vec<AreaLit> areaDetectors;

		int qhead;


		bool propagate(vec<Lit> & conflict);
		void buildAreaGEQReason(T area, vec<Lit> & conflict);
		void buildAreaLTReason(T area,vec<Lit> & conflict);
		void buildPointContainedReason(const Point<D,T> & p,vec<Lit> & conflict);
		void buildPointNotContainedReason(const Point<D,T> & p,vec<Lit> & conflict);
		void buildPointOnHullOrDisabledReason(Var pointVar,const Point<D,T> & p, vec<Lit> & conflict);
		void buildPointNotOnHullOrDisabledReason(Var pointVar,const Point<D,T> & p, vec<Lit> & conflict);

		void buildReason(Lit p, vec<Lit> & reason, CRef marker);
		bool checkSatisfied();
		Lit decide();
		void addAreaDetectorLit(T areaGreaterEqThan, Var v);
		void addPointContainmentLit(Point<D,T> p,Var outerVar);
		void addPointOnHullLit(Var pointVar,Var outerVar);
		ConvexHullDetector(int detectorID,PointSet<D,T> & under, PointSet<D,T> & over, GeometryTheorySolver<D,T> * outer,  double seed=1);
		virtual ~ConvexHullDetector(){

		}
private:
		void findFarPoints(Line<D,T> & testline, bool includeCollinearPoints,std::vector<Point<D,T> > & test_set,std::vector<Point<D,T> > & min_set,ConvexPolygon<D,T> &hull){
			test_set.clear();
			T hullside =0;
			for(int i =0;hullside==0 && i<hull.getVertices().size();i++){
				hullside =testline.whichSide(hull.getVertices()[i]);
			}
			assert(hullside!=0);
			for(int i = 0;i<over.size();i++){
				if(!over.pointEnabled(i)){
					int side = testline.whichSide(over[i]);
					if(side!=hullside && (includeCollinearPoints ||side!=0)){
						test_set.push_back(over[i]);
						if(hasSet && test_set.size()>=min_set.size()){
							return;//shortcut
						}
					}
				}
			}
			if(test_set.size()<min_set.size() || !hasSet){
				hasSet=true;
				min_set = test_set;
				//test_set.copyTo(min_set);
			}
		}

		void findContainingTriangle2d_helper(ConvexPolygon<2,T> & polygon, int first_vertex, int last_vertex,const Point<2,T> & point, ConvexPolygon<2,T> & triangle_out){
			//recurse on this segment of the polygon, finding a triangle that contains the point.
			//precondition: the point is contained in this convex segment of the polygon

			//Noah's algorithm: pick 3 vertices in the polygon. 2 of them are adjacent, and the third is arbitrary (but should probably be the index that is farthest in both directions from the adjacent vertices)
			 //Check if they contain the point; if they do, return them.
			 //Else, check which of the two sides of the triangle with the non-adjacent vertex the point is on. Recurse on that sub polygon.

			 //When recursing, 2 of the three vertices are already selected (they are the vertices from the existing triangle), so we only have to pick one more vertex.
			 //Since we already know that the point isn't on the other side of those two vertices, we only have to check two sides in the case where the point is not contained.
			assert(first_vertex!=last_vertex);
			assert(polygon.containsInRange(point,first_vertex,last_vertex));
			triangle_out.clear();
			std::vector<Point<2,T>> & polygon_vertices = polygon.getVertices();
			Point<2,T> & a = polygon_vertices[first_vertex];
			Point<2,T> & b = polygon_vertices[last_vertex];
			int mid_point = 0;
			if(first_vertex<last_vertex){
				mid_point = (last_vertex-first_vertex)/2 + first_vertex;
			}else{
				mid_point = (first_vertex-last_vertex)/2 + last_vertex;
			}
			Point<2,T> & c = polygon_vertices[mid_point];
			triangle_out.addVertexUnchecked(a);
			triangle_out.addVertexUnchecked(b);
			triangle_out.addVertexUnchecked(c);

			if(triangle_out.contains(point))
				return;//we are done
			else{
				Line<2,T> testLine(a,c);
				assert(testLine.whichSide(point)!= 0);//else we would have already found the point

				if(testLine.whichSide(point)!= testLine.whichSide(b)){
					 findContainingTriangle2d_helper(polygon,first_vertex,mid_point,point,triangle_out);
				}else{
#ifndef NDEBUG
					Line<2,T> dbgLine(b,c);
					assert(dbgLine.whichSide(point)!= 0);//else we would have already found the point
					assert(dbgLine.whichSide(point)!= dbgLine.whichSide(a));
#endif
					findContainingTriangle2d_helper(polygon,mid_point,last_vertex,point,triangle_out);
				}
			}

		}

		void findContainingTriangle2d( ConvexPolygon<D,T> & polygon,const Point<D,T> & point, ConvexPolygon<D,T> & triangle_out){
			assert(polygon.contains(point));

			findContainingTriangle2d_helper(polygon,0,polygon.getVertices().size()-1,point,triangle_out);
		}

		/*

		void findContainingTriangle2d( ConvexPolygon<D,T> & polygon,  Point<2,T> & point, ConvexPolygon<2,T> & triangle_out){
			 //Noah's algorithm: pick 3 equidistant (in terms of index in the list of vertices) vertices.
			 //Check if they contain the point; if they do, return them.
			 //Else, check which of the three sides of the triangle the point is on. Recurse on that sub polygon.

			 //In this implementation, when recursing, 2 of the three vertices are already selected (they are the vertices from the existing triangle), so we only have to pick one more vertex.
			 //Since we already know that the point isn't on the other side of those two vertices, we only have to check two sides in the case where the point is not contained.
			 vec<Point<2,T>> & polygon_vertices = polygon.getVertices();
			 int skip = polygon_vertices.size()/3;
			 Point<2,T> & a = w[0];
			 Point<2,T> & b = w[skip];
			 Point<2,T> & c = w[2*skip];

			 triangle_out.clear();
			 triangle_out.addVertex(a);triangle_out.addVertex(b); triangle_out.addVertex(c);
			 if(triangle_out.contains(point)){
				 return;
			 }
			 //Check which of the three sides of the triangle contains the point
			 Line<2,T> testLine;
			 testLine.a = a;
			 testLine.b = b;
			 if(testLine.whichSide(point)!= testLine.whichSide(c)){
				 findContainingTriangle2d_helper(polygon_vertices,0,skip,point,triangle_out);
				 return;
			 }
			 testLine.a = b;
			 testLine.b = c;
			 if(testLine.whichSide(point)!= testLine.whichSide(c)){
				 findContainingTriangle2d_helper(polygon_vertices,skip,2*skip,point,triangle_out);
				 return;
			 }
			 findContainingTriangle2d_helper(polygon_vertices,2*skip,0,point,triangle_out);
		}*/
};

template<unsigned int D, class T>
ConvexHullDetector<D,T>::ConvexHullDetector(int detectorID,PointSet<D,T> & under, PointSet<D,T> & over, GeometryTheorySolver<D,T> * outer,double seed):
GeometryDetector(detectorID),outer(outer),under(under),over(over),rnd_seed(seed){

	point_contained_marker=outer->newReasonMarker(getID());
	point_not_contained_marker=outer->newReasonMarker(getID());

	point_on_hull_marker=outer->newReasonMarker(getID());
	point_not_on_hull_marker=outer->newReasonMarker(getID());
	area_geq_marker=outer->newReasonMarker(getID());
	area_not_geq_marker=outer->newReasonMarker(getID());

	if(hullAlg== ConvexHullAlg::ALG_QUICKHULL){
		over_hull = new QuickConvexHull<D,T>(over);
		under_hull = new QuickConvexHull<D,T>(under);
	}else if (hullAlg== ConvexHullAlg::ALG_MONOTONE_HULL){

		over_hull = new MonotoneConvexHull<D,T>(over);
		under_hull = new MonotoneConvexHull<D,T>(under);

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
void ConvexHullDetector<D,T>::addAreaDetectorLit(T areaGreaterEqThan, Var outerVar){
	Var v = outer->newVar(outerVar,getID());
	Lit l = mkLit(v,false);

	areaDetectors.push({areaGreaterEqThan,l});

}


template<unsigned int D, class T>
void ConvexHullDetector<D,T>::buildReason(Lit p, vec<Lit> & reason, CRef marker){
		reason.push(p);
		if(marker==area_geq_marker){

			for(auto & a:areaDetectors){
				if(var(a.l)==var(p)){
					buildAreaGEQReason(a.areaGreaterEqThan, reason);
					return;
				}
			}


		}else if(marker==area_not_geq_marker){
			for(auto & a:areaDetectors){
				if(var(a.l)==var(p)){
					buildAreaLTReason(a.areaGreaterEqThan, reason);
					return;
				}
			}

		}else if(marker==point_contained_marker){
			for(auto & a:pointContainedLits){
				if(var(a.l)==var(p)){
					buildPointContainedReason(a.p, reason);
					return;
				}
			}


		}else if(marker==point_not_contained_marker){
			for(auto & a:pointContainedLits){
				if(var(a.l)==var(p)){
					buildPointNotContainedReason(a.p, reason);
					return;
				}
			}



		}else if(marker==point_on_hull_marker){
			for(auto & a:pointOnHullLits){
				if(var(a.l)==var(p)){
					buildPointOnHullOrDisabledReason(a.pointVar,a.p, reason);
					return;
				}
			}


		}else if(marker==point_not_on_hull_marker){
			for(auto & a:pointOnHullLits){
				if(var(a.l)==var(p)){
					buildPointNotOnHullOrDisabledReason(a.pointVar,a.p, reason);
					return;
				}
			}


		}else{
			assert(false);
		}
		assert(false);
}

template<unsigned int D, class T>
bool ConvexHullDetector<D,T>::propagate(vec<Lit> & conflict){

	static int iter = 0;
	if(++iter==9){
		int a=1;
	}

		over_hull->update();

		under_hull->update();

		if(areaDetectors.size()){

			T over_area = over_hull->getHull().getArea();
			T under_area = under_hull->getHull().getArea();
			double tover =  over_area.get_d();
			double tunder = under_area.get_d();
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
						buildAreaGEQReason(under_area,conflict);
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
						buildAreaLTReason(over_area,conflict);
						return false;
					}
				}


			}
		}
		ConvexPolygon<D,T> & p_over = over_hull->getHull();
		ConvexPolygon<D,T> & p_under = under_hull->getHull();
		//If we are making many queries, it is probably worth it to pre-process the polygon and then make the queries.
		for(int i =0;i<pointContainedLits.size();i++){
			Point<D,T> & point = pointContainedLits[i].p;
			Lit l = pointContainedLits[i].l;


			if(p_under.contains(point)){
				//l is true
				if(outer->value(l)==l_True){
					//do nothing
				}else if(outer->value(l)==l_Undef){

					outer->enqueue(l,point_contained_marker) ;
				}else if (outer->value(l)==l_False){
					conflict.push(l);
					buildPointContainedReason(point,conflict);
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
					buildPointNotContainedReason(point,conflict);
					return false;
				}
			}
		}
		if(pointOnHullLits.size()){
#ifndef NDEBUG
			for(bool b:under_hull_members)
				assert(!b);
			for(bool b:over_hull_members)
				assert(!b);
#endif
			under_hull_members.growTo(over.size());
			over_hull_members.growTo(under.size());
			for(auto & p:p_under){
				under_hull_members[p.getID()]=true;
			}
			for(auto & p:p_over){
				over_hull_members[p.getID()]=true;
			}
			for(int i =0;i<pointOnHullLits.size();i++){
				Point<D,T> & point = pointOnHullLits[i].p;
				Lit l = pointOnHullLits[i].l;

				//check if a point is a member of the hull

				bool under_member=under_hull_members[point.getID()];
				bool over_member=over_hull_members[point.getID()];

				if(!under.pointEnabled(point.getID()) ||under_member){
					//l is true
					if(outer->value(l)==l_True){
						//do nothing
					}else if(outer->value(l)==l_Undef){

						outer->enqueue(l,point_contained_marker) ;
					}else if (outer->value(l)==l_False){
						conflict.push(l);
						buildPointOnHullOrDisabledReason( pointOnHullLits[i].pointVar, point,conflict);
						return false;
					}
				}else if (!over_member){
					l=~l;
					//l is true
					if(outer->value(l)==l_True){
						//do nothing
					}else if(outer->value(l)==l_Undef){
						outer->enqueue(l,point_not_contained_marker) ;
					}else if (outer->value(l)==l_False){
						conflict.push(l);
						buildPointNotOnHullOrDisabledReason(pointOnHullLits[i].pointVar,point,conflict);
						return false;
					}
				}
			}
			for(auto & p:p_under){
				under_hull_members[p.getID()]=false;
			}
			for(auto & p:p_over){
				over_hull_members[p.getID()]=false;
			}
		}
		return true;
	}
template<unsigned int D, class T>
bool ConvexHullDetector<D,T>::checkSatisfied(){

	std::vector<Point<D,T>> enabled_points;
/*	for (auto & p:over.getEnabledPoints(enabled_points)){
		printf("(%f, %f)",p.x,p.y);
	}
	printf("\n");*/
	MonotoneConvexHull<D,T> cv(over);
	T expectOver = over_hull->getHull().getArea();
	T expectUnder =under_hull->getHull().getArea();
	T area = cv.getHull().getArea();

	assert(equal_epsilon(area, expectOver));
	assert(equal_epsilon(area, expectUnder));
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
	//the reason that the area is greater or equal to the current value is the set of points in the convex hull (all of which are enabled).
	for(auto & p:under_hull->getHull()){
		int pID = p.getID();
		assert(under.pointEnabled(pID));
		conflict.push(mkLit(outer->getPointVar(pID),false));

	}
}
template<unsigned int D, class T>
void ConvexHullDetector<D,T>::buildAreaLTReason(T area,vec<Lit> & conflict){
	//the reason that the area is less than some value is that some point that is OUTSIDE the convex hull is not enabled.
	for(int i = 0;i<over.size();i++){
		if(!over.pointEnabled(i) && ! over_hull->getHull().contains(over[i])){
			conflict.push(mkLit(outer->getPointVar(i),false));
		}
	}
}
template<unsigned int D, class T>
void ConvexHullDetector<D,T>::buildPointContainedReason(const Point<D,T> & s, vec<Lit> & conflict){

}

template<unsigned int D, class T>
void ConvexHullDetector<D,T>::buildPointNotContainedReason(const Point<D,T> & s, vec<Lit> & conflict){

}
template<unsigned int D, class T>
void ConvexHullDetector<D,T>::buildPointOnHullOrDisabledReason(Var pointVar,const Point<D,T> & p, vec<Lit> & conflict){

}

template<unsigned int D, class T>
void ConvexHullDetector<D,T>::buildPointNotOnHullOrDisabledReason(Var pointVar,const Point<D,T> & p, vec<Lit> & conflict){

}

template<>
void ConvexHullDetector<2,double>::buildPointOnHullOrDisabledReason(Var pointVar,const Point<2,double> & p, vec<Lit> & conflict);
template< >
void ConvexHullDetector<2,double>::buildPointNotOnHullOrDisabledReason(Var pointVar,const Point<2,double> & p, vec<Lit> & conflict);
template<>
void ConvexHullDetector<2,double>::buildPointContainedReason(const Point<2,double> & s,vec<Lit> & conflict);
template<>
void ConvexHullDetector<2,double>::buildPointNotContainedReason(const Point<2,double> & s, vec<Lit> & conflict);


template<>
void ConvexHullDetector<2,mpq_class>::buildPointOnHullOrDisabledReason(Var pointVar,const Point<2,mpq_class> & p, vec<Lit> & conflict);
template< >
void ConvexHullDetector<2,mpq_class>::buildPointNotOnHullOrDisabledReason(Var pointVar,const Point<2,mpq_class> & p, vec<Lit> & conflict);
template<>
void ConvexHullDetector<2,mpq_class>::buildPointContainedReason(const Point<2,mpq_class> & s,vec<Lit> & conflict);
template<>
void ConvexHullDetector<2,mpq_class>::buildPointNotContainedReason(const Point<2,mpq_class> & s, vec<Lit> & conflict);

#endif
