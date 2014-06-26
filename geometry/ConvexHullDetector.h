
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
		void buildPointContainedReason(const Point<D,T> & p,vec<Lit> & conflict);
		void buildPointNotContainedReason(const Point<D,T> & p,vec<Lit> & conflict);


		void buildReason(Lit p, vec<Lit> & reason, CRef marker);
		bool checkSatisfied();
		Lit decide();
		void addAreaDetectorLit(double areaGreaterEqThan, Var v);
		void addPointContainmentLit(Point<D,T> p,Var outerVar);

		ConvexHullDetector(int detectorID,PointSet<D,T> & under, PointSet<D,T> & over, GeometryTheorySolver<D,T> * outer,  double seed=1);
		virtual ~ConvexHullDetector(){

		}
private:
		void findFarPoints(Line<D,T> & testline,vec<Point<D,T> > & test_set,vec<Point<D,T> > & min_set,ConvexPolygon<D,T> &hull){
			test_set.clear();
			T hullside =0;
			for(int i =0;hullside==0 && i<hull.getVerticeS().size();i++){
				hullside =testline.whichSide()(hull.getVertices()[i]);
			}
			assert(hullside!=0);
			for(int i = 0;i<over.size();i++){
				if(!over.pointEnabled(i)){
					if(testline.whichSide(over[i])!=hullside){
						test_set.push(over[i]);
						if(hasSet && test_set.size()>=min_set.size()){
							return;//shortcut
						}
					}
				}
			}
			if(test_set.size()<min_set.size() || !hasSet){
				hasSet=true;
				test_set.copyTo(min_set);
			}
		}

		void findContainingTriangle2d_helper(ConvexPolygon<D,T> & polygon, int first_vertex, int last_vertex, Point<D,T> & point, ConvexPolygon<D,T> & triangle_out){
			//recurse on this segment of the polygon, finding a triangle that contains the point.
			//precondition: the point is contained in this convex segment of the polygon

			//Noah's algorithm: pick 3 vertices in the polygon. 2 of them are adjacent, and the third is arbitrary (but should probably be the index that is farthest in both directions from the adjacent vertices)
			 //Check if they contain the point; if they do, return them.
			 //Else, check which of the two sides of the triangle with the non-adjacent vertex the point is on. Recurse on that sub polygon.

			 //When recursing, 2 of the three vertices are already selected (they are the vertices from the existing triangle), so we only have to pick one more vertex.
			 //Since we already know that the point isn't on the other side of those two vertices, we only have to check two sides in the case where the point is not contained.
			assert(first_vertex!=last_vertex);
			assert(polygon.contains(point,first_vertex,last_vertex));
			triangle_out.clear();
			vec<Point<2,T>> & polygon_vertices = polygon.getVertices();
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
					 findContainingTriangle2d_helper(polygon_vertices,first_vertex,mid_point,point,triangle_out);
				}else{
#ifndef NDEBUG
					Line<2,T> dbgLine(b,c);
					assert(dbgLine.whichSide(point)!= 0);//else we would have already found the point
					assert(dbgLine.whichSide(point)!= dbgLine.whichSide(a));
#endif
					findContainingTriangle2d_helper(polygon_vertices,mid_point,last_vertex,point,triangle_out);
				}
			}

		}

		void findContainingTriangle2d( ConvexPolygon<D,T> & polygon,  Point<2,T> & point, ConvexPolygon<2,T> & triangle_out){
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

	vec<Point<D,T>> enabled_points;
/*	for (auto & p:over.getEnabledPoints(enabled_points)){
		printf("(%f, %f)",p.x,p.y);
	}
	printf("\n");*/
	MonotoneConvexHull<D,T> cv(over);

	T area = cv.getHull().getArea();
	T expectOver = over_hull->getHull().getArea();
	T expectUnder =under_hull->getHull().getArea();
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
void ConvexHullDetector<D,T>::buildPointContainedReason(const Point<D,T> & s,vec<Lit> & conflict){
	//A point is contained because there exists three enabled points that form a containing triangle
	//(might have to handle the edge case of a 2 point hull depending on how we interpret such a hull)
	//its not immediately clear how we would pick which 3 points to include, given that there may be a choice
	ConvexPolygon<D,T> &hull = under_hull->getHull();
	assert(hull.contains(s));

	if(D==2){
		ConvexPolygon<D,T> triangle;
		findContainingTriangle2d(hull,s,triangle);
		assert(triangle.contains(s));
		for(auto & p: triangle){
			int id = p.getID();
			Var v = outer->getPointVar(id);
			assert(outer->value(v)==l_True);
			conflict.push(mkLit(v,true));
		}
	}else{
		assert(false);//not yet implemented
	}

/*
	 //Set the first vertex of the containing triangle to be (arbitrarily) point 0.
	 //then walk around the hull until we find a vertex that is on the other side of the point - that will be vertex 2.
	 //then we continue on to find vertex 3

	 //note: this can also compute the area (which is the sum of p2[0]*p1[1] - p1[0]*p2[1]); could potentially combine these...
	 for(int i = 1;i<w.size();i++){
		 Point<2,double> p1 = i>0 ? (w[i-1]-s):((w.last()-s));
		 Point<2,double> p2 = w[i]-s;
		 bool contained = (p2[0]*p1[1] - p1[0]*p2[1]) >0;
		 if(!contained){
			 return false;
		 }
	 }*/
}

template<unsigned int D,class T>
void ConvexHullDetector<D,T>::buildPointNotContainedReason(const Point<D,T> & s, vec<Lit> & conflict){
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
	if(D==2){
		Line<D,T> testline;


		testline.a = s;
		ConvexPolygon<D,T> &hull = over_hull->getHull();
		Line<D,T> top_line;
		Line<D,T> bottom_line;
		top_line.a = s;
		bottom_line.a=s;
		top_line.b = hull.getVertices()[0];
		bottom_line.b = hull.getVertices().last();
		/*PointSet<D,T> pt;
		pt.addPoint(s,0);
		pt.setPointEnabled(0,true);
		if(over.pointEnabled(i)){
			int id = pt.addPoint(over[i]);
			pt.setPointEnabled(id,true)
		}
		MonotoneConvexHull<D,T> mt(pt);
		*/
		vec<Point<D,T> > min_set;
		vec<Point<D,T> > test_set;
		hasSet=false;

		for(int i = 0;i<over.size();i++){
			Point<D,T> & p = over[i];
			if(!over.pointEnabled(i)){
				if(!hull.contains(p)){
					testline.b = p;

					//now check to see if test line intersects the hull (it can still intersect on the otherside
					if(!testline.intersects(hull)){
						//then this is a point outside the hull that we can try testing.
						findFarPoints(testline,test_set,min_set,hull);
					}
				}
			}else{
				assert(hull.contains(p));
			}
		}
		//find the extreme points on the hull itself
		for(auto & p:hull.getVertices()){
			if(top_line.whichSide(p)>=0){
				top_line.b=p;
			}

			if(bottom_line.whichSide(p)<=0){
				bottom_line.b=p;
			}
		}
		findFarPoints(top_line,test_set,min_set,hull);
		findFarPoints(bottom_line,test_set,min_set,hull);

	}else{
		assert(false);//not yet implemented
	}

}

#endif
