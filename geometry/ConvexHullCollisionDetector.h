#ifndef CONVEX_HULL_COLLISION_DETECTOR_H_
#define CONVEX_HULL_COLLISION_DETECTOR_H_

#include "core/SolverTypes.h"
#include "PointSet.h"
#include "GeometryDetector.h"
#include "core/Config.h"
#include "ConvexHull.h"
#include "MonotoneConvexHull.h"
#include "QuickConvexHull.h"
#include "Polygon.h"
#include "Line.h"
#include "ConvexHullDetector.h"

using namespace Minisat;



template<unsigned int D, class T>
class ConvexHullCollisionDetector:public GeometryDetector{
		bool hasSet=false;
public:
		GeometryTheorySolver<D,T> * outer;
		//int within;
		PointSet<D,T> & under1;
		PointSet<D,T> & over1;
		PointSet<D,T> & under2;
		PointSet<D,T> & over2;
		//The idea here is to avoid recomputing convex hulls if they are already being computed elsewhere.
		//It would be cleaner to just recompute the hulls here.

		//Down the line, this should really expand to handle more than just pairwise comparisons, and apply bounding box detection as well...
		ConvexHull<D,T> & hull1_under;
		ConvexHull<D,T> & hull1_over;

		ConvexHull<D,T> & hull2_under;
		ConvexHull<D,T> & hull2_over;

		double rnd_seed;
		//Intersection means they don't just collide, but also pass through each other
		CRef hulls_intersect_marker;
		CRef hulls_not_intersect_marker;

		//Collision includes the case where one or more points are exactly shared between the two hulls
		CRef hulls_collision_marker;
		CRef hulls_not_collision_marker;

		Lit collisionLit=lit_Undef;
		Lit intersectionLit = lit_Undef;

		bool propagate(vec<Lit> & conflict);
		void buildCollisionReason(vec<Lit> & conflict);
		void buildNotCollisionReason(vec<Lit> & conflict);

		void buildReason(Lit p, vec<Lit> & reason, CRef marker);
		bool checkSatisfied();

		void addCollisionDetectorLit(Var v, bool require_intersection);

		ConvexHullCollisionDetector(int detectorID,GeometryTheorySolver<D,T> * outer,PointSet<D,T> & under1, PointSet<D,T> & over1,PointSet<D,T> & under2, PointSet<D,T> & over2,ConvexHullDetector<D,T> & hull1_under,ConvexHullDetector<D,T> & hull1_over,ConvexHullDetector<D,T> & hull2_under,ConvexHullDetector<D,T> & hull2_over, double seed=1);
private:
		bool findSeparatingAxis(ConvexPolygon<D,T> & h1, ConvexPolygon<D,T> & h2, PointSet<2,mpq_class> & pointset1, PointSet<2,mpq_class> & pointset2,  std::vector<std::pair<Point<2,mpq_class> ,mpq_class>>  &projection_out1,std::vector<std::pair<Point<2,mpq_class> ,mpq_class>>  &projection_out2);


		void findContainingTriangle2d_helper(ConvexPolygon<2,T> & polygon, int first_vertex, int last_vertex,const Point<2,T> & point, NConvexPolygon<2,T> & triangle_out){
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
			ConvexPolygon<2,T> & polygon_vertices =(ConvexPolygon<2,T> &) *this;
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

		void findContainingTriangle2d( ConvexPolygon<D,T> & polygon,const Point<D,T> & point, NConvexPolygon<D,T> & triangle_out){
			assert(polygon.contains(point));

			findContainingTriangle2d_helper(polygon,0,polygon.size()-1,point,triangle_out);
		}

};

template<unsigned int D, class T>
ConvexHullCollisionDetector<D,T>::ConvexHullCollisionDetector(int detectorID,GeometryTheorySolver<D,T> * outer,PointSet<D,T> & under1, PointSet<D,T> & over1,PointSet<D,T> & under2, PointSet<D,T> & over2,ConvexHullDetector<D,T> & hull1_under,ConvexHullDetector<D,T> & hull1_over,ConvexHullDetector<D,T> & hull2_under,ConvexHullDetector<D,T> & hull2_over, double seed):
GeometryDetector(detectorID),outer(outer),under1(under1),over1(over1),under2(under2),over2(over2),hull1_under(hull1_under),hull1_over(hull1_over),hull2_under(hull2_under),hull2_over(hull2_over),rnd_seed(seed)
{
	hulls_intersect_marker = outer->newReasonMarker(detectorID);
	hulls_not_intersect_marker = outer->newReasonMarker(detectorID);

	hulls_collision_marker = outer->newReasonMarker(detectorID);
	hulls_not_collision_marker = outer->newReasonMarker(detectorID);


}
template<unsigned int D, class T>
void ConvexHullCollisionDetector<D,T>::addCollisionDetectorLit(Var outerVar, bool require_intersection){
	if(require_intersection){
		if(this->intersectionLit ==lit_Undef){
			Var v = outer->newVar(outerVar,getID());
			Lit l = mkLit(v,false);
			intersectionLit=l;
		}else{
			outer->makeEqualInSolver(mkLit(outerVar),outer->toSolver(intersectionLit));
		}
	}else{
		if(collisionLit==lit_Undef){
			Var v = outer->newVar(outerVar,getID());
			Lit l = mkLit(v,false);
			collisionLit=l;
		}else{
			outer->makeEqualInSolver(mkLit(outerVar),outer->toSolver(collisionLit));
		}
	}
}
template<unsigned int D, class T>
bool ConvexHullCollisionDetector<D,T>::propagate(vec<Lit> & conflict){

	//check whether the two convex hulls intersect by checking their
	auto & h1_under = hull1_under.getHull();
	auto & h1_over = hull1_over.getHull();
	auto & h2_under = hull2_under.getHull();
	auto & h2_over = hull2_over.getHull();

	if(h1_under.intersects(h2_under)){
		Lit l = collisionLit;
		if(outer->value(l)==l_True){
			//do nothing
		}else if(outer->value(l)==l_Undef){
			outer->enqueue(l,hulls_collision_marker) ;
		}else if (outer->value(l)==l_False){
			conflict.push(l);
			buildCollisionReason(conflict);
			return false;
		}
	}else if (not h1_under.intersects(h2_under)){
		Lit l = ~collisionLit;
		if(outer->value(l)==l_True){
			//do nothing
		}else if(outer->value(l)==l_Undef){
			outer->enqueue(l,hulls_not_collision_marker) ;
		}else if (outer->value(l)==l_False){
			conflict.push(l);
			buildNotCollisionReason(conflict);
			return false;
		}
	}

	return false;
}
template<unsigned int D, class T>
void ConvexHullCollisionDetector<D,T>::buildReason(Lit p, vec<Lit> & reason, CRef marker){
	reason.push(p);
	if(marker==hulls_collision_marker){
		buildCollisionReason(reason);
	}else if(marker==hulls_not_collision_marker){
		buildNotCollisionReason(reason);
	}else if(marker==hulls_intersect_marker){

	}else if(marker==hulls_not_intersect_marker){


	}else{
		assert(false);
	}
}

template<unsigned int D, class T>
bool ConvexHullCollisionDetector<D,T>::checkSatisfied(){
	auto & h1_under = hull1_under.getHull();
	auto & h1_over = hull1_over.getHull();
	auto & h2_under = hull2_under.getHull();
	auto & h2_over = hull2_over.getHull();
	if (collisionLit!=lit_Undef){
		if (outer->value(collisionLit)==l_True){
			if(! h1_under.intersects(h2_under)){
				return false;
			}
		}else if (outer->value(collisionLit)==l_False){
			if(h1_over.intersects(h2_over)){
				return false;
			}
		}
	}
	return true;
}

template<>
void ConvexHullCollisionDetector<2, mpq_class>::buildCollisionReason(vec<Lit> & conflict);

template<>
void ConvexHullCollisionDetector<2, mpq_class>::buildNotCollisionReason(vec<Lit> & conflict);

template<>
void ConvexHullCollisionDetector<2, double>::buildCollisionReason(vec<Lit> & conflict);

template<>
void ConvexHullCollisionDetector<2, double>::buildNotCollisionReason(vec<Lit> & conflict);

template<>
bool ConvexHullCollisionDetector<2, mpq_class>::findSeparatingAxis(ConvexPolygon<2, mpq_class> & h1, ConvexPolygon<2, mpq_class> & h2, PointSet<2,mpq_class> & pointset1, PointSet<2,mpq_class> & pointset2,  std::vector<std::pair<Point<2,mpq_class> ,mpq_class>>  &projection_out1,std::vector<std::pair<Point<2,mpq_class> ,mpq_class>>  &projection_out2);


#endif
