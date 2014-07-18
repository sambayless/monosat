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
#include <vector>
using namespace Minisat;



template<unsigned int D, class T>
class ConvexHullCollisionDetector:public GeometryDetector{
		bool hasSet=false;
public:
		GeometryTheorySolver<D,T> * outer;
		//int within;

		std::vector<PointSet<D,T>> & under_sets;
		std::vector<PointSet<D,T>> & over_sets;
		std::vector<ConvexHullDetector<D,T>*> & convexHullDetectors;


		double rnd_seed;
		//Intersection means they don't just collide, but also pass through each other


		//Collision includes the case where one or more points are exactly shared between the two hulls
		CRef hulls_collision_marker;
		CRef hulls_not_collision_marker;


		struct CollisionLit{
			int pointSet1;
			int pointSet2;
			Lit inclusiveLit;
			Lit exclusiveLit;
		};
		std::vector<CollisionLit> collisionLits;

		std::vector<bool> has_collision_detector;
		std::vector<int> collision_pointsets;

		bool propagate(vec<Lit> & conflict);
		void buildCollisionReason(vec<Lit> & conflict, int pointSet1, int pointSet2, bool inclusive);
		void buildNotCollisionReason(vec<Lit> & conflict, int pointSet1, int pointSet2, bool inclusive);

		void buildReason(Lit p, vec<Lit> & reason, CRef marker);
		bool checkSatisfied();

		void addCollisionDetectorLit(int pointSet1, int pointSet2,Var outerVar, bool inclusive);
		ConvexHullCollisionDetector(int detectorID,GeometryTheorySolver<D,T> * outer,std::vector<PointSet<D,T>> & under_sets, std::vector<PointSet<D,T>> & over_sets,std::vector<ConvexHullDetector<D,T>*> & convexHullDetectors, double seed);

private:

		void buildCollisionReason2d(vec<Lit> & conflict, int pointSet1, int pointSet2, bool inclusive);
		void buildNotCollisionReason2d(vec<Lit> & conflict, int pointSet1, int pointSet2, bool inclusive);
		bool findSeparatingAxis(ConvexPolygon<D,T> & h1, ConvexPolygon<D,T> & h2, PointSet<D,T> & pointset1, PointSet<D,T> & pointset2,  std::vector<std::pair<Point<D,T> ,T>>  &projection_out1,std::vector<std::pair<Point<D,T> ,T>>  &projection_out2, bool inclusive){
			if(D==2){
				return findSeparatingAxis2d((ConvexPolygon<2,T> &)h1,(ConvexPolygon<2,T> &)h2,(PointSet<2,T> &)pointset1,(PointSet<2,T> &)pointset2,(std::vector<std::pair<Point<2,T> ,T>>  &)projection_out1,(std::vector<std::pair<Point<2,T> ,T>>  &)projection_out2,inclusive);
			}else{
				assert(false);
			}
			return false;
		}

		bool findSeparatingAxis2d(ConvexPolygon<2,T> & h1, ConvexPolygon<2,T> & h2, PointSet<2,T> & pointset1, PointSet<2,T> & pointset2,  std::vector<std::pair<Point<2,T> ,T>>  &projection_out1,std::vector<std::pair<Point<2,T> ,T>>  &projection_out2, bool inclusive);


		void findContainingTriangle2d_helper(ConvexPolygon<2,T> & polygon, int first_vertex, int last_vertex,const Point<2,T> & point, NConvexPolygon<2,T> & triangle_out, bool inclusive){
			//recurse on this segment of the polygon, finding a triangle that contains the point.
			//precondition: the point is contained in this convex segment of the polygon

			//Noah's algorithm: pick 3 vertices in the polygon. 2 of them are adjacent, and the third is arbitrary (but should probably be the index that is farthest in both directions from the adjacent vertices)
			 //Check if they contain the point; if they do, return them.
			 //Else, check which of the two sides of the triangle with the non-adjacent vertex the point is on. Recurse on that sub polygon.

			 //When recursing, 2 of the three vertices are already selected (they are the vertices from the existing triangle), so we only have to pick one more vertex.
			 //Since we already know that the point isn't on the other side of those two vertices, we only have to check two sides in the case where the point is not contained.
			assert(first_vertex!=last_vertex);
			assert(polygon.containsInRange(point,first_vertex,last_vertex, inclusive));
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

			if(triangle_out.contains(point,inclusive))
				return;//we are done
			else{
				Line<2,T> testLine(a,c);
				assert(testLine.whichSide(point)!= 0);//else we would have already found the point

				if(testLine.whichSide(point)!= testLine.whichSide(b)){
					 findContainingTriangle2d_helper(polygon,first_vertex,mid_point,point,triangle_out,inclusive);
				}else{
#ifndef NDEBUG
					Line<2,T> dbgLine(b,c);
					assert(dbgLine.whichSide(point)!= 0);//else we would have already found the point
					assert(dbgLine.whichSide(point)!= dbgLine.whichSide(a));
#endif
					findContainingTriangle2d_helper(polygon,mid_point,last_vertex,point,triangle_out,inclusive);
				}
			}

		}

		void findContainingTriangle2d( ConvexPolygon<D,T> & polygon,const Point<D,T> & point, NConvexPolygon<D,T> & triangle_out, bool inclusive){
			assert(polygon.contains(point, inclusive));

			findContainingTriangle2d_helper(polygon,0,polygon.size()-1,point,triangle_out, inclusive);
		}

};

template<unsigned int D, class T>
ConvexHullCollisionDetector<D,T>::ConvexHullCollisionDetector(int detectorID,GeometryTheorySolver<D,T> * outer,std::vector<PointSet<D,T>> & under_sets, std::vector<PointSet<D,T>> & over_sets,std::vector<ConvexHullDetector<D,T>*> & convexHullDetectors, double seed):
GeometryDetector(detectorID),outer(outer),under_sets(under_sets),over_sets(over_sets),convexHullDetectors(convexHullDetectors),rnd_seed(seed)
{
	hulls_collision_marker = outer->newReasonMarker(detectorID);
	hulls_not_collision_marker = outer->newReasonMarker(detectorID);


}
template<unsigned int D, class T>
void ConvexHullCollisionDetector<D,T>::addCollisionDetectorLit(int pointSet1, int pointSet2, Var outerVar,bool inclusive){
	if(pointSet1==pointSet2){
		outer->makeTrueInSolver(mkLit(outerVar));
		return;
	}else if (pointSet1>pointSet2){
		std::swap(pointSet1,pointSet2);
	}

	for(auto & p:collisionLits){
		if (p.pointSet1 == pointSet1 && p.pointSet2 ==pointSet2){
			if(inclusive && p.inclusiveLit!=lit_Undef){
				outer->makeEqualInSolver(mkLit(outerVar),outer->toSolver(p.inclusiveLit));
			}else if (!inclusive && p.exclusiveLit!=lit_Undef){
				outer->makeEqualInSolver(mkLit(outerVar),outer->toSolver(p.exclusiveLit));
			}else if (inclusive && p.inclusiveLit==lit_Undef){
				Var v = outer->newVar(outerVar,getID());
				Lit l = mkLit(v,false);
				p.inclusiveLit=l;
			}else if (!inclusive && p.exclusiveLit==lit_Undef){
				Var v = outer->newVar(outerVar,getID());
				Lit l = mkLit(v,false);
				p.exclusiveLit=l;
			}else{
				assert(false);
			}
			return;
		}
	}
	Var v = outer->newVar(outerVar,getID());
	Lit l = mkLit(v,false);
	if(inclusive){
		collisionLits.push_back({pointSet1,pointSet2,l,lit_Undef});
	}else{
		collisionLits.push_back({pointSet1,pointSet2,lit_Undef,l});
	}
	assert(pointSet2>pointSet1);

	while(has_collision_detector.size()<=pointSet2)
		has_collision_detector.push_back(false);
	if(!has_collision_detector[pointSet1]){
		collision_pointsets.push_back(pointSet1);
	}
	if(!has_collision_detector[pointSet2]){
		collision_pointsets.push_back(pointSet2);
	}
}
template<unsigned int D, class T>
bool ConvexHullCollisionDetector<D,T>::propagate(vec<Lit> & conflict){

	//for now, I'm just iterating through each possible collision and doing pairwise checks. In the future, it would be nice to make this more efficient.
	for(auto & c:collisionLits){
		//this is really ugly.
		auto & h1_under = convexHullDetectors[c.pointSet1]->getConvexHull(false)->getHull();
		auto & h1_over =  convexHullDetectors[c.pointSet1]->getConvexHull(true)->getHull();
		auto & h2_under = convexHullDetectors[c.pointSet2]->getConvexHull(false)->getHull();
		auto & h2_over = convexHullDetectors[c.pointSet2]->getConvexHull(true)->getHull();
/*		cout<<h1_under << "\n";
		cout << h2_under<<"\n";
		cout<<h1_over << "\n";
		cout << h2_over<<"\n";*/
		if(c.inclusiveLit!=lit_Undef){
			if(h1_under.intersects(h2_under,true)){
				Lit l = c.inclusiveLit;
				if(outer->value(l)==l_True){
					//do nothing
				}else if(outer->value(l)==l_Undef){
					outer->enqueue(l,hulls_collision_marker) ;
				}else if (outer->value(l)==l_False){
					conflict.push(l);
					buildCollisionReason(conflict,c.pointSet1,c.pointSet2,true);
					return false;
				}
				//If the hulls intersect inclusively, they MAY also intersect exclusively, so we need to check that too (could probably combine much of the calculations in these checks in the future)
				if (c.exclusiveLit !=lit_Undef){
					if(h1_under.intersects(h2_under,false)){
						Lit l = c.exclusiveLit;
						if(outer->value(l)==l_True){
							//do nothing
						}else if(outer->value(l)==l_Undef){
							outer->enqueue(l,hulls_collision_marker) ;
						}else if (outer->value(l)==l_False){
							conflict.push(l);
							buildCollisionReason(conflict,c.pointSet1,c.pointSet2,false);
							return false;
						}
					}else if (! h1_over.intersects(h2_over,false)){
						//If the hulls DID intersect inclusively, the we need to check if they also intersect exclusively.
						Lit l = ~c.exclusiveLit;
						if(outer->value(l)==l_True){
							//do nothing
						}else if(outer->value(l)==l_Undef){
							outer->enqueue(l,hulls_not_collision_marker) ;
						}else if (outer->value(l)==l_False){
							conflict.push(l);
							buildNotCollisionReason(conflict,c.pointSet1,c.pointSet2,false);
							return false;
						}
					}
				}
			}else if (c.exclusiveLit!=lit_Undef && ! h1_over.intersects(h2_over,false)){
				//If the hulls DID intersect inclusively, the we need to check if they also intersect exclusively.
				Lit l = ~c.exclusiveLit;
				if(outer->value(l)==l_True){
					//do nothing
				}else if(outer->value(l)==l_Undef){
					outer->enqueue(l,hulls_not_collision_marker) ;
				}else if (outer->value(l)==l_False){
					conflict.push(l);
					buildNotCollisionReason(conflict,c.pointSet1,c.pointSet2,false);
					return false;
				}
				//If the hulls do not intersect exclusively, they may still intersect inclusively, so we need to check that.
				if (!h1_over.intersects(h2_over,true)){
					Lit l = ~c.inclusiveLit;
					if(outer->value(l)==l_True){
						//do nothing
					}else if(outer->value(l)==l_Undef){
						outer->enqueue(l,hulls_not_collision_marker) ;
					}else if (outer->value(l)==l_False){
						conflict.push(l);
						buildNotCollisionReason(conflict,c.pointSet1,c.pointSet2,true);
						return false;
					}
				}
			}else if (c.exclusiveLit==lit_Undef && ! h1_over.intersects(h2_over,true)){
				Lit l = ~c.inclusiveLit;
				if(outer->value(l)==l_True){
					//do nothing
				}else if(outer->value(l)==l_Undef){
					outer->enqueue(l,hulls_not_collision_marker) ;
				}else if (outer->value(l)==l_False){
					conflict.push(l);
					buildNotCollisionReason(conflict,c.pointSet1,c.pointSet2,true);
					return false;
				}
			}

		}else if (c.exclusiveLit !=lit_Undef){
			if(h1_under.intersects(h2_under,false)){
				Lit l = c.exclusiveLit;
				if(outer->value(l)==l_True){
					//do nothing
				}else if(outer->value(l)==l_Undef){
					outer->enqueue(l,hulls_collision_marker) ;
				}else if (outer->value(l)==l_False){
					conflict.push(l);
					buildCollisionReason(conflict,c.pointSet1,c.pointSet2,false);
					return false;
				}
			}else if (! h1_over.intersects(h2_over,false)){
				Lit l = ~c.exclusiveLit;
				if(outer->value(l)==l_True){
					//do nothing
				}else if(outer->value(l)==l_Undef){
					outer->enqueue(l,hulls_not_collision_marker) ;
				}else if (outer->value(l)==l_False){
					conflict.push(l);
					buildNotCollisionReason(conflict,c.pointSet1,c.pointSet2,false);
					return false;
				}
			}
		}
	}
	return true;
}
template<unsigned int D, class T>
void ConvexHullCollisionDetector<D,T>::buildReason(Lit p, vec<Lit> & reason, CRef marker){
	reason.push(p);

	for(auto & a:collisionLits){
		if(var(a.inclusiveLit)==var(p)){
			buildCollisionReason(reason,a.pointSet1,a.pointSet2,true);
			return;
		}else if (var(a.exclusiveLit)==var(p)){
			buildNotCollisionReason(reason,a.pointSet1,a.pointSet2,false);
			return;
		}
	}
	assert(false);
}
template<unsigned int D, class T>
void ConvexHullCollisionDetector<D,T>::buildCollisionReason(vec<Lit> & conflict,int pointSet1,int pointSet2, bool inclusive){
	if(D==2){
		buildCollisionReason2d(conflict,pointSet1, pointSet2,inclusive);
	}
}

template<unsigned int D, class T>
void ConvexHullCollisionDetector<D,T>::buildNotCollisionReason(vec<Lit> & conflict,int pointSet1,int pointSet2, bool inclusive){
	if(D==2){
		buildNotCollisionReason2d(conflict,pointSet1, pointSet2,inclusive);
	}
}

template<unsigned int D, class T>
void ConvexHullCollisionDetector<D,T>::buildCollisionReason2d(vec<Lit> & conflict,int pointSet1,int pointSet2, bool inclusive){
	//If the two polygons intersect, there are two possible cases (these are not mutually exclusive)
	//1) 1 polygon has at least one point contained in the other polygon. Learn that either that point must be disabled, or at least
	//		one point from its containing triangle must be disabled (as in the point containment theory).
	//2) There exists an edge (or, more generally, a line segment between any two vertices, as opposed to just an edge) from each poygon,
	//		such that the two edges intersect. Learn that one of the end points must be disabled.
	//Note that if we are NOT considering edge points to to collide, then these conditions can miss collisions, and we have to add the following case.
	//3) There exist three vertices that are shared between the two hulls.



	ConvexPolygon<2, T> & h1 =(ConvexPolygon<2, T> &) convexHullDetectors[pointSet1]->getConvexHull(false)->getHull();
	ConvexPolygon<2, T> & h2 = (ConvexPolygon<2, T> &) convexHullDetectors[pointSet2]->getConvexHull(false)->getHull();

	assert(h1.intersects(h2,inclusive));
	assert(h2.intersects(h1,inclusive));
	//so, first, check each edge segment to see if there is an intersecting line.
	//I'm doing this first, on the hypothesis that these line intersection constraints may form better learned clauses.
	//(since cases 1 and 2 overlap, we have to make a choice about which to favour...)

	//it is probably possible to improve on this quadratic time search...
	for(int i = 0; i<h1.size();i++){
		Point<2,T> & prev = h1[i-1];
		Point<2,T> & p = h1[i];
		LineSegment<2,T> edge1(prev,p);
		for(int j = 0;  j<h2.size();j++){
			Point<2,T> & prev2 = h2[j-1];
			Point<2,T> & p2 = h2[j];
			LineSegment<2,T> edge2(prev2,p2);
			if(edge1.intersects(edge2,inclusive)){
				//learn that one of the endpoints of these two intersecting lines must be disabled
				conflict.push(~mkLit(outer->getPointVar(prev.getID())));
				conflict.push(~mkLit(outer->getPointVar(p.getID())));
				conflict.push(~mkLit(outer->getPointVar(prev2.getID())));
				conflict.push(~mkLit(outer->getPointVar(p2.getID())));
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
			triangle.clear();
			//findContainingTriangle2d(h2,h1[i],triangle,inclusive);
			h2.contains(h1[i],triangle,inclusive);
			conflict.push(~mkLit(outer->getPointVar(h1[i].getID())));
			for(auto & p: triangle){
				int id = p.getID();
				Var v = outer->getPointVar(id);
				assert(outer->value(v)==l_True);
				conflict.push(mkLit(v,true));
			}
			return;
		}
	}

	for(int i = 0; i<h2.size();i++){
		if(h1.contains(h2[i],inclusive)){
			static NConvexPolygon<2,T> triangle;
			triangle.clear();
			//findContainingTriangle2d(h1,h2[i],triangle,inclusive);
			h1.contains(h2[i],triangle,inclusive);
			conflict.push(~mkLit(outer->getPointVar(h2[i].getID())));
			for(auto & p: triangle){
				int id = p.getID();
				Var v = outer->getPointVar(id);
				assert(outer->value(v)==l_True);
				conflict.push(mkLit(v,true));
			}
			return;
		}
	}
	assert(!inclusive);
	int count=0;
	//this probably doesn't have to be quadratic
	for(int i = 0; i<h1.size();i++){
		Point<2,T> & p = h1[i];

		for(int j = 0;  j<h2.size();j++){
			Point<2,T> & p2 = h2[j];
			if(p==p2){
				{
					int id = p.getID();
					Var v = outer->getPointVar(id);
					assert(outer->value(v)==l_True);
					conflict.push(mkLit(v,true));
				}
				{
					int id = p2.getID();
					Var v = outer->getPointVar(id);
					assert(outer->value(v)==l_True);
					conflict.push(mkLit(v,true));
				}
				if(count++==3){
					return;
				}
			}
		}
	}
	assert(false);
}

template<unsigned int D, class T>
void ConvexHullCollisionDetector<D,T>::buildNotCollisionReason2d(vec<Lit> & conflict,int pointSet1,int pointSet2, bool inclusive){
	ConvexPolygon<2, T> & h1 = (ConvexPolygon<2, T> &) convexHullDetectors[pointSet1]->getConvexHull(true)->getHull();
	ConvexPolygon<2, T> & h2 = (ConvexPolygon<2, T> &) convexHullDetectors[pointSet2]->getConvexHull(true)->getHull();
	PointSet<D,T> & over1 = over_sets[pointSet1];
	PointSet<D,T> & over2 = over_sets[pointSet2];
	assert(!h1.intersects(h2, inclusive));
	if(h1.size()==0){
		//then the reason that the shapes don't collide is that at least one point in each of them must be enabled
		//can we improve on this?
		for(int i = 0;i<over1.size();i++){
			if(!over1.pointEnabled(i)){
				Lit l = mkLit( outer->getPointVar(over1[i].getID()));
				assert(outer->value(l)==l_False);
				conflict.push(l);
			}
		}
		return;
	}else if (h2.size()==0){
		for(int i = 0;i<over2.size();i++){
			if(!over2.pointEnabled(i)){
				Lit l = mkLit( outer->getPointVar(over2[i].getID()));
				assert(outer->value(l)==l_False);
				conflict.push(l);
			}
		}
		return;
	}else if (h1.size()==1 && h2.size()==1){
		assert(h1[0]!=h2[0]);//identical points are considered to collide
		//can we improve on this?
		for(int i = 0;i<over1.size();i++){
			if(!over1.pointEnabled(i)){
				Lit l = mkLit( outer->getPointVar(over1[i].getID()));
				assert(outer->value(l)==l_False);
				conflict.push(l);
			}
		}
		for(int i = 0;i<over2.size();i++){
			if(!over2.pointEnabled(i)){
				Lit l = mkLit( outer->getPointVar(over2[i].getID()));
				assert(outer->value(l)==l_False);
				conflict.push(l);
			}
		}
		return;
	}

	//the reason that the two convex hulls do NOT collide is that there exists a separating axis between them.
	//Find that axis, and then find all the positions of the disabled points of both polygons along that axis.
	//let h1 be left of h2.
	//Then either some disabled point in h1 that is to the right of the rightmost point in h1 must be enabled,
	//or some disabled point in h2 that is to the left of the leftmost point in h2 must be enabled.

	//Note: It may be possible to improve on this analysis!
	std::vector<std::pair<Point<2,T> ,T>>  projection1;
	std::vector<std::pair<Point<2,T> ,T>>  projection2;
	bool check = findSeparatingAxis(h1,h2,over1,over2, projection1,projection2,inclusive);
	assert(check);
	 T leftmost1 = numeric<T>::infinity();
	 T rightmost1 = -numeric<T>::infinity();

	for(auto & p:projection1){
		int pID = p.first.getID();
		int pointset = outer->getPointset(pID);

		int pointsetIndex = outer->getPointsetIndex(pID);

		if(over1.pointEnabled(pointsetIndex)){
			if(p.second<leftmost1)
				leftmost1 = p.second;
			if(p.second>rightmost1)
				rightmost1 = p.second;
		}
	}

	 T leftmost2 = numeric<T>::infinity();
	 T rightmost2 = -numeric<T>::infinity();

	for(auto & p:projection2){
		int pID = p.first.getID();
		int pointset = outer->getPointset(pID);
		int pointsetIndex = outer->getPointsetIndex(pID);

		if(over2.pointEnabled(pointsetIndex)){
		if(p.second<leftmost2)
			leftmost2 = p.second;
		if(p.second>rightmost2)
			rightmost2 = p.second;
		}
	}


	bool h1_is_left;
	bool h1_is_right;
	if(inclusive){
		assert(rightmost1<leftmost2 || rightmost2<leftmost1);
		h1_is_left=rightmost1<leftmost2;
		h1_is_right=rightmost1>leftmost2;
	}else{
		assert(rightmost1<=leftmost2 || rightmost2<=leftmost1);
		h1_is_left=rightmost1<=leftmost2;
		h1_is_right=rightmost1>=leftmost2;
	}

	bool h2_is_left;
	bool h2_is_right;
	if(inclusive){
		h2_is_left=rightmost1>leftmost2;
		h2_is_right=rightmost1<leftmost2;
	}else{
		h2_is_left=rightmost1>=leftmost2;
		h2_is_right=rightmost1<=leftmost2;
	}

	for(auto & p:projection1){
		int pID = p.first.getID();
		int pointset = outer->getPointset(pID);
		int pointsetIndex = outer->getPointsetIndex(pID);

		if(!over1.pointEnabled(pointsetIndex)){
			if(h1_is_left && p.second>rightmost1){
				//then if we enable this point, the two hulls will move closer to each other.
				//we can probably improve on this, by only considering points that are also >= the rightmost _disabled_ point of pointset2...
				Lit l = mkLit(outer->getPointVar(pID));
				assert(outer->value(l)==l_False);
				conflict.push(l);
			}
			if(h1_is_right && p.second<leftmost1){
				Lit l = mkLit(outer->getPointVar(pID));
				assert(outer->value(l)==l_False);
				conflict.push(l);
			}
		}
	}
	for(auto & p:projection2){
		int pID = p.first.getID();
		int pointset = outer->getPointset(pID);
		int pointsetIndex = outer->getPointsetIndex(pID);

		if(!over2.pointEnabled(pointsetIndex)){
			if(h2_is_left && p.second>rightmost2){
				//then if we enable this point, the two hulls will move closer to each other.
				//we can probably improve on this, by only considering points that are also >= the rightmost _disabled_ point of pointset2...
				Lit l = mkLit(outer->getPointVar(pID));
				assert(outer->value(l)==l_False);
				conflict.push(l);
			}
			if(h2_is_right && p.second<leftmost2){
				Lit l = mkLit(outer->getPointVar(pID));
				assert(outer->value(l)==l_False);
				conflict.push(l);
			}
		}
	}
}
template<unsigned int D, class T>
bool ConvexHullCollisionDetector<D, T>::findSeparatingAxis2d(ConvexPolygon<2, T> & hull1, ConvexPolygon<2, T> & hull2, PointSet<2,T> & pointset1_t, PointSet<2,T> & pointset2_t, std::vector<std::pair<Point<2,T> ,T>>  &projection_out1_t,std::vector<std::pair<Point<2,T> ,T>>  &projection_out2_t, bool inclusive){
	ConvexPolygon<2, T> & h1 = (hull1.size()<=hull2.size())?hull1:hull2;
	ConvexPolygon<2, T> & h2 = (hull1.size()<=hull2.size())?hull2:hull1;
	std::vector<std::pair<Point<2,T> ,T>>  &projection_out1 = (hull1.size()<=hull2.size())?projection_out1_t:projection_out2_t;
	std::vector<std::pair<Point<2,T> ,T>>  &projection_out2 = (hull1.size()<=hull2.size())?projection_out2_t:projection_out1_t;
	PointSet<2,T> & pointset1 =  (hull1.size()<=hull2.size())?pointset1_t:pointset2_t;
	PointSet<2,T> & pointset2 =  (hull1.size()<=hull2.size())?pointset2_t:pointset1_t;


	projection_out1.clear();
	projection_out2.clear();
	if(h1.size()==0 || h2.size()==0){
		return false;
	}else if (h1.size()==1 && h2.size()==1){
		 Point<2,T> un_normalized_normal = h2[0]-h1[0];
		 //now place all the remaining (disabled) points on the axis as well
		 for(int i = 0;i<pointset1.size();i++){
				 auto & p = pointset1[i];
				 T projection = un_normalized_normal.dot(p);
				 projection_out1.push_back({p,projection});
		 }
		 for(int i = 0;i<pointset2.size();i++){
				 auto & p = pointset2[i];
				 T projection = un_normalized_normal.dot(p);
				 projection_out2.push_back({p,projection});
		 }
		return true;
	}


	 //Separating Axis Theorem for collision detection between two convex polygons
	 //loop through each edge in _each_ polygon and project both polygons onto that edge's normal.
	 //If any of the projections are non-intersection, then these don't collide; else, they do collide
	 if(h1.size()>1){
		for(int i = 0;i<h1.size();i++){
		 auto & p = h1[i];
		 auto & prev = (h1)[i-1];
		 Point<2,T> edge = p-prev;
		 Point<2,T> un_normalized_normal(-edge.y, edge.x);
		 projection_out1.clear();
		 projection_out2.clear();
		 //now project both polygons onto to this normal and see if they overlap, by finding the minimum and maximum distances
		 //Note that since we are NOT normalizing the normal vector, the projection is distorted along that vector
		 //(this still allows us to check overlaps, but means that the minimum distance found between the two shapes may be incorrect)
		 T left = numeric<T>::infinity();
		 T right = -numeric<T>::infinity();
		 for (auto & p:h1){
			 T projection = un_normalized_normal.dot(p);
			 projection_out1.push_back({p,projection});
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
		 for (auto & p:h2){
				 T projection = un_normalized_normal.dot(p);
				 projection_out2.push_back({p,projection});
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
			 //now place all the remaining (disabled) points on the axis as well
			 for(int i = 0;i<pointset1.size();i++){
				 if(!pointset1.pointEnabled(i)){
					 auto & p = pointset1[i];
					 T projection = un_normalized_normal.dot(p);
					 projection_out1.push_back({p,projection});
				 }
			 }
			 for(int i = 0;i<pointset2.size();i++){
				 if(!pointset2.pointEnabled(i)){
					 auto & p = pointset2[i];
					 T projection = un_normalized_normal.dot(p);
					 projection_out2.push_back({p,projection});
				 }
			 }
			 return true;
		 }
	 }
	}
	if(h2.size()>1){
	 //now test the axis produced by the other polygon
	 for(int i = 0;i<h2.size();i++){
		 auto & p = h2[i];
		 auto & prev = h2[i-1];
		 Point<2,T> edge = p-prev;
		 Point<2,T> un_normalized_normal(-edge.y, edge.x);
		 projection_out1.clear();
		 projection_out2.clear();
		 T left = numeric<T>::infinity();
		 T right = -numeric<T>::infinity();
		 for (auto & p:h2){
			 T projection = un_normalized_normal.dot(p);
			 projection_out2.push_back({p,projection});
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
		 for (auto & p:h1){
			 T projection = un_normalized_normal.dot(p);
			 projection_out1.push_back({p,projection});
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
			 for(int i = 0;i<pointset1.size();i++){
				 if(!pointset1.pointEnabled(i)){
					 auto & p = pointset1[i];
					 T projection = un_normalized_normal.dot(p);
					 projection_out1.push_back({p,projection});
				 }
			 }
			 for(int i = 0;i<pointset2.size();i++){
				 if(!pointset2.pointEnabled(i)){
					 auto & p = pointset2[i];
					 T projection = un_normalized_normal.dot(p);
					 projection_out2.push_back({p,projection});
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


template<unsigned int D, class T>
bool ConvexHullCollisionDetector<D,T>::checkSatisfied(){
	for(auto & c:collisionLits){
		//this is really ugly.
		auto & h1_under = convexHullDetectors[c.pointSet1]->getConvexHull(false)->getHull();
		auto & h1_over =  convexHullDetectors[c.pointSet1]->getConvexHull(true)->getHull();
		auto & h2_under = convexHullDetectors[c.pointSet2]->getConvexHull(false)->getHull();
		auto & h2_over = convexHullDetectors[c.pointSet2]->getConvexHull(true)->getHull();

		if (c.inclusiveLit!=lit_Undef){
			if (outer->value(c.inclusiveLit)==l_True){
				if(! h1_over.intersects(h2_over,true)){
					return false;
				}
			}else if (outer->value(c.inclusiveLit)==l_False){
				if(h1_under.intersects(h2_under,true)){
					cout<<h1_over << "\n";
					cout << h2_over<<"\n";
					return false;
				}
			}
		}
		if (c.exclusiveLit!=lit_Undef){
			if (outer->value(c.exclusiveLit)==l_True){
				if(! h1_over.intersects(h2_over,false)){
					cout<<h1_under << "\n";
					cout << h2_under<<"\n";
					return false;
				}
			}else if (outer->value(c.exclusiveLit)==l_False){
				if(h1_under.intersects(h2_under,false)){
					return false;
				}
			}
		}
	}
	return true;
}


#endif
