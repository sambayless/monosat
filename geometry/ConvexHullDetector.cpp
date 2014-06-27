#include "ConvexHullDetector.h"
#include "GeometryTheory.h"

template<>
void ConvexHullDetector<2,double>::buildPointContainedReason(const Point<2,double> & s,vec<Lit> & conflict){
	//A point is contained because there exists three enabled points that form a containing triangle
	//(might have to handle the edge case of a 2 point hull depending on how we interpret such a hull)
	//its not immediately clear how we would pick which 3 points to include, given that there may be a choice
	ConvexPolygon<2,double> &hull = under_hull->getHull();
	assert(hull.contains(s));


	ConvexPolygon<2,double> triangle;
	findContainingTriangle2d(hull,s,triangle);
	assert(triangle.contains(s));
	for(auto & p: triangle){
		int id = p.getID();
		Var v = outer->getPointVar(id);
		assert(outer->value(v)==l_True);
		conflict.push(mkLit(v,true));
	}

}

template<>
void ConvexHullDetector<2,double>::buildPointNotContainedReason(const Point<2,double> & s, vec<Lit> & conflict){
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

		Line<2,double> testline;


		testline.a = s;
		ConvexPolygon<2,double> &hull = over_hull->getHull();

		//edge cases:
		if(hull.size()<3){
			//then we just report that _some_ disabled vertex must be enabled
			for(int i = 0;i<over.size();i++){
				Point<2,double> & p = over[i];
				if(!over.pointEnabled(i)){
					conflict.push(mkLit( outer->getPointVar(p.getID()),false));
				}
			}
			return;
		}

		Line<2,double> top_line;
		Line<2,double> bottom_line;
		top_line.a = s;
		bottom_line.a=s;
		top_line.b = hull.getVertices()[0];
		bottom_line.b = hull.getVertices().back();

		std::vector<Point<2,double> > min_set;
		std::vector<Point<2,double> > test_set;
		hasSet=false;

		for(int i = 0;i<over.size();i++){
			Point<2,double> & p = over[i];
			if(!over.pointEnabled(i)){
				if(!hull.contains(p)){
					testline.b = p;

					//now check to see if test line intersects the hull (it can still intersect on the otherside
					if(!testline.intersects(hull)){
						//then this is a point outside the hull that we can try testing.
						findFarPoints(testline,true,test_set,min_set,hull);
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
		findFarPoints(top_line,true,test_set,min_set,hull);
		findFarPoints(bottom_line,true,test_set,min_set,hull);

		for(int i = 0;i<min_set.size();i++){
			Point<2,double> & p = min_set[i];
			assert(!over.pointEnabled(p.getID()));
			Lit l = mkLit( outer->getPointVar(p.getID()),false);
			assert(outer->value(l)==l_False);
			conflict.push(l);

		}

}

template<>
void ConvexHullDetector<2,double>::buildPointOnHullOrDisabledReason(Var pointVar,const Point<2,double> & s, vec<Lit> & conflict){
	//If the point IS on the hull (or disabled), it is either because it is currently disabled...

	assert(outer->getPointVar(point.getID())==pointVar);
	Lit enabledLit = mkLit(pointVar,false);

	if(outer->value(enabledLit)==l_False){
		conflict.push(enabledLit);
		return;
	}

	//or because there is some line we can form that passes through this point, such that no enabled point is on the far side of the line, and the entire hull (except for this point)
	//is on the near side.
	//Note that as we are forcing the convex hull to contain all collinear points on the boundary,
	//this code is handles the case where enabled points are exactly on this line in the opposite way as  'buildPointNotContainedReason' does


	Line<2,double> testline;

	testline.a = s;
	ConvexPolygon<2,double> &hull = over_hull->getHull();

	//edge cases:
	if(hull.size()<3){
		//then we just report that _some_ disabled vertex must be enabled
		for(int i = 0;i<over.size();i++){
			Point<2,double> & p = over[i];
			if(!over.pointEnabled(i)){
				conflict.push(mkLit( outer->getPointVar(p.getID()),false));
			}
		}
		return;
	}

		Line<2,double> top_line;
		Line<2,double> bottom_line;
		top_line.a = s;
		bottom_line.a=s;
		top_line.b = hull.getVertices()[0];
		bottom_line.b = hull.getVertices().back();

		std::vector<Point<2,double> > min_set;
		std::vector<Point<2,double> > test_set;
		hasSet=false;

		for(int i = 0;i<over.size();i++){
			Point<2,double> & p = over[i];
			if(!over.pointEnabled(i)){
				if(!hull.contains(p)){
					testline.b = p;

					//now check to see if test line intersects the hull (it can still intersect on the otherside
					if(!testline.intersects(hull)){
						//then this is a point outside the hull that we can try testing.
						//Note that as we are forcing the convex hull to contain all collinear points on the boundary,
						//this code is handles the case where enabled points are exactly on this line in the opposite way as  'buildPointNotContainedReason' does

						findFarPoints(testline,false,test_set,min_set,hull);
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
		findFarPoints(top_line,true,test_set,min_set,hull);
		findFarPoints(bottom_line,true,test_set,min_set,hull);

		for(int i = 0;i<min_set.size();i++){
			Point<2,double> & p = min_set[i];
			assert(!over.pointEnabled(p.getID()));
			Lit l = mkLit( outer->getPointVar(p.getID()),false);
			assert(outer->value(l)==l_False);
			conflict.push(l);

		}


}

template< >
void ConvexHullDetector<2,double>::buildPointNotOnHullOrDisabledReason(Var pointVar,const Point<2,double> & s, vec<Lit> & conflict){
	//If the point is not (on the hull or disabled), it is because it is both currently enabled...
	assert(outer->getPointVar(s.getID())==pointVar);
	Lit enabledLit = mkLit(pointVar,false);
	assert(outer->getValue(enabledLit)==l_True);
	conflict.push(~enabledLit);
	//AND, it is because there exist three other enabled points that form a triangle that contain this point:
	ConvexPolygon<2,double> &hull = under_hull->getHull();
#ifndef NDEBUG
	for(auto & p:hull)
		assert(p.getID()!=s.getID());
#endif
	ConvexPolygon<2,double> triangle;
	findContainingTriangle2d(hull,s,triangle);
	assert(triangle.contains(s));
	for(auto & p: triangle){
		int id = p.getID();
		Var v = outer->getPointVar(id);
		assert(outer->value(v)==l_True);
		conflict.push(mkLit(v,true));
	}


}

template<>
void ConvexHullDetector<2,mpq_class>::buildPointContainedReason(const Point<2,mpq_class> & s,vec<Lit> & conflict){
	assert(false);
}

template<>
void ConvexHullDetector<2,mpq_class>::buildPointNotContainedReason(const Point<2,mpq_class> & s, vec<Lit> & conflict){
	assert(false);
}

template<>
void ConvexHullDetector<2,mpq_class>::buildPointOnHullOrDisabledReason(Var pointVar,const Point<2,mpq_class> & s, vec<Lit> & conflict){
	assert(false);
}

template< >
void ConvexHullDetector<2,mpq_class>::buildPointNotOnHullOrDisabledReason(Var pointVar,const Point<2,mpq_class> & s, vec<Lit> & conflict){
	assert(false);
}
