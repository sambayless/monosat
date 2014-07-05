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

	assert(outer->getPointVar(s.getID())==pointVar);
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
	assert(outer->value(enabledLit)==l_True);
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
void ConvexHullDetector<2,double>::buildConvexIntersectsReason(ConvexPolygon<2,double> & polygon,vec<Lit> & conflict){
	//If the two polygons intersect, there are two possible cases (these are not mutually exclusive)
	//1) 1 polygon has at least one point contained in the other polygon. Learn that either that point must be disabled, or at least
	//		one point from its containing triangle must be disabled (as in the point containment theory).
	//2) There exists an edge (or, more generally, a line segment between any two vertices, as opposed to just an edge) from each poygon,
	//		such that the two edges intersect. Learn that one of the end points must be disabled.
	ConvexPolygon<2, double> & h1 = under_hull->getHull();
	ConvexPolygon<2, double> & h2 = polygon;


	//so, first, check each edge segment to see if there is an intersecting line.
	//I'm doing this first, on the hypothesis that these line intersection constraints may form better learned clauses.
	//(since cases 1 and 2 overlap, we have to make a choice about which to favour...)

	//it is probably possible to improve on this quadratic time search...
	for(int i = 0; i<h1.size();i++){
		Point<2,double> & prev = h1[i-1];
		Point<2,double> & p = h1[i];
		LineSegment<2,double> edge1(prev,p);
		for(int j = 0;  j<h2.size();j++){
			Point<2,double> & prev2 = h1[i-1];
			Point<2,double> & p2 = h1[i];
			LineSegment<2,double> edge2(prev2,p2);
			if(edge1.intersects(edge2)){
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
		if(h2.contains(h1[i])){
			ConvexPolygon<2,double> triangle;
			findContainingTriangle2d(h2,h1[i],triangle);
			assert(triangle.contains(h1[i]));
			conflict.push(~mkLit(outer->getPointVar(h1[i].getID())));
			return;
		}
	}

	for(int i = 0; i<h2.size();i++){
		if(h1.contains(h2[i])){
			ConvexPolygon<2,double> triangle;
			findContainingTriangle2d(h1,h2[i],triangle);
			assert(triangle.contains(h2[i]));
			for(auto & p: triangle){
				int id = p.getID();
				Var v = outer->getPointVar(id);
				assert(outer->value(v)==l_True);
				conflict.push(mkLit(v,true));
			}
			return;
		}
	}
}

template<>
void ConvexHullDetector<2,double>::buildConvexNotIntersectsReason(ConvexPolygon<2,double> & polygon,vec<Lit> & conflict){
	ConvexPolygon<2, double> & h1 = over_hull->getHull();
	ConvexPolygon<2, double> & h2 = polygon;


	if(h1.size()==0){
		//then the reason that the shapes don't collide is that at least one point in each of them must be enabled
		//can we improve on this?
		for(int i = 0;i<h1.size();i++){
			if(!over.pointEnabled(i)){
				Lit l = mkLit( outer->getPointVar(over[i].getID()));
				assert(outer->value(l)==l_False);
				conflict.push(l);
			}
		}
	}else if (h2.size()==0){
		return;
	}else if (h2.size()==1){
		buildPointNotContainedReason(h2[0], conflict);
	}else if (h1.size()==1 ){
		assert(h1[0]!=h2[0]);//identical points are considered to collide
	}

	//the reason that the two convex hulls do NOT collide is that there exists a separating axis between them.
	//Find that axis, and then find all the positions of the disabled points of both polygons along that axis.
	//let h1 be left of h2.
	//Then either some disabled point in h1 that is to the right of the rightmost point in h1 must be enabled,
	//or some disabled point in h2 that is to the left of the leftmost point in h2 must be enabled.

	//Note: It may be possible to improve on this analysis!
	vec<std::pair<Point<2,double> ,double>>  projection;
	vec<std::pair<Point<2,double> ,double>>  projection2;
	findSeparatingAxis(h1,h2,over, projection,projection2);

	double leftmost1 = std::numeric_limits<double>::max();
	 double rightmost1 = std::numeric_limits<double>::min();

	for(auto & p:projection){
		int pID = p.first.getID();
		int pointset = outer->getPointset(pID);
		assert(pointset==over.getID());
		int pointsetIndex = outer->getPointsetIndex(pID);
		if(over.pointEnabled(pointsetIndex)){
			if(p.second<leftmost1)
				leftmost1 = p.second;
			if(p.second>rightmost1)
				rightmost1 = p.second;
		}
	}

	double leftmost2 = std::numeric_limits<double>::max();
	double rightmost2 = std::numeric_limits<double>::min();

	for(auto & p:projection2){
		int pID = p.first.getID();
		int pointset = outer->getPointset(pID);
		int pointsetIndex = outer->getPointsetIndex(pID);


		if(p.second<leftmost2)
			leftmost2 = p.second;
		if(p.second>rightmost2)
			rightmost2 = p.second;

	}

	assert(rightmost1<leftmost2 || rightmost2<leftmost1);
	bool h1_is_left=rightmost1<leftmost2;

	for(auto & p:projection){
		int pID = p.first.getID();
		int pointset = outer->getPointset(pID);
		int pointsetIndex = outer->getPointsetIndex(pID);
		assert(pointset==over.getID());
		if(!over.pointEnabled(pointsetIndex)){
			if(h1_is_left ? (p.second>rightmost1):((p.second<leftmost1))){
				//then if we enable this point, the two hulls will move closer to each other.
				//we can probably improve on this, by only considering points that are also >= the rightmost _disabled_ point of pointset2...

				Lit l = mkLit(outer->getPointVar(pID));
				assert(outer->value(l)==l_False);
				conflict.push(l);
			}
		}
	}

}

template<>
void ConvexHullDetector<2,mpq_class>::buildPointContainedReason(const Point<2,mpq_class> & s,vec<Lit> & conflict){
	//A point is contained because there exists three enabled points that form a containing triangle
	//(might have to handle the edge case of a 2 point hull depending on how we interpret such a hull)
	//its not immediately clear how we would pick which 3 points to include, given that there may be a choice
	ConvexPolygon<2,mpq_class> &hull = under_hull->getHull();
	assert(hull.contains(s));


	ConvexPolygon<2,mpq_class> triangle;
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
void ConvexHullDetector<2,mpq_class>::buildPointNotContainedReason(const Point<2,mpq_class> & s, vec<Lit> & conflict){
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

			Line<2,mpq_class> testline;


			testline.a = s;
			ConvexPolygon<2,mpq_class> &hull = over_hull->getHull();

			//edge cases:
			if(hull.size()<3){
				//then we just report that _some_ disabled vertex must be enabled
				for(int i = 0;i<over.size();i++){
					Point<2,mpq_class> & p = over[i];
					if(!over.pointEnabled(i)){
						conflict.push(mkLit( outer->getPointVar(p.getID()),false));
					}
				}
				return;
			}

			Line<2,mpq_class> top_line;
			Line<2,mpq_class> bottom_line;
			top_line.a = s;
			bottom_line.a=s;
			top_line.b = hull.getVertices()[0];
			bottom_line.b = hull.getVertices().back();

			std::vector<Point<2,mpq_class> > min_set;
			std::vector<Point<2,mpq_class> > test_set;
			hasSet=false;

			for(int i = 0;i<over.size();i++){
				Point<2,mpq_class> & p = over[i];
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
				Point<2,mpq_class> & p = min_set[i];
				assert(!over.pointEnabled(p.getID()));
				Lit l = mkLit( outer->getPointVar(p.getID()),false);
				assert(outer->value(l)==l_False);
				conflict.push(l);

			}
}

/*
template<>
void ConvexHullDetector<2,mpq_class>::buildLineIntersectsReason(LineSegment<2,mpq_class> & s, vec<Lit> & conflict){
	//If the two polygons intersect, there are two possible cases (these are not mutually exclusive)
	//1) 1 polygon has at least one point contained in the other polygon. Learn that either that point must be disabled, or at least
	//		one point from its containing triangle must be disabled (as in the point containment theory).
	//2) There exists an edge (or, more generally, a line segment between any two vertices, as opposed to just an edge) from each poygon,
	//		such that the two edges intersect. Learn that one of the end points must be disabled.
	ConvexPolygon<2, mpq_class> & hull = under_hull->getHull();

	//First, check to see if either end point is contained within the convex hull
	if(hull.contains( s.a)){
		//find a containing triangle around one of these points
		ConvexPolygon<2,mpq_class> triangle;
		findContainingTriangle2d(hull,s.a,triangle);
		assert(triangle.contains(s.a));
		for(auto & p: triangle){
			int id = p.getID();
			Var v = outer->getPointVar(id);
			assert(outer->value(v)==l_True);
			conflict.push(mkLit(v,true));
		}
	}else if  (hull.contains(s.b)){
		ConvexPolygon<2,mpq_class> triangle;
		findContainingTriangle2d(hull,s.b,triangle);
		assert(triangle.contains(s.b));
		for(auto & p: triangle){
			int id = p.getID();
			Var v = outer->getPointVar(id);
			assert(outer->value(v)==l_True);
			conflict.push(mkLit(v,true));
		}
	}

	//If neither endpoint is contained, then it may still be the case that the line segment intersects the hull (if it passes through at least one edge of the hull).
	//So, iterate through the edges and find one that intersects the line

	for(int i = 0; i<hull.size();i++){
		Point<2,mpq_class> & prev = hull[i-1];
		Point<2,mpq_class> & p = hull[i];
		LineSegment<2,mpq_class> edge(prev,p);

		if(edge.intersects(s)){
			//learn that one of the endpoints of the intersecting line must be disabled
			conflict.push(~mkLit(outer->getPointVar(prev.getID())));
			conflict.push(~mkLit(outer->getPointVar(p.getID())));

			return;
		}

	}
}
template<>
void ConvexHullDetector<2,mpq_class>::buildLineNotIntersectsReason( LineSegment<2,mpq_class> & s, vec<Lit> & conflict){
		ConvexPolygon<2, mpq_class> & hull = under_hull->getHull();


		if(hull.size()==0){
			//then the reason that the shapes don't collide is that at least one of them must be enabled
			//can we improve on this?
			for(int i = 0;i<over.size();i++){
				if(!over.pointEnabled(i)){
					Lit l = mkLit( outer->getPointVar(over[i].getID()));
					assert(outer->value(l)==l_False);
					conflict.push(l);
				}
			}
		}else if (hull.size()==1){
			//can we improve on this?
			for(int i = 0;i<over.size();i++){
				if(!over.pointEnabled(i)){
					Lit l = mkLit( outer->getPointVar(over[i].getID()));
					assert(outer->value(l)==l_False);
					conflict.push(l);
				}
			}
		}

		//the reason that the two convex hulls do NOT collide is that there exists a separating axis between them.
		//Find that axis, and then find all the positions of the disabled points of both polygons along that axis.
		//let hull be left of h2.
		//Then either some disabled point in hull that is to the right of the rightmost point in hull must be enabled,
		//or some disabled point in h2 that is to the left of the leftmost point in h2 must be enabled.

		//Note: It may be possible to improve on this analysis!
		vec<std::pair<Point<2,mpq_class> ,mpq_class>>  projection1;
		findSeparatingAxis(hull,s,over,projection1);

		 mpq_class leftmost1 = std::numeric_limits<mpq_class>::max();
		 mpq_class rightmost1 = std::numeric_limits<mpq_class>::min();

		for(auto & p:projection1){
			int pID = p.first.getID();
			int pointset = outer->getPointset(pID);
			assert(pointset==over1.getID());
			int pointsetIndex = outer->getPointsetIndex(pID);
			if(over1.pointEnabled(pointsetIndex)){
				if(p.second<leftmost1)
					leftmost1 = p.second;
				if(p.second>rightmost1)
					rightmost1 = p.second;
			}
		}

		 mpq_class leftmost2 = std::numeric_limits<mpq_class>::max();
		 mpq_class rightmost2 = std::numeric_limits<mpq_class>::min();

		for(auto & p:projection2){
			int pID = p.first.getID();
			int pointset = outer->getPointset(pID);
			int pointsetIndex = outer->getPointsetIndex(pID);
			assert(pointset==over2.getID());
			if(over2.pointEnabled(pointsetIndex)){
			if(p.second<leftmost2)
				leftmost2 = p.second;
			if(p.second>rightmost2)
				rightmost2 = p.second;
			}
		}

		assert(rightmost1<leftmost2 || rightmost2<leftmost1);
		bool hull_is_left=rightmost1<leftmost2;

		for(auto & p:projection1){
			int pID = p.first.getID();
			int pointset = outer->getPointset(pID);
			int pointsetIndex = outer->getPointsetIndex(pID);
			assert(pointset==over1.getID());
			if(!over1.pointEnabled(pointsetIndex)){
				if(hull_is_left ? (p.second>rightmost1):((p.second<leftmost1))){
					//then if we enable this point, the two hulls will move closer to each other.
					//we can probably improve on this, by only considering points that are also >= the rightmost _disabled_ point of pointset2...

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
			assert(pointset==over2.getID());
			if(!over2.pointEnabled(pointsetIndex)){
				if(hull_is_left ? (p.second<leftmost2):((p.second>rightmost2))){
					//then if we enable this point, the two hulls will move closer to each other.
					//we can probably improve on this, by only considering points that are also >= the rightmost _disabled_ point of pointset2...
					Lit l = mkLit(outer->getPointVar(pID));
					assert(outer->value(l)==l_False);
					conflict.push(l);
				}
			}
		}

};
template<>
void ConvexHullDetector<2,double>::buildLineIntersectsReason( LineSegment<2,double> & s, vec<Lit> & conflict){
	//If the two polygons intersect, there are two possible cases (these are not mutually exclusive)
	//1) 1 polygon has at least one point contained in the other polygon. Learn that either that point must be disabled, or at least
	//		one point from its containing triangle must be disabled (as in the point containment theory).
	//2) There exists an edge (or, more generally, a line segment between any two vertices, as opposed to just an edge) from each poygon,
	//		such that the two edges intersect. Learn that one of the end points must be disabled.
	ConvexPolygon<2, double> & hull = under_hull->getHull();

	//First, check to see if either end point is contained within the convex hull
	if(hull.contains( s.a)){
		//find a containing triangle around one of these points
		ConvexPolygon<2,double> triangle;
		findContainingTriangle2d(hull,s.a,triangle);
		assert(triangle.contains(s.a));
		for(auto & p: triangle){
			int id = p.getID();
			Var v = outer->getPointVar(id);
			assert(outer->value(v)==l_True);
			conflict.push(mkLit(v,true));
		}
	}else if  (hull.contains(s.b)){
		ConvexPolygon<2,double> triangle;
		findContainingTriangle2d(hull,s.b,triangle);
		assert(triangle.contains(s.b));
		for(auto & p: triangle){
			int id = p.getID();
			Var v = outer->getPointVar(id);
			assert(outer->value(v)==l_True);
			conflict.push(mkLit(v,true));
		}
	}

	//If neither endpoint is contained, then it may still be the case that the line segment intersects the hull (if it passes through at least one edge of the hull).
	//So, iterate through the edges and find one that intersects the line

	for(int i = 0; i<hull.size();i++){
		Point<2,double> & prev = hull[i-1];
		Point<2,double> & p = hull[i];
		LineSegment<2,double> edge(prev,p);

		if(edge.intersects(s)){
			//learn that one of the endpoints of the intersecting line must be disabled
			conflict.push(~mkLit(outer->getPointVar(prev.getID())));
			conflict.push(~mkLit(outer->getPointVar(p.getID())));

			return;
		}

	}
}
template<>
void ConvexHullDetector<2,double>::buildLineNotIntersectsReason( LineSegment<2,double> & s, vec<Lit> & conflict){

}
*/


template<>
void ConvexHullDetector<2,mpq_class>::buildConvexIntersectsReason(ConvexPolygon<2,mpq_class> & polygon,vec<Lit> & conflict){
	//If the two polygons intersect, there are two possible cases (these are not mutually exclusive)
	//1) 1 polygon has at least one point contained in the other polygon. Learn that either that point must be disabled, or at least
	//		one point from its containing triangle must be disabled (as in the point containment theory).
	//2) There exists an edge (or, more generally, a line segment between any two vertices, as opposed to just an edge) from each poygon,
	//		such that the two edges intersect. Learn that one of the end points must be disabled.
	ConvexPolygon<2, mpq_class> & h1 = under_hull->getHull();
	ConvexPolygon<2, mpq_class> & h2 = polygon;


	//so, first, check each edge segment to see if there is an intersecting line.
	//I'm doing this first, on the hypothesis that these line intersection constraints may form better learned clauses.
	//(since cases 1 and 2 overlap, we have to make a choice about which to favour...)

	//it is probably possible to improve on this quadratic time search...
	for(int i = 0; i<h1.size();i++){
		Point<2,mpq_class> & prev = h1[i-1];
		Point<2,mpq_class> & p = h1[i];
		LineSegment<2,mpq_class> edge1(prev,p);
		for(int j = 0;  j<h2.size();j++){
			Point<2,mpq_class> & prev2 = h1[i-1];
			Point<2,mpq_class> & p2 = h1[i];
			LineSegment<2,mpq_class> edge2(prev2,p2);
			if(edge1.intersects(edge2)){
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
		if(h2.contains(h1[i])){
			ConvexPolygon<2,mpq_class> triangle;
			findContainingTriangle2d(h2,h1[i],triangle);
			assert(triangle.contains(h1[i]));
			conflict.push(~mkLit(outer->getPointVar(h1[i].getID())));
			return;
		}
	}

	for(int i = 0; i<h2.size();i++){
		if(h1.contains(h2[i])){
			ConvexPolygon<2,mpq_class> triangle;
			findContainingTriangle2d(h1,h2[i],triangle);
			assert(triangle.contains(h2[i]));
			for(auto & p: triangle){
				int id = p.getID();
				Var v = outer->getPointVar(id);
				assert(outer->value(v)==l_True);
				conflict.push(mkLit(v,true));
			}
			return;
		}
	}
}

template<>
void ConvexHullDetector<2,mpq_class>::buildConvexNotIntersectsReason(ConvexPolygon<2,mpq_class> & polygon,vec<Lit> & conflict){
	ConvexPolygon<2, mpq_class> & h1 = over_hull->getHull();
	ConvexPolygon<2, mpq_class> & h2 = polygon;


	if(h1.size()==0){
		//then the reason that the shapes don't collide is that at least one point in each of them must be enabled
		//can we improve on this?
		for(int i = 0;i<h1.size();i++){
			if(!over.pointEnabled(i)){
				Lit l = mkLit( outer->getPointVar(over[i].getID()));
				assert(outer->value(l)==l_False);
				conflict.push(l);
			}
		}
	}else if (h2.size()==0){
		return;
	}else if (h2.size()==1){
		buildPointNotContainedReason(h2[0], conflict);
	}else if (h1.size()==1 ){
		assert(h1[0]!=h2[0]);//identical points are considered to collide
	}

	//the reason that the two convex hulls do NOT collide is that there exists a separating axis between them.
	//Find that axis, and then find all the positions of the disabled points of both polygons along that axis.
	//let h1 be left of h2.
	//Then either some disabled point in h1 that is to the right of the rightmost point in h1 must be enabled,
	//or some disabled point in h2 that is to the left of the leftmost point in h2 must be enabled.

	//Note: It may be possible to improve on this analysis!
	vec<std::pair<Point<2,mpq_class> ,mpq_class>>  projection;
	vec<std::pair<Point<2,mpq_class> ,mpq_class>>  projection2;
	findSeparatingAxis(h1,h2,over, projection,projection2);

	 mpq_class leftmost1 = std::numeric_limits<mpq_class>::max();
	 mpq_class rightmost1 = std::numeric_limits<mpq_class>::min();

	for(auto & p:projection){
		int pID = p.first.getID();
		int pointset = outer->getPointset(pID);
		assert(pointset==over.getID());
		int pointsetIndex = outer->getPointsetIndex(pID);
		if(over.pointEnabled(pointsetIndex)){
			if(p.second<leftmost1)
				leftmost1 = p.second;
			if(p.second>rightmost1)
				rightmost1 = p.second;
		}
	}

	 mpq_class leftmost2 = std::numeric_limits<mpq_class>::max();
	 mpq_class rightmost2 = std::numeric_limits<mpq_class>::min();

	for(auto & p:projection2){
		int pID = p.first.getID();
		int pointset = outer->getPointset(pID);
		int pointsetIndex = outer->getPointsetIndex(pID);


		if(p.second<leftmost2)
			leftmost2 = p.second;
		if(p.second>rightmost2)
			rightmost2 = p.second;

	}

	assert(rightmost1<leftmost2 || rightmost2<leftmost1);
	bool h1_is_left=rightmost1<leftmost2;

	for(auto & p:projection){
		int pID = p.first.getID();
		int pointset = outer->getPointset(pID);
		int pointsetIndex = outer->getPointsetIndex(pID);
		assert(pointset==over.getID());
		if(!over.pointEnabled(pointsetIndex)){
			if(h1_is_left ? (p.second>rightmost1):((p.second<leftmost1))){
				//then if we enable this point, the two hulls will move closer to each other.
				//we can probably improve on this, by only considering points that are also >= the rightmost _disabled_ point of pointset2...

				Lit l = mkLit(outer->getPointVar(pID));
				assert(outer->value(l)==l_False);
				conflict.push(l);
			}
		}
	}

}

template<>
bool ConvexHullDetector<2, mpq_class>::findSeparatingAxis(ConvexPolygon<2, mpq_class> & hull1, ConvexPolygon<2, mpq_class> & hull2, PointSet<2,mpq_class> & pointset1, vec<std::pair<Point<2,mpq_class> ,mpq_class>>  &projection_out1, vec<std::pair<Point<2,mpq_class> ,mpq_class>>  &projection_out2){
	ConvexPolygon<2, mpq_class> & h1 = (hull1.size()<=hull2.size())?hull1:hull2;
	ConvexPolygon<2, mpq_class> & h2 = (hull1.size()<=hull2.size())?hull2:hull1;

	projection_out1.clear();
	projection_out2.clear();
	if(h1.size()==0 || h2.size()==0){
		return false;
	}else if (h1.size()==1 && h2.size()==1){
		 Point<2,mpq_class> un_normalized_normal = h2[0]-h1[0];
		 //now place all the remaining (disabled) points on the axis as well
		 for(int i = 0;i<pointset1.size();i++){
			 auto & p = pointset1[i];
			 mpq_class projection = un_normalized_normal.dot(p);
			 projection_out1.push({p,projection});
		 }
		 for(int i = 0;i<h2.size();i++){
			 auto & p = h2[i];
			 mpq_class projection = un_normalized_normal.dot(p);
			 projection_out2.push({p,projection});
		 }
		return true;
	}


	 std::vector<Point<2,mpq_class> > &  w = h1.getVertices();

	 //Separating Axis Theorem for collision detection between two convex polygons
	 //loop through each edge in _each_ polygon and project both polygons onto that edge's normal.
	 //If any of the projections are non-intersection, then these don't collide; else, they do collide
	 for(int i = 0;i<h1.size();i++){
		 auto & p = h1[i];
		 auto & prev = (h1)[i-1];
		 Point<2,mpq_class> edge = p-prev;
		 Point<2,mpq_class> un_normalized_normal(-edge.y, edge.x);
		 projection_out1.clear();
		 projection_out2.clear();
		 //now project both polygons onto to this normal and see if they overlap, by finding the minimum and maximum distances
		 //Note that since we are NOT normalizing the normal vector, the projection is distorted along that vector
		 //(this still allows us to check overlaps, but means that the minimum distance found between the two shapes may be incorrect)
		 mpq_class left = std::numeric_limits<mpq_class>::max();
		 mpq_class right = std::numeric_limits<mpq_class>::min();
		 for (auto & p:h1){
			 mpq_class projection = un_normalized_normal.dot(p);
			 projection_out1.push({p,projection});
			 if (projection < left) {
				  left = projection;
			 } else if (projection > right) {
				  right = projection;
			 }
		 }
		 bool overlaps = false;
		 for (auto & p:h2){
			 mpq_class projection = un_normalized_normal.dot(p);
			 projection_out2.push({p,projection});
			 if (projection >= left && projection <= right ) {
				 overlaps=true;
				 break;
			 }
		 }
		 if(!overlaps){
			 //now place all the remaining (disabled) points on the axis as well
			 for(int i = 0;i<pointset1.size();i++){
				 if(!pointset1.pointEnabled(i)){
					 auto & p = pointset1[i];
					 mpq_class projection = un_normalized_normal.dot(p);
					 projection_out1.push({p,projection});
				 }
			 }

			 return true;
		 }
	 }

	 //now test the axis produced by the other polygon
	 for(int i = 0;i<h2.size();i++){
		 auto & p = (h1)[i];
		 auto & prev = (h1)[i-1];
		 Point<2,mpq_class> edge = p-prev;
		 Point<2,mpq_class> un_normalized_normal(-edge.y, edge.x);
		 projection_out1.clear();
		 projection_out2.clear();
		 mpq_class left = std::numeric_limits<mpq_class>::max();
		 mpq_class right = std::numeric_limits<mpq_class>::min();
		 for (auto & p:h1){
			 mpq_class projection = un_normalized_normal.dot(p);
			 projection_out1.push({p,projection});
			 if (projection < left) {
				  left = projection;
			 } else if (projection > right) {
				  right = projection;
			 }
		 }
		 bool overlaps = false;
		 for (auto & p:h2){
			 mpq_class projection = un_normalized_normal.dot(p);
			 projection_out2.push({p,projection});
			 if (projection >= left && projection <= right ) {
				 overlaps=true;
				 break;
			 }
		 }
		 if(!overlaps){
			 for(int i = 0;i<pointset1.size();i++){
				 if(!pointset1.pointEnabled(i)){
					 auto & p = pointset1[i];
					 mpq_class projection = un_normalized_normal.dot(p);
					 projection_out1.push({p,projection});
				 }
			 }
			 return true;
		 }
	 }
	 projection_out1.clear();
	 projection_out2.clear();
	 return false;


}

template<>
bool ConvexHullDetector<2, double>::findSeparatingAxis(ConvexPolygon<2, double> & hull1, ConvexPolygon<2, double> & hull2, PointSet<2,double> & pointset1, vec<std::pair<Point<2,double> ,double>>  &projection_out1, vec<std::pair<Point<2,double> ,double>>  &projection_out2){
	ConvexPolygon<2, double> & h1 = (hull1.size()<=hull2.size())?hull1:hull2;
		ConvexPolygon<2, double> & h2 = (hull1.size()<=hull2.size())?hull2:hull1;

		projection_out1.clear();
		projection_out2.clear();
		if(h1.size()==0 || h2.size()==0){
			return false;
		}else if (h1.size()==1 && h2.size()==1){
			 Point<2,double> un_normalized_normal = h2[0]-h1[0];
			 //now place all the remaining (disabled) points on the axis as well
			 for(int i = 0;i<pointset1.size();i++){
				 auto & p = pointset1[i];
				 double projection = un_normalized_normal.dot(p);
				 projection_out1.push({p,projection});
			 }
			 for(int i = 0;i<h2.size();i++){
				 auto & p = h2[i];
				 double projection = un_normalized_normal.dot(p);
				 projection_out2.push({p,projection});
			 }
			return true;
		}


		 std::vector<Point<2,double> > &  w = h1.getVertices();

		 //Separating Axis Theorem for collision detection between two convex polygons
		 //loop through each edge in _each_ polygon and project both polygons onto that edge's normal.
		 //If any of the projections are non-intersection, then these don't collide; else, they do collide
		 for(int i = 0;i<h1.size();i++){
			 auto & p = h1[i];
			 auto & prev = (h1)[i-1];
			 Point<2,double> edge = p-prev;
			 Point<2,double> un_normalized_normal(-edge.y, edge.x);
			 projection_out1.clear();
			 projection_out2.clear();
			 //now project both polygons onto to this normal and see if they overlap, by finding the minimum and maximum distances
			 //Note that since we are NOT normalizing the normal vector, the projection is distorted along that vector
			 //(this still allows us to check overlaps, but means that the minimum distance found between the two shapes may be incorrect)
			 double left = std::numeric_limits<double>::max();
			 double right = std::numeric_limits<double>::min();
			 for (auto & p:h1){
				 double projection = un_normalized_normal.dot(p);
				 projection_out1.push({p,projection});
				 if (projection < left) {
					  left = projection;
				 } else if (projection > right) {
					  right = projection;
				 }
			 }
			 bool overlaps = false;
			 for (auto & p:h2){
				 double projection = un_normalized_normal.dot(p);
				 projection_out2.push({p,projection});
				 if (projection >= left && projection <= right ) {
					 overlaps=true;
					 break;
				 }
			 }
			 if(!overlaps){
				 //now place all the remaining (disabled) points on the axis as well
				 for(int i = 0;i<pointset1.size();i++){
					 if(!pointset1.pointEnabled(i)){
						 auto & p = pointset1[i];
						 double projection = un_normalized_normal.dot(p);
						 projection_out1.push({p,projection});
					 }
				 }

				 return true;
			 }
		 }

		 //now test the axis produced by the other polygon
		 for(int i = 0;i<h2.size();i++){
			 auto & p = (h1)[i];
			 auto & prev = (h1)[i-1];
			 Point<2,double> edge = p-prev;
			 Point<2,double> un_normalized_normal(-edge.y, edge.x);
			 projection_out1.clear();
			 projection_out2.clear();
			 double left = std::numeric_limits<double>::max();
			 double right = std::numeric_limits<double>::min();
			 for (auto & p:h1){
				 double projection = un_normalized_normal.dot(p);
				 projection_out1.push({p,projection});
				 if (projection < left) {
					  left = projection;
				 } else if (projection > right) {
					  right = projection;
				 }
			 }
			 bool overlaps = false;
			 for (auto & p:h2){
				 double projection = un_normalized_normal.dot(p);
				 projection_out2.push({p,projection});
				 if (projection >= left && projection <= right ) {
					 overlaps=true;
					 break;
				 }
			 }
			 if(!overlaps){
				 for(int i = 0;i<pointset1.size();i++){
					 if(!pointset1.pointEnabled(i)){
						 auto & p = pointset1[i];
						 double projection = un_normalized_normal.dot(p);
						 projection_out1.push({p,projection});
					 }
				 }
				 return true;
			 }
		 }
		 projection_out1.clear();
		 projection_out2.clear();
		 return false;


}

template<>
void ConvexHullDetector<2,mpq_class>::buildPointOnHullOrDisabledReason(Var pointVar,const Point<2,mpq_class> & s, vec<Lit> & conflict){
	//If the point IS on the hull (or disabled), it is either because it is currently disabled...

	assert(outer->getPointVar(s.getID())==pointVar);
	Lit enabledLit = mkLit(pointVar,false);

	if(outer->value(enabledLit)==l_False){
		conflict.push(enabledLit);
		return;
	}

	//or because there is some line we can form that passes through this point, such that no enabled point is on the far side of the line, and the entire hull (except for this point)
	//is on the near side.
	//Note that as we are forcing the convex hull to contain all collinear points on the boundary,
	//this code is handles the case where enabled points are exactly on this line in the opposite way as  'buildPointNotContainedReason' does


	Line<2,mpq_class> testline;

	testline.a = s;
	ConvexPolygon<2,mpq_class> &hull = over_hull->getHull();

	//edge cases:
	if(hull.size()<3){
		//then we just report that _some_ disabled vertex must be enabled
		for(int i = 0;i<over.size();i++){
			Point<2,mpq_class> & p = over[i];
			if(!over.pointEnabled(i)){
				conflict.push(mkLit( outer->getPointVar(p.getID()),false));
			}
		}
		return;
	}

		Line<2,mpq_class> top_line;
		Line<2,mpq_class> bottom_line;
		top_line.a = s;
		bottom_line.a=s;
		top_line.b = hull.getVertices()[0];
		bottom_line.b = hull.getVertices().back();

		std::vector<Point<2,mpq_class> > min_set;
		std::vector<Point<2,mpq_class> > test_set;
		hasSet=false;

		for(int i = 0;i<over.size();i++){
			Point<2,mpq_class> & p = over[i];
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
			Point<2,mpq_class> & p = min_set[i];
			assert(!over.pointEnabled(p.getID()));
			Lit l = mkLit( outer->getPointVar(p.getID()),false);
			assert(outer->value(l)==l_False);
			conflict.push(l);

		}

}

template< >
void ConvexHullDetector<2,mpq_class>::buildPointNotOnHullOrDisabledReason(Var pointVar,const Point<2,mpq_class> & s, vec<Lit> & conflict){
	//If the point is not (on the hull or disabled), it is because it is both currently enabled...
	assert(outer->getPointVar(s.getID())==pointVar);
	Lit enabledLit = mkLit(pointVar,false);
	assert(outer->value(enabledLit)==l_True);
	conflict.push(~enabledLit);
	//AND, it is because there exist three other enabled points that form a triangle that contain this point:
	ConvexPolygon<2,mpq_class> &hull = under_hull->getHull();
#ifndef NDEBUG
	for(auto & p:hull)
		assert(p.getID()!=s.getID());
#endif
	ConvexPolygon<2,mpq_class> triangle;
	findContainingTriangle2d(hull,s,triangle);
	assert(triangle.contains(s));
	for(auto & p: triangle){
		int id = p.getID();
		Var v = outer->getPointVar(id);
		assert(outer->value(v)==l_True);
		conflict.push(mkLit(v,true));
	}
}
