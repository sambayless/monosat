#include "ConvexHullCollisionDetector.h"
#include "LineSegment.h"
#include "GeometryTheory.h"

template<>
void ConvexHullCollisionDetector<2, mpq_class>::buildCollisionReason(vec<Lit> & conflict){

	//If the two polygons intersect, there are two possible cases (these are not mutually exclusive)
	//1) 1 polygon has at least one point contained in the other polygon. Learn that either that point must be disabled, or at least
	//		one point from its containing triangle must be disabled (as in the point containment theory).
	//2) There exists an edge (or, more generally, a line segment between any two vertices, as opposed to just an edge) from each poygon,
	//		such that the two edges intersect. Learn that one of the end points must be disabled.
	ConvexPolygon<2, mpq_class> & h1 = hull1_under.getHull();
	ConvexPolygon<2, mpq_class> & h2 = hull2_under.getHull();


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
		if(h2.contains(h1[i])){
			static ConvexPolygon<2,mpq_class> triangle;
			findContainingTriangle2d(h2,h1[i],triangle);
			assert(triangle.contains(h1[i]));
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
		if(h1.contains(h2[i])){
			static ConvexPolygon<2,mpq_class> triangle;
			findContainingTriangle2d(h1,h2[i],triangle);
			assert(triangle.contains(h2[i]));
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

}

template<>
void ConvexHullCollisionDetector<2, mpq_class>::buildNotCollisionReason(vec<Lit> & conflict){

	ConvexPolygon<2, mpq_class> & h1 = hull1_over.getHull();
	ConvexPolygon<2, mpq_class> & h2 = hull2_over.getHull();


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
	}else if (h2.size()==0){
		for(int i = 0;i<over2.size();i++){
			if(!over2.pointEnabled(i)){
				Lit l = mkLit( outer->getPointVar(over2[i].getID()));
				assert(outer->value(l)==l_False);
				conflict.push(l);
			}
		}
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
	}

	//the reason that the two convex hulls do NOT collide is that there exists a separating axis between them.
	//Find that axis, and then find all the positions of the disabled points of both polygons along that axis.
	//let h1 be left of h2.
	//Then either some disabled point in h1 that is to the right of the rightmost point in h1 must be enabled,
	//or some disabled point in h2 that is to the left of the leftmost point in h2 must be enabled.

	//Note: It may be possible to improve on this analysis!
	vec<std::pair<Point<2,mpq_class> ,mpq_class>>  projection1;
	vec<std::pair<Point<2,mpq_class> ,mpq_class>>  projection2;
	findSeparatingAxis(h1,h2,over1,over2, projection1,projection2);

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
	bool h1_is_left=rightmost1<leftmost2;

	for(auto & p:projection1){
		int pID = p.first.getID();
		int pointset = outer->getPointset(pID);
		int pointsetIndex = outer->getPointsetIndex(pID);
		assert(pointset==over1.getID());
		if(!over1.pointEnabled(pointsetIndex)){
			if(h1_is_left ? (p.second>rightmost1):((p.second<leftmost1))){
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
			if(h1_is_left ? (p.second<leftmost2):((p.second>rightmost2))){
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
void ConvexHullCollisionDetector<2, double>::buildCollisionReason(vec<Lit> & conflict){

}

template<>
void ConvexHullCollisionDetector<2, double>::buildNotCollisionReason(vec<Lit> & conflict){

}
template<>
bool ConvexHullCollisionDetector<2, mpq_class>::findSeparatingAxis(ConvexPolygon<2, mpq_class> & hull1, ConvexPolygon<2, mpq_class> & hull2, PointSet<2,mpq_class> & pointset1, PointSet<2,mpq_class> & pointset2, vec<std::pair<Point<2,mpq_class> ,mpq_class>>  &projection_out1_t,vec<std::pair<Point<2,mpq_class> ,mpq_class>>  &projection_out2_t){
	ConvexPolygon<2, mpq_class> & h1 = (hull1.size()<=hull2.size())?hull1:hull2;
	ConvexPolygon<2, mpq_class> & h2 = (hull1.size()<=hull2.size())?hull2:hull1;
	vec<std::pair<Point<2,mpq_class> ,mpq_class>>  &projection_out1 = (hull1.size()<=hull2.size())?projection_out1_t:projection_out2_t;
	vec<std::pair<Point<2,mpq_class> ,mpq_class>>  &projection_out2 = (hull1.size()<=hull2.size())?projection_out2_t:projection_out1_t;

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
		 for(int i = 0;i<pointset2.size();i++){
				 auto & p = pointset2[i];
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
			 for(int i = 0;i<pointset2.size();i++){
				 if(!pointset2.pointEnabled(i)){
					 auto & p = pointset2[i];
					 mpq_class projection = un_normalized_normal.dot(p);
					 projection_out2.push({p,projection});
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
			 for(int i = 0;i<pointset2.size();i++){
				 if(!pointset2.pointEnabled(i)){
					 auto & p = pointset2[i];
					 mpq_class projection = un_normalized_normal.dot(p);
					 projection_out2.push({p,projection});
				 }
			 }
			 return true;
		 }
	 }
	 projection_out1.clear();
	 projection_out2.clear();
	 return false;


}
