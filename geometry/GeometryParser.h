/****************************************************************************************[Solver.h]
 The MIT License (MIT)

 Copyright (c) 2014, Sam Bayless

 Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
 associated documentation files (the "Software"), to deal in the Software without restriction,
 including without limitation the rights to use, copy, modify, merge, publish, distribute,
 sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all copies or
 substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
 NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
 OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 **************************************************************************************************/

#ifndef GEOMETRY_PARSER_H_
#define GEOMETRY_PARSER_H_

#include <stdio.h>
#include <gmpxx.h>
#include "utils/ParseUtils.h"
#include "graph/GraphParser.h"
#include "geometry/GeometryTypes.h"
#include "core/SolverTypes.h"
#include "geometry/GeometryTheory.h"
#include "core/Dimacs.h"
#include "core/Config.h"
#include "mtl/Vec.h"
#include <vector>

namespace Monosat {

// GEOMETRY Parser:
template<class B, class Solver, class T = double>
class GeometryParser: public Parser<B, Solver> {
	using Parser<B, Solver>::mapVar;
	//vec<GeometryTheorySolver<1,mpq_class>*> space_1D;
	GeometryTheorySolver<2, T> *space_2D = nullptr;
	//vec<GeometryTheorySolver<3,double>*> space_3D;
	vec<char> tmp_str;
	struct ParsePoint {
		Var var;
		vec<T> position;
	};

	vec<vec<ParsePoint>> pointsets;

	struct ConvexHullArea {
		int pointsetID;
		T area;
		Var v;
	};

	vec<ConvexHullArea> convex_hull_areas;

	/*
	 struct ConvexHullPointContained{
	 int pointsetID;
	 ParsePoint point;

	 };

	 vec<ConvexHullPointContained> convex_hull_point_containments;
	 */

	/*	struct ConvexHullLineIntersection{
	 int pointsetID;
	 Var var;
	 ParsePoint p;
	 ParsePoint q;
	 };

	 vec<ConvexHullLineIntersection> convex_hull_line_intersections;*/

	struct ConvexHullPolygonIntersection {
		int pointsetID;
		bool inclusive;
		Var var;
		vec<ParsePoint> points;
	};

	vec<ConvexHullPolygonIntersection> convex_hull_polygon_intersections;

	struct ConvexHullPoint {
		int pointsetID;
		Var pointVar;
		Var pointOnHull;
	};
	vec<ConvexHullPoint> convex_hull_points;

	struct ConvexHullsIntersection {
		int pointsetID1;
		int pointsetID2;
		Var var;
		bool inclusive;
	};
	vec<ConvexHullsIntersection> convex_hulls_intersect;

	void parsePoint(B &in, int d, ParsePoint & point) {
		//for now, points are given as 32-bit integers. This is less convenient than, say, floats, but makes parsing much easier.
		//in the future, we could also support inputting arbitrary precision integers... but for now, we aren't going to.
		
		//is there a better way to read in an arbitrary-sized vector of points?
		for (int i = 0; i < d; i++) {
			
			int n = parseInt(in);
			//double p = parseDouble(in,tmp_str);
			point.position.push((T) n);
		}
	}
	
	void readPoint(B& in, Solver& S) {
		
		int pointsetID = parseInt(in);
		
		Var v = parseInt(in)-1;
		v = mapVar(S,v);
		int d = parseInt(in);
		pointsets.growTo(pointsetID + 1);
		pointsets[pointsetID].push();
		//stringstream ss(in);
		ParsePoint & point = pointsets[pointsetID].last();
		parsePoint(in, d, point);
		
		point.var = v;

	}
	
	void readConvexHullArea(B& in, Solver& S) {
		
		//hull_area_lt pointsetID area var
		int pointset = parseInt(in); //ID of the hull
		Var var = parseInt(in) - 1;
		var = mapVar(S,var);
		stringstream ss(in);
		T area;
		ss >> area;
		//T area =  parseDouble(in,tmp_str);
		convex_hull_areas.push( { pointset, area, var });
	}
	void readConvexHullIntersectsPolygon(B& in, Solver& S, bool inclusive) {
		//convex_hull_intersects_polygon pointsetID var NumPoints D p1 p2 ... pD q1 q2 ... qD ...
		
		convex_hull_polygon_intersections.push();
		convex_hull_polygon_intersections.last().inclusive = inclusive;
		
		int pointsetID = parseInt(in);
		Var v = parseInt(in) - 1;
		v = mapVar(S,v);
		convex_hull_polygon_intersections.last().pointsetID = pointsetID;
		convex_hull_polygon_intersections.last().var = v;
		int numPoints = parseInt(in);
		
		int d = parseInt(in);
		for (int i = 0; i < numPoints; i++) {
			convex_hull_polygon_intersections.last().points.push();
			ParsePoint & p = convex_hull_polygon_intersections.last().points.last();
			parsePoint(in, d, p);
		}
	}
	/*void readConvexHullIntersectsLine(B& in, Solver& S){
	 //convex_hull_intersects_line pointsetID var D p1 p2 ... qD q1 q2 ... qD

	 convex_hull_line_intersections.push();
	 ParsePoint & p = convex_hull_line_intersections.last().p;
	 ParsePoint & q = convex_hull_line_intersections.last().q;
	 int pointsetID = parseInt(in);
	 int v = parseInt(in)-1;
	 convex_hull_line_intersections.last().pointsetID = pointsetID;
	 convex_hull_line_intersections.last().var = v;
	 int d = parseInt(in);
	 parsePoint(in,d,p);
	 parsePoint(in,d,q);

	 //stringstream ss(in);
	 for(int i = 0;i<d;i++){
	 //skipWhitespace(in);
	 long n = parseInt(in);
	 //double p = parseDouble(in,tmp_str);
	 point.position.push((T)n);
	 }
	 for(int i = 0;i<d;i++){

	 ss>>p;
	 point.position.push(p);
	 }



	 }*/
	/*
	 void readConvexHullPointContained(B& in, Solver& S) {

	 //hull_point_contained pointsetID var D p1 p2 ... pD

	 convex_hull_point_containments.push();
	 ParsePoint & point = convex_hull_point_containments.last().point;
	 int pointsetID = parseInt(in);
	 int v = parseInt(in)-1;
	 int d = parseInt(in);
	 parsePoint(in,d,point);

	 //stringstream ss(in);
	 for(int i = 0;i<d;i++){
	 //skipWhitespace(in);
	 long n = parseInt(in);
	 //double p = parseDouble(in,tmp_str);
	 point.position.push((T)n);
	 }
	 for(int i = 0;i<d;i++){

	 ss>>p;
	 point.position.push(p);
	 }

	 point.var = v;
	 convex_hull_point_containments.last().pointsetID = pointsetID;
	 }*/
	void readConvexHullPointOnHull(B& in, Solver& S) {
		
		//point_on_hull pointsetID pointVar var
		
		int pointsetID = parseInt(in);
		int pointVar = parseInt(in) - 1;
		Var var = parseInt(in) - 1;
		var = mapVar(S,var);
		convex_hull_points.push( { pointsetID, pointVar, var });
	}
	
	void readConvexHullsIntersect(B& in, Solver& S, bool inclusive) {
		
		//hulls_intersect pointsetID1 pointsetID2 var
		//v is true iff the two convex hulls intersect
		int pointsetID1 = parseInt(in);
		int pointsetID2 = parseInt(in);
		Var v = parseInt(in) - 1;
		v = mapVar(S,v);
		convex_hulls_intersect.push( { pointsetID1, pointsetID2, v, inclusive });
		
	}
	
public:
	
	GeometryParser():Parser<B, Solver>("Geometry"){

	}


	bool parseLine(B& in, Solver& S) {
		
		skipWhitespace(in);
		if (*in == EOF) {
			return false;
		} else if (match(in, "convex_hull_area_gt")) {
			readConvexHullArea(in, S);
		}/*else if (match(in,"convex_hull_containment")){
		 readConvexHullPointContained(in,S);
		 }else if (match(in,"convex_hull_intersects_line")){
		 readConvexHullIntersectsLine(in,S);
		 }*/
		else if (match(in, "convex_hull_collides_polygon")) {
			readConvexHullIntersectsPolygon(in, S, true);
		} else if (match(in, "convex_hull_intersects_polygon")) {
			readConvexHullIntersectsPolygon(in, S, false);
		} else if (match(in, "point_on_convex_hull")) {
			readConvexHullPointOnHull(in, S);
		} else if (match(in, "convex_hulls_intersect")) {
			readConvexHullsIntersect(in, S, false);
		} else if (match(in, "convex_hulls_collide")) {
			readConvexHullsIntersect(in, S, true);
		} else if (match(in, "heightmap_volume")) {
			
		} else if (match(in, "euclidian_steiner_tree_weight")) {
			
		} else if (match(in, "rectilinear_steiner_tree_weight")) {
			
		} else if (match(in, "point")) {
			//add a point to a point set
			readPoint(in, S);
			
		} else {
			return false;
		}
		return true;
	}
	
	void implementConstraints(Solver & S) {
		//build point sets in their appropriate spaces.
		//for now, we only support up to 3 dimensions
		vec<int> pointsetDim;
		for (int i = 0; i < pointsets.size(); i++) {
			vec<ParsePoint> & pointset = pointsets[i];
			if (pointset.size() == 0)
				continue;
			ParsePoint & firstP = pointset[0];
			int D = firstP.position.size();
			
			if (D == 1) {
				/* space_1D.growTo(i+1,nullptr);
				 if(!space_1D[i]){
				 space_1D[i] = new GeometryTheorySolver<1,double>(&S);
				 S.addTheory(space_1D[i]);
				 }*/
			} else if (D == 2) {
				
				if (!space_2D) {
					space_2D = new GeometryTheorySolver<2, T>(&S);
					S.addTheory(space_2D);
				}
			}/*else if (D==3){
			 space_3D.growTo(i+1,nullptr);
			 if(!space_3D[i]){
			 space_3D[i] = new GeometryTheorySolver<3,double>(&S);
			 S.addTheory(space_3D[i]);
			 }
			 }*/else {
				
				parse_errorf("Only points of dimension 1 and 2 currently supported (found point %d of dimension %d), aborting\n",
						i, D);

			}
			
			pointsetDim.growTo(i + 1, -1);
			pointsetDim[i] = D;
			for (ParsePoint & p : pointset) {
				if (p.position.size() != D) {
					parse_errorf("All points in a pointset must have the same dimensionality\n");

				}
				if (D == 1) {
					/* Point<1,double> pnt(p.position);
					 space_1D[i]->newPoint(pnt,p.var);*/
				} else if (D == 2) {
					Point<2, T> pnt(p.position);
					space_2D->newPoint(i, pnt, p.var);
				}/*else if(D==3){
				 Point<3,double> pnt(p.position);
				 space_3D[i]->newPoint(pnt,p.var);
				 }*/else {
					assert(false);
				}
			}
		}
		
		for (auto & c : convex_hull_areas) {
			if (c.pointsetID >= pointsetDim.size() || c.pointsetID < 0 || pointsetDim[c.pointsetID] < 0) {
				parse_errorf("Bad pointsetID %d\n", c.pointsetID);
			}
			if (c.pointsetID < 0 || c.pointsetID >= pointsetDim.size() || pointsetDim[c.pointsetID] < 0) {
				parse_errorf( "Pointset %d is undefined, aborting!", c.pointsetID);

			}
			int D = pointsetDim[c.pointsetID];
			
			if (D == 1) {
				// space_1D[c.pointsetID]->convexHullArea(c.area,c.v);
			} else if (D == 2) {
				space_2D->convexHullArea(c.pointsetID, c.area, c.v);
			}/*else if(D==3){
			 space_3D[c.pointsetID]->convexHullArea(c.area,c.v);
			 }*/else {
				assert(false);
			}
			
		}
		for (auto & c : convex_hull_polygon_intersections) {
			if (c.pointsetID >= pointsetDim.size() || c.pointsetID < 0 || pointsetDim[c.pointsetID] < 0) {
				parse_errorf( "Bad pointsetID %d\n", c.pointsetID);
			}
			if (c.pointsetID < 0 || c.pointsetID >= pointsetDim.size() || pointsetDim[c.pointsetID] < 0) {
				parse_errorf( "Pointset %d is undefined, aborting!", c.pointsetID);

			}
			int D = pointsetDim[c.pointsetID];
			
			if (D == 1) {
				
				// space_1D[c.pointsetID]->convexHullContains(pnt, c.point.var);
			} else if (D == 2) {
				
				std::vector<Point<2, T> > p2;
				for (ParsePoint &p : c.points) {
					p2.push_back(p.position);
				}
				space_2D->convexHullIntersectsPolygon(c.pointsetID, p2, c.var, c.inclusive);
			} else if (D == 3) {
				
			} else {
				assert(false);
			}
		}
		/*	 for (auto & c: convex_hull_point_containments){
		 if(c.pointsetID>=pointsetDim.size() || c.pointsetID<0 || pointsetDim[c.pointsetID]<0){
		 parse_errorf("Bad pointsetID %d\n", c.pointsetID);
		 }
		 int D = pointsetDim[c.pointsetID];

		 if(D==1){
		 Point<1,T> pnt(c.point.position);
		 // space_1D[c.pointsetID]->convexHullContains(pnt, c.point.var);
		 }else if(D==2){
		 Point<2,T> pnt(c.point.position);
		 space_2D->convexHullContains(c.pointsetID,pnt, c.point.var);
		 }else if(D==3){
		 Point<3,double> pnt(c.point.position);
		 space_3D[c.pointsetID]->convexHullContains(pnt, c.point.var);
		 }else{
		 assert(false);
		 }
		 }*/

		/*	 for (auto & c:convex_hull_line_intersections){
		 if(c.pointsetID>=pointsetDim.size() || c.pointsetID<0 || pointsetDim[c.pointsetID]<0){
		 parse_errorf("Bad pointsetID %d\n", c.pointsetID);
		 }
		 int D = pointsetDim[c.pointsetID];

		 if(D==1){
		 // Point<1,T> pnt(c.point.position);
		 // space_1D[c.pointsetID]->convexHullContains(pnt, c.point.var);
		 }else if(D==2){
		 Point<2,T> p(c.p.position);
		 Point<2,T> q(c.q.position);
		 space_2D->convexHullIntersectsLine(c.pointsetID,p,q,c.var);
		 }else if(D==3){
		 Point<3,double> pnt(c.point.position);
		 space_3D[c.pointsetID]->convexHullContains(pnt, c.point.var);
		 }else{
		 assert(false);
		 }
		 }*/

		for (auto & c : convex_hull_points) {
			if (c.pointsetID >= pointsetDim.size() || c.pointsetID < 0 || pointsetDim[c.pointsetID] < 0) {
				parse_errorf( "Bad pointsetID %d\n", c.pointsetID);
			}
			if (c.pointsetID < 0 || c.pointsetID >= pointsetDim.size() || pointsetDim[c.pointsetID] < 0) {
				parse_errorf( "Pointset %d is undefined, aborting!", c.pointsetID);

			}
			parse_errorf("point_on_convex_hull not currently supported, aborting.\n");

			/* int D = pointsetDim[c.pointsetID];

			 if(D==1){
			 //space_1D[c.pointsetID]->convexHullContains(c.pointVar,c.pointOnHull);
			 }else if(D==2){
			 if(!S.hasTheory(c.pointVar)){
			 parse_errorf("Variable %d is not a point variable, aborting!",c.pointVar);

			 }
			 Var theoryVar = S.getTheoryVar(c.pointVar);
			 int pointID = space_2D->getPointID(theoryVar);
			 int pointset = space_2D->getPointset(pointID);
			 assert(pointset == c.pointsetID);
			 int pointIndex = space_2D->getPointsetIndex(pointID);

			 assert(pointIndex>=0);
			 space_2D->pointOnHull(c.pointsetID,pointIndex,c.pointOnHull);
			 }else if(D==3){
			 space_3D[c.pointsetID]->convexHullContains(c.pointVar,c.pointOnHull);
			 }else{
			 assert(false);
			 }*/

		}
		
		for (auto & c : convex_hulls_intersect) {
			if (c.pointsetID1 < 0 || c.pointsetID1 >= pointsetDim.size() || c.pointsetID1 < 0
					|| pointsetDim[c.pointsetID1] < 0) {
				parse_errorf( "Bad pointsetID %d\n", c.pointsetID1);

			}
			if (c.pointsetID2 < 0 || c.pointsetID2 >= pointsetDim.size() || c.pointsetID2 < 0
					|| pointsetDim[c.pointsetID2] < 0) {
				parse_errorf( "Bad pointsetID %d\n", c.pointsetID2);

			}
			int D = pointsetDim[c.pointsetID1];
			
			if (pointsetDim[c.pointsetID2] != D) {
				parse_errorf(
						"Cannot intersect convex hulls in different dimensions (%d has dimension %d, while %d has dimension %d)\n",
						c.pointsetID1, D, c.pointsetID2, pointsetDim[c.pointsetID2]);
			}
			if (D == 1) {
				
			} else if (D == 2) {
				
				space_2D->convexHullsIntersect(c.pointsetID1, c.pointsetID2, c.var, c.inclusive);
			} else {
				assert(false);
			}
			
		}
		
	}
	
};

//=================================================================================================
}
;

#endif /* GRAPH_PARSER_H_ */
