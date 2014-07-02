/*
 * graph_parser.h
 *
 *  Created on: 2013-06-08
 *      Author: sam
 */

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



namespace Minisat {


// GEOMETRY Parser:
template<class B, class Solver, class T=double>
class GeometryParser:public Parser<B,Solver>{
	//vec<GeometryTheorySolver<1,mpq_class>*> space_1D;
	vec<GeometryTheorySolver<2,T>*> space_2D;
	//vec<GeometryTheorySolver<3,double>*> space_3D;
	vec<char> tmp_str;
	struct ParsePoint{
		Var var;
		vec<T> position;
	};

	vec<vec<ParsePoint>> pointsets;

	struct ConvexHullArea{
		int pointsetID;
		T area;
		Var v;
	};

	vec<ConvexHullArea> convex_hull_areas;

	struct ConvexHullPointContained{
		int pointsetID;
		ParsePoint point;

	};

	vec<ConvexHullPointContained> convex_hull_point_containments;

	struct ConvexHullPoint{
		int pointsetID;
		Var pointVar;
		Var pointOnHull;
	};
	vec<ConvexHullPoint> convex_hull_points;

	struct ConvexHullsIntersection{
		int pointsetID1;
		int pointsetID2;
		Var pointOnHull;
	};
	vec<ConvexHullsIntersection> convex_hulls_intersect;

	void parsePoint(B &in, int d, ParsePoint & point){
		//for now, points are given as 32-bit integers. This is less convenient than, say, floats, but makes parsing much easier.
		//in the future, we could also support inputting arbitrary precision integers... but for now, we aren't going to.

		//is there a better way to read in an arbitrary-sized vector of points?
		for(int i = 0;i<d;i++){

			int n = parseInt(in);
			//double p = parseDouble(in,tmp_str);
			point.position.push((T)n);
		}
	}


	void readPoint(B& in, Solver& S ) {

		int pointsetID = parseInt(in);

		int v = parseInt(in);
		int d = parseInt(in);
		pointsets.growTo(pointsetID+1);
		pointsets[pointsetID].push();
		//stringstream ss(in);
		ParsePoint & point = pointsets[pointsetID].last();
		parsePoint(in,d,point);

		point.var = v-1;
		while ( point.var >= S.nVars()) S.newVar();
	}

void readConvexHullArea(B& in, Solver& S) {

    //hull_area_lt pointsetID area var
	int pointset = parseInt(in); //ID of the hull
	int var = parseInt(in)-1;
	stringstream ss(in);
	T area;
	ss>>area;
	//T area =  parseDouble(in,tmp_str);
	convex_hull_areas.push({pointset,area,var});
}

void readConvexHullPointContained(B& in, Solver& S) {

    //hull_point_contained pointsetID D p1 p2 ... pD var

    convex_hull_point_containments.push();
    ParsePoint & point = convex_hull_point_containments.last().point;
	int pointsetID = parseInt(in);
	int v = parseInt(in)-1;
	int d = parseInt(in);
	parsePoint(in,d,point);

	//stringstream ss(in);
/*	for(int i = 0;i<d;i++){
		//skipWhitespace(in);
		long n = parseInt(in);
		//double p = parseDouble(in,tmp_str);
		point.position.push((T)n);
	}*/
/*	for(int i = 0;i<d;i++){

		ss>>p;
		point.position.push(p);
	}*/

	point.var = v;
	convex_hull_point_containments.last().pointsetID = pointsetID;
}
void readConvexHullPointOnHull(B& in, Solver& S) {

    //point_on_hull pointsetID pointVar var

	int pointsetID = parseInt(in);
	int pointVar = parseInt(in)-1;
	int var = parseInt(in)-1;
	convex_hull_points.push({pointsetID,pointVar, var});
}


void readConvexHullsIntersect(B& in, Solver& S) {

    //hulls_intersect pointsetID1 pointsetID2 var
	//v is true iff the two convex hulls intersect
	int pointsetID1 = parseInt(in);
	int pointsetID2 = parseInt(in);
	int v = parseInt(in)-1;
	convex_hulls_intersect.push({pointsetID1,pointsetID2,v});


}

public:

 bool parseLine(B& in, Solver& S){

		skipWhitespace(in);
		if (*in == EOF)
			return false;
		else if (match(in,"point")){
			//add a point to a point set
			readPoint(in,S);

		}else if (match(in,"convex_hull_area_gt")){
			readConvexHullArea(in,S);
		}else if (match(in,"convex_hull_containment")){
			readConvexHullPointContained(in,S);
		}else if (match(in,"point_on_convex_hull")){
			readConvexHullPointContained(in,S);
		}else if (match(in,"convex_hulls_intersect")){
			readConvexHullsIntersect(in,S);
		}else if (match(in,"heightmap_volume")){

		}else if (match(in, "euclidian_steiner_tree_weight")){

		}else if (match(in, "rectilinear_steiner_tree_weight")){

		}else{
			return false;
		}
		return true;
 }

 void implementConstraints(Solver & S){
	 //build point sets in their appropriate spaces.
	 //for now, we only support up to 3 dimensions
	 vec<int> pointsetDim;
	 for(int i = 0;i<pointsets.size();i++){
		 vec<ParsePoint> & pointset = pointsets[i];
		 if(pointset.size()==0)
			 continue;
		 ParsePoint & firstP = pointset[0];
		 int D = firstP.position.size();

		 if(D==1){
			/* space_1D.growTo(i+1,nullptr);
			 if(!space_1D[i]){
				 space_1D[i] = new GeometryTheorySolver<1,double>(&S);
				 S.addTheory(space_1D[i]);
			 }*/
		 }else if (D==2){
			 space_2D.growTo(i+1,nullptr);
			 if(!space_2D[i]){
				 space_2D[i] = new GeometryTheorySolver<2,T>(&S);
				 S.addTheory(space_2D[i]);
			 }
		 }/*else if (D==3){
			 space_3D.growTo(i+1,nullptr);
			 if(!space_3D[i]){
				 space_3D[i] = new GeometryTheorySolver<3,double>(&S);
				 S.addTheory(space_3D[i]);
			 }
		 }*/else{

			 fprintf(stderr,"Only points of dimension 1 and 2 currently supported (found point %d of dimension %d), aborting\n", i,D);
			 exit(3);
		 }

		 pointsetDim.growTo(i+1,-1);
		 pointsetDim[i]=D;
		 for(ParsePoint & p:pointset){
			 if(p.position.size()!=D){
				 fprintf(stderr,"All points in a pointset must have the same dimensionality\n");
				 exit(3);
			 }
			 if(D==1){
				/* Point<1,double> pnt(p.position);
				 space_1D[i]->newPoint(pnt,p.var);*/
			 }else if(D==2){
				 Point<2,T> pnt(p.position);
				 space_2D[i]->newPoint(pnt,p.var);
			 }/*else if(D==3){
				 Point<3,double> pnt(p.position);
				 space_3D[i]->newPoint(pnt,p.var);
			 }*/else{
				 assert(false);
			 }
		 }
	 }

	 for (auto & c: convex_hull_areas){
		 if(c.pointsetID>=pointsetDim.size() || c.pointsetID<0 || pointsetDim[c.pointsetID]<0){
			 fprintf(stderr,"Bad pointsetID %d\n", c.pointsetID);
		 }
		 int D = pointsetDim[c.pointsetID];
		 if(D==1){
			// space_1D[c.pointsetID]->convexHullArea(c.area,c.v);
		 }else if(D==2){
			 space_2D[c.pointsetID]->convexHullArea(c.area,c.v);
		 }/*else if(D==3){
			 space_3D[c.pointsetID]->convexHullArea(c.area,c.v);
		 }*/else{
			 assert(false);
		 }

	 }
	 for (auto & c: convex_hull_point_containments){
		 if(c.pointsetID>=pointsetDim.size() || c.pointsetID<0 || pointsetDim[c.pointsetID]<0){
			 fprintf(stderr,"Bad pointsetID %d\n", c.pointsetID);
		 }
		 int D = pointsetDim[c.pointsetID];

		 if(D==1){
			 Point<1,T> pnt(c.point.position);
			// space_1D[c.pointsetID]->convexHullContains(pnt, c.point.var);
		 }else if(D==2){
			 Point<2,T> pnt(c.point.position);
			 space_2D[c.pointsetID]->convexHullContains(pnt, c.point.var);
		 }/*else if(D==3){
			 Point<3,double> pnt(c.point.position);
			 space_3D[c.pointsetID]->convexHullContains(pnt, c.point.var);
		 }*/else{
			 assert(false);
		 }

	 }

	 for (auto & c: convex_hull_points){
		 if(c.pointsetID>=pointsetDim.size() || c.pointsetID<0 || pointsetDim[c.pointsetID]<0){
			 fprintf(stderr,"Bad pointsetID %d\n", c.pointsetID);
		 }
		 int D = pointsetDim[c.pointsetID];
		 if(D==1){
			 //space_1D[c.pointsetID]->convexHullContains(c.pointVar,c.pointOnHull);
		 }else if(D==2){
			 space_2D[c.pointsetID]->convexHullContains(c.pointVar,c.pointOnHull);
		 }/*else if(D==3){
			 space_3D[c.pointsetID]->convexHullContains(c.pointVar,c.pointOnHull);
		 }*/else{
			 assert(false);
		 }

	 }

	 for (auto & c:convex_hulls_intersect){
		 if(c.pointsetID1>=pointsetDim.size() || c.pointsetID1<0 || pointsetDim[c.pointsetID1]<0){
			 fprintf(stderr,"Bad pointsetID %d\n", c.pointsetID1);
		 }
		 if(c.pointsetID2>=pointsetDim.size() || c.pointsetID2<0 || pointsetDim[c.pointsetID2]<0){
			 fprintf(stderr,"Bad pointsetID %d\n", c.pointsetID2);
		 }
		 int D = pointsetDim[c.pointsetID1];
		 if( pointsetDim[c.pointsetID2] != D){
			 fprintf(stderr,"Cannot intersect convex hulls in different dimensions (%d has dimension %d, while %d has dimension %d)\n",c.pointsetID1,D, c.pointsetID2, pointsetDim[c.pointsetID2]);
		 }
		 if(D==1){

		 }else if(D==2){

			 space_2D[c.pointsetID]->convexHullsIntersect(c.pointsetID1,c.pointsetID2,c.var);
		 }else{
			 assert(false);
		 }


	 }

 }

};

//=================================================================================================
};



#endif /* GRAPH_PARSER_H_ */
