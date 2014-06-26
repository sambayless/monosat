/*
 * graph_parser.h
 *
 *  Created on: 2013-06-08
 *      Author: sam
 */

#ifndef GEOMETRY_PARSER_H_
#define GEOMETRY_PARSER_H_


#include <stdio.h>

#include "utils/ParseUtils.h"
#include "graph/GraphParser.h"
#include "geometry/GeometryTypes.h"
#include "core/SolverTypes.h"
#include "geometry/GeometryTheory.h"
#include "core/Dimacs.h"
#include "core/Config.h"
namespace Minisat {


// GEOMETRY Parser:
template<class B, class Solver>
class GeometryParser:public Parser<B,Solver>{
	vec<GeometryTheorySolver<1,double>*> space_1D;
	vec<GeometryTheorySolver<2,double>*> space_2D;
	//vec<GeometryTheorySolver<3,double>*> space_3D;
	vec<char> tmp_str;
	struct ParsePoint{
		Var var;
		vec<double> position;
	};

	vec<vec<ParsePoint>> pointsets;

	struct ConvexHullArea{
		int pointsetID;
		double area;
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
	void readPoint(B& in, Solver& S ) {

		int pointsetID = parseInt(in);
		int d = parseInt(in);
		pointsets.growTo(pointsetID+1);
		pointsets[pointsetID].push();
		ParsePoint & point = pointsets[pointsetID].last();
		for(int i = 0;i<d;i++){
			double p = parseDouble(in,tmp_str);
			point.position.push(p);
		}
		int v = parseInt(in);
		point.var = v-1;
		while ( point.var >= S.nVars()) S.newVar();
	}

void readConvexHullArea(B& in, Solver& S) {

    //hull_area_lt pointsetID area var
	int pointset = parseInt(in); //ID of the hull
	double area =  parseDouble(in,tmp_str);
	int var = parseInt(in)-1;

	convex_hull_areas.push({pointset,area,var});
}
void readConvexHullPointContained(B& in, Solver& S) {

    //hull_point_contained pointsetID D p1 p2 ... pD var

    convex_hull_point_containments.push();
    ParsePoint & point = convex_hull_point_containments.last().point;
	int pointsetID = parseInt(in);
	int d = parseInt(in);
	for(int i = 0;i<d;i++){
		int p = parseDouble(in,tmp_str);
		point.position.push(p);
	}
	int v = parseInt(in)-1;
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
			 space_1D.growTo(i+1,nullptr);
			 if(!space_1D[i]){
				 space_1D[i] = new GeometryTheorySolver<1,double>(&S);
				 S.addTheory(space_1D[i]);
			 }
		 }else if (D==2){
			 space_2D.growTo(i+1,nullptr);
			 if(!space_2D[i]){
				 space_2D[i] = new GeometryTheorySolver<2,double>(&S);
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
				 Point<1,double> pnt(p.position);
				 space_1D[i]->newPoint(pnt,p.var);
			 }else if(D==2){
				 Point<2,double> pnt(p.position);
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
			 space_1D[c.pointsetID]->convexHullArea(c.area,c.v);
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
			 Point<1,double> pnt(c.point.position);
			 space_1D[c.pointsetID]->convexHullContains(pnt, c.point.var);
		 }else if(D==2){
			 Point<2,double> pnt(c.point.position);
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
			 space_1D[c.pointsetID]->convexHullContains(c.pointVar,c.pointOnHull);
		 }else if(D==2){
			 space_2D[c.pointsetID]->convexHullContains(c.pointVar,c.pointOnHull);
		 }/*else if(D==3){
			 space_3D[c.pointsetID]->convexHullContains(c.pointVar,c.pointOnHull);
		 }*/else{
			 assert(false);
		 }

	 }
 }

};

//=================================================================================================
};



#endif /* GRAPH_PARSER_H_ */
