/*
 * graph_parser.h
 *
 *  Created on: 2013-06-08
 *      Author: sam
 */

#ifndef GRAPH_PARSER_H_
#define GRAPH_PARSER_H_


#include <stdio.h>

#include "utils/ParseUtils.h"
#include "graph/GraphParser.h"
#include "core/SolverTypes.h"
#include "geometry/GeometryTheory.h"
#include "core/Dimacs.h"
#include "core/Config.h"
namespace Minisat {


// GEOMETRY Parser:
template<class B, class Solver>
class GeometryParser:public Parser<B,Solver>{

	GeometryTheorySolver * space_2D=nullptr;
	GeometryTheorySolver * space_3D=nullptr;


void readConvexHull(B& in, Solver& S, GeometryTheorySolver* G) {


    if(!eagerMatch(in,"hull")){
    	 printf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
    }


	int hullID = parseInt(in); //ID of the hull
	int d = parseInt(in);
	G->createConvexHull(hullID,d);

}

void readHullPoint(B& in, Solver& S,GeometryTheorySolver* G,vec<int> & point ) {
	if(opt_ignore_graph){
		skipLine(in);
		return;
	}
	 if(!eagerMatch(in,"hullpoint")){
    	printf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
    }
    ++in;
    point.clear();
        int hullID = parseInt(in);
        int d = parseInt(in);
        int hullVar = parseInt(in)-1;
        for(int i = 0;i<d;i++){
        	int p = parseInt(in);
        	point.push(p);
        }

        G->addHullPoint(hullID,mkLit(hullVar,false), point);

}
void readConvexHullPointContained(B& in, Solver& S,GeometryTheorySolver* G,vec<int> & point ) {
	if(opt_ignore_graph){
		skipLine(in);
		return;
	}

    ++in;
    point.clear();
        int pID = parseInt(in);
        int d = parseInt(in);
        int containedVar = parseInt(in)-1;
        for(int i = 0;i<d;i++){
        	int p = parseInt(in);
        	point.push(p);
        }

        G->addPointContainmentLit(pID,mkLit(containedVar,false), point);

}
public:

 bool parseLine(B& in, Solver& S){

		skipWhitespace(in);
		if (*in == EOF)
			return false;
		else if (*in == 'p'){
			return false;
		}else if (match(in,"point")){
			//add a point to a point set
		}else if (match(in,"convex_hull_area_lt")){

		}else if (match(in,"convex_hull_containment")){

		}else if (match(in,"heightmap_volume")){

		}else if (match(in, "euclidian_steiner_tree_weight")){

		}else if (match(in, "rectilinear_steiner_tree_weight")){

		}
		return false;
 }

 void implementConstraints(Solver & S){

 }

};

//=================================================================================================
};



#endif /* GRAPH_PARSER_H_ */
