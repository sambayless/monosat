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
#include "core/SolverTypes.h"
#include "geometry/GeometryTheory.h"

#include "core/Config.h"
namespace Minisat {

//=================================================================================================
// GEOMETRY Parser:




template<class B, class Solver>
static void readConvexHull(B& in, Solver& S, GeometryTheorySolver* G) {
	if(opt_ignore_graph){
		skipLine(in);
		return;
	}

    if(!eagerMatch(in,"hull")){
    	 printf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
    }


	int hullID = parseInt(in); //ID of the hull
	int d = parseInt(in);
	G->createConvexHull(hullID,d);

}

template<class B, class Solver>
static void readHullPoint(B& in, Solver& S,GeometryTheorySolver* G,vec<int> & point ) {
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
template<class B, class Solver>
static void readPointContained(B& in, Solver& S,GeometryTheorySolver* G,vec<int> & point ) {
	if(opt_ignore_graph){
		skipLine(in);
		return;
	}
	 if(!eagerMatch(in,"hullpoint")){
    	printf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
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

template<class B, class Solver>
static void readWeightedEdge(B& in, Solver& S, vec<GraphTheory*> & graphs) {
	if(opt_ignore_graph){
		skipLine(in);
		return;
	}
    if(*in != 'w'){
    	printf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
    }
    ++in;

        int graphID = parseInt(in);
        int from = parseInt(in);
        int to=parseInt(in);
        int edgeVar = parseInt(in)-1;
        int weight = parseInt(in);
        /*if(edgeVar==-1){
        	edgeVar=edge_var++-1;
        }*/
        if(graphID <0 || graphID>=graphs.size() || !graphs[graphID]){
        	printf("PARSE ERROR! Undeclared graph identifier %d for edge %d\n",graphID, edgeVar), exit(3);
        }
        if(edgeVar<0){
        	printf("PARSE ERROR! Edge variables must be >=0, was %d\n", edgeVar), exit(3);
        }
        GraphTheory * graph = graphs[graphID];
        while (edgeVar >= S.nVars()) S.newVar();
        graph->newEdge(from,to,edgeVar,weight);

}

template<class B, class Solver>
static void readReach(B& in, Solver& S, vec<GraphTheory*> & graphs) {
	if(opt_ignore_graph){
		skipLine(in);
		return;
	}
	//r g u w var is a reach querry: var is true if can u reach w in graph g, false otherwise
    if(*in != 'r'){
    	printf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
    }
    ++in;

        int graphID = parseInt(in);
        int from = parseInt(in);
       // int steps = parseInt(in);
         int to=parseInt(in);
        int reachVar = parseInt(in)-1;
        if(graphID <0 || graphID>=graphs.size() || !graphs[graphID]){
        	printf("PARSE ERROR! Undeclared graph identifier %d for edge %d\n",graphID, reachVar), exit(3);
        }
        if(reachVar<0){
        	printf("PARSE ERROR! Edge variables must be >=0, was %d\n", reachVar), exit(3);
        }
        GraphTheory * graph = graphs[graphID];
        while (reachVar+graph->nNodes() >= S.nVars()) S.newVar();
        graph->reaches(from,to,reachVar);


}

template<class B, class Solver>
static void readDistance(B& in, Solver& S, vec<GraphTheory*> & graphs) {
	if(opt_ignore_graph){
		skipLine(in);
		return;
	}
	//d g u w var dist is a reach querry: var is true if can u reach w in graph g, false otherwise
    if(*in != 'd'){
    	printf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
    }
    ++in;

        int graphID = parseInt(in);
        int from = parseInt(in);
         int to=parseInt(in);
        int reachVar = parseInt(in)-1;
        int steps = parseInt(in);
        if(graphID <0 || graphID>=graphs.size() || !graphs[graphID]){
        	printf("PARSE ERROR! Undeclared graph identifier %d for edge %d\n",graphID, reachVar), exit(3);
        }
        if(reachVar<0){
        	printf("PARSE ERROR! Edge variables must be >=0, was %d\n", reachVar), exit(3);
        }
        GraphTheory * graph = graphs[graphID];
        while (reachVar+graph->nNodes() >= S.nVars()) S.newVar();
        graph->reaches(from,to,reachVar,steps);


}



template<class B, class Solver>
static void readMinSpanningTreeConstraint(B& in, Solver& S, vec<GraphTheory*> & graphs) {
	if(opt_ignore_graph){
		skipLine(in);
		return;
	}
	//m g maxweight var is a minimum spanning tree weight constraint: var is true if the mst of the tree is <= maxweight
    if(*in != 'm'){
    	printf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
    }
    ++in;

        int graphID = parseInt(in);
        int maxweight = parseInt(in);
        int reachVar = parseInt(in)-1;
        if(graphID <0 || graphID>=graphs.size() || !graphs[graphID]){
        	printf("PARSE ERROR! Undeclared graph identifier %d for edge %d\n",graphID, reachVar), exit(3);
        }
        if(reachVar<0){
        	printf("PARSE ERROR! Edge variables must be >=0, was %d\n", reachVar), exit(3);
        }
        GraphTheory * graph = graphs[graphID];
        while (reachVar+graph->nNodes() >= S.nVars()) S.newVar();
        graph->minimumSpanningTree(reachVar,maxweight);


}
template<class B, class Solver>
static void readMinSpanningTreeEdgeConstraint(B& in, Solver& S, vec<GraphTheory*> & graphs) {
	if(opt_ignore_graph){
		skipLine(in);
		return;
	}
	//x g u v var is a minimum spanning tree EDGE constraint: var is true if the UNIQUE mst of the tree includes the UNDIRECTED edge from u to v. (directed graphs will be treated as undirected here)
	//Note that if there are multiple edges from u to v, this constraint is true if any of them are in the mst.
    if(*in != 'x'){
    	printf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
    }
    ++in;

        int graphID = parseInt(in);
        int edgeVar = parseInt(in)-1;
        //int from = parseInt(in);
        //int to = parseInt(in);
        int reachVar = parseInt(in)-1;
        if(graphID <0 || graphID>=graphs.size() || !graphs[graphID]){
        	printf("PARSE ERROR! Undeclared graph identifier %d for edge %d\n",graphID, reachVar), exit(3);
        }
        if(reachVar<0){
        	printf("PARSE ERROR! Edge variables must be >=0, was %d\n", reachVar), exit(3);
        }
        GraphTheory * graph = graphs[graphID];
        while (reachVar+graph->nNodes() >= S.nVars()) S.newVar();
        graph->edgeInMinimumSpanningTree(edgeVar,reachVar);


}


template<class B, class Solver>
static void readMaxFlowConstraint(B& in, Solver& S, vec<GraphTheory*> & graphs) {
	if(opt_ignore_graph){
		skipLine(in);
		return;
	}
	//f g s t flow var is a max flow constraint, true if the max flow is >= var
    if(*in != 'f'){
    	printf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
    }
    ++in;

        int graphID = parseInt(in);
        int  s= parseInt(in);
        int t = parseInt(in);
        int flow = parseInt(in);
        int reachVar = parseInt(in)-1;
        if(graphID <0 || graphID>=graphs.size() || !graphs[graphID]){
        	printf("PARSE ERROR! Undeclared graph identifier %d for edge %d\n",graphID, reachVar), exit(3);
        }
        if(reachVar<0){
        	printf("PARSE ERROR! Edge variables must be >=0, was %d\n", reachVar), exit(3);
        }
        GraphTheory * graph = graphs[graphID];
        while (reachVar+graph->nNodes() >= S.nVars()) S.newVar();
        graph->maxFlow(s,t,flow,reachVar);


}




template<class B, class Solver>
static void readMinConnectedComponentsConstraint(B& in, Solver& S, vec<GraphTheory*> & graphs) {
	if(opt_ignore_graph){
		skipLine(in);
		return;
	}
	//n g min_components var is a minimum connected component count constraint, true if the max flow is >= var
    if(*in != 'n'){
    	printf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
    }
    ++in;

        int graphID = parseInt(in);
        int min_components = parseInt(in);
        int reachVar = parseInt(in)-1;
        if(graphID <0 || graphID>=graphs.size() || !graphs[graphID]){
        	printf("PARSE ERROR! Undeclared graph identifier %d for edge %d\n",graphID, reachVar), exit(3);
        }
        if(reachVar<0){
        	printf("PARSE ERROR! Edge variables must be >=0, was %d\n", reachVar), exit(3);
        }
        GraphTheory * graph = graphs[graphID];
        while (reachVar+graph->nNodes() >= S.nVars()) S.newVar();
        graph->minConnectedComponents(min_components,reachVar);
}

template<class B, class Solver>
static void parse_GRAPH_main(B& in, Solver& S, vec<std::pair<int,std::string> > * symbols=NULL ) {
	vec<GraphTheory*> graphs;
	vec<Lit> lits;
    int vars    = 0;
    int clauses = 0;
    int cnt     = 0;
    std::string symbol;
    for (;;){
        skipWhitespace(in);
        if (*in == EOF) break;
        else if (*in == 'p'){
        	skipLine(in);

        } else if (*in == 'c'){
        	if(symbols && eagerMatch(in,"c var ")){
        		//this is a variable symbol map
        		skipWhitespace(in);
        		int v = parseInt(in);
        		if(v<=0){
        			printf("PARSE ERROR! Variables must be positive: %c\n", *in), exit(3);
        		}

        		v--; //subtract one to get the variable id

        		symbol.clear();
        		skipWhitespace(in);
        		while(*in != '\n' && ! isWhitespace(*in)){
        			symbol.push_back(*in);
        			++in;
        		}
        		if(symbol.size()==0){
        			printf("PARSE ERROR! Empty symbol: %c\n", *in), exit(3);
        		}
        		symbols->push();
        		symbols->last().first=v;
        		symbols->last().second=symbol;
        		skipLine(in);
        	}else{
        		//just a comment
        		skipLine(in);
        	}
        } else if (*in == 'p'){
            skipLine(in);
        }else if (*in == 'g'){
        	++in;
        	skipWhitespace(in);
        	if(*in=='d'){
        		//for now, only digraphs are supported
        		 readDiGraph(in, S,graphs);
        	}else{
        		printf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
        	}
        }else if (*in == 'e'){
            cnt++;
            readEdge(in, S,graphs);
        }else if (*in == 'w'){
            cnt++;
            readWeightedEdge(in, S,graphs);
        }else if (*in == 'r'){
            readReach(in, S,graphs);
        }else if (*in == 'd'){
            readDistance(in, S,graphs);
        }else if (*in == 'm'){
            readMinSpanningTreeConstraint(in, S,graphs);
        }else if (*in == 'x'){
            readMinSpanningTreeEdgeConstraint(in, S,graphs);
        }else if (*in == 'f'){
            readMaxFlowConstraint(in, S,graphs);
        }else if (*in == 'n'){
            readMinConnectedComponentsConstraint(in, S,graphs);
        }else{
            cnt++;
            readClause(in, S, lits);
            S.addClause_(lits);
        }

    }
    for(int i = 0;i<graphs.size();i++){
    	if(graphs[i])
    	   graphs[i]->implementConstraints();
    }

}

// Inserts problem into solver.
//
template<class Solver>
static void parse_GRAPH(gzFile input_stream, Solver& S, vec<std::pair<int,std::string> > * symbols=NULL) {
    StreamBuffer in(input_stream);
    parse_GRAPH_main(in, S,symbols);

}

//=================================================================================================
}



#endif /* GRAPH_PARSER_H_ */
