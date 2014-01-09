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
#include "graph/GraphTheory.h"
#include "graph/TestGraph.h"
#include "core/Config.h"
namespace Minisat {

//=================================================================================================
// GRAPH Parser:

//Each graph file has the following (ascii) format:

//p graph\n
//c anything
//g digraph n e id
//e g u w var
//r g u w var
//...

//first line is just magic
//following lines are either g digraph nodes edges id
//or edge specifiers, which are of the form:
//e graphID fromNode toNode literal
//If lit is 0 or 1, this is a constant edge (0 for false, 1 for true)
//r g u w var is a reach querry: var is true if can u reach w in graph g, false otherwise

template<class B, class Solver>
static void readDiGraph(B& in, Solver& S, vec<GraphTheory*> & graphs) {
    int     g, n,e, ev;
    if(!eagerMatch(in,"digraph")){
    	 printf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
    }


        n = parseInt(in);//num nodes
        e = parseInt(in);//num edges (I'm ignoring this currently)
      //  ev = parseInt(in);//the variable of the first graph edge.
        g=parseInt(in);//id of the graph
        GraphTheory * graph = NULL;
        if(opt_graph)
        	graph= new GraphTheorySolver(&S);
        else
        	graph= new TestGraph(&S);
        graph->newNodes(n);
        graphs.growTo(g+1);
        graphs[g]=graph;
        S.addTheory(graph);
      //  return ev;
}

template<class B, class Solver>
static void readEdge(B& in, Solver& S, vec<GraphTheory*> & graphs) {

    if(*in != 'e'){
    	printf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
    }
    ++in;

        int graphID = parseInt(in);
        int from = parseInt(in);
        int to=parseInt(in);
        int edgeVar = parseInt(in)-1;
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
        graph->newEdge(from,to,edgeVar);

}

template<class B, class Solver>
static void readReach(B& in, Solver& S, vec<GraphTheory*> & graphs) {
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
            //if (eagerMatch(in, "p graph")){
                //vars    = parseInt(in);
                //clauses = parseInt(in);
                // SATRACE'06 hack
                // if (clauses > 4000000)
                //     S.eliminate(true);
           // }else{
           //     printf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
           // }
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
        }else if (*in == 'r'){
            readReach(in, S,graphs);
        }else if (*in == 'd'){
            readDistance(in, S,graphs);
        }else{
            cnt++;
            readClause(in, S, lits);
            S.addClause_(lits);
        }

    }

}

// Inserts problem into solver.
//
template<class Solver>
static void parse_GRAPH(gzFile input_stream, Solver& S, vec<std::pair<int,std::string> > * symbols=NULL) {
    StreamBuffer in(input_stream);
    parse_GRAPH_main(in, S,symbols); }

//=================================================================================================
}



#endif /* GRAPH_PARSER_H_ */
