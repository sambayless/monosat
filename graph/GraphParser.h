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
#include "pb/PbTheory.h"
#include <set>
#include <string>
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
//r g u w var is a reach query: var is true if can u reach w in graph g, false otherwise

template<class B, class Solver>
static void readDiGraph(B& in, Solver& S, vec<GraphTheory*> & graphs) {
	if(opt_ignore_theories){
		skipLine(in);
		return;
	}

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
        	graph= new GraphTheorySolver(&S,g);
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
	if(opt_ignore_theories){
		skipLine(in);
		return;
	}
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
static void readWeightedEdge(B& in, Solver& S, vec<GraphTheory*> & graphs) {
	if(opt_ignore_theories){
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
static void readConnect(B& in, Solver& S, vec<GraphTheory*> & graphs) {
	if(opt_ignore_theories){
		skipLine(in);
		return;
	}
	//u g u w var is an undirected reachability query: var is true if can u reach w in graph g, false otherwise
    if(*in != 'u'){
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
        while (reachVar >= S.nVars()) S.newVar();
        graph->connects(from,to,reachVar);


}


template<class B, class Solver>
static void readReach(B& in, Solver& S, vec<GraphTheory*> & graphs) {
	if(opt_ignore_theories){
		skipLine(in);
		return;
	}
	//r g u w var is a reach query: var is true if can u reach w in graph g, false otherwise
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
        while (reachVar >= S.nVars()) S.newVar();
        graph->reaches(from,to,reachVar);


}

template<class B, class Solver>
static void readDistance(B& in, Solver& S, vec<GraphTheory*> & graphs) {
	if(opt_ignore_theories){
		skipLine(in);
		return;
	}
	//d g u w var dist is a reach query: var is true if can u reach w in graph g, false otherwise
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
        while (reachVar >= S.nVars()) S.newVar();
        graph->reaches(from,to,reachVar,steps);


}



template<class B, class Solver>
static void readMinSpanningTreeConstraint(B& in, Solver& S, vec<GraphTheory*> & graphs) {
	if(opt_ignore_theories){
		skipLine(in);
		return;
	}
	//m g maxweight var is a minimum spanning tree weight constraint: var is true if the mst of the graph is <= maxweight
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
        while (reachVar >= S.nVars()) S.newVar();
        graph->minimumSpanningTree(reachVar,maxweight);


}
template<class B, class Solver>
static void readMinSpanningTreeEdgeConstraint(B& in, Solver& S, vec<GraphTheory*> & graphs) {
	if(opt_ignore_theories){
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
        while (reachVar >= S.nVars()) S.newVar();
        graph->edgeInMinimumSpanningTree(edgeVar,reachVar);


}


template<class B, class Solver>
static void readMaxFlowConstraint(B& in, Solver& S, vec<GraphTheory*> & graphs) {
	if(opt_ignore_theories){
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
        while (reachVar >= S.nVars()) S.newVar();
        graph->maxFlow(s,t,flow,reachVar);


}




template<class B, class Solver>
static void readMinConnectedComponentsConstraint(B& in, Solver& S, vec<GraphTheory*> & graphs) {
	if(opt_ignore_theories){
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
        while (reachVar >= S.nVars()) S.newVar();
        graph->minConnectedComponents(min_components,reachVar);
}

struct SteinerStruct{
	int id;
	vec<std::pair<int, Var> > terminals;
	vec<std::pair<int,Var> > weight_constraints;
	SteinerStruct(int id):id(id){

	}
};

/*template<class B, class Solver>
static void readSteinerTreeDeclaration(B& in, Solver& S, vec<GraphTheory*> & graphs, vec<vec<SteinerStruct*>> steiners) {
	if(opt_ignore_theories){
		skipLine(in);
		return;
	}
	//steiner_tree graphID steinerID
	//Creates a steiner tree with the given id. Steiner constraints can then refer to that steinerID.
	if(!eagerMatch(in,"steiner_tree")){
    	printf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
    }


        int graphID = parseInt(in);
        int steinerID = parseInt(in);

        if(graphID <0 || graphID>=graphs.size() || !graphs[graphID]){
         	printf("PARSE ERROR! Undeclared graph identifier %d\n",graphID), exit(3);
         }

        steiners.growTo(graphs.size());
        steiners[graphID].growTo(steinerID+1);
        if(steiners[graphID][steinerID]!=0){
        	printf("PARSE ERROR! Multiple declarations of steiner tree %d for graph %d\n",steinerID, graphID), exit(3);
        }
        steiners[graphID][steinerID] = new SteinerStruct{steinerID};
}*/
template<class B, class Solver>
static void readSteinerTreeTerminal(B& in, Solver& S, vec<GraphTheory*> & graphs, vec<vec<SteinerStruct*>>& steiners) {
	if(opt_ignore_theories){
		skipLine(in);
		return;
	}
	//steiner_terminal graphID steinerID node var
	//Declares that the specified node in graphID is a terminal in the steiner tree (conditioned on var being assigned to true)



	int graphID = parseInt(in);
	int steinerID = parseInt(in);
	int node = parseInt(in);
	int var = parseInt(in)-1;
	if(graphID <0 || graphID>=graphs.size() || !graphs[graphID]){
		printf("PARSE ERROR! Undeclared graph identifier %d \n",graphID), exit(3);
	}
    steiners.growTo(graphs.size());
    steiners[graphID].growTo(steinerID+1);
	if(!steiners[graphID][steinerID]){
		steiners[graphID][steinerID] = new SteinerStruct(steinerID);
	}
	steiners[graphID][steinerID]->terminals.push({node,var});

}
template<class B, class Solver>
static void readSteinerTreeConstraint(B& in, Solver& S, vec<GraphTheory*> & graphs, vec<vec<SteinerStruct*>>& steiners) {
	if(opt_ignore_theories){
		skipLine(in);
		return;
	}
	//steiner_minweight graphID steinerID maxweight var is a minimum spanning tree weight constraint: var is true iff the minimum steiner tree of the graph is <= maxweight

	int graphID = parseInt(in);
	int steinerID = parseInt(in);
	int maxweight = parseInt(in);
	int var = parseInt(in)-1;
	if(graphID <0 || graphID>=graphs.size() || !graphs[graphID]){
		printf("PARSE ERROR! Undeclared graph identifier %d \n",graphID), exit(3);
	}
    steiners.growTo(graphs.size());
    steiners[graphID].growTo(steinerID+1);
	if(!steiners[graphID][steinerID]){
		steiners[graphID][steinerID] = new SteinerStruct{steinerID};
	}
	steiners[graphID][steinerID]->weight_constraints.push({maxweight,var});
}

template<class B, class Solver>
static void readPB(B & in,vec<Lit> & lits, vec<int> & weights, Solver & S, PbTheory * pb){
	if(opt_ignore_theories){
		skipLine(in);
		return;
	}
	//pb constraints are in this form:
	//first a size integer, >= 1
	//then, a list of literals
	//then, either a 0, or another size integer followed by an optional same-length list of weights
	//next, an operator: '<','<=','=','>','>=', '!=' (the last one is non-standard...)
	//followed by a comparison integer
	//next, EITHER a 0, or a 1, or a 2.
	//If a 0, then this is an unconditionally true constraint.
	//If a 1 or 2, then there must be one last literal.
	//If 1, then this is a one sided constraint, which is enforced if the final literal is true, and otherwise has no effect.
	//If 2, then this is a two sided constraint: if the final literal is false, then the condtion is enforced to be false.

	//o <size> lit1 lit2 ... [0 | <size>weight1 weight2 ...] 'op' total 2 headLiteral
	 if(!eagerMatch(in,"o")){
		 printf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
	 }
	 lits.clear();
	 weights.clear();
	 int size =  parseInt(in);
	 if(size<=0){
		 printf("PARSE ERROR! Empty PB clause\n"), exit(3);
	 }

	 for(int i = 0;i<size;i++){
		 int parsed_lit = parseInt(in);
		 if (parsed_lit == 0) break;
		 int v = abs(parsed_lit)-1;
		 while (v >= S.nVars()) S.newVar();
		 lits.push( (parsed_lit > 0) ? mkLit(v) : ~mkLit(v) );
	 }

	 int wsize = parseInt(in);
	 if(wsize != 0 && wsize !=size){
		 printf("PARSE ERROR! Number of weights must either be the same as the size of the clause, or 0.\n"), exit(3);
	 }
	 for(int i = 0;i<wsize;i++){
		 int parsed_weight = parseInt(in);
		 weights.push(parsed_weight);
	 }
	 if(wsize==0){
		 for(int i = 0;i<size;i++)
			 weights.push(1);
	 }
	 skipWhitespace(in);
	 PbTheory::PbType op;
	 //read the operator:
	 if(*in == '<'){
		 ++in;
		 if(*in=='='){
			 ++in;
			 op = PbTheory::PbType::LE;
		 }else{
			 op = PbTheory::PbType::LT;
		 }
	 }else if (*in=='>'){
		 ++in;
		 if(*in=='='){
			 ++in;
			 op = PbTheory::PbType::GE;
		 }else{
			 op = PbTheory::PbType::GT;
		 }
	 }else if (*in == '!'){
		 ++in;
		 if(*in !='='){
			 printf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
		 }
		 ++in;
		 op = PbTheory::PbType::NE;
	 }else if (*in == '='){
		 ++in;
		 op = PbTheory::PbType::EQ;
	 }
	 int comparison = parseInt(in);
	 int type = parseInt(in);
	 Lit head = lit_Undef;
	 bool oneSided =false;
	 if(type==0){
		 //done

	 }else {
		 if(type==1){
			 oneSided=true;
		 }
		 int parsed_lit = parseInt(in);

		 int v = abs(parsed_lit)-1;
		 while (v >= S.nVars()) S.newVar();
		 head =  (parsed_lit > 0) ? mkLit(v) : ~mkLit(v) ;
	 }
	 assert(lits.size()==weights.size());
	 pb->addConstraint(lits,weights,comparison,head,op,oneSided?PbTheory::ConstraintSide::Upper: PbTheory::ConstraintSide::Both);
}

template<class B, class Solver>
static void parse_GRAPH_main(B& in, Solver& S, vec<std::pair<int,std::string> > * symbols=NULL ) {
	vec<GraphTheory*> graphs;
	vec<vec<SteinerStruct*>> steiners;
	PbTheory * pbtheory=nullptr;
	vec<Lit> lits;
	vec<int> weights;
    int vars    = 0;
    int clauses = 0;
    int cnt     = 0;
    //std::set<std::string> used_symbols;
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
     /*   		if(symbols && used_symbols.count(symbol)){
        			printf("PARSE ERROR! Duplicated symbol: %c\n", *symbol.c_str()), exit(3);
        		}
        		used_symbols.insert(symbol);*/

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
        }else if (*in == 'u'){
            readConnect(in, S,graphs);
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
        }else if (*in == 's'){
        	if(!eagerMatch(in,"steiner_")){
        		printf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
        	}
        	if(*in=='t'){
        		if(!eagerMatch(in,"terminal")){
        			printf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
				}
        		readSteinerTreeTerminal(in,S,graphs, steiners);
        	}else if (*in=='m'){
          		if(!eagerMatch(in,"minweight")){
          			printf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
    				}
        		readSteinerTreeConstraint(in,S,graphs, steiners);
        	}
        }else if (*in == 'o'){
        	//Pseudoboolean constraint: o is for opb... can't use p, unfortunately...
        	if(!pbtheory){
        		pbtheory = new PbTheory(&S);
        	    S.addTheory(pbtheory);
        	}
        	readPB(in,lits,weights,S,pbtheory);
        }else{
            cnt++;
            readClause(in, S, lits);
            S.addClause_(lits);
        }

    }
    //clear any unmapped symbols (have to do this _before_ implementing constraints, which may introduce new variables)
    if(symbols){
		int i,j=0;
		for(i = 0;i<symbols->size();i++ ){
			std::pair<int,std::string> p = (*symbols)[i];
			Var v = p.first;
			if(v<=S.nVars()){
				//keep this symbol
				(*symbols)[j++] = (*symbols)[i];
			}
		}
		symbols->shrink(i-j);
    }

    for(int gid = 0;gid<steiners.size();gid++){
    	for(auto & steiner:steiners[gid]){
    		if(steiner){
				graphs[gid]->addSteinerTree(steiner->terminals,steiner->id);
				for(auto & weight : steiner->weight_constraints){
					graphs[gid]->addSteinerWeightConstraint(steiner->id, weight.first,weight.second);
				}
				delete(steiner);
    		}
    	}
    }

    for(int i = 0;i<graphs.size();i++){
    	if(graphs[i])
    	   graphs[i]->implementConstraints();
    }
    if(pbtheory)
    	pbtheory->implementConstraints();
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
