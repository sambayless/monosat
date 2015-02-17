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

#ifndef GRAPH_PARSER_H_
#define GRAPH_PARSER_H_

#include <stdio.h>

#include "utils/ParseUtils.h"
#include "core/SolverTypes.h"
#include "graph/GraphTheory.h"
#include "bv/BVTheorySolver.h"
#include "core/Config.h"
#include "pb/PbTheory.h"
#include "core/Dimacs.h"
#include <gmpxx.h>
#include <set>
#include <string>
#include <sstream>
namespace Monosat {

struct SteinerStruct {
	int id;
	vec<std::pair<int, Var> > terminals;
	vec<std::pair<int, Var> > weight_constraints;
	SteinerStruct(int id) :
			id(id) {
		
	}
};

//=================================================================================================
// GRAPH Parser:
template<class B, class Solver>
class GraphParser: public Parser<B, Solver> {
	bool precise;
	BVTheorySolver<long>*& bvTheory;
	vec<GraphTheorySolver<long>*> graphs;
	vec<GraphTheorySolver<double>*> graphs_float;
	vec<GraphTheorySolver<mpq_class>*> graphs_rational;

	enum class GraphType {
		INTEGER, FLOAT, RATIONAL
	};

	vec<vec<SteinerStruct*>> steiners;
	PbTheory * pbtheory = nullptr;

	vec<int> weights;
	vec<Lit> lits;
	int count = 0;
	vec<char> tmp;
	
	struct BVEdge{
		int graphID;
		int from;
		int to;
		int edgeVar;
		int bvID;
	};
	vec<BVEdge> bvedges;

	void readDiGraph(B& in, GraphType graph_type, Solver& S) {
		if (opt_ignore_theories) {
			skipLine(in);
			return;
		}
		
		if (!precise) {
			if (graph_type == GraphType::RATIONAL) {
				graph_type = GraphType::FLOAT;
			}
		}
		
		int g, n, e, ev;
		
		n = parseInt(in); //num nodes
		e = parseInt(in); //num edges (I'm ignoring this currently)
		//  ev = parseInt(in);//the variable of the first graph edge.
		g = parseInt(in);  //id of the graph
		graphs.growTo(g + 1);
		graphs_float.growTo(g + 1);
		graphs_rational.growTo(g + 1);
		if (graph_type == GraphType::INTEGER) {
			GraphTheorySolver<long> *graph = new GraphTheorySolver<long>(&S, g);
			graph->newNodes(n);
			graphs[g] = graph;
			S.addTheory(graph);
		} else if (graph_type == GraphType::FLOAT) {
			GraphTheorySolver<double> *graph = new GraphTheorySolver<double>(&S, g);
			graph->newNodes(n);
			graphs_float[g] = graph;
			S.addTheory(graph);
		} else if (graph_type == GraphType::RATIONAL) {
			GraphTheorySolver<mpq_class> *graph = new GraphTheorySolver<mpq_class>(&S, g);
			graph->newNodes(n);
			graphs_rational[g] = graph;
			S.addTheory(graph);
		}
		//  return ev;
	}

	void readEdgeBV(B& in, Solver& S) {
		if (opt_ignore_theories) {
			skipLine(in);
			return;
		}

		++in;

		int graphID = parseInt(in);
		int from = parseInt(in);
		int to = parseInt(in);
		int edgeVar = parseInt(in) - 1;

		if (graphID < 0 || graphID >= graphs.size()) {
			printf("PARSE ERROR! Undeclared graph identifier %d for edge %d\n", graphID, edgeVar), exit(1);
		}
		if (edgeVar < 0) {
			printf("PARSE ERROR! Edge variables must be >=0, was %d\n", edgeVar), exit(1);
		}
		while (edgeVar >= S.nVars())
			S.newVar();

		skipWhitespace(in);
		int bvID = parseInt(in);
/*		static vec<Var> bv;
		bv.clear();

		while(*in != "\n"){
			skipWhitespace(in);
			int bvVar = parseInt(in)-1;
			bv.push(bvVar);
		}*/

		if (graphs[graphID]) {
			bvedges.push({graphID,from, to, edgeVar, bvID});
			//graphs[graphID]->newEdgeBV(from, to, edgeVar, bvID);
		}/* else if (graphs_float[graphID]) {
			graphs_float[graphID]->newEdge(from, to, edgeVar,bvID);
		} else if (graphs_rational[graphID]) {
			graphs_rational[graphID]->newEdge(from, to, edgeVar,bvID);
		} */else {
			printf("PARSE ERROR! Undeclared graph identifier %d for edge %d\n", graphID, edgeVar), exit(1);
			exit(1);
		}
	}
	
	void readEdge(B& in, Solver& S) {
		if (opt_ignore_theories) {
			skipLine(in);
			return;
		}
		
		++in;

		int graphID = parseInt(in);
		int from = parseInt(in);
		int to = parseInt(in);
		int edgeVar = parseInt(in) - 1;
		
		if (graphID < 0 || graphID >= graphs.size()) {
			printf("PARSE ERROR! Undeclared graph identifier %d for edge %d\n", graphID, edgeVar), exit(1);
		}
		if (edgeVar < 0) {
			printf("PARSE ERROR! Edge variables must be >=0, was %d\n", edgeVar), exit(1);
		}
		while (edgeVar >= S.nVars())
			S.newVar();
		
		skipWhitespace(in);
		if(*in=='\n' || *in==0){
			//this is an unweighted edge
			
			if (graphs[graphID]) {
				graphs[graphID]->newEdge(from, to, edgeVar);
			} else if (graphs_float[graphID]) {
				graphs_float[graphID]->newEdge(from, to, edgeVar);
			} else if (graphs_rational[graphID]) {
				graphs_rational[graphID]->newEdge(from, to, edgeVar);
			} else {
				printf("PARSE ERROR! Undeclared graph identifier %d for edge %d\n", graphID, edgeVar), exit(1);
				exit(1);
			}
			
		}else{

			if (graphs[graphID]) {
				int weight = parseInt(in);
				graphs[graphID]->newEdge(from, to, edgeVar, weight);
			} else if (graphs_float[graphID]) {
				//float can be either a plain integer, or a rational in the form '123/456', or a floating point in decimal format
				double weight = parseDouble(in, tmp);
				skipWhitespace(in);
				if (*in == '/') {
					++in;
					double denom = parseDouble(in, tmp);
					weight /= denom;
				}
				graphs_float[graphID]->newEdge(from, to, edgeVar, weight);
			} else if (graphs_rational[graphID]) {
				//rational can be either a plain integer, or a rational in the form '123/456', or a floating point value
				std::stringstream ss;
				skipWhitespace(in);

				while (*in != '\n') {
					ss << (*in);
					++in;
				}

				//first, try to interpret this string as a double:
				try {
					double value = std::stod(ss.str());
					mpq_class weight(value);
					weight.canonicalize();
					graphs_rational[graphID]->newEdge(from, to, edgeVar, weight);
				} catch (std::exception& e) {
					//if that fails, attempt to read it in directly as an fraction:
					mpq_class weight(ss.str());
					weight.canonicalize();
					graphs_rational[graphID]->newEdge(from, to, edgeVar, weight);
				}

			} else {
				printf("PARSE ERROR! Undeclared graph identifier %d for edge %d\n", graphID, edgeVar), exit(1);
				exit(1);
			}
		}
	}

	void readReach(B& in, Solver& S) {
		if (opt_ignore_theories) {
			skipLine(in);
			return;
		}
		//reach grachID u w var is a reach query: var is true if can u reach w in graph g, false otherwise
		
		++in;
		
		int graphID = parseInt(in);
		int from = parseInt(in);
		// int steps = parseInt(in);
		int to = parseInt(in);
		int reachVar = parseInt(in) - 1;
		if (graphID < 0 || graphID >= graphs.size()) {
			printf("PARSE ERROR! Undeclared graph identifier %d for edge %d\n", graphID, reachVar), exit(1);
		}
		if (reachVar < 0) {
			printf("PARSE ERROR! Edge variables must be >=0, was %d\n", reachVar), exit(1);
		}
		
		while (reachVar >= S.nVars())
			S.newVar();
		
		if (graphs[graphID]) {
			graphs[graphID]->reaches(from, to, reachVar);
		} else if (graphs_float[graphID]) {
			graphs_float[graphID]->reaches(from, to, reachVar);
		} else if (graphs_rational[graphID]) {
			graphs_rational[graphID]->reaches(from, to, reachVar);
		} else {
			printf("PARSE ERROR! Undeclared graph identifier %d\n", graphID), exit(1);
			exit(1);
		}
	}
	
	void readDistance(B& in, Solver& S, bool leq = false) {
		if (opt_ignore_theories) {
			skipLine(in);
			return;
		}
		//distance_lt grachID u w var dist is a reach query: var is true if can u reach w in graph g, false otherwise
		
		++in;
		
		int graphID = parseInt(in);
		int from = parseInt(in);
		int to = parseInt(in);
		int reachVar = parseInt(in) - 1;
		int steps = parseInt(in);
		if (graphID < 0 || graphID >= graphs.size()) {
			printf("PARSE ERROR! Undeclared graph identifier %d for edge %d\n", graphID, reachVar), exit(1);
		}
		if (reachVar < 0) {
			printf("PARSE ERROR! Edge variables must be >=0, was %d\n", reachVar), exit(1);
		}
		
		if (graphs[graphID]) {
			graphs[graphID]->reaches(from, to, reachVar, steps);
		} else if (graphs_float[graphID]) {
			graphs_float[graphID]->reaches(from, to, reachVar, steps);
		} else if (graphs_rational[graphID]) {
			graphs_rational[graphID]->reaches(from, to, reachVar, steps);
		} else {
			printf("PARSE ERROR! Undeclared graph identifier %d\n", graphID), exit(1);
			exit(1);
		}
	}
	void readShortestPath(B& in, Solver& S, bool leq = false) {
		if (opt_ignore_theories) {
			skipLine(in);
			return;
		}
		//distance_lt grachID u w var dist is a reach query: var is true if can u reach w in graph g, false otherwise
		
		++in;
		static vec<char> tmp;
		int graphID = parseInt(in);
		int from = parseInt(in);
		int to = parseInt(in);
		int reachVar = parseInt(in) - 1;
		
		if (graphID < 0 || graphID >= graphs.size()) {
			printf("PARSE ERROR! Undeclared graph identifier %d for edge %d\n", graphID, reachVar), exit(1);
		}
		if (reachVar < 0) {
			printf("PARSE ERROR! Edge variables must be >=0, was %d\n", reachVar), exit(1);
		}
		
		if (graphs[graphID]) {
			int weight = parseInt(in);
			graphs[graphID]->distance(from, to, reachVar, weight);
		} else if (graphs_float[graphID]) {
			double weight = parseDouble(in, tmp);
			skipWhitespace(in);
			if (*in == '/') {
				++in;
				double denom = parseDouble(in, tmp);
				weight /= denom;
			}
			graphs_float[graphID]->distance(from, to, reachVar, weight);
		} else if (graphs_rational[graphID]) {
			std::stringstream ss;
			skipWhitespace(in);
			//rational can be either a plain integer, or a rational in the form '123/456', or a floating point value
			while (*in != '\n') {
				ss << (*in);
				++in;
			}
			//first, try to interpret this string as a double:
			try {
				double value = std::stod(ss.str());
				mpq_class weight(value);
				weight.canonicalize();
				graphs_rational[graphID]->distance(from, to, reachVar, weight);
			} catch (std::exception& e) {
				//if that fails, attempt to read it in directly as an fraction:
				mpq_class weight(ss.str());
				weight.canonicalize();
				graphs_rational[graphID]->distance(from, to, reachVar, weight);
			}
			
		} else {
			printf("PARSE ERROR! Undeclared graph identifier %d\n", graphID), exit(1);
			exit(1);
		}
	}

	void readAcyclic(B& in, Solver& S, bool directed = false) {
			if (opt_ignore_theories) {
				skipLine(in);
				return;
			}
			//distance_lt grachID u w var dist is a reach query: var is true if can u reach w in graph g, false otherwise

			++in;
			static vec<char> tmp;
			int graphID = parseInt(in);
			int reachVar = parseInt(in) - 1;

			if (graphID < 0 || graphID >= graphs.size()) {
				printf("PARSE ERROR! Undeclared graph identifier %d for edge %d\n", graphID, reachVar), exit(1);
			}
			if (reachVar < 0) {
				printf("PARSE ERROR! Variables must be >=0, was %d\n", reachVar), exit(1);
			}
			if (graphs[graphID]) {
				graphs[graphID]->acyclic(reachVar,directed);
			} else if (graphs_float[graphID]) {
				graphs_float[graphID]->acyclic(reachVar,directed);
			} else if (graphs_rational[graphID]) {
				graphs_rational[graphID]->acyclic(reachVar,directed);
			} else {
				printf("PARSE ERROR! Undeclared graph identifier %d\n", graphID), exit(1);
				exit(1);
			}
		}

	void readMinSpanningTreeConstraint(B& in, Solver& S, bool inclusive) {
		if (opt_ignore_theories) {
			skipLine(in);
			return;
		}

		++in;
		
		int graphID = parseInt(in);

		int reachVar = parseInt(in) - 1;


		if (graphID < 0 || graphID >= graphs.size()) {
			printf("PARSE ERROR! Undeclared graph identifier %d for edge %d\n", graphID, reachVar), exit(1);
		}
		if (reachVar < 0) {
			printf("PARSE ERROR! Edge variables must be >=0, was %d\n", reachVar), exit(1);
		}
		
		while (reachVar >= S.nVars())
			S.newVar();
		if (graphs[graphID]) {
			long maxweight = parseInt(in);
			graphs[graphID]->minimumSpanningTree(reachVar, maxweight,inclusive);
		} else if (graphs_float[graphID]) {
			//float can be either a plain integer, or a rational (interpreted at floating point precision) in the form '123/456', or a floating point in decimal format
			double maxweight = parseDouble(in, tmp);
			skipWhitespace(in);
			if (*in == '/') {
				++in;
				double denom = parseDouble(in, tmp);
				maxweight /= denom;
			}
			graphs_float[graphID]->minimumSpanningTree(reachVar, maxweight,inclusive);
		} else if (graphs_rational[graphID]) {
			//rational can be either a plain integer, or a rational in the form '123/456', or a floating point value
			std::stringstream ss;
			skipWhitespace(in);

			while (*in != '\n') {
				ss << (*in);
				++in;
			}

			//first, try to interpret this string as a double:
			try {
				double value = std::stod(ss.str());
				mpq_class maxweight(value);
				maxweight.canonicalize();
				graphs_rational[graphID]->minimumSpanningTree(reachVar, maxweight,inclusive);
			} catch (std::exception& e) {
				//if that fails, attempt to read it in directly as an fraction:
				mpq_class maxweight(ss.str());
				maxweight.canonicalize();
				graphs_rational[graphID]->minimumSpanningTree(reachVar, maxweight,inclusive);
			}

		} else {
			printf("PARSE ERROR! Undeclared graph identifier %d\n", graphID), exit(1);
			exit(1);
		}

	}
	void readMinSpanningTreeEdgeConstraint(B& in, Solver& S) {
		if (opt_ignore_theories) {
			skipLine(in);
			return;
		}
		//mst_edge grachID u v var is a minimum spanning tree EDGE constraint: var is true if the UNIQUE mst of the tree includes the UNDIRECTED edge from u to v. (directed graphs will be treated as undirected here)
		//Note that if there are multiple edges from u to v, this constraint is true if any of them are in the mst.
		
		++in;
		
		int graphID = parseInt(in);
		int edgeVar = parseInt(in) - 1;
		//int from = parseInt(in);
		//int to = parseInt(in);
		int reachVar = parseInt(in) - 1;
		if (graphID < 0 || graphID >= graphs.size()) {
			printf("PARSE ERROR! Undeclared graph identifier %d for edge %d\n", graphID, reachVar), exit(1);
		}
		if (reachVar < 0) {
			printf("PARSE ERROR! Edge variables must be >=0, was %d\n", reachVar), exit(1);
		}
		
		while (reachVar >= S.nVars())
			S.newVar();
		
		if (graphs[graphID]) {
			graphs[graphID]->edgeInMinimumSpanningTree(edgeVar, reachVar);
		} else if (graphs_float[graphID]) {
			graphs_float[graphID]->edgeInMinimumSpanningTree(edgeVar, reachVar);
		} else if (graphs_rational[graphID]) {
			graphs_rational[graphID]->edgeInMinimumSpanningTree(edgeVar, reachVar);
		} else {
			printf("PARSE ERROR! Undeclared graph identifier %d\n", graphID), exit(1);
			exit(1);
		}
	}
	void readOldMaxFlowConstraint(B& in, Solver& S) {
		if (opt_ignore_theories) {
			skipLine(in);
			return;
		}

		++in;

		int graphID = parseInt(in);
		int s = parseInt(in);
		int t = parseInt(in);
		long flow = parseInt(in);
		int reachVar = parseInt(in) - 1; //note: maximum flow constraint format has been changed since the paper. The order of reachVar and flow after the paper, to allow for non-integer flow constraints.


		if (graphID < 0 || graphID >= graphs.size()) {
			printf("PARSE ERROR! Undeclared graph identifier %d for edge %d\n", graphID, reachVar), exit(1);
		}
		if (reachVar < 0) {
			printf("PARSE ERROR! Edge variables must be >=0, was %d\n", reachVar), exit(1);
		}

		while (reachVar >= S.nVars())
			S.newVar();

		if (graphs[graphID]) {
			graphs[graphID]->maxFlow(s, t, flow, reachVar);
		} else if (graphs_float[graphID]) {
			graphs_float[graphID]->maxFlow(s, t, flow, reachVar);
		} else if (graphs_rational[graphID]) {
			graphs_rational[graphID]->maxFlow(s, t, flow, reachVar);
		} else {
			printf("PARSE ERROR! Undeclared graph identifier %d\n", graphID), exit(1);
			exit(1);
		}
	}
	void readMaxFlowConstraint(B& in, Solver& S, bool inclusive) {
		if (opt_ignore_theories) {
			skipLine(in);
			return;
		}

		++in;

		int graphID = parseInt(in);
		int s = parseInt(in);
		int t = parseInt(in);
		int reachVar = parseInt(in) - 1; //note: switched the order of reachVar and flow after the paper, to allow for non-integer flow constraints in the future...


		if (graphID < 0 || graphID >= graphs.size()) {
			printf("PARSE ERROR! Undeclared graph identifier %d for edge %d\n", graphID, reachVar), exit(1);
		}
		if (reachVar < 0) {
			printf("PARSE ERROR! Edge variables must be >=0, was %d\n", reachVar), exit(1);
		}

		while (reachVar >= S.nVars())
			S.newVar();

		//parse the flow constraint appropriately for the type of the graph:

		if (graphs[graphID]) {
			long flow = parseInt(in);
			graphs[graphID]->maxFlow(s, t, flow, reachVar,inclusive);
		} else if (graphs_float[graphID]) {
			//float can be either a plain integer, or a rational (interpreted at floating point precision) in the form '123/456', or a floating point in decimal format
			double flow = parseDouble(in, tmp);
			skipWhitespace(in);
			if (*in == '/') {
				++in;
				double denom = parseDouble(in, tmp);
				flow /= denom;
			}
			graphs_float[graphID]->maxFlow(s, t, flow, reachVar,inclusive);
		} else if (graphs_rational[graphID]) {
			//rational can be either a plain integer, or a rational in the form '123/456', or a floating point value
			std::stringstream ss;
			skipWhitespace(in);

			while (*in != '\n') {
				ss << (*in);
				++in;
			}

			//first, try to interpret this string as a double:
			try {
				double value = std::stod(ss.str());
				mpq_class flow(value);
				flow.canonicalize();
				graphs_rational[graphID]->maxFlow(s, t, flow, reachVar,inclusive);
			} catch (std::exception& e) {
				//if that fails, attempt to read it in directly as an fraction:
				mpq_class flow(ss.str());
				flow.canonicalize();
				graphs_rational[graphID]->maxFlow(s, t, flow, reachVar,inclusive);
			}

		} else {
			printf("PARSE ERROR! Undeclared graph identifier %d\n", graphID), exit(1);
			exit(1);
		}

	}
	
	void readMinConnectedComponentsConstraint(B& in, Solver& S) {
		if (opt_ignore_theories) {
			skipLine(in);
			return;
		}
		//connected_component_count_lt grachID min_components var is a minimum connected component count constraint, true if the max flow is >= var
		
		++in;
		
		int graphID = parseInt(in);
		int min_components = parseInt(in);
		int reachVar = parseInt(in) - 1;
		if (graphID < 0 || graphID >= graphs.size()) {
			printf("PARSE ERROR! Undeclared graph identifier %d for edge %d\n", graphID, reachVar), exit(1);
		}
		if (reachVar < 0) {
			printf("PARSE ERROR! Edge variables must be >=0, was %d\n", reachVar), exit(1);
		}
		
		while (reachVar >= S.nVars())
			S.newVar();
		
		if (graphs[graphID]) {
			graphs[graphID]->minConnectedComponents(min_components, reachVar);
		} else if (graphs_float[graphID]) {
			graphs_float[graphID]->minConnectedComponents(min_components, reachVar);
		} else if (graphs_rational[graphID]) {
			graphs_rational[graphID]->minConnectedComponents(min_components, reachVar);
		} else {
			printf("PARSE ERROR! Undeclared graph identifier %d\n", graphID), exit(1);
			exit(1);
		}
	}
	
	/* void readSteinerTreeDeclaration(B& in, Solver& S, vec<GraphTheorySolver*> & graphs, vec<vec<SteinerStruct*>> steiners) {
	 if(opt_ignore_theories){
	 skipLine(in);
	 return;
	 }
	 //steiner_tree graphID steinerID
	 //Creates a steiner tree with the given id. Steiner constraints can then refer to that steinerID.
	 if(!eagerMatch(in,"steiner_tree")){
	 printf("PARSE ERROR! Unexpected char: %c\n", *in), exit(1);
	 }


	 int graphID = parseInt(in);
	 int steinerID = parseInt(in);

	 if(graphID <0 || graphID>=graphs.size()){
	 printf("PARSE ERROR! Undeclared graph identifier %d\n",graphID), exit(1);
	 }

	 steiners.growTo(graphs.size());
	 steiners[graphID].growTo(steinerID+1);
	 if(steiners[graphID][steinerID]!=0){
	 printf("PARSE ERROR! Multiple declarations of steiner tree %d for graph %d\n",steinerID, graphID), exit(1);
	 }
	 steiners[graphID][steinerID] = new SteinerStruct{steinerID};
	 }*/
	void readSteinerTreeTerminal(B& in, Solver& S) {
		if (opt_ignore_theories) {
			skipLine(in);
			return;
		}
		//steiner_terminal graphID steinerID node var
		//Declares that the specified node in graphID is a terminal in the steiner tree (conditioned on var being assigned to true)
		
		int graphID = parseInt(in);
		int steinerID = parseInt(in);
		int node = parseInt(in);
		int var = parseInt(in) - 1;
		if (graphID < 0 || graphID >= graphs.size()) {
			printf("PARSE ERROR! Undeclared graph identifier %d \n", graphID), exit(1);
		}
		steiners.growTo(graphs.size());
		steiners[graphID].growTo(steinerID + 1);
		if (!steiners[graphID][steinerID]) {
			steiners[graphID][steinerID] = new SteinerStruct(steinerID);
		}
		steiners[graphID][steinerID]->terminals.push( { node, var });
		
	}
	void readSteinerTreeConstraint(B& in, Solver& S) {
		if (opt_ignore_theories) {
			skipLine(in);
			return;
		}
		//steiner_minweight graphID steinerID maxweight var is a minimum spanning tree weight constraint: var is true iff the minimum steiner tree of the graph is <= maxweight
		
		int graphID = parseInt(in);
		int steinerID = parseInt(in);
		int maxweight = parseInt(in);
		int var = parseInt(in) - 1;
		if (graphID < 0 || graphID >= graphs.size()) {
			printf("PARSE ERROR! Undeclared graph identifier %d \n", graphID), exit(1);
		}
		steiners.growTo(graphs.size());
		steiners[graphID].growTo(steinerID + 1);
		if (!steiners[graphID][steinerID]) {
			steiners[graphID][steinerID] = new SteinerStruct { steinerID };
		}
		steiners[graphID][steinerID]->weight_constraints.push( { maxweight, var });
	}
	
	void readBV(B& in, Solver& S){
		if (opt_ignore_theories) {
			skipLine(in);
			return;
		}

	}

	void readPB(B & in, vec<Lit> & lits, vec<int> & weights, Solver & S, PbTheory * pb) {
		if (opt_ignore_theories) {
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
		
		//pb_lt <size> lit1 lit2 ... [0 | <size>weight1 weight2 ...] 'op' total 2 headLiteral
		
		lits.clear();
		weights.clear();
		int size = parseInt(in);
		if (size <= 0) {
			printf("PARSE ERROR! Empty PB clause\n"), exit(1);
		}
		
		for (int i = 0; i < size; i++) {
			int parsed_lit = parseInt(in);
			if (parsed_lit == 0)
				break;
			int v = abs(parsed_lit) - 1;
			while (v >= S.nVars())
				S.newVar();
			lits.push((parsed_lit > 0) ? mkLit(v) : ~mkLit(v));
		}
		
		int wsize = parseInt(in);
		if (wsize != 0 && wsize != size) {
			printf("PARSE ERROR! Number of weights must either be the same as the size of the clause, or 0.\n"), exit(1);
		}
		for (int i = 0; i < wsize; i++) {
			int parsed_weight = parseInt(in);
			weights.push(parsed_weight);
		}
		if (wsize == 0) {
			for (int i = 0; i < size; i++)
				weights.push(1);
		}
		skipWhitespace(in);
		PbTheory::PbType op = PbTheory::PbType::EQ;
		//read the operator:
		if (*in == '<') {
			++in;
			if (*in == '=') {
				++in;
				op = PbTheory::PbType::LE;
			} else {
				op = PbTheory::PbType::LT;
			}
		} else if (*in == '>') {
			++in;
			if (*in == '=') {
				++in;
				op = PbTheory::PbType::GE;
			} else {
				op = PbTheory::PbType::GT;
			}
		} else if (*in == '!') {
			++in;
			if (*in != '=') {
				printf("PARSE ERROR! Unexpected char: %c\n", *in), exit(1);
			}
			++in;
			op = PbTheory::PbType::NE;
		} else if (*in == '=') {
			++in;
			op = PbTheory::PbType::EQ;
		}else{
			printf("PARSE ERROR! Bad PB constraint\n"), exit(1);
		}
		int comparison = parseInt(in);
		int type = parseInt(in);
		Lit head = lit_Undef;
		bool oneSided = false;
		if (type == 0) {
			//done
			
		} else {
			if (type == 1) {
				oneSided = true;
			}
			int parsed_lit = parseInt(in);
			
			int v = abs(parsed_lit) - 1;
			while (v >= S.nVars())
				S.newVar();
			head = (parsed_lit > 0) ? mkLit(v) : ~mkLit(v);
		}
		assert(lits.size() == weights.size());
		pb->addConstraint(lits, weights, comparison, head, op,
				oneSided ? PbTheory::ConstraintSide::Upper : PbTheory::ConstraintSide::Both);
	}
	
public:
	GraphParser(bool precise, BVTheorySolver<long>*& bvTheory) :
			precise(precise),bvTheory(bvTheory) {
		
	}
	bool parseLine(B& in, Solver& S) {
		
		skipWhitespace(in);
		if (*in == EOF)
			return false;
		else if (*in == 'c') {
			//just a comment
			return false;
		} else if (match(in, "digraph")) {
			skipWhitespace(in);
			if (match(in, "int")) {
				readDiGraph(in, GraphType::INTEGER, S);
			} else if (match(in, "float")) {
				readDiGraph(in, GraphType::FLOAT, S);
			} else if (match(in, "rational")) {
				readDiGraph(in, GraphType::RATIONAL, S);
			} else {
				//assume the graph is integer
				readDiGraph(in, GraphType::INTEGER, S);
			}
			
			skipWhitespace(in);
			//if(*in=='d'){
			//for now, only digraphs are supported
			
			//}else{
			//	printf("PARSE ERROR! Unexpected char: %c\n", *in), exit(1);
			//}
			return true;
		}else if (match(in, "edge_bv")) {
			count++;
			readEdgeBV(in, S);
			return true;
		} else if (match(in, "edge")) {
			count++;
			readEdge(in, S);
			return true;
		} else if (match(in, "weighted_edge")) {
			count++;
			readEdge(in, S);
			return true;
		}else if (match(in, "reach")) {
			readReach(in, S);
			return true;
		} else if (match(in, "distance_lt")) {
			readDistance(in, S);
			return true;
		} else if (match(in, "distance_leq")) {
			readDistance(in, S, true);
			return true;
		} else if (match(in, "weighted_distance_lt")) {
			readShortestPath(in, S);
			return true;
		} else if (match(in, "weighted_distance_leq")) {
			readShortestPath(in, S, true);
			return true;
		}else if (match(in, "mst_weight_leq")) {
			readMinSpanningTreeConstraint(in, S,true);
			return true;
		}else if (match(in, "mst_weight_lt")) {
			readMinSpanningTreeConstraint(in, S,false);
			return true;
		}else if (match(in, "mst_edge")) {
			readMinSpanningTreeEdgeConstraint(in, S);
			return true;
		} else if (match(in, "maximum_flow_geq")) {
			readMaxFlowConstraint(in, S,true);
			return true;
		} else if (match(in, "maximum_flow_gt")) {
			readMaxFlowConstraint(in, S,false);
			return true;
		} else if (match(in, "max_flow_gt")) {
			//for compatibility with an old file format (don't use this!)
			readOldMaxFlowConstraint(in, S);
			return true;
		}else if (match(in, "connected_component_count_lt")) {
			readMinConnectedComponentsConstraint(in, S);
			return true;
		}else if (match(in, "acyclic")) {
			//A _directed_ acyclic graph constraint
			readAcyclic(in, S,true);
			return true;
		}else if (match(in, "forest")) {
			//An _undirected_ acyclic graph constraint
			readAcyclic(in, S,false);
			return true;
		}else if (match(in, "pb_lt")) {

			if (!pbtheory) {
				pbtheory = new PbTheory(&S);
				S.addTheory(pbtheory);
			}
			readPB(in, lits, weights, S, pbtheory);
			return true;
		}
		return false;
	}
	
	void implementConstraints(Solver & S) {
		for (int i = 0; i < graphs.size(); i++) {
			if(graphs[i])
				graphs[i]->setBVTheory(bvTheory);
		}
	/*	for (int i = 0; i < graphs_float.size(); i++) {
			if(graphs_float[i])
				graphs_float[i]->setComparator(comparison);
		}
		for (int i = 0; i < graphs_rational.size(); i++) {
			if(graphs_rational[i])
				graphs_rational[i]->setComparator(comparison);
		}*/
		//not really implemented, yet!
	/*	for (int gid = 0; gid < steiners.size(); gid++) {
			for (auto & steiner : steiners[gid]) {
				if (steiner) {
					graphs[gid]->steinerTree(steiner->terminals, steiner->id);
					for (auto & weight : steiner->weight_constraints) {
						graphs[gid]->addSteinerWeightConstraint(steiner->id, weight.first, weight.second);
					}
					delete (steiner);
				}
			}
		}*/
		for (auto & e:bvedges){
			graphs[e.graphID]->newEdgeBV(e.from, e.to, e.edgeVar, e.bvID);
		}
		for (int i = 0; i < graphs.size(); i++) {
			if (graphs[i])
				graphs[i]->implementConstraints();
		}
		for (int i = 0; i < graphs_float.size(); i++) {
			if (graphs_float[i])
				graphs_float[i]->implementConstraints();
		}
		for (int i = 0; i < graphs_rational.size(); i++) {
			if (graphs_rational[i])
				graphs_rational[i]->implementConstraints();
		}
		if (pbtheory)
			pbtheory->implementConstraints();
	}

	
};

//=================================================================================================
}
;

#endif /* GRAPH_PARSER_H_ */
