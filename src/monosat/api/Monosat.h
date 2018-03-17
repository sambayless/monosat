/****************************************************************************************[Solver.h]
 The MIT License (MIT)

 Copyright (c) 2016, Sam Bayless

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

#ifndef MONOSAT_H_
#define MONOSAT_H_
//Monosat library interface, in C
#include <stdint.h>
#ifdef __cplusplus
#include "monosat/utils/ParseUtils.h"
#include "monosat/utils/Options.h"
#include "monosat/core/Solver.h"
#include "monosat/simp/SimpSolver.h"
#include "monosat/graph/GraphTheory.h"
#include "monosat/geometry/GeometryTheory.h"
#include "monosat/fsm/FSMTheory.h"
#include "monosat/pb/PbTheory.h"
#include "monosat/amo/AMOTheory.h"
#include "Monosat.h"
#include "monosat/core/Dimacs.h"
#include "monosat/bv/BVParser.h"
#include "monosat/graph/GraphParser.h"
#include "monosat/amo/AMOParser.h"
#include "monosat/core/Optimize.h"
#include "monosat/pb/PbSolver.h"
#include "monosat/routing/FlowRouter.h"
extern "C"
{
typedef Monosat::SimpSolver *  SolverPtr;
typedef Monosat::GraphTheorySolver<int64_t> * GraphTheorySolver_long;
typedef Monosat::GraphTheorySolver<double>*  GraphTheorySolver_double;
typedef Monosat::BVTheorySolver<int64_t>* BVTheoryPtr;
typedef Monosat::FSMTheorySolver * FSMTheorySolverPtr;
typedef Monosat::FlowRouter<int64_t> * FlowRouterPtr;
typedef int64_t Weight;
#else
typedef void * SolverPtr;
typedef void *  BVTheoryPtr;
typedef void *  GraphTheorySolver_long;
typedef void *  GraphTheorySolver_double;
typedef void *  FSMTheorySolverPtr;
typedef void *  FlowRouterPtr;
typedef int Var;
typedef int64_t Weight;
#endif


  inline int varToLit(int var, bool negated){
	  return toInt(mkLit(var,negated));
  }

 inline  int litToVar(int lit){
	  return lit/2;
  }
 const char * getVersion();
  SolverPtr newSolver();
  SolverPtr newSolver_arg(const char*argv);
#ifndef JNA
  //Java Native Access sometimes has problems with string arrays
  SolverPtr newSolver_args(int argc, char**argv);
#endif

  void deleteSolver (SolverPtr S);
  //If set, dump constraints to this file (as they are asserted in the solver)
  void setOutputFile(SolverPtr S,const char * output);
  void readGNF(SolverPtr S, const char  * filename);

  bool solve(SolverPtr S);
  bool solveAssumptions(SolverPtr S,int * assumptions, int n_assumptions);
  //Solve under assumptions, and also minimize a set of BVs (in order of precedence)
  //bool solveAssumptions_MinBVs(SolverPtr S,int * assumptions, int n_assumptions, int * minimize_bvs, int n_minimize_bvs);

  //Sets the (approximate) time limit in seconds before returning l_Undef from solveLimited; ignored by solve(). Set to <0 to disable time limit.
  void setTimeLimit(SolverPtr S,int seconds);
  //Sets the (approximate) memory limit in megabytes before returning l_Undef from solveLimited; ignored by solve(). Set to <0 to disable memory limit.
  void setMemoryLimit(SolverPtr S,int mb);
  //Sets the maximum number of (additional) conflicts allowed in the solver before returning l_Undef from solveLimited; ignored by solve(). Set to <0 to disable conflict limit.
  void setConflictLimit(SolverPtr S,int num_conflicts);
  //Sets the maximum number of (additional) propagations allowed in the solver before returning l_Undef from solveLimited; ignored by solve(). Set to <0 to disable propagation limit.
  void setPropagationLimit(SolverPtr S,int num_propagations);

  //Returns 0 for satisfiable, 1 for proved unsatisfiable, 2 for failed to find a solution (within any resource limits that have been set)
  int solveLimited(SolverPtr S);
  //Returns 0 for satisfiable, 1 for proved unsatisfiable, 2 for failed to find a solution (within any resource limits that have been set)
  int solveAssumptionsLimited(SolverPtr S,int * assumptions, int n_assumptions);

  //Solve under assumptions, and also minimize a set of BVs (in order of precedence)
  //Returns 0 for satisfiable, 1 for proved unsatisfiable, 2 for failed to find a solution (within any resource limits that have been set)
  //int solveAssumptionsLimited_MinBVs(SolverPtr S,int * assumptions, int n_assumptions, int * minimize_bvs, int n_minimize_bvs);

  bool lastSolutionWasOptimal(SolverPtr S);

  //If the last solution was unsat, then this get the 'conflict clause' produced by the solver (a subset of the assumptions which are sufficient to cause the instance to be UNSAT).
  //Fills the given pointer with the first max_store_size literals of the conflict clause, and returns the number of literals in the conflict clause. Set store_clause to null and max_store_size to 0 to find the size of the conflict clause
  //Returns -1 if the solver has no conflict clause from the most recent solve() call (because that call was not UNSAT)
  int getConflictClause(SolverPtr S, int * store_clause, int max_store_size);

  void backtrack(SolverPtr S);
  int newVar(SolverPtr S);
  void setDecisionVar(SolverPtr S,int var,bool decidable);
  bool isDecisionVar(SolverPtr S,int var);

  //Static, lexicographic heuristic. Larger values are higher priority (decided first); default priority is 0
  void setDecisionPriority(SolverPtr S,int var, int priority);
  int getDecisionPriority(SolverPtr S,int var);
  // Which polarity the decision heuristic should use for a variable (by default).
  void setDecisionPolarity(SolverPtr S,Var v, bool b);
  bool getDecisionPolarity(SolverPtr S,Var v);

  //The solver will (sometimes) instantiate an arbitrary true literal for use as a constant.
  //Call this method to a) force that literal to be instantiate, and b) get that literal.
  int true_lit(SolverPtr S);

  //Prevents this literal from being simplified by the preprocessor
  bool disallowLiteralSimplification(SolverPtr S, int lit);

  //permanently disable SAT-based preprocessing in this solver
  void disablePreprocessing(SolverPtr S);
  int nVars(SolverPtr S);
  int nClauses(SolverPtr S);
  int nBitvectors(SolverPtr S,BVTheoryPtr bv);

  bool addClause(SolverPtr S,int * lits, int n_lits);
  bool addUnitClause(SolverPtr S,int lit);
  bool addBinaryClause(SolverPtr S,int lit1, int lit2);
  bool addTertiaryClause(SolverPtr S,int lit1, int lit2, int lit3);

  //remove any optimization objectives from the solver
  void clearOptimizationObjectives(SolverPtr S);

  void maximizeBV(SolverPtr S,  BVTheoryPtr bv, int bvID);
  void minimizeBV(SolverPtr S,  BVTheoryPtr bv, int bvID);
  void maximizeLits(SolverPtr S, int * lits, int n_lits);
  void minimizeLits(SolverPtr S, int * lits, int n_lits);
  void maximizeWeightedLits(SolverPtr S, int * lits, int * weights, int n_lits);
  void minimizeWeightedLits(SolverPtr S, int * lits, int * weights, int n_lits);

//theory interface for bitvectors
  BVTheoryPtr initBVTheory(SolverPtr S);
  int newBitvector_const(SolverPtr S, BVTheoryPtr bv, int bvWidth, Weight constval);
  int newBitvector_anon(SolverPtr S, BVTheoryPtr bv, int bvWidth);
  int newBitvector(SolverPtr S, BVTheoryPtr bv, int * bits, int n_bits);
  int bv_width(SolverPtr S, BVTheoryPtr  bv,int bvID);
  int newBVComparison_const_lt(SolverPtr S, BVTheoryPtr bv, int bvID, Weight weight);
  int newBVComparison_bv_lt(SolverPtr S, BVTheoryPtr bv, int bvID, int compareID);
  int newBVComparison_const_leq(SolverPtr S, BVTheoryPtr bv, int bvID, Weight weight);
  int newBVComparison_bv_leq(SolverPtr S, BVTheoryPtr bv, int bvID, int compareID);
  int newBVComparison_const_gt(SolverPtr S, BVTheoryPtr bv, int bvID, Weight weight);
  int newBVComparison_bv_gt(SolverPtr S, BVTheoryPtr bv, int bvID, int compareID);
  int newBVComparison_const_geq(SolverPtr S, BVTheoryPtr bv, int bvID, Weight weight);
  int newBVComparison_bv_geq(SolverPtr S, BVTheoryPtr bv, int bvID, int compareID);

  //Convert the specified bitvector, as well as any other bitvectors in its cone of influence, into pure CNF
  void bv_bitblast(SolverPtr S, BVTheoryPtr bv,int bvID);

  void bv_concat( SolverPtr S, BVTheoryPtr bv,int aID, int bID, int resultID);
  void bv_slice( SolverPtr S, BVTheoryPtr bv,int aID, int lower, int upper, int resultID);
  void bv_not( SolverPtr S, BVTheoryPtr bv,int bvaID, int bvResultID);
  void bv_and( SolverPtr S, BVTheoryPtr bv,int bvaID, int bvbID, int bvResultID);
  void bv_nand( SolverPtr S, BVTheoryPtr bv,int bvaID, int bvbID, int bvResultID);
  void bv_or( SolverPtr S, BVTheoryPtr bv,int bvaID, int bvbID, int bvResultID);
  void bv_nor( SolverPtr S, BVTheoryPtr bv,int bvaID, int bvbID, int bvResultID);
  void bv_xor( SolverPtr S, BVTheoryPtr bv,int bvaID, int bvbID, int bvResultID);
  void bv_xnor( SolverPtr S, BVTheoryPtr bv,int bvaID, int bvbID, int bvResultID);

  void bv_ite( SolverPtr S, BVTheoryPtr bv, int condition_lit,int bvThenID, int bvElseID, int bvResultID);

  void bv_addition( SolverPtr S, BVTheoryPtr bv, int bvID1, int bvID2, int resultID);
  void bv_subtraction( SolverPtr S, BVTheoryPtr bv, int bvID1, int bvID2, int resultID);
void bv_multiply(SolverPtr S, BVTheoryPtr bv, int bvID1, int bvID2, int resultID);
void bv_divide(SolverPtr S, BVTheoryPtr bv, int bvID1,  int bvID2, int resultID);
  void bv_min(SolverPtr S, BVTheoryPtr bv,  int* args,int n_args,int resultID);
  void bv_max(SolverPtr S, BVTheoryPtr bv,  int* args,int n_args,int resultID);
  void bv_popcount(SolverPtr S, BVTheoryPtr bv,  int* args,int n_args, int resultID);
void bv_unary(SolverPtr S, BVTheoryPtr bv, int * args, int n_args, int resultID);

  //simple at-most-one constraint: asserts that at most one of the set of variables (NOT LITERALS) may be true.
  //for small numbers of variables, consider using a direct CNF encoding instead
  void at_most_one(SolverPtr S, int * vars, int n_vars);

  void assertPB_lt(SolverPtr S, int rhs, int n_args, int * literals, int * coefficients);
  void assertPB_leq(SolverPtr S, int rhs, int n_args, int * literals, int * coefficients);
  void assertPB_eq(SolverPtr S, int rhs, int n_args, int * literals, int * coefficients);
  void assertPB_geq(SolverPtr S, int rhs, int n_args, int * literals, int * coefficients);
  void assertPB_gt(SolverPtr S, int rhs, int n_args, int * literals, int * coefficients);
  //Convert any pb constraints in the solver into cnf (will be called automatically before solve())
  void flushPB(SolverPtr S);
//theory interface for graphs

  GraphTheorySolver_long newGraph(SolverPtr S);

  int newNode(SolverPtr S,GraphTheorySolver_long G);
  int newEdge(SolverPtr S, GraphTheorySolver_long G,int from,int  to,  Weight weight);
  int newEdge_double(SolverPtr S, GraphTheorySolver_double G,int from,int  to,  double weight);
  int newEdge_bv(SolverPtr S, GraphTheorySolver_long G,int from,int  to, int bvID);
  int nNodes(SolverPtr S,GraphTheorySolver_long G);
  int nEdges(SolverPtr S,GraphTheorySolver_long G);
  int reaches(SolverPtr S,GraphTheorySolver_long G,int from, int to);
  int shortestPathUnweighted_lt_const(SolverPtr S,GraphTheorySolver_long G,int from, int to, int steps);
  int shortestPathUnweighted_leq_const(SolverPtr S,GraphTheorySolver_long G,int from, int to, int steps);
  int shortestPath_lt_const(SolverPtr S,GraphTheorySolver_long G,int from, int to, Weight dist);
  int shortestPath_leq_const(SolverPtr S,GraphTheorySolver_long G,int from, int to, Weight dist);
  int shortestPath_lt_bv(SolverPtr S,GraphTheorySolver_long G,int from, int to, int bvID);
  int shortestPath_leq_bv(SolverPtr S,GraphTheorySolver_long G,int from, int to, int bvID);
  int maximumFlow_geq(SolverPtr S,GraphTheorySolver_long G,int source, int sink, Weight weight);
  int maximumFlow_gt(SolverPtr S,GraphTheorySolver_long G,int source, int sink, Weight weight);
  int maximumFlow_geq_bv(SolverPtr S,GraphTheorySolver_long G,int source, int sink, int bvID);
  int maximumFlow_gt_bv(SolverPtr S,GraphTheorySolver_long G,int source, int sink, int bvID);
  int minimumSpanningTree_leq(SolverPtr S,GraphTheorySolver_long G, Weight weight);
  int minimumSpanningTree_lt(SolverPtr S,GraphTheorySolver_long G,int source, int sink, Weight weight);
  int acyclic_undirected(SolverPtr S,GraphTheorySolver_long G);
  int acyclic_directed(SolverPtr S,GraphTheorySolver_long G);
  void newEdgeSet(SolverPtr S,GraphTheorySolver_long G,int * edges, int n_edges, bool enforceEdgeAssignment);

  //this enables a heuristic on this graph, from the RUC paper, which sets assigned edges to zero weight, to encourage edge-reuse in solutions
  void graph_setAssignEdgesToWeight(SolverPtr S,GraphTheorySolver_long G, int64_t weight);
  //flow routing interface

  FlowRouterPtr createFlowRouting(SolverPtr S,GraphTheorySolver_long G, int sourceNode,int destNode,int maxflowLit);
  void addRoutingNet(SolverPtr S,GraphTheorySolver_long G, FlowRouterPtr router, int disabledEdge, int n_members, int * edge_lits, int * reach_lits);


  //theory interface for finite state machines
  FSMTheorySolverPtr initFSMTheory(SolverPtr S);
  int newFSM(SolverPtr S,FSMTheorySolverPtr fsmTheory, int inputAlphabet, int outputAlphabet);
  int newState(SolverPtr S,FSMTheorySolverPtr fsmTheory, int fsmID);

  int newTransition(SolverPtr S,FSMTheorySolverPtr fsmTheory, int fsmID, int fromNode, int toNode,int inputLabel, int outputLabel);
  int newString(SolverPtr S,FSMTheorySolverPtr fsmTheory, int * str,int len);
  int fsmAcceptsString(SolverPtr S,FSMTheorySolverPtr fsmTheory, int fsmID, int startNode, int acceptNode,int stringID);
  int fsmCompositionAccepts(Monosat::SimpSolver * S, Monosat::FSMTheorySolver *  fsmTheory,   int fsmGeneratorID,int fsmAcceptorID, int gen_startNode, int gen_acceptNode, int acceptor_startNode, int acceptor_acceptNode, int stringID);

  //model query
  //For a given literal (not variable!), returns 0 for true, 1 for false, 2 for unassigned.
  int getModel_Literal(SolverPtr S,int lit);
  //Get an assignment to a bitvector in the model. The model may find a range of satisfying assignments to the bitvector;
  //If getMaximumValue is true, this function returns the maximum satisfying assignment to the bitvector in the model; else it returns the smallest.
  Weight getModel_BV(SolverPtr S, BVTheoryPtr bv, int bvID, bool getMaximumValue);
  //graph queries:
  //maxflow_literal is the literal (not variable!) that is the atom for the maximum flow query
  Weight getModel_MaxFlow(SolverPtr S,GraphTheorySolver_long G,int maxflow_literal);
  //maxflow_literal is the literal (not variable!) that is the atom for the maximum flow query
  Weight getModel_EdgeFlow(SolverPtr S,GraphTheorySolver_long G,int maxflow_literal, int edgeLit);
  Weight getModel_AcyclicEdgeFlow(SolverPtr S,GraphTheorySolver_long G,int maxflow_literal, int edgeLit);

  Weight getModel_MinimumSpanningTreeWeight(SolverPtr S,GraphTheorySolver_long G,int spanning_tree_literal);
  int getModel_Path_Nodes_Length(SolverPtr S,GraphTheorySolver_long G,int reach_or_distance_literal);
  int getModel_Path_Nodes(SolverPtr S,GraphTheorySolver_long G,int reach_or_distance_literal, int store_length, int * store);

  int getModel_Path_EdgeLits_Length(SolverPtr S,GraphTheorySolver_long G,int reach_or_distance_literal);
  int getModel_Path_EdgeLits(SolverPtr S,GraphTheorySolver_long G,int reach_or_distance_literal, int store_length, int * store);


#ifdef __cplusplus
}
#endif

#endif /* MONOSAT_H_ */
