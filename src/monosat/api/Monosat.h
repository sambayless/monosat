/**************************************************************************************************
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
#include <stdio.h>
#ifdef __cplusplus
#include "monosat/utils/ParseUtils.h"
#include "monosat/utils/Options.h"
#include "monosat/core/Solver.h"
#include "monosat/simp/SimpSolver.h"
#include "monosat/graph/GraphTheory.h"
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
#include <stdbool.h>
typedef void * SolverPtr;
typedef void *  BVTheoryPtr;
typedef void *  GraphTheorySolver_long;
typedef void *  GraphTheorySolver_double;
typedef void *  FSMTheorySolverPtr;
typedef void *  FlowRouterPtr;
typedef int Var;
typedef int64_t Weight;
#endif


  int varToLit(int var, bool negated);

  inline int litToVar(int lit){
	  return lit/2;
  }
 const char * getVersion(void);
  SolverPtr newSolver(void);
  SolverPtr newSolver_arg(const char*argv);
#ifndef JNA
  //Java Native Access sometimes has problems with string arrays
  SolverPtr newSolver_args(int argc, char**argv);
#endif

  void deleteSolver (SolverPtr S);
  //Return true if the solver has not yet proven a the formula UNSAT
  bool ok(SolverPtr S);
  //If set, dump constraints to this file (as they are asserted in the solver)
  void setOutputFile(SolverPtr S,const char * output);
  //Load a GNF, and run any embedded solve/optimize calls
  void readGNF(SolverPtr S, const char  * filename);
  //Load a GNF, but ignore any embedded solve/optimize calls
  void loadGNF(SolverPtr S, const char  * filename);
  //flush constraints to file
  void flushFile (SolverPtr S);
  //stop writing constraints to file, and close the file (if any)
  //this will be called automatically if the solver is deleted
  void closeFile (SolverPtr S);

  bool solve(SolverPtr S);
  bool solveAssumptions(SolverPtr S,int * assumptions, int n_assumptions);
  //Solve under assumptions, and also minimize a set of BVs (in order of precedence)
  //bool solveAssumptions_MinBVs(SolverPtr S,int * assumptions, int n_assumptions, int * minimize_bvs, int n_minimize_bvs);

  //Sets the (approximate) time limit in seconds before returning l_Undef from solveLimited; ignored by solve(). Set to <0 to disable time limit.
  void setTimeLimit(SolverPtr S,int seconds);

  //Sets the maximum number of (additional) conflicts allowed in the solver before returning l_Undef from solveLimited; ignored by solve(). Set to <0 to disable conflict limit.
  void setConflictLimit(SolverPtr S,int num_conflicts);
  //Sets the maximum number of (additional) propagations allowed in the solver before returning l_Undef from solveLimited; ignored by solve(). Set to <0 to disable propagation limit.
  void setPropagationLimit(SolverPtr S,int num_propagations);

  uint64_t nConflicts(SolverPtr S);

  uint64_t nPropagations(SolverPtr S);


  //Returns 0 for satisfiable, 1 for proved unsatisfiable, 2 for failed to find a solution (within any resource limits that have been set)
  int solveLimited(SolverPtr S);
  //Returns 0 for satisfiable, 1 for proved unsatisfiable, 2 for failed to find a solution (within any resource limits that have been set)
  int solveAssumptionsLimited(SolverPtr S,int * assumptions, int n_assumptions);

  //Solve under assumptions, and also minimize a set of BVs (in order of precedence)
  //Returns 0 for satisfiable, 1 for proved unsatisfiable, 2 for failed to find a solution (within any resource limits that have been set)
  //int solveAssumptionsLimited_MinBVs(SolverPtr S,int * assumptions, int n_assumptions, int * minimize_bvs, int n_minimize_bvs);

  //If the solver is operating with limits on the number of conflicts, memory usage, or runtime (see the setXXXXXLimit() methods), then optimization queries
  //are not guaranteed to produce an optimal (or locally optimal) solution.
  //This function returns false if the last optimization call to the solver produced a non-optimal result due to such a limit, and true otherwise.
  bool lastSolutionWasOptimal(SolverPtr S);

  //If the last solution was UNSAT, and assumptions were enforced in the solver, then this get the 'conflict clause' produced by the solver (a subset of the assumptions which are sufficient to cause the instance to be UNSAT).
  //Fills the given pointer with the first max_store_size literals of the conflict clause, and returns the number of literals in the conflict clause. Set store_clause to null and max_store_size to 0 to find the size of the conflict clause
  //Returns -1 if the solver has no conflict clause from the most recent solve() call (because that call was not UNSAT)
  int getConflictClause(SolverPtr S, int * store_clause, int max_store_size);

  //Reduce the given set of mutually UNSAT assumptions into a (locally) minimal-sized set of assumptions that are still UNSAT.
  //If the assumptions are not UNSAT, this method does returns -1.
  //Else, it returns the number of literals in the unsat core, which it will also store in the first number_literals
  //entries in unsat_assumptions
  //Note: This method will make repeated (and potentially expensive) calls to the SAT solver to attempt to remove literals from the
  //set of assumptions.
  //Note also that if any of the setXXXXXLimits() are applied to the solver, this may not produce a locally minimal unsat core.
  int minimizeUnsatCore(SolverPtr S, int * unsat_assumptions, int n_lits);

  //After UNSAT solve calls with assumptions, the solver will find a 'conflict clause' consisting of a subset of the assumptions
  //which are sufficient to make the solver UNSAT (see getConflictClause).
  //Normally, the conflict clause is produced as a side effect of proving the query unsat, with the solver only removing literals
  //from the conflict clause on a best effort basis.
  //This method will make repeated (and potentially expensive) calls to the SAT solver to attempt to remove further literals from
  //the conflict clause.
  //Afterward, the conflict clause can be obtained using getConflictClause().
  //NOTE: this function may be expensive, is not required to get a conflict clause; getConflictClause() can be used after any unsat call with assumptions, even
  //without calling minimizeConflictClause().
  //Note also that if any of the setXXXXXLimits() are applied to the solver, this may not produce a locally minimal conflict clause.
  void minimizeConflictClause(SolverPtr S);

  //Restart the solver after a solution (this will be done automatically in most circumstances, and so is not usually called)
  //This only undoes the decisions made internally by the solver, without removing its variables or constraints.
  void backtrack(SolverPtr S);

  //Create a new variable in the sat solver, and return that variable.
  //Note: Most functions expect a literal, not a variable.
  //To convert a variable into a literal, use varToLit(variable)
  int newVar(SolverPtr S);

  //Create a new named variable, only if this name is unique and valid (or if it is empty).
  //If the name is not empty, and is non-unique or invalid, an exception will be thrown and no variable will be created at all.
  int newNamedVar(SolverPtr S,const char  * varname);

  //Associate a unique name with this variable.
  //varname must consist of printable ascii characters, and may not contain a newline.
  //If varname is null or length-0, then this will remove any existing name.
  //If varname is non-unique, or if this variable already has a name, throw an excpetion
  void setVariableName(SolverPtr S, int variable, const char  * varname);
  bool variableHasName(SolverPtr S, int variable);
  bool hasVariableWithName(SolverPtr S, const char * name);
  const char * getVariableName(SolverPtr S, int variable);
  int getVariable(SolverPtr S, const char * varname);

  //Get the nth named variable (in the order that the named variables were assigned names)
  int getNamedVariableN(SolverPtr S, int n);
  //Get the number of named variables in the solver
  int nNamedVariables(SolverPtr S);

  //Release this literal back to the sat solver, so that its variable can be eventually reused (after the next backtrack to 0).
  //The literal will be assigned to true in this process.
  void releaseLiteral(SolverPtr S, int literal);
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
  bool disallowLiteralSimplification(SolverPtr S, int var);

  //permanently disable SAT-based preprocessing in this solver
  void disablePreprocessing(SolverPtr S);
  int nVars(SolverPtr S);
  int nClauses(SolverPtr S);
  int nLearnedClauses(SolverPtr S);
  int nBitvectors(SolverPtr S,BVTheoryPtr bv);

  bool addClause(SolverPtr S,int * lits, int n_lits);
  bool addUnitClause(SolverPtr S,int lit);
  bool addBinaryClause(SolverPtr S,int lit1, int lit2);
  bool addTertiaryClause(SolverPtr S,int lit1, int lit2, int lit3);
  //Add n_pairs binary clauses: (first_args[0] OR second_args[0]) AND (first_args[1] OR second_args[1])...
  void addBinaryClauses(SolverPtr S,int * first_args,int * second_args, int n_pairs);
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

  void setBitvectorName(SolverPtr S, BVTheoryPtr bv, int bvID, const char * name);
  const char * getBitvectorName(SolverPtr S, BVTheoryPtr bv, int bvID);
  int getBitvector(SolverPtr S, BVTheoryPtr bv, const char * name);

  //True if this bit vector has a non-empty name.
  bool bitvectorHasName(SolverPtr S, BVTheoryPtr bv, int bvID);
  //True if there exists a bitvector with the given name.
  bool hasBitvectorWithName(SolverPtr S, BVTheoryPtr bv, const char * name);


  //Get the number of named bitvectors in the solver
  int nNamedBitvectors(SolverPtr S, BVTheoryPtr bv);
  //Get the nth named bitvector in the solver (in the order they were named)
  int getNamedBitvectorN(SolverPtr S, BVTheoryPtr bv,int n);

  int bv_width(SolverPtr S, BVTheoryPtr  bv,int bvID);
  int bv_nBits(SolverPtr S, BVTheoryPtr  bv,int bvID);
  int bv_bit(SolverPtr S, BVTheoryPtr  bv,int bvID, int bit);



  int newBVComparison_const_lt(SolverPtr S, BVTheoryPtr bv, int bvID, Weight weight);
  int newBVComparison_bv_lt(SolverPtr S, BVTheoryPtr bv, int bvID, int compareID);
  int newBVComparison_const_leq(SolverPtr S, BVTheoryPtr bv, int bvID, Weight weight);
  int newBVComparison_bv_leq(SolverPtr S, BVTheoryPtr bv, int bvID, int compareID);
  int newBVComparison_const_gt(SolverPtr S, BVTheoryPtr bv, int bvID, Weight weight);
  int newBVComparison_bv_gt(SolverPtr S, BVTheoryPtr bv, int bvID, int compareID);
  int newBVComparison_const_geq(SolverPtr S, BVTheoryPtr bv, int bvID, Weight weight);
  int newBVComparison_bv_geq(SolverPtr S, BVTheoryPtr bv, int bvID, int compareID);
  int newBVComparison_const_eq(SolverPtr S, BVTheoryPtr bv, int bvID, Weight weight);
  int newBVComparison_bv_eq(SolverPtr S, BVTheoryPtr bv, int bvID, int compareID);
  int newBVComparison_const_neq(SolverPtr S, BVTheoryPtr bv, int bvID, Weight weight);
  int newBVComparison_bv_neq(SolverPtr S, BVTheoryPtr bv, int bvID, int compareID);

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
  //Create a new graph with the given name (or an unamed graph, if name is null or empty).
  //If name is non-null and not empty, then it must be unique.
  GraphTheorySolver_long newGraph_Named(SolverPtr S, const char * name, int bitwidth);
  //If there exists a graph in the solver with the given name, return a pointer to it.
  GraphTheorySolver_long getGraph(SolverPtr S, const char * name);

  const char * getGraphName(SolverPtr S, GraphTheorySolver_long G);


  int getGraphWidth(SolverPtr S, GraphTheorySolver_long G);
  int newNode(SolverPtr S,GraphTheorySolver_long G);
  int newNode_Named(SolverPtr S,GraphTheorySolver_long G,  const char * name);
  bool hasNamedNode(SolverPtr S,GraphTheorySolver_long G,  const char * name);
  const char * getNodeName(SolverPtr S, GraphTheorySolver_long G, int nodeID);
  int newEdge(SolverPtr S, GraphTheorySolver_long G,int from,int  to,  Weight weight);
  int newEdge_double(SolverPtr S, GraphTheorySolver_double G,int from,int  to,  double weight);
  int newEdge_bv(SolverPtr S, GraphTheorySolver_long G,int from,int  to, int bvID);
  int nNodes(SolverPtr S,GraphTheorySolver_long G);
  int nEdges(SolverPtr S,GraphTheorySolver_long G);
  int getEdgeLiteralN(SolverPtr S,GraphTheorySolver_long G, int n);
  int getEdge_to(SolverPtr S,GraphTheorySolver_long G, int edgeLit);
  int getEdge_from(SolverPtr S,GraphTheorySolver_long G, int edgeLit);
  Weight getEdge_weight_const(SolverPtr S,GraphTheorySolver_long G, int edgeLit);
  int getEdge_weight_bv(SolverPtr S,GraphTheorySolver_long G, int edgeLit);
  bool edgeHasBVWeight(SolverPtr S,GraphTheorySolver_long G, int edgeLit);

  int reaches(SolverPtr S,GraphTheorySolver_long G,int from, int to);
  int reachesBackward(SolverPtr S,GraphTheorySolver_long G,int from, int to);
  int onPath(SolverPtr S,GraphTheorySolver_long G,int nodeOnPath,int from, int to);
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
  int minimumSpanningTree_lt(SolverPtr S,GraphTheorySolver_long G, Weight weight);
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
  int fsmCompositionAccepts(SolverPtr S, FSMTheorySolverPtr  fsmTheory,   int fsmGeneratorID,int fsmAcceptorID, int gen_startNode, int gen_acceptNode, int acceptor_startNode, int acceptor_acceptNode, int stringID);

  //model query
  //Returns true if the solver has a model (in which case it is safe to query the model), false otherwise
  bool hasModel(SolverPtr S);
  //For a given literal (not variable!), returns 0 for true, 1 for false, 2 for unassigned.
  int getModel_Literal(SolverPtr S,int lit);
  //Check if a literal is known to be a constant by the solver (eg, if it has an assignment at decision level 0).
  //For a given literal (not variable!), returns 0 for if the literal is a true constant, 1 if it is a false constant,
  //or 2 is it not assigned at level 0.
  int getConstantModel_Literal(SolverPtr S,int lit);
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
