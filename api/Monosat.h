/*
 * MonoSAT.h
 *
 *  Created on: May 21, 2015
 *      Author: sam
 */

#ifndef MONOSAT_H_
#define MONOSAT_H_
//Monosat library interface
#ifdef __cplusplus
#include <stdint.h>

extern "C"
{
typedef Monosat::SimpSolver *  SolverPtr;
typedef Monosat::GraphTheorySolver<int64_t> * GraphTheorySolver_long;
typedef Monosat::GraphTheorySolver<double>*  GraphTheorySolver_double;
typedef Monosat::BVTheorySolver<int64_t>* BVTheoryPtr;
typedef Monosat::FSMTheorySolver * FSMTheorySolverPtr;
#else
typedef void * SolverPtr;
typedef void*  BVTheoryPtr;
typedef void*  GraphTheorySolver_long;
typedef void*  GraphTheorySolver_double;
typedef void*  FSMTheorySolverPtr;
typedef int Var;
#endif


  int varToLit(int var, bool negated){
	  return toInt(mkLit(var,negated));
  }

  int litToVar(int lit){
	  return lit/2;
  }
  SolverPtr newSolver();
  SolverPtr newSolver_arg(char*argv);
#ifndef JNA
  //Java Native Access sometimes has problems with string arrays
  SolverPtr newSolver_args(int argc, char**argv);
#endif

  void deleteSolver (SolverPtr S);

  void readGNF(SolverPtr S, const char  * filename);

  bool solve(SolverPtr S);
  bool solveAssumptions(SolverPtr S,int * assumptions, int n_assumptions);
  //Solve under assumptions, and also minimize a set of BVs (in order of precedence)
  bool solveAssumptions_MinBVs(SolverPtr S,int * assumptions, int n_assumptions, int * minimize_bvs, int n_minimize_bvs);

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
  int solveAssumptionsLimited_MinBVs(SolverPtr S,int * assumptions, int n_assumptions, int * minimize_bvs, int n_minimize_bvs);



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

  int nVars(SolverPtr S);
  int nClauses(SolverPtr S);
  int nBitvectors(SolverPtr S,BVTheoryPtr bv);

  bool addClause(SolverPtr S,int * lits, int n_lits);
  bool addUnitClause(SolverPtr S,int lit);
  bool addBinaryClause(SolverPtr S,int lit1, int lit2);
  bool addTertiaryClause(SolverPtr S,int lit1, int lit2, int lit3);
  //theory interface for bitvectors
  BVTheoryPtr initBVTheory(SolverPtr S);
  int newBitvector_const(SolverPtr S, BVTheoryPtr bv, int bvWidth, int64_t constval);
  int newBitvector_anon(SolverPtr S, BVTheoryPtr bv, int bvWidth);
  int newBitvector(SolverPtr S, BVTheoryPtr bv, int * bits, int n_bits);
  int bv_width(SolverPtr S, BVTheoryPtr  bv,int bvID);
  int newBVComparison_const_lt(SolverPtr S, BVTheoryPtr bv, int bvID, int64_t weight);
  int newBVComparison_bv_lt(SolverPtr S, BVTheoryPtr bv, int bvID, int compareID);
  int newBVComparison_const_leq(SolverPtr S, BVTheoryPtr bv, int bvID, int64_t weight);
  int newBVComparison_bv_leq(SolverPtr S, BVTheoryPtr bv, int bvID, int compareID);
  int newBVComparison_const_gt(SolverPtr S, BVTheoryPtr bv, int bvID, int64_t weight);
  int newBVComparison_bv_gt(SolverPtr S, BVTheoryPtr bv, int bvID, int compareID);
  int newBVComparison_const_geq(SolverPtr S, BVTheoryPtr bv, int bvID, int64_t weight);
  int newBVComparison_bv_geq(SolverPtr S, BVTheoryPtr bv, int bvID, int compareID);


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
  void bv_min(SolverPtr S, BVTheoryPtr bv,  int* args,int n_args,int resultID);
  void bv_max(SolverPtr S, BVTheoryPtr bv,  int* args,int n_args,int resultID);
  void bv_popcount(SolverPtr S, BVTheoryPtr bv,  int* args,int n_args, int resultID);

  //simple at-most-one constraint: asserts that at most one of the set of variables (NOT LITERALS) may be true.
  //for small numbers of variables, consider using a direct CNF encoding instead
  void at_most_one(SolverPtr S, int * vars, int n_vars);

  //theory interface for graphs

  GraphTheorySolver_long newGraph(SolverPtr S);

  int newNode(SolverPtr S,GraphTheorySolver_long G);
  int newEdge(SolverPtr S, GraphTheorySolver_long G,int from,int  to,  int64_t weight);
  int newEdge_double(SolverPtr S, GraphTheorySolver_double G,int from,int  to,  double weight);
  int newEdge_bv(SolverPtr S, GraphTheorySolver_long G,int from,int  to, int bvID);
  int reaches(SolverPtr S,GraphTheorySolver_long G,int from, int to);
  int shortestPathUnweighted_lt_const(SolverPtr S,GraphTheorySolver_long G,int from, int to, int steps);
  int shortestPathUnweighted_leq_const(SolverPtr S,GraphTheorySolver_long G,int from, int to, int steps);
  int shortestPath_lt_const(SolverPtr S,GraphTheorySolver_long G,int from, int to, int64_t dist);
  int shortestPath_leq_const(SolverPtr S,GraphTheorySolver_long G,int from, int to, int64_t dist);
  int shortestPath_lt_bv(SolverPtr S,GraphTheorySolver_long G,int from, int to, int bvID);
  int shortestPath_leq_bv(SolverPtr S,GraphTheorySolver_long G,int from, int to, int bvID);
  int maximumFlow_geq(SolverPtr S,GraphTheorySolver_long G,int source, int sink, int64_t weight);
  int maximumFlow_gt(SolverPtr S,GraphTheorySolver_long G,int source, int sink, int64_t weight);
  int maximumFlow_geq_bv(SolverPtr S,GraphTheorySolver_long G,int source, int sink, int bvID);
  int maximumFlow_gt_bv(SolverPtr S,GraphTheorySolver_long G,int source, int sink, int bvID);
  int minimumSpanningTree_leq(SolverPtr S,GraphTheorySolver_long G, int64_t weight);
  int minimumSpanningTree_lt(SolverPtr S,GraphTheorySolver_long G,int source, int sink, int64_t weight);
  int acyclic_undirected(SolverPtr S,GraphTheorySolver_long G);
  int acyclic_directed(SolverPtr S,GraphTheorySolver_long G);

  //theory interface for finite state machines
  FSMTheorySolverPtr initFSMTheory(SolverPtr S);
  int newFSM(SolverPtr S,FSMTheorySolverPtr fsmTheory, int inputAlphabet, int outputAlphabet);
  int newState(SolverPtr S,FSMTheorySolverPtr fsmTheory, int fsmID);

  int newTransition(SolverPtr S,FSMTheorySolverPtr fsmTheory, int fsmID, int fromNode, int toNode,int inputLabel, int outputLabel);
  int newString(SolverPtr S,FSMTheorySolverPtr fsmTheory, int * str,int len);
  int fsmAcceptsString(SolverPtr S,FSMTheorySolverPtr fsmTheory, int fsmID, int startNode, int acceptNode,int stringID);

  //model query
  //For a given literal (not variable!), returns 0 for true, 1 for false, 2 for unassigned.
  int getModel_Literal(SolverPtr S,int lit);
  //Get an assignment to a bitvector in the model. The model may find a range of satisfying assignments to the bitvector;
  //If getMaximumValue is true, this function returns the maximum satisfying assignment to the bitvector in the model; else it returns the smallest.
  int64_t getModel_BV(SolverPtr S, BVTheoryPtr bv, int bvID, bool getMaximumValue);
  //graph queries:
  //maxflow_literal is the literal (not variable!) that is the atom for the maximum flow query
  int64_t getModel_MaxFlow(SolverPtr S,GraphTheorySolver_long G,int maxflow_literal);
  //maxflow_literal is the literal (not variable!) that is the atom for the maximum flow query
  int64_t getModel_EdgeFlow(SolverPtr S,GraphTheorySolver_long G,int maxflow_literal, int edgeLit);
  int64_t getModel_AcyclicEdgeFlow(SolverPtr S,GraphTheorySolver_long G,int maxflow_literal, int edgeLit);

  int64_t getModel_MinimumSpanningTreeWeight(SolverPtr S,GraphTheorySolver_long G,int spanning_tree_literal);
/*
  //Returns the number of nodes in the path length for this reachability or shortest path literal (1+number of edges)
  int getModel_PathLength(SolverPtr S,GraphTheorySolver_long G,int reach_or_shortest_path_lit);
  void getModel_Path(SolverPtr S,GraphTheorySolver_long G,int reach_or_shortest_path_lit, int * store_path);
  void getModel_PathByEdgeLit(SolverPtr S,GraphTheorySolver_long G,int reach_or_shortest_path_lit, int * store_path);
*/


#ifdef __cplusplus
}
#endif

#endif /* MONOSAT_H_ */
