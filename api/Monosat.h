/*
 * MonoSAT.h
 *
 *  Created on: May 21, 2015
 *      Author: sam
 */

#ifndef MONOSAT_H_
#define MONOSAT_H_
//Monosat library interface


extern "C"  //Tells the compile to use C-linkage for the next scope.
{



  void * newSolver(int argc, char**argv);


  void deleteSolver (Monosat::SimpSolver * S);

  bool solve(Monosat::SimpSolver * S);
  bool solveAssumption(Monosat::SimpSolver * S,int * assumptions, int n_assumptions);
  void backtrack(Monosat::SimpSolver * S);
  int newVar(Monosat::SimpSolver * S);

  bool addClause(Monosat::SimpSolver * S,int * assumptions, int n_assumptions);
  bool addUnitClause(Monosat::SimpSolver * S,int lit);
  bool addBinaryClause(Monosat::SimpSolver * S,int lit1, int lit2);
  bool addTertiaryClause(Monosat::SimpSolver * S,int lit1, int lit2, int lit3);
  //theory interface for bitvectors
  void * initBVTheory(Monosat::SimpSolver * S);
  int newBitvector_const(Monosat::SimpSolver * S, Monosat::BVTheorySolver<long> * bv, int bvWidth, long constval);
  int newBitvector(Monosat::SimpSolver * S, Monosat::BVTheorySolver<long> * bv, int * bits, int n_bits);
  int newBVComparison_const_lt(Monosat::SimpSolver * S, Monosat::BVTheorySolver<long> * bv, int bvID, long weight);
  int newBVComparison_bv_lt(Monosat::SimpSolver * S, Monosat::BVTheorySolver<long> * bv, int bvID, int compareID);
  int newBVComparison_const_leq(Monosat::SimpSolver * S, Monosat::BVTheorySolver<long> * bv, int bvID, long weight);
  int newBVComparison_bv_leq(Monosat::SimpSolver * S, Monosat::BVTheorySolver<long> * bv, int bvID, int compareID);
  int newBVComparison_const_gt(Monosat::SimpSolver * S, Monosat::BVTheorySolver<long> * bv, int bvID, long weight);
  int newBVComparison_bv_gt(Monosat::SimpSolver * S, Monosat::BVTheorySolver<long> * bv, int bvID, int compareID);
  int newBVComparison_const_geq(Monosat::SimpSolver * S, Monosat::BVTheorySolver<long> * bv, int bvID, long weight);
  int newBVComparison_bv_geq(Monosat::SimpSolver * S, Monosat::BVTheorySolver<long> * bv, int bvID, int compareID);

  void bv_addition( Monosat::SimpSolver * S, Monosat::BVTheorySolver<long> * bv, int bvID1, int bvID2, int resultID);

  //theory interface for graphs

  void * newGraph(Monosat::SimpSolver * S);

  int newNode(Monosat::SimpSolver * S,Monosat::GraphTheorySolver<long> *G);
  int newEdge(Monosat::SimpSolver * S, Monosat::GraphTheorySolver<long> *G,int from,int  to,  long weight);
  int newEdge_bv(Monosat::SimpSolver * S, Monosat::GraphTheorySolver<long> *G,int from,int  to, int bvID);
  int reaches(Monosat::SimpSolver * S,Monosat::GraphTheorySolver<long> *G,int from, int to);
  int shortestPathUnweighted_lt_const(Monosat::SimpSolver * S,Monosat::GraphTheorySolver<long> *G,int from, int to, int steps);
  int shortestPathUnweighted_leq_const(Monosat::SimpSolver * S,Monosat::GraphTheorySolver<long> *G,int from, int to, int steps);
  int shortestPath_lt_const(Monosat::SimpSolver * S,Monosat::GraphTheorySolver<long> *G,int from, int to, long dist);
  int shortestPath_leq_const(Monosat::SimpSolver * S,Monosat::GraphTheorySolver<long> *G,int from, int to, long dist);
  int shortestPath_lt_bv(Monosat::SimpSolver * S,Monosat::GraphTheorySolver<long> *G,int from, int to, int bvID);
  int shortestPath_leq_bv(Monosat::SimpSolver * S,Monosat::GraphTheorySolver<long> *G,int from, int to, int bvID);
  int maximumFlow_geq(Monosat::SimpSolver * S,Monosat::GraphTheorySolver<long> *G,int source, int sink, long weight);
  int maximumFlow_gt(Monosat::SimpSolver * S,Monosat::GraphTheorySolver<long> *G,int source, int sink, long weight);
  int maximumFlow_geq_bv(Monosat::SimpSolver * S,Monosat::GraphTheorySolver<long> *G,int source, int sink, int bvID);
  int maximumFlow_gt_bv(Monosat::SimpSolver * S,Monosat::GraphTheorySolver<long> *G,int source, int sink, int bvID);
  int minimumSpanningTree_leq(Monosat::SimpSolver * S,Monosat::GraphTheorySolver<long> *G, long weight);
  int minimumSpanningTree_lt(Monosat::SimpSolver * S,Monosat::GraphTheorySolver<long> *G,int source, int sink, long weight);
  int acyclic_undirected(Monosat::SimpSolver * S,Monosat::GraphTheorySolver<long> *G);
  int acyclic_directed(Monosat::SimpSolver * S,Monosat::GraphTheorySolver<long> *G);

  //model query
  //For a given literal (not variable!), returns 0 for true, 1 for false, 2 for unassigned.
  int getModel_Literal(Monosat::SimpSolver * S,int lit);
  long getModel_BV(Monosat::SimpSolver * S, Monosat::BVTheorySolver<long> * bv, int bvID);
  //graph queries:
  //maxflow_literal is the literal (not variable!) that is the atom for the maximum flow query
  long getModel_MaxFlow(Monosat::SimpSolver * S,Monosat::GraphTheorySolver<long> *G,int maxflow_literal);
  //maxflow_literal is the literal (not variable!) that is the atom for the maximum flow query
  long getModel_EdgeFlow(Monosat::SimpSolver * S,Monosat::GraphTheorySolver<long> *G,int maxflow_literal, int edgeLit);
  long getModel_MinimumSpanningTreeWeight(Monosat::SimpSolver * S,Monosat::GraphTheorySolver<long> *G,int spanning_tree_literal);

}
#endif /* MONOSAT_H_ */
