/*
The MIT License (MIT)

Copyright (c) 2018, Sam Bayless

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
*/

package monosat;

import java.nio.IntBuffer;
import java.util.Vector;

/**
 * Low-level JNI mapping of the C interface to MonoSAT. This class is not intended for end users,
 * and its interface may change without warning if the native library changes. See instead
 * monosat.Solver.
 */
final class MonosatJNI { // package level access specifier
  static {
    System.loadLibrary("monosat");
  }

  // Get the version string from the native library.
  public static native String getVersion();

  // literal/variable transformations
  public static native int varToLit(int var, boolean negated);

  public static native int litToVar(int lit);

  // create/destroy solver pointers
  public static native long newSolver();

  public static native long newSolver(String args);

  public static native void deleteSolver(long solverPtr);

  public static native boolean ok(long solverPtr);

  // If set, dump constraints to this file (as they are asserted in the solver)
  public static native void setOutputFile(long solverPtr, String filename);

  // load constraints from a file, and execute any embedded solve/optimize calls
  public static native void readGNF(long solverPtr, String filename);

  // load constraints from a file, but ignore any embedded solve/optimize calls
  public static native void loadGNF(long solverPtr, String filename);

  //flush constraints to file (if an output file is open)
  public static native void flushFile (long solverPtr);
  //stop writing constraints to file, and close the file (if any)
  //this will be called automatically if the solver is deleted
  public static native void closeFile (long solverPtr);



    /**
     * Adds a name to a literal. The name must consist of printable ascii characters,
     * must not include `~`, and must be unique.
     * If the variable previously had a name, it will have multiple names after this call.
     * When a literal is assigned a name, its opposite polarity is also assigned that name,
     * prefixed by a '~'.
     */
    public static native void addLiteralName(long solverPtr, int literal, String name);

    // True if the specified variable has the given name
    public static native boolean literalHasName(long solverPtr, int literal, String name);

    /*
     * Retrieve a literal's name. If the literal has more than one name, then name index
     * selects which name to return (in the order they were assigned).
     * Otherwise, it should be set to 0.
     * If the literal has no names, or if nameIndex is out of range, this returns
     * the empty string.
     */
    public static native String getLiteralName(long solverPtr, int literal, int nameIndex);

    /**
     * Return the number of names associated with a literal.
     */
    public static native int literalNameCount(long solverPtr, int literal);

    public static native boolean hasLiteralWithName(long solverPtr, String name);

    //get the number of named literals in the solver
    public static native int nNamedLiterals(long solverPtr);
    //Get the nth named literal
    public static native int getNamedLiteralN(long solverPtr, int n);

    public static native int getLiteral(long solverPtr, String name);

    /**
     * Adds a name to a literal. The name must consist of printable ascii characters,
     * must not include `~`, and must be unique.
     * If the variable previously had a name, it will have multiple names after this call.
     */
    @Deprecated
    public static native void addVariableName(long solverPtr, int literal, String name);

    // True if the specified variable has the given name
    @Deprecated
    public static native boolean variableHasName(long solverPtr, int variable, String name);


    /*
     * Get a variable's name. If the variable has more than one name, then name index
     * selects which name to return (in the order they were assigned).
     * Otherwise, it should be set to 0.
     * If the variable has no names, or if nameIndex is out of range, this returns
     * the empty string.
     */
    @Deprecated
    public static native String getVariableName(long solverPtr, int variable, int nameIndex);


    //Return the number of names associated with a variable.
    @Deprecated
    public static native int variableNameCount(long solverPtr, int variable);
    @Deprecated
    public static native boolean hasVariableWithName(long solverPtr, String name);

    //get the number of named variables in the solver
    @Deprecated
    public static native int nNamedVariables(long solverPtr);
    //Get the nth named variable
    @Deprecated
    public static native int getNamedVariableN(long solverPtr, int n);

    @Deprecated
    public static native int getVariable(long solverPtr, String name);

  // basic solver functions
  public static native boolean solve(long solverPtr);

  public static native boolean solveAssumptions(
      long solverPtr, IntBuffer assumps, int n_assumptions);

  // Sets the (approximate) time limit in seconds before returning l_Undef from solveLimited;
  // ignored by solve(). Set to <0 to disable time limit.
  public static native void setTimeLimit(long solverPtr, int seconds);


  // Sets the maximum number of (additional) conflicts allowed in the solver before returning
  // l_Undef from solveLimited; ignored by solve(). Set to <0 to disable conflict limit.
  public static native void setConflictLimit(long solverPtr, int num_conflicts);

  //Number of conflicts in the solver
  public static native long nConflicts(long solverPtr);

  // Sets the maximum number of (additional) propagation rounds allowed in the solver before
  // returning l_Undef from
  // solveLimited; ignored by solve(). Set to <0 to disable propagation limit.
  public static native void setPropagationLimit(long solverPtr, int num_propagations);

  //number of propagation rounds in the solver
  public static native long nPropagations(long solverPtr);

  // Returns 0 for satisfiable, 1 for proved unsatisfiable, 2 for failed to find a solution (within
  // any resource limits that have been set)
  public static native int solveLimited(long solverPtr);

  // Returns 0 for satisfiable, 1 for proved unsatisfiable, 2 for failed to find a solution (within
  // any resource limits that have been set)
  public static native int solveAssumptionsLimited(
      long solverPtr, IntBuffer assumptions, int n_assumptions);

  public static native boolean lastSolutionWasOptimal(long solverPtr);

  // If the last solution was unsat, then this get the 'conflict clause' produced by the solver (a
  // subset of the assumptions which are sufficient to cause the instance to be UNSAT).
  // Fills the given pointer with the first max_store_size literals of the conflict clause, and
  // returns the number of literals in the conflict clause. Set store_clause to null and
  // max_store_size to 0 to find the size of the conflict clause
  // Returns -1 if the solver has no conflict clause from the most recent solve() call (because that
  // call was not UNSAT)
  public static native int getConflictClause(long solverPtr, IntBuffer store_clause, int length);

  // Reduce the given set of mutually UNSAT assumptions into a (locally) minimal-sized set of
  // assumptions that are still UNSAT.
  // If the assumptions are not UNSAT, this method does returns -1.
  // Else, it returns the number of literals in the unsat core, which it will also store in the
  // first number_literals
  // entries in unsat_assumptions
  // Note: This method will make repeated (and potentially expensive) calls to the SAT solver to
  // attempt to remove literals from the
  // set of assumptions.
  // Note also that if any of the setXXXXXLimits() are applied to the solver, this may not produce a
  // locally minimal unsat core.
  public static native int minimizeUnsatCore(
      long solverPtr, IntBuffer unsat_assumptions, int length);

  // After UNSAT solve calls with assumptions, the solver will find a 'conflict clause' consisting
  // of a subset of the assumptions
  // which are sufficient to make the solver UNSAT (see getConflictClause).
  // Normally, the conflict clause is produced as a side effect of proving the query unsat, with the
  // solver only removing literals
  // from the conflict clause on a best effort basis.
  // This method will make repeated (and potentially expensive) calls to the SAT solver to attempt
  // to remove further literals from
  // the conflict clause.
  // Afterward, the conflict clause can be obtained using getConflictClause().
  // NOTE: this function may be expensive, is not required to get a conflict clause;
  // getConflictClause() can be used after any unsat call with assumptions, even
  // without calling minimizeConflictClause().
  // Note also that if any of the setXXXXXLimits() are applied to the solver, this may not produce a
  // locally minimal conflict clause.
  public static native void minimizeConflictClause(long solverPtr);

  public static native void backtrack(long solverPtr);

  public static native int newVar(long solverPtr);

  public static native int newNamedVar(long solverPtr, String name);

  public static native void releaseLiteral(long solverPtr, int literal);

  public static native void setDecisionVar(long solverPtr, int var, boolean decidable);

  public static native boolean isDecisionVar(long solverPtr, int var);

  // Static, lexicographic heuristic. Larger values are higher priority (decided first); default
  // priority is 0
  public static native void setDecisionPriority(long solverPtr, int var, int priority);

  public static native int getDecisionPriority(long solverPtr, int var);

  // Which polarity the decision heuristic should use for a variable (by default).
  public static native void setDecisionPolarity(long solverPtr, int v, boolean b);

  public static native boolean getDecisionPolarity(long solverPtr, int v);

  // The solver will (sometimes) instantiate an arbitrary true literal for use as a constant.
  // Call this method to a) force that literal to be instantiate, and b) get that literal.
  public static native int true_lit(long solverPtr);

  // Prevents this literal from being simplified by the preprocessor
  public static native boolean disallowLiteralSimplification(long solverPtr, int var);

  // permanently disable SAT-based pre-processing in this solver
  public static native void disablePreprocessing(long solverPtr);

  public static native int nVars(long solverPtr);

  public static native int nClauses(long solverPtr);

  public static native int nLearnedClauses(long solverPtr);

  public static native int nBitvectors(long solverPtr, long bvPtr);

  public static native boolean addClause(long solverPtr, IntBuffer lits, int n_lits);

  public static native boolean addUnitClause(long solverPtr, int lit);

  public static native boolean addBinaryClause(long solverPtr, int lit1, int lit2);

  public static native boolean addTertiaryClause(long solverPtr, int lit1, int lit2, int lit3);

  // remove any optimization objectives from the solver
  public static native void clearOptimizationObjectives(long solverPtr);

  public static native void maximizeBV(long solverPtr, long bvPtr, int bvID);

  public static native void minimizeBV(long solverPtr, long bvPtr, int bvID);

  public static native void maximizeLits(long solverPtr, IntBuffer lits, int n_lits);

  public static native void minimizeLits(long solverPtr, IntBuffer lits, int n_lits);

  public static native void maximizeWeightedLits(
      long solverPtr, IntBuffer lits, IntBuffer weights, int n_lits);

  public static native void minimizeWeightedLits(
      long solverPtr, IntBuffer lits, IntBuffer weights, int n_lits);

  // theory interface for bitvectors
  public static native long initBVTheory(long solverPtr);

  public static native int newBitvector_const(
      long solverPtr, long bvPtr, int bvWidth, long constVal);

  public static native int newBitvector_anon(long solverPtr, long bvPtr, int bvWidth);

  public static native int newBitvector_lazy(long solverPtr, long bvPtr, IntBuffer bits, int n_bits);

  public static native int newBitvector(long solverPtr, long bvPtr, IntBuffer bits, int n_bits);

  public static native void setBitvectorName(long solverPtr, long bvPtr, int bvID, String name);

  public static native boolean bitvectorHasName(long solverPtr, long bvPtr, int bvID);
  public static native boolean hasBitvectorWithName(long solverPtr, long bvPtr, String name);

  public static native String getBitvectorName(long solverPtr, long bvPtr, int bvID, int nameIndex);

  public static native int getBitvectorNameCount(long solverPtr, long bvPtr, int bvID);

  public static native int getBitvectorWidth(long solverPtr, long bvPtr, int bvID);
  // Number of defined literals in the BV, which may be 0 or the bitvector width
  public static native int nBitvectorBits(long solverPtr, long bvPtr, int bvID);

  public static native int getBitvector(long solverPtr, long bvPtr, String name);

  //get the number of named bitvector in the solver
  public static native int nNamedBitvectors(long solverPtr,long bvPtr);
  //Get the nth named bitvector
  public static native int getNamedBitvectorN(long solverPtr,long bvPtr, int n);

  public static native int getBitvectorBit(long solverPtr, long bvPtr, int bvID, int bit);

  public static native int newBVComparison_const_lt(
      long solverPtr, long bvPtr, int bvID, long constVal);

  public static native int newBVComparison_bv_lt(
      long solverPtr, long bvPtr, int bvID, int compareID);

  public static native int newBVComparison_const_leq(
      long solverPtr, long bvPtr, int bvID, long constVal);

  public static native int newBVComparison_bv_leq(
      long solverPtr, long bvPtr, int bvID, int compareID);

  public static native int newBVComparison_const_gt(
      long solverPtr, long bvPtr, int bvID, long constVal);

  public static native int newBVComparison_bv_gt(
      long solverPtr, long bvPtr, int bvID, int compareID);

  public static native int newBVComparison_const_geq(
      long solverPtr, long bvPtr, int bvID, long constVal);

  public static native int newBVComparison_bv_geq(
      long solverPtr, long bvPtr, int bvID, int compareID);

  public static native int newBVComparison_const_eq(
      long solverPtr, long bvPtr, int bvID, long constVal);

  public static native int newBVComparison_bv_eq(
      long solverPtr, long bvPtr, int bvID, int compareID);

  public static native int newBVComparison_const_neq(
      long solverPtr, long bvPtr, int bvID, long constVal);

  public static native int newBVComparison_bv_neq(
      long solverPtr, long bvPtr, int bvID, int compareID);

  // Convert the specified bitvector, as well as any other bitvectors in its cone of influence, into
  // pure CNF
  public static native void bv_bitblast(long solverPtr, long bvPtr, int bvID);

  public static native void bv_concat(long solverPtr, long bvPtr, int aID, int bID, int resultID);

  public static native void bv_slice(
      long solverPtr, long bvPtr, int aID, int lower, int upper, int resultID);

  public static native void bv_not(long solverPtr, long bvPtr, int bvaID, int bvResultID);

  public static native void bv_and(
      long solverPtr, long bvPtr, int bvaID, int bvbID, int bvResultID);

  public static native void bv_nand(
      long solverPtr, long bvPtr, int bvaID, int bvbID, int bvResultID);

  public static native void bv_or(long solverPtr, long bvPtr, int bvaID, int bvbID, int bvResultID);

  public static native void bv_nor(
      long solverPtr, long bvPtr, int bvaID, int bvbID, int bvResultID);

  public static native void bv_xor(
      long solverPtr, long bvPtr, int bvaID, int bvbID, int bvResultID);

  public static native void bv_xnor(
      long solverPtr, long bvPtr, int bvaID, int bvbID, int bvResultID);

  public static native void bv_ite(
      long solverPtr, long bvPtr, int condition_lit, int bvThenID, int bvElseID, int bvResultID);

  public static native void bv_addition(
      long solverPtr, long bvPtr, int bvID1, int bvID2, int resultID);

  public static native void bv_subtraction(
      long solverPtr, long bvPtr, int bvID1, int bvID2, int resultID);

  public static native void bv_multiply(
      long solverPtr, long bvPtr, int bvID1, int bvID2, int resultID);

  public static native void bv_divide(
      long solverPtr, long bvPtr, int bvID1, int bvID2, int resultID);

  public static native void bv_min(
      long solverPtr, long bvPtr, IntBuffer args, int n_args, int resultID);

  public static native void bv_max(
      long solverPtr, long bvPtr, IntBuffer args, int n_args, int resultID);

  public static native void bv_popcount(
      long solverPtr, long bvPtr, IntBuffer args, int n_args, int resultID);

  public static native void bv_unary(
      long solverPtr, long bvPtr, IntBuffer args, int n_args, int resultID);

  // simple at-most-one constraint: asserts that at most one of the set of lit
  // may be true.
  // for small numbers of variables, consider using a direct CNF encoding instead
  public static native void at_most_one_lit(long solverPtr, IntBuffer lits, int n_lits);

  public static native void assertPB_lt(
      long solverPtr, int rhs, int n_args, IntBuffer literals, IntBuffer coefficients);

  public static native void assertPB_leq(
      long solverPtr, int rhs, int n_args, IntBuffer literals, IntBuffer coefficients);

  public static native void assertPB_eq(
      long solverPtr, int rhs, int n_args, IntBuffer literals, IntBuffer coefficients);

  public static native void assertPB_geq(
      long solverPtr, int rhs, int n_args, IntBuffer literals, IntBuffer coefficients);

  public static native void assertPB_gt(
      long solverPtr, int rhs, int n_args, IntBuffer literals, IntBuffer coefficients);

  // Convert any pb constraints in the solver into cnf (will be called automatically before solve())
  public static native void flushPB(long solverPtr);

  // theory interface for graphs
  public static native long newGraph(long solverPtr);

  public static native long newGraph_Named(long solverPtr, String name, int bitwidth);

  public static native String getGraphName(long solverPtr, long graphPtr);

  public static native int getGraphWidth(long solverPtr, long graphPtr);
  /**
   * Get a pointer to a graph with this name, that is already defined in the solver (or return 0, if
   * there is no such graph)
   *
   * @param solverPtr
   * @param name
   * @return a pointer to a graph with this name, that is already defined in the solver (or return
   *     0, if there is no such * graph)
   */
  public static native long getGraph(long solverPtr, String name);

  public static native int newNode(long solverPtr, long graphPtr);

  public static native int newNode_Named(long solverPtr, long graphPtr, String name);

  public static native boolean hasNamedNode(long solverPtr, long graphPtr, String name);

  public static native String getNodeName(long solverPtr, long graphPtr, int node);

  public static native int newEdge(
      long solverPtr, long graphPtr, int from, int to, long constweight);

  public static native int newEdge_bv(long solverPtr, long graphPtr, int from, int to, int bvID);

  public static native int nNodes(long solverPtr, long graphPtr);

  public static native int nEdges(long solverPtr, long graphPtr);

  public static native int getEdgeLiteralN(long solverPtr, long graphPtr, int n);

  public static native int getEdge_to(long solverPtr, long graphPtr, int edgeLit);

  public static native int getEdge_from(long solverPtr, long graphPtr, int edgeLit);

  public static native long getEdge_weight_const(long solverPtr, long graphPtr, int edgeLit);

  public static native int getEdge_weight_bv(long solverPtr, long graphPtr, int edgeLit);

  public static native boolean edgeHasBVWeight(long solverPtr, long graphPtr, int edgeLit);
  
  public static native int reaches(long solverPtr, long graphPtr, int from, int to);

  public static native int reachesBackward(long solverPtr, long graphPtr, int from, int to);

  public static native int onPath(long solverPtr, long graphPtr, int nodeOnPath, int from, int to);

  public static native int shortestPathUnweighted_lt_const(
      long solverPtr, long graphPtr, int from, int to, int steps);

  public static native int shortestPathUnweighted_leq_const(
      long solverPtr, long graphPtr, int from, int to, int steps);

  public static native int shortestPath_lt_const(
      long solverPtr, long graphPtr, int from, int to, long dist);

  public static native int shortestPath_leq_const(
      long solverPtr, long graphPtr, int from, int to, long dist);

  public static native int shortestPath_lt_bv(
      long solverPtr, long graphPtr, int from, int to, int bvID);

  public static native int shortestPath_leq_bv(
      long solverPtr, long graphPtr, int from, int to, int bvID);

  public static native int maximumFlow_geq(
      long solverPtr, long graphPtr, int source, int sink, long constweight);

  public static native int maximumFlow_gt(
      long solverPtr, long graphPtr, int source, int sink, long constweight);

  public static native int maximumFlow_geq_bv(
      long solverPtr, long graphPtr, int source, int sink, int bvID);

  public static native int maximumFlow_gt_bv(
      long solverPtr, long graphPtr, int source, int sink, int bvID);

  public static native int minimumSpanningTree_leq(long solverPtr, long graphPtr, long constweight);

  public static native int minimumSpanningTree_lt(
      long solverPtr, long graphPtr, int source, int sink, long constweight);

  public static native int acyclic_undirected(long solverPtr, long graphPtr);

  public static native int acyclic_directed(long solverPtr, long graphPtr);

  public static native void newEdgeSet(
      long solverPtr, long graphPtr, IntBuffer edges, int n_edges, boolean enforceEdgeAssignment);

  // this enables a heuristic on this graph, from the RUC paper, which sets assigned edges to zero
  // long, to encourage edge-reuse in solutions
  public static native void graph_setAssignEdgesToWeight(
      long solverPtr, long graphPtr, long weight);
  // flow routing interface

  public static native long createFlowRouting(
      long solverPtr, long graphPtr, int sourceNode, int destNode, int maxflowLit);

  public static native void addRoutingNet(
      long solverPtr,
      long graphPtr,
      long routerPtr,
      int disabledEdge,
      int n_members,
      IntBuffer edge_lits,
      IntBuffer reach_lits);

  // model query
  // Returns true if the solver has a model (in which case it is safe to query the model), false
  // otherwise
  public static native boolean hasModel(long solverPtr);

  // For a given literal (not variable!), returns 0 for true, 1 for false, 2 for unassigned.
  public static native int getModel_Literal(long solverPtr, int lit);

  // Check if a literal is known to be a constant by the solver (eg, if it has an assignment at
  // decision level 0).
  // For a given literal (not variable!), returns 0 for if the literal is a true constant, 1 if it
  // is a false constant,
  // or 2 is it not assigned at level 0.
  public static native int getConstantModel_Literal(long solverPtr, int lit);

  // Get an assignment to a bitvector in the model. The model may find a range of satisfying
  // assignments to the bitvector;
  // If getMaximumValue is true, this function returns the maximum satisfying assignment to the
  // bitvector in the model; else it returns the smallest.
  public static native long getModel_BV(
      long solverPtr, long bvPtr, int bvID, boolean getMaximumValue);

  // graph queries:
  // maxflow_literal is the literal (not variable!) that is the atom for the maximum flow query
  public static native long getModel_MaxFlow(long solverPtr, long graphPtr, int maxflow_literal);

  // maxflow_literal is the literal (not variable!) that is the atom for the maximum flow query
  public static native long getModel_EdgeFlow(
      long solverPtr, long graphPtr, int maxflow_literal, int edgeLit);

  public static native long getModel_AcyclicEdgeFlow(
      long solverPtr, long graphPtr, int maxflow_literal, int edgeLit);

  public static native long getModel_MinimumSpanningTreeWeight(
      long solverPtr, long graphPtr, int spanning_tree_literal);

  public static native int getModel_Path_Nodes_Length(
      long solverPtr, long graphPtr, int reach_or_distance_literal);

  public static native int getModel_Path_Nodes(
      long solverPtr,
      long graphPtr,
      int reach_or_distance_literal,
      int store_length,
      IntBuffer store);

  public static native int getModel_Path_EdgeLits_Length(
      long solverPtr, long graphPtr, int reach_or_distance_literal);

  public static native int getModel_Path_EdgeLits(
      long solverPtr,
      long graphPtr,
      int reach_or_distance_literal,
      int store_length,
      IntBuffer store);

  // Logic building methods
  // Note: Many of these methods are capitalized, defying Java naming conventions.
  // This is because the lower case versions (and, or...) are keywords in C++, and so cannot be used
  // in the JNI bridge.
  public static native int And_(long solverPtr, int lit_a, int lit_b, int lit_out);

  public static native int Ands_(long solverPtr, IntBuffer lits, int n_lits, int lit_out);

  public static native void AssertImpliesAnd_(
      long solverPtr, int implies, IntBuffer lits, int n_lits, int lit_out);

  public static native void AssertImpliesAnd(long solverPtr, int implies, IntBuffer lits, int n_lits);

  public static native int Ands(long solverPtr, IntBuffer lits, int n_lits);

  public static native int And(long solverPtr, int lit_a, int lit_b);

  public static native int Or_(long solverPtr, int lit_a, int lit_b, int lit_out);

  public static native int Ors_(long solverPtr, IntBuffer lits, int n_lits, int lit_out);

  // If this gate is true, then all of lits must be true.
  // But if this gate is false, lits may be true or false.
  public static native int ImpliesAnd(long solverPtr, IntBuffer lits, int n_lits, int lit_out);

  // If this gate is true, then at least one of lits must be true.
  // But if this gate is false, lits may be true or false.
  public static native int ImpliesOr(long solverPtr, IntBuffer lits, int n_lits);

  public static native int ImpliesOr_(long solverPtr, IntBuffer lits, int n_lits, int lit_out);

  // This is an OR condition that holds only if implies is true
  public static native void AssertImpliesOr_(
      long solverPtr, int implies, IntBuffer lits, int n_lits, int lit_out);

    public static native void AssertImpliesOr(
            long solverPtr, int implies, IntBuffer lits, int n_lits);

  public static native int Ors(long solverPtr, IntBuffer lits, int n_lits);

  public static native int Or(long solverPtr, int lit_a, int lit_b);

  public static native int Nors(long solverPtr, IntBuffer lits, int n_lits);

  public static native int Nor(long solverPtr, int lit_a, int lit_b);

  public static native int Nands(long solverPtr, IntBuffer lits, int n_lits);

  public static native int Nand(long solverPtr, int lit_a, int lit_b);

  public static native int Xors(long solverPtr, IntBuffer lits, int n_lits);

  public static native int Xor(long solverPtr, int lit_a, int lit_b);

  public static native int Xnors(long solverPtr, IntBuffer lits, int n_lits);

  public static native int Xnor(long solverPtr, int lit_a, int lit_b);

  public static native int Implies(long solverPtr, int lit_a, int lit_b);

  public static native int Implies_(long solverPtr, int lit_a, int lit_b, int lit_out);

  public static native int Ite(long solverPtr, int lit_cond, int lit_thn, int lit_els);

  public static native int Ite_(
      long solverPtr, int lit_cond, int lit_thn, int lit_els, int lit_out);

  public static native int Add(
      long solverPtr, IntBuffer lits_a, IntBuffer lits_b, int n_lits, IntBuffer lits_out);

  public static native int Add_(
      long solverPtr,
      IntBuffer lits_a,
      IntBuffer lits_b,
      int n_lits,
      IntBuffer lits_out,
      int carry_lit);

  public static native int Subtract(
      long solverPtr, IntBuffer lits_a, IntBuffer lits_b, int n_lits, IntBuffer lits_out);

  public static native int Subtract_(
      long solverPtr,
      IntBuffer lits_a,
      IntBuffer lits_b,
      int n_lits,
      IntBuffer lits_out,
      int borrow_lit);

  // perform two's complement negation (long solverPtr,invert bits and add 1)
  // IntBuffer lits_out must have enough room to hold n_lits
  public static native void Negate(long solverPtr, IntBuffer lits, int n_lits, IntBuffer lits_out);

  // perform two's complement negation (long solverPtr,invert bits and add 1)
  public static native void Negate_(long solverPtr, IntBuffer lits, int n_lits, IntBuffer lits_out);

  public static native void Assert(long solverPtr, int lit);

  public static native void AssertOrTertiary(long solverPtr, int lit_a, int lit_b, int lit_c);

  public static native void AssertOrs(long solverPtr, IntBuffer lits, int n_lits);

  public static native void AssertOr(long solverPtr, int lit_a, int lit_b);

  public static native void AssertNands(long solverPtr, IntBuffer lits, int n_lits);

  public static native void AssertNand(long solverPtr, int lit_a, int lit_b);

  public static native void AssertAnds(long solverPtr, IntBuffer lits, int n_lits);

  public static native void AssertAnd(long solverPtr, int lit_a, int lit_b);

  public static native void AssertNors(long solverPtr, IntBuffer lits, int n_lits);

  public static native void AssertNor(long solverPtr, int lit_a, int lit_b);

  public static native void AssertXor(long solverPtr, int lit_a, int lit_b);

  public static native void AssertXors(long solverPtr, IntBuffer lits, int n_lits);

  public static native void AssertXnors(long solverPtr, IntBuffer lits, int n_lits);

  public static native void AssertXnor(long solverPtr, int lit_a, int lit_b);

  public static native void AssertImplies(long solverPtr, int lit_a, int lit_b);

  public static native void AssertEqual(long solverPtr, int lit_a, int lit_b);

  public static native void AssertAllSame(long solverPtr, IntBuffer lits, int n_lits);

  public static native int Equals(long solverPtr, IntBuffer A_lits, IntBuffer B_lits, int n_lits);

  public static native int LEQ(long solverPtr, IntBuffer A_lits, IntBuffer B_lits, int n_lits);

  public static native int LT(long solverPtr, IntBuffer A_lits, IntBuffer B_lits, int n_lits);

  public static native void AssertEquals(
      long solverPtr, IntBuffer A_lits, IntBuffer B_lits, int n_lits);

  public static native void AssertLEQ(
      long solverPtr, IntBuffer A_lits, IntBuffer B_lits, int n_lits);

  public static native void AssertLT(
      long solverPtr, IntBuffer A_lits, IntBuffer B_lits, int n_lits);

  // uses n^2 binary clauses to create a simple at-most-one constraint.
  // if you have more than 20 or so literals, strongly consider using a pseudo-Boolean constraint
  // solver instead
  public static native void AssertAMO(long solverPtr, IntBuffer lits, int n_lits);

  // uses n^2 binary clauses to create a simple exactly-one-constraint.
  // if you have more than 20 or so literals, strongly consider using a pseudo-Boolean constraint
  // solver instead
  public static native void AssertExactlyOne(long solverPtr, IntBuffer lits, int n_lits);

  /** Simple test driver, useful for debugging native library linking issues. */
  @SuppressWarnings("unchecked")
  public static void main(String[] args) {
    System.out.println("Loading native libraries:");
    try {
      java.lang.reflect.Field LIBRARIES = ClassLoader.class.getDeclaredField("loadedLibraryNames");
      LIBRARIES.setAccessible(true);
      Vector<String> libraries = (Vector<String>) LIBRARIES.get(ClassLoader.getSystemClassLoader());
      for (String s : libraries) System.out.println(s);

      String version = MonosatJNI.getVersion(); // invoke the native method
      System.out.println("Loaded MonoSAT: " + version);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  /**
   * Check whether a string is a valid monosat ID. If not, throw an exception; if it is, return that
   * string.
   *
   * @param name The string to check
   * @return str, unaltered.
   * @throws IllegalArgumentException If the string contains illegal characters
   */
  protected static String validID(String name) {
    if (name.chars()
        .allMatch(c -> (c < 128 && !Character.isISOControl(c) && !Character.isWhitespace(c)))) {
      return name;
    } else {
      throw new IllegalArgumentException(
          "IDs are restricted to printable ASCII characeters, and "
              + "may not include whitespace or newlines: \""
              + name
              + "\" is not a valid name.");
    }
  }
}
