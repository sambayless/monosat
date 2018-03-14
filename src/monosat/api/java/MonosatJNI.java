import java.util.ArrayList;
/**
 * Low-level JNI mapping of the C interface to MonoSAT.
 */
public class MonosatJNI {
    static {
        System.loadLibrary("monosat");
    }

    public native String getVersion();

    //literal/variable transformations
    public native int varToLit(int var, boolean negated);
    public native int litToVar(int lit);

    //create/destroy solver pointers
    public native long newSolver();
    public native long newSolver(String args);
    public native long newSolver(ArrayList<String> args);
    public native void deleteSolver(long solverPtr);

    //If set, dump constraints to this file (as they are asserted in the solver)
    public native void setOutputFile(long solverPtr, String filename);
    //load constraints from a file
    public native void readGNF(long solverPtr, String filename);

    //basic solver functions
    public native boolean solve(long solverPtr);
    public native boolean solveAssumptions(long solverPtr, int [] assumps);

    //Sets the (approximate) time limit in seconds before returning l_Undef from solveLimited; ignored by solve(). Set to <0 to disable time limit.
    public native void setTimeLimit(long solverPtr,int seconds);
    //Sets the (approximate) memory limit in megabytes before returning l_Undef from solveLimited; ignored by solve(). Set to <0 to disable memory limit.
    public native void setMemoryLimit(long solverPtr,int mb);
    //Sets the maximum number of (additional) conflicts allowed in the solver before returning l_Undef from solveLimited; ignored by solve(). Set to <0 to disable conflict limit.
    public native void setConflictLimit(long solverPtr,int num_conflicts);
    //Sets the maximum number of (additional) propagations allowed in the solver before returning l_Undef from solveLimited; ignored by solve(). Set to <0 to disable propagation limit.
    public native void setPropagationLimit(long solverPtr,int num_propagations);

    //Returns 0 for satisfiable, 1 for proved unsatisfiable, 2 for failed to find a solution (within any resource limits that have been set)
    public native int solveLimited(long solverPtr);
    //Returns 0 for satisfiable, 1 for proved unsatisfiable, 2 for failed to find a solution (within any resource limits that have been set)
    public native int solveAssumptionsLimited(long solverPtr,int [] assumptions);


    public native boolean lastSolutionWasOptimal(long solverPtr);

    //If the last solution was unsat, then this get the 'conflict clause' produced by the solver (a subset of the assumptions which are sufficient to cause the instance to be UNSAT).
    //Fills the given pointer with the first max_store_size literals of the conflict clause, and returns the number of literals in the conflict clause. Set store_clause to null and max_store_size to 0 to find the size of the conflict clause
    //Returns -1 if the solver has no conflict clause from the most recent solve() call (because that call was not UNSAT)
    public native int getConflictClause(long solverPtr, int [] store_clause, int max_store_size);

    public native void backtrack(long solverPtr);
    public native int newVar(long solverPtr);
    public native void setDecisionVar(long solverPtr,int var,boolean decidable);
    public native boolean isDecisionVar(long solverPtr,int var);

    //Static, lexicographic heuristic. Larger values are higher priority (decided first); default priority is 0
    public native void setDecisionPriority(long solverPtr,int var, int priority);
    public native int getDecisionPriority(long solverPtr,int var);
    // Which polarity the decision heuristic should use for a variable (by default).
    public native void setDecisionPolarity(long solverPtr,int v, boolean b);
    public native boolean getDecisionPolarity(long solverPtr,int v);

    //The solver will (sometimes) instantiate an arbitrary true literal for use as a constant.
    //Call this method to a) force that literal to be instantiate, and b) get that literal.
    public native int true_lit(long solverPtr);

    //Prevents this literal from being simplified by the preprocessor
    public native boolean disallowLiteralSimplification(long solverPtr, int lit);

    //permanently disable SAT-based preprocessing in this solver
    public native void disablePreprocessing(long solverPtr);
    public native int nVars(long solverPtr);
    public native int nClauses(long solverPtr);
    public native int nBitvectors(long solverPtr,long bvPtr);

    public native boolean addClause(long solverPtr,int [] lits, int n_lits);
    public native boolean addUnitClause(long solverPtr,int lit);
    public native boolean addBinaryClause(long solverPtr,int lit1, int lit2);
    public native boolean addTertiaryClause(long solverPtr,int lit1, int lit2, int lit3);

    //remove any optimization objectives from the solver
    public native void clearOptimizationObjectives(long solverPtr);

    public native void maximizeBV(long solverPtr,  long bvPtr, int bvID);
    public native void minimizeBV(long solverPtr,  long bvPtr, int bvID);
    public native void maximizeLits(long solverPtr, int [] lits, int n_lits);
    public native void minimizeLits(long solverPtr, int [] lits, int n_lits);
    public native void maximizeWeightedLits(long solverPtr, int [] lits, int [] weights, int n_lits);
    public native void minimizeWeightedLits(long solverPtr, int [] lits, int [] weights, int n_lits);

    //theory interface for bitvectors
    public native long initBVTheory(long solverPtr);
    public native int newBitvector_const(long solverPtr, long bvPtr, int bvWidth, long constval);
    public native int newBitvector_anon(long solverPtr, long bvPtr, int bvWidth);
    public native int newBitvector(long solverPtr, long bvPtr, int [] bits, int n_bits);
    public native int bv_width(long solverPtr, long bvPtr,int bvID);
    public native int newBVComparison_const_lt(long solverPtr, long bvPtr, int bvID, long constval);
    public native int newBVComparison_bv_lt(long solverPtr, long bvPtr, int bvID, int compareID);
    public native int newBVComparison_const_leq(long solverPtr, long bvPtr, int bvID, long constval);
    public native int newBVComparison_bv_leq(long solverPtr, long bvPtr, int bvID, int compareID);
    public native int newBVComparison_const_gt(long solverPtr, long bvPtr, int bvID, long constval);
    public native int newBVComparison_bv_gt(long solverPtr, long bvPtr, int bvID, int compareID);
    public native int newBVComparison_const_geq(long solverPtr, long bvPtr, int bvID, long constval);
    public native int newBVComparison_bv_geq(long solverPtr, long bvPtr, int bvID, int compareID);

    //Convert the specified bitvector, as well as any other bitvectors in its cone of influence, into pure CNF
    public native void bv_bitblast(long solverPtr, long bvPtr,int bvID);

    public native void bv_concat( long solverPtr, long bvPtr,int aID, int bID, int resultID);
    public native void bv_slice( long solverPtr, long bvPtr,int aID, int lower, int upper, int resultID);
    public native void bv_not( long solverPtr, long bvPtr,int bvaID, int bvResultID);
    public native void bv_and( long solverPtr, long bvPtr,int bvaID, int bvbID, int bvResultID);
    public native void bv_nand( long solverPtr, long bvPtr,int bvaID, int bvbID, int bvResultID);
    public native void bv_or( long solverPtr, long bvPtr,int bvaID, int bvbID, int bvResultID);
    public native void bv_nor( long solverPtr, long bvPtr,int bvaID, int bvbID, int bvResultID);
    public native void bv_xor( long solverPtr, long bvPtr,int bvaID, int bvbID, int bvResultID);
    public native void bv_xnor( long solverPtr, long bvPtr,int bvaID, int bvbID, int bvResultID);

    public native void bv_ite( long solverPtr, long bvPtr, int condition_lit,int bvThenID, int bvElseID, int bvResultID);

    public native void bv_addition( long solverPtr, long bvPtr, int bvID1, int bvID2, int resultID);
    public native void bv_subtraction( long solverPtr, long bvPtr, int bvID1, int bvID2, int resultID);
    public native void bv_multiply(long solverPtr, long bvPtr, int bvID1, int bvID2, int resultID);
    public native void bv_divide(long solverPtr, long bvPtr, int bvID1,  int bvID2, int resultID);
    public native void bv_min(long solverPtr, long bvPtr,  int[] args,int n_args,int resultID);
    public native void bv_max(long solverPtr, long bvPtr,  int[] args,int n_args,int resultID);
    public native void bv_popcount(long solverPtr, long bvPtr,  int[] args,int n_args, int resultID);
    public native void bv_unary(long solverPtr, long bvPtr, int [] args, int n_args, int resultID);

    //simple at-most-one constraint: asserts that at most one of the set of variables (NOT LITERALS) may be true.
    //for small numbers of variables, consider using a direct CNF encoding instead
    public native void at_most_one(long solverPtr, int [] vars, int n_vars);

    public native void assertPB_lt(long solverPtr, int rhs, int n_args, int [] literals, int [] coefficients);
    public native void assertPB_leq(long solverPtr, int rhs, int n_args, int [] literals, int [] coefficients);
    public native void assertPB_eq(long solverPtr, int rhs, int n_args, int [] literals, int [] coefficients);
    public native void assertPB_geq(long solverPtr, int rhs, int n_args, int [] literals, int [] coefficients);
    public native void assertPB_gt(long solverPtr, int rhs, int n_args, int [] literals, int [] coefficients);
    //Convert any pb constraints in the solver into cnf (will be called automatically before solve())
    public native void flushPB(long solverPtr);
    //theory interface for graphs
    public native long newGraph(long solverPtr);

    public native int newNode(long solverPtr,long graphPtr);
    public native int newEdge(long solverPtr, long graphPtr,int from,int  to,  long constweight);

    public native int newEdge_bv(long solverPtr, long graphPtr,int from,int  to, int bvID);
    public native int nNodes(long solverPtr,long graphPtr);
    public native int nEdges(long solverPtr,long graphPtr);
    public native int reaches(long solverPtr,long graphPtr,int from, int to);
    public native int shortestPathUnweighted_lt_const(long solverPtr,long graphPtr,int from, int to, int steps);
    public native int shortestPathUnweighted_leq_const(long solverPtr,long graphPtr,int from, int to, int steps);
    public native int shortestPath_lt_const(long solverPtr,long graphPtr,int from, int to, long dist);
    public native int shortestPath_leq_const(long solverPtr,long graphPtr,int from, int to, long dist);
    public native int shortestPath_lt_bv(long solverPtr,long graphPtr,int from, int to, int bvID);
    public native int shortestPath_leq_bv(long solverPtr,long graphPtr,int from, int to, int bvID);
    public native int maximumFlow_geq(long solverPtr,long graphPtr,int source, int sink, long constweight);
    public native int maximumFlow_gt(long solverPtr,long graphPtr,int source, int sink, long constweight);
    public native int maximumFlow_geq_bv(long solverPtr,long graphPtr,int source, int sink, int bvID);
    public native int maximumFlow_gt_bv(long solverPtr,long graphPtr,int source, int sink, int bvID);
    public native int minimumSpanningTree_leq(long solverPtr,long graphPtr, long constweight);
    public native int minimumSpanningTree_lt(long solverPtr,long graphPtr,int source, int sink, long constweight);
    public native int acyclic_undirected(long solverPtr,long graphPtr);
    public native int acyclic_directed(long solverPtr,long graphPtr);
    public native void newEdgeSet(long solverPtr,long graphPtr,int [] edges, int n_edges, boolean enforceEdgeAssignment);

    //this enables a heuristic on this graph, from the RUC paper, which sets assigned edges to zero long, to encourage edge-reuse in solutions
    public native void graph_setAssignEdgesToWeight(long solverPtr,long graphPtr, long weight);
    //flow routing interface

    public native long createFlowRouting(long solverPtr,long graphPtr, int sourceNode,int destNode,int maxflowLit);
    public native void addRoutingNet(long solverPtr,long graphPtr, long routerPtr, int disabledEdge, int n_members, int [] edge_lits, int [] reach_lits);

    //model query
    //For a given literal (not variable!), returns 0 for true, 1 for false, 2 for unassigned.
    public native int getModel_Literal(long solverPtr,int lit);
    //Get an assignment to a bitvector in the model. The model may find a range of satisfying assignments to the bitvector;
    //If getMaximumValue is true, this function returns the maximum satisfying assignment to the bitvector in the model; else it returns the smallest.
    public native long getModel_BV(long solverPtr, long bvPtr, int bvID, boolean getMaximumValue);
    //graph queries:
    //maxflow_literal is the literal (not variable!) that is the atom for the maximum flow query
    public native long getModel_MaxFlow(long solverPtr,long graphPtr,int maxflow_literal);
    //maxflow_literal is the literal (not variable!) that is the atom for the maximum flow query
    public native long getModel_EdgeFlow(long solverPtr,long graphPtr,int maxflow_literal, int edgeLit);
    public native long getModel_AcyclicEdgeFlow(long solverPtr,long graphPtr,int maxflow_literal, int edgeLit);

    public native long getModel_MinimumSpanningTreeWeight(long solverPtr,long graphPtr,int spanning_tree_literal);
    public native int getModel_Path_Nodes_Length(long solverPtr,long graphPtr,int reach_or_distance_literal);
    public native int getModel_Path_Nodes(long solverPtr,long graphPtr,int reach_or_distance_literal, int store_length, int [] store);

    public native int getModel_Path_EdgeLits_Length(long solverPtr,long graphPtr,int reach_or_distance_literal);
    public native int getModel_Path_EdgeLits(long solverPtr,long graphPtr,int reach_or_distance_literal, int store_length, int [] store);



    // Test Driver
    public static void main(String[] args) {
        String version =  new MonosatJNI().getVersion();  // invoke the native method
        System.out.println(version);
    }


}
