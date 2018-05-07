/***************************************************************************************[Solver.cc]
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
 **************************************************************************************************/

package monosat;

import java.util.ArrayList;
import java.util.Vector;
import java.nio.IntBuffer;

/**
 * Low-level JNI mapping of the C interface to MonoSAT.
 * This class is not intended for end users, and its interface
 * may change without warning if the native library changes.
 * See instead monosat.Solver.
 */
final class MonosatJNI {//package level access specifier
    static {
        System.loadLibrary("monosat");
    }

    //Get the version string from the native library.
    public native static String getVersion();

    //literal/variable transformations
    public native static int varToLit(int var, boolean negated);

    public native static int litToVar(int lit);

    //create/destroy solver pointers
    public native static long newSolver();

    public native static long newSolver(String args);

    public native static void deleteSolver(long solverPtr);

    //If set, dump constraints to this file (as they are asserted in the solver)
    public native static void setOutputFile(long solverPtr, String filename);

    //load constraints from a file
    public native static void readGNF(long solverPtr, String filename);

    //basic solver functions
    public native static boolean solve(long solverPtr);

    public native static boolean solveAssumptions(long solverPtr, IntBuffer assumps, int n_assumptions);

    //Sets the (approximate) time limit in seconds before returning l_Undef from solveLimited; ignored by solve(). Set to <0 to disable time limit.
    public native static void setTimeLimit(long solverPtr, int seconds);

    //Sets the (approximate) memory limit in megabytes before returning l_Undef from solveLimited; ignored by solve(). Set to <0 to disable memory limit.
    public native static void setMemoryLimit(long solverPtr, int mb);

    //Sets the maximum number of (additional) conflicts allowed in the solver before returning l_Undef from solveLimited; ignored by solve(). Set to <0 to disable conflict limit.
    public native static void setConflictLimit(long solverPtr, int num_conflicts);

    //Sets the maximum number of (additional) propagations allowed in the solver before returning l_Undef from solveLimited; ignored by solve(). Set to <0 to disable propagation limit.
    public native static void setPropagationLimit(long solverPtr, int num_propagations);

    //Returns 0 for satisfiable, 1 for proved unsatisfiable, 2 for failed to find a solution (within any resource limits that have been set)
    public native static int solveLimited(long solverPtr);

    //Returns 0 for satisfiable, 1 for proved unsatisfiable, 2 for failed to find a solution (within any resource limits that have been set)
    public native static int solveAssumptionsLimited(long solverPtr, IntBuffer assumptions, int n_assumptions);

    public native static boolean lastSolutionWasOptimal(long solverPtr);

    //If the last solution was unsat, then this get the 'conflict clause' produced by the solver (a subset of the assumptions which are sufficient to cause the instance to be UNSAT).
    //Fills the given pointer with the first max_store_size literals of the conflict clause, and returns the number of literals in the conflict clause. Set store_clause to null and max_store_size to 0 to find the size of the conflict clause
    //Returns -1 if the solver has no conflict clause from the most recent solve() call (because that call was not UNSAT)
    public native static int getConflictClause(long solverPtr, IntBuffer store_clause, int length);

  //Reduce the given set of mutually UNSAT assumptions into a (locally) minimal-sized set of assumptions that are still UNSAT.
  //If the assumptions are not UNSAT, this method does returns -1.
  //Else, it returns the number of literals in the unsat core, which it will also store in the first number_literals
  //entries in unsat_assumptions
  //Note: This method will make repeated (and potentially expensive) calls to the SAT solver to attempt to remove literals from the
  //set of assumptions.
  //Note also that if any of the setXXXXXLimits() are applied to the solver, this may not produce a locally minimal unsat core.
    public native static int minimizeUnsatCore(long solverPtr, IntBuffer unsat_assumptions, int length);

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
    public native static void minimizeConflictClause(long solverPtr);

    public native static void backtrack(long solverPtr);

    public native static int newVar(long solverPtr);

    public native static void releaseLiteral(long solverPtr, int literal);

    public native static void setDecisionVar(long solverPtr, int var, boolean decidable);

    public native static boolean isDecisionVar(long solverPtr, int var);

    //Static, lexicographic heuristic. Larger values are higher priority (decided first); default priority is 0
    public native static void setDecisionPriority(long solverPtr, int var, int priority);

    public native static int getDecisionPriority(long solverPtr, int var);

    // Which polarity the decision heuristic should use for a variable (by default).
    public native static void setDecisionPolarity(long solverPtr, int v, boolean b);

    public native static boolean getDecisionPolarity(long solverPtr, int v);

    //The solver will (sometimes) instantiate an arbitrary true literal for use as a constant.
    //Call this method to a) force that literal to be instantiate, and b) get that literal.
    public native static int true_lit(long solverPtr);

    //Prevents this literal from being simplified by the preprocessor
    public native static boolean disallowLiteralSimplification(long solverPtr, int var);

    //permanently disable SAT-based preprocessing in this solver
    public native static void disablePreprocessing(long solverPtr);

    public native static int nVars(long solverPtr);

    public native static int nClauses(long solverPtr);

    public native static int nBitvectors(long solverPtr, long bvPtr);

    public native static boolean addClause(long solverPtr, IntBuffer lits, int n_lits);

    public native static boolean addUnitClause(long solverPtr, int lit);

    public native static boolean addBinaryClause(long solverPtr, int lit1, int lit2);

    public native static boolean addTertiaryClause(long solverPtr, int lit1, int lit2, int lit3);

    //remove any optimization objectives from the solver
    public native static void clearOptimizationObjectives(long solverPtr);

    public native static void maximizeBV(long solverPtr, long bvPtr, int bvID);

    public native static void minimizeBV(long solverPtr, long bvPtr, int bvID);

    public native static void maximizeLits(long solverPtr, IntBuffer lits, int n_lits);

    public native static void minimizeLits(long solverPtr, IntBuffer lits, int n_lits);

    public native static void maximizeWeightedLits(long solverPtr, IntBuffer lits, IntBuffer weights, int n_lits);

    public native static void minimizeWeightedLits(long solverPtr, IntBuffer lits, IntBuffer weights, int n_lits);

    //theory interface for bitvectors
    public native static long initBVTheory(long solverPtr);

    public native static int newBitvector_const(long solverPtr, long bvPtr, int bvWidth, long constval);

    public native static int newBitvector_anon(long solverPtr, long bvPtr, int bvWidth);

    public native static int newBitvector(long solverPtr, long bvPtr, IntBuffer bits, int n_bits);

    public native static int bv_width(long solverPtr, long bvPtr, int bvID);

    public native static int newBVComparison_const_lt(long solverPtr, long bvPtr, int bvID, long constval);

    public native static int newBVComparison_bv_lt(long solverPtr, long bvPtr, int bvID, int compareID);

    public native static int newBVComparison_const_leq(long solverPtr, long bvPtr, int bvID, long constval);

    public native static int newBVComparison_bv_leq(long solverPtr, long bvPtr, int bvID, int compareID);

    public native static int newBVComparison_const_gt(long solverPtr, long bvPtr, int bvID, long constval);

    public native static int newBVComparison_bv_gt(long solverPtr, long bvPtr, int bvID, int compareID);

    public native static int newBVComparison_const_geq(long solverPtr, long bvPtr, int bvID, long constval);

    public native static int newBVComparison_bv_geq(long solverPtr, long bvPtr, int bvID, int compareID);

    public native static int newBVComparison_const_eq(long solverPtr, long bvPtr, int bvID, long constval);

    public native static int newBVComparison_bv_eq(long solverPtr, long bvPtr, int bvID, int compareID);

    public native static int newBVComparison_const_neq(long solverPtr, long bvPtr, int bvID, long constval);

    public native static int newBVComparison_bv_neq(long solverPtr, long bvPtr, int bvID, int compareID);

    //Convert the specified bitvector, as well as any other bitvectors in its cone of influence, into pure CNF
    public native static void bv_bitblast(long solverPtr, long bvPtr, int bvID);

    public native static void bv_concat(long solverPtr, long bvPtr, int aID, int bID, int resultID);

    public native static void bv_slice(long solverPtr, long bvPtr, int aID, int lower, int upper, int resultID);

    public native static void bv_not(long solverPtr, long bvPtr, int bvaID, int bvResultID);

    public native static void bv_and(long solverPtr, long bvPtr, int bvaID, int bvbID, int bvResultID);

    public native static void bv_nand(long solverPtr, long bvPtr, int bvaID, int bvbID, int bvResultID);

    public native static void bv_or(long solverPtr, long bvPtr, int bvaID, int bvbID, int bvResultID);

    public native static void bv_nor(long solverPtr, long bvPtr, int bvaID, int bvbID, int bvResultID);

    public native static void bv_xor(long solverPtr, long bvPtr, int bvaID, int bvbID, int bvResultID);

    public native static void bv_xnor(long solverPtr, long bvPtr, int bvaID, int bvbID, int bvResultID);

    public native static void bv_ite(long solverPtr, long bvPtr, int condition_lit, int bvThenID, int bvElseID, int bvResultID);

    public native static void bv_addition(long solverPtr, long bvPtr, int bvID1, int bvID2, int resultID);

    public native static void bv_subtraction(long solverPtr, long bvPtr, int bvID1, int bvID2, int resultID);

    public native static void bv_multiply(long solverPtr, long bvPtr, int bvID1, int bvID2, int resultID);

    public native static void bv_divide(long solverPtr, long bvPtr, int bvID1, int bvID2, int resultID);

    public native static void bv_min(long solverPtr, long bvPtr, IntBuffer args, int n_args, int resultID);

    public native static void bv_max(long solverPtr, long bvPtr, IntBuffer args, int n_args, int resultID);

    public native static void bv_popcount(long solverPtr, long bvPtr, IntBuffer args, int n_args, int resultID);

    public native static void bv_unary(long solverPtr, long bvPtr, IntBuffer args, int n_args, int resultID);

    //simple at-most-one constraint: asserts that at most one of the set of variables (NOT LITERALS) may be true.
    //for small numbers of variables, consider using a direct CNF encoding instead
    public native static void at_most_one(long solverPtr, IntBuffer vars, int n_vars);

    public native static void assertPB_lt(long solverPtr, int rhs, int n_args, IntBuffer literals, IntBuffer coefficients);

    public native static void assertPB_leq(long solverPtr, int rhs, int n_args, IntBuffer literals, IntBuffer coefficients);

    public native static void assertPB_eq(long solverPtr, int rhs, int n_args, IntBuffer literals, IntBuffer coefficients);

    public native static void assertPB_geq(long solverPtr, int rhs, int n_args, IntBuffer literals, IntBuffer coefficients);

    public native static void assertPB_gt(long solverPtr, int rhs, int n_args, IntBuffer literals, IntBuffer coefficients);

    //Convert any pb constraints in the solver into cnf (will be called automatically before solve())
    public native static void flushPB(long solverPtr);

    //theory interface for graphs
    public native static long newGraph(long solverPtr);

    public native static int newNode(long solverPtr, long graphPtr);

    public native static int newEdge(long solverPtr, long graphPtr, int from, int to, long constweight);

    public native static int newEdge_bv(long solverPtr, long graphPtr, int from, int to, int bvID);

    public native static int nNodes(long solverPtr, long graphPtr);

    public native static int nEdges(long solverPtr, long graphPtr);

    public native static int reaches(long solverPtr, long graphPtr, int from, int to);

    public native static int reachesBackward(long solverPtr, long graphPtr, int from, int to);

    public native static int onPath(long solverPtr, long graphPtr,int nodeOnPath,int from, int to);

    public native static int shortestPathUnweighted_lt_const(long solverPtr, long graphPtr, int from, int to, int steps);

    public native static int shortestPathUnweighted_leq_const(long solverPtr, long graphPtr, int from, int to, int steps);

    public native static int shortestPath_lt_const(long solverPtr, long graphPtr, int from, int to, long dist);

    public native static int shortestPath_leq_const(long solverPtr, long graphPtr, int from, int to, long dist);

    public native static int shortestPath_lt_bv(long solverPtr, long graphPtr, int from, int to, int bvID);

    public native static int shortestPath_leq_bv(long solverPtr, long graphPtr, int from, int to, int bvID);

    public native static int maximumFlow_geq(long solverPtr, long graphPtr, int source, int sink, long constweight);

    public native static int maximumFlow_gt(long solverPtr, long graphPtr, int source, int sink, long constweight);

    public native static int maximumFlow_geq_bv(long solverPtr, long graphPtr, int source, int sink, int bvID);

    public native static int maximumFlow_gt_bv(long solverPtr, long graphPtr, int source, int sink, int bvID);

    public native static int minimumSpanningTree_leq(long solverPtr, long graphPtr, long constweight);

    public native static int minimumSpanningTree_lt(long solverPtr, long graphPtr, int source, int sink, long constweight);

    public native static int acyclic_undirected(long solverPtr, long graphPtr);

    public native static int acyclic_directed(long solverPtr, long graphPtr);

    public native static void newEdgeSet(long solverPtr, long graphPtr, IntBuffer edges, int n_edges, boolean enforceEdgeAssignment);

    //this enables a heuristic on this graph, from the RUC paper, which sets assigned edges to zero long, to encourage edge-reuse in solutions
    public native static void graph_setAssignEdgesToWeight(long solverPtr, long graphPtr, long weight);
    //flow routing interface

    public native static long createFlowRouting(long solverPtr, long graphPtr, int sourceNode, int destNode, int maxflowLit);

    public native static void addRoutingNet(long solverPtr, long graphPtr, long routerPtr, int disabledEdge, int n_members, IntBuffer edge_lits, IntBuffer reach_lits);

    //model query
    //Returns true if the solver has a model (in which case it is safe to query the model), false otherwise
    public native static boolean hasModel(long solverPtr);

    //For a given literal (not variable!), returns 0 for true, 1 for false, 2 for unassigned.
    public native static int getModel_Literal(long solverPtr, int lit);

    //Check if a literal is known to be a constant by the solver (eg, if it has an assignment at decision level 0).
    //For a given literal (not variable!), returns 0 for if the literal is a true constant, 1 if it is a false constant,
    //or 2 is it not assigned at level 0.
    public native static int getConstantModel_Literal(long solverPtr, int lit);

    //Get an assignment to a bitvector in the model. The model may find a range of satisfying assignments to the bitvector;
    //If getMaximumValue is true, this function returns the maximum satisfying assignment to the bitvector in the model; else it returns the smallest.
    public native static long getModel_BV(long solverPtr, long bvPtr, int bvID, boolean getMaximumValue);

    //graph queries:
    //maxflow_literal is the literal (not variable!) that is the atom for the maximum flow query
    public native static long getModel_MaxFlow(long solverPtr, long graphPtr, int maxflow_literal);

    //maxflow_literal is the literal (not variable!) that is the atom for the maximum flow query
    public native static long getModel_EdgeFlow(long solverPtr, long graphPtr, int maxflow_literal, int edgeLit);

    public native static long getModel_AcyclicEdgeFlow(long solverPtr, long graphPtr, int maxflow_literal, int edgeLit);

    public native static long getModel_MinimumSpanningTreeWeight(long solverPtr, long graphPtr, int spanning_tree_literal);

    public native static int getModel_Path_Nodes_Length(long solverPtr, long graphPtr, int reach_or_distance_literal);

    public native static int getModel_Path_Nodes(long solverPtr, long graphPtr, int reach_or_distance_literal, int store_length, IntBuffer store);

    public native static int getModel_Path_EdgeLits_Length(long solverPtr, long graphPtr, int reach_or_distance_literal);

    public native static int getModel_Path_EdgeLits(long solverPtr, long graphPtr, int reach_or_distance_literal, int store_length, IntBuffer store);

    //Logic building methods
    //Note: Many of these methods are capitalized, defying Java naming conventions.
    //This is because the lower case versions (and, or...) are keywords in C++, and so cannot be used in the JNI bridge.
    public native static int And_(long solverPtr, int lit_a, int lit_b, int lit_out);

    public native static int Ands_(long solverPtr, IntBuffer lits, int n_lits, int lit_out);

    public native static void AssertImpliesAnd_(long solverPtr, int implies, IntBuffer lits, int n_lits, int lit_out);

    public native static int Ands(long solverPtr, IntBuffer lits, int n_lits);

    public native static int And(long solverPtr, int lit_a, int lit_b);

    public native static int Or_(long solverPtr, int lit_a, int lit_b, int lit_out);

    public native static int Ors_(long solverPtr, IntBuffer lits, int n_lits, int lit_out);

    //If this gate is true, then all of vals must be true.
    //But if this gate is false, vals may be true or false.
    public native static int ImpliesAnd(long solverPtr, IntBuffer lits, int n_lits, int lit_out);

    //If this gate is true, then at least one of vals must be true.
    //But if this gate is false, vals may be true or false.
    public native static int ImpliesOr(long solverPtr, IntBuffer lits, int n_lits);

    public native static int ImpliesOr_(long solverPtr, IntBuffer lits, int n_lits, int lit_out);

    //This is an OR condition that holds only if implies is true
    public native static void AssertImpliesOr_(long solverPtr, int implies, IntBuffer lits, int n_lits, int lit_out);

    public native static int Ors(long solverPtr, IntBuffer lits, int n_lits);

    public native static int Or(long solverPtr, int lit_a, int lit_b);

    public native static int Nors(long solverPtr, IntBuffer lits, int n_lits);

    public native static int Nor(long solverPtr, int lit_a, int lit_b);

    public native static int Nands(long solverPtr, IntBuffer lits, int n_lits);

    public native static int Nand(long solverPtr, int lit_a, int lit_b);

    public native static int Xors(long solverPtr, IntBuffer lits, int n_lits);

    public native static int Xor(long solverPtr, int lit_a, int lit_b);


    public native static int Xnors(long solverPtr, IntBuffer lits, int n_lits);

    public native static int Xnor(long solverPtr, int lit_a, int lit_b);

    public native static int Implies(long solverPtr, int lit_a, int lit_b);

    public native static int Implies_(long solverPtr, int lit_a, int lit_b, int lit_out);

    public native static int Ite(long solverPtr, int lit_cond, int lit_thn, int lit_els);

    public native static int Ite_(long solverPtr, int lit_cond, int lit_thn, int lit_els, int lit_out);

    public native static int Add(long solverPtr, IntBuffer lits_a, IntBuffer lits_b, int n_lits, IntBuffer lits_out);

    public native static int Add_(long solverPtr, IntBuffer lits_a, IntBuffer lits_b, int n_lits, IntBuffer lits_out, int carry_lit);

    public native static int Subtract(long solverPtr, IntBuffer lits_a, IntBuffer lits_b, int n_lits, IntBuffer lits_out);

    public native static int Subtract_(long solverPtr, IntBuffer lits_a, IntBuffer lits_b, int n_lits, IntBuffer lits_out, int borrow_lit);

    //perform two's complement negation (long solverPtr,invert bits and add 1)
//IntBuffer lits_out must have enough room to hold n_lits
    public native static void Negate(long solverPtr, IntBuffer lits, int n_lits, IntBuffer lits_out);

    //perform two's complement negation (long solverPtr,invert bits and add 1)
    public native static void Negate_(long solverPtr, IntBuffer lits, int n_lits, IntBuffer lits_out);

    public native static void Assert(long solverPtr, int lit);

    public native static void AssertOrTertiary(long solverPtr, int lit_a, int lit_b, int lit_c);

    public native static void AssertOrs(long solverPtr, IntBuffer lits, int n_lits);

    public native static void AssertOr(long solverPtr, int lit_a, int lit_b);

    public native static void AssertNands(long solverPtr, IntBuffer lits, int n_lits);

    public native static void AssertNand(long solverPtr, int lit_a, int lit_b);

    public native static void AssertAnds(long solverPtr, IntBuffer lits, int n_lits);

    public native static void AssertAnd(long solverPtr, int lit_a, int lit_b);

    public native static void AssertNors(long solverPtr, IntBuffer lits, int n_lits);

    public native static void AssertNor(long solverPtr, int lit_a, int lit_b);

    public native static void AssertXor(long solverPtr, int lit_a, int lit_b);

    public native static void AssertXors(long solverPtr, IntBuffer lits, int n_lits);

    public native static void AssertXnors(long solverPtr, IntBuffer lits, int n_lits);

    public native static void AssertXnor(long solverPtr, int lit_a, int lit_b);

    public native static void AssertImplies(long solverPtr, int lit_a, int lit_b);

    public native static void AssertEqual(long solverPtr, int lit_a, int lit_b);

    public native static void AssertAllSame(long solverPtr, IntBuffer lits, int n_lits);

    public native static int Equals(long solverPtr, IntBuffer A_lits, IntBuffer B_lits, int n_lits);

    public native static int LEQ(long solverPtr, IntBuffer A_lits, IntBuffer B_lits, int n_lits);

    public native static int LT(long solverPtr, IntBuffer A_lits, IntBuffer B_lits, int n_lits);

    public native static void AssertEquals(long solverPtr, IntBuffer A_lits, IntBuffer B_lits, int n_lits);

    public native static void AssertLEQ(long solverPtr, IntBuffer A_lits, IntBuffer B_lits, int n_lits);

    public native static void AssertLT(long solverPtr, IntBuffer A_lits, IntBuffer B_lits, int n_lits);

    //uses n^2 binary clauses to create a simple at-most-one constraint.
//if you have more than 20 or so literals, strongly consider using a pseudo-Boolean constraint solver instead
    public native static void AssertAMO(long solverPtr, IntBuffer lits, int n_lits);

    //uses n^2 binary clauses to create a simple exactly-one-constraint.
//if you have more than 20 or so literals, strongly consider using a pseudo-Boolean constraint solver instead
    public native static void AssertExactlyOne(long solverPtr, IntBuffer lits, int n_lits);

   /**
   * Simple test driver, useful for debugging native library linking issues.
    */
    @SuppressWarnings("unchecked")
    public static void main(String[] args) {
        System.out.println("Loading native libraries:");
        try {
            java.lang.reflect.Field LIBRARIES = ClassLoader.class.getDeclaredField("loadedLibraryNames");
            LIBRARIES.setAccessible(true);
            Vector<String> libraries = (Vector<String>) LIBRARIES.get(ClassLoader.getSystemClassLoader());
            for (String s : libraries)
                System.out.println(s);

            String version = MonosatJNI.getVersion();  // invoke the native method
            System.out.println("Loaded MonoSAT: " + version);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }


}
