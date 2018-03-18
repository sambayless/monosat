from libc cimport bool
from libc.stdint cimport  int64_t
cdef extern from "monosat/api/Monosat.h":

    ctypedef void* SolverPtr

    ctypedef void* BVTheoryPtr

    ctypedef void* GraphTheorySolver_long

    ctypedef void* GraphTheorySolver_double

    ctypedef void* FSMTheorySolverPtr

    ctypedef void* FlowRouterPtr

    ctypedef int Var

    ctypedef int64_t Weight

    char* getVersion()

    SolverPtr newSolver()

    SolverPtr newSolver_arg(char* argv)

    

    void deleteSolver(SolverPtr S)

    void setOutputFile(SolverPtr S, char* output)

    void readGNF(SolverPtr S, char* filename)

    bool solve(SolverPtr S)

    bool solveAssumptions(SolverPtr S, int64_t assumptions, int n_assumptions)

    void setTimeLimit(SolverPtr S, int seconds)

    void setMemoryLimit(SolverPtr S, int mb)

    void setConflictLimit(SolverPtr S, int num_conflicts)

    void setPropagationLimit(SolverPtr S, int num_propagations)

    int solveLimited(SolverPtr S)

    int solveAssumptionsLimited(SolverPtr S, int64_t assumptions, int n_assumptions)

    bool lastSolutionWasOptimal(SolverPtr S)

    int getConflictClause(SolverPtr S, int64_t store_clause, int max_store_size)

    void backtrack(SolverPtr S)

    int newVar(SolverPtr S)

    void setDecisionVar(SolverPtr S, int var, bool decidable)

    bool isDecisionVar(SolverPtr S, int var)

    void setDecisionPriority(SolverPtr S, int var, int priority)

    int getDecisionPriority(SolverPtr S, int var)

    void setDecisionPolarity(SolverPtr S, Var v, bool b)

    bool getDecisionPolarity(SolverPtr S, Var v)

    int true_lit(SolverPtr S)

    bool disallowLiteralSimplification(SolverPtr S, int lit)

    void disablePreprocessing(SolverPtr S)

    int nVars(SolverPtr S)

    int nClauses(SolverPtr S)

    int nBitvectors(SolverPtr S, BVTheoryPtr bv)

    bool addClause(SolverPtr S, int64_t lits, int n_lits)

    bool addUnitClause(SolverPtr S, int lit)

    bool addBinaryClause(SolverPtr S, int lit1, int lit2)

    bool addTertiaryClause(SolverPtr S, int lit1, int lit2, int lit3)

    void addBinaryClauses(SolverPtr S, int64_t first_args, int64_t second_args, int n_pairs)

    void clearOptimizationObjectives(SolverPtr S)

    void maximizeBV(SolverPtr S, BVTheoryPtr bv, int bvID)

    void minimizeBV(SolverPtr S, BVTheoryPtr bv, int bvID)

    void maximizeLits(SolverPtr S, int64_t lits, int n_lits)

    void minimizeLits(SolverPtr S, int64_t lits, int n_lits)

    void maximizeWeightedLits(SolverPtr S, int64_t lits, int64_t weights, int n_lits)

    void minimizeWeightedLits(SolverPtr S, int64_t lits, int64_t weights, int n_lits)

    BVTheoryPtr initBVTheory(SolverPtr S)

    int newBitvector_const(SolverPtr S, BVTheoryPtr bv, int bvWidth, Weight constval)

    int newBitvector_anon(SolverPtr S, BVTheoryPtr bv, int bvWidth)

    int newBitvector(SolverPtr S, BVTheoryPtr bv, int64_t bits, int n_bits)

    int bv_width(SolverPtr S, BVTheoryPtr bv, int bvID)

    int newBVComparison_const_lt(SolverPtr S, BVTheoryPtr bv, int bvID, Weight weight)

    int newBVComparison_bv_lt(SolverPtr S, BVTheoryPtr bv, int bvID, int compareID)

    int newBVComparison_const_leq(SolverPtr S, BVTheoryPtr bv, int bvID, Weight weight)

    int newBVComparison_bv_leq(SolverPtr S, BVTheoryPtr bv, int bvID, int compareID)

    int newBVComparison_const_gt(SolverPtr S, BVTheoryPtr bv, int bvID, Weight weight)

    int newBVComparison_bv_gt(SolverPtr S, BVTheoryPtr bv, int bvID, int compareID)

    int newBVComparison_const_geq(SolverPtr S, BVTheoryPtr bv, int bvID, Weight weight)

    int newBVComparison_bv_geq(SolverPtr S, BVTheoryPtr bv, int bvID, int compareID)

    void bv_bitblast(SolverPtr S, BVTheoryPtr bv, int bvID)

    void bv_concat(SolverPtr S, BVTheoryPtr bv, int aID, int bID, int resultID)

    void bv_slice(SolverPtr S, BVTheoryPtr bv, int aID, int lower, int upper, int resultID)

    void bv_not(SolverPtr S, BVTheoryPtr bv, int bvaID, int bvResultID)

    void bv_and(SolverPtr S, BVTheoryPtr bv, int bvaID, int bvbID, int bvResultID)

    void bv_nand(SolverPtr S, BVTheoryPtr bv, int bvaID, int bvbID, int bvResultID)

    void bv_or(SolverPtr S, BVTheoryPtr bv, int bvaID, int bvbID, int bvResultID)

    void bv_nor(SolverPtr S, BVTheoryPtr bv, int bvaID, int bvbID, int bvResultID)

    void bv_xor(SolverPtr S, BVTheoryPtr bv, int bvaID, int bvbID, int bvResultID)

    void bv_xnor(SolverPtr S, BVTheoryPtr bv, int bvaID, int bvbID, int bvResultID)

    void bv_ite(SolverPtr S, BVTheoryPtr bv, int condition_lit, int bvThenID, int bvElseID, int bvResultID)

    void bv_addition(SolverPtr S, BVTheoryPtr bv, int bvID1, int bvID2, int resultID)

    void bv_subtraction(SolverPtr S, BVTheoryPtr bv, int bvID1, int bvID2, int resultID)

    void bv_multiply(SolverPtr S, BVTheoryPtr bv, int bvID1, int bvID2, int resultID)

    void bv_divide(SolverPtr S, BVTheoryPtr bv, int bvID1, int bvID2, int resultID)

    void bv_min(SolverPtr S, BVTheoryPtr bv, int64_t args, int n_args, int resultID)

    void bv_max(SolverPtr S, BVTheoryPtr bv, int64_t args, int n_args, int resultID)

    void bv_popcount(SolverPtr S, BVTheoryPtr bv, int64_t args, int n_args, int resultID)

    void bv_unary(SolverPtr S, BVTheoryPtr bv, int64_t args, int n_args, int resultID)

    void at_most_one(SolverPtr S, int64_t vars, int n_vars)

    void assertPB_lt(SolverPtr S, int rhs, int n_args, int64_t literals, int64_t coefficients)

    void assertPB_leq(SolverPtr S, int rhs, int n_args, int64_t literals, int64_t coefficients)

    void assertPB_eq(SolverPtr S, int rhs, int n_args, int64_t literals, int64_t coefficients)

    void assertPB_geq(SolverPtr S, int rhs, int n_args, int64_t literals, int64_t coefficients)

    void assertPB_gt(SolverPtr S, int rhs, int n_args, int64_t literals, int64_t coefficients)

    void flushPB(SolverPtr S)

    GraphTheorySolver_long newGraph(SolverPtr S)

    int newNode(SolverPtr S, GraphTheorySolver_long G)

    int newEdge(SolverPtr S, GraphTheorySolver_long G, int _from, int to, Weight weight)

    int newEdge_double(SolverPtr S, GraphTheorySolver_double G, int _from, int to, double weight)

    int newEdge_bv(SolverPtr S, GraphTheorySolver_long G, int _from, int to, int bvID)

    int nNodes(SolverPtr S, GraphTheorySolver_long G)

    int nEdges(SolverPtr S, GraphTheorySolver_long G)

    int reaches(SolverPtr S, GraphTheorySolver_long G, int _from, int to)

    int shortestPathUnweighted_lt_const(SolverPtr S, GraphTheorySolver_long G, int _from, int to, int steps)

    int shortestPathUnweighted_leq_const(SolverPtr S, GraphTheorySolver_long G, int _from, int to, int steps)

    int shortestPath_lt_const(SolverPtr S, GraphTheorySolver_long G, int _from, int to, Weight dist)

    int shortestPath_leq_const(SolverPtr S, GraphTheorySolver_long G, int _from, int to, Weight dist)

    int shortestPath_lt_bv(SolverPtr S, GraphTheorySolver_long G, int _from, int to, int bvID)

    int shortestPath_leq_bv(SolverPtr S, GraphTheorySolver_long G, int _from, int to, int bvID)

    int maximumFlow_geq(SolverPtr S, GraphTheorySolver_long G, int source, int sink, Weight weight)

    int maximumFlow_gt(SolverPtr S, GraphTheorySolver_long G, int source, int sink, Weight weight)

    int maximumFlow_geq_bv(SolverPtr S, GraphTheorySolver_long G, int source, int sink, int bvID)

    int maximumFlow_gt_bv(SolverPtr S, GraphTheorySolver_long G, int source, int sink, int bvID)

    int minimumSpanningTree_leq(SolverPtr S, GraphTheorySolver_long G, Weight weight)

    int minimumSpanningTree_lt(SolverPtr S, GraphTheorySolver_long G, int source, int sink, Weight weight)

    int acyclic_undirected(SolverPtr S, GraphTheorySolver_long G)

    int acyclic_directed(SolverPtr S, GraphTheorySolver_long G)

    void newEdgeSet(SolverPtr S, GraphTheorySolver_long G, int64_t edges, int n_edges, bool enforceEdgeAssignment)

    void graph_setAssignEdgesToWeight(SolverPtr S, GraphTheorySolver_long G, int64_t weight)

    FlowRouterPtr createFlowRouting(SolverPtr S, GraphTheorySolver_long G, int sourceNode, int destNode, int maxflowLit)

    void addRoutingNet(SolverPtr S, GraphTheorySolver_long G, FlowRouterPtr router, int disabledEdge, int n_members, int64_t edge_lits, int64_t reach_lits)

    FSMTheorySolverPtr initFSMTheory(SolverPtr S)

    int newFSM(SolverPtr S, FSMTheorySolverPtr fsmTheory, int inputAlphabet, int outputAlphabet)

    int newState(SolverPtr S, FSMTheorySolverPtr fsmTheory, int fsmID)

    int newTransition(SolverPtr S, FSMTheorySolverPtr fsmTheory, int fsmID, int _fromNode, int toNode, int inputLabel, int outputLabel)

    int newString(SolverPtr S, FSMTheorySolverPtr fsmTheory, int64_t str, int len)

    int fsmAcceptsString(SolverPtr S, FSMTheorySolverPtr fsmTheory, int fsmID, int startNode, int acceptNode, int stringID)

    int fsmCompositionAccepts(SolverPtr S, FSMTheorySolverPtr fsmTheory, int fsmGeneratorID, int fsmAcceptorID, int gen_startNode, int gen_acceptNode, int acceptor_startNode, int acceptor_acceptNode, int stringID)

    int getModel_Literal(SolverPtr S, int lit)

    Weight getModel_BV(SolverPtr S, BVTheoryPtr bv, int bvID, bool getMaximumValue)

    Weight getModel_MaxFlow(SolverPtr S, GraphTheorySolver_long G, int maxflow_literal)

    Weight getModel_EdgeFlow(SolverPtr S, GraphTheorySolver_long G, int maxflow_literal, int edgeLit)

    Weight getModel_AcyclicEdgeFlow(SolverPtr S, GraphTheorySolver_long G, int maxflow_literal, int edgeLit)

    Weight getModel_MinimumSpanningTreeWeight(SolverPtr S, GraphTheorySolver_long G, int spanning_tree_literal)

    int getModel_Path_Nodes_Length(SolverPtr S, GraphTheorySolver_long G, int reach_or_distance_literal)

    int getModel_Path_Nodes(SolverPtr S, GraphTheorySolver_long G, int reach_or_distance_literal, int store_length, int64_t store)

    int getModel_Path_EdgeLits_Length(SolverPtr S, GraphTheorySolver_long G, int reach_or_distance_literal)

    int getModel_Path_EdgeLits(SolverPtr S, GraphTheorySolver_long G, int reach_or_distance_literal, int store_length, int64_t store)
