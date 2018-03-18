#cython: c_string_encoding=ascii  # for cython>=0.19
#cython: embedsignature=False
from  libcpp.string  cimport string as libcpp_string
from  libcpp.string  cimport string as libcpp_utf8_string
from  libcpp.set     cimport set as libcpp_set
from  libcpp.vector  cimport vector as libcpp_vector
from  libcpp.pair    cimport pair as libcpp_pair
from  libcpp.map     cimport map  as libcpp_map
from  libcpp cimport bool
from  libc.string cimport const_char
from cython.operator cimport dereference as deref, preincrement as inc, address as address


from monosat cimport BVTheoryPtr
from monosat cimport FSMTheorySolverPtr
from monosat cimport FlowRouterPtr
from monosat cimport GraphTheorySolver_double
from monosat cimport GraphTheorySolver_long
from monosat cimport SolverPtr
from monosat cimport Var
from monosat cimport Weight
from monosat cimport acyclic_directed as _acyclic_directed_monosat
from monosat cimport acyclic_undirected as _acyclic_undirected_monosat
from monosat cimport addBinaryClause as _addBinaryClause_monosat
from monosat cimport addBinaryClauses as _addBinaryClauses_monosat
from monosat cimport addClause as _addClause_monosat
from monosat cimport addRoutingNet as _addRoutingNet_monosat
from monosat cimport addTertiaryClause as _addTertiaryClause_monosat
from monosat cimport addUnitClause as _addUnitClause_monosat
from monosat cimport assertPB_eq as _assertPB_eq_monosat
from monosat cimport assertPB_geq as _assertPB_geq_monosat
from monosat cimport assertPB_gt as _assertPB_gt_monosat
from monosat cimport assertPB_leq as _assertPB_leq_monosat
from monosat cimport assertPB_lt as _assertPB_lt_monosat
from monosat cimport at_most_one as _at_most_one_monosat
from monosat cimport backtrack as _backtrack_monosat
from monosat cimport bv_addition as _bv_addition_monosat
from monosat cimport bv_and as _bv_and_monosat
from monosat cimport bv_bitblast as _bv_bitblast_monosat
from monosat cimport bv_concat as _bv_concat_monosat
from monosat cimport bv_divide as _bv_divide_monosat
from monosat cimport bv_ite as _bv_ite_monosat
from monosat cimport bv_max as _bv_max_monosat
from monosat cimport bv_min as _bv_min_monosat
from monosat cimport bv_multiply as _bv_multiply_monosat
from monosat cimport bv_nand as _bv_nand_monosat
from monosat cimport bv_nor as _bv_nor_monosat
from monosat cimport bv_not as _bv_not_monosat
from monosat cimport bv_or as _bv_or_monosat
from monosat cimport bv_popcount as _bv_popcount_monosat
from monosat cimport bv_slice as _bv_slice_monosat
from monosat cimport bv_subtraction as _bv_subtraction_monosat
from monosat cimport bv_unary as _bv_unary_monosat
from monosat cimport bv_width as _bv_width_monosat
from monosat cimport bv_xnor as _bv_xnor_monosat
from monosat cimport bv_xor as _bv_xor_monosat
from monosat cimport clearOptimizationObjectives as _clearOptimizationObjectives_monosat
from monosat cimport createFlowRouting as _createFlowRouting_monosat
from monosat cimport deleteSolver as _deleteSolver_monosat
from monosat cimport disablePreprocessing as _disablePreprocessing_monosat
from monosat cimport disallowLiteralSimplification as _disallowLiteralSimplification_monosat
from monosat cimport flushPB as _flushPB_monosat
from monosat cimport fsmAcceptsString as _fsmAcceptsString_monosat
from monosat cimport fsmCompositionAccepts as _fsmCompositionAccepts_monosat
from monosat cimport getConflictClause as _getConflictClause_monosat
from monosat cimport getDecisionPolarity as _getDecisionPolarity_monosat
from monosat cimport getDecisionPriority as _getDecisionPriority_monosat
from monosat cimport getModel_AcyclicEdgeFlow as _getModel_AcyclicEdgeFlow_monosat
from monosat cimport getModel_BV as _getModel_BV_monosat
from monosat cimport getModel_EdgeFlow as _getModel_EdgeFlow_monosat
from monosat cimport getModel_Literal as _getModel_Literal_monosat
from monosat cimport getModel_MaxFlow as _getModel_MaxFlow_monosat
from monosat cimport getModel_MinimumSpanningTreeWeight as _getModel_MinimumSpanningTreeWeight_monosat
from monosat cimport getModel_Path_EdgeLits as _getModel_Path_EdgeLits_monosat
from monosat cimport getModel_Path_EdgeLits_Length as _getModel_Path_EdgeLits_Length_monosat
from monosat cimport getModel_Path_Nodes as _getModel_Path_Nodes_monosat
from monosat cimport getModel_Path_Nodes_Length as _getModel_Path_Nodes_Length_monosat
from monosat cimport getVersion as _getVersion_monosat
from monosat cimport graph_setAssignEdgesToWeight as _graph_setAssignEdgesToWeight_monosat
from monosat cimport initBVTheory as _initBVTheory_monosat
from monosat cimport initFSMTheory as _initFSMTheory_monosat
from monosat cimport isDecisionVar as _isDecisionVar_monosat
from monosat cimport lastSolutionWasOptimal as _lastSolutionWasOptimal_monosat
from monosat cimport maximizeBV as _maximizeBV_monosat
from monosat cimport maximizeLits as _maximizeLits_monosat
from monosat cimport maximizeWeightedLits as _maximizeWeightedLits_monosat
from monosat cimport maximumFlow_geq as _maximumFlow_geq_monosat
from monosat cimport maximumFlow_geq_bv as _maximumFlow_geq_bv_monosat
from monosat cimport maximumFlow_gt as _maximumFlow_gt_monosat
from monosat cimport maximumFlow_gt_bv as _maximumFlow_gt_bv_monosat
from monosat cimport minimizeBV as _minimizeBV_monosat
from monosat cimport minimizeLits as _minimizeLits_monosat
from monosat cimport minimizeWeightedLits as _minimizeWeightedLits_monosat
from monosat cimport minimumSpanningTree_leq as _minimumSpanningTree_leq_monosat
from monosat cimport minimumSpanningTree_lt as _minimumSpanningTree_lt_monosat
from monosat cimport nBitvectors as _nBitvectors_monosat
from monosat cimport nClauses as _nClauses_monosat
from monosat cimport nEdges as _nEdges_monosat
from monosat cimport nNodes as _nNodes_monosat
from monosat cimport nVars as _nVars_monosat
from monosat cimport newBVComparison_bv_geq as _newBVComparison_bv_geq_monosat
from monosat cimport newBVComparison_bv_gt as _newBVComparison_bv_gt_monosat
from monosat cimport newBVComparison_bv_leq as _newBVComparison_bv_leq_monosat
from monosat cimport newBVComparison_bv_lt as _newBVComparison_bv_lt_monosat
from monosat cimport newBVComparison_const_geq as _newBVComparison_const_geq_monosat
from monosat cimport newBVComparison_const_gt as _newBVComparison_const_gt_monosat
from monosat cimport newBVComparison_const_leq as _newBVComparison_const_leq_monosat
from monosat cimport newBVComparison_const_lt as _newBVComparison_const_lt_monosat
from monosat cimport newBitvector as _newBitvector_monosat
from monosat cimport newBitvector_anon as _newBitvector_anon_monosat
from monosat cimport newBitvector_const as _newBitvector_const_monosat
from monosat cimport newEdge as _newEdge_monosat
from monosat cimport newEdgeSet as _newEdgeSet_monosat
from monosat cimport newEdge_bv as _newEdge_bv_monosat
from monosat cimport newEdge_double as _newEdge_double_monosat
from monosat cimport newFSM as _newFSM_monosat
from monosat cimport newGraph as _newGraph_monosat
from monosat cimport newNode as _newNode_monosat
from monosat cimport newSolver as _newSolver_monosat
from monosat cimport newSolver_arg as _newSolver_arg_monosat
from monosat cimport newState as _newState_monosat
from monosat cimport newString as _newString_monosat
from monosat cimport newTransition as _newTransition_monosat
from monosat cimport newVar as _newVar_monosat
from monosat cimport reaches as _reaches_monosat
from monosat cimport readGNF as _readGNF_monosat
from monosat cimport setConflictLimit as _setConflictLimit_monosat
from monosat cimport setDecisionPolarity as _setDecisionPolarity_monosat
from monosat cimport setDecisionPriority as _setDecisionPriority_monosat
from monosat cimport setDecisionVar as _setDecisionVar_monosat
from monosat cimport setMemoryLimit as _setMemoryLimit_monosat
from monosat cimport setOutputFile as _setOutputFile_monosat
from monosat cimport setPropagationLimit as _setPropagationLimit_monosat
from monosat cimport setTimeLimit as _setTimeLimit_monosat
from monosat cimport shortestPathUnweighted_leq_const as _shortestPathUnweighted_leq_const_monosat
from monosat cimport shortestPathUnweighted_lt_const as _shortestPathUnweighted_lt_const_monosat
from monosat cimport shortestPath_leq_bv as _shortestPath_leq_bv_monosat
from monosat cimport shortestPath_leq_const as _shortestPath_leq_const_monosat
from monosat cimport shortestPath_lt_bv as _shortestPath_lt_bv_monosat
from monosat cimport shortestPath_lt_const as _shortestPath_lt_const_monosat
from monosat cimport solve as _solve_monosat
from monosat cimport solveAssumptions as _solveAssumptions_monosat
from monosat cimport solveAssumptionsLimited as _solveAssumptionsLimited_monosat
from monosat cimport solveLimited as _solveLimited_monosat
from monosat cimport true_lit as _true_lit_monosat


def acyclic_directed( S ,  G ):
    """Cython signature: int acyclic_directed(void* S, void* G)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(G, (int, long)), 'arg G wrong type'


    cdef int _r = _acyclic_directed_monosat((<void*>S), (<void*>G))
    py_result = <int>_r
    return py_result

def acyclic_undirected( S ,  G ):
    """Cython signature: int acyclic_undirected(void* S, void* G)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(G, (int, long)), 'arg G wrong type'


    cdef int _r = _acyclic_undirected_monosat((<void*>S), (<void*>G))
    py_result = <int>_r
    return py_result

def addBinaryClause( S ,  lit1 ,  lit2 ):
    """Cython signature: bool addBinaryClause(void* S, int lit1, int lit2)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(lit1, (int, long)), 'arg lit1 wrong type'
    assert isinstance(lit2, (int, long)), 'arg lit2 wrong type'



    cdef bool _r = _addBinaryClause_monosat((<void*>S), (<int>lit1), (<int>lit2))
    py_result = <bool>_r
    return py_result

def addBinaryClauses( S ,  first_args ,  second_args ,  n_pairs ):
    """Cython signature: void addBinaryClauses(void* S, int64_t first_args, int64_t second_args, int n_pairs)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(first_args, (int, long)), 'arg first_args wrong type'
    assert isinstance(second_args, (int, long)), 'arg second_args wrong type'
    assert isinstance(n_pairs, (int, long)), 'arg n_pairs wrong type'




    _addBinaryClauses_monosat((<void*>S), (<int64_t>first_args), (<int64_t>second_args), (<int>n_pairs))

def addClause( S ,  lits ,  n_lits ):
    """Cython signature: bool addClause(void* S, int64_t lits, int n_lits)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(lits, (int, long)), 'arg lits wrong type'
    assert isinstance(n_lits, (int, long)), 'arg n_lits wrong type'



    cdef bool _r = _addClause_monosat((<void*>S), (<int64_t>lits), (<int>n_lits))
    py_result = <bool>_r
    return py_result

def addRoutingNet( S ,  G ,  router ,  disabledEdge ,  n_members ,  edge_lits ,  reach_lits ):
    """Cython signature: void addRoutingNet(void* S, void* G, void* router, int disabledEdge, int n_members, int64_t edge_lits, int64_t reach_lits)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(G, (int, long)), 'arg G wrong type'
    assert isinstance(router, (int, long)), 'arg router wrong type'
    assert isinstance(disabledEdge, (int, long)), 'arg disabledEdge wrong type'
    assert isinstance(n_members, (int, long)), 'arg n_members wrong type'
    assert isinstance(edge_lits, (int, long)), 'arg edge_lits wrong type'
    assert isinstance(reach_lits, (int, long)), 'arg reach_lits wrong type'







    _addRoutingNet_monosat((<void*>S), (<void*>G), (<void*>router), (<int>disabledEdge), (<int>n_members), (<int64_t>edge_lits), (<int64_t>reach_lits))

def addTertiaryClause( S ,  lit1 ,  lit2 ,  lit3 ):
    """Cython signature: bool addTertiaryClause(void* S, int lit1, int lit2, int lit3)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(lit1, (int, long)), 'arg lit1 wrong type'
    assert isinstance(lit2, (int, long)), 'arg lit2 wrong type'
    assert isinstance(lit3, (int, long)), 'arg lit3 wrong type'




    cdef bool _r = _addTertiaryClause_monosat((<void*>S), (<int>lit1), (<int>lit2), (<int>lit3))
    py_result = <bool>_r
    return py_result

def addUnitClause( S ,  lit ):
    """Cython signature: bool addUnitClause(void* S, int lit)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(lit, (int, long)), 'arg lit wrong type'


    cdef bool _r = _addUnitClause_monosat((<void*>S), (<int>lit))
    py_result = <bool>_r
    return py_result

def assertPB_eq( S ,  rhs ,  n_args ,  literals ,  coefficients ):
    """Cython signature: void assertPB_eq(void* S, int rhs, int n_args, int64_t literals, int64_t coefficients)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(rhs, (int, long)), 'arg rhs wrong type'
    assert isinstance(n_args, (int, long)), 'arg n_args wrong type'
    assert isinstance(literals, (int, long)), 'arg literals wrong type'
    assert isinstance(coefficients, (int, long)), 'arg coefficients wrong type'





    _assertPB_eq_monosat((<void*>S), (<int>rhs), (<int>n_args), (<int64_t>literals), (<int64_t>coefficients))

def assertPB_geq( S ,  rhs ,  n_args ,  literals ,  coefficients ):
    """Cython signature: void assertPB_geq(void* S, int rhs, int n_args, int64_t literals, int64_t coefficients)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(rhs, (int, long)), 'arg rhs wrong type'
    assert isinstance(n_args, (int, long)), 'arg n_args wrong type'
    assert isinstance(literals, (int, long)), 'arg literals wrong type'
    assert isinstance(coefficients, (int, long)), 'arg coefficients wrong type'





    _assertPB_geq_monosat((<void*>S), (<int>rhs), (<int>n_args), (<int64_t>literals), (<int64_t>coefficients))

def assertPB_gt( S ,  rhs ,  n_args ,  literals ,  coefficients ):
    """Cython signature: void assertPB_gt(void* S, int rhs, int n_args, int64_t literals, int64_t coefficients)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(rhs, (int, long)), 'arg rhs wrong type'
    assert isinstance(n_args, (int, long)), 'arg n_args wrong type'
    assert isinstance(literals, (int, long)), 'arg literals wrong type'
    assert isinstance(coefficients, (int, long)), 'arg coefficients wrong type'





    _assertPB_gt_monosat((<void*>S), (<int>rhs), (<int>n_args), (<int64_t>literals), (<int64_t>coefficients))

def assertPB_leq( S ,  rhs ,  n_args ,  literals ,  coefficients ):
    """Cython signature: void assertPB_leq(void* S, int rhs, int n_args, int64_t literals, int64_t coefficients)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(rhs, (int, long)), 'arg rhs wrong type'
    assert isinstance(n_args, (int, long)), 'arg n_args wrong type'
    assert isinstance(literals, (int, long)), 'arg literals wrong type'
    assert isinstance(coefficients, (int, long)), 'arg coefficients wrong type'





    _assertPB_leq_monosat((<void*>S), (<int>rhs), (<int>n_args), (<int64_t>literals), (<int64_t>coefficients))

def assertPB_lt( S ,  rhs ,  n_args ,  literals ,  coefficients ):
    """Cython signature: void assertPB_lt(void* S, int rhs, int n_args, int64_t literals, int64_t coefficients)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(rhs, (int, long)), 'arg rhs wrong type'
    assert isinstance(n_args, (int, long)), 'arg n_args wrong type'
    assert isinstance(literals, (int, long)), 'arg literals wrong type'
    assert isinstance(coefficients, (int, long)), 'arg coefficients wrong type'





    _assertPB_lt_monosat((<void*>S), (<int>rhs), (<int>n_args), (<int64_t>literals), (<int64_t>coefficients))

def at_most_one( S ,  vars ,  n_vars ):
    """Cython signature: void at_most_one(void* S, int64_t vars, int n_vars)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(vars, (int, long)), 'arg vars wrong type'
    assert isinstance(n_vars, (int, long)), 'arg n_vars wrong type'



    _at_most_one_monosat((<void*>S), (<int64_t>vars), (<int>n_vars))

def backtrack( S ):
    """Cython signature: void backtrack(void* S)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'

    _backtrack_monosat((<void*>S))

def bv_addition( S ,  bv ,  bvID1 ,  bvID2 ,  resultID ):
    """Cython signature: void bv_addition(void* S, void* bv, int bvID1, int bvID2, int resultID)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(bv, (int, long)), 'arg bv wrong type'
    assert isinstance(bvID1, (int, long)), 'arg bvID1 wrong type'
    assert isinstance(bvID2, (int, long)), 'arg bvID2 wrong type'
    assert isinstance(resultID, (int, long)), 'arg resultID wrong type'





    _bv_addition_monosat((<void*>S), (<void*>bv), (<int>bvID1), (<int>bvID2), (<int>resultID))

def bv_and( S ,  bv ,  bvaID ,  bvbID ,  bvResultID ):
    """Cython signature: void bv_and(void* S, void* bv, int bvaID, int bvbID, int bvResultID)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(bv, (int, long)), 'arg bv wrong type'
    assert isinstance(bvaID, (int, long)), 'arg bvaID wrong type'
    assert isinstance(bvbID, (int, long)), 'arg bvbID wrong type'
    assert isinstance(bvResultID, (int, long)), 'arg bvResultID wrong type'





    _bv_and_monosat((<void*>S), (<void*>bv), (<int>bvaID), (<int>bvbID), (<int>bvResultID))

def bv_bitblast( S ,  bv ,  bvID ):
    """Cython signature: void bv_bitblast(void* S, void* bv, int bvID)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(bv, (int, long)), 'arg bv wrong type'
    assert isinstance(bvID, (int, long)), 'arg bvID wrong type'



    _bv_bitblast_monosat((<void*>S), (<void*>bv), (<int>bvID))

def bv_concat( S ,  bv ,  aID ,  bID ,  resultID ):
    """Cython signature: void bv_concat(void* S, void* bv, int aID, int bID, int resultID)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(bv, (int, long)), 'arg bv wrong type'
    assert isinstance(aID, (int, long)), 'arg aID wrong type'
    assert isinstance(bID, (int, long)), 'arg bID wrong type'
    assert isinstance(resultID, (int, long)), 'arg resultID wrong type'





    _bv_concat_monosat((<void*>S), (<void*>bv), (<int>aID), (<int>bID), (<int>resultID))

def bv_divide( S ,  bv ,  bvID1 ,  bvID2 ,  resultID ):
    """Cython signature: void bv_divide(void* S, void* bv, int bvID1, int bvID2, int resultID)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(bv, (int, long)), 'arg bv wrong type'
    assert isinstance(bvID1, (int, long)), 'arg bvID1 wrong type'
    assert isinstance(bvID2, (int, long)), 'arg bvID2 wrong type'
    assert isinstance(resultID, (int, long)), 'arg resultID wrong type'





    _bv_divide_monosat((<void*>S), (<void*>bv), (<int>bvID1), (<int>bvID2), (<int>resultID))

def bv_ite( S ,  bv ,  condition_lit ,  bvThenID ,  bvElseID ,  bvResultID ):
    """Cython signature: void bv_ite(void* S, void* bv, int condition_lit, int bvThenID, int bvElseID, int bvResultID)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(bv, (int, long)), 'arg bv wrong type'
    assert isinstance(condition_lit, (int, long)), 'arg condition_lit wrong type'
    assert isinstance(bvThenID, (int, long)), 'arg bvThenID wrong type'
    assert isinstance(bvElseID, (int, long)), 'arg bvElseID wrong type'
    assert isinstance(bvResultID, (int, long)), 'arg bvResultID wrong type'






    _bv_ite_monosat((<void*>S), (<void*>bv), (<int>condition_lit), (<int>bvThenID), (<int>bvElseID), (<int>bvResultID))

def bv_max( S ,  bv ,  args ,  n_args ,  resultID ):
    """Cython signature: void bv_max(void* S, void* bv, int64_t args, int n_args, int resultID)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(bv, (int, long)), 'arg bv wrong type'
    assert isinstance(args, (int, long)), 'arg args wrong type'
    assert isinstance(n_args, (int, long)), 'arg n_args wrong type'
    assert isinstance(resultID, (int, long)), 'arg resultID wrong type'





    _bv_max_monosat((<void*>S), (<void*>bv), (<int64_t>args), (<int>n_args), (<int>resultID))

def bv_min( S ,  bv ,  args ,  n_args ,  resultID ):
    """Cython signature: void bv_min(void* S, void* bv, int64_t args, int n_args, int resultID)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(bv, (int, long)), 'arg bv wrong type'
    assert isinstance(args, (int, long)), 'arg args wrong type'
    assert isinstance(n_args, (int, long)), 'arg n_args wrong type'
    assert isinstance(resultID, (int, long)), 'arg resultID wrong type'





    _bv_min_monosat((<void*>S), (<void*>bv), (<int64_t>args), (<int>n_args), (<int>resultID))

def bv_multiply( S ,  bv ,  bvID1 ,  bvID2 ,  resultID ):
    """Cython signature: void bv_multiply(void* S, void* bv, int bvID1, int bvID2, int resultID)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(bv, (int, long)), 'arg bv wrong type'
    assert isinstance(bvID1, (int, long)), 'arg bvID1 wrong type'
    assert isinstance(bvID2, (int, long)), 'arg bvID2 wrong type'
    assert isinstance(resultID, (int, long)), 'arg resultID wrong type'





    _bv_multiply_monosat((<void*>S), (<void*>bv), (<int>bvID1), (<int>bvID2), (<int>resultID))

def bv_nand( S ,  bv ,  bvaID ,  bvbID ,  bvResultID ):
    """Cython signature: void bv_nand(void* S, void* bv, int bvaID, int bvbID, int bvResultID)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(bv, (int, long)), 'arg bv wrong type'
    assert isinstance(bvaID, (int, long)), 'arg bvaID wrong type'
    assert isinstance(bvbID, (int, long)), 'arg bvbID wrong type'
    assert isinstance(bvResultID, (int, long)), 'arg bvResultID wrong type'





    _bv_nand_monosat((<void*>S), (<void*>bv), (<int>bvaID), (<int>bvbID), (<int>bvResultID))

def bv_nor( S ,  bv ,  bvaID ,  bvbID ,  bvResultID ):
    """Cython signature: void bv_nor(void* S, void* bv, int bvaID, int bvbID, int bvResultID)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(bv, (int, long)), 'arg bv wrong type'
    assert isinstance(bvaID, (int, long)), 'arg bvaID wrong type'
    assert isinstance(bvbID, (int, long)), 'arg bvbID wrong type'
    assert isinstance(bvResultID, (int, long)), 'arg bvResultID wrong type'





    _bv_nor_monosat((<void*>S), (<void*>bv), (<int>bvaID), (<int>bvbID), (<int>bvResultID))

def bv_not( S ,  bv ,  bvaID ,  bvResultID ):
    """Cython signature: void bv_not(void* S, void* bv, int bvaID, int bvResultID)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(bv, (int, long)), 'arg bv wrong type'
    assert isinstance(bvaID, (int, long)), 'arg bvaID wrong type'
    assert isinstance(bvResultID, (int, long)), 'arg bvResultID wrong type'




    _bv_not_monosat((<void*>S), (<void*>bv), (<int>bvaID), (<int>bvResultID))

def bv_or( S ,  bv ,  bvaID ,  bvbID ,  bvResultID ):
    """Cython signature: void bv_or(void* S, void* bv, int bvaID, int bvbID, int bvResultID)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(bv, (int, long)), 'arg bv wrong type'
    assert isinstance(bvaID, (int, long)), 'arg bvaID wrong type'
    assert isinstance(bvbID, (int, long)), 'arg bvbID wrong type'
    assert isinstance(bvResultID, (int, long)), 'arg bvResultID wrong type'





    _bv_or_monosat((<void*>S), (<void*>bv), (<int>bvaID), (<int>bvbID), (<int>bvResultID))

def bv_popcount( S ,  bv ,  args ,  n_args ,  resultID ):
    """Cython signature: void bv_popcount(void* S, void* bv, int64_t args, int n_args, int resultID)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(bv, (int, long)), 'arg bv wrong type'
    assert isinstance(args, (int, long)), 'arg args wrong type'
    assert isinstance(n_args, (int, long)), 'arg n_args wrong type'
    assert isinstance(resultID, (int, long)), 'arg resultID wrong type'





    _bv_popcount_monosat((<void*>S), (<void*>bv), (<int64_t>args), (<int>n_args), (<int>resultID))

def bv_slice( S ,  bv ,  aID ,  lower ,  upper ,  resultID ):
    """Cython signature: void bv_slice(void* S, void* bv, int aID, int lower, int upper, int resultID)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(bv, (int, long)), 'arg bv wrong type'
    assert isinstance(aID, (int, long)), 'arg aID wrong type'
    assert isinstance(lower, (int, long)), 'arg lower wrong type'
    assert isinstance(upper, (int, long)), 'arg upper wrong type'
    assert isinstance(resultID, (int, long)), 'arg resultID wrong type'






    _bv_slice_monosat((<void*>S), (<void*>bv), (<int>aID), (<int>lower), (<int>upper), (<int>resultID))

def bv_subtraction( S ,  bv ,  bvID1 ,  bvID2 ,  resultID ):
    """Cython signature: void bv_subtraction(void* S, void* bv, int bvID1, int bvID2, int resultID)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(bv, (int, long)), 'arg bv wrong type'
    assert isinstance(bvID1, (int, long)), 'arg bvID1 wrong type'
    assert isinstance(bvID2, (int, long)), 'arg bvID2 wrong type'
    assert isinstance(resultID, (int, long)), 'arg resultID wrong type'





    _bv_subtraction_monosat((<void*>S), (<void*>bv), (<int>bvID1), (<int>bvID2), (<int>resultID))

def bv_unary( S ,  bv ,  args ,  n_args ,  resultID ):
    """Cython signature: void bv_unary(void* S, void* bv, int64_t args, int n_args, int resultID)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(bv, (int, long)), 'arg bv wrong type'
    assert isinstance(args, (int, long)), 'arg args wrong type'
    assert isinstance(n_args, (int, long)), 'arg n_args wrong type'
    assert isinstance(resultID, (int, long)), 'arg resultID wrong type'





    _bv_unary_monosat((<void*>S), (<void*>bv), (<int64_t>args), (<int>n_args), (<int>resultID))

def bv_width( S ,  bv ,  bvID ):
    """Cython signature: int bv_width(void* S, void* bv, int bvID)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(bv, (int, long)), 'arg bv wrong type'
    assert isinstance(bvID, (int, long)), 'arg bvID wrong type'



    cdef int _r = _bv_width_monosat((<void*>S), (<void*>bv), (<int>bvID))
    py_result = <int>_r
    return py_result

def bv_xnor( S ,  bv ,  bvaID ,  bvbID ,  bvResultID ):
    """Cython signature: void bv_xnor(void* S, void* bv, int bvaID, int bvbID, int bvResultID)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(bv, (int, long)), 'arg bv wrong type'
    assert isinstance(bvaID, (int, long)), 'arg bvaID wrong type'
    assert isinstance(bvbID, (int, long)), 'arg bvbID wrong type'
    assert isinstance(bvResultID, (int, long)), 'arg bvResultID wrong type'





    _bv_xnor_monosat((<void*>S), (<void*>bv), (<int>bvaID), (<int>bvbID), (<int>bvResultID))

def bv_xor( S ,  bv ,  bvaID ,  bvbID ,  bvResultID ):
    """Cython signature: void bv_xor(void* S, void* bv, int bvaID, int bvbID, int bvResultID)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(bv, (int, long)), 'arg bv wrong type'
    assert isinstance(bvaID, (int, long)), 'arg bvaID wrong type'
    assert isinstance(bvbID, (int, long)), 'arg bvbID wrong type'
    assert isinstance(bvResultID, (int, long)), 'arg bvResultID wrong type'





    _bv_xor_monosat((<void*>S), (<void*>bv), (<int>bvaID), (<int>bvbID), (<int>bvResultID))

def clearOptimizationObjectives( S ):
    """Cython signature: void clearOptimizationObjectives(void* S)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'

    _clearOptimizationObjectives_monosat((<void*>S))

def createFlowRouting( S ,  G ,  sourceNode ,  destNode ,  maxflowLit ):
    """Cython signature: void* createFlowRouting(void* S, void* G, int sourceNode, int destNode, int maxflowLit)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(G, (int, long)), 'arg G wrong type'
    assert isinstance(sourceNode, (int, long)), 'arg sourceNode wrong type'
    assert isinstance(destNode, (int, long)), 'arg destNode wrong type'
    assert isinstance(maxflowLit, (int, long)), 'arg maxflowLit wrong type'





    cdef void* _r = _createFlowRouting_monosat((<void*>S), (<void*>G), (<int>sourceNode), (<int>destNode), (<int>maxflowLit))
    py_result = <object>_r
    return py_result

def deleteSolver( S ):
    """Cython signature: void deleteSolver(void* S)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'

    _deleteSolver_monosat((<void*>S))

def disablePreprocessing( S ):
    """Cython signature: void disablePreprocessing(void* S)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'

    _disablePreprocessing_monosat((<void*>S))

def disallowLiteralSimplification( S ,  lit ):
    """Cython signature: bool disallowLiteralSimplification(void* S, int lit)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(lit, (int, long)), 'arg lit wrong type'


    cdef bool _r = _disallowLiteralSimplification_monosat((<void*>S), (<int>lit))
    py_result = <bool>_r
    return py_result

def flushPB( S ):
    """Cython signature: void flushPB(void* S)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'

    _flushPB_monosat((<void*>S))

def fsmAcceptsString( S ,  fsmTheory ,  fsmID ,  startNode ,  acceptNode ,  stringID ):
    """Cython signature: int fsmAcceptsString(void* S, void* fsmTheory, int fsmID, int startNode, int acceptNode, int stringID)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(fsmTheory, (int, long)), 'arg fsmTheory wrong type'
    assert isinstance(fsmID, (int, long)), 'arg fsmID wrong type'
    assert isinstance(startNode, (int, long)), 'arg startNode wrong type'
    assert isinstance(acceptNode, (int, long)), 'arg acceptNode wrong type'
    assert isinstance(stringID, (int, long)), 'arg stringID wrong type'






    cdef int _r = _fsmAcceptsString_monosat((<void*>S), (<void*>fsmTheory), (<int>fsmID), (<int>startNode), (<int>acceptNode), (<int>stringID))
    py_result = <int>_r
    return py_result

def fsmCompositionAccepts( S ,  fsmTheory ,  fsmGeneratorID ,  fsmAcceptorID ,  gen_startNode ,  gen_acceptNode ,  acceptor_startNode ,  acceptor_acceptNode ,  stringID ):
    """Cython signature: int fsmCompositionAccepts(void* S, void* fsmTheory, int fsmGeneratorID, int fsmAcceptorID, int gen_startNode, int gen_acceptNode, int acceptor_startNode, int acceptor_acceptNode, int stringID)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(fsmTheory, (int, long)), 'arg fsmTheory wrong type'
    assert isinstance(fsmGeneratorID, (int, long)), 'arg fsmGeneratorID wrong type'
    assert isinstance(fsmAcceptorID, (int, long)), 'arg fsmAcceptorID wrong type'
    assert isinstance(gen_startNode, (int, long)), 'arg gen_startNode wrong type'
    assert isinstance(gen_acceptNode, (int, long)), 'arg gen_acceptNode wrong type'
    assert isinstance(acceptor_startNode, (int, long)), 'arg acceptor_startNode wrong type'
    assert isinstance(acceptor_acceptNode, (int, long)), 'arg acceptor_acceptNode wrong type'
    assert isinstance(stringID, (int, long)), 'arg stringID wrong type'









    cdef int _r = _fsmCompositionAccepts_monosat((<void*>S), (<void*>fsmTheory), (<int>fsmGeneratorID), (<int>fsmAcceptorID), (<int>gen_startNode), (<int>gen_acceptNode), (<int>acceptor_startNode), (<int>acceptor_acceptNode), (<int>stringID))
    py_result = <int>_r
    return py_result

def getConflictClause( S ,  store_clause ,  max_store_size ):
    """Cython signature: int getConflictClause(void* S, int64_t store_clause, int max_store_size)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(store_clause, (int, long)), 'arg store_clause wrong type'
    assert isinstance(max_store_size, (int, long)), 'arg max_store_size wrong type'



    cdef int _r = _getConflictClause_monosat((<void*>S), (<int64_t>store_clause), (<int>max_store_size))
    py_result = <int>_r
    return py_result

def getDecisionPolarity( S ,  v ):
    """Cython signature: bool getDecisionPolarity(void* S, int v)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(v, (int, long)), 'arg v wrong type'


    cdef bool _r = _getDecisionPolarity_monosat((<void*>S), (<int>v))
    py_result = <bool>_r
    return py_result

def getDecisionPriority( S ,  var ):
    """Cython signature: int getDecisionPriority(void* S, int var)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(var, (int, long)), 'arg var wrong type'


    cdef int _r = _getDecisionPriority_monosat((<void*>S), (<int>var))
    py_result = <int>_r
    return py_result

def getModel_AcyclicEdgeFlow( S ,  G ,  maxflow_literal ,  edgeLit ):
    """Cython signature: int64_t getModel_AcyclicEdgeFlow(void* S, void* G, int maxflow_literal, int edgeLit)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(G, (int, long)), 'arg G wrong type'
    assert isinstance(maxflow_literal, (int, long)), 'arg maxflow_literal wrong type'
    assert isinstance(edgeLit, (int, long)), 'arg edgeLit wrong type'




    cdef int64_t _r = _getModel_AcyclicEdgeFlow_monosat((<void*>S), (<void*>G), (<int>maxflow_literal), (<int>edgeLit))
    py_result = <int64_t>_r
    return py_result

def getModel_BV( S ,  bv ,  bvID ,  getMaximumValue ):
    """Cython signature: int64_t getModel_BV(void* S, void* bv, int bvID, bool getMaximumValue)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(bv, (int, long)), 'arg bv wrong type'
    assert isinstance(bvID, (int, long)), 'arg bvID wrong type'
    assert isinstance(getMaximumValue, (int, long)), 'arg getMaximumValue wrong type'




    cdef int64_t _r = _getModel_BV_monosat((<void*>S), (<void*>bv), (<int>bvID), (<bool>getMaximumValue))
    py_result = <int64_t>_r
    return py_result

def getModel_EdgeFlow( S ,  G ,  maxflow_literal ,  edgeLit ):
    """Cython signature: int64_t getModel_EdgeFlow(void* S, void* G, int maxflow_literal, int edgeLit)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(G, (int, long)), 'arg G wrong type'
    assert isinstance(maxflow_literal, (int, long)), 'arg maxflow_literal wrong type'
    assert isinstance(edgeLit, (int, long)), 'arg edgeLit wrong type'




    cdef int64_t _r = _getModel_EdgeFlow_monosat((<void*>S), (<void*>G), (<int>maxflow_literal), (<int>edgeLit))
    py_result = <int64_t>_r
    return py_result

def getModel_Literal( S ,  lit ):
    """Cython signature: int getModel_Literal(void* S, int lit)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(lit, (int, long)), 'arg lit wrong type'


    cdef int _r = _getModel_Literal_monosat((<void*>S), (<int>lit))
    py_result = <int>_r
    return py_result

def getModel_MaxFlow( S ,  G ,  maxflow_literal ):
    """Cython signature: int64_t getModel_MaxFlow(void* S, void* G, int maxflow_literal)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(G, (int, long)), 'arg G wrong type'
    assert isinstance(maxflow_literal, (int, long)), 'arg maxflow_literal wrong type'



    cdef int64_t _r = _getModel_MaxFlow_monosat((<void*>S), (<void*>G), (<int>maxflow_literal))
    py_result = <int64_t>_r
    return py_result

def getModel_MinimumSpanningTreeWeight( S ,  G ,  spanning_tree_literal ):
    """Cython signature: int64_t getModel_MinimumSpanningTreeWeight(void* S, void* G, int spanning_tree_literal)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(G, (int, long)), 'arg G wrong type'
    assert isinstance(spanning_tree_literal, (int, long)), 'arg spanning_tree_literal wrong type'



    cdef int64_t _r = _getModel_MinimumSpanningTreeWeight_monosat((<void*>S), (<void*>G), (<int>spanning_tree_literal))
    py_result = <int64_t>_r
    return py_result

def getModel_Path_EdgeLits( S ,  G ,  reach_or_distance_literal ,  store_length ,  store ):
    """Cython signature: int getModel_Path_EdgeLits(void* S, void* G, int reach_or_distance_literal, int store_length, int64_t store)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(G, (int, long)), 'arg G wrong type'
    assert isinstance(reach_or_distance_literal, (int, long)), 'arg reach_or_distance_literal wrong type'
    assert isinstance(store_length, (int, long)), 'arg store_length wrong type'
    assert isinstance(store, (int, long)), 'arg store wrong type'





    cdef int _r = _getModel_Path_EdgeLits_monosat((<void*>S), (<void*>G), (<int>reach_or_distance_literal), (<int>store_length), (<int64_t>store))
    py_result = <int>_r
    return py_result

def getModel_Path_EdgeLits_Length( S ,  G ,  reach_or_distance_literal ):
    """Cython signature: int getModel_Path_EdgeLits_Length(void* S, void* G, int reach_or_distance_literal)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(G, (int, long)), 'arg G wrong type'
    assert isinstance(reach_or_distance_literal, (int, long)), 'arg reach_or_distance_literal wrong type'



    cdef int _r = _getModel_Path_EdgeLits_Length_monosat((<void*>S), (<void*>G), (<int>reach_or_distance_literal))
    py_result = <int>_r
    return py_result

def getModel_Path_Nodes( S ,  G ,  reach_or_distance_literal ,  store_length ,  store ):
    """Cython signature: int getModel_Path_Nodes(void* S, void* G, int reach_or_distance_literal, int store_length, int64_t store)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(G, (int, long)), 'arg G wrong type'
    assert isinstance(reach_or_distance_literal, (int, long)), 'arg reach_or_distance_literal wrong type'
    assert isinstance(store_length, (int, long)), 'arg store_length wrong type'
    assert isinstance(store, (int, long)), 'arg store wrong type'





    cdef int _r = _getModel_Path_Nodes_monosat((<void*>S), (<void*>G), (<int>reach_or_distance_literal), (<int>store_length), (<int64_t>store))
    py_result = <int>_r
    return py_result

def getModel_Path_Nodes_Length( S ,  G ,  reach_or_distance_literal ):
    """Cython signature: int getModel_Path_Nodes_Length(void* S, void* G, int reach_or_distance_literal)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(G, (int, long)), 'arg G wrong type'
    assert isinstance(reach_or_distance_literal, (int, long)), 'arg reach_or_distance_literal wrong type'



    cdef int _r = _getModel_Path_Nodes_Length_monosat((<void*>S), (<void*>G), (<int>reach_or_distance_literal))
    py_result = <int>_r
    return py_result

def getVersion():
    """Cython signature: char * getVersion()"""
    cdef char  * _r = <char *>_getVersion_monosat()
    py_result = <char *>(_r)
    return py_result

def graph_setAssignEdgesToWeight( S ,  G ,  weight ):
    """Cython signature: void graph_setAssignEdgesToWeight(void* S, void* G, int64_t weight)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(G, (int, long)), 'arg G wrong type'
    assert isinstance(weight, (int, long)), 'arg weight wrong type'



    _graph_setAssignEdgesToWeight_monosat((<void*>S), (<void*>G), (<int64_t>weight))

def initBVTheory( S ):
    """Cython signature: void* initBVTheory(void* S)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'

    cdef void* _r = _initBVTheory_monosat((<void*>S))
    py_result = <object>_r
    return py_result

def initFSMTheory( S ):
    """Cython signature: void* initFSMTheory(void* S)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'

    cdef void* _r = _initFSMTheory_monosat((<void*>S))
    py_result = <object>_r
    return py_result

def isDecisionVar( S ,  var ):
    """Cython signature: bool isDecisionVar(void* S, int var)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(var, (int, long)), 'arg var wrong type'


    cdef bool _r = _isDecisionVar_monosat((<void*>S), (<int>var))
    py_result = <bool>_r
    return py_result

def lastSolutionWasOptimal( S ):
    """Cython signature: bool lastSolutionWasOptimal(void* S)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'

    cdef bool _r = _lastSolutionWasOptimal_monosat((<void*>S))
    py_result = <bool>_r
    return py_result

def maximizeBV( S ,  bv ,  bvID ):
    """Cython signature: void maximizeBV(void* S, void* bv, int bvID)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(bv, (int, long)), 'arg bv wrong type'
    assert isinstance(bvID, (int, long)), 'arg bvID wrong type'



    _maximizeBV_monosat((<void*>S), (<void*>bv), (<int>bvID))

def maximizeLits( S ,  lits ,  n_lits ):
    """Cython signature: void maximizeLits(void* S, int64_t lits, int n_lits)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(lits, (int, long)), 'arg lits wrong type'
    assert isinstance(n_lits, (int, long)), 'arg n_lits wrong type'



    _maximizeLits_monosat((<void*>S), (<int64_t>lits), (<int>n_lits))

def maximizeWeightedLits( S ,  lits ,  weights ,  n_lits ):
    """Cython signature: void maximizeWeightedLits(void* S, int64_t lits, int64_t weights, int n_lits)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(lits, (int, long)), 'arg lits wrong type'
    assert isinstance(weights, (int, long)), 'arg weights wrong type'
    assert isinstance(n_lits, (int, long)), 'arg n_lits wrong type'




    _maximizeWeightedLits_monosat((<void*>S), (<int64_t>lits), (<int64_t>weights), (<int>n_lits))

def maximumFlow_geq( S ,  G ,  source ,  sink ,  weight ):
    """Cython signature: int maximumFlow_geq(void* S, void* G, int source, int sink, int64_t weight)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(G, (int, long)), 'arg G wrong type'
    assert isinstance(source, (int, long)), 'arg source wrong type'
    assert isinstance(sink, (int, long)), 'arg sink wrong type'
    assert isinstance(weight, (int, long)), 'arg weight wrong type'





    cdef int _r = _maximumFlow_geq_monosat((<void*>S), (<void*>G), (<int>source), (<int>sink), (<int64_t>weight))
    py_result = <int>_r
    return py_result

def maximumFlow_geq_bv( S ,  G ,  source ,  sink ,  bvID ):
    """Cython signature: int maximumFlow_geq_bv(void* S, void* G, int source, int sink, int bvID)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(G, (int, long)), 'arg G wrong type'
    assert isinstance(source, (int, long)), 'arg source wrong type'
    assert isinstance(sink, (int, long)), 'arg sink wrong type'
    assert isinstance(bvID, (int, long)), 'arg bvID wrong type'





    cdef int _r = _maximumFlow_geq_bv_monosat((<void*>S), (<void*>G), (<int>source), (<int>sink), (<int>bvID))
    py_result = <int>_r
    return py_result

def maximumFlow_gt( S ,  G ,  source ,  sink ,  weight ):
    """Cython signature: int maximumFlow_gt(void* S, void* G, int source, int sink, int64_t weight)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(G, (int, long)), 'arg G wrong type'
    assert isinstance(source, (int, long)), 'arg source wrong type'
    assert isinstance(sink, (int, long)), 'arg sink wrong type'
    assert isinstance(weight, (int, long)), 'arg weight wrong type'





    cdef int _r = _maximumFlow_gt_monosat((<void*>S), (<void*>G), (<int>source), (<int>sink), (<int64_t>weight))
    py_result = <int>_r
    return py_result

def maximumFlow_gt_bv( S ,  G ,  source ,  sink ,  bvID ):
    """Cython signature: int maximumFlow_gt_bv(void* S, void* G, int source, int sink, int bvID)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(G, (int, long)), 'arg G wrong type'
    assert isinstance(source, (int, long)), 'arg source wrong type'
    assert isinstance(sink, (int, long)), 'arg sink wrong type'
    assert isinstance(bvID, (int, long)), 'arg bvID wrong type'





    cdef int _r = _maximumFlow_gt_bv_monosat((<void*>S), (<void*>G), (<int>source), (<int>sink), (<int>bvID))
    py_result = <int>_r
    return py_result

def minimizeBV( S ,  bv ,  bvID ):
    """Cython signature: void minimizeBV(void* S, void* bv, int bvID)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(bv, (int, long)), 'arg bv wrong type'
    assert isinstance(bvID, (int, long)), 'arg bvID wrong type'



    _minimizeBV_monosat((<void*>S), (<void*>bv), (<int>bvID))

def minimizeLits( S ,  lits ,  n_lits ):
    """Cython signature: void minimizeLits(void* S, int64_t lits, int n_lits)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(lits, (int, long)), 'arg lits wrong type'
    assert isinstance(n_lits, (int, long)), 'arg n_lits wrong type'



    _minimizeLits_monosat((<void*>S), (<int64_t>lits), (<int>n_lits))

def minimizeWeightedLits( S ,  lits ,  weights ,  n_lits ):
    """Cython signature: void minimizeWeightedLits(void* S, int64_t lits, int64_t weights, int n_lits)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(lits, (int, long)), 'arg lits wrong type'
    assert isinstance(weights, (int, long)), 'arg weights wrong type'
    assert isinstance(n_lits, (int, long)), 'arg n_lits wrong type'




    _minimizeWeightedLits_monosat((<void*>S), (<int64_t>lits), (<int64_t>weights), (<int>n_lits))

def minimumSpanningTree_leq( S ,  G ,  weight ):
    """Cython signature: int minimumSpanningTree_leq(void* S, void* G, int64_t weight)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(G, (int, long)), 'arg G wrong type'
    assert isinstance(weight, (int, long)), 'arg weight wrong type'



    cdef int _r = _minimumSpanningTree_leq_monosat((<void*>S), (<void*>G), (<int64_t>weight))
    py_result = <int>_r
    return py_result

def minimumSpanningTree_lt( S ,  G ,  source ,  sink ,  weight ):
    """Cython signature: int minimumSpanningTree_lt(void* S, void* G, int source, int sink, int64_t weight)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(G, (int, long)), 'arg G wrong type'
    assert isinstance(source, (int, long)), 'arg source wrong type'
    assert isinstance(sink, (int, long)), 'arg sink wrong type'
    assert isinstance(weight, (int, long)), 'arg weight wrong type'





    cdef int _r = _minimumSpanningTree_lt_monosat((<void*>S), (<void*>G), (<int>source), (<int>sink), (<int64_t>weight))
    py_result = <int>_r
    return py_result

def nBitvectors( S ,  bv ):
    """Cython signature: int nBitvectors(void* S, void* bv)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(bv, (int, long)), 'arg bv wrong type'


    cdef int _r = _nBitvectors_monosat((<void*>S), (<void*>bv))
    py_result = <int>_r
    return py_result

def nClauses( S ):
    """Cython signature: int nClauses(void* S)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'

    cdef int _r = _nClauses_monosat((<void*>S))
    py_result = <int>_r
    return py_result

def nEdges( S ,  G ):
    """Cython signature: int nEdges(void* S, void* G)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(G, (int, long)), 'arg G wrong type'


    cdef int _r = _nEdges_monosat((<void*>S), (<void*>G))
    py_result = <int>_r
    return py_result

def nNodes( S ,  G ):
    """Cython signature: int nNodes(void* S, void* G)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(G, (int, long)), 'arg G wrong type'


    cdef int _r = _nNodes_monosat((<void*>S), (<void*>G))
    py_result = <int>_r
    return py_result

def nVars( S ):
    """Cython signature: int nVars(void* S)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'

    cdef int _r = _nVars_monosat((<void*>S))
    py_result = <int>_r
    return py_result

def newBVComparison_bv_geq( S ,  bv ,  bvID ,  compareID ):
    """Cython signature: int newBVComparison_bv_geq(void* S, void* bv, int bvID, int compareID)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(bv, (int, long)), 'arg bv wrong type'
    assert isinstance(bvID, (int, long)), 'arg bvID wrong type'
    assert isinstance(compareID, (int, long)), 'arg compareID wrong type'




    cdef int _r = _newBVComparison_bv_geq_monosat((<void*>S), (<void*>bv), (<int>bvID), (<int>compareID))
    py_result = <int>_r
    return py_result

def newBVComparison_bv_gt( S ,  bv ,  bvID ,  compareID ):
    """Cython signature: int newBVComparison_bv_gt(void* S, void* bv, int bvID, int compareID)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(bv, (int, long)), 'arg bv wrong type'
    assert isinstance(bvID, (int, long)), 'arg bvID wrong type'
    assert isinstance(compareID, (int, long)), 'arg compareID wrong type'




    cdef int _r = _newBVComparison_bv_gt_monosat((<void*>S), (<void*>bv), (<int>bvID), (<int>compareID))
    py_result = <int>_r
    return py_result

def newBVComparison_bv_leq( S ,  bv ,  bvID ,  compareID ):
    """Cython signature: int newBVComparison_bv_leq(void* S, void* bv, int bvID, int compareID)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(bv, (int, long)), 'arg bv wrong type'
    assert isinstance(bvID, (int, long)), 'arg bvID wrong type'
    assert isinstance(compareID, (int, long)), 'arg compareID wrong type'




    cdef int _r = _newBVComparison_bv_leq_monosat((<void*>S), (<void*>bv), (<int>bvID), (<int>compareID))
    py_result = <int>_r
    return py_result

def newBVComparison_bv_lt( S ,  bv ,  bvID ,  compareID ):
    """Cython signature: int newBVComparison_bv_lt(void* S, void* bv, int bvID, int compareID)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(bv, (int, long)), 'arg bv wrong type'
    assert isinstance(bvID, (int, long)), 'arg bvID wrong type'
    assert isinstance(compareID, (int, long)), 'arg compareID wrong type'




    cdef int _r = _newBVComparison_bv_lt_monosat((<void*>S), (<void*>bv), (<int>bvID), (<int>compareID))
    py_result = <int>_r
    return py_result

def newBVComparison_const_geq( S ,  bv ,  bvID ,  weight ):
    """Cython signature: int newBVComparison_const_geq(void* S, void* bv, int bvID, int64_t weight)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(bv, (int, long)), 'arg bv wrong type'
    assert isinstance(bvID, (int, long)), 'arg bvID wrong type'
    assert isinstance(weight, (int, long)), 'arg weight wrong type'




    cdef int _r = _newBVComparison_const_geq_monosat((<void*>S), (<void*>bv), (<int>bvID), (<int64_t>weight))
    py_result = <int>_r
    return py_result

def newBVComparison_const_gt( S ,  bv ,  bvID ,  weight ):
    """Cython signature: int newBVComparison_const_gt(void* S, void* bv, int bvID, int64_t weight)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(bv, (int, long)), 'arg bv wrong type'
    assert isinstance(bvID, (int, long)), 'arg bvID wrong type'
    assert isinstance(weight, (int, long)), 'arg weight wrong type'




    cdef int _r = _newBVComparison_const_gt_monosat((<void*>S), (<void*>bv), (<int>bvID), (<int64_t>weight))
    py_result = <int>_r
    return py_result

def newBVComparison_const_leq( S ,  bv ,  bvID ,  weight ):
    """Cython signature: int newBVComparison_const_leq(void* S, void* bv, int bvID, int64_t weight)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(bv, (int, long)), 'arg bv wrong type'
    assert isinstance(bvID, (int, long)), 'arg bvID wrong type'
    assert isinstance(weight, (int, long)), 'arg weight wrong type'




    cdef int _r = _newBVComparison_const_leq_monosat((<void*>S), (<void*>bv), (<int>bvID), (<int64_t>weight))
    py_result = <int>_r
    return py_result

def newBVComparison_const_lt( S ,  bv ,  bvID ,  weight ):
    """Cython signature: int newBVComparison_const_lt(void* S, void* bv, int bvID, int64_t weight)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(bv, (int, long)), 'arg bv wrong type'
    assert isinstance(bvID, (int, long)), 'arg bvID wrong type'
    assert isinstance(weight, (int, long)), 'arg weight wrong type'




    cdef int _r = _newBVComparison_const_lt_monosat((<void*>S), (<void*>bv), (<int>bvID), (<int64_t>weight))
    py_result = <int>_r
    return py_result

def newBitvector( S ,  bv ,  bits ,  n_bits ):
    """Cython signature: int newBitvector(void* S, void* bv, int64_t bits, int n_bits)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(bv, (int, long)), 'arg bv wrong type'
    assert isinstance(bits, (int, long)), 'arg bits wrong type'
    assert isinstance(n_bits, (int, long)), 'arg n_bits wrong type'




    cdef int _r = _newBitvector_monosat((<void*>S), (<void*>bv), (<int64_t>bits), (<int>n_bits))
    py_result = <int>_r
    return py_result

def newBitvector_anon( S ,  bv ,  bvWidth ):
    """Cython signature: int newBitvector_anon(void* S, void* bv, int bvWidth)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(bv, (int, long)), 'arg bv wrong type'
    assert isinstance(bvWidth, (int, long)), 'arg bvWidth wrong type'



    cdef int _r = _newBitvector_anon_monosat((<void*>S), (<void*>bv), (<int>bvWidth))
    py_result = <int>_r
    return py_result

def newBitvector_const( S ,  bv ,  bvWidth ,  constval ):
    """Cython signature: int newBitvector_const(void* S, void* bv, int bvWidth, int64_t constval)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(bv, (int, long)), 'arg bv wrong type'
    assert isinstance(bvWidth, (int, long)), 'arg bvWidth wrong type'
    assert isinstance(constval, (int, long)), 'arg constval wrong type'




    cdef int _r = _newBitvector_const_monosat((<void*>S), (<void*>bv), (<int>bvWidth), (<int64_t>constval))
    py_result = <int>_r
    return py_result

def newEdge( S ,  G ,  _from ,  to ,  weight ):
    """Cython signature: int newEdge(void* S, void* G, int _from, int to, int64_t weight)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(G, (int, long)), 'arg G wrong type'
    assert isinstance(_from, (int, long)), 'arg _from wrong type'
    assert isinstance(to, (int, long)), 'arg to wrong type'
    assert isinstance(weight, (int, long)), 'arg weight wrong type'





    cdef int _r = _newEdge_monosat((<void*>S), (<void*>G), (<int>_from), (<int>to), (<int64_t>weight))
    py_result = <int>_r
    return py_result

def newEdgeSet( S ,  G ,  edges ,  n_edges ,  enforceEdgeAssignment ):
    """Cython signature: void newEdgeSet(void* S, void* G, int64_t edges, int n_edges, bool enforceEdgeAssignment)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(G, (int, long)), 'arg G wrong type'
    assert isinstance(edges, (int, long)), 'arg edges wrong type'
    assert isinstance(n_edges, (int, long)), 'arg n_edges wrong type'
    assert isinstance(enforceEdgeAssignment, (int, long)), 'arg enforceEdgeAssignment wrong type'





    _newEdgeSet_monosat((<void*>S), (<void*>G), (<int64_t>edges), (<int>n_edges), (<bool>enforceEdgeAssignment))

def newEdge_bv( S ,  G ,  _from ,  to ,  bvID ):
    """Cython signature: int newEdge_bv(void* S, void* G, int _from, int to, int bvID)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(G, (int, long)), 'arg G wrong type'
    assert isinstance(_from, (int, long)), 'arg _from wrong type'
    assert isinstance(to, (int, long)), 'arg to wrong type'
    assert isinstance(bvID, (int, long)), 'arg bvID wrong type'





    cdef int _r = _newEdge_bv_monosat((<void*>S), (<void*>G), (<int>_from), (<int>to), (<int>bvID))
    py_result = <int>_r
    return py_result

def newEdge_double( S ,  G ,  _from ,  to , double weight ):
    """Cython signature: int newEdge_double(void* S, void* G, int _from, int to, double weight)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(G, (int, long)), 'arg G wrong type'
    assert isinstance(_from, (int, long)), 'arg _from wrong type'
    assert isinstance(to, (int, long)), 'arg to wrong type'
    assert isinstance(weight, float), 'arg weight wrong type'





    cdef int _r = _newEdge_double_monosat((<void*>S), (<void*>G), (<int>_from), (<int>to), (<double>weight))
    py_result = <int>_r
    return py_result

def newFSM( S ,  fsmTheory ,  inputAlphabet ,  outputAlphabet ):
    """Cython signature: int newFSM(void* S, void* fsmTheory, int inputAlphabet, int outputAlphabet)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(fsmTheory, (int, long)), 'arg fsmTheory wrong type'
    assert isinstance(inputAlphabet, (int, long)), 'arg inputAlphabet wrong type'
    assert isinstance(outputAlphabet, (int, long)), 'arg outputAlphabet wrong type'




    cdef int _r = _newFSM_monosat((<void*>S), (<void*>fsmTheory), (<int>inputAlphabet), (<int>outputAlphabet))
    py_result = <int>_r
    return py_result

def newGraph( S ):
    """Cython signature: void* newGraph(void* S)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'

    cdef void* _r = _newGraph_monosat((<void*>S))
    py_result = <object>_r
    return py_result

def newNode( S ,  G ):
    """Cython signature: int newNode(void* S, void* G)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(G, (int, long)), 'arg G wrong type'


    cdef int _r = _newNode_monosat((<void*>S), (<void*>G))
    py_result = <int>_r
    return py_result

def newSolver():
    """Cython signature: void* newSolver()"""
    cdef void* _r = _newSolver_monosat()
    if _r is NULL:
        raise MemoryError()
    py_result = <object>_r
    return py_result

def newSolver_arg(bytes argv ):
    """Cython signature: void* newSolver_arg(char * argv)"""
    assert isinstance(argv, bytes), 'arg argv wrong type'

    cdef void* _r = _newSolver_arg_monosat((<char *>argv))
    if _r is NULL:
        raise MemoryError()
    py_result = <object>_r
    return py_result

def newState( S ,  fsmTheory ,  fsmID ):
    """Cython signature: int newState(void* S, void* fsmTheory, int fsmID)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(fsmTheory, (int, long)), 'arg fsmTheory wrong type'
    assert isinstance(fsmID, (int, long)), 'arg fsmID wrong type'



    cdef int _r = _newState_monosat((<void*>S), (<void*>fsmTheory), (<int>fsmID))
    py_result = <int>_r
    return py_result

def newString( S ,  fsmTheory ,  str ,  len ):
    """Cython signature: int newString(void* S, void* fsmTheory, int64_t str, int len)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(fsmTheory, (int, long)), 'arg fsmTheory wrong type'
    assert isinstance(str, (int, long)), 'arg str wrong type'
    assert isinstance(len, (int, long)), 'arg len wrong type'




    cdef int _r = _newString_monosat((<void*>S), (<void*>fsmTheory), (<int64_t>str), (<int>len))
    py_result = <int>_r
    return py_result

def newTransition( S ,  fsmTheory ,  fsmID ,  _fromNode ,  toNode ,  inputLabel ,  outputLabel ):
    """Cython signature: int newTransition(void* S, void* fsmTheory, int fsmID, int _fromNode, int toNode, int inputLabel, int outputLabel)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(fsmTheory, (int, long)), 'arg fsmTheory wrong type'
    assert isinstance(fsmID, (int, long)), 'arg fsmID wrong type'
    assert isinstance(_fromNode, (int, long)), 'arg _fromNode wrong type'
    assert isinstance(toNode, (int, long)), 'arg toNode wrong type'
    assert isinstance(inputLabel, (int, long)), 'arg inputLabel wrong type'
    assert isinstance(outputLabel, (int, long)), 'arg outputLabel wrong type'







    cdef int _r = _newTransition_monosat((<void*>S), (<void*>fsmTheory), (<int>fsmID), (<int>_fromNode), (<int>toNode), (<int>inputLabel), (<int>outputLabel))
    py_result = <int>_r
    return py_result

def newVar( S ):
    """Cython signature: int newVar(void* S)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'

    cdef int _r = _newVar_monosat((<void*>S))
    py_result = <int>_r
    return py_result

def reaches( S ,  G ,  _from ,  to ):
    """Cython signature: int reaches(void* S, void* G, int _from, int to)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(G, (int, long)), 'arg G wrong type'
    assert isinstance(_from, (int, long)), 'arg _from wrong type'
    assert isinstance(to, (int, long)), 'arg to wrong type'




    cdef int _r = _reaches_monosat((<void*>S), (<void*>G), (<int>_from), (<int>to))
    py_result = <int>_r
    return py_result

def readGNF( S , bytes filename ):
    """Cython signature: void readGNF(void* S, char * filename)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(filename, bytes), 'arg filename wrong type'


    _readGNF_monosat((<void*>S), (<char *>filename))

def setConflictLimit( S ,  num_conflicts ):
    """Cython signature: void setConflictLimit(void* S, int num_conflicts)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(num_conflicts, (int, long)), 'arg num_conflicts wrong type'


    _setConflictLimit_monosat((<void*>S), (<int>num_conflicts))

def setDecisionPolarity( S ,  v ,  b ):
    """Cython signature: void setDecisionPolarity(void* S, int v, bool b)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(v, (int, long)), 'arg v wrong type'
    assert isinstance(b, (int, long)), 'arg b wrong type'



    _setDecisionPolarity_monosat((<void*>S), (<int>v), (<bool>b))

def setDecisionPriority( S ,  var ,  priority ):
    """Cython signature: void setDecisionPriority(void* S, int var, int priority)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(var, (int, long)), 'arg var wrong type'
    assert isinstance(priority, (int, long)), 'arg priority wrong type'



    _setDecisionPriority_monosat((<void*>S), (<int>var), (<int>priority))

def setDecisionVar( S ,  var ,  decidable ):
    """Cython signature: void setDecisionVar(void* S, int var, bool decidable)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(var, (int, long)), 'arg var wrong type'
    assert isinstance(decidable, (int, long)), 'arg decidable wrong type'



    _setDecisionVar_monosat((<void*>S), (<int>var), (<bool>decidable))

def setMemoryLimit( S ,  mb ):
    """Cython signature: void setMemoryLimit(void* S, int mb)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(mb, (int, long)), 'arg mb wrong type'


    _setMemoryLimit_monosat((<void*>S), (<int>mb))

def setOutputFile( S , bytes output ):
    """Cython signature: void setOutputFile(void* S, char * output)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(output, bytes), 'arg output wrong type'


    _setOutputFile_monosat((<void*>S), (<char *>output))

def setPropagationLimit( S ,  num_propagations ):
    """Cython signature: void setPropagationLimit(void* S, int num_propagations)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(num_propagations, (int, long)), 'arg num_propagations wrong type'


    _setPropagationLimit_monosat((<void*>S), (<int>num_propagations))

def setTimeLimit( S ,  seconds ):
    """Cython signature: void setTimeLimit(void* S, int seconds)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(seconds, (int, long)), 'arg seconds wrong type'


    _setTimeLimit_monosat((<void*>S), (<int>seconds))

def shortestPathUnweighted_leq_const( S ,  G ,  _from ,  to ,  steps ):
    """Cython signature: int shortestPathUnweighted_leq_const(void* S, void* G, int _from, int to, int steps)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(G, (int, long)), 'arg G wrong type'
    assert isinstance(_from, (int, long)), 'arg _from wrong type'
    assert isinstance(to, (int, long)), 'arg to wrong type'
    assert isinstance(steps, (int, long)), 'arg steps wrong type'





    cdef int _r = _shortestPathUnweighted_leq_const_monosat((<void*>S), (<void*>G), (<int>_from), (<int>to), (<int>steps))
    py_result = <int>_r
    return py_result

def shortestPathUnweighted_lt_const( S ,  G ,  _from ,  to ,  steps ):
    """Cython signature: int shortestPathUnweighted_lt_const(void* S, void* G, int _from, int to, int steps)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(G, (int, long)), 'arg G wrong type'
    assert isinstance(_from, (int, long)), 'arg _from wrong type'
    assert isinstance(to, (int, long)), 'arg to wrong type'
    assert isinstance(steps, (int, long)), 'arg steps wrong type'





    cdef int _r = _shortestPathUnweighted_lt_const_monosat((<void*>S), (<void*>G), (<int>_from), (<int>to), (<int>steps))
    py_result = <int>_r
    return py_result

def shortestPath_leq_bv( S ,  G ,  _from ,  to ,  bvID ):
    """Cython signature: int shortestPath_leq_bv(void* S, void* G, int _from, int to, int bvID)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(G, (int, long)), 'arg G wrong type'
    assert isinstance(_from, (int, long)), 'arg _from wrong type'
    assert isinstance(to, (int, long)), 'arg to wrong type'
    assert isinstance(bvID, (int, long)), 'arg bvID wrong type'





    cdef int _r = _shortestPath_leq_bv_monosat((<void*>S), (<void*>G), (<int>_from), (<int>to), (<int>bvID))
    py_result = <int>_r
    return py_result

def shortestPath_leq_const( S ,  G ,  _from ,  to ,  dist ):
    """Cython signature: int shortestPath_leq_const(void* S, void* G, int _from, int to, int64_t dist)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(G, (int, long)), 'arg G wrong type'
    assert isinstance(_from, (int, long)), 'arg _from wrong type'
    assert isinstance(to, (int, long)), 'arg to wrong type'
    assert isinstance(dist, (int, long)), 'arg dist wrong type'





    cdef int _r = _shortestPath_leq_const_monosat((<void*>S), (<void*>G), (<int>_from), (<int>to), (<int64_t>dist))
    py_result = <int>_r
    return py_result

def shortestPath_lt_bv( S ,  G ,  _from ,  to ,  bvID ):
    """Cython signature: int shortestPath_lt_bv(void* S, void* G, int _from, int to, int bvID)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(G, (int, long)), 'arg G wrong type'
    assert isinstance(_from, (int, long)), 'arg _from wrong type'
    assert isinstance(to, (int, long)), 'arg to wrong type'
    assert isinstance(bvID, (int, long)), 'arg bvID wrong type'





    cdef int _r = _shortestPath_lt_bv_monosat((<void*>S), (<void*>G), (<int>_from), (<int>to), (<int>bvID))
    py_result = <int>_r
    return py_result

def shortestPath_lt_const( S ,  G ,  _from ,  to ,  dist ):
    """Cython signature: int shortestPath_lt_const(void* S, void* G, int _from, int to, int64_t dist)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(G, (int, long)), 'arg G wrong type'
    assert isinstance(_from, (int, long)), 'arg _from wrong type'
    assert isinstance(to, (int, long)), 'arg to wrong type'
    assert isinstance(dist, (int, long)), 'arg dist wrong type'





    cdef int _r = _shortestPath_lt_const_monosat((<void*>S), (<void*>G), (<int>_from), (<int>to), (<int64_t>dist))
    py_result = <int>_r
    return py_result

def solve( S ):
    """Cython signature: bool solve(void* S)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'

    cdef bool _r = _solve_monosat((<void*>S))
    py_result = <bool>_r
    return py_result

def solveAssumptions( S ,  assumptions ,  n_assumptions ):
    """Cython signature: bool solveAssumptions(void* S, int64_t assumptions, int n_assumptions)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(assumptions, (int, long)), 'arg assumptions wrong type'
    assert isinstance(n_assumptions, (int, long)), 'arg n_assumptions wrong type'



    cdef bool _r = _solveAssumptions_monosat((<void*>S), (<int64_t>assumptions), (<int>n_assumptions))
    py_result = <bool>_r
    return py_result

def solveAssumptionsLimited( S ,  assumptions ,  n_assumptions ):
    """Cython signature: int solveAssumptionsLimited(void* S, int64_t assumptions, int n_assumptions)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'
    assert isinstance(assumptions, (int, long)), 'arg assumptions wrong type'
    assert isinstance(n_assumptions, (int, long)), 'arg n_assumptions wrong type'



    cdef int _r = _solveAssumptionsLimited_monosat((<void*>S), (<int64_t>assumptions), (<int>n_assumptions))
    py_result = <int>_r
    return py_result

def solveLimited( S ):
    """Cython signature: int solveLimited(void* S)"""
    assert isinstance(S, (int, long)), 'arg S wrong type'

    cdef int _r = _solveLimited_monosat((<void*>S))
    py_result = <int>_r
    return py_result

def true_lit( S ):
    """Cython signature: int true_lit(void* S)"""


    cdef int _r = _true_lit_monosat((<void*>S))
    py_result = <int>_r
    return py_result 
