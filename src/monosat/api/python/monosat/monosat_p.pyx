#cython: c_string_encoding=ascii  # for cython>=0.19
#cython: embedsignature=False
from cpython cimport array
import array
cimport cpython.pycapsule as pycapsule
from libc.stdint cimport  int64_t
from libc.stdint cimport  uint64_t
from libc.stdint cimport  uintptr_t
from  libc.string cimport const_char
from cython.operator cimport dereference as deref, preincrement as inc, address as address


from monosat_header cimport BVTheoryPtr

from monosat_header cimport FlowRouterPtr
from monosat_header cimport GraphTheorySolver_double
from monosat_header cimport GraphTheorySolver_long
from monosat_header cimport SolverPtr
from monosat_header cimport Var
from monosat_header cimport Weight
from monosat_header cimport acyclic_directed as _acyclic_directed_monosat
from monosat_header cimport acyclic_undirected as _acyclic_undirected_monosat
from monosat_header cimport addBinaryClause as _addBinaryClause_monosat
from monosat_header cimport addBinaryClauses as _addBinaryClauses_monosat
from monosat_header cimport addClause as _addClause_monosat
from monosat_header cimport addRoutingNet as _addRoutingNet_monosat
from monosat_header cimport addTertiaryClause as _addTertiaryClause_monosat
from monosat_header cimport addUnitClause as _addUnitClause_monosat
from monosat_header cimport assertPB_eq as _assertPB_eq_monosat
from monosat_header cimport assertPB_geq as _assertPB_geq_monosat
from monosat_header cimport assertPB_gt as _assertPB_gt_monosat
from monosat_header cimport assertPB_leq as _assertPB_leq_monosat
from monosat_header cimport assertPB_lt as _assertPB_lt_monosat
from monosat_header cimport at_most_one as _at_most_one_monosat
from monosat_header cimport backtrack as _backtrack_monosat
from monosat_header cimport bv_addition as _bv_addition_monosat
from monosat_header cimport bv_and as _bv_and_monosat
from monosat_header cimport bv_bitblast as _bv_bitblast_monosat
from monosat_header cimport bv_concat as _bv_concat_monosat
from monosat_header cimport bv_divide as _bv_divide_monosat
from monosat_header cimport bv_ite as _bv_ite_monosat
from monosat_header cimport bv_max as _bv_max_monosat
from monosat_header cimport bv_min as _bv_min_monosat
from monosat_header cimport bv_multiply as _bv_multiply_monosat
from monosat_header cimport bv_nand as _bv_nand_monosat
from monosat_header cimport bv_nor as _bv_nor_monosat
from monosat_header cimport bv_not as _bv_not_monosat
from monosat_header cimport bv_or as _bv_or_monosat
from monosat_header cimport bv_popcount as _bv_popcount_monosat
from monosat_header cimport bv_slice as _bv_slice_monosat
from monosat_header cimport bv_subtraction as _bv_subtraction_monosat
from monosat_header cimport bv_unary as _bv_unary_monosat
from monosat_header cimport bv_width as _bv_width_monosat
from monosat_header cimport bv_xnor as _bv_xnor_monosat
from monosat_header cimport bv_xor as _bv_xor_monosat
from monosat_header cimport clearOptimizationObjectives as _clearOptimizationObjectives_monosat
from monosat_header cimport createFlowRouting as _createFlowRouting_monosat
from monosat_header cimport deleteSolver as _deleteSolver_monosat
from monosat_header cimport disablePreprocessing as _disablePreprocessing_monosat
from monosat_header cimport disallowLiteralSimplification as _disallowLiteralSimplification_monosat
from monosat_header cimport flushPB as _flushPB_monosat

from monosat_header cimport getConflictClause as _getConflictClause_monosat
from monosat_header cimport minimizeUnsatCore as _minimizeUnsatCore_monosat
from monosat_header cimport minimizeConflictClause as _minimizeConflictClause_monosat

from monosat_header cimport getDecisionPolarity as _getDecisionPolarity_monosat
from monosat_header cimport getDecisionPriority as _getDecisionPriority_monosat
from monosat_header cimport getModel_AcyclicEdgeFlow as _getModel_AcyclicEdgeFlow_monosat
from monosat_header cimport getModel_BV as _getModel_BV_monosat
from monosat_header cimport getModel_EdgeFlow as _getModel_EdgeFlow_monosat
from monosat_header cimport getModel_Literal as _getModel_Literal_monosat
from monosat_header cimport getModel_MaxFlow as _getModel_MaxFlow_monosat
from monosat_header cimport getModel_MinimumSpanningTreeWeight as _getModel_MinimumSpanningTreeWeight_monosat
from monosat_header cimport getModel_Path_EdgeLits as _getModel_Path_EdgeLits_monosat
from monosat_header cimport getModel_Path_EdgeLits_Length as _getModel_Path_EdgeLits_Length_monosat
from monosat_header cimport getModel_Path_Nodes as _getModel_Path_Nodes_monosat
from monosat_header cimport getModel_Path_Nodes_Length as _getModel_Path_Nodes_Length_monosat
from monosat_header cimport getVersion as _getVersion_monosat
from monosat_header cimport graph_setAssignEdgesToWeight as _graph_setAssignEdgesToWeight_monosat
from monosat_header cimport initBVTheory as _initBVTheory_monosat

from monosat_header cimport isDecisionVar as _isDecisionVar_monosat
from monosat_header cimport lastSolutionWasOptimal as _lastSolutionWasOptimal_monosat
from monosat_header cimport maximizeBV as _maximizeBV_monosat
from monosat_header cimport maximizeLits as _maximizeLits_monosat
from monosat_header cimport maximizeWeightedLits as _maximizeWeightedLits_monosat
from monosat_header cimport maximumFlow_geq as _maximumFlow_geq_monosat
from monosat_header cimport maximumFlow_geq_bv as _maximumFlow_geq_bv_monosat
from monosat_header cimport maximumFlow_gt as _maximumFlow_gt_monosat
from monosat_header cimport maximumFlow_gt_bv as _maximumFlow_gt_bv_monosat
from monosat_header cimport minimizeBV as _minimizeBV_monosat
from monosat_header cimport minimizeLits as _minimizeLits_monosat
from monosat_header cimport minimizeWeightedLits as _minimizeWeightedLits_monosat
from monosat_header cimport minimumSpanningTree_leq as _minimumSpanningTree_leq_monosat
from monosat_header cimport minimumSpanningTree_lt as _minimumSpanningTree_lt_monosat
from monosat_header cimport nBitvectors as _nBitvectors_monosat
from monosat_header cimport nClauses as _nClauses_monosat
from monosat_header cimport nEdges as _nEdges_monosat
from monosat_header cimport nNodes as _nNodes_monosat
from monosat_header cimport nVars as _nVars_monosat
from monosat_header cimport newBVComparison_bv_geq as _newBVComparison_bv_geq_monosat
from monosat_header cimport newBVComparison_bv_gt as _newBVComparison_bv_gt_monosat
from monosat_header cimport newBVComparison_bv_leq as _newBVComparison_bv_leq_monosat
from monosat_header cimport newBVComparison_bv_lt as _newBVComparison_bv_lt_monosat
from monosat_header cimport newBVComparison_const_geq as _newBVComparison_const_geq_monosat
from monosat_header cimport newBVComparison_const_gt as _newBVComparison_const_gt_monosat
from monosat_header cimport newBVComparison_const_leq as _newBVComparison_const_leq_monosat
from monosat_header cimport newBVComparison_const_lt as _newBVComparison_const_lt_monosat
from monosat_header cimport newBitvector as _newBitvector_monosat
from monosat_header cimport newBitvector_anon as _newBitvector_anon_monosat
from monosat_header cimport newBitvector_const as _newBitvector_const_monosat
from monosat_header cimport newEdge as _newEdge_monosat
from monosat_header cimport newEdgeSet as _newEdgeSet_monosat
from monosat_header cimport newEdge_bv as _newEdge_bv_monosat
from monosat_header cimport newEdge_double as _newEdge_double_monosat

from monosat_header cimport newGraph as _newGraph_monosat
from monosat_header cimport newNode as _newNode_monosat
from monosat_header cimport newSolver as _newSolver_monosat
from monosat_header cimport newSolver_arg as _newSolver_arg_monosat
from monosat_header cimport newState as _newState_monosat
from monosat_header cimport newString as _newString_monosat
from monosat_header cimport newTransition as _newTransition_monosat
from monosat_header cimport newVar as _newVar_monosat
from monosat_header cimport reaches as _reaches_monosat
from monosat_header cimport reachesBackward as _reachesBackward_monosat
from monosat_header cimport onPath as _onPath_monosat
from monosat_header cimport readGNF as _readGNF_monosat
from monosat_header cimport loadGNF as _loadGNF_monosat
from monosat_header cimport setConflictLimit as _setConflictLimit_monosat
from monosat_header cimport setDecisionPolarity as _setDecisionPolarity_monosat
from monosat_header cimport setDecisionPriority as _setDecisionPriority_monosat
from monosat_header cimport setDecisionVar as _setDecisionVar_monosat
from monosat_header cimport setOutputFile as _setOutputFile_monosat
from monosat_header cimport setPropagationLimit as _setPropagationLimit_monosat
from monosat_header cimport setTimeLimit as _setTimeLimit_monosat
from monosat_header cimport shortestPathUnweighted_leq_const as _shortestPathUnweighted_leq_const_monosat
from monosat_header cimport shortestPathUnweighted_lt_const as _shortestPathUnweighted_lt_const_monosat
from monosat_header cimport shortestPath_leq_bv as _shortestPath_leq_bv_monosat
from monosat_header cimport shortestPath_leq_const as _shortestPath_leq_const_monosat
from monosat_header cimport shortestPath_lt_bv as _shortestPath_lt_bv_monosat
from monosat_header cimport shortestPath_lt_const as _shortestPath_lt_const_monosat
from monosat_header cimport connectedComponents_geq_const as _connectedComponents_geq_const_monosat

from monosat_header cimport solve as _solve_monosat
from monosat_header cimport solveAssumptions as _solveAssumptions_monosat
from monosat_header cimport solveAssumptionsLimited as _solveAssumptionsLimited_monosat
from monosat_header cimport solveLimited as _solveLimited_monosat
from monosat_header cimport true_lit as _true_lit_monosat


def acyclic_directed( S ,  G ):
    """Cython signature: int acyclic_directed(void* S, void* G)"""

    


    cdef int _r = _acyclic_directed_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(G,NULL)))
    py_result = <int>_r
    return py_result

def acyclic_undirected( S ,  G ):
    """Cython signature: int acyclic_undirected(void* S, void* G)"""

    


    cdef int _r = _acyclic_undirected_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(G,NULL)))
    py_result = <int>_r
    return py_result

def addBinaryClause( S ,  lit1 ,  lit2 ):
    """Cython signature: bint addBinaryClause(void* S, int lit1, int lit2)"""

    assert isinstance(lit1, (int, long)), 'arg lit1 wrong type'
    assert isinstance(lit2, (int, long)), 'arg lit2 wrong type'



    cdef bint _r = _addBinaryClause_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<int>lit1), (<int>lit2))
    py_result = <bint>_r
    return py_result

def addBinaryClauses( S ,  first_args ,  second_args ,  n_pairs ):
    """Cython signature: void addBinaryClauses(void* S, int* first_args, int* second_args, int n_pairs)"""


    cdef array.array a =  array.array('i', first_args)
    cdef array.array b =  array.array('i', second_args)

    _addBinaryClauses_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<int*>a.data.as_ints),(<int*>b.data.as_ints), (<int>n_pairs))

def addClause( S ,  lits ,  n_lits ):
    """Cython signature: bint addClause(void* S, int* lits, int n_lits)"""

    assert isinstance(lits, list), 'arg lits wrong type'
    assert isinstance(n_lits, (int, long)), 'arg n_lits wrong type'


    cdef array.array a =  array.array('i', lits)
    cdef bint _r = _addClause_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<int*>a.data.as_ints), (<int>n_lits))
    py_result = <bint>_r
    return py_result

def addRoutingNet( S ,  G ,  router ,  disabledEdge ,  n_members ,  edge_lits ,  reach_lits ):
    """Cython signature: void addRoutingNet(void* S, void* G, void* router, int disabledEdge, int n_members, int* edge_lits, int* reach_lits)"""

    
    assert isinstance(router, (int, long)), 'arg router wrong type'
    assert isinstance(disabledEdge, (int, long)), 'arg disabledEdge wrong type'
    assert isinstance(n_members, (int, long)), 'arg n_members wrong type'
    assert isinstance(edge_lits, (int, long)), 'arg edge_lits wrong type'
    assert isinstance(reach_lits, (int, long)), 'arg reach_lits wrong type'


    cdef array.array a =  array.array('i', edge_lits)
    cdef array.array b =  array.array('i', reach_lits)
    _addRoutingNet_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(G,NULL)), (<void*>router), (<int>disabledEdge), (<int>n_members), (<int*>a.data.as_ints), (<int*>b.data.as_ints))

def addTertiaryClause( S ,  lit1 ,  lit2 ,  lit3 ):
    """Cython signature: bint addTertiaryClause(void* S, int lit1, int lit2, int lit3)"""

    assert isinstance(lit1, (int, long)), 'arg lit1 wrong type'
    assert isinstance(lit2, (int, long)), 'arg lit2 wrong type'
    assert isinstance(lit3, (int, long)), 'arg lit3 wrong type'




    cdef bint _r = _addTertiaryClause_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<int>lit1), (<int>lit2), (<int>lit3))
    py_result = <bint>_r
    return py_result

def addUnitClause( S ,  lit ):
    """Cython signature: bint addUnitClause(void* S, int lit)"""

    assert isinstance(lit, (int, long)), 'arg lit wrong type'


    cdef bint _r = _addUnitClause_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<int>lit))
    py_result = <bint>_r
    return py_result

def assertPB_eq( S ,  rhs ,  n_args ,  literals ,  coefficients ):
    """Cython signature: void assertPB_eq(void* S, int rhs, int n_args, int* literals, int* coefficients)"""

    assert isinstance(rhs, (int, long)), 'arg rhs wrong type'
    assert isinstance(n_args, (int, long)), 'arg n_args wrong type'
    assert isinstance(literals, (int, long)), 'arg literals wrong type'
    assert isinstance(coefficients, (int, long)), 'arg coefficients wrong type'

    cdef array.array a =  array.array('i', literals)
    cdef array.array b =  array.array('i', coefficients)



    _assertPB_eq_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<int>rhs), (<int>n_args), (<int*>a.data.as_ints), (<int*>b.data.as_ints))

def assertPB_geq( S ,  rhs ,  n_args ,  literals ,  coefficients ):
    """Cython signature: void assertPB_geq(void* S, int rhs, int n_args, int* literals, int* coefficients)"""

    assert isinstance(rhs, (int, long)), 'arg rhs wrong type'
    assert isinstance(n_args, (int, long)), 'arg n_args wrong type'
    assert isinstance(literals, (int, long)), 'arg literals wrong type'
    assert isinstance(coefficients, (int, long)), 'arg coefficients wrong type'

    cdef array.array a =  array.array('i', literals)
    cdef array.array b =  array.array('i', coefficients)



    _assertPB_geq_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<int>rhs), (<int>n_args), (<int*>a.data.as_ints), (<int*>b.data.as_ints))

def assertPB_gt( S ,  rhs ,  n_args ,  literals ,  coefficients ):
    """Cython signature: void assertPB_gt(void* S, int rhs, int n_args, int* literals, int* coefficients)"""

    assert isinstance(rhs, (int, long)), 'arg rhs wrong type'
    assert isinstance(n_args, (int, long)), 'arg n_args wrong type'
    assert isinstance(literals, (int, long)), 'arg literals wrong type'
    assert isinstance(coefficients, (int, long)), 'arg coefficients wrong type'


    cdef array.array a =  array.array('i', literals)
    cdef array.array b =  array.array('i', coefficients)



    _assertPB_gt_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<int>rhs), (<int>n_args), (<int*>a.data.as_ints), (<int*>b.data.as_ints))

def assertPB_leq( S ,  rhs ,  n_args ,  literals ,  coefficients ):
    """Cython signature: void assertPB_leq(void* S, int rhs, int n_args, int* literals, int* coefficients)"""

    assert isinstance(rhs, (int, long)), 'arg rhs wrong type'
    assert isinstance(n_args, (int, long)), 'arg n_args wrong type'
    assert isinstance(literals, (int, long)), 'arg literals wrong type'
    assert isinstance(coefficients, (int, long)), 'arg coefficients wrong type'


    cdef array.array a =  array.array('i', literals)
    cdef array.array b =  array.array('i', coefficients)



    _assertPB_leq_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<int>rhs), (<int>n_args),  (<int*>a.data.as_ints), (<int*>b.data.as_ints))

def assertPB_lt( S ,  rhs ,  n_args ,  literals ,  coefficients ):
    """Cython signature: void assertPB_lt(void* S, int rhs, int n_args, int* literals, int* coefficients)"""

    assert isinstance(rhs, (int, long)), 'arg rhs wrong type'
    assert isinstance(n_args, (int, long)), 'arg n_args wrong type'
    assert isinstance(literals, (int, long)), 'arg literals wrong type'
    assert isinstance(coefficients, (int, long)), 'arg coefficients wrong type'


    cdef array.array a =  array.array('i', literals)
    cdef array.array b =  array.array('i', coefficients)



    _assertPB_lt_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<int>rhs), (<int>n_args), (<int*>a.data.as_ints), (<int*>b.data.as_ints))

def at_most_one( S ,  vars ,  n_vars ):
    """Cython signature: void at_most_one(void* S, int* vars, int n_vars)"""

    assert isinstance(vars, list), 'arg vars wrong type'
    assert isinstance(n_vars, (int, long)), 'arg n_vars wrong type'

    cdef array.array a =  array.array('i', vars)

    _at_most_one_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)),  (<int*>a.data.as_ints), (<int>n_vars))

def backtrack( S ):
    """Cython signature: void backtrack(void* S)"""


    _backtrack_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)))

def bv_addition( S ,  bv ,  bvID1 ,  bvID2 ,  resultID ):
    """Cython signature: void bv_addition(void* S, void* bv, int bvID1, int bvID2, int resultID)"""

    
    assert isinstance(bvID1, (int, long)), 'arg bvID1 wrong type'
    assert isinstance(bvID2, (int, long)), 'arg bvID2 wrong type'
    assert isinstance(resultID, (int, long)), 'arg resultID wrong type'





    _bv_addition_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(bv,NULL)), (<int>bvID1), (<int>bvID2), (<int>resultID))

def bv_and( S ,  bv ,  bvaID ,  bvbID ,  bvResultID ):
    """Cython signature: void bv_and(void* S, void* bv, int bvaID, int bvbID, int bvResultID)"""

    
    assert isinstance(bvaID, (int, long)), 'arg bvaID wrong type'
    assert isinstance(bvbID, (int, long)), 'arg bvbID wrong type'
    assert isinstance(bvResultID, (int, long)), 'arg bvResultID wrong type'





    _bv_and_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(bv,NULL)), (<int>bvaID), (<int>bvbID), (<int>bvResultID))

def bv_bitblast( S ,  bv ,  bvID ):
    """Cython signature: void bv_bitblast(void* S, void* bv, int bvID)"""

    
    



    _bv_bitblast_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(bv,NULL)), (<int>bvID))

def bv_concat( S ,  bv ,  aID ,  bID ,  resultID ):
    """Cython signature: void bv_concat(void* S, void* bv, int aID, int bID, int resultID)"""

    
    assert isinstance(aID, (int, long)), 'arg aID wrong type'
    assert isinstance(bID, (int, long)), 'arg bID wrong type'
    assert isinstance(resultID, (int, long)), 'arg resultID wrong type'





    _bv_concat_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(bv,NULL)), (<int>aID), (<int>bID), (<int>resultID))

def bv_divide( S ,  bv ,  bvID1 ,  bvID2 ,  resultID ):
    """Cython signature: void bv_divide(void* S, void* bv, int bvID1, int bvID2, int resultID)"""

    
    assert isinstance(bvID1, (int, long)), 'arg bvID1 wrong type'
    assert isinstance(bvID2, (int, long)), 'arg bvID2 wrong type'
    assert isinstance(resultID, (int, long)), 'arg resultID wrong type'





    _bv_divide_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(bv,NULL)), (<int>bvID1), (<int>bvID2), (<int>resultID))

def bv_ite( S ,  bv ,  condition_lit ,  bvThenID ,  bvElseID ,  bvResultID ):
    """Cython signature: void bv_ite(void* S, void* bv, int condition_lit, int bvThenID, int bvElseID, int bvResultID)"""

    
    assert isinstance(condition_lit, (int, long)), 'arg condition_lit wrong type'
    assert isinstance(bvThenID, (int, long)), 'arg bvThenID wrong type'
    assert isinstance(bvElseID, (int, long)), 'arg bvElseID wrong type'
    assert isinstance(bvResultID, (int, long)), 'arg bvResultID wrong type'






    _bv_ite_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(bv,NULL)), (<int>condition_lit), (<int>bvThenID), (<int>bvElseID), (<int>bvResultID))

def bv_max( S ,  bv ,  args ,  n_args ,  resultID ):
    """Cython signature: void bv_max(void* S, void* bv, int* args, int n_args, int resultID)"""

    
    assert isinstance(args, (int, long)), 'arg args wrong type'
    assert isinstance(n_args, (int, long)), 'arg n_args wrong type'
    assert isinstance(resultID, (int, long)), 'arg resultID wrong type'


    cdef array.array a =  array.array('i', args)


    _bv_max_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(bv,NULL)),(<int*>a.data.as_ints), (<int>n_args), (<int>resultID))

def bv_min( S ,  bv ,  args ,  n_args ,  resultID ):
    """Cython signature: void bv_min(void* S, void* bv, int* args, int n_args, int resultID)"""

    
    assert isinstance(args, (int, long)), 'arg args wrong type'
    assert isinstance(n_args, (int, long)), 'arg n_args wrong type'
    assert isinstance(resultID, (int, long)), 'arg resultID wrong type'


    cdef array.array a =  array.array('i', args)


    _bv_min_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(bv,NULL)), (<int*>a.data.as_ints), (<int>n_args), (<int>resultID))

def bv_multiply( S ,  bv ,  bvID1 ,  bvID2 ,  resultID ):
    """Cython signature: void bv_multiply(void* S, void* bv, int bvID1, int bvID2, int resultID)"""

    
    assert isinstance(bvID1, (int, long)), 'arg bvID1 wrong type'
    assert isinstance(bvID2, (int, long)), 'arg bvID2 wrong type'
    assert isinstance(resultID, (int, long)), 'arg resultID wrong type'





    _bv_multiply_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(bv,NULL)), (<int>bvID1), (<int>bvID2), (<int>resultID))

def bv_nand( S ,  bv ,  bvaID ,  bvbID ,  bvResultID ):
    """Cython signature: void bv_nand(void* S, void* bv, int bvaID, int bvbID, int bvResultID)"""

    
    assert isinstance(bvaID, (int, long)), 'arg bvaID wrong type'
    assert isinstance(bvbID, (int, long)), 'arg bvbID wrong type'
    assert isinstance(bvResultID, (int, long)), 'arg bvResultID wrong type'





    _bv_nand_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(bv,NULL)), (<int>bvaID), (<int>bvbID), (<int>bvResultID))

def bv_nor( S ,  bv ,  bvaID ,  bvbID ,  bvResultID ):
    """Cython signature: void bv_nor(void* S, void* bv, int bvaID, int bvbID, int bvResultID)"""

    
    assert isinstance(bvaID, (int, long)), 'arg bvaID wrong type'
    assert isinstance(bvbID, (int, long)), 'arg bvbID wrong type'
    assert isinstance(bvResultID, (int, long)), 'arg bvResultID wrong type'





    _bv_nor_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(bv,NULL)), (<int>bvaID), (<int>bvbID), (<int>bvResultID))

def bv_not( S ,  bv ,  bvaID ,  bvResultID ):
    """Cython signature: void bv_not(void* S, void* bv, int bvaID, int bvResultID)"""

    
    assert isinstance(bvaID, (int, long)), 'arg bvaID wrong type'
    assert isinstance(bvResultID, (int, long)), 'arg bvResultID wrong type'




    _bv_not_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(bv,NULL)), (<int>bvaID), (<int>bvResultID))

def bv_or( S ,  bv ,  bvaID ,  bvbID ,  bvResultID ):
    """Cython signature: void bv_or(void* S, void* bv, int bvaID, int bvbID, int bvResultID)"""

    
    assert isinstance(bvaID, (int, long)), 'arg bvaID wrong type'
    assert isinstance(bvbID, (int, long)), 'arg bvbID wrong type'
    assert isinstance(bvResultID, (int, long)), 'arg bvResultID wrong type'





    _bv_or_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(bv,NULL)), (<int>bvaID), (<int>bvbID), (<int>bvResultID))

def bv_popcount( S ,  bv ,  args ,  n_args ,  resultID ):
    """Cython signature: void bv_popcount(void* S, void* bv, int* args, int n_args, int resultID)"""

    
    assert isinstance(args, (int, long)), 'arg args wrong type'
    assert isinstance(n_args, (int, long)), 'arg n_args wrong type'
    assert isinstance(resultID, (int, long)), 'arg resultID wrong type'


    cdef array.array a =  array.array('i', args)


    _bv_popcount_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(bv,NULL)),(<int*>a.data.as_ints), (<int>n_args), (<int>resultID))

def bv_slice( S ,  bv ,  aID ,  lower ,  upper ,  resultID ):
    """Cython signature: void bv_slice(void* S, void* bv, int aID, int lower, int upper, int resultID)"""

    
    assert isinstance(aID, (int, long)), 'arg aID wrong type'
    assert isinstance(lower, (int, long)), 'arg lower wrong type'
    assert isinstance(upper, (int, long)), 'arg upper wrong type'
    assert isinstance(resultID, (int, long)), 'arg resultID wrong type'






    _bv_slice_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(bv,NULL)), (<int>aID), (<int>lower), (<int>upper), (<int>resultID))

def bv_subtraction( S ,  bv ,  bvID1 ,  bvID2 ,  resultID ):
    """Cython signature: void bv_subtraction(void* S, void* bv, int bvID1, int bvID2, int resultID)"""

    
    assert isinstance(bvID1, (int, long)), 'arg bvID1 wrong type'
    assert isinstance(bvID2, (int, long)), 'arg bvID2 wrong type'
    assert isinstance(resultID, (int, long)), 'arg resultID wrong type'





    _bv_subtraction_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(bv,NULL)), (<int>bvID1), (<int>bvID2), (<int>resultID))

def bv_unary( S ,  bv ,  args ,  n_args ,  resultID ):
    """Cython signature: void bv_unary(void* S, void* bv, int* args, int n_args, int resultID)"""

    
    assert isinstance(args, (int, long)), 'arg args wrong type'
    assert isinstance(n_args, (int, long)), 'arg n_args wrong type'
    assert isinstance(resultID, (int, long)), 'arg resultID wrong type'


    cdef array.array a =  array.array('i', args)


    _bv_unary_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(bv,NULL)), (<int*>a.data.as_ints), (<int>n_args), (<int>resultID))

def bv_width( S ,  bv ,  bvID ):
    """Cython signature: int bv_width(void* S, void* bv, int bvID)"""

    
    



    cdef int _r = _bv_width_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(bv,NULL)), (<int>bvID))
    py_result = <int>_r
    return py_result

def bv_xnor( S ,  bv ,  bvaID ,  bvbID ,  bvResultID ):
    """Cython signature: void bv_xnor(void* S, void* bv, int bvaID, int bvbID, int bvResultID)"""

    
    assert isinstance(bvaID, (int, long)), 'arg bvaID wrong type'
    assert isinstance(bvbID, (int, long)), 'arg bvbID wrong type'
    assert isinstance(bvResultID, (int, long)), 'arg bvResultID wrong type'





    _bv_xnor_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(bv,NULL)), (<int>bvaID), (<int>bvbID), (<int>bvResultID))

def bv_xor( S ,  bv ,  bvaID ,  bvbID ,  bvResultID ):
    """Cython signature: void bv_xor(void* S, void* bv, int bvaID, int bvbID, int bvResultID)"""

    
    assert isinstance(bvaID, (int, long)), 'arg bvaID wrong type'
    assert isinstance(bvbID, (int, long)), 'arg bvbID wrong type'
    assert isinstance(bvResultID, (int, long)), 'arg bvResultID wrong type'





    _bv_xor_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(bv,NULL)), (<int>bvaID), (<int>bvbID), (<int>bvResultID))

def clearOptimizationObjectives( S ):
    """Cython signature: void clearOptimizationObjectives(void* S)"""


    _clearOptimizationObjectives_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)))

def createFlowRouting( S ,  G ,  sourceNode ,  destNode ,  maxflowLit ):
    """Cython signature: void* createFlowRouting(void* S, void* G, int sourceNode, int destNode, int maxflowLit)"""

    
    assert isinstance(sourceNode, (int, long)), 'arg sourceNode wrong type'
    assert isinstance(destNode, (int, long)), 'arg destNode wrong type'
    assert isinstance(maxflowLit, (int, long)), 'arg maxflowLit wrong type'





    cdef void* _r = _createFlowRouting_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(G,NULL)), (<int>sourceNode), (<int>destNode), (<int>maxflowLit))
    py_result =   pycapsule.PyCapsule_New(_r, NULL, NULL)
    return py_result

def deleteSolver( S ):
    """Cython signature: void deleteSolver(void* S)"""


    _deleteSolver_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)))

def disablePreprocessing( S ):
    """Cython signature: void disablePreprocessing(void* S)"""


    _disablePreprocessing_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)))

def disallowLiteralSimplification( S ,  lit ):
    """Cython signature: bint disallowLiteralSimplification(void* S, int lit)"""

    assert isinstance(lit, (int, long)), 'arg lit wrong type'


    cdef bint _r = _disallowLiteralSimplification_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<int>lit))
    py_result = <bint>_r
    return py_result

def flushPB( S ):
    """Cython signature: void flushPB(void* S)"""


    _flushPB_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)))


def getConflictClause( S ,  store_clause ,  max_store_size ):
    """Cython signature: int getConflictClause(void* S, int* store_clause, int max_store_size)"""

    assert isinstance(store_clause, list), 'arg store_clause wrong type'
    assert isinstance(max_store_size, (int, long)), 'arg max_store_size wrong type'

    cdef array.array a =  array.array('i', store_clause)

    cdef int _r = _getConflictClause_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)),(<int*>a.data.as_ints), (<int>max_store_size))
    py_result = <int>_r
    return py_result

def minimizeUnsatCore( S ,  assumptions ,  n_lits ):
    """Cython signature: int minimizeUnsatCore(void* S, int* assumptions, int n_lits)"""

    assert isinstance(assumptions, list), 'arg assumptions wrong type'
    assert isinstance(n_lits, (int, long)), 'arg n_lits wrong type'

    cdef array.array a =  array.array('i', assumptions)

    cdef int _r = _minimizeUnsatCore_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)),(<int*>a.data.as_ints), (<int>n_lits))
    py_result = <int>_r
    return py_result

def minimizeConflictClause( S ):
    """Cython signature: int minimizeConflictClause(void* S)"""
    _minimizeConflictClause_monosat(<void*>pycapsule.PyCapsule_GetPointer(S,NULL));

def getDecisionPolarity( S ,  v ):
    """Cython signature: bint getDecisionPolarity(void* S, int v)"""

    assert isinstance(v, (int, long)), 'arg v wrong type'


    cdef bint _r = _getDecisionPolarity_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<int>v))
    py_result = <bint>_r
    return py_result

def getDecisionPriority( S ,  var ):
    """Cython signature: int getDecisionPriority(void* S, int var)"""

    assert isinstance(var, (int, long)), 'arg var wrong type'


    cdef int _r = _getDecisionPriority_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<int>var))
    py_result = <int>_r
    return py_result

def getModel_AcyclicEdgeFlow( S ,  G ,  maxflow_literal ,  edgeLit ):
    """Cython signature: int64_t getModel_AcyclicEdgeFlow(void* S, void* G, int maxflow_literal, int edgeLit)"""

    
    assert isinstance(maxflow_literal, (int, long)), 'arg maxflow_literal wrong type'
    assert isinstance(edgeLit, (int, long)), 'arg edgeLit wrong type'




    cdef int64_t _r = _getModel_AcyclicEdgeFlow_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(G,NULL)), (<int>maxflow_literal), (<int>edgeLit))
    py_result = <int64_t>_r
    return py_result

def getModel_BV( S ,  bv ,  bvID ,  getMaximumValue ):
    """Cython signature: int64_t getModel_BV(void* S, void* bv, int bvID, bint getMaximumValue)"""

    
    
    assert isinstance(getMaximumValue, (int, long)), 'arg getMaximumValue wrong type'




    cdef int64_t _r = _getModel_BV_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(bv,NULL)), (<int>bvID), (<bint>getMaximumValue))
    py_result = <int64_t>_r
    return py_result

def getModel_EdgeFlow( S ,  G ,  maxflow_literal ,  edgeLit ):
    """Cython signature: int64_t getModel_EdgeFlow(void* S, void* G, int maxflow_literal, int edgeLit)"""

    
    assert isinstance(maxflow_literal, (int, long)), 'arg maxflow_literal wrong type'
    assert isinstance(edgeLit, (int, long)), 'arg edgeLit wrong type'




    cdef int64_t _r = _getModel_EdgeFlow_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(G,NULL)), (<int>maxflow_literal), (<int>edgeLit))
    py_result = <int64_t>_r
    return py_result

def getModel_Literal( S ,  lit ):
    """Cython signature: int getModel_Literal(void* S, int lit)"""

    assert isinstance(lit, (int, long)), 'arg lit wrong type'


    cdef int _r = _getModel_Literal_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<int>lit))
    py_result = <int>_r
    return py_result

def getModel_MaxFlow( S ,  G ,  maxflow_literal ):
    """Cython signature: int64_t getModel_MaxFlow(void* S, void* G, int maxflow_literal)"""

    
    assert isinstance(maxflow_literal, (int, long)), 'arg maxflow_literal wrong type'



    cdef int64_t _r = _getModel_MaxFlow_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(G,NULL)), (<int>maxflow_literal))
    py_result = <int64_t>_r
    return py_result

def getModel_MinimumSpanningTreeWeight( S ,  G ,  spanning_tree_literal ):
    """Cython signature: int64_t getModel_MinimumSpanningTreeWeight(void* S, void* G, int spanning_tree_literal)"""

    
    assert isinstance(spanning_tree_literal, (int, long)), 'arg spanning_tree_literal wrong type'



    cdef int64_t _r = _getModel_MinimumSpanningTreeWeight_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(G,NULL)), (<int>spanning_tree_literal))
    py_result = <int64_t>_r
    return py_result

def getModel_Path_EdgeLits( S ,  G ,  reach_or_distance_literal ,  store_length ,  store ):
    """Cython signature: int getModel_Path_EdgeLits(void* S, void* G, int reach_or_distance_literal, int store_length, int* store)"""

    
    assert isinstance(reach_or_distance_literal, (int, long)), 'arg reach_or_distance_literal wrong type'
    assert isinstance(store_length,(int,long)), 'arg store_length wrong type'
    assert isinstance(store, list), 'arg store wrong type'


    initializer = list(range(store_length))
    assert(len(initializer)==store_length)
    cdef array.array a =  array.array('i', initializer) #is there a better way to initialize an array.array to length n?
    cdef int _r = _getModel_Path_EdgeLits_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(G,NULL)), (<int>reach_or_distance_literal), (<int>store_length), (<int*>a.data.as_ints))
    store.clear()
    for lit in a:
        store.append(lit)
    assert(len(store)==store_length)
    py_result = <int>_r
    return py_result

def getModel_Path_EdgeLits_Length( S ,  G ,  reach_or_distance_literal ):
    """Cython signature: int getModel_Path_EdgeLits_Length(void* S, void* G, int reach_or_distance_literal)"""

    
    assert isinstance(reach_or_distance_literal, (int, long)), 'arg reach_or_distance_literal wrong type'
    cdef int _r = _getModel_Path_EdgeLits_Length_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(G,NULL)), (<int>reach_or_distance_literal))
    py_result = <int>_r
    return py_result

def getModel_Path_Nodes( S ,  G ,  reach_or_distance_literal ,  store_length ,  store ):
    """Cython signature: int getModel_Path_Nodes(void* S, void* G, int reach_or_distance_literal, int store_length, int* store)"""

    
    assert isinstance(reach_or_distance_literal, (int, long)), 'arg reach_or_distance_literal wrong type'
    assert isinstance(store_length, (int,long)), 'arg store_length wrong type'
    assert isinstance(store, list), 'arg store wrong type'


    initializer = list(range(store_length))
    assert(len(initializer)==store_length)
    cdef array.array a =  array.array('i', initializer) #is there a better way to initialize an array.array to length n?

    cdef int _r = _getModel_Path_Nodes_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(G,NULL)), (<int>reach_or_distance_literal), (<int>store_length),(<int*>a.data.as_ints))
    store.clear()
    for node in a:
        store.append(node)
    assert(len(store)==store_length)
    py_result = <int>_r
    return py_result

def getModel_Path_Nodes_Length( S ,  G ,  reach_or_distance_literal ):
    """Cython signature: int getModel_Path_Nodes_Length(void* S, void* G, int reach_or_distance_literal)"""

    
    assert isinstance(reach_or_distance_literal, (int, long)), 'arg reach_or_distance_literal wrong type'



    cdef int _r = _getModel_Path_Nodes_Length_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(G,NULL)), (<int>reach_or_distance_literal))
    py_result = <int>_r
    return py_result

def getVersion():
    """Cython signature: char * getVersion()"""
    cdef char  * _r = <char *>_getVersion_monosat()
    py_result = <char *>(_r)
    return py_result

def graph_setAssignEdgesToWeight( S ,  G ,  weight ):
    """Cython signature: void graph_setAssignEdgesToWeight(void* S, void* G, int64_t weight)"""

    
    assert isinstance(weight, (int, long)), 'arg weight wrong type'



    _graph_setAssignEdgesToWeight_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(G,NULL)), (<int64_t>weight))

def initBVTheory( S ):
    """Cython signature: void* initBVTheory(void* S)"""
    #

    #print "{0:x}".format(<uintptr_t>(<void*>pycapsule.PyCapsule_GetPointer(S,NULL)))


    cdef void* _r = _initBVTheory_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)))
    py_result =   pycapsule.PyCapsule_New(_r, NULL, NULL)
    return py_result


def isDecisionVar( S ,  var ):
    """Cython signature: bint isDecisionVar(void* S, int var)"""

    assert isinstance(var, (int, long)), 'arg var wrong type'


    cdef bint _r = _isDecisionVar_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<int>var))
    py_result = <bint>_r
    return py_result

def lastSolutionWasOptimal( S ):
    """Cython signature: bint lastSolutionWasOptimal(void* S)"""


    cdef bint _r = _lastSolutionWasOptimal_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)))
    py_result = <bint>_r
    return py_result

def maximizeBV( S ,  bv ,  bvID ):
    """Cython signature: void maximizeBV(void* S, void* bv, int bvID)"""

   
    _maximizeBV_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(bv,NULL)), (<int>bvID))

def maximizeLits( S ,  lits ,  n_lits ):
    """Cython signature: void maximizeLits(void* S, int* lits, int n_lits)"""

    assert isinstance(lits, list), 'arg lits wrong type'
    assert isinstance(n_lits, (int, long)), 'arg n_lits wrong type'

    cdef array.array a =  array.array('i', lits)

    _maximizeLits_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<int*>a.data.as_ints), (<int>n_lits))

def maximizeWeightedLits( S ,  lits ,  weights ,  n_lits ):
    """Cython signature: void maximizeWeightedLits(void* S, int* lits, int* weights, int n_lits)"""

    assert isinstance(lits, list), 'arg lits wrong type'
    assert isinstance(weights, (int, long)), 'arg weights wrong type'
    assert isinstance(n_lits, (int, long)), 'arg n_lits wrong type'

    cdef array.array a =  array.array('i', lits)
    cdef array.array b =  array.array('i', weights)

    _maximizeWeightedLits_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<int*>a.data.as_ints), (<int*>b.data.as_ints), (<int>n_lits))

def maximumFlow_geq( S ,  G ,  source ,  sink ,  weight ):
    """Cython signature: int maximumFlow_geq(void* S, void* G, int source, int sink, int64_t weight)"""

    
    assert isinstance(source, (int, long)), 'arg source wrong type'
    assert isinstance(sink, (int, long)), 'arg sink wrong type'
    assert isinstance(weight, (int, long)), 'arg weight wrong type'





    cdef int _r = _maximumFlow_geq_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(G,NULL)), (<int>source), (<int>sink), (<int64_t>weight))
    py_result = <int>_r
    return py_result

def maximumFlow_geq_bv( S ,  G ,  source ,  sink ,  bvID ):
    """Cython signature: int maximumFlow_geq_bv(void* S, void* G, int source, int sink, int bvID)"""

    
    assert isinstance(source, (int, long)), 'arg source wrong type'
    assert isinstance(sink, (int, long)), 'arg sink wrong type'
    





    cdef int _r = _maximumFlow_geq_bv_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(G,NULL)), (<int>source), (<int>sink), (<int>bvID))
    py_result = <int>_r
    return py_result

def maximumFlow_gt( S ,  G ,  source ,  sink ,  weight ):
    """Cython signature: int maximumFlow_gt(void* S, void* G, int source, int sink, int64_t weight)"""

    
    assert isinstance(source, (int, long)), 'arg source wrong type'
    assert isinstance(sink, (int, long)), 'arg sink wrong type'
    assert isinstance(weight, (int, long)), 'arg weight wrong type'





    cdef int _r = _maximumFlow_gt_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(G,NULL)), (<int>source), (<int>sink), (<int64_t>weight))
    py_result = <int>_r
    return py_result

def maximumFlow_gt_bv( S ,  G ,  source ,  sink ,  bvID ):
    """Cython signature: int maximumFlow_gt_bv(void* S, void* G, int source, int sink, int bvID)"""

    
    assert isinstance(source, (int, long)), 'arg source wrong type'
    assert isinstance(sink, (int, long)), 'arg sink wrong type'
    





    cdef int _r = _maximumFlow_gt_bv_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(G,NULL)), (<int>source), (<int>sink), (<int>bvID))
    py_result = <int>_r
    return py_result

def minimizeBV( S ,  bv ,  bvID ):
    """Cython signature: void minimizeBV(void* S, void* bv, int bvID)"""

    
    



    _minimizeBV_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(bv,NULL)), (<int>bvID))

def minimizeLits( S ,  lits ,  n_lits ):
    """Cython signature: void minimizeLits(void* S, int* lits, int n_lits)"""

    assert isinstance(lits, list), 'arg lits wrong type'
    assert isinstance(n_lits, (int, long)), 'arg n_lits wrong type'

    cdef array.array a =  array.array('i', lits)


    _minimizeLits_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<int*>a.data.as_ints), (<int>n_lits))

def minimizeWeightedLits( S ,  lits ,  weights ,  n_lits ):
    """Cython signature: void minimizeWeightedLits(void* S, int* lits, int* weights, int n_lits)"""

    assert isinstance(lits, list), 'arg lits wrong type'
    assert isinstance(weights, (int, long)), 'arg weights wrong type'
    assert isinstance(n_lits, (int, long)), 'arg n_lits wrong type'

    cdef array.array a =  array.array('i', lits)
    cdef array.array b =  array.array('i', weights)


    _minimizeWeightedLits_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<int*>a.data.as_ints), (<int*>b.data.as_ints), (<int>n_lits))

def minimumSpanningTree_leq( S ,  G ,  weight ):
    """Cython signature: int minimumSpanningTree_leq(void* S, void* G, int64_t weight)"""

    
    assert isinstance(weight, (int, long)), 'arg weight wrong type'



    cdef int _r = _minimumSpanningTree_leq_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(G,NULL)), (<int64_t>weight))
    py_result = <int>_r
    return py_result

def minimumSpanningTree_lt( S ,  G ,  weight ):
    """Cython signature: int minimumSpanningTree_lt(void* S, void* G, int64_t weight)"""
    assert isinstance(weight, (int, long)), 'arg weight wrong type'

    cdef int _r = _minimumSpanningTree_lt_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(G,NULL)), (<int64_t>weight))
    py_result = <int>_r
    return py_result

def nBitvectors( S ,  bv ):
    """Cython signature: int nBitvectors(void* S, void* bv)"""

    


    cdef int _r = _nBitvectors_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(bv,NULL)))
    py_result = <int>_r
    return py_result

def nClauses( S ):
    """Cython signature: int nClauses(void* S)"""


    cdef int _r = _nClauses_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)))
    py_result = <int>_r
    return py_result

def nEdges( S ,  G ):
    """Cython signature: int nEdges(void* S, void* G)"""

    


    cdef int _r = _nEdges_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(G,NULL)))
    py_result = <int>_r
    return py_result

def nNodes( S ,  G ):
    """Cython signature: int nNodes(void* S, void* G)"""

    


    cdef int _r = _nNodes_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(G,NULL)))
    py_result = <int>_r
    return py_result

def nVars( S ):
    """Cython signature: int nVars(void* S)"""


    cdef int _r = _nVars_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)))
    py_result = <int>_r
    return py_result

def newBVComparison_bv_geq( S ,  bv ,  bvID ,  compareID ):
    """Cython signature: int newBVComparison_bv_geq(void* S, void* bv, int bvID, int compareID)"""

    
    
    assert isinstance(compareID, (int, long)), 'arg compareID wrong type'




    cdef int _r = _newBVComparison_bv_geq_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(bv,NULL)), (<int>bvID), (<int>compareID))
    py_result = <int>_r
    return py_result

def newBVComparison_bv_gt( S ,  bv ,  bvID ,  compareID ):
    """Cython signature: int newBVComparison_bv_gt(void* S, void* bv, int bvID, int compareID)"""

    
    
    assert isinstance(compareID, (int, long)), 'arg compareID wrong type'




    cdef int _r = _newBVComparison_bv_gt_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(bv,NULL)), (<int>bvID), (<int>compareID))
    py_result = <int>_r
    return py_result

def newBVComparison_bv_leq( S ,  bv ,  bvID ,  compareID ):
    """Cython signature: int newBVComparison_bv_leq(void* S, void* bv, int bvID, int compareID)"""

    
    
    assert isinstance(compareID, (int, long)), 'arg compareID wrong type'




    cdef int _r = _newBVComparison_bv_leq_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(bv,NULL)), (<int>bvID), (<int>compareID))
    py_result = <int>_r
    return py_result

def newBVComparison_bv_lt( S ,  bv ,  bvID ,  compareID ):
    """Cython signature: int newBVComparison_bv_lt(void* S, void* bv, int bvID, int compareID)"""

    
    
    assert isinstance(compareID, (int, long)), 'arg compareID wrong type'




    cdef int _r = _newBVComparison_bv_lt_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(bv,NULL)), (<int>bvID), (<int>compareID))
    py_result = <int>_r
    return py_result

def newBVComparison_const_geq( S ,  bv ,  bvID ,  weight ):
    """Cython signature: int newBVComparison_const_geq(void* S, void* bv, int bvID, int64_t weight)"""

    
    
    assert isinstance(weight, (int, long)), 'arg weight wrong type'




    cdef int _r = _newBVComparison_const_geq_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(bv,NULL)), (<int>bvID), (<int64_t>weight))
    py_result = <int>_r
    return py_result

def newBVComparison_const_gt( S ,  bv ,  bvID ,  weight ):
    """Cython signature: int newBVComparison_const_gt(void* S, void* bv, int bvID, int64_t weight)"""

    
    
    assert isinstance(weight, (int, long)), 'arg weight wrong type'




    cdef int _r = _newBVComparison_const_gt_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(bv,NULL)), (<int>bvID), (<int64_t>weight))
    py_result = <int>_r
    return py_result

def newBVComparison_const_leq( S ,  bv ,  bvID ,  weight ):
    """Cython signature: int newBVComparison_const_leq(void* S, void* bv, int bvID, int64_t weight)"""

    
    
    assert isinstance(weight, (int, long)), 'arg weight wrong type'




    cdef int _r = _newBVComparison_const_leq_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(bv,NULL)), (<int>bvID), (<int64_t>weight))
    py_result = <int>_r
    return py_result

def newBVComparison_const_lt( S ,  bv ,  bvID ,  weight ):
    """Cython signature: int newBVComparison_const_lt(void* S, void* bv, int bvID, int64_t weight)"""

    
    
    assert isinstance(weight, (int, long)), 'arg weight wrong type'




    cdef int _r = _newBVComparison_const_lt_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(bv,NULL)), (<int>bvID), (<int64_t>weight))
    py_result = <int>_r
    return py_result

def newBitvector( S ,  bv ,  bits ,  n_bits ):
    """Cython signature: int newBitvector(void* S, void* bv, int* bits, int n_bits)"""

    
    assert isinstance(bits, list), 'arg bits wrong type'
    assert isinstance(n_bits, (int, long)), 'arg n_bits wrong type'


    cdef array.array a =  array.array('i', bits)

    cdef int _r = _newBitvector_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(bv,NULL)), (<int*>a.data.as_ints), (<int>n_bits))
    py_result = <int>_r
    return py_result

def newBitvector_anon( S ,  bv ,  bvWidth ):
    """Cython signature: int newBitvector_anon(void* S, void* bv, int bvWidth)"""

    
    assert isinstance(bvWidth, (int, long)), 'arg bvWidth wrong type'



    cdef int _r = _newBitvector_anon_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(bv,NULL)), (<int>bvWidth))
    py_result = <int>_r
    return py_result

def newBitvector_const( S ,  bv ,  bvWidth ,  constval ):
    """Cython signature: int newBitvector_const(void* S, void* bv, int bvWidth, int64_t constval)"""

    
    assert isinstance(bvWidth, (int, long)), 'arg bvWidth wrong type'
    assert isinstance(constval, (int, long)), 'arg constval wrong type'




    cdef int _r = _newBitvector_const_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(bv,NULL)), (<int>bvWidth), (<int64_t>constval))
    py_result = <int>_r
    return py_result

def newEdge( S ,  G ,  _from ,  to ,  weight ):
    """Cython signature: int newEdge(void* S, void* G, int _from, int to, int64_t weight)"""

    
    assert isinstance(_from, (int, long)), 'arg _from wrong type'
    assert isinstance(to, (int, long)), 'arg to wrong type'
    assert isinstance(weight, (int, long)), 'arg weight wrong type'





    cdef int _r = _newEdge_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(G,NULL)), (<int>_from), (<int>to), (<int64_t>weight))
    py_result = <int>_r
    return py_result

def newEdgeSet( S ,  G ,  edges ,  n_edges ,  enforceEdgeAssignment ):
    """Cython signature: void newEdgeSet(void* S, void* G, int* edges, int n_edges, bint enforceEdgeAssignment)"""
    cdef array.array a =  array.array('i', edges)
    _newEdgeSet_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(G,NULL)), (<int*>a.data.as_ints), (<int>n_edges), (<bint>enforceEdgeAssignment))

def newEdge_bv( S ,  G ,  _from ,  to ,  bvID ):
    """Cython signature: int newEdge_bv(void* S, void* G, int _from, int to, int bvID)"""

    
    assert isinstance(_from, (int, long)), 'arg _from wrong type'
    assert isinstance(to, (int, long)), 'arg to wrong type'
    





    cdef int _r = _newEdge_bv_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(G,NULL)), (<int>_from), (<int>to), (<int>bvID))
    py_result = <int>_r
    return py_result

def newEdge_double( S ,  G ,  _from ,  to , double weight ):
    """Cython signature: int newEdge_double(void* S, void* G, int _from, int to, double weight)"""

    
    assert isinstance(_from, (int, long)), 'arg _from wrong type'
    assert isinstance(to, (int, long)), 'arg to wrong type'
    assert isinstance(weight, float), 'arg weight wrong type'





    cdef int _r = _newEdge_double_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(G,NULL)), (<int>_from), (<int>to), (<double>weight))
    py_result = <int>_r
    return py_result


def newGraph( S ):
    """Cython signature: void* newGraph(void* S)"""


    cdef void* _r = _newGraph_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)))
    py_result =   pycapsule.PyCapsule_New(_r, NULL, NULL)
    return py_result

def newNode( S ,  G ):
    """Cython signature: int newNode(void* S, void* G)"""

    


    cdef int _r = _newNode_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(G,NULL)))
    py_result = <int>_r
    return py_result

def newSolver():
    """Cython signature: void* newSolver()"""
    cdef void* _r = _newSolver_monosat()
    if _r is NULL:
        raise MemoryError()
    py_result =   pycapsule.PyCapsule_New(_r, NULL, NULL)
    return py_result

def newSolver_arg(bytes argv ):
    """Cython signature: void* newSolver_arg(char * argv)"""
    assert isinstance(argv, bytes), 'arg argv wrong type'

    cdef void* _r = _newSolver_arg_monosat((<char *>argv))
    #print "{0:x}".format(<uintptr_t>(_r))
    if _r is NULL:
        raise MemoryError()
    py_result =   pycapsule.PyCapsule_New(_r, NULL, NULL)

    return py_result


def newVar( S ):
    """Cython signature: int newVar(void* S)"""
    

    cdef int _r = _newVar_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)))
    py_result = <int>_r
    return py_result

def reaches( S ,  G ,  _from ,  to ):
    """Cython signature: int reaches(void* S, void* G, int _from, int to)"""
    
    
    assert isinstance(_from, (int, long)), 'arg _from wrong type'
    assert isinstance(to, (int, long)), 'arg to wrong type'
    cdef int _r = _reaches_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(G,NULL)), (<int>_from), (<int>to))
    py_result = <int>_r
    return py_result

def reachesBackward( S ,  G ,  _from ,  to ):
    """Cython signature: int reaches(void* S, void* G, int _from, int to)"""


    assert isinstance(_from, (int, long)), 'arg _from wrong type'
    assert isinstance(to, (int, long)), 'arg to wrong type'
    cdef int _r = _reachesBackward_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(G,NULL)), (<int>_from), (<int>to))
    py_result = <int>_r
    return py_result


def onPath( S ,  G , nodeOnPath, _from ,  to ):
    """Cython signature: int reaches(void* S, void* G, int nodeOnPath,int _from, int to)"""


    assert isinstance(_from, (int, long)), 'arg _from wrong type'
    assert isinstance(to, (int, long)), 'arg to wrong type'
    assert isinstance(nodeOnPath, (int, long)), 'arg to wrong type'
    cdef int _r = _onPath_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(G,NULL)), (<int>nodeOnPath), (<int>_from), (<int>to))
    py_result = <int>_r
    return py_result

def readGNF( S , bytes filename ):
    """Cython signature: void readGNF(void* S, char * filename)"""
    assert isinstance(filename, bytes), 'arg filename wrong type'
    _readGNF_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<char *>filename))

def loadGNF( S , bytes filename ):
    """Cython signature: void loadGNF(void* S, char * filename)"""
    assert isinstance(filename, bytes), 'arg filename wrong type'
    _loadGNF_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<char *>filename))

def setConflictLimit( S ,  num_conflicts ):
    """Cython signature: void setConflictLimit(void* S, int num_conflicts)"""
    
    assert isinstance(num_conflicts, (int, long)), 'arg num_conflicts wrong type'


    _setConflictLimit_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<int>num_conflicts))

def setDecisionPolarity( S ,  v ,  b ):
    """Cython signature: void setDecisionPolarity(void* S, int v, bint b)"""
    
    assert isinstance(v, (int, long)), 'arg v wrong type'
    assert isinstance(b, (int, long)), 'arg b wrong type'



    _setDecisionPolarity_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<int>v), (<bint>b))

def setDecisionPriority( S ,  var ,  priority ):
    """Cython signature: void setDecisionPriority(void* S, int var, int priority)"""
    
    assert isinstance(var, (int, long)), 'arg var wrong type'
    assert isinstance(priority, (int, long)), 'arg priority wrong type'



    _setDecisionPriority_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<int>var), (<int>priority))

def setDecisionVar( S ,  var ,  decidable ):
    """Cython signature: void setDecisionVar(void* S, int var, bint decidable)"""
    
    assert isinstance(var, (int, long)), 'arg var wrong type'
    assert isinstance(decidable, (int, long)), 'arg decidable wrong type'



    _setDecisionVar_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<int>var), (<bint>decidable))

def setOutputFile( S , bytes output ):
    """Cython signature: void setOutputFile(void* S, char * output)"""
    
    assert isinstance(output, bytes), 'arg output wrong type'


    _setOutputFile_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<char *>output))

def setPropagationLimit( S ,  num_propagations ):
    """Cython signature: void setPropagationLimit(void* S, int num_propagations)"""
    
    assert isinstance(num_propagations, (int, long)), 'arg num_propagations wrong type'


    _setPropagationLimit_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<int>num_propagations))

def setTimeLimit( S ,  seconds ):
    """Cython signature: void setTimeLimit(void* S, int seconds)"""
    
    assert isinstance(seconds, (int, long)), 'arg seconds wrong type'


    _setTimeLimit_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<int>seconds))

def shortestPathUnweighted_leq_const( S ,  G ,  _from ,  to ,  steps ):
    """Cython signature: int shortestPathUnweighted_leq_const(void* S, void* G, int _from, int to, int steps)"""
    
    
    assert isinstance(_from, (int, long)), 'arg _from wrong type'
    assert isinstance(to, (int, long)), 'arg to wrong type'
    assert isinstance(steps, (int, long)), 'arg steps wrong type'





    cdef int _r = _shortestPathUnweighted_leq_const_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(G,NULL)), (<int>_from), (<int>to), (<int>steps))
    py_result = <int>_r
    return py_result

def shortestPathUnweighted_lt_const( S ,  G ,  _from ,  to ,  steps ):
    """Cython signature: int shortestPathUnweighted_lt_const(void* S, void* G, int _from, int to, int steps)"""
    
    
    assert isinstance(_from, (int, long)), 'arg _from wrong type'
    assert isinstance(to, (int, long)), 'arg to wrong type'
    assert isinstance(steps, (int, long)), 'arg steps wrong type'





    cdef int _r = _shortestPathUnweighted_lt_const_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(G,NULL)), (<int>_from), (<int>to), (<int>steps))
    py_result = <int>_r
    return py_result

def shortestPath_leq_bv( S ,  G ,  _from ,  to ,  bvID ):
    """Cython signature: int shortestPath_leq_bv(void* S, void* G, int _from, int to, int bvID)"""
    
    
    assert isinstance(_from, (int, long)), 'arg _from wrong type'
    assert isinstance(to, (int, long)), 'arg to wrong type'
    





    cdef int _r = _shortestPath_leq_bv_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(G,NULL)), (<int>_from), (<int>to), (<int>bvID))
    py_result = <int>_r
    return py_result

def shortestPath_leq_const( S ,  G ,  _from ,  to ,  dist ):
    """Cython signature: int shortestPath_leq_const(void* S, void* G, int _from, int to, int64_t dist)"""
    
    
    assert isinstance(_from, (int, long)), 'arg _from wrong type'
    assert isinstance(to, (int, long)), 'arg to wrong type'
    assert isinstance(dist, (int, long)), 'arg dist wrong type'





    cdef int _r = _shortestPath_leq_const_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(G,NULL)), (<int>_from), (<int>to), (<int64_t>dist))
    py_result = <int>_r
    return py_result

def shortestPath_lt_bv( S ,  G ,  _from ,  to ,  bvID ):
    """Cython signature: int shortestPath_lt_bv(void* S, void* G, int _from, int to, int bvID)"""
    
    
    assert isinstance(_from, (int, long)), 'arg _from wrong type'
    assert isinstance(to, (int, long)), 'arg to wrong type'
    





    cdef int _r = _shortestPath_lt_bv_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(G,NULL)), (<int>_from), (<int>to), (<int>bvID))
    py_result = <int>_r
    return py_result

def shortestPath_lt_const( S ,  G ,  _from ,  to ,  dist ):
    """Cython signature: int shortestPath_lt_const(void* S, void* G, int _from, int to, int64_t dist)"""
    
    
    assert isinstance(_from, (int, long)), 'arg _from wrong type'
    assert isinstance(to, (int, long)), 'arg to wrong type'
    assert isinstance(dist, (int, long)), 'arg dist wrong type'





    cdef int _r = _shortestPath_lt_const_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(G,NULL)), (<int>_from), (<int>to), (<int64_t>dist))
    py_result = <int>_r
    return py_result

def connectedComponents_geq_const( S ,  G , components):
    """Cython signature: int connectedComponents_geq_const(void* S, void* G, int components)"""

    assert isinstance(components, (int, long)), 'arg _from wrong type'

    cdef int _r = _connectedComponents_geq_const_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<void*>pycapsule.PyCapsule_GetPointer(G,NULL)), (<int>components))
    py_result = <int>_r
    return py_result

def solve( S ):
    """Cython signature: bint solve(void* S)"""
    

    cdef bint _r = _solve_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)))
    py_result = <bint>_r
    return py_result

def solveAssumptions( S ,  assumptions,  n_assumptions):
    """Cython signature: bint solveAssumptions(void* S, int* assumptions, int n_assumptions)"""
    

    cdef array.array a =  array.array('i', assumptions)
    cdef bint _r = _solveAssumptions_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), <int*>a.data.as_ints, (<int>n_assumptions))
    py_result = <bint>_r
    return py_result

def solveAssumptionsLimited( S ,  assumptions ,  n_assumptions ):
    """Cython signature: int solveAssumptionsLimited(void* S, int* assumptions, int n_assumptions)"""
    
    assert isinstance(assumptions, list), 'arg assumptions wrong type'
    assert isinstance(n_assumptions, (int, long)), 'arg n_assumptions wrong type'


    cdef array.array a =  array.array('i', assumptions)
    cdef int _r = _solveAssumptionsLimited_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)), (<int*>a.data.as_ints), (<int>n_assumptions))
    py_result = <int>_r
    return py_result

def solveLimited( S ):
    """Cython signature: int solveLimited(void* S)"""
    

    cdef int _r = _solveLimited_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)))
    py_result = <int>_r
    return py_result

def true_lit( S ):
    """Cython signature: int true_lit(void* S)"""


    cdef int _r = _true_lit_monosat((<void*>pycapsule.PyCapsule_GetPointer(S,NULL)))
    py_result = <int>_r
    return py_result 
