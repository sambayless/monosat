# The MIT License (MIT)
#
# Copyright (c) 2014, Sam Bayless
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
# associated documentation files (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge, publish, distribute,
# sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or
# substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
# NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
# OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
from __future__ import division
from __future__ import print_function

# Python interface to MonoSAT
# Includes _optional_ support for cython, otherwise falling back on ctypes

import os
import platform
from monosat.singleton import Singleton
from os import path
from enum import Enum
import warnings

module_path = os.path.abspath(path.dirname(__file__))
# If use_cython is true, then monosat will attempt to load the faster cython libraries first, falling back on ctypes if
# the cython interface is not available
# Set use_cython to false to force use of the ctypes api, rather than cython
use_cython = False
try:
    if not use_cython:
        raise Exception("Load ctypes")
    # attempt to load the (faster) Cython interface to monosat.
    # the must have been compiled previously, using setup.py
    # import pyximport
    # pyximport.install()
    import monosat.monosat_p

    print("Using cython interface")


    # cython doesn't use these conversion functions, but ctypes does, so define them as passthroughs if cython is used
    def c_int(x):
        return x


    def c_int64(x):
        return x


    def c_int64(x):
        return x


    def c_bvID(x):
        return x


    def c_bvID_p(x):
        return x


    def c_solver_p(x):
        return x


    def c_uintp(x):
        return x


    def c_intp(x):
        return x


    def c_var(x):
        return x


    def c_var_p(x):
        return x


    def c_bool(x):
        return x


    def c_bv_p(x):
        return x


    def c_fsm_theory_p(x):
        return x


    def c_fsm_p(x):
        return x


    def c_graph_p(x):
        return x


    def c_literal(x):
        return x


    def c_literal_p(x):
        return x


    def c_var(x):
        return x


    def c_var_p(x):
        return x


except:
    use_cython = False
    # Fall back on the cyypes interface
    from ctypes import *

    library_monosat = "libmonosat.so"
    if platform.system() == "Windows":
        library_monosat = "libmonosat.dll"
    elif platform.system() == "Darwin":
        library_monosat = "libmonosat.dylib"
    try:
        try:
            # First try loading monosat from the python library directory
            _monosat_c = cdll.LoadLibrary(os.path.join(module_path, library_monosat))
        except Exception as e:
            # then fall back to loading from the system path
            _monosat_c = cdll.LoadLibrary(library_monosat)
    except Exception as e:
        print("Unable to load libmonosat dynamic library")
        raise e

    # libc = CDLL('libc.so.6')
    null_ptr = POINTER(c_int)()

    c_uint_p = POINTER(c_int)
    c_int_p = POINTER(c_int)

    # Note: longs specified as 64bit on all platforms.
    c_long_p = POINTER(c_int64)
    null_ptr64 = POINTER(c_int64)()

    c_solver_p = POINTER(c_int64)
    c_graph_p = POINTER(c_int64)
    c_bv_p = POINTER(c_int64)
    c_fsm_theory_p = POINTER(c_int64)
    c_fsm_p = POINTER(c_int64)

    c_literal = c_int
    c_literal_p = c_int_p
    c_bvID = c_int
    c_bvID_p = c_int_p
    c_var = c_int
    c_var_p = c_int_p


def dimacs(l):
    assert isinstance(l, int)
    if l & 1:
        return -(l // 2 + 1)
    else:
        return l // 2 + 1


class Ineq(Enum):
    LT = -2
    LEQ = -1
    EQ = 0
    GEQ = 1
    GT = 2


class Solver:

    def __init__(self, monosat_c, arguments=None, output_file=None):
        self.elapsed_time = 0
        if arguments is None:
            arguments = []
        elif isinstance(arguments, str):
            # split this up by whitespace
            arguments = arguments.split()
        self.arguments = arguments
        arr = "monosat " + " ".join(arguments)

        self._ptr = monosat_c.newSolver_arg(bytes(arr, "utf-8"))
        if output_file is not None and len(output_file) > 0:
            monosat_c.setOutputFile(self._ptr, c_char_p(output_file.encode("ascii")))
        self.bvtheory = monosat_c.initBVTheory(self._ptr)
        self.arguments = None
        self.symbolmap = dict()
        self.graphs = []

        self._true = None

    def delete(self):
        Monosat().monosat_c.deleteSolver(self._ptr)
        self._ptr = None
        if Monosat().solver is self:
            Monosat().solver = None


# very simple, low-level python to Monosat's C interface.
class Monosat(metaclass=Singleton):

    def __init__(self):
        self._managers = dict()
        self.elapsed_time = 0
        if use_cython:
            self.monosat_c = monosat.monosat_p
        else:
            self.monosat_c = _monosat_c
            self._int_array = (c_int * (1024))()
            self._int_array2 = (c_int * (1024))()
            self._long_array = (c_int64 * (1024))()
            # Set the return types for each function

            self.monosat_c.getVersion.argtypes = []
            self.monosat_c.getVersion.restype = c_char_p

            self.monosat_c.newSolver.argtypes = []
            self.monosat_c.newSolver.restype = c_solver_p

            self.monosat_c.newSolver_arg.argtypes = [c_char_p]
            self.monosat_c.newSolver_arg.restype = c_solver_p

            self.monosat_c.newSolver_args.argtypes = [c_int, POINTER(c_char_p)]
            self.monosat_c.newSolver_args.restype = c_solver_p

            self.monosat_c.deleteSolver.argtypes = [c_solver_p]

            self.monosat_c.readGNF.argtypes = [c_solver_p, c_char_p]
            self.monosat_c.loadGNF.argtypes = [c_solver_p, c_char_p]

            self.monosat_c.solve.argtypes = [c_solver_p]
            self.monosat_c.solve.restype = c_bool

            self.monosat_c.solveAssumptions.argtypes = [c_solver_p, c_literal_p, c_int]
            self.monosat_c.solveAssumptions.restype = c_bool

            self.monosat_c.solveLimited.argtypes = [c_solver_p]
            self.monosat_c.solveLimited.restype = c_int

            self.monosat_c.solveAssumptionsLimited.argtypes = [
                c_solver_p,
                c_literal_p,
                c_int,
            ]
            self.monosat_c.solveAssumptionsLimited.restype = c_int

            self.monosat_c.lastSolutionWasOptimal.argtypes = [c_solver_p]
            self.monosat_c.lastSolutionWasOptimal.restype = c_bool

            self.monosat_c.getConflictClause.argtypes = [c_solver_p, c_int_p, c_int]
            self.monosat_c.getConflictClause.restype = c_int

            self.monosat_c.minimizeUnsatCore.argtypes = [c_solver_p, c_int_p, c_int]
            self.monosat_c.minimizeUnsatCore.restype = c_int

            self.monosat_c.minimizeConflictClause.argtypes = [c_solver_p]

            self.monosat_c.setTimeLimit.argtypes = [c_solver_p, c_int]

            self.monosat_c.setConflictLimit.argtypes = [c_solver_p, c_int]
            self.monosat_c.setPropagationLimit.argtypes = [c_solver_p, c_int]

            self.monosat_c.backtrack.argtypes = [c_solver_p]

            self.monosat_c.newVar.argtypes = [c_solver_p]
            self.monosat_c.newVar.restype = c_int

            self.monosat_c.setDecisionVar.argtypes = [c_solver_p, c_var, c_bool]

            self.monosat_c.isDecisionVar.argtypes = [c_solver_p, c_var]
            self.monosat_c.isDecisionVar.restype = c_bool

            self.monosat_c.setDecisionPriority.argtypes = [c_solver_p, c_var, c_int]

            self.monosat_c.getDecisionPriority.argtypes = [c_solver_p, c_var]
            self.monosat_c.getDecisionPriority.restype = c_int

            self.monosat_c.setDecisionPolarity.argtypes = [c_solver_p, c_var, c_bool]

            self.monosat_c.getDecisionPolarity.argtypes = [c_solver_p, c_var]
            self.monosat_c.getDecisionPolarity.restype = c_bool

            self.monosat_c.disallowLiteralSimplification.argtypes = [
                c_solver_p,
                c_literal,
            ]
            self.monosat_c.disallowLiteralSimplification.restype = c_bool

            self.monosat_c.addClause.argtypes = [c_solver_p, c_literal_p, c_int]
            self.monosat_c.addClause.restype = c_bool

            self.monosat_c.addUnitClause.argtypes = [c_solver_p, c_literal]
            self.monosat_c.addUnitClause.restype = c_bool

            self.monosat_c.addBinaryClause.argtypes = [c_solver_p, c_literal, c_literal]
            self.monosat_c.addBinaryClause.restype = c_bool

            self.monosat_c.addBinaryClauses.argtypes = [
                c_solver_p,
                c_int_p,
                c_int_p,
                c_int,
            ]

            self.monosat_c.addTertiaryClause.argtypes = [
                c_solver_p,
                c_literal,
                c_literal,
                c_literal,
            ]
            self.monosat_c.addTertiaryClause.restype = c_bool

            self.monosat_c.true_lit.argtypes = [c_solver_p]
            self.monosat_c.true_lit.restype = c_int

            self.monosat_c.clearOptimizationObjectives.argtypes = [c_solver_p]
            self.monosat_c.maximizeBV.argtypes = [c_solver_p, c_bv_p, c_int]
            self.monosat_c.minimizeBV.argtypes = [c_solver_p, c_bv_p, c_int]
            self.monosat_c.maximizeLits.argtypes = [c_solver_p, c_int_p, c_int]
            self.monosat_c.minimizeLits.argtypes = [c_solver_p, c_int_p, c_int]
            self.monosat_c.maximizeWeightedLits.argtypes = [
                c_solver_p,
                c_int_p,
                c_int_p,
                c_int,
            ]
            self.monosat_c.minimizeWeightedLits.argtypes = [
                c_solver_p,
                c_int_p,
                c_int_p,
                c_int,
            ]

            self.monosat_c.at_most_one.argtypes = [c_solver_p, c_var_p, c_int]
            self.monosat_c.assertPB_lt.argtypes = [
                c_solver_p,
                c_int,
                c_int,
                c_int_p,
                c_int_p,
            ]
            self.monosat_c.assertPB_leq.argtypes = [
                c_solver_p,
                c_int,
                c_int,
                c_int_p,
                c_int_p,
            ]
            self.monosat_c.assertPB_eq.argtypes = [
                c_solver_p,
                c_int,
                c_int,
                c_int_p,
                c_int_p,
            ]
            self.monosat_c.assertPB_geq.argtypes = [
                c_solver_p,
                c_int,
                c_int,
                c_int_p,
                c_int_p,
            ]
            self.monosat_c.assertPB_gt.argtypes = [
                c_solver_p,
                c_int,
                c_int,
                c_int_p,
                c_int_p,
            ]
            self.monosat_c.flushPB.argtypes = [c_solver_p]

            self.monosat_c.initBVTheory.argtypes = [c_solver_p]
            self.monosat_c.initBVTheory.restype = c_bv_p

            self.monosat_c.newBitvector_anon.argtypes = [c_solver_p, c_bv_p, c_int]
            self.monosat_c.newBitvector_anon.restype = c_bvID

            self.monosat_c.newBitvector_const.argtypes = [
                c_solver_p,
                c_bv_p,
                c_int,
                c_int64,
            ]
            self.monosat_c.newBitvector_const.restype = c_bvID

            self.monosat_c.newBitvector.argtypes = [c_solver_p, c_bv_p, c_var_p, c_int]
            self.monosat_c.newBitvector.restype = c_bvID

            self.monosat_c.nBitvectors.argtypes = [c_solver_p, c_bv_p]
            self.monosat_c.nBitvectors.restype = c_int
            self.monosat_c.newBVComparison_const_lt.argtypes = [
                c_solver_p,
                c_bv_p,
                c_bvID,
                c_int64,
            ]
            self.monosat_c.newBVComparison_const_lt.restype = c_literal

            self.monosat_c.newBVComparison_bv_lt.argtypes = [
                c_solver_p,
                c_bv_p,
                c_bvID,
                c_bvID,
            ]
            self.monosat_c.newBVComparison_bv_lt.restype = c_literal

            self.monosat_c.newBVComparison_const_leq.argtypes = [
                c_solver_p,
                c_bv_p,
                c_bvID,
                c_int64,
            ]
            self.monosat_c.newBVComparison_const_leq.restype = c_literal

            self.monosat_c.newBVComparison_bv_leq.argtypes = [
                c_solver_p,
                c_bv_p,
                c_bvID,
                c_bvID,
            ]
            self.monosat_c.newBVComparison_bv_leq.restype = c_literal

            self.monosat_c.newBVComparison_const_gt.argtypes = [
                c_solver_p,
                c_bv_p,
                c_bvID,
                c_int64,
            ]
            self.monosat_c.newBVComparison_const_gt.restype = c_literal

            self.monosat_c.newBVComparison_bv_gt.argtypes = [
                c_solver_p,
                c_bv_p,
                c_bvID,
                c_bvID,
            ]
            self.monosat_c.newBVComparison_bv_gt.restype = c_literal

            self.monosat_c.newBVComparison_const_geq.argtypes = [
                c_solver_p,
                c_bv_p,
                c_bvID,
                c_int64,
            ]
            self.monosat_c.newBVComparison_const_geq.restype = c_literal

            self.monosat_c.newBVComparison_bv_geq.argtypes = [
                c_solver_p,
                c_bv_p,
                c_bvID,
                c_bvID,
            ]
            self.monosat_c.newBVComparison_bv_geq.restype = c_literal

            self.monosat_c.bv_addition.argtypes = [
                c_solver_p,
                c_bv_p,
                c_bvID,
                c_bvID,
                c_bvID,
            ]
            self.monosat_c.bv_subtraction.argtypes = [
                c_solver_p,
                c_bv_p,
                c_bvID,
                c_bvID,
                c_bvID,
            ]

            self.monosat_c.bv_multiply.argtypes = [
                c_solver_p,
                c_bv_p,
                c_bvID,
                c_bvID,
                c_bvID,
            ]
            self.monosat_c.bv_divide.argtypes = [
                c_solver_p,
                c_bv_p,
                c_bvID,
                c_bvID,
                c_bvID,
            ]

            self.monosat_c.bv_ite.argtypes = [
                c_solver_p,
                c_bv_p,
                c_literal,
                c_bvID,
                c_bvID,
                c_bvID,
            ]

            self.monosat_c.bv_min.argtypes = [
                c_solver_p,
                c_bv_p,
                c_bvID_p,
                c_int,
                c_bvID,
            ]
            self.monosat_c.bv_max.argtypes = [
                c_solver_p,
                c_bv_p,
                c_bvID_p,
                c_int,
                c_bvID,
            ]

            self.monosat_c.bv_not.argtypes = [c_solver_p, c_bv_p, c_bvID, c_bvID]
            self.monosat_c.bv_and.argtypes = [
                c_solver_p,
                c_bv_p,
                c_bvID,
                c_bvID,
                c_bvID,
            ]
            self.monosat_c.bv_nand.argtypes = [
                c_solver_p,
                c_bv_p,
                c_bvID,
                c_bvID,
                c_bvID,
            ]
            self.monosat_c.bv_or.argtypes = [c_solver_p, c_bv_p, c_bvID, c_bvID, c_bvID]
            self.monosat_c.bv_nor.argtypes = [
                c_solver_p,
                c_bv_p,
                c_bvID,
                c_bvID,
                c_bvID,
            ]
            self.monosat_c.bv_xor.argtypes = [
                c_solver_p,
                c_bv_p,
                c_bvID,
                c_bvID,
                c_bvID,
            ]
            self.monosat_c.bv_xnor.argtypes = [
                c_solver_p,
                c_bv_p,
                c_bvID,
                c_bvID,
                c_bvID,
            ]

            self.monosat_c.bv_concat.argtypes = [
                c_solver_p,
                c_bv_p,
                c_bvID,
                c_bvID,
                c_bvID,
            ]
            self.monosat_c.bv_popcount.argtypes = [
                c_solver_p,
                c_bv_p,
                c_literal_p,
                c_int,
                c_bvID,
            ]
            self.monosat_c.bv_unary.argtypes = [
                c_solver_p,
                c_bv_p,
                c_literal_p,
                c_int,
                c_bvID,
            ]

            self.monosat_c.bv_bitblast.argtypes = [c_solver_p, c_bv_p, c_bvID]

            self.monosat_c.bv_slice.argtypes = [
                c_solver_p,
                c_bv_p,
                c_bvID,
                c_int,
                c_int,
                c_bvID,
            ]

            self.monosat_c.newGraph.argtypes = [c_solver_p]
            self.monosat_c.newGraph.restype = c_graph_p

            self.monosat_c.newNode.argtypes = [c_solver_p, c_graph_p]
            self.monosat_c.newNode.restype = c_int

            self.monosat_c.nNodes.argtypes = [c_solver_p, c_graph_p]
            self.monosat_c.nNodes.restype = c_int

            self.monosat_c.nEdges.argtypes = [c_solver_p, c_graph_p]
            self.monosat_c.nEdges.restype = c_int

            self.monosat_c.newEdge.argtypes = [
                c_solver_p,
                c_graph_p,
                c_int,
                c_int,
                c_int64,
            ]
            self.monosat_c.newEdge.restype = c_literal

            self.monosat_c.newEdge_double.argtypes = [
                c_solver_p,
                c_graph_p,
                c_int,
                c_int,
                c_double,
            ]
            self.monosat_c.newEdge_double.restype = c_literal

            self.monosat_c.newEdge_bv.argtypes = [
                c_solver_p,
                c_graph_p,
                c_int,
                c_int,
                c_bvID,
            ]
            self.monosat_c.newEdge_bv.restype = c_literal

            self.monosat_c.reaches.argtypes = [c_solver_p, c_graph_p, c_int, c_int]
            self.monosat_c.reaches.restype = c_literal

            self.monosat_c.reachesBackward.argtypes = [
                c_solver_p,
                c_graph_p,
                c_int,
                c_int,
            ]
            self.monosat_c.reachesBackward.restype = c_literal

            self.monosat_c.onPath.argtypes = [
                c_solver_p,
                c_graph_p,
                c_int,
                c_int,
                c_int,
            ]
            self.monosat_c.onPath.restype = c_literal

            self.monosat_c.shortestPathUnweighted_lt_const.argtypes = [
                c_solver_p,
                c_graph_p,
                c_int,
                c_int,
                c_int64,
            ]
            self.monosat_c.shortestPathUnweighted_lt_const.restype = c_literal

            self.monosat_c.shortestPathUnweighted_leq_const.argtypes = [
                c_solver_p,
                c_graph_p,
                c_int,
                c_int,
                c_int64,
            ]
            self.monosat_c.shortestPathUnweighted_leq_const.restype = c_literal

            self.monosat_c.shortestPath_lt_const.argtypes = [
                c_solver_p,
                c_graph_p,
                c_int,
                c_int,
                c_int64,
            ]
            self.monosat_c.shortestPath_lt_const.restype = c_literal

            self.monosat_c.shortestPath_leq_const.argtypes = [
                c_solver_p,
                c_graph_p,
                c_int,
                c_int,
                c_int64,
            ]
            self.monosat_c.shortestPath_leq_const.restype = c_literal

            self.monosat_c.shortestPath_lt_bv.argtypes = [
                c_solver_p,
                c_graph_p,
                c_int,
                c_int,
                c_bvID,
            ]
            self.monosat_c.shortestPath_lt_bv.restype = c_literal

            self.monosat_c.shortestPath_leq_bv.argtypes = [
                c_solver_p,
                c_graph_p,
                c_int,
                c_int,
                c_bvID,
            ]
            self.monosat_c.shortestPath_leq_bv.restype = c_literal

            self.monosat_c.nVars.argtypes = [c_solver_p]
            self.monosat_c.nVars.restype = c_int

            self.monosat_c.nClauses.argtypes = [c_solver_p]
            self.monosat_c.nClauses.restype = c_int

            self.monosat_c.maximumFlow_geq.argtypes = [
                c_solver_p,
                c_graph_p,
                c_int,
                c_int,
                c_int64,
            ]
            self.monosat_c.maximumFlow_geq.restype = c_literal

            self.monosat_c.maximumFlow_gt.argtypes = [
                c_solver_p,
                c_graph_p,
                c_int,
                c_int,
                c_int64,
            ]
            self.monosat_c.maximumFlow_gt.restype = c_literal

            self.monosat_c.maximumFlow_geq_bv.argtypes = [
                c_solver_p,
                c_graph_p,
                c_int,
                c_int,
                c_bvID,
            ]
            self.monosat_c.maximumFlow_geq_bv.restype = c_literal

            self.monosat_c.maximumFlow_gt_bv.argtypes = [
                c_solver_p,
                c_graph_p,
                c_int,
                c_int,
                c_bvID,
            ]
            self.monosat_c.maximumFlow_gt_bv.restype = c_literal

            self.monosat_c.minimumSpanningTree_leq.argtypes = [
                c_solver_p,
                c_graph_p,
                c_int64,
            ]
            self.monosat_c.minimumSpanningTree_leq.restype = c_literal

            self.monosat_c.minimumSpanningTree_lt.argtypes = [
                c_solver_p,
                c_graph_p,
                c_int64,
            ]
            self.monosat_c.minimumSpanningTree_lt.restype = c_literal

            self.monosat_c.acyclic_undirected.argtypes = [c_solver_p, c_graph_p]
            self.monosat_c.acyclic_undirected.restype = c_literal

            self.monosat_c.acyclic_directed.argtypes = [c_solver_p, c_graph_p]
            self.monosat_c.acyclic_directed.restype = c_literal

            self.monosat_c.connectedComponents_geq_const.argtypes = [
                c_solver_p,
                c_graph_p,
                c_int
            ]
            self.monosat_c.connectedComponents_geq_const.restype = c_literal

            self.monosat_c.newEdgeSet.argtypes = [
                c_solver_p,
                c_graph_p,
                c_literal_p,
                c_int,
                c_bool,
            ]

            self.monosat_c.initFSMTheory.argtypes = [c_solver_p]
            self.monosat_c.initFSMTheory.restype = c_fsm_theory_p

            self.monosat_c.newFSM.argtypes = [c_solver_p, c_fsm_theory_p, c_int, c_int]
            self.monosat_c.newFSM.restype = c_int

            self.monosat_c.newState.argtypes = [c_solver_p, c_fsm_theory_p, c_int]
            self.monosat_c.newState.restype = c_int

            self.monosat_c.newTransition.argtypes = [
                c_solver_p,
                c_fsm_theory_p,
                c_int,
                c_int,
                c_int,
                c_int,
                c_int,
            ]
            self.monosat_c.newTransition.restype = c_int

            self.monosat_c.newString.argtypes = [
                c_solver_p,
                c_fsm_theory_p,
                c_int_p,
                c_int,
            ]
            self.monosat_c.newString.restype = c_int

            self.monosat_c.fsmAcceptsString.argtypes = [
                c_solver_p,
                c_fsm_theory_p,
                c_int,
                c_int,
                c_int,
                c_int,
            ]
            self.monosat_c.fsmAcceptsString.restype = c_int

            self.monosat_c.fsmCompositionAccepts.argtypes = [
                c_solver_p,
                c_fsm_theory_p,
                c_int,
                c_int,
                c_int,
                c_int,
                c_int,
                c_int,
                c_int
            ]
            self.monosat_c.fsmCompositionAccepts.restype = c_int

            self.monosat_c.getModel_Literal.argtypes = [c_solver_p, c_literal]
            self.monosat_c.getModel_Literal.restype = c_int

            self.monosat_c.getModel_BV.argtypes = [c_solver_p, c_bv_p, c_bvID, c_bool]
            self.monosat_c.getModel_BV.restype = c_int64

            self.monosat_c.getModel_MaxFlow.argtypes = [
                c_solver_p,
                c_graph_p,
                c_literal,
            ]
            self.monosat_c.getModel_MaxFlow.restype = c_int64

            self.monosat_c.getModel_EdgeFlow.argtypes = [
                c_solver_p,
                c_graph_p,
                c_literal,
                c_literal,
            ]
            self.monosat_c.getModel_EdgeFlow.restype = c_int64

            self.monosat_c.getModel_AcyclicEdgeFlow.argtypes = [
                c_solver_p,
                c_graph_p,
                c_literal,
                c_literal,
            ]
            self.monosat_c.getModel_AcyclicEdgeFlow.restype = c_int64

            self.monosat_c.getModel_MinimumSpanningTreeWeight.argtypes = [
                c_solver_p,
                c_graph_p,
                c_literal,
            ]
            self.monosat_c.getModel_MinimumSpanningTreeWeight.restype = c_int64

            self.monosat_c.getModel_Path_Nodes_Length.argtypes = [
                c_solver_p,
                c_graph_p,
                c_literal,
            ]
            self.monosat_c.getModel_Path_Nodes_Length.restype = c_int

            self.monosat_c.getModel_Path_Nodes.argtypes = [
                c_solver_p,
                c_graph_p,
                c_literal,
                c_int,
                c_int_p,
            ]
            self.monosat_c.getModel_Path_Nodes.restype = c_int

            self.monosat_c.getModel_Path_EdgeLits_Length.argtypes = [
                c_solver_p,
                c_graph_p,
                c_literal,
            ]
            self.monosat_c.getModel_Path_EdgeLits_Length.restype = c_int

            self.monosat_c.getModel_Path_EdgeLits.argtypes = [
                c_solver_p,
                c_graph_p,
                c_literal,
                c_int,
                c_int_p,
            ]
            self.monosat_c.getModel_Path_EdgeLits.restype = c_int

            self.monosat_c.createFlowRouting.argtypes = [
                c_solver_p,
                c_graph_p,
                c_int,
                c_int,
                c_literal,
            ]
            self.monosat_c.createFlowRouting.restype = c_int64

            self.monosat_c.addRoutingNet.argtypes = [
                c_solver_p,
                c_graph_p,
                c_int64,
                c_literal,
                c_int,
                c_literal_p,
                c_literal_p,
            ]

            self.monosat_c.graph_setAssignEdgesToWeight.argtypes = [
                c_solver_p,
                c_graph_p,
                c_int64,
            ]

        self.solver = None
        self.newSolver()
        # For many (but not all) instances, the following settings may give good performance:
        # self.newSolver("-verb=0 -verb-time=0 -rnd-theory-freq=0.99 -no-decide-bv-intrinsic  -decide-bv-bitwise  -decide-graph-bv -decide-theories -no-decide-graph-rnd   -lazy-maxflow-decisions -conflict-min-cut -conflict-min-cut-maxflow -reach-underapprox-cnf -check-solution ")

    def _getManagers(self):
        solver = self.getSolver()
        if solver not in self._managers:
            self._managers[solver] = dict()
        return self._managers[solver]

    def getSolver(self):
        return self.solver

    def setSolver(self, solver):
        self.solver = solver

    def getVersion(self):
        return self.monosat_c.getVersion()

    def init(self, arguments=None, output_file=None):
        warnings.warn(
            "Monosat().init() is deprecated, use Monosat().newSolver() instead",
            DeprecationWarning,
        )
        return self.newSolver(arguments, output_file)

    def setOutputFile(self, output_file=None):
        warnings.warn(
            "Monosat().setOutputFile() is deprecated, use Monosat().newSolver(output_file=filename) instead",
            DeprecationWarning,
        )
        self.newSolver(arguments=self.solver.arguments, output_file=output_file)

    def newSolver(self, arguments=None, output_file=None, define_true=True):
        if self.solver is not None:
            self.solver.delete()
        self.solver = Solver(self.monosat_c, arguments, output_file)
        if define_true:
            self.solver._true = self.getTrue()
        return self.solver

    def loadConstraints(self, filename, process_solve_statements=False):
        if process_solve_statements:
            self.monosat_c.readGNF(self.solver._ptr, c_char_p(filename.encode("ascii")))
        else:
            self.monosat_c.loadGNF(self.solver._ptr, c_char_p(filename.encode("ascii")))

    def getEmptyIntArray(self, length):
        if length > len(self._int_array):
            self._int_array = (c_int * length)()
        return self._int_array

    def getIntArray(self, nums):
        if use_cython:
            return nums
        if len(nums) > len(self._int_array):
            self._int_array = (c_int * len(nums))()
        for i, n in enumerate(nums):
            self._int_array[i] = c_int(n)
        return self._int_array

    def getIntArray2(self, nums):
        if use_cython:
            return nums
        if len(nums) > len(self._int_array2):
            self._int_array2 = (c_int * len(nums))()
        for i, n in enumerate(nums):
            self._int_array2[i] = c_int(n)
        return self._int_array2

    def getLongArray(self, nums):
        if use_cython:
            return nums
        if len(nums) > len(self._long_array):
            self._long_array = (c_int64 * len(nums))()
        for i, n in enumerate(nums):
            self._long_array[i] = c_int64(n)
        return self._long_array

    def intArrayToList(self, array_pointer, length):
        ret = []
        for i in range(length):
            ret.append(array_pointer[i])
        return ret

    def setTimeLimit(self, seconds):
        if seconds is None or seconds < 0:
            self.monosat_c.setTimeLimit(self.solver._ptr, -1)
        else:
            self.monosat_c.setTimeLimit(self.solver._ptr, seconds)

    def setConflictLimit(self, conflicts):
        if conflicts is None or conflicts < 0:
            self.monosat_c.setConflictLimit(self.solver._ptr, -1)
        else:
            self.monosat_c.setConflictLimit(self.solver._ptr, conflicts)

    def setPropagationLimit(self, propagations):
        if propagations is None or propagations < 0:
            self.monosat_c.setPropagationLimit(self.solver._ptr, -1)
        else:
            self.monosat_c.setPropagationLimit(self.solver._ptr, propagations)

    def lastSolutionWasOptimal(self):
        return self.monosat_c.lastSolutionWasOptimal(self.solver._ptr)

    def getConflictClause(self):
        conflict_size = self.monosat_c.getConflictClause(self.solver._ptr, null_ptr, 0)
        if conflict_size < 0:
            return None  # No conflict clause in the solver

        conflict_ptr = self.getEmptyIntArray(conflict_size)
        length = self.monosat_c.getConflictClause(
            self.solver._ptr, conflict_ptr, conflict_size
        )
        if length != conflict_size:
            raise RuntimeError("Error reading conflict clause")

        return self.intArrayToList(conflict_ptr, conflict_size)

    def minimizeUnsatCore(self, assumptions):
        lp = self.getIntArray(assumptions)

        core_size = self.monosat_c.minimizeUnsatCore(
            self.solver._ptr, lp, len(assumptions)
        )
        assert core_size >= 0
        assert core_size <= len(assumptions)
        return self.intArrayToList(lp, core_size)

    def minimizeConflictClause(self):
        self.monosat_c.minimizeConflictClause(self.solver._ptr)

    def solve(self, assumptions=None):
        self.backtrack()
        if assumptions is None:
            assumptions = []

        lp = self.getIntArray(assumptions)

        return self.monosat_c.solveAssumptions(self.solver._ptr, lp, len(assumptions))

    def solveLimited(self, assumptions=None):
        self.backtrack()
        if assumptions is None:
            assumptions = []

        lp = self.getIntArray(assumptions)

        r = self.monosat_c.solveAssumptionsLimited(
            self.solver._ptr, lp, len(assumptions)
        )

        if r == 0:
            return True
        elif r == 1:
            return False
        else:
            assert r == 2
            return None

    def backtrack(self):
        return self.monosat_c.backtrack(self.solver._ptr)

    def addUnitClause(self, clause):

        if isinstance(clause, int):
            self.monosat_c.addUnitClause(self.solver._ptr, clause)
        elif len(clause) == 1:
            self.monosat_c.addUnitClause(self.solver._ptr, clause[0])
        else:
            assert False

    def addBinaryClause(self, l0, l1):
        self.monosat_c.addBinaryClause(self.solver._ptr, l0, l1)

    # clauses should be a list of pairs of literals, each of which will be added as a binary clause
    def addBinaryClauses(self, clauses):

        if len(clauses) == 1:
            self.addBinaryClause(clauses[0][0], clauses[0][1])
        else:
            firsts, seconds = zip(*clauses)

            first_p = self.getIntArray(firsts)
            second_p = self.getIntArray2(seconds)
            self.monosat_c.addBinaryClauses(
                self.solver._ptr, first_p, second_p, len(clauses)
            )

    def addTertiaryClause(self, l0, l1, l2):
        self.monosat_c.addTertiaryClause(self.solver._ptr, l0, l1, l2)

    def addClause(self, clause):
        if isinstance(clause, int):
            self.addUnitClause(clause)
        elif len(clause) == 1:
            self.addUnitClause(clause[0])
        elif len(clause) == 2:
            self.addBinaryClause(clause[0], clause[1])
        elif len(clause) == 3:
            self.addTertiaryClause(clause[0], clause[1], clause[2])
        else:
            lp = self.getIntArray(clause)
            self.monosat_c.addClause(self.solver._ptr, lp, len(clause))

    def clearOptimizationObjectives(self):
        self.monosat_c.clearOptimizationObjectives(self.solver._ptr)

    def maximizeBV(self, bvID):
        self.monosat_c.maximizeBV(self.solver._ptr, self.solver.bvtheory, c_int(bvID))

    def minimizeBV(self, bvID):
        self.monosat_c.minimizeBV(self.solver._ptr, self.solver.bvtheory, c_int(bvID))

    def maximizeLits(self, lits):
        lp = self.getIntArray(lits)
        self.monosat_c.maximizeLits(self.solver._ptr, lp, len(lits))

    def minimizeLits(self, lits):
        lp = self.getIntArray(lits)
        self.monosat_c.minimizeLits(self.solver._ptr, lp, len(lits))

    def maximizeWeightedLits(self, lits, weights):
        lp = self.getIntArray(lits)
        lp2 = self.getIntArray2(weights)
        self.monosat_c.maximizeWeightedLits(self.solver._ptr, lp, lp2, len(lits))

    def minimizeWeightedLits(self, lits, weights):
        lp = self.getIntArray(lits)
        lp2 = self.getIntArray2(weights)
        self.monosat_c.minimizeWeightedLits(self.solver._ptr, lp, lp2, len(lits))

    # convenience code to and together to lits and return a new
    def addAnd(self, lit1, lit2):
        out = self.newLit()
        self.addTertiaryClause(out, self.Not(lit1), self.Not(lit2))
        self.addBinaryClause(self.Not(out), lit1)
        self.addBinaryClause(self.Not(out), lit2)
        return out

    def newLit(self, allow_simplification=False):
        varID = int(self.monosat_c.newVar(self.solver._ptr))
        if not allow_simplification:
            self.monosat_c.disallowLiteralSimplification(self.solver._ptr, varID * 2)
        return varID * 2

    def setDecisionVar(self, v, b):
        self.monosat_c.setDecisionVar(self.solver._ptr, c_var(v), c_bool(b))

    def isDecisionVar(self, var):
        return self.monosat_c.isDecisionVar(self.solver._ptr, c_var(v))

    def setDecisionPriority(self, v, priority):
        self.monosat_c.setDecisionPriority(self.solver._ptr, c_var(v), c_int(priority))

    def getDecisionPriority(self, v):
        self.monosat_c.getDecisionPriority(self.solver._ptr, c_var(v))

    def setDecisionPolarity(self, v, b):
        self.monosat_c.setDecisionPolarity(self.solver._ptr, c_var(v), c_bool(b))

    def getDecisionPolarity(self, var):
        return self.monosat_c.getDecisionPolarity(self.solver._ptr, c_var(v))

    def nVars(self):
        return self.monosat_c.nVars(self.solver._ptr)

    def nClauses(self):
        return self.monosat_c.nClauses(self.solver._ptr)

    def true(self):
        return self.solver._true

    def false(self):
        return self.Not(self.solver._true)

    def Not(self, lit):
        return lit ^ 1

    def isPositive(self, lit):
        return (~lit) & 1

    def isNegative(self, lit):
        return lit & 1

    def setSymbol(self, lit, symbol):
        pass

    def getSymbol(self, lit):
        pass

    def AssertAtMostOne(self, clause):
        self.backtrack()
        newclause = []
        for l in clause:
            if self.isPositive(l):
                newclause.append(l)
            else:
                # Create a new, positive literal and assert that it is equal to the old, negative one.
                l2 = self.newLit()
                self.addBinaryClause(self.Not(l), l2)
                self.addBinaryClause(self.Not(l2), l)
                newclause.append(l2)

        lp = self.getIntArray(newclause)
        for i in range(len(newclause)):
            l = lp[i]
            assert self.isPositive(l)
            assert l % 2 == 0  # all of the literals must be positive
            lp[i] = l // 2

        self.monosat_c.at_most_one(self.solver._ptr, lp, len(newclause))

    def pbOpToStr(self, op):
        if op == Ineq.LT:
            return "<"
        elif op == Ineq.LEQ:
            return "<="
        elif op == Ineq.EQ:
            return "=="
        elif op == Ineq.GEQ:
            return ">="
        elif op == Ineq.GT:
            return ">"

    def AssertPB(self, lits, coefs, op, rhs):
        self.backtrack()
        lp = self.getIntArray(lits)
        lp2 = self.getIntArray2(coefs)
        crhs = c_int(rhs)
        assert len(lits) == len(coefs)
        if op == Ineq.LT:
            self.monosat_c.assertPB_lt(self.solver._ptr, crhs, len(lits), lp, lp2)
        elif op == Ineq.LEQ:
            self.monosat_c.assertPB_leq(self.solver._ptr, crhs, len(lits), lp, lp2)
        elif op == Ineq.EQ:
            self.monosat_c.assertPB_eq(self.solver._ptr, crhs, len(lits), lp, lp2)
        elif op == Ineq.GEQ:
            self.monosat_c.assertPB_geq(self.solver._ptr, crhs, len(lits), lp, lp2)
        elif op == Ineq.GT:
            self.monosat_c.assertPB_gt(self.solver._ptr, crhs, len(lits), lp, lp2)

    # Convert any psuedo-boolean constraints into cnf in the solver (will be called automatically during solving)
    def flushPB(self):
        self.monosat_c.flushPB(self.solver._ptr)

    # def preprocess(self,disable_future_preprocessing=False):
    #    self.monosat_c.preprocess(disable_future_preprocessing)

    # bv interface
    def newBitvector_anon(self, width):
        self.backtrack()
        bvID = self.monosat_c.newBitvector_anon(
            self.solver._ptr, self.solver.bvtheory, width
        )
        return bvID

    def newBitvector_const(self, width, val):
        self.backtrack()
        bvID = self.monosat_c.newBitvector_const(
            self.solver._ptr, self.solver.bvtheory, width, c_int64(val)
        )
        return bvID

    def nBitvectors(self):
        return self.monosat_c.nBitvectors(self.solver._ptr, self.solver.bvtheory)

    def newBitvector(self, bits):
        self.backtrack()
        arr = self.getIntArray(bits)
        bvID = self.monosat_c.newBitvector(
            self.solver._ptr, self.solver.bvtheory, arr, c_int(len(bits))
        )
        self._num_bvs = bvID + 1
        return bvID

    def getTrue(self):
        self.backtrack()
        if self.solver._true is None:
            l = self.monosat_c.true_lit(self.solver._ptr)
            self.addUnitClause(
                l
            )  # This isn't technically required by the solver; only doing it to ensure the unit clause gets recorded in the GNF...
            self.solver._true = l
        return self.solver._true

    def newBVComparison_const_lt(self, bvID, val):
        self.backtrack()
        l = self.monosat_c.newBVComparison_const_lt(
            self.solver._ptr, self.solver.bvtheory, c_bvID(bvID), c_int64(val)
        )
        return l

    def newBVComparison_bv_lt(self, bvID1, bvID2):
        self.backtrack()
        l = self.monosat_c.newBVComparison_bv_lt(
            self.solver._ptr, self.solver.bvtheory, c_bvID(bvID1), c_bvID(bvID2)
        )
        return l

    def newBVComparison_const_leq(self, bvID, val):
        self.backtrack()
        l = self.monosat_c.newBVComparison_const_leq(
            self.solver._ptr, self.solver.bvtheory, c_bvID(bvID), c_int64(val)
        )
        return l

    def newBVComparison_bv_leq(self, bvID1, bvID2):
        self.backtrack()
        l = self.monosat_c.newBVComparison_bv_leq(
            self.solver._ptr, self.solver.bvtheory, c_bvID(bvID1), c_bvID(bvID2)
        )
        return l

    def newBVComparison_const_gt(self, bvID, val):
        self.backtrack()
        l = self.monosat_c.newBVComparison_const_gt(
            self.solver._ptr, self.solver.bvtheory, c_bvID(bvID), c_int64(val)
        )
        return l

    def newBVComparison_bv_gt(self, bvID1, bvID2):
        self.backtrack()
        l = self.monosat_c.newBVComparison_bv_gt(
            self.solver._ptr, self.solver.bvtheory, c_bvID(bvID1), c_bvID(bvID2)
        )
        return l

    def newBVComparison_const_geq(self, bvID, val):
        self.backtrack()
        l = self.monosat_c.newBVComparison_const_geq(
            self.solver._ptr, self.solver.bvtheory, c_bvID(bvID), c_int64(val)
        )
        return l

    def newBVComparison_bv_geq(self, bvID1, bvID2):
        self.backtrack()
        l = self.monosat_c.newBVComparison_bv_geq(
            self.solver._ptr, self.solver.bvtheory, c_bvID(bvID1), c_bvID(bvID2)
        )
        return l

    def bv_addition(self, aID, bID, resultID):
        self.backtrack()
        self.monosat_c.bv_addition(
            self.solver._ptr,
            self.solver.bvtheory,
            c_bvID(aID),
            c_bvID(bID),
            c_bvID(resultID),
        )

    def bv_subtraction(self, aID, bID, resultID):
        self.backtrack()
        self.monosat_c.bv_subtraction(
            self.solver._ptr,
            self.solver.bvtheory,
            c_bvID(aID),
            c_bvID(bID),
            c_bvID(resultID),
        )

    def bv_multiply(self, aID, bID, resultID):
        self.backtrack()
        self.monosat_c.bv_multiply(
            self.solver._ptr,
            self.solver.bvtheory,
            c_bvID(aID),
            c_bvID(bID),
            c_bvID(resultID),
        )

    def bv_divide(self, aID, bID, resultID):
        self.backtrack()
        self.monosat_c.bv_divide(
            self.solver._ptr,
            self.solver.bvtheory,
            c_bvID(aID),
            c_bvID(bID),
            c_bvID(resultID),
        )

    def bv_ite(self, condition_lit, thnID, elsID, resultID):
        self.backtrack()
        self.monosat_c.bv_ite(
            self.solver._ptr,
            self.solver.bvtheory,
            condition_lit,
            c_bvID(thnID),
            c_bvID(elsID),
            c_bvID(resultID),
        )

    def bv_min(self, args, resultID):
        self.backtrack()
        lp = self.getIntArray(args)
        self.monosat_c.bv_min(
            self.solver._ptr, self.solver.bvtheory, lp, len(args), c_bvID(resultID)
        )

    def bv_max(self, args, resultID):
        self.backtrack()
        lp = self.getIntArray(args)
        self.monosat_c.bv_max(
            self.solver._ptr, self.solver.bvtheory, lp, len(args), c_bvID(resultID)
        )

    def bv_not(self, aID, resultID):
        self.backtrack()
        self.monosat_c.bv_not(
            self.solver._ptr, self.solver.bvtheory, c_bvID(aID), c_bvID(resultID)
        )

    def bv_and(self, aID, bID, resultID):
        self.backtrack()
        self.monosat_c.bv_and(
            self.solver._ptr,
            self.solver.bvtheory,
            c_bvID(aID),
            c_bvID(bID),
            c_bvID(resultID),
        )

    def bv_nand(self, aID, bID, resultID):
        self.backtrack()
        self.monosat_c.bv_nand(
            self.solver._ptr,
            self.solver.bvtheory,
            c_bvID(aID),
            c_bvID(bID),
            c_bvID(resultID),
        )

    def bv_or(self, aID, bID, resultID):
        self.backtrack()
        self.monosat_c.bv_or(
            self.solver._ptr,
            self.solver.bvtheory,
            c_bvID(aID),
            c_bvID(bID),
            c_bvID(resultID),
        )

    def bv_nor(self, aID, bID, resultID):
        self.backtrack()
        self.monosat_c.bv_nor(
            self.solver._ptr,
            self.solver.bvtheory,
            c_bvID(aID),
            c_bvID(bID),
            c_bvID(resultID),
        )

    def bv_xor(self, aID, bID, resultID):
        self.backtrack()
        self.monosat_c.bv_xor(
            self.solver._ptr,
            self.solver.bvtheory,
            c_bvID(aID),
            c_bvID(bID),
            c_bvID(resultID),
        )

    def bv_xnor(self, aID, bID, resultID):
        self.backtrack()
        self.monosat_c.bv_xnor(
            self.solver._ptr,
            self.solver.bvtheory,
            c_bvID(aID),
            c_bvID(bID),
            c_bvID(resultID),
        )

    def bv_concat(self, aID, bID, resultID):
        self.backtrack()
        self.monosat_c.bv_concat(
            self.solver._ptr,
            self.solver.bvtheory,
            c_bvID(aID),
            c_bvID(bID),
            c_bvID(resultID),
        )

    def bv_bitblast(self, bvID):
        self.backtrack()
        self.monosat_c.bv_bitblast(self.solver._ptr, self.solver.bvtheory, c_bvID(bvID))

    def bv_slice(self, aID, lower, upper, resultID):
        self.backtrack()
        self.monosat_c.bv_slice(
            self.solver._ptr,
            self.solver.bvtheory,
            c_bvID(aID),
            lower,
            upper,
            c_bvID(resultID),
        )

    def bv_popcount(self, args, resultID):
        self.backtrack()
        newargs = []
        for l in args:
            if self.isPositive(l):
                # newargs.append(l)
                l2 = self.newLit()
                self.addBinaryClause(self.Not(l), l2)
                self.addBinaryClause(self.Not(l2), l)
                newargs.append(l2)
            else:
                # Create a new, positive literal and assert that it is equal to the old, negative one.
                l2 = self.newLit()
                self.addBinaryClause(self.Not(l), l2)
                self.addBinaryClause(self.Not(l2), l)
                newargs.append(l2)

        lp = self.getIntArray(newargs)
        self.monosat_c.bv_popcount(
            self.solver._ptr, self.solver.bvtheory, lp, len(newargs), c_bvID(resultID)
        )

    def bv_unary(self, args, resultID):
        self.backtrack()
        newargs = []
        for l in args:
            if self.isPositive(l):
                # newargs.append(l)
                l2 = self.newLit()
                self.addBinaryClause(self.Not(l), l2)
                self.addBinaryClause(self.Not(l2), l)
                newargs.append(l2)
            else:
                # Create a new, positive literal and assert that it is equal to the old, negative one.
                l2 = self.newLit()
                self.addBinaryClause(self.Not(l), l2)
                self.addBinaryClause(self.Not(l2), l)
                newargs.append(l2)

        lp = self.getIntArray(newargs)
        self.monosat_c.bv_unary(
            self.solver._ptr, self.solver.bvtheory, lp, len(newargs), c_bvID(resultID)
        )

    # Monosat fsm interface

    def newFSM(self, in_labels, out_labels):
        return self.monosat_c.newFSM(self.solver._ptr, null_ptr64, in_labels, out_labels)

    def newState(self, fsm_id):
        return self.monosat_c.newState(self.solver._ptr, null_ptr64, fsm_id)

    def newTransition(self, fsm_id, from_state, to_state, in_label, out_label):
        return self.monosat_c.newTransition(self.solver._ptr, null_ptr64, fsm_id, from_state, to_state, in_label, out_label)

    # A string is an array of positive integers
    def newString(self, string_int_array):
        lp = self.getIntArray(string_int_array)
        return self.monosat_c.newString(self.solver._ptr, null_ptr64, lp, len(string_int_array))

    def fsmAcceptsString(self, fsm_id, starting_state, accepting_state, strID):
        return self.monosat_c.fsmAcceptsString(
            self.solver._ptr, null_ptr64, fsm_id, starting_state, accepting_state, strID
        )

    def fsmCompositionAccepts(
            self,
            fsm_generator_id,
            fsm_acceptor_id,
            gen_starting_state,
            gen_accepting_state,
            accept_starting_state,
            accept_accepting_state,
            strID,
    ):
        return self.monosat_c.fsmCompositionAccepts(
            self.solver._ptr,
            null_ptr64,
            fsm_generator_id,
            fsm_acceptor_id,
            gen_starting_state,
            gen_accepting_state,
            accept_starting_state,
            accept_accepting_state,
            strID,
        )

    # Monosat graph interface

    def newGraph(self):
        self.backtrack()
        g = self.monosat_c.newGraph(self.solver._ptr)
        gid = len(self.solver.graphs)
        self.solver.graphs.append(g)

        return g

    def getGraph(self, id):
        return self.solver.graphs[id]

    def nGraphs(self):
        return len(self.solver.graphs)

    def newNode(self, graph):
        self.backtrack()
        return self.monosat_c.newNode(self.solver._ptr, graph)

    def nNodes(self, graph):
        return self.monosat_c.nNodes(self.solver._ptr, graph)

    def nEdges(self, graph):
        return self.monosat_c.nEdges(self.solver._ptr, graph)

    def newEdge(self, graph, u, v, weight):
        # self.backtrack()
        l = self.monosat_c.newEdge(
            self.solver._ptr, graph, c_int(u), c_int(v), c_int64(weight)
        )
        return l

    def newEdge_double(self, graph, u, v, weight):
        self.backtrack()
        l = self.monosat_c.newEdge_double(
            self.solver._ptr, graph, c_int(u), c_int(v), c_double(weight)
        )
        return l

    def newEdge_bv(self, graph, u, v, bvID):
        self.backtrack()
        l = self.monosat_c.newEdge_bv(
            self.solver._ptr, graph, c_int(u), c_int(v), c_bvID(bvID)
        )
        return l

    def newEdgeSet(self, graph, edges, enforceEdgeAssignments=True):
        self.backtrack()
        lp = self.getIntArray(edges)
        self.monosat_c.newEdgeSet(
            self.solver._ptr, graph, lp, len(edges), c_bool(enforceEdgeAssignments)
        )

    def assignWeightsTo(self, graph, weight):
        self.monosat_c.graph_setAssignEdgesToWeight(
            self.solver._ptr, graph, c_int64(weight)
        )

    def enforceRouting(self, graph, source, destination, nets, maxflowlit):
        r_ptr = self.monosat_c.createFlowRouting(
            self.solver._ptr, graph, c_int(source), c_int(destination), maxflowlit
        )
        for dest_edge_lits, net_reach_lits, disabled_edge_lit in nets:
            lp = self.getIntArray(dest_edge_lits)
            lp2 = self.getIntArray2(net_reach_lits)
            self.monosat_c.addRoutingNet(
                self.solver._ptr,
                graph,
                r_ptr,
                disabled_edge_lit,
                len(dest_edge_lits),
                lp,
                lp2,
            )

    def reaches(self, graph, u, v):
        self.checkNode(graph, u)
        self.checkNode(graph, v)
        self.backtrack()
        l = self.monosat_c.reaches(self.solver._ptr, graph, c_int(u), c_int(v))
        return l

    def reachesBackward(self, graph, u, v):
        self.checkNode(graph, u)
        self.checkNode(graph, v)
        self.backtrack()
        l = self.monosat_c.reachesBackward(self.solver._ptr, graph, c_int(u), c_int(v))
        return l

    def onPath(self, graph, nodeOnPath, u, v):
        self.checkNode(graph, u)
        self.checkNode(graph, v)
        self.checkNode(graph, nodeOnPath)
        self.backtrack()
        l = self.monosat_c.onPath(
            self.solver._ptr, graph, c_int(nodeOnPath), c_int(u), c_int(v)
        )
        return l

    def shortestPathUnweighted_lt_const(self, graph, u, v, dist):
        self.checkNode(graph, u)
        self.checkNode(graph, v)
        self.backtrack()
        l = self.monosat_c.shortestPathUnweighted_lt_const(
            self.solver._ptr, graph, c_int(u), c_int(v), c_int64(dist)
        )
        return l

    def shortestPathUnweighted_leq_const(self, graph, u, v, dist):
        self.checkNode(graph, u)
        self.checkNode(graph, v)
        self.backtrack()
        l = self.monosat_c.shortestPathUnweighted_leq_const(
            self.solver._ptr, graph, c_int(u), c_int(v), c_int64(dist)
        )
        return l

    def shortestPath_lt_const(self, graph, u, v, dist):
        self.checkNode(graph, u)
        self.checkNode(graph, v)
        self.backtrack()
        l = self.monosat_c.shortestPath_lt_const(
            self.solver._ptr, graph, c_int(u), c_int(v), c_int64(dist)
        )
        return l

    def shortestPath_leq_const(self, graph, u, v, dist):
        self.checkNode(graph, u)
        self.checkNode(graph, v)
        self.backtrack()
        l = self.monosat_c.shortestPath_leq_const(
            self.solver._ptr, graph, c_int(u), c_int(v), c_int64(dist)
        )
        return l

    def shortestPath_lt_bv(self, graph, u, v, bvID):
        self.checkNode(graph, u)
        self.checkNode(graph, v)
        self.checkBV(bvID)
        self.backtrack()
        l = self.monosat_c.shortestPath_lt_bv(
            self.solver._ptr, graph, c_int(u), c_int(v), c_bvID(bvID)
        )
        return l

    def shortestPath_leq_bv(self, graph, u, v, bvID):
        self.checkNode(graph, u)
        self.checkNode(graph, v)
        self.checkBV(bvID)
        self.backtrack()
        l = self.monosat_c.shortestPath_leq_bv(
            self.solver._ptr, graph, c_int(u), c_int(v), c_bvID(bvID)
        )
        return l

    def maximumFlow_geq(self, graph, s, t, flow):
        self.checkNode(graph, s)
        self.checkNode(graph, t)
        self.backtrack()
        l = self.monosat_c.maximumFlow_geq(
            self.solver._ptr, graph, c_int(s), c_int(t), c_int64(flow)
        )
        return l

    def maximumFlow_gt(self, graph, s, t, flow):
        self.checkNode(graph, s)
        self.checkNode(graph, t)
        self.backtrack()
        l = self.monosat_c.maximumFlow_gt(
            self.solver._ptr, graph, c_int(s), c_int(t), c_int64(flow)
        )
        return l

    def maximumFlow_geq_bv(self, graph, s, t, bvID):
        self.checkNode(graph, s)
        self.checkNode(graph, t)
        self.checkBV(bvID)
        self.backtrack()
        l = self.monosat_c.maximumFlow_geq_bv(
            self.solver._ptr, graph, c_int(s), c_int(t), c_bvID(bvID)
        )
        return l

    def maximumFlow_gt_bv(self, graph, s, t, bvID):
        self.checkNode(graph, s)
        self.checkNode(graph, t)
        self.checkBV(bvID)
        self.backtrack()
        l = self.monosat_c.maximumFlow_gt_bv(
            self.solver._ptr, graph, c_int(s), c_int(t), c_bvID(bvID)
        )
        return l

    def minimumSpanningTree_leq(self, graph, weight):
        self.backtrack()
        l = self.monosat_c.minimumSpanningTree_leq(
            self.solver._ptr, graph, c_int64(weight)
        )
        return l

    def minimumSpanningTree_lt(self, graph, weight):
        self.backtrack()
        l = self.monosat_c.minimumSpanningTree_lt(
            self.solver._ptr, graph, c_int64(weight)
        )
        return l

    def connectedComponents_geq_const(self, graph, components):
        self.backtrack()
        l = self.monosat_c.connectedComponents_geq_const(self.solver._ptr, graph, c_int(components))
        return l

    def acyclic_undirected(self, graph):
        self.backtrack()
        l = self.monosat_c.acyclic_undirected(self.solver._ptr, graph)
        return l

    def acyclic_directed(self, graph):
        self.backtrack()
        l = self.monosat_c.acyclic_directed(self.solver._ptr, graph)
        return l

    # 0 = true, 1=false, 2=unassigned
    def getModel_Literal(self, lit):
        self.checkLit(lit)
        return self.monosat_c.getModel_Literal(self.solver._ptr, lit)

    def getModel_BV(self, bvID, getMaximumValue=False):
        self.checkBV(bvID)
        return self.monosat_c.getModel_BV(
            self.solver._ptr,
            self.solver.bvtheory,
            c_bvID(bvID),
            c_bool(getMaximumValue),
        )

    def getModel_MaxFlow(self, graph, flowlit):
        self.checkLit(flowlit)
        return self.monosat_c.getModel_MaxFlow(self.solver._ptr, graph, flowlit)

    def getModel_EdgeFlow(self, graph, flowlit, edgelit):
        self.checkLit(flowlit)
        self.checkLit(edgelit)
        return self.monosat_c.getModel_EdgeFlow(
            self.solver._ptr, graph, flowlit, edgelit
        )

    def getModel_AcyclicEdgeFlow(self, graph, flowlit, edgelit):
        self.checkLit(flowlit)
        self.checkLit(edgelit)
        return self.monosat_c.getModel_AcyclicEdgeFlow(
            self.solver._ptr, graph, flowlit, edgelit
        )

    def getModel_MinimumSpanningTreeWeight(self, graph, mstlit):
        self.checkLit(mstlit)
        return self.monosat_c.getModel_MinimumSpanningTreeWeight(
            self.solver._ptr, graph, mstlit
        )

    def getModel_Path_Nodes(self, graph, reach_or_distance_lit):
        self.checkLit(reach_or_distance_lit)
        arg_length = self.monosat_c.getModel_Path_Nodes_Length(
            self.solver._ptr, graph, reach_or_distance_lit
        )
        if arg_length < 0:
            return None
        elif arg_length == 0:
            return []
        if use_cython:
            path = []
            l = self.monosat_c.getModel_Path_Nodes(
                self.solver._ptr, graph, reach_or_distance_lit, arg_length, path
            )
            if l != arg_length:
                raise RuntimeError("Error reading path model")
            return path
        else:
            path = list(range(arg_length))
            path_pointer = self.getIntArray(path)
            l = self.monosat_c.getModel_Path_Nodes(
                self.solver._ptr, graph, reach_or_distance_lit, arg_length, path_pointer
            )
            if l != arg_length:
                raise RuntimeError("Error reading path model")

            return self.intArrayToList(path_pointer, arg_length)

    def getModel_Path_EdgeLits(self, graph, reach_or_distance_lit):
        self.checkLit(reach_or_distance_lit)
        arg_length = self.monosat_c.getModel_Path_EdgeLits_Length(
            self.solver._ptr, graph, reach_or_distance_lit
        )
        if arg_length < 0:
            return None
        elif arg_length == 0:
            return []
        if use_cython:
            path = []
            l = self.monosat_c.getModel_Path_EdgeLits(
                self.solver._ptr, graph, reach_or_distance_lit, arg_length, path
            )
            if l != arg_length:
                raise RuntimeError("Error reading path model")
            return path
        else:
            path = list(range(arg_length))
            path_pointer = self.getIntArray(path)
            l = self.monosat_c.getModel_Path_EdgeLits(
                self.solver._ptr, graph, reach_or_distance_lit, arg_length, path_pointer
            )
            if l != arg_length:
                raise RuntimeError("Error reading path model")

            return self.intArrayToList(path_pointer, arg_length)

    def checkLit(self, l):
        if l < 0:
            raise RuntimeError("Bad literal %d" % (l))
        v = l // 2
        if v >= self.nVars():
            raise RuntimeError(
                "Bad variable %d (largest variable is %d)" % (v, self.nVars() - 1)
            )

    def checkBV(self, bvID):
        if bvID < 0 or bvID >= self.nBitvectors():
            raise RuntimeError("Bitvector %d does not exist" % (bvID))

    def checkNode(self, g, s):

        if s < 0 or s >= self.nNodes(g):
            raise RuntimeError(
                "Node %d does not exist in this graph (the largest node is %d)"
                % (s, self.nNodes(g) - 1)
            )
