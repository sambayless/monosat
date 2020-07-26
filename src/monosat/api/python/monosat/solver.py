import time
from monosat.bvtheory import BitVector
from monosat.logic import *
from monosat.monosat_c import Monosat
from monosat.pbtheory import PBManager


def FoundOptimal():
    return Monosat().lastSolutionWasOptimal()


class SolveException(Exception):
    pass


def Solve(
        *assumptions,
        preprocessing=True,
        bvs_to_minimize=None,
        time_limit_seconds=None,
        conflict_limit=None
):
    WriteConstraints()
    if time_limit_seconds is None or time_limit_seconds <= 0:
        time_limit_seconds = -1

    if conflict_limit is None or conflict_limit <= 0:
        conflict_limit = -1

    Monosat().setTimeLimit(time_limit_seconds)

    Monosat().setConflictLimit(conflict_limit)

    # if preprocessing:
    #    Monosat().preprocess();
    # print("Solving in Monosat...")
    t = time.process_time()

    if isinstance(assumptions, Var):
        assumptions = [assumptions]
    elif assumptions is None:
        assumptions = []
    elif len(assumptions) == 1 and not isinstance(assumptions[0], Var):
        try:
            if len(assumptions[0]) != None:
                assumptions = assumptions[0]
        except:
            pass

    if isinstance(bvs_to_minimize, BitVector):
        bvs_to_minimize = [bvs_to_minimize]
    elif bvs_to_minimize is None:
        bvs_to_minimize = []

    for bv in bvs_to_minimize:
        bvID = bv.getID()
        Monosat().minimizeBV(bvID)

    r = Monosat().solveLimited([x.getLit() for x in assumptions])
    if r is None:
        raise SolveException(
            "MonoSAT aborted before solving (possibly do to a time or memory limit)"
        )
    Monosat().elapsed_time += time.process_time() - t
    found_optimal = Monosat().lastSolutionWasOptimal()
    if r is None:
        raise SolveException(
            "MonoSAT aborted before solving (possibly due to a time or memory limit)"
        )
    elif r and not found_optimal:
        print(
            "MonoSAT found a satisfying solution, but it might not be optimal (due to a time or memory limit)"
        )

    return r


# If the most recent solve() call was UNSAT, returns a (negated) subset of the assumptions, at least one of which must be
# true in any satisfying solution.
# If minimize is true, the solver will expend effort trying to find a locally minimla subest
def getConflictClause(minimize=False):
    if minimize:
        Monosat().minimizeConflictClause()
    conf_clause = Monosat().getConflictClause()
    if conf_clause is None:
        return None
    else:
        vars = []
        for v in conf_clause:
            vars.append(Var(v))
        return vars


def minimizeUnsatCore(assumptions):
    assumption_ints = [a.lit for a in assumptions]
    conf_clause = Monosat().minimizeUnsatCore(assumption_ints)
    assert conf_clause is not None
    vars = []
    for v in conf_clause:
        vars.append(Var(v))
    return vars


# optimization support
def clearOptimizationObjectives():
    Monosat().clearOptimizationObjectives()


def maximize(bitvector_or_literals, weights=None):
    if isinstance(bitvector_or_literals, Var):
        arg = [bitvector_or_literals]
    if isinstance(weights, int):
        weights = [weights]

    if isinstance(bitvector_or_literals, BitVector):
        Monosat().maximizeBV(bitvector_or_literals.getID())
    else:
        lit_ints = [l.getLit() for l in bitvector_or_literals]
        if weights is None:
            Monosat().maximizeLits(lit_ints)
        else:
            Monosat().maximizeWeightedLits(lit_ints, weights)


def minimize(bitvector_or_literals, weights=None):
    if isinstance(bitvector_or_literals, Var):
        bitvector_or_literals = [bitvector_or_literals]
    if isinstance(weights, int):
        weights = [weights]

    if isinstance(bitvector_or_literals, BitVector.getID()):
        Monosat().minimizeBV(bitvector_or_literals)
    else:
        lit_ints = [l.getLit() for l in bitvector_or_literals]
        if weights is None:
            Monosat().minimizeLits(lit_ints)
        else:
            Monosat().minimizeWeightedLits(lit_ints, weights)


def WriteConstraints():
    _writePBCosntraints()


def _writePBCosntraints():
    # write any pseudoboolean constraints
    if not PBManager().hasConstraints():
        return

    t = time.process_time()
    pbmgr = PBManager()
    pbmgr.flush()
    d = time.process_time()
    PBManager().elapsed_time += d - t
