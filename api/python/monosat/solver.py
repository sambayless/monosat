from monosat.monosat_c import Monosat
from monosat.logic import *
from monosat.bvtheory import BitVector
from monosat.graphtheory import Graph
from monosat.pbtheory import PBManager
import time
_monosat = Monosat()
_pbm = PBManager()
_pbm.elapsed_time=0
_monosat.elapsed_time=0


def Solve(assumptions=None, preprocessing=True):
    WriteConstraints()
        
    #if preprocessing:
    #    _monosat.preprocess();
    print("Solving in Monosat...")
    t = time.clock()
    if assumptions is not None:
        r= _monosat.solveAssumptions([x.getLit() for x in assumptions])
    else:        
        r= _monosat.solve()
    _monosat.elapsed_time +=  time.clock()-t
    return r

def WriteConstraints():
    
    _writePBCosntraints()

def _writePBCosntraints():
    #write any pseudoboolean constraints
    if not _pbm.hasConstraints():
        return
    

    t = time.clock()
    pbmgr = PBManager()
    pbmgr.flush();
    d = time.clock()
    _pbm.elapsed_time += d-t