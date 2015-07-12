from monosat.monosat_c import Monosat
from monosat.logic import *
from monosat.bvtheory import bv
from monosat.graphtheory import Graph
from monosat.pbtheory import PBManager
import time
_monosat = Monosat()
_pbm = PBManager()
_pbm.elapsed_time=0
_monosat.elapsed_time=0
def solve(assumptions=None, preprocessing=True):
    #first, write any pseudoboolean constraints
    _writePBCosntraints()
        
    #if preprocessing:
    #    _monosat.preprocess();
    print("Solving in Monosat...")
    t = time.clock()
    if assumptions is not None:
        r= _monosat.solveAssumptions((x.getLit() for x in assumptions))
    else:        
        r= _monosat.solve()
    _monosat.elapsed_time +=  time.clock()-t
    return r

def _writePBCosntraints():

    if not _pbm.hasConstraints():
        return
    

    t = time.clock()
    pbmgr = PBManager()
    pbmgr.flush();
    d = time.clock()
    _pbm.elapsed_time += d-t