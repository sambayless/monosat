from monosat.monosat_c import Monosat
from monosat.logic import *
from monosat.bvtheory import BitVector
from monosat.graphtheory import Graph
from monosat.pbtheory import PBManager
import time



def Solve(assumptions=None, preprocessing=True):
    WriteConstraints()
        
    #if preprocessing:
    #    Monosat().preprocess();
    print("Solving in Monosat...")
    t = time.clock()
    if assumptions is not None:
        r= Monosat().solveAssumptions([x.getLit() for x in assumptions])
    else:        
        r= Monosat().solve()
    Monosat().elapsed_time +=  time.clock()-t
    return r

def WriteConstraints():
    
    _writePBCosntraints()

def _writePBCosntraints():
    #write any pseudoboolean constraints
    if not PBManager().hasConstraints():
        return
    

    t = time.clock()
    pbmgr = PBManager()
    pbmgr.flush();
    d = time.clock()
    PBManager().elapsed_time += d-t