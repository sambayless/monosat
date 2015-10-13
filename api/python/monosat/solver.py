from monosat.monosat_c import Monosat
from monosat.logic import *
from monosat.bvtheory import BitVector
from monosat.graphtheory import Graph
from monosat.pbtheory import PBManager
import time



def Solve(assumptions=None, preprocessing=True,bvs_to_minimize=None):
    WriteConstraints()
        
    #if preprocessing:
    #    Monosat().preprocess();
    print("Solving in Monosat...")
    t = time.clock()
    if assumptions is not None or bvs_to_minimize is not None:
        if isinstance(assumptions,Var):
            assumptions=[assumptions]
        elif assumptions is None:
            assumptions = []
            
        if isinstance(bvs_to_minimize,BitVector):
            bvs_to_minimize=[bvs_to_minimize]
        elif bvs_to_minimize is None:
            bvs_to_minimize=[]
        
            
        r= Monosat().solveAssumptions([x.getLit() for x in assumptions],[bv.getID() for bv in bvs_to_minimize])
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