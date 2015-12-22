from monosat.monosat_c import Monosat
from monosat.logic import *
from monosat.bvtheory import BitVector
from monosat.graphtheory import Graph
from monosat.pbtheory import PBManager
import time



def Solve(assumptions=None, preprocessing=True,bvs_to_minimize=None,time_limit_seconds=None, memory_limit_mb=None,conflict_limit=None):
    WriteConstraints()
    if time_limit_seconds is None or time_limit_seconds <=0:
        time_limit_seconds=-1
    if memory_limit_mb is None or memory_limit_mb <=0:
        memory_limit_mb=-1    
    if conflict_limit is None or conflict_limit <=0:
        conflict_limit=-1

    Monosat().setTimeLimit(time_limit_seconds)
    Monosat().setMemoryLimit(memory_limit_mb)
    Monosat().setConflictLimit(conflict_limit)
     
    #if preprocessing:
    #    Monosat().preprocess();
    print("Solving in Monosat...")
    t = time.clock()

    if isinstance(assumptions,Var):
        assumptions=[assumptions]
    elif assumptions is None:
        assumptions = []
            
    if isinstance(bvs_to_minimize,BitVector):
        bvs_to_minimize=[bvs_to_minimize]
    elif bvs_to_minimize is None:
        bvs_to_minimize=[]
        
  
    r= Monosat().solveLimited([x.getLit() for x in assumptions],[bv.getID() for bv in bvs_to_minimize])     
    if r is None:
        raise RuntimeError("MonoSAT aborted before solving (possibly do to a time or memory limit)")
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