
from monosat import *

vm_bv_width = 4


#Monosat().newSolver("-theory-prop-during-simp -theory-prop-during-fast-simp -opt-conflict-limit=1000 -opt-time-limit=-1 -verb-time=4  -verb=1 -no-rnd-theory-order -theory-order-vsids -decide-opt-lits -binary-search -vsids-both -rnd-theory-freq=0.990000 -vsids-balance=1.000000 -no-decide-bv-intrinsic  -decide-bv-bitwise  -decide-graph-bv -decide-theories -no-decide-graph-rnd   -lazy-maxflow-decisions -no-conflict-min-cut -conflict-min-cut-maxflow -reach-underapprox-cnf -check-solution -no-theory-prop-during-assumps")



used_cores = BitVector(vm_bv_width)
used_ram = BitVector(vm_bv_width)

sum_cores = BitVector(vm_bv_width,0);

fsum = sum_cores +used_cores
Assert(fsum <=4)
fsum.bitblast()


sum_ram = BitVector(vm_bv_width,0);

fsum = sum_ram +used_ram
Assert(fsum <=4)
fsum.bitblast()


if   Solve():
    print("SATISFIABLE")
else:
    print("UNSAT")
