from monosat import *    

import functools
import math
import os
import random
import random
import sys


print("begin encode");

seed = random.randint(1,100000)

random.seed(seed)
print("RandomSeed=" + str(seed))


filename="/tmp/test_bv"

filename = filename +  ".gnf"
Monosat().init("-decide-graph-bv -no-decide-theories -no-decide-graph-rnd   -lazy-maxflow-decisions -conflict-min-cut -conflict-min-cut-maxflow -reach-underapprox-cnf ")
Monosat().setOutputFile(open(filename,'w'))


bv0 = BitVector(4)
bv1 = BitVector(4)
bv2 = BitVector(4)
bv3 = BitVector(4)




g = Graph()

nodes = []
for n in range(4):
    nodes.append(g.addNode())
    
Assert(g.addEdge(0,1,bv0))
Assert(g.addEdge(0,2,bv1))
Assert(g.addEdge(1,2,bv2))
Assert(g.addEdge(2,3,bv3))

bv4 = BitVector(4)
bv5 = BitVector(4)

#Assert(g.distance_leq(0,3,bv5))
#Assert(Not(g.distance_lt(0,3,bv5)))

#Assert(bv5<6)
Assert(bv4 >3)

Assert(bv2>1)
Assert(bv3>1)

Assert(bv5<7)

mf = g.maxFlowGreaterOrEqualTo(1,3,bv4)
Assert(mf)
Assert(Not(g.maxFlowGreaterOrEqualTo(1,3,bv5)))

Assert(bv4+bv5==10)
#Assert(bv4+bv5==11)
#Assert(bv6>4)
bvs=[bv0,bv1,bv2,bv3,bv4,bv5]
#Assert(bv6>6)
#Assert(bv6<10)

result=Solve()
   
print(str(result))
if result:
    print(bv1.value())
    print(bv2.value())
    print(bv3.value())
    print(bv4.value())
    print(bv5.value())
    print(g.getMaxFlow(mf))


