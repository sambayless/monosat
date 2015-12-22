    
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





g = Graph()

nodes = []
for n in range(4):
    nodes.append(g.addNode())
    

a = g.addEdge(0,1)
b = g.addEdge(1,2)
c = g.addEdge(2,3)
d = g.addEdge(0,3)
AssertExactlyOne((a,b,c,d))
g.newEdgeSet((a,b,c,d))

e = Var()
Assert(e==Not(c))


mf = g.maxFlowGreaterOrEqualTo(0,3,1)
Assert(mf)

result=Solve()
   
print("Result: " + str(result))
if result:
    print(a.value())
    print(b.value())
    print(c.value())
    print(d.value())
    print(g.getMaxFlow(mf))
    for (v,w,var,wt) in g.getAllEdges():
        #if  var in model:
        if var.value():   
            print("%d->%d"%(v,w))


print("Done")
