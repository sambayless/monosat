from monosat import *

from networkx.algorithms.dag import is_directed_acyclic_graph
from networkx.relabel import convert_node_labels_to_integers
import functools
import math
import networkx as nx
import networkx.algorithms.components
import os
import random
import random
import sys

print("begin encode");

seed = random.randint(1,100000)

random.seed(seed)
print("RandomSeed=" + str(seed))

size =6





g = Graph()

nodes = []
nxg = nx.gnp_random_graph(9,0.35,127,True)
nxg = convert_node_labels_to_integers(nxg)
for n in nxg:
    nodes.append(g.addNode())
nedges=0

edges =[]
for u,v in nxg.edges():
    print("%d -> %d"%(u,v))
    edges.append(g.addEdge(nodes[u],nodes[v]))



Assert(g.acyclic(True))
Assert(g.reaches(0,8))    
Assert(g.reaches(8,0))    
#This is unsat, because the second two constraints force a cycle, while the first prohibits it

result = Solve()
print("Result is " + str(result))
assert(result==False)

