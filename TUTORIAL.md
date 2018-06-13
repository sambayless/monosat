# MonoSAT Python Tutorial

This is a brief introduction to MonoSAT's Z3-inspired Python 3 library, which you can use to
conveniently construct and solve formulas with MonoSAT. 
You can find the corresponding code, along with many other examples, in [tutorial.py].

Before going any further, see the installation instructions for the Python library in [README].
Also, be warned that this library has only been tested with Python 3.3+, and may not work on earlier
versions (in particular, Python 2 may not work at all).

Using MonoSAT in Python is as simple as:
```py
#Import the MonoSAT library
from monosat import *

#Create two Boolean variables:
a = Var() 
b = Var() 
c = Or(a, Not(b)) 

#Add a unit clause to the solver, asserting that variable c must be true
Assert(c)

result = Solve() #Solve the instance in MonoSAT, return either True if the instance is SAT, and False if it is UNSAT
if result:
	print("SAT")
else:
	print("UNSAT")
```
	
You can continue making further assertions, creating new variables, and making incremental calls to the solver:
```py
d = Var()
Assert(Implies(d, Or(a, b)))
Assert(d)

result = Solve()
if result:
	print("SAT")
	print("a: " + str(a.value())) 
	print("b: " + str(b.value()))
	print("c: " + str(c.value()))
	print("d: " + str(d.value()))
else:
	print("UNSAT")
```

MonoSAT also supports an alternative syntax using Python's bitwise operators:
```py
Assert(~(a & b))
```

There is no way to remove assertions from MonoSAT yet, however, you can use the assumption mechanism to
temporarily assert that a variable must be true (or false):
```py
result = Solve([b])
```

If in the previous call, MonoSAT was only UNSAT under an assumption, the solver can still be used in subsequent calls:
```py
result = Solve([~b])
```

### Theory Support
 
Now, onto the interesting stuff. 
MonoSAT has support for several useful theories, including some common ones (Bitvectors, Cardinality constraints), 
and some uncommon ones - especially, graph predicates such as reachability, shortest paths, maximum flows, and minimum spanning tree length.
In fact, MonoSAT has support for many more theories from other domains, including finite state machines, but the graph theory is the most well supported, currently.

Constructing a graph in MonoSAT is as easy as:
```py
g = Graph()

#Create three nodes
n1 = g.addNode()
n2 = g.addNode()
n3 = g.addNode()

#Add three directed edges to the graph
e1 = g.addEdge(n1, n2) 
e2 = g.addEdge(n2, n3) 
e3 = g.addEdge(n1, n3)
```

e1, e2, and e3 are *symbolic edges*, meaning that the edge (n1,n2) is included in G if and only if the
theory atom e1 is assigned to True by MonoSAT.
You can use e1,e2, and e3 just like variables in MonoSAT, and in that way control which edges are in the graph
using arbitrary Boolean logic:
```py
Assert(Not(And(e1, e2, e3)))
Assert(Or(e1, e3))

#You can even mix these symbolic edge variables with other logic from MonoSAT
Assert(Implies(c, e1)) 
 
#Once you have created a graph and some edges, you can assert graph properties about that graph:
#For example, you can assert that node n3 must be reachable from node n1, in g
Assert(g.reaches(n1, n3))
```

Graph predicates are 'double sided', so you can also assert that they are false, in order to 
prevent one node from reaching another:
```py
Assert(Not(g.reaches(n2, n1)))

#You can also mix graph predicates in with arbitrary logic, just like variables and edges
Assert(Or(~b, ~g.reaches(n1, n2)))
```

Edges can also have weights, represented as fixed-width, bounded bitvectors.
(By bounded bitvectors, we mean that every bitvector in MonoSAT is asserted to 
be in the range [0, Max], and can never overflow/underflow.)
```py
#create a bitvector of width 4
bv1 = BitVector(4)
bv2 = BitVector(4)
bv3 = BitVector(4)
```

BitVectors support addition and comparisons to constants, but do not yet directly support negatives 
or subtraction (the bitvectors are unsigned):
```py
Assert(bv1+bv2 <= 7)

Assert(bv1 + bv3 >= 3)
Assert(bv1 >= 2)

result = Solve()
if result:
	print("SAT")
	print("bv1: " + str(bv1.value())) 
	print("bv2: " + str(bv2.value()))
	print("bv3: " + str(bv3.value()))
else:
	print("UNSAT")
```

When creating an edge, you can use bitvectors (or Python ints) as edge weights (otherwise, by default, every edge has weight '1'):
```py
#Create a new graph
g2 = Graph()
#Create three nodes
n4 = g2.addNode()
n5 = g2.addNode()
n6 = g2.addNode()

#Add three weighted edges to the graph
e4 = g2.addEdge(n4, n5, bv1) 
e5 = g2.addEdge(n5, n6, bv2) 
e6 = g2.addEdge(n4, n6, bv3)
```

MonoSAT supports several useful graph predicates in addition to reachability, including:
Shortest path constraints:
```py
#Assert that the distance from n1 to n3 is less or equal to 1 (edges have default weights of 1)
Assert(g2.distance_leq(n4, n6, 3)) 

#You can also use BitVectors in the arguments of graph predicates:
bv4 = BitVector(4)
Assert(Not(g2.distance_lt(n4, n6, bv4)))
Assert(bv4 == (bv1 + bv2))
```

MonoSAT also features highly optimized support for maximum flow constraints, allowing for comparisons against either a python integer, or a bitvector:
```py
Assert(g2.maxFlowGreaterOrEqualTo(n4, n6, 3))

bv5 = BitVector(4)
Assert(g2.maxFlowGreaterOrEqualTo(n4, n6, bv5))
```

Just like with reachability and shortest path constraints, these maximum flow predicates are two sided
so you can assert that the maximum flow must be less than a given bitvector, or you can include the
maximum flow predicate as part of arbitrary Boolean logic:
```py
Assert(Or(~c, ~g2.maxFlowGreaterOrEqualTo(n4, n6, bv5 + 1)))
```
	
MonoSAT also features support for minimum spanning tree constraints (in undirected graphs):
```py
g3 = Graph()
n7 = g3.addNode()
n8 = g3.addNode()
n9 = g3.addNode()

#Add three weighted, undirected edges to the graph
e7 = g3.addUndirectedEdge(n7, n8, 1) 
e8 = g3.addUndirectedEdge(n8, n9, 2) 
e9 = g3.addUndirectedEdge(n7, n9, 4)

Assert(g3.minimumSpanningTreeLessEq(3))
Assert(~g3.minimumSpanningTreeLessEq(1))
```

(Watch out, though: minimum spanning tree constraints don't support bitvectors yet.)


[tutorial.py]: examples/python/tutorial.py
[README]: README.md
