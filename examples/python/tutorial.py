# MonoSAT Python Tutorial

#This is a brief introduction to MonoSAT's Z3-inspired Python 3 library, which you can use to
#conveniently construct and solve formulas with MonoSAT.

#Before going any further, see the installation instructions for the Python library in [README].
#Also, be warned that this library has only been tested with Python 3.3+, and may not work on earlier
#versions (in particular, Python 2 may not work at all).

#Using MonoSAT in Python is as simple as:

#Import the MonoSAT library
from monosat import *

#Create two Boolean variables:
a = Var() 
b = Var()

# c is true if a is true or b is false, and false otherwise
c = Or(a, Not(b))

#Add a unit clause to the solver, asserting that variable c must be true
Assert(c)

#Solve the instance in MonoSAT, returning either True if the instance is SAT, and False if it is UNSAT
result = Solve()
if result:
	print("SAT")
	#After a satisfiable call to Solve(), you can query the assignments given by the solver to
	#individual variables using v.value()
	print("a: " + str(a.value())) 
	print("b: " + str(b.value()))
	print("c: " + str(c.value()))
else:
	print("UNSAT")
	
# After a solve call, you can continue making further assertions, creating new variables,
# and making incremental calls to the solver

d= Var()
Assert(Implies(d, Or(a,b)))

# There are also assertion forms for the common logic constructions, which are slightly more efficient than creating a
# new literal and asserting it to true. An equivalent way to accomplish the above would have been:
AssertImplies(d, Or(a,b))

# Note that d does not yet have an assignment in the solver, and so d.value() will return None until the next solve call
print("Variable 'd' is unassigned, and so has value " + str(d.value()))

result = Solve()
if result:
	print("SAT")
	print("a: " + str(a.value())) 
	print("b: " + str(b.value()))
	print("c: " + str(c.value()))
	print("d: " + str(d.value())) # now d is assigned
else:
	print("UNSAT")

#You can use the '~' operator to apply negation, the same way as Not()
Assert(~(And(a, b)))
result = Solve()
if result:
	print("SAT")
	print("a: " + str(a.value())) 
	print("b: " + str(b.value()))
	print("c: " + str(c.value()))
	print("d: " + str(d.value()))
else:
	print("UNSAT")

#There is no way to remove assertions from MonoSAT yet, however, you can use assumptions to
#temporarily assert that a variable must be true (or false):

result = Solve([b])
if result:
	print("SAT")
	print("a: " + str(a.value())) 
	print("b: " + str(b.value()))
	print("c: " + str(c.value()))
	print("d: " + str(d.value()))
else:
	print("Temporarily UNSAT, under the assumption that 'b' is true")

#If in the previous call, MonoSAT was only UNSAT under an assumption, the solver can still be used in subsequent calls:
result = Solve([~b])
if result:
	print("SAT")
	print("a: " + str(a.value())) 
	print("b: " + str(b.value()))
	print("c: " + str(c.value()))
	print("d: " + str(d.value()))
else:
	print("UNSAT (under the assumption that 'b' is False)")

### Theory Support

# Now, onto the interesting stuff.
# In addition to Boolean logic, MonoSAT supports an extensive theory of finite graphs, including
# support for many common graph predicates such as reachability, shortest paths, maximum flows, acyclicity, and
# minimum spanning trees.
# MonoSAT also has support for BitVectors and Cardinality/Pseudo-Boolean constraints.

#Constructing a graph in MonoSAT is as easy as:
g = Graph()

#Create three nodes
n0 = g.addNode()
n1 = g.addNode()
n2 = g.addNode()

#Add three directed edges to the graph.
#You can also create undirected edges, using g.addUndirectedEdge().
e0 = g.addEdge(n0,n1) 
e1 = g.addEdge(n1,n2) 
e2 = g.addEdge(n0,n2)

#e0, e1, and e2 are *symbolic edges*, meaning that the edge (n0,n1) is included in G if and only if the
#theory atom e0 is assigned to True by MonoSAT.
#You can use e0,e1, and e2 just like variables in MonoSAT, and in that way control which edges are in the graph
#using arbitrary Boolean logic:

AssertNand(e0,e1,e2) # This is logically equivalent to Assert(Not(And(e0,e1,e2)))
AssertOr(e0,e2)

#You can even mix these symbolic edge variables with other logic from MonoSAT
AssertImplies(c, e0)
 
#Once you have created a graph and some edges, you can assert graph properties about that graph.
#For example, you can assert that node n2 must be reachable from node n0, in g
Assert(g.reaches(n0,n2))

result = Solve()
if result:
	print("SAT")
	print("e0: " + str(e0.value())) 
	print("e1: " + str(e1.value()))
	print("e2: " + str(e2.value()))

else:
	print("UNSAT")
	
#Graph predicates are 'double sided', so you can also assert that they are false, in order to 
#prevent one node from reaching another:
Assert(Not(g.reaches(n1,n0)))

#You can also mix graph predicates in with arbitrary logic, just like variables and edges
Assert(Or(~b, ~g.reaches(n0,n1)))

result = Solve()
if result:
	print("SAT")
	print("e0: " + str(e0.value())) 
	print("e1: " + str(e1.value()))
	print("e2: " + str(e2.value()))

else:
	print("UNSAT")

#Edges can also have weights, represented as fixed-width, bounded bitvectors.
#(By bounded bitvectors, we mean that every bitvector in MonoSAT is asserted to 
#be in the range [0, Max], and can never overflow/underflow)

#create a bitvector of width 4
bv0 = BitVector(4)
bv1 = BitVector(4)
bv2 = BitVector(4)


# BitVectors support addition, subtraction, and comparisons, but do not yet directly support
# negative values (the bitvectors are unsigned).
Assert(bv0+bv1 <= 7)

Assert(bv0 + bv2 >= bv1)
Assert(bv0 >= 2)

result = Solve()
if result:
	print("SAT")
	print("bv0: " + str(bv0.value())) 
	print("bv1: " + str(bv1.value()))
	print("bv2: " + str(bv2.value()))

else:
	print("UNSAT")

#When creating an edge, you can use bitvectors (or Python ints) as edge weights (otherwise, by default, every edge has weight '1'):
#Create a new graph
g2 = Graph()
#Create three nodes
n4 = g2.addNode()
n5 = g2.addNode()
n6 = g2.addNode()

#Add three weighted edges to the graph
#Weights may be bitvectors, or integer constants.
e3 = g2.addEdge(n4,n5, bv0) 
e4 = g2.addEdge(n5,n6, bv1) 
e5 = g2.addEdge(n4,n6, bv2)

#MonoSAT supports several useful graph predicates in addition to reachability, including:
#Shortest path constraints:
#Assert that the distance from n0 to n2 is less or equal to 3 (edges have default weights of 1)
Assert(g2.distance_leq(n4,n6,3)) 

#You can also use BitVectors in the arguments of graph predicates:
bv3 = BitVector(4)
Assert(Not(g2.distance_lt(n4,n6,bv3)))
Assert(bv3 == (bv0 + bv1))

result = Solve()
if result:
	print("SAT")
	print("e3: " + str(e3.value())) 
	print("e4: " + str(e4.value()))
	print("e5: " + str(e5.value()))
	print("bv0: " + str(bv0.value())) 
	print("bv1: " + str(bv1.value()))
	print("bv2: " + str(bv2.value()))
	print("bv3: " + str(bv3.value()))
else:
	print("UNSAT")
	
	
#MonoSAT also features highly optimized support for maximum flow constraints, allowing for comparisons against either a python integer, or a bitvector:
Assert(g2.maxFlowGreaterOrEqualTo(n4,n6,3))

bv4 = BitVector(4)
Assert(g2.maxFlowGreaterOrEqualTo(n4,n6,bv4))

#Just like with reachability and distance constraints, these maximum flow predicates are two sided
#so you can assert that the maximum flow must be less than a given bitvector, or you can include the
#maximum flow predicate as part of arbitrary Boolean logic 
Assert(Or(~c,~g2.maxFlowGreaterOrEqualTo(n4,n6,bv4+1)))

result = Solve()
if result:
	print("SAT")
	print("e3: " + str(e3.value())) 
	print("e4: " + str(e4.value()))
	print("e5: " + str(e5.value()))
	print("bv0: " + str(bv0.value())) 
	print("bv1: " + str(bv1.value()))
	print("bv2: " + str(bv2.value()))
	print("bv4: " + str(bv4.value()))
else:
	print("UNSAT")
	
result = Solve([bv4==4])
if result:
	print("SAT")
	print("e3: " + str(e3.value())) 
	print("e4: " + str(e4.value()))
	print("e5: " + str(e5.value()))
	print("bv0: " + str(bv0.value())) 
	print("bv1: " + str(bv1.value()))
	print("bv2: " + str(bv2.value()))
	print("bv4: " + str(bv4.value()))
else:
	print("UNSAT")

result = Solve([bv4>4, bv4<7])
if result:
	print("SAT")
	print("e3: " + str(e3.value())) 
	print("e4: " + str(e4.value()))
	print("e5: " + str(e5.value()))
	print("bv0: " + str(bv0.value())) 
	print("bv1: " + str(bv1.value()))
	print("bv2: " + str(bv2.value()))
	print("bv4: " + str(bv4.value()))
else:
	print("UNSAT")
	
#MonoSAT also features good support for minimum spanning tree constraints (in undirected graphs):
g3 = Graph()
n7 = g3.addNode()
n8 = g3.addNode()
n9 = g3.addNode()

#Add three weighted, undirected edges to the graph
e6 = g3.addUndirectedEdge(n7,n8, 1) 
e7 = g3.addUndirectedEdge(n8,n9, 2) 
e8 = g3.addUndirectedEdge(n7,n9, 4)

Assert(g3.minimumSpanningTreeLessEq(3))
Assert(~g3.minimumSpanningTreeLessEq(1))

result = Solve()
if result:
	print("SAT")
	print("e6: " + str(e6.value())) 
	print("e7: " + str(e7.value()))
	print("e8: " + str(e8.value()))
else:
	print("UNSAT")

#(Minimum spanning tree constraints don't support bitvectors yet, but they could in the future)
