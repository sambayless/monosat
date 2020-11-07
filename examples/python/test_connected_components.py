from monosat import *

g= Graph()
e1 = g.addEdge(0,1)
e2 = g.addEdge(1,2)
e3 = g.addEdge(0,2)
Assert(Or(e1,e2, e3))

# Test Connected component constraints
# Note that the graph is interpreted as undirected for this constraint - so an edge connecting two nodes, in either
# direction, makes them part of one component.

# Note also that the node indices start at '0'
assert(g.numNodes()==3)

# Because at least one edge above must be true, and there are 3 nodes, there can be at most 2 connected components,
# so the 'g.connectedComponentsGreaterOrEqualTo(3)' cannot be true
atleast3components = g.connectedComponentsGreaterOrEqualTo(3)
assert(not Solve(atleast3components))
# but the negated version should be sat
assert(Solve (Not(atleast3components)))

# cannot have less than 1 connected component

atleast1component = g.connectedComponentsGreaterOrEqualTo(1)
assert(Solve(atleast1component))
assert(not Solve(Not(atleast1component)))

# There can, however, be at least 2 components

atleast2components = g.connectedComponentsGreaterOrEqualTo(2)
assert(Solve(atleast2components))

# and we can also have less than 2 connected components
assert(Solve(Not(atleast2components)))

# But if we assert 2 of the edges false, then we cannot have fewer than 2 connected components
assert(not Solve(Not(atleast2components), Not(e1), Not(e2)))
# If we enforce at least 2 components, this becomes SAT again
assert(Solve(atleast2components, Not(e1), Not(e2)))

# make this constraint permanent
g.AssertConnectedComponentsGreaterOrEqualTo(2)
# Still SAT
assert(Solve())

# Still SAT
g.AssertConnectedComponentsEqualTo(2)
assert(Solve())

# But no longer SAT after this constraint
g.AssertConnectedComponentsLessOrEqualto(1)
assert(not Solve())