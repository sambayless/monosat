from monosat import *
g= Graph()
e1 = g.addEdge(1,2)
e2 = g.addEdge(2,3)
e3 = g.addEdge(1,3)
Assert(Not(And(e1,e3)))
Assert(~g.reaches(1,3))
Assert(~g.reaches(1,3))

try:
    Assert(~g.reaches(1,4))
    # should not be reachable, as node 4 doesn't exist, so the above should raise an exception
    assert(false)
except RuntimeError as e:
    pass

try:
    Assert(~g.reaches(1,4))
    # should not be reachable, as node 4 doesn't exist, so the above should raise an exception
    assert(false)
except RuntimeError as e:
    pass

result = Solve()
print(result)


     