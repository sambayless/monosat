from monosat import *
g= Graph()
e1 = g.addEdge(1,2) 
e2 = g.addEdge(2,3)
e3 = g.addEdge(1,3)
Assert(Not(And(e1,e3)))

Assert(Or(g.reaches(1,3), g.distance_leq(1,3,2)))

result = Solve()
print(result)