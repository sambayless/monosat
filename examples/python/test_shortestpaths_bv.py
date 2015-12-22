from monosat import *
g= Graph()
bv1 = BitVector(4)
bv2 = BitVector(4)
bv3 = BitVector(4)
e1 = g.addEdge(1,2,bv1) 
e2 = g.addEdge(2,3,bv2)

Assert(g.distance_leq(1,3,bv3))
Assert(Not(g.distance_lt(1,3,bv3)))

Assert((bv1 + bv3)==9)

result = Solve()
if(result):
  print(result)
  print(e1.value())
  print(e2.value())
  print(bv1.value())
  print(bv2.value())
  print(bv3.value())