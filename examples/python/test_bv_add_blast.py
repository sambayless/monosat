import random
from monosat import *

monosat.bvtheory.BVManager().bitblast_addition=True

print("begin encode");

seed = random.randint(1,100000)

random.seed(seed)
print("RandomSeed=" + str(seed))



bv1 = BitVector(4)
bv2 = BitVector(4)
bv3 = BitVector(4)
bv4 = BitVector(4)


g = Graph()

nodes = []
for n in range(4):
    nodes.append(g.addNode())
    
Assert(g.addEdge(0,1,bv1))
Assert(g.addEdge(0,2,bv2))
Assert(g.addEdge(1,2,bv3))
Assert(g.addEdge(2,3,bv4))

bv5 = BitVector(4)
bv6 = BitVector(4)

Assert(g.distance_leq(0,3,bv5))
Assert(Not(g.distance_lt(0,3,bv5)))

Assert(bv5<6)
Assert(bv5 > 4)

#Assert(bv6>=5)

Assert(g.distance_leq(1,3,bv6))
Assert(Not(g.distance_lt(1,3,bv6)))

Assert(bv2==10)
Assert(bv5+bv6<11)
Assert(bv5+bv6>8)
#Assert(bv6>6)
#Assert(bv6<10)

result=Solve()
print("Result is " + str(result))

if result:

    print(bv1.value())
    print(bv2.value())
    print(bv3.value())
    print(bv4.value())
    print(bv5.value())
    print(bv6.value())

print("Done!\n")


