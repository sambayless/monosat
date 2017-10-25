import random
from monosat import *

print("begin encode");
Monosat().setOutputFile("/tmp/bv.gnf")
seed = random.randint(1,100000)

random.seed(seed)
print("RandomSeed=" + str(seed))




bv = BitVector(32)
assert(Solve())
Assert(bv<=0)
Assert(bv>8)

result=Solve()
print("Result is " + str(result))
assert(result==False)

