import random
from monosat import *

print("begin encode");

seed = random.randint(1,100000)

random.seed(seed)
print("RandomSeed=" + str(seed))


a = Var()
b= Var()

bv1 = BitVector(4)
bv2 = BitVector(4)
bv3 = BitVector(4)
bv4 = BitVector(4)

Assert(bv1<=5)
Assert(bv2 >=10)

Assert(bv2==bv3)

AssertIf(a, bv1==bv2, bv3==bv4)
AssertIf(Not(a),  bv1==bv2, bv3==bv4)

result =Solve()

print("Result is " + str(result))
assert(result==False)

