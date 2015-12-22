from monosat import *

import functools
import math
import os
import random
import random
import sys

print("begin encode");

seed = random.randint(1,100000)

random.seed(seed)
print("RandomSeed=" + str(seed))


bv1 = BitVector(4)
bv2 = BitVector(4)
bv3 = BitVector(4)


Assert(bv1==4)
Assert(bv2 == 3)
Assert(bv3==1)
bv4 = Max(bv1,bv2,bv3)

result =Solve()

print("Result is " + str(result))
assert(result==True)
print(bv1.value())
print(bv2.value())
print(bv3.value())

print(bv4.value())
assert(bv4.value()==4)

