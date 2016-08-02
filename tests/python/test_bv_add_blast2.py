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
bv4 = BitVector(4)



bv5 = BitVector(4)
bv6 = BitVector(4)

Assert(bv5>2)
Assert(bv6 > 4)


s = bv5+bv6
Assert(s<10)
s.bitblast();

#Assert(bv6>6)
#Assert(bv6<10)

result =Solve()

print("Result is " + str(result))
print(bv5.value())
print(bv6.value())
print(s.value());
assert(result==True)

