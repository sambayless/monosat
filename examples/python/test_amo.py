from monosat import *

import functools
import math
import os
import random
import random
import sys


print("Testing MonoSAT's simple at-most-one constraint theory")
print("For small numbers of variables, try also pseudo-Boolean constraint encodings.")
vars=[]
for v in range(10):
    vars.append(Var())
    
AssertAtMostOne(vars)

result=Solve()
print("Result is " + str(result))
assert(result==True)
for v in vars:
    print(v.value())

AssertOr(vars)
result=Solve()
print("Result is " + str(result))
assert(result==True)
for v in vars:
    print(v.value())
    
Assert(vars[0])
Assert(vars[1])
result=Solve()
print("Result is " + str(result))
assert(result==False)
