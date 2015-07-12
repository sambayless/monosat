from monosat import *

import functools
import math
import os
import random
import random
import sys


print("Note: currently, pseudo-Boolean encodings rely on Minisat+, which must be installed and on the path (as 'minisat+') in order for MonoSAT to use it.")
vars=[]
for v in range(10):
    vars.append(Var())

AssertEqualPB(vars,2);
AssertLessEqPB(vars, 4)
AssertGreaterThanPB(vars, 1)

result=Solve()
print("Result is " + str(result))
assert(result==True)

