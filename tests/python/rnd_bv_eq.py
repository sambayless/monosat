import functools
import math
from monosat import *
import os
from random import shuffle
import random
import random
import sys
import itertools

filename=None
if __name__ == "__main__":
    seed = random.randint(1,100000)

    
    if len(sys.argv)>1:
        filename=sys.argv[1]
    if len(sys.argv)>2:
        seed=int(sys.argv[2])
        

Monosat.newSolver(output_file=filename)
print("begin encode");

random.seed(seed)
print("RandomSeed=" + str(seed))

width=4


bvs = []
nbvs= 10




for i in range(nbvs):
    bvs.append(BitVector(width))

pairs=[ (a,b) for a,b in itertools.combinations(range(nbvs), 2) if a!=b]
neqs = len(pairs)
print(len(pairs))
shuffle(pairs)
eq_pairs =  pairs
print(eq_pairs)

eqs=[]
for a,b in eq_pairs:
    assert(a!=b)
    eqs.append(bvs[a]==bvs[b])

AssertEqualPB(eqs, len(eqs)//2);

vals = []
for bv in bvs:
    n=random.randint(0,1<<width-1)
    vals.append(bv==n)

AssertEqualPB(vals, len(bvs)//2+1)

result =Solve()

print("Result is " + str(result))

if result:
    sys.exit(10)
else:
    sys.exit(20)
