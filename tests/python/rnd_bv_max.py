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


print("begin encode");

Monosat().newSolver(filename)

random.seed(seed)
print("RandomSeed=" + str(seed))

width=4

bvs = []
nbvs= 20
nsubsets=5

for i in range(nbvs):
    bvs.append(BitVector(width, random.randint(0,1<<width-1)))

subsets=[]
selecteds=[]
for i in range(nsubsets):
    subset=[]
    selected=[]
    selecteds.append(selected)
    subsets.append(subset)
    for bv in bvs:
        select=Var()
        selected.append(select)
        subset.append(Ite(select,bv,BitVector(width,0))) #1<<width-1
    AssertEqualPB(selected, nbvs//2)
    
maxs = []
for subset in subsets:
    maxs.append(Max(subset))
    
for a,b in itertools.combinations(maxs,2):
    Assert(a!=b)
    
result =Solve()

print("Result is " + str(result))

if result:
    print(" ".join(str(bv.value()) for bv in bvs))
    for subset in selecteds:
        print(" ".join(str(s.value()) for s in subset))
    for m in maxs:
        print(str(m.value()))
        
    sys.exit(10)
else:
    sys.exit(20)
