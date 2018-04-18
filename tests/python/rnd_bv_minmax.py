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
Monosat.newSolver(output_file=filename)
#seed = random.randint(1,100000) # 3538

random.seed(seed)
print("RandomSeed=" + str(seed))

width=4

bvs = []
nbvs= 15
nsubsets=6
vals = []
for i in range(nbvs):
    r = random.randint(0,(1<<width)-1)
    bvs.append(BitVector(width, r))
    vals.append(r)
print(str(vals))


max_subsets=[]
min_subsets=[]
selecteds=[]
for i in range(nsubsets):
    max_subset=[]
    min_subset=[]
    selected=[]
    selecteds.append(selected)
    max_subsets.append(max_subset)
    min_subsets.append(min_subset)
    for bv in bvs:
        select=Var()
        selected.append(select)
        max_subset.append(Ite(select,bv,BitVector(width,0))) #1<<width-1
        min_subset.append(Ite(select,bv,BitVector(width,(1<<width)-1))) #1<<width-1
    AssertEqualPB(selected, nbvs//2)
    
maxs = []
for subset in max_subsets:
    maxs.append(Max(subset))

mins = []
for subset in min_subsets:
    mins.append(Min(subset))

for a,b in itertools.combinations(maxs,2):
    Assert(a!=b)
    
for a,b in itertools.combinations(mins,2):
    Assert(a!=b)
    
result =Solve()

print("Result is " + str(result))

if result:
    print(" ".join(str(bv.value()) for bv in bvs))
    for subset in selecteds:
        print(" ".join(str(s.value()) for s in subset))
    print("maxs:")
    for m in maxs:
        print(str(m.value()))
    print("mins:")
    for m in mins:
        print(str(m.value()))
    print("done")
    sys.exit(10)
else:
    sys.exit(20)
