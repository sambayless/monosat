from monosat import *
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

Monosat().newSolver(filename)

print("begin encode");

random.seed(seed)
print("RandomSeed=" + str(seed))
#
width=4


nbvs= 15
nsubsets=8

selects = []
for n in range(nbvs):
    selects.append(true() if random.random()>0.5 else false())

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
    weights = []
    for v in selects:
        weights.append(random.randint(-10,10))
        select = Var()
        selected.append(select)
        max_subset.append(If(select,v,false())) #1<<width-1
        #min_subset.append(If(select,Not(v),false())) #1<<width-1

    weightval = sum(weights[0:len(weights)//2])
    random.shuffle(weights)
    AssertEqualPB(selected, weightval,weights)



maxs = []
for subset in max_subsets:
    maxs.append(PopCount(subset,method="BV",return_bv=True))

#mins = []
#for subset in min_subsets:
#    mins.append(PopCount(subset,method="BV",return_bv=True))

for a,b in itertools.combinations(maxs,2):
    Assert(a!=b)
"""
for a,b in itertools.combinations(mins,2):
    Assert(a!=b)"""

result =Solve()

print("Result is " + str(result))

if result:
    print(" ".join(str(s.value()) for s in selects))
    for subset in selecteds:
        print(" ".join(str(s.value()) for s in subset))
    print("maxs:")
    for m in maxs:
        print(str(m.value()))
    #print("mins:")
    #for m in mins:
    #    print(str(m.value()))
    print("done")
    sys.exit(10)
else:
    sys.exit(20)






Monosat().setOutputFile("/tmp/test.gnf")
vars=[]
for v in range(10):
    vars.append(Var())

AssertEqualPB(vars,2);
AssertLessEqPB(vars, 4)
AssertGreaterThanPB(vars, 1)

weights = list(range(10))

maximize(vars,weights)
result=Solve()
sum = 0
for i,v in enumerate(vars):
    print(v.value())
    if (v.value()):
        sum+=weights[i]
print("Result is " + str(result))
assert(result==True)

