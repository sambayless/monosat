from monosat import *

import functools
import math
import os
import random
import random
import sys

print("Generate random 0-1 diophantine equation")

seed = random.randint(1,100000)

random.seed(seed)
print("RandomSeed=" + str(seed))

bitwidth=8
n_equations=3
n_vars=3
max_co = 10

vars=[]
for i in range(n_vars):
    vars.append(BitVector(bitwidth))

for i in range(n_equations):
    coefficients=[]
    for c in range(n_vars):
        coefficients.append(random.randint(0,1))
    val = random.randint(0,max_co*n_vars)
    rhs = BitVector(bitwidth, val)
    sum = BitVector(bitwidth,0)
    
    for j in range(len(coefficients)):
        if coefficients[j]>0:
            sum = sum + ( vars[j])
    Assert(sum==rhs)
    print(str(i) + ": " + str(coefficients) + " =  " + str(val) )
     
result =Solve()

print("Result is " + str(result))

if(result):

    for i in range(len(vars)):
        print(str(i) + "=  " + str(vars[i].value()))


