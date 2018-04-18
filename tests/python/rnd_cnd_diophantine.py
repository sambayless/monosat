from monosat import *

import functools
import math
import os
import random
import random
import sys

def run(seed,filename):
    #-debug-learnts=/tmp/testlearnts4.gnf -debug-bounds -debug-analysis
    Monosat().newSolver("  -verb=0 -verb-time=0 -rnd-theory-freq=0.99 -no-decide-bv-intrinsic  -decide-bv-bitwise  -decide-graph-bv -decide-theories -no-decide-graph-rnd   -lazy-maxflow-decisions -conflict-min-cut -conflict-min-cut-maxflow -reach-underapprox-cnf -check-solution ",output_file=filename)
    if filename is not None:
        print("Writing to %s"%(filename))
    
    
    random.seed(seed)
    print("RandomSeed=" + str(seed))
    
    bitwidth=8
    n_equations=3
    n_vars=4
    max_co = 10
    
    vars=[]
    for i in range(n_vars):
        vars.append(BitVector(bitwidth))
    
    coef_choice = []
    for i in range(n_vars):
        coef_choice.append(Var())
    
    for i in range(n_equations):
        coefficients=[]
        for c in range(n_vars):
            coefficients.append(random.randint(0,3))
        alt_coefficients=[]
        for c in range(n_vars):
            alt_coefficients.append(random.randint(0,3))
        val = random.randint(0,max_co*n_vars)
        rhs = BitVector(bitwidth, val)
        sum = BitVector(bitwidth,0)
        
        for j in range(len(coefficients)):
            opt_sum1 = BitVector(bitwidth,0) 
            opt_sum2 = BitVector(bitwidth,0) 
            for k in range( coefficients[j]):          
                opt_sum1= opt_sum1 + ( vars[j]) #simulation multiplication for now
            for k in range( alt_coefficients[j]):          
                opt_sum2= opt_sum2 + ( vars[j]) #simulation multiplication for now
            sum+=Ite(coef_choice[j],opt_sum1,opt_sum2)
                
        Assert(sum==rhs)
        print(str(i) + ": " + str(coefficients) + " =  " + str(val) )
         
    result =Solve()
    
    print("Result is " + str(result))
    
    if(result):
    
        for i in range(len(vars)):
            print(str(i) + "=  " + str(vars[i].value()))
    return result
    
if __name__ == "__main__":
    seed = random.randint(1,100000)
    filename=None
    
    if len(sys.argv)>1:
        filename=sys.argv[1]
    if len(sys.argv)>2:
        seed=int(sys.argv[2])
        
    r =run(seed,filename)
    if r:
        sys.exit(10)
    else:
        sys.exit(20)