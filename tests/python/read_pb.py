#!/usr/bin/env python3
from __future__ import division, print_function
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
import bz2
#note: for testing purposes, not recommended for practical pb solving!
Monosat().newSolver("-verb=1 -pb-verb=1",output_file="/tmp/monosat_pb.gnf")



if __name__ == "__main__":
    seed = random.randint(1,100000)
    filename=None

    if len(sys.argv)>1:
        filename=sys.argv[1]
    if len(sys.argv)>2:
        seed=int(sys.argv[2])
random.seed(seed)

if filename is None:
    print("Usage: python3 read_pb.py <filename.opb>")
    sys.exit(0)

used_vars = dict()

goal = []

print("begin encode");
for line in (bz2.open(filename,"rt") if filename.endswith("bz2") else open(filename)):
    print(line, file=sys.stderr)
    line = line.strip();
    if line.startswith("*"):
        continue
    if(line.endswith(";")):
        line = line[0:-1]
    #this is a goal line
    words = line.split()
    is_objective = False
    objective = None
    if words[0] == "min:" or words[0] == "max:":
        is_objective = True
        objective=words[0]
        if words[0]=="min:":
            objective="min"
        else:
            objective="max"
        words = words[1:]


    const = None
    comp = None
    rhs = None
    lits = []
    pbs = []
    weights = []
    for w in words:

        if w.startswith("<") or w.startswith(">") or w.startswith("=") or w.startswith("!"):
            #this is a comparison
            if is_objective:
                raise Exception("Comparison in goal line")
            comp = w
        else:
            #try to convert w to a number
            try:
                assert (const is None)
                const = int(w)
                #this is a constant
                if comp is not None:
                    rhs = const
                    const=None
                    break

            except:
                neg = False
                if w[0]=="~":
                    neg = true()
                    w = w[1:]
                assert(w[0]=="x")
                varnum = int(w[1:])
                assert(varnum>0)
                if varnum not in used_vars:
                    used_vars[varnum] = Var("x%d"%(varnum)).getLit()

                var = Var(used_vars[varnum])
                if const is None:
                    const = 1
                lit = var if not neg else Not(var)
                lits.append(lit)
                pbs.append(varnum if not neg else -varnum)
                weights.append(const)
                const = None

    assert(const is None)
    assert(len(weights) == len(lits))
    #print(lits)
    #print(weights)
    #print(comp)
    #print(rhs)

    if not is_objective:
        assert(comp is not None)
        assert(rhs is not None)
        if (comp =="<"):
            AssertLessThanPB(lits,rhs,weights);
        elif (comp ==">"):
            AssertGreaterThanPB(lits,rhs,weights);
        elif (comp =="<="):
            AssertLessEqPB(lits,rhs,weights);
        elif (comp ==">="):
            AssertGreaterEqPB(lits,rhs,weights);
        elif (comp.startswith("=")):
            AssertEqualPB(lits,rhs,weights);
        else:
            raise Exception("Unknown operator " + str(comp))
    else:
        goal = [objective,lits,pbs,weights]
        assert(comp is None)
        assert(rhs is None)
        assert(objective is not None)
        if objective == "min":
            minimize(lits,weights)
        elif objective=="max":
            maximize(lits,weights)
        else:
            raise Exception("Unknown objective " + str(objective))
#print("RandomSeed=" + str(seed))

result=Solve()


print("Result is " + str(result))

if result:
    witness = []
    model = dict()
    witstr = "v ";
    for pb,l in used_vars.items():
        lit = Var(l)

        if lit.value():
            witstr+="x%d "%(pb)
            witness.append(pb)
            model[pb]=True
            model[-pb]=False
        elif lit.value()==False:
            witstr+="-x%d "%(pb)
            witness.append(-pb)
            model[pb]=False
            model[-pb]=True
        else:
            print("No full model for x%d"%(pb), file=sys.stderr)
            assert(False)
    witstr.rstrip();
    goalval = 0
    for pb,w in zip(goal[2],goal[3]):
        if pb not in model or -pb not in model:
            print("No model for x%d"%(pb), file=sys.stderr)
        assert(pb in model and -pb in model)

        if model[pb]:
            goalval+=w
    print("c Optimal solution: %d"%(goalval))
    print("s OPTIMUM FOUND")
    print("o %d"%(goalval))
    print(witstr)
else:
    print("s UNSATISFIABLE")
