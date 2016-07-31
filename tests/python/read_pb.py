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
Monosat().init("-pb-verb=1")



if __name__ == "__main__":
    seed = random.randint(1,100000)
    filename=None

    if len(sys.argv)>1:
        filename=sys.argv[1]
    if len(sys.argv)>2:
        seed=int(sys.argv[2])
random.seed(seed)

used_vars = dict()

print("begin encode");
for line in bz2.open(filename,"rt"):
    print(line)
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
                if varnum not in used_vars:
                    used_vars[varnum] = Var("x%d"%(varnum))

                var = used_vars[varnum]
                if const is None:
                    const = 1
                lit = var if not neg else Not(var)
                lits.append(lit)
                weights.append(const)
                const = None

    assert(const is None)
    assert(len(weights) == len(lits))
    print(lits)
    print(weights)
    print(comp)
    print(rhs)

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
        assert(comp is None)
        assert(rhs is None)
        assert(objective is not None)
        if objective == "min":
            minimize(lits,weights)
        elif objective=="max":
            maximize(lits,weights)
        else:
            raise Exception("Unknown objective " + str(objective))
print("RandomSeed=" + str(seed))

result=Solve()

print("Result is " + str(result))
assert(result==True)

