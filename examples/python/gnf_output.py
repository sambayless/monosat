# Demonstrate writing and loading GNF files

from monosat import *
Monosat().setOutputFile("tutorial.gnf")
a = Var()
b = Var()
c = Or(a, Not(b))

Assert(c)
Assert(Not(c))

#trivial contradition, so result should be false
result = Solve()
print(result)

Monosat().init("") #clear and reinitialize the solver
#solver now has no constraints, so it should be sat
print(Solve())

# read in the previously saved constraints (these can also be read in on the command line, eg ./monosat tutorial.gnf
Monosat().readGNF("tutorial.gnf")

result2 = Solve()
print(result2)
assert(result==result2)