# Demonstrate writing and loading GNF files

from monosat import *
#To write to an output file, create a new solver and pass it the name of a file to write to
#All constraints created in that new solver will be copied to that file
Monosat().newSolver(output_file="example.gnf")
a = Var()
b = Var()
c = Or(a, Not(b))

Assert(c)
Assert(Not(c))

#trivial contradition, so result should be UNSAT
result = Solve()
print(result)

Monosat().newSolver() #create a new solver, this time without any output file set.
#solver now has no constraints, so it should be SAT
print(Solve())

# read in the previously saved constraints (these can also be read in on the command line, eg ./monosat tutorial.gnf
Monosat().readGNF("example.gnf")

result2 = Solve()
print(result2)
assert(result==result2)