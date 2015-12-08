from monosat import *

S1 = Monosat().newSolver()
S2 = Monosat().newSolver()

Monosat().setSolver(S1)

a = Var() 
b = Var() 
c = Or(a, Not(b)) 

Assert(c)

result = Solve()
print(result)
Monosat().setSolver(S2)
a2 = Var() 
b2 = Var() 
c2 = Or(a2, Not(b2)) 

Assert(c2)
Assert(~c2)

result = Solve()
print(result)

Monosat().setSolver(S1)
result = Solve()
print(result)

print(str(a))
print(str(a2))

found_error=False
try:
    Assert(a2)
except RuntimeError as e:
    #correct, there should be an error here
    found_error=True    
assert(found_error)

Monosat().setSolver(S2)
found_error=False
try:
    Assert(a)
except RuntimeError as e:
    #correct, there should be an error here
    found_error=True    
assert(found_error)

print("Done")
    