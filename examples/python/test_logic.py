from monosat import *

a = Var() 
b = Var() 
c = Or(a, Not(b)) 

Assert(c)

result = Solve()
print(result)