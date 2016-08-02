from monosat import *

vars=[]
for v in range(10):
    vars.append(Var())

AssertEqualPB(vars,2);
AssertLessEqPB(vars, 4)
AssertGreaterThanPB(vars, 1)

result=Solve()
for v in vars:
    print(v.value())
print("Result is " + str(result))
assert(result==True)

