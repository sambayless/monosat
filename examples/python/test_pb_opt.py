from monosat import *

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

