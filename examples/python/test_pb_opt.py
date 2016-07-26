from monosat import *
Monosat().setOutputFile("/tmp/test.gnf")
vars=[]
for v in range(10):
    vars.append(Var())

AssertEqualPB(vars,2);
AssertLessEqPB(vars, 4)
AssertGreaterThanPB(vars, 1)

weights = list(range(10))

maximize(vars,weights)
result=Solve()
for v in vars:
    print(v.value())
print("Result is " + str(result))
assert(result==True)

