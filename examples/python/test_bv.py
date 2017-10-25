from monosat import *

bv = BitVector(32)
assert(Solve())
Assert(bv<=0)
Assert(bv>8)

result=Solve()
print("Result is " + str(result))
assert(result==False)

