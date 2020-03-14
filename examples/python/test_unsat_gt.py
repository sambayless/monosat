from monosat import *

bv1 = BitVector(4)
bv2 = BitVector(4)
assert(Solve())
Assert(bv1==15)
assert(Solve())
Assert(bv2>bv1)
assert(not Solve())



