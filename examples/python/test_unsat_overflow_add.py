from monosat import *

bv1 = BitVector(4)
assert(Solve())
# 15 is the largest value representable in a 4 bit unsigned bitvector
Assert(bv1==15)
# In MonoSAT, bitvectors cannot overflow
bv2 = bv1 + 1
# so bv2 = 15 + 1 will be UNSAT (whereas in other SMT solvers this might be SAT)
assert(not Solve())





