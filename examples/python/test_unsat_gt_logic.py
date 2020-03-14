from monosat import *

bv1 = BitVector(4)
bv2 = BitVector(4)
assert(Solve())
Assert(bv1==15)
assert(Solve())
# The Logic.py Arithmetic operators do not apply to bitvectors
# this method will throw an exception
try:
    GreaterThan(bv2, bv1)
    assert(False,"Should throw exception because Logic.GreaterThan does not apply to bitvectors")
except Exception:
    # Exception expected
    print("Logic.GreaterThan does not apply to bitvectors")




