from monosat import *

bv1 = BitVector(4)
bv2 = BitVector(4)
bv3 = bv1.concat(bv2)

Assert(bv1<=5)
Assert(bv1>0)
Assert(bv2 >=10)

result =Solve()
assert(result)
for v in reversed(bv1.bits()):
    print("1" if v.value() else 0,end='')
print(" = " + str(bv1.value()) +"\n")

for v in reversed(bv2.bits()):
    print("1" if v.value() else 0,end='')
print(" = " + str(bv2.value()) +"\n")

for v in reversed(bv3.bits()):
    print("1" if v.value() else 0,end='')
print(" = " + str(bv3.value()) +"\n")