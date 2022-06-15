# The MIT License (MIT)
#
# Copyright (c) 2014, Sam Bayless
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
# associated documentation files (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge, publish, distribute,
# sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or
# substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
# NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
# OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import monosat.monosat_c
import sys
from collections.abc import Iterable
from monosat.logic import *
from monosat.manager import Manager
from monosat.monosat_c import Monosat

debug = False


# Collects a set of graphs to encode together into a formula
class BVManager(metaclass=Manager):

    def __init__(self):
        self.bvs = []
        self.aux_bvs = []
        self.comparisons = []
        self.consts = dict()
        self.ites = []
        self._monosat = monosat.monosat_c.Monosat()
        self.bitblast_addition = False
        self.bitblast_addition_shadow = False

    def _setConstant(self, val, width, bv):
        self.consts[(val, width)] = bv

    def _hasConstant(self, val, width):
        return (val, width) in self.consts

    def _getConstant(self, val, width):
        return self.consts[(val, width)]

    def clear(self):
        self.bvs = []
        self.aux_bvs = []
        self.comparisons = []
        self.ites = []
        self.consts = []

    # Each "string" must actually be a list of positive integers
    def Bv(self, width=8, const_value=None):

        return BitVector(self, width, const_value)

    def hasBitVector(self, id):
        return id >= 0 and id < len(self.bvs) and self.bvs[id] is not None

    def getBitVector(self, id):
        return self.bvs[id]

    def getBitVectors(self):
        return self.bvs

    def Ite(self, i, t, e):
        assert isinstance(i, Var)
        assert isinstance(t, BitVector)
        assert isinstance(e, BitVector)
        assert t.width() == e.width()
        result = self.Bv(t.width())
        self._monosat.bv_ite(i.getLit(), t.getID(), e.getID(), result.getID())
        return result

    def Max(self, *bvs):
        width = None
        for bv in bvs:
            if not isinstance(bv, BitVector):
                raise TypeError("Arguments of Max must be bitvectors")
            if width is None:
                width = bv.width()
            if bv.width() != width:
                raise ValueError("All arguments of Max must have same bitwidth")

        return BitVector(self, width, "max", bvs)

    def Min(self, *bvs):
        width = None
        for bv in bvs:
            if not isinstance(bv, BitVector):
                raise TypeError("Arguments of Min must be bitvectors")
            if width is None:
                width = bv.width()
            if bv.width() != width:
                raise ValueError("All arguments of Min must have same bitwidth")

        return BitVector(self, width, "min", bvs)

    def write(self, f):
        for bv in self.bvs:
            bv.write(f)
        for i, bv in enumerate(self.aux_bvs):
            bv.pid = len(self.bvs) + i
            bv.write(f)
        for (i, t, e, r) in self.ites:
            f.write(
                "bv_Ite %d %d %d %d\n" % (dimacs(i), t.getID(), e.getID(), r.getID())
            )


def _bv_Max(*bvs):
    return BVManager().Max(*bvs)


def _bv_Min(*bvs):
    return BVManager().Min(*bvs)


def _bv_Ite(i, t, e):
    return BVManager().Ite(i, t, e)


# def Bv(width, const_value=None):
#    return BVManager().Bv(width,const_value)


def _checkBVs(bvs):
    for bv in bvs:
        assert isinstance(bv, BitVector)
        if bv.mgr._solver != Monosat().getSolver():
            raise RuntimeError(
                "Bitvector %s does not belong to current solver, aborting (use setSolver() to correct this)"
                % (str(bv))
            )


class BitVector:

    def __init__(self, mgr, width=None, op=None, args=None):
        assigned_bits = None
        if isinstance(mgr, int):
            # Shift the arguments over 1
            args = op
            op = width
            width = mgr
            mgr = BVManager()
        elif isinstance(mgr, Iterable):
            assigned_bits = list(mgr)
            mgr = BVManager()
            # Build this bitvector from a list of elements
            assert op is None
            assert args is None
            if width is not None:
                assert width == len(assigned_bits)
            width = len(assigned_bits)

        assert mgr._solver == Monosat().getSolver()
        assert width is not None

        self.mgr = mgr
        self._bv = None
        self.symbol = None

        self._width = width
        self._constant = None

        self.pid = None

        if args is None and isinstance(op, (int, float)):
            val = int(op)
            originalval = val

            if val < 0:
                val += 1 << width
                # val=0
                # print("Warning: negative bitvectors not yet supported, setting to 0", file=sys.stderr)

            if val >= (1 << width):
                val = (1 << width) - 1
                print(
                    "Warning: value %d is too large to represent with a width-%d bitvector, setting to %d"
                    % (originalval, width, val),
                    file=sys.stderr,
                )

            self._constant = val

            # if(val>0):
            #    print("Warning: value %d is too large to represent with a width-%d bitvector, setting to %d"%(originalval,width, (1<<width)-1), file=sys.stderr)
            op = None
            args = None
            if not mgr._hasConstant(val, width):
                self.pid = mgr._monosat.newBitvector_const(width, val)
                mgr._setConstant(val, width, self)
            else:
                self.pid = mgr._getConstant(val, width).pid
            self._bv = [None] * width
            # fill bv with constants, for convenience elsewhere
            for i in range(width - 1, -1, -1):
                v = 1 << i
                if val >= v:
                    val -= v
                    self._bv[i] = true()
                else:
                    self._bv[i] = false()
        elif (args is None and op == "anon") or op == "~":
            # create an anomymous bitvector (has no literals)
            self.pid = mgr._monosat.newBitvector_anon(width)
        else:
            self._bv = []
            if assigned_bits is None:
                for _ in range(width):
                    self._bv.append(Var())
            else:
                # Can't do this, because Monosat doesn't support multiple theories on the same variable
                # self._bv = assigned_bits
                for i in range(width):
                    v = Var()
                    self._bv.append(v)
                    AssertEq(v, assigned_bits[i])

            # arr = (c_int*width)()
            # for i,v in enumerate(self._bv):
            #    arr[i]=c_int(v.getLit()//2)
            bits = [v.getLit() // 2 for v in self._bv]

            self.pid = mgr._monosat.newBitvector(bits)

        if op == "+":
            assert len(args) == 2
            _checkBVs((self, args[0], args[1]))
            if not mgr.bitblast_addition:
                mgr._monosat.bv_addition(args[0].getID(), args[1].getID(), self.getID())
            if mgr.bitblast_addition or mgr.bitblast_addition_shadow:
                carry = false()
                for i, (a, b, out) in enumerate(zip(args[0], args[1], self)):
                    r, carry2 = Add(a, b, carry)
                    AssertEq(out, r)
                    carry = carry2
                Assert(Not(carry))  # disallow overflow.
        elif op == "-":
            _checkBVs((self, args[0], args[1]))
            # mgr._monosat.bv_addition(self.getID(), args[1].getID(), args[0].getID())
            mgr._monosat.bv_subtraction(args[0].getID(), args[1].getID(), self.getID())
        elif op == "*":
            mgr._monosat.bv_multiply(args[0].getID(), args[1].getID(), self.getID())
        elif op == "/":
            mgr._monosat.bv_divide(args[0].getID(), args[1].getID(), self.getID())
        elif op == "~":
            _checkBVs((self, args[0]))
            mgr._monosat.bv_not(args[0].getID(), self.getID())
        elif op == "&":
            _checkBVs((self, args[0]))
            mgr._monosat.bv_and(args[0].getID(), args[1].getID(), self.getID())
        elif op == "~&":
            _checkBVs((self, args[0]))
            mgr._monosat.bv_nand(args[0].getID(), args[1].getID(), self.getID())
        elif op == "|":
            _checkBVs((self, args[0]))
            mgr._monosat.bv_or(args[0].getID(), args[1].getID(), self.getID())
        elif op == "~|":
            _checkBVs((self, args[0]))
            mgr._monosat.bv_nor(args[0].getID(), args[1].getID(), self.getID())
        elif op == "^":
            _checkBVs((self, args[0]))
            mgr._monosat.bv_xor(args[0].getID(), args[1].getID(), self.getID())
        elif op == "~^":
            _checkBVs((self, args[0]))
            mgr._monosat.bv_xnor(args[0].getID(), args[1].getID(), self.getID())
        elif op == "slice":
            _checkBVs((self, args[0]))
            mgr._monosat.bv_slice(args[0].getID(), args[1], args[2], self.getID())
        elif op == "concat":
            _checkBVs((self, args[0], args[1]))
            mgr._monosat.bv_concat(args[0].getID(), args[1].getID(), self.getID())
        elif op == "min":
            _checkBVs((self,))
            _checkBVs((args))
            mgr._monosat.bv_min([x.getID() for x in args], self.getID())
        elif op == "max":
            _checkBVs((self,))
            _checkBVs((args))
            mgr._monosat.bv_max([x.getID() for x in args], self.getID())
        elif op == "popcount":
            mgr._monosat.bv_popcount((l.getLit() for l in args), self.getID())

    def isConst(self):
        return self._constant is not None

    def __repr__(self):
        if self.isConst():
            return "bv%d=" % (self.pid) + str(self._constant)
        else:
            return "bv%d" % (self.pid)

    def checkValue(self, val):
        if val < 0 or val >= 1 << self.width():
            print(
                "Error: value %d is too large to represent with a width-%d bitvector"
                % (val, self.width()),
                file=sys.stderr,
            )
            assert False

    def value(self):
        return self.mgr._monosat.getModel_BV(self.pid)

    def bitblast(self):
        self.mgr._monosat.bv_bitblast(self.pid)

    def getSymbol(self):
        return self.symbol

    def setSymbol(self, s):
        self.symbol = str(s)

    """                
    def __eq__(self,other):
        if other and isinstance(other,self.__class__):
            return other.getID()==self.getID()
        return False"""

    def __hash__(self):
        return self.getID()

    def width(self):
        return self._width

    def getID(self):
        return self.pid

    def __getitem__(self, index):
        if isinstance(index, slice):
            assert index.step is None or index.step == 1
            lower = index.start if index.start else 0
            upper = index.stop if index.stop else self.width() - 1
            assert upper > lower
            assert lower >= 0
            assert upper <= self.width()
            return BitVector(self.mgr, upper - lower, "slice", (self, lower, upper - 1))
        return self._bv[index]

    def bits(self):
        return list(self._bv)

    def __len__(self):
        return self._width

    def __add__(self, other):
        if not isinstance(other, BitVector):
            other = BitVector(self.mgr, self.width(), other)
        return BitVector(self.mgr, self.width(), "+", (self, other))

    __radd__ = __add__

    def __sub__(self, other):
        if not isinstance(other, BitVector):
            other = BitVector(self.mgr, self.width(), other)
        return BitVector(self.mgr, self.width(), "-", (self, other))

    __rsub__ = __sub__

    def __mul__(self, other):
        if not isinstance(other, BitVector):
            other = BitVector(self.mgr, self.width(), other)
        return BitVector(self.mgr, self.width(), "*", (self, other))

    __rmul__ = __mul__

    def __div__(self, other):
        if not isinstance(other, BitVector):
            other = BitVector(self.mgr, self.width(), other)
        return BitVector(self.mgr, self.width(), "/", (self, other))

    __rdiv__ = __div__

    def lt(self, compareTo):
        if isinstance(compareTo, BitVector):
            return Var(
                self.mgr._monosat.newBVComparison_bv_lt(self.getID(), compareTo.getID())
            )
        else:
            self.checkValue(int(compareTo))
            return Var(
                self.mgr._monosat.newBVComparison_const_lt(self.getID(), int(compareTo))
            )

    def leq(self, compareTo):
        if isinstance(compareTo, BitVector):
            return Var(
                self.mgr._monosat.newBVComparison_bv_leq(
                    self.getID(), compareTo.getID()
                )
            )
        else:
            self.checkValue(int(compareTo))
            return Var(
                self.mgr._monosat.newBVComparison_const_leq(
                    self.getID(), int(compareTo)
                )
            )

    def gt(self, compareTo):
        if isinstance(compareTo, BitVector):
            return Var(
                self.mgr._monosat.newBVComparison_bv_gt(self.getID(), compareTo.getID())
            )
        else:
            self.checkValue(int(compareTo))
            return Var(
                self.mgr._monosat.newBVComparison_const_gt(self.getID(), int(compareTo))
            )
        #    compareTo = BitVector(self.mgr,self.width(),compareTo)

    def geq(self, compareTo):
        if isinstance(compareTo, BitVector):
            return Var(
                self.mgr._monosat.newBVComparison_bv_geq(
                    self.getID(), compareTo.getID()
                )
            )
        else:
            self.checkValue(int(compareTo))
            return Var(
                self.mgr._monosat.newBVComparison_const_geq(
                    self.getID(), int(compareTo)
                )
            )

    def eq(self, compareTo):
        # if not isinstance(compareTo, BitVector):
        #    compareTo = BitVector(self.mgr,self.width(),compareTo)
        return And(self.leq(compareTo), self.geq(compareTo))

    def neq(self, compareTo):
        # if not isinstance(compareTo, BitVector):
        #    compareTo = BitVector(self.mgr,self.width(),compareTo)
        return Nand(self.leq(compareTo), self.geq(compareTo))

    def __lt__(self, other):
        return self.lt(other)

    def __le__(self, other):
        return self.leq(other)

    def __gt__(self, other):
        return self.gt(other)

    def __ge__(self, other):
        return self.geq(other)

    def __eq__(self, other):
        return self.eq(other)

    def __ne__(self, other):
        return self.neq(other)

    def Not(self):
        return BitVector(self.mgr, self.width(), "~", (self,))

    def __invert__(self):
        return self.Not()

    # bitwise operators
    def And(self, other):
        return BitVector(self.mgr, self.width(), "&", (self, other))

    def Nand(self, other):
        return BitVector(self.mgr, self.width(), "~&", (self, other))

    def Or(self, other):
        return BitVector(self.mgr, self.width(), "|", (self, other))

    def Nor(self, other):
        return BitVector(self.mgr, self.width(), "~|", (self, other))

    def Xor(self, other):
        return BitVector(self.mgr, self.width(), "^", (self, other))

    def Xnor(self, other):
        return BitVector(self.mgr, self.width(), "~^", (self, other))

    def __and__(self, other):
        return self.And(other)

    def __or__(self, other):
        return self.Or(other)

    def __xor__(self, other):
        return self.Xor(other)

    """def __getslice__(self,lower,upper):
        assert(upper>lower)
        assert(lower>0)
        assert(upper<=self.width())
        return BitVector(self.mgr,upper-lower,'slice',(self,lower,upper)) """

    def concat(self, other):
        return BitVector(
            self.mgr, self.width() + other.width(), "concat", (self, other)
        )

    def write(self, f):

        if self._constant is not None:
            f.write(
                "bv const %d %d " % (self.pid, self.width())
                + str(self._constant)
                + "\n"
            )

        else:
            f.write("bv %d %d" % (self.pid, self.width()))
            for var in self._bv:
                f.write(" %d" % (var.getInputLiteral()))
            f.write("\n")
        if self.symbol is not None:
            f.write("c bv %d " % (self.pid) + str(self.symbol) + "\n")

        if self.op == "+":
            f.write("bv + %d" % (self.pid))
            assert len(self.args) == 2
            for arg in self.args:
                # if isinstance(arg,BitVector):
                f.write(" %d" % (arg.getID()))
                # else:
                #    f. write(" " +str(arg))
            f.write("\n")

        for (v, compareTo, op) in self._comparisons:
            # if (isinstance(compareTo,BitVector)):
            f.write(
                "bv "
                + str(op)
                + " %d %d  " % (v.getInputLiteral(), self.pid)
                + " "
                + str(compareTo.getID())
                + "\n"
            )
            # else:
            #    f.write("bv " + str(op) + " %d bv%d  "%(v.getInputLiteral(), self.pid) + str(compareTo) + "\n")
