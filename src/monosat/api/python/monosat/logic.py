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

from collections.abc import Iterable
import math
import warnings

import monosat.monosat_c
import numbers

_monosat = monosat.monosat_c.Monosat()


def _checkLits(vars):
    for v in vars:
        assert isinstance(v, Var)
        if v._solver != _monosat.getSolver():
            raise RuntimeError(
                "Variable %s does not belong to current solver, aborting (use setSolver() to correct this)"
                % (str(v))
            )


def getSymbols():
    return _monosat.symbolmap


class Var:

    def __init__(self, symbol=None, allow_simplification=False):
        # Warning: if allow_simplification is set to true, then the variable may be eliminated by the solver at the next Solve();
        # After that point, the variable would no longer be safe to use in clauses or constraints in subsequent calls to Solve().
        # Use carefully!
        self._solver = _monosat.getSolver()
        if isinstance(symbol, bool):
            self.lit = _monosat.true() if symbol else _monosat.false()
        elif isinstance(symbol, int):
            self.lit = symbol
        elif isinstance(symbol, Var):
            self.lit = symbol.lit
        else:
            self.lit = _monosat.newLit(allow_simplification)
            if symbol is not None:
                _monosat.setSymbol(self.lit, str(symbol))

    def __repr__(self):
        if self.isConstTrue():
            return "T"
        elif self.isConstFalse():
            return "F"
        else:
            return str(self.getLit())

    def getSymbol(self):
        return self.symbol

    def setDecisionPriority(self, priority):
        _monosat.setDecisionPriority(self.getVar(), priority)

    def getDecisionPriority(self):
        return _monosat.getDecisionPriority(self.getVar())

    def setDecisionPolarity(self, b):
        _monosat.setDecisionPolarity(self.getVar(), b)

    def getDecisionPolarity(self):
        return _monosat.getDecisionPolarity(self.getVar())

    def setDecisionVar(self, decidable):
        _monosat.setDecisionVar(self.getVar(), decidable)

    def isDecisionVar(self):
        return _monosat.isDecisionVar(self.getVar())

    def setSymbol(self, s):
        pass
        # if self.symbol is not None and self.symbol in _monosat.symbolmap:
        #    del _monosat.symbolmap[self.symbol]
        # if s is not None:
        #    self.symbol = str(s)
        #    _monosat.symbolmap[self.symbol]=self

    def __hash__(self):
        return self.lit

    def __eq__(self, other):
        return isinstance(other, Var) and self.lit == other.lit

    def getLit(self):
        return self.lit

    def getVar(self):
        return self.lit // 2

    def value(self):
        # 0 = true, 1=false, 2=unassigned
        # return _monosat.getModel_Literal(self.lit) ==0 #Treat unassigned as false
        v = _monosat.getModel_Literal(self.lit)
        if v == 0:
            return True
        elif v == 1:
            return False
        else:
            return None

    def isConst(self):
        return self.isConstTrue() or self.isConstFalse()

    def isConstTrue(self):
        return self.getLit() == _monosat.true()

    def isConstFalse(self):
        return self.getLit() == _monosat.false()

    def sign(self):
        return (self.lit % 2) != 0

    def And(self, other):
        o = VAR(other)
        if self.isConstFalse() or o.isConstFalse():
            return false()

        if self.isConstTrue():
            return o
        if o.isConstTrue():
            return self

        v = Var()
        _checkLits((self, o, v))
        _monosat.addTertiaryClause(
            v.getLit(), _monosat.Not(self.getLit()), _monosat.Not(o.getLit())
        )
        _monosat.addBinaryClause(_monosat.Not(v.getLit()), self.getLit())
        _monosat.addBinaryClause(_monosat.Not(v.getLit()), o.getLit())
        return v
        # return Var(_monosat.addAnd( self.getLit(),o.getLit()))

    def Or(self, other):
        o = VAR(other)
        if self.isConstTrue() or o.isConstTrue():
            return true()

        if self.isConstFalse():
            return o
        if o.isConstFalse():
            return self
        v = Var()
        _checkLits((self, o, v))
        _monosat.addTertiaryClause(_monosat.Not(v.getLit()), self.getLit(), o.getLit())
        _monosat.addBinaryClause(v.getLit(), _monosat.Not(self.getLit()))
        _monosat.addBinaryClause(v.getLit(), _monosat.Not(o.getLit()))

        return v
        # return Var()
        # return Var(_monosat.Not( _monosat.addAnd(  _monosat.Not( self.getLit()),_monosat.Not(o.getLit()))))

    def Nor(self, other):
        return ~(self.Or(other))
        """o=VAR(other)
        if(self.isConstTrue() or o.isConstTrue()):
            return false()
        
        if(self.isConstFalse()):
            return o.Not();
        if(o.isConstFalse()):
            return self.Not();"""

        # return Var( _monosat.addAnd(  _monosat.Not( self.getLit()),_monosat.Not(o.getLit())))

    def Nand(self, other):
        return ~(self.And(other))
        """o=VAR(other)
        if(self.isConstFalse() or o.isConstFalse()):
            return true()
        
        if(self.isConstTrue()):
            return o.Not();
        if(o.isConstTrue()):
            return self.Not();
        return Var(_monosat.Not(_monosat.addAnd( self.getLit(),o.getLit())))"""

    # Keep this operator overload
    def __invert__(self):
        return self.Not()

    def Not(self):
        return Var(_monosat.Not(self.getLit()))

    # No longer overloading the __==__ operator, as that caused problems with dictionaries/hash functions.
    def Eq(self, other):
        return self.Xnor(other)

    def Neq(self, other):
        return self.Xor(other)

    # Note: No longer overloading the __xor__ operator, as that caused problems with dictionaries/hash functions.
    def Xor(self, other):
        o = VAR(other)

        if self.isConst() and o.isConst():
            return true() if self.getLit() != o.getLit() else false()

        if self.isConstTrue():
            return o.Not()

        if self.isConstFalse():
            return o

        if o.isConstTrue():
            return self.Not()

        if o.isConstFalse():
            return self

        v = Var()
        _checkLits((self, o, v))
        # If both inputs 0, output must be 0.
        _monosat.addTertiaryClause(_monosat.Not(v.getLit()), self.getLit(), o.getLit())

        # If one in put 1, the other input 0, output must be 1.
        _monosat.addTertiaryClause(v.getLit(), _monosat.Not(self.getLit()), o.getLit())

        # If one in put 1, the other input 0, output must be 1.
        _monosat.addTertiaryClause(v.getLit(), self.getLit(), _monosat.Not(o.getLit()))

        # If both inputs 1, output must be 0.
        _monosat.addTertiaryClause(
            _monosat.Not(v.getLit()),
            _monosat.Not(self.getLit()),
            _monosat.Not(o.getLit()),
        )

        # a=_monosat.addAnd(self.getLit(),_monosat.Not( o.getLit()))
        # b=_monosat.addAnd(_monosat.Not(self.getLit()),o.getLit())

        # return Var(_monosat.Not(_monosat.addAnd(_monosat.Not( a),_monosat.Not(b))))
        return v

    def Xnor(self, other):
        return ~(self.Xor(other))
        """ o=VAR(other)
        
        if(self.isConst() and o.isConst()):
            return true() if self.getLit() == o.getLit() else false()
        
        if(self.isConstFalse()):
            return o.Not();
        
        if(self.isConstTrue()):
            return o;
        
        if(o.isConstFalse()):
            return self.Not();
        
        if(o.isConstTrue()):
            return self;
        
        
        a=_monosat.addAnd(self.getLit(),_monosat.Not( o.getLit()))
        b=_monosat.addAnd(_monosat.Not(self.getLit()),o.getLit())
        
        return Var(_monosat.addAnd(_monosat.Not( a),_monosat.Not(b)))        """

    def Implies(self, other):
        o = VAR(other)
        _checkLits((self, o))
        # return ~(self.And(~other))
        return Var(
            _monosat.Not(_monosat.addAnd(self.getLit(), _monosat.Not(o.getLit())))
        )


def true():
    return Var(_monosat.true())


def false():
    return Var(_monosat.false())


# for internal use only
def _addClause(args):
    clause = [VAR(x) for x in args]
    _checkLits(clause)
    _monosat.addClause([x.getLit() for x in clause])


def _addSafeClause(args):
    _checkLits(args)
    _monosat.addClause([x.getLit() for x in args])


def IsBoolVar(a):
    from monosat.bvtheory import BitVector

    if a is None:
        return False
    if isinstance(a, Var):
        return True
    elif isinstance(a, bool):
        # Yes, I mean this to be an explicit check against true, false literals, and not against falsy-ness.
        return True
    elif isinstance(a, BitVector) or isinstance(a, numbers.Integral):
        return False
    else:
        return False


def VAR(a):
    from monosat.bvtheory import BitVector

    if a is None:
        raise TypeError("Cannot convert None to symbolic Booleans")
    if isinstance(a, Var):
        return a
    elif isinstance(a, bool):
        # Yes, I mean this to be an explicit check against true, false literals, and not against falsy-ness.
        if a:
            return true()
        else:
            return false()
    elif isinstance(a, BitVector) or isinstance(a, numbers.Integral):
        raise TypeError(
            "Cannot use BitVectors as arguments to functions requiring symbolic Booleans"
        )
    else:
        raise TypeError("Cannot convert to symbolic Booleans")


def Ite(i, t, e):
    # Is there a cleaner way to do this?
    from monosat.bvtheory import BitVector

    if IsBoolVar(t) and IsBoolVar(e):
        return _boolean_Ite(i, t, e)
    if isinstance(t, BitVector):
        bi = VAR(i)
        if isinstance(e, BitVector):
            pass
        elif isinstance(e, numbers.Integral):
            e = monosat.BitVector(t.width(), int(e))
        else:
            raise TypeError(
                "Cannot mix BitVectors and symbolic Booleans in If-Then-Else arguments"
            )
        return monosat.bvtheory._bv_Ite(bi, t, e)
    elif isinstance(t, numbers.Integral) and isinstance(e, BitVector):
        bi = VAR(i)
        t = monosat.BitVector(t.width(), int(e))
        return monosat.bvtheory._bv_Ite(bi, t, e)
    elif isinstance(t, numbers.Integral) and isinstance(e, numbers.Integral):
        # Cannot (easily) derive the width of the bitvector if both arguments are ints
        raise TypeError(
            "At least one of the (Then,Else) arguments to If-Then-Else must be a symbolic Boolean or a BitVector (they cannot both be Python primitives)"
        )
    else:
        raise TypeError(
            "Arguments to If-Then-Else must either be symbolic Booleans, or BitVectors"
        )


# This needs to work for other types of then/else args, too
def _boolean_Ite(i, t, e):
    ni = VAR(i)
    nt = VAR(t)

    if e is None:
        e = true()  # Not false!
    ne = VAR(e)

    if isTrue(ni):
        return t
    elif isFalse(ni):
        return e
    _checkLits((ni, nt, ne))
    l = _monosat.Not(_monosat.addAnd(ni.getLit(), _monosat.Not(nt.getLit())))
    r = _monosat.Not(
        _monosat.addAnd(_monosat.Not(ni.getLit()), _monosat.Not(ne.getLit()))
    )

    return Var(_monosat.addAnd(l, r))


def If(condition, thn, els=None):
    return Ite(condition, thn, els)


def And(*args):
    if len(args) == 0:
        return false()
    elif len(args) == 1:
        if isinstance(args[0], Iterable):
            return And(*args[0])
        return VAR(args[0])
    else:
        out = Var()
        clause = []
        for a in args:
            clause.append(_monosat.Not(VAR(a).getLit()))
            _monosat.addBinaryClause(a.getLit(), _monosat.Not(out.getLit()))
        clause.append(out.getLit())
        _monosat.addClause(clause)
        return out


def Or(*args):
    if len(args) == 0:
        return false()
    elif len(args) == 1:
        if isinstance(args[0], Iterable):
            return Or(*args[0])
        return VAR(args[0])
    else:
        out = Var()
        clause = []
        for a in args:
            clause.append(VAR(a).getLit())
            _monosat.addBinaryClause(_monosat.Not(a.getLit()), out.getLit())
        clause.append(_monosat.Not(out.getLit()))
        _monosat.addClause(clause)
        return out


def Nand(*args):
    return And(*args).Not()


def Nor(*args):
    return Or(*args).Not()


def Not(a):
    return VAR(a).Not()


def Xor(*args):
    if len(args) == 0:
        return false()
    elif len(args) == 1:
        if isinstance(args[0], Iterable):
            return Xor(*args[0])
        return VAR(args[0])
    else:
        a = VAR(args[0])
        for i in range(1, len(args)):
            a = a.Xor(VAR(args[i]))
        return a

    # return VAR(a).Xor(VAR(b))


def Xnor(*args):
    return Xor(*args).Not()


# def Xnor(a,b):
#    return VAR(a).Xnor(VAR(b))

# a implies b
def Implies(a, b):
    return Or(Not(VAR(a)), VAR(b))  # VAR(a) & Not(b)


def Eq(a, b):
    return Xnor(a, b)


def Neq(a, b):
    return Xor(a, b)


# Asserting versions of these, to avoid allocating extra literals...
def AssertTrue(*args):
    for arg in args:
        Assert(arg)


def AssertFalse(*args):
    for arg in args:
        Assert(~arg)


def AssertOr(*args):
    if len(args) == 0:
        Assert(false())
    elif len(args) == 1:
        if isinstance(args[0], Iterable):
            _addClause(args[0])
    else:
        _addClause(args)


def AssertNor(*args):
    if len(args) == 0:
        Assert(true())
    elif len(args) == 1:
        if isinstance(args[0], Iterable):
            AssertNor(*args[0])
    else:
        # AssertAnd((v.Not() for v in args))
        for v in args:
            Assert(VAR(v).Not())


def AssertAnd(*args):
    if len(args) == 0:
        Assert(false())
    elif len(args) == 1:
        if isinstance(args[0], Iterable):
            AssertAnd(*args[0])
    else:
        for v in args:
            Assert(v)


def AssertNand(*args):
    if len(args) == 0:
        Assert(true())
    elif len(args) == 1:
        if isinstance(args[0], Iterable):
            AssertNand(*args[0])
    else:
        AssertOr(*[VAR(v).Not() for v in args])


def AssertXor(*args):
    if len(args) == 0:
        return false()
    elif len(args) == 1:
        if isinstance(args[0], Iterable):
            AssertXor(*args[0])
    elif len(args) == 2:
        _addSafeClause((VAR(args[0]), VAR(args[1])))
        _addSafeClause((~VAR(args[0]), ~VAR(args[1])))
    else:
        # This can probably be done better...
        a = VAR(args[0])
        for i in range(1, len(args)):
            a = a.Xor(VAR(args[i]))
        Assert(a)


def AssertXnor(*args):
    if len(args) == 0:
        return true()
    elif len(args) == 1:
        if isinstance(args[0], Iterable):
            AssertXnor(*args[0])
    elif len(args) == 2:
        _addSafeClause((~VAR(args[0]), VAR(args[1])))
        _addSafeClause((VAR(args[0]), ~VAR(args[1])))
    else:
        # This can probably be done better...
        a = VAR(args[0])
        for i in range(1, len(args)):
            a = a.Xor(VAR(args[i]))
        Assert(~a)


# a implies b
def AssertImplies(a, b):
    _addSafeClause((~VAR(a), VAR(b)))
    # return Not(VAR(a)).Or(  VAR(b)) #VAR(a).And( Not(b))


def AssertEachImplies(pairs):
    _monosat.addBinaryClauses([(x.getLit(), y.getLit()) for x, y in pairs])


def AssertIff(a, b):
    AssertXnor((a, b))


def AssertEq(*args):
    if len(args) == 1 and isinstance(args[0], Iterable):
        return AssertEq(*args[0])
    if len(args) == 2:
        AssertXnor(args)
    elif len(args) > 2:
        AssertOr(And(args), Nor(args))
    else:
        print(
            "Warning - equality constraints have no affect on arguments of length %d"
            % (len(args))
        )
        pass


def AssertNeq(a, b):
    AssertXor((a, b))
    # return Xor(a,b)


def AssertClause(clause):
    _addClause(clause)
    # _monosat.addClause([a.getLit() for a in clause])


def Assert(a):
    if a is False or a is True:
        raise Exception(
            "Error: asserted python Boolean literal (not symbolic variable)"
        )

    a = VAR(a)
    if a.isConstFalse():
        warnings.warn(
            "Asserted a trivial contradiction",
            UserWarning,
        )
        # print("Warning: asserted constant false variable")
    _checkLits((a,))
    _monosat.addUnitClause(a.getLit())


def AssertIf(If, Thn, Els=None):
    if Els:
        AssertAnd(Or(Not(If), Thn), Or(If, Els))
    else:
        AssertOr(Not(If), Thn)


def HalfAdd(a, b):
    aV = VAR(a)
    bV = VAR(b)

    sum = a.Xor(b)
    carry = a.And(b)

    return sum, carry


# Add these bits (or arrays of bits)
def Add(a, b, c=False):
    if isinstance(a, Iterable):
        return _AddArray(a, b, c)
    aV = VAR(a)
    bV = VAR(b)
    cV = VAR(c)

    if cV.isConstFalse():
        sum = a.Xor(b)
        carry = a.And(b)
    elif cV.isConstTrue():
        sum = (a.Xor(b)).Not()
        carry = a.Or(
            b
        )  # Because carry is true, if either a or b, then we will output a carry
    else:
        sum = (a.Xor(b)).Xor(c)
        carry = (a.And(b)).Or((c.And((a.Xor(b)))))

    return sum, carry


# Subtract these bits (or arrays of bits) using 2's complement
def Subtract(a, b):
    if isinstance(a, Iterable):
        return _SubtractArray(a, b)
    aV = VAR(a)
    bV = Not(VAR(b))

    sum = (a.Xor(b)).Not()
    carry = a.Or(
        b
    )  # Because carry is true, if either a or b, then we will output a carry

    return sum, carry


def AddOne(array, bit=True):
    if isFalse(bit):
        return list(array)
    a1 = list(array)
    output = []
    sum, carry = HalfAdd(a1[0], bit)
    output.append(sum)
    for i in range(1, len(a1)):
        sum, carry = HalfAdd(a1[i], carry)
        output.append(sum)
    return output, carry


def _numberAsArray(number):
    if isinstance(number, (int, float, complex)):
        val1 = number
        n = 0
        num = []
        while 1 << n <= val1:
            if (1 << n & val1) == 0:
                num.append(false())
            else:
                num.append(true())
            n += 1
        return num
    else:
        raise (Exception("Cannot convert to binary: " + str(number)))


def _AddArray(array1, array2, carry=False):
    if not isinstance(array1, Iterable):
        array1 = _numberAsArray(array1)
    if not isinstance(array2, Iterable):
        array2 = _numberAsArray(array2)

    a1 = list(array1)
    a2 = list(array2)
    while len(a1) < len(a2):
        a1.append(false())
    while len(a2) < len(a1):
        a2.append(false())

    output = []
    for i in range(0, len(a1)):
        sum, carry = Add(a1[i], a2[i], carry)
        output.append(sum)

    output.append(carry)
    return output, carry
    # return output+carry


def _SubtractArray(array1, array2):
    """returns array1 - array2
    """
    a1 = list(array1)
    a2 = list(array2)
    while len(a1) < len(a2):
        a1.append(false())
    while len(a2) < len(a1):
        a2.append(false())
    carry = true()
    output = []
    for i in range(0, len(a1)):
        sum, carry = Add(a1[i], Not(a2[i]))
        output.append(sum)

    return output, carry


def Min(*args):
    from monosat.bvtheory import BitVector

    if len(args) == 1 and isinstance(args[0], Iterable):
        return Min(*args[0])
    if len(args) > 0 and isinstance(args[0], BitVector):
        return monosat.bvtheory._bv_Min(*args)
    else:
        assert len(args) == 2  # In the future, support more arguments...
        return _Min(*args)


def Max(*args):
    from monosat.bvtheory import BitVector

    if len(args) == 1 and isinstance(args[0], Iterable):
        return Max(*args[0])

    if len(args) > 0 and isinstance(args[0], BitVector):
        return monosat.bvtheory._bv_Max(*args)
    else:
        assert len(args) == 2  # In the future, support more arguments...
        return _Max(*args)


def _Min(array1, array2):
    a = list(array1)
    b = list(array2)
    while len(a) < len(b):
        a.append(false())
    while len(b) < len(a):
        b.append(false())

    a_less_eq = LessEq(a, b)
    output = []
    for (a_bit, b_bit) in zip(a, b):
        output.append(If(a_less_eq, a_bit, b_bit))
    return output


def _Max(array1, array2):
    a = list(array1)
    b = list(array2)
    while len(a) < len(b):
        a.append(false())
    while len(b) < len(a):
        b.append(false())

    a_less_eq = LessEq(a, b)
    output = []
    for (a_bit, b_bit) in zip(a, b):
        output.append(If(a_less_eq, b_bit, a_bit))
    return output


# True IFF num1 is < num2, <= num2
# Note: This method operates on lists of Boolean literals. It does not support BitVector arguments.
# Use BitVector.geq/leq/gt/lt to compare BitVectors.
def _LessOrEqual(num1, num2):
    from monosat.bvtheory import BitVector
    if isinstance(num1, BitVector) or isinstance(num2, BitVector):
        warnings.warn(
            "Logic.LessEq/GreaterEq/LessThan/GreaterThan "
            "do not apply to bitvectors, use BitVector.geq/leq/gt/lt instead",
            UserWarning,
        )
        raise Exception("Logic._LessOrEqual does not apply to bitvector arguments (use BitVector.geq/leq/gt/lt instead)")

    if isinstance(num1, (int, float, complex)) and isinstance(
            num2, (int, float, complex)
    ):
        return (num1 < num2, num1 <= num2)

    if isinstance(num1, (int, float, complex)):
        num1 = _numberAsArray(num1)
    if isinstance(num2, (int, float, complex)):
        num2 = _numberAsArray(num2)

    allconst2 = True
    if not isinstance(num2, (int, float, complex)):
        for v in num2:
            aV = VAR(v)
            if not aV.isConst():
                allconst2 = False
                break

    if allconst2:
        a1 = list(num1)
        if isinstance(num2, (int, float, complex)):
            val2 = num2
            n = 0
            num2 = []
            while 1 << n <= val2:
                if (1 << n & val2) == 0:
                    num2.append(false())
                else:
                    num2.append(true())
                n += 1
        else:
            val2 = 0
            for i in range(len(num2)):
                aV = VAR(num2[i])
                assert aV.isConst()
                if aV.isConstTrue():
                    val2 += 1 << i

        if val2 > math.pow(2, len(num1)):
            return (true(), true())

        a2 = list(num2)
        while len(a1) < len(a2):
            a1.append(false())
        while len(a2) < len(a1):
            a2.append(false())

        lt = false()
        gt = false()
        for i in range(len(a1) - 1, -1, -1):
            a = VAR(a1[i])
            b = VAR(a2[i])
            gt = gt.Or(And(Not(lt), a, Not(b)))
            lt = gt.Or(And(Not(gt), Not(a), b))

        return lt, Or(lt, Not(gt))

    a1 = list(num1)
    a2 = list(num2)
    while len(a1) < len(a2):
        a1.append(false())
    while len(a2) < len(a1):
        a2.append(false())
    lt = false()
    gt = false()
    for i in range(len(a1) - 1, -1, -1):
        a = VAR(a1[i])
        b = VAR(a2[i])
        gt = gt.Or(And(Not(lt), a, Not(b)))
        lt = gt.Or(And(Not(gt), Not(a), b))

    return lt, Or(lt, Not(gt))


def isTrue(v):
    return VAR(v).isConstTrue()


def isFalse(v):
    return VAR(v).isConstFalse()


def isConst(v):
    return VAR(v).isConst()


def LessThan(num1, num2):
    return _LessOrEqual(num1, num2)[0]


def LessEq(num1, num2):
    return _LessOrEqual(num1, num2)[1]


def Eq(*args):
    return Equal(*args)


def Equal(*args):
    if len(args) == 1 and isinstance(args, Iterable):
        args = args[0]

    if len(args) < 2:
        return true()
    elif len(args) == 2:
        return Equal2_(args[0], args[1])
    from monosat.bvtheory import BitVector

    all_vars = True
    all_bv = True
    assert len(args) >= 2
    for a in args:
        if not isinstance(a, Var):
            all_vars = False
        if not isinstance(a, BitVector):
            all_bv = False
    assert not (all_bv and all_vars)

    if all_vars:
        allequal = true()
        for i in range(len(args) - 1):
            ai = VAR(args[i])
            bi = VAR(args[i + 1])
            allequal = allequal.And(ai.Xnor(bi))
        return allequal
    elif all_bv:
        allequal = true()
        for i in range(len(args) - 1):
            ai = args[i]
            bi = args[i + 1]
            assert isinstance(ai, BitVector)
            assert isinstance(bi, BitVector)
            allequal = allequal.And(ai == bi)
        return allequal
    else:
        assert False


def Equal2_(num1, num2):
    if isinstance(num1, (bool, int, float, complex)) and isinstance(
            num2, (bool, int, float, complex)
    ):
        return num1 == num2
    elif isinstance(num1, Var) and isinstance(num2, Var):
        return num1.Xnor(num2)
    elif isinstance(num1, (bool, int, float, complex)):
        return Equal(num2, num1)
    elif isinstance(num2, (bool, int, float, complex)):
        n = 0
        t = []
        while 1 << n <= num2:
            if (1 << n & num2) == 0:
                t.append(false())
            else:
                t.append(true())
            n += 1

        num2 = t

    a = list(num1)
    b = list(num2)
    while len(a) < len(b):
        a.append(false())
    while len(b) < len(a):
        b.append(false())

    allequal = true()
    for i in range(len(a)):
        ai = VAR(a[i])
        bi = VAR(b[i])
        allequal = allequal.And(ai.Xnor(bi))
    return allequal


def GreaterThan(num1, num2):
    return LessThan(num2, num1)


def GreaterEq(num1, num2):
    return LessEq(num2, num1)


def PopCount(vars, **kwargs):
    method = "TREE" if "method" not in kwargs else kwargs["method"]
    if method == "CSA":
        return _PopCountHarleySeal(vars)
    elif method == "TREE":
        return _PopCountPairs(vars)
    elif method == "TABLE":
        return _PopCountTable(vars, kwargs["max"] if "max" in kwargs else None)
    elif method == "UNARY":
        # Warning: this will return the population count in a unary, not binary, vector, of length len(vars).
        return _PopCountUnary(vars)
    elif method == "BVSUM":
        # Warning: this will return the population count in a unary, not binary, vector, of length len(vars).
        return _PopCountBVSum(
            vars,
            kwargs["tree_addition"] if "tree_addition" in kwargs else False,
            kwargs["return_bv"] if "return_bv" in kwargs else False,
        )
    elif method == "BV":
        # Warning: this will return the population count in a unary, not binary, vector, of length len(vars).
        return _PopCountBV(
            vars, kwargs["return_bv"] if "return_bv" in kwargs else False
        )

    elif method == "NAIVE":
        return _PopCountNaive(vars)
    else:
        raise Exception("Unknown popcount method " + str(method))


def _next_power_of_2(x):
    # from http://stackoverflow.com/a/14267557
    return 2 ** (x - 1).bit_length()


def _grouper(n, iterable, padvalue=None):
    from itertools import chain, repeat

    "grouper(3, 'abcdefg', 'x') --> ('a','b','c'), ('d','e','f'), ('g','x','x')"
    return zip(*[chain(iterable, repeat(padvalue, n - 1))] * n)


def _PopCountUnary(vars):
    """
    vars = [false() for v in vars]
    for v in vars:
        
    
    
    #Add some redundant constraints, so that the sovler doesn't need to work to learn them:
    any_vars=Or(vars)
    any_outs = Or(output)
    AssertEq(any_vars,any_outs)
    
    return output;  """
    pass


def _PopCountBV(vars, retBitVector=False):
    from monosat.bvtheory import BitVector

    bvwidth = math.ceil(math.log(len(vars), 2))
    sm = BitVector(bvwidth, "popcount", vars)
    # Add some redundant constraints, so that the sovler doesn't need to work to learn them:
    output = sm.bits()
    any_vars = Or(vars)
    any_outs = Or(output)
    AssertEq(any_vars, any_outs)

    if retBitVector:
        return sm
    else:
        return output


def _PopCountBVSum(vars, tree_addition=False, retBitVector=False):
    from monosat.bvtheory import BitVector

    bvwidth = math.ceil(math.log(len(vars), 2))
    if not tree_addition:
        sm = bv(bvwidth, 0)
        for v in vars[i]:
            sm = If(v, sm + 1, sm)

    else:
        elements = [If(v, BitVector(bvwidth, 1), BitVector(bvwidth, 0)) for v in vars]
        while len(elements) > 1:
            if len(elements) % 2 == 1:
                elements.append(bv(bvwidth, 0))

            next_elements = []
            for a, b in zip(elements[::2], elements[1::2]):
                next_elements.append(a + b)
            elements = next_elements
        assert len(elements) == 1
        sm = elements[0]

    # Add some redundant constraints, so that the sovler doesn't need to work to learn them:
    output = sm.bits()
    any_vars = Or(vars)
    any_outs = Or(output)
    AssertEq(any_vars, any_outs)

    if retBitVector:
        return sm
    else:
        return output


def _PopCountTable(vars, max_value=None):
    if max_value is None:
        max_value = len(vars)

    if max_value > 10:
        print(
            "Warning: building a pop-count table for more than 10 values (was %d) may require a large amount of memory!"
            % (max_value),
            file=sys.stderr,
        )

    output = []
    maxwidth = math.ceil(math.log(len(vars), 2))
    for v in maxwidth:
        output.append(Var())

    for n in range((1 << len(vars)) - 1):
        pop = bin(n).count("1")
        if pop <= max_value:  # make this faster in the future!!
            a = true()
            for i, b in enumerate(bin(n)):
                v = vars[i]
                if b == "0":
                    a = a.And(Not(v))
                else:
                    a = a.And(v)
            for i in range(maxwidth):
                bit = 1 << i
                if pop & bit:
                    AssertImplies(a, output[i])
                else:
                    AssertImplies(a, Not(output[i]))

    # Add some redundant constraints, so that the sovler doesn't need to work to learn them:
    any_vars = Or(vars)
    any_outs = Or(output)
    AssertEq(any_vars, any_outs)

    return output


def _PopCountPairs(vars):
    # From Hackers Delight, 2nd ed.
    # required_len = _next_power_of_2(len(vars))
    # while( len(vars)<required_len):
    #    vars = vars+false()
    maxwidth = math.ceil(math.log(len(vars), 2)) + 1
    total = [false()] * maxwidth  # is it safe to make this maxwidth-1?
    # ones=[false()]*len(vars)
    ones = false()
    size = 2
    # Split the vars into adjacent pairs, sum each pair, the repeat.
    sums = []

    for a, b in _grouper(2, vars, None):
        sums.append((a, b))

    while len(sums) > 1:
        next_sums = []
        for a, b in _grouper(2, sums, None):
            assert a is not None
            if b is not None:
                s, carry = _AddArray(a, b)
                # carry will always be false
            else:
                s = a
            next_sums.append(s)
        sums = next_sums

        # assert(ignore.isConstFalse())
    assert len(sums) == 1
    output = sums[0]
    # Add some redundant constraints, so that the sovler doesn't need to work to learn them:
    any_vars = Or(vars)
    any_outs = Or(output)
    AssertEq(any_vars, any_outs)

    # Can't do this, unless vars is a power of 2...
    # all_vars=And(vars)
    # all_outs=And(total2)
    # AssertEq(all_vars,all_outs)
    return output


def _CSA(a, b, c):
    assert isinstance(a, Iterable) == isinstance(b, Iterable)
    assert isinstance(a, Iterable) == isinstance(c, Iterable)
    if isinstance(a, Iterable):
        assert len(a) == len(b)
        assert len(a) == len(c)
        ps = []
        sc = []
        for i in range(len(a)):
            psi, sci = _SCA(a[i], b[i], c[i])
            ps.append(psi)
            sc.append(sci)
        return ps, sc

    # from wikipedia
    # ps = a xor b xor c
    # sc = and(a,  b) or and(a,c) or and(b,c)
    ps = Xor(a, b, c)
    sc = Or(And(a, b), And(a, c), And(b, c))
    return ps, sc


def _PopCountHarleySeal(vars):
    # From Hackers Delight, 2nd ed.
    if len(vars) % 2 == 1:
        vars = vars + false()
    maxwidth = math.ceil(math.log(len(vars), 2)) + 1
    total = [false()] * maxwidth  # is it safe to make this maxwidth-1?
    # ones=[false()]*len(vars)
    ones = false()
    for i in range(0, len(vars), 2):
        twos, ones = _CSA(ones, vars[i], vars[i + 1])
        total, ignore = AddOne(total, twos)  # _AddArray(total,)
        # assert(ignore.isConstFalse())
    # total = 2*total + pop(ones) - operating at the bit level, this is equivalent to left shifting total, and putting 'ones' in the 1's column
    total2 = [ones]
    total2.extend(total)

    # Add some redundant constraints, so that the sovler doesn't need to work to learn them:
    any_vars = Or(vars)
    any_outs = Or(total2)
    AssertEq(any_vars, any_outs)

    # Can't do this, unless vars is a power of 2...
    # all_vars=And(vars)
    # all_outs=And(total2)
    # AssertEq(all_vars,all_outs)
    return total2


"""def _PopCountHarleySeal4_3(vars):

    while len(vars)%4!=0:
        vars = vars+false()
    
    sums=[] 
    for i in range(0,len(vars),4):
        sum,carry =  Add( vars[i+1],vars[i+2],vars[i+3])
        sum2,z0 = HalfAdd(vars[i],carry) 
        z2,z1 = Add(sum,sum2)
        sums.append([z0,z1,z2])"""


def _PopCountNaive(vars):
    if isinstance(vars, (int, float, complex)):
        return bin(vars).count("1")

    allconst = True
    for v in vars:
        aV = VAR(v)
        if not aV.isConst():
            allconst = False
            break

    if allconst:
        count = 0
        for v in vars:
            aV = VAR(v)
            if v.isConstTrue():
                count += 1
        n = 0
        output = []
        while 1 << n <= count:
            if (1 << n & count) == 0:
                output.append(false())
            else:
                output.append(true())
            n += 1
        return output

    maxwidth = math.ceil(math.log(len(vars), 2)) + 2
    output = []
    for i in range(maxwidth):
        output.append(false())
    # simple, suboptimal adder

    for v in vars:
        a = VAR(v)
        if not a.isConstFalse():
            output, carry = Add(output, [a])

    return output

    while len(vars) > 1:
        next = []
        next_carry = []
        for i in range(0, len(vars), 3):
            assert i < len(vars)
            a = vars[i]
            b = false()
            c = false()
            if i + 1 < len(vars):
                b = vars[i + 1]
            if i + 2 < len(vars):
                c = vars[i + 2]
            sum, carry = Add(a, b, c)
            next.append(sum)
            next.append(carry)
        vars = next


# Returns true if exactly one variable in the array is true.
def OneHot(*arrayVars):
    sum, carry = false(), false()
    eversum = false()
    evercarry = false()
    for x in arrayVars:
        sum, carry = Add(sum, x, carry)
        eversum = Or(eversum, sum)
        evercarry = Or(evercarry, carry)

    return And(eversum, ~evercarry)


def PopLE(constant, *arrayVars):
    return PopLT(constant + 1, *arrayVars)


# True if the population count of the vars in the array is exactly equal to the (integer) constant
# Note: this is an UNSIGNED comparison


def PopEq(compareTo, *arrayVars):
    if len(arrayVars) == 1 and isinstance(arrayVars[0], Iterable):
        arrayVars = arrayVars[0]
    # build up a tree of add operations
    if isinstance(compareTo, (bool, int, float, complex)):
        if compareTo < 0:
            print(
                "Warning: Passed negative constant to popcount methods, which expect unsigned arguments"
            )
            return false()
            # shouldn't really come up...
        elif compareTo == 0:
            return Not(Or(*arrayVars))
        elif compareTo == 1:
            return OneHot(*arrayVars)
        elif compareTo == len(arrayVars):
            return And(*arrayVars)
        elif compareTo > len(arrayVars):
            return false()

    count = PopCount(arrayVars)
    return Equal(count, compareTo)


# Note: this is an UNSIGNED comparison
def PopLT(compareTo, *arrayVars):
    if len(arrayVars) == 1 and isinstance(arrayVars[0], Iterable):
        arrayVars = arrayVars[0]
    # build up a tree of add operations
    if isinstance(compareTo, (bool, int, float, complex)):
        if compareTo < 0:
            print(
                "Warning: Passed negative constant to popcount methods, which expect unsigned arguments"
            )
            return false()
            # shouldn't really come up...
        elif compareTo == 0:
            return false()
        elif compareTo == 1:
            return Not(Or(arrayVars))
        elif compareTo == 2:
            return OneHot(*arrayVars)
        elif compareTo == len(arrayVars):
            return Not(And(arrayVars))
        elif compareTo > len(arrayVars):
            return true()
    count = PopCount(arrayVars)
    return LessThan(count, compareTo)


# Note: this is an UNSIGNED comparison
def PopGE(constant, *arrayVars):
    return PopGT(constant - 1, *arrayVars)


# True if the population count of the vars in the array is exactly equal to the (integer) constant
# Note: this is an UNSIGNED comparison
def PopGT(constant, *arrayVars):
    return Not(PopLE(constant, *arrayVars))
