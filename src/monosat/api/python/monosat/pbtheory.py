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

import os
import shutil
import sys
import time
import monosat.monosat_c
from monosat.logic import *
from monosat.manager import Manager
from monosat.monosat_c import Monosat, dimacs, Ineq
from tempfile import NamedTemporaryFile

debug = False


# Collects a set of graphs to encode together into a formula
class PBManager(metaclass=Manager):

    def setPB(self, pb):
        self.pb = pb

    def __init__(self):
        self.pb = MonosatPB()  # MinisatPlus()
        self.import_time = 0
        self.elapsed_time = 0

    def clear(self):
        if self.pb:
            self.pb.clear()

    def AssertLessThanPB(self, clause, val, weights=None):
        self.AssertPB(clause, val, "<", weights)
        # self.constraints.append((clause,val,'<',weights))

    def AssertGreaterThanPB(self, clause, val, weights=None):
        self.AssertPB(clause, val, ">", weights)
        # self.constraints.append((clause,val,'>',weights))

    def AssertLessEqPB(self, clause, val, weights=None):
        self.AssertPB(clause, val, "<=", weights)
        # self.constraints.append((clause,val,'<=',weights))

    def AssertGreaterEqPB(self, clause, val, weights=None):
        self.AssertPB(clause, val, ">=", weights)
        # self.constraints.append((clause,val,'>=',weights))

    def AssertEqualPB(self, clause, val, weights=None):
        # self.constraints.append((clause,val,'=',weights))
        self.AssertPB(clause, val, "=", weights)

    def AssertNotEqualPB(self, clause, val, weights=None):
        # self.constraints.append((clause,val,'=',weights))
        self.AssertPB(clause, val, "!=", weights)

    def AssertRangePB(self, clause, lowerBound, upperBound, weights=None):
        self.pb.AssertRangePB(clause, lowerBound, upperBound, weights)

    def AssertPB(self, clause, val, constraint, weights=None):
        self.pb.AssertPB(clause, val, constraint, weights)

    def twoSidedRangePB(self, clause, lowerBound, upperBound, weights=None):
        return self.pb.twoSidedRangePB(clause, lowerBound, upperBound, weights)

    def twoSidedPB(self, clause, val, constraint, weights=None, condition=None):
        return self.pb.twoSidedPB(clause, val, constraint, weights, condition)

    def conditionalPB(self, clause, val, constraint, weights=None, condition=None):
        return self.pb.conditionalPB(clause, val, constraint, weights, condition)

    def conditionalRangePB(
            self, clause, lowerBound, upperBound, weights=None, condition=None
    ):
        return self.pb.conditionalRangePB(
            clause, lowerBound, upperBound, weights, condition
        )

    def LessThanPB(self, clause, val, weights=None, condition=None):
        return self.conditionalPB(clause, val, "<", weights, condition)

    def GreaterThanPB(self, clause, val, weights=None, condition=None):
        # self.conditional_constraints.append((clause,val,condition,'>',weights))
        return self.conditionalPB(clause, val, ">", weights, condition)

    def LessEqPB(self, clause, val, weights=None, condition=None):
        # self.conditional_constraints.append((clause,val,condition,'<=',weights))
        return self.conditionalPB(clause, val, "<=", weights, condition)

    def GreaterEqPB(self, clause, val, weights=None, condition=None):
        # self.conditional_constraints.append((clause,val,condition,'>=',weights))
        return self.conditionalPB(clause, val, ">=", weights, condition)

    def EqualPB(self, clause, val, weights=None, condition=None):
        # self.conditional_constraints.append((clause,val,condition,'=',weights))
        return self.conditionalPB(clause, val, "=", weights, condition)

    def AssertAtMostOne(self, clause):
        new_clause = []
        # create fresh variables, so that theory vars can be used as arguments here.
        for l in clause:
            l2 = Var()
            new_clause.append(l2)
            AssertEq(l2, l)

        Monosat().AssertAtMostOne([l.getLit() for l in new_clause])

    def AssertExactlyOne(self, clause):
        AssertClause(clause)
        self.AssertAtMostOne(clause)

    def hasConstraints(self):
        return self.pb.hasConstraints()

    def write(self):
        self.pb.write()

    def flush(self):
        self.pb.write()
        self.pb.clear()


class MonosatPB:

    def __init__(self):
        self._monosat = monosat.monosat_c.Monosat()

    def clear(self):
        pass

    def AssertLessThanPB(self, clause, val, weights=None):
        self.AssertPB(clause, val, "<", weights)
        # self.constraints.append((clause,val,'<',weights))

    def AssertGreaterThanPB(self, clause, val, weights=None):
        self.AssertPB(clause, val, ">", weights)
        # self.constraints.append((clause,val,'>',weights))

    def AssertLessEqPB(self, clause, val, weights=None):
        self.AssertPB(clause, val, "<=", weights)
        # self.constraints.append((clause,val,'<=',weights))

    def AssertGreaterEqPB(self, clause, val, weights=None):
        self.AssertPB(clause, val, ">=", weights)
        # self.constraints.append((clause,val,'>=',weights))

    def AssertEqualPB(self, clause, val, weights=None):
        # self.constraints.append((clause,val,'=',weights))
        self.AssertPB(clause, val, "=", weights)

    def AssertRangePB(self, clause, lowerBound, upperBound, weights=None):
        if lowerBound == upperBound:
            self.AssertEqualPB(clause, lowerBound, weights)
        else:
            self.AssertGreaterEqPB(clause, lowerBound, weights)
            self.AssertLessEqPB(clause, upperBound, weights)
            # self.constraints.append((clause,lowerBound,'>=',weights))
            # self.constraints.append((clause,upperBound,'<=',weights))

    def getIneq(self, constraint):
        if constraint == "<":
            return Ineq.LT
        elif constraint == "<=":
            return Ineq.LEQ
        elif constraint == "=" or constraint == "==":
            return Ineq.EQ
        elif constraint == ">=":
            return Ineq.GEQ
        elif constraint == ">":
            return Ineq.GEQ
        else:
            raise Exception("Unknown operator " + str(constraint))

    def AssertPB(self, clause, val, constraint, weights=None):
        if constraint == "!=" or constraint == "<>":
            self.AssertNotEqualPB(clause, val, weights)
            return
        nclause = [l.getLit() for l in clause]
        nweights = weights
        if nweights is None:
            nweights = []
        nweights = list(nweights)
        while len(nweights) < len(nclause):
            nweights.append(1)

        # need all the variables to be in positive polarity...
        """
        for i,l in enumerate(clause):
            if weights is not None:
                w = weights[i]
            else:
                w=1
            if w==0:
                continue
            if l.isConstTrue():
                val-=w
                continue
            elif l.isConstFalse():
                continue

            if Monosat().isPositive(l.getLit()):
                nclause.append(l.getLit());
                nweights.append(w)
            else:
                nclause.append(Not(l).getLit())
                nweights.append(-w)
                val-=w
        """

        self._monosat.AssertPB(nclause, nweights, self.getIneq(constraint), val)

    def _negate(self, constraint):
        if constraint == "<":
            return ">="
        elif constraint == "<=":
            return ">"
        elif constraint == ">":
            return "<="
        elif constraint == ">=":
            return "<"
        elif constraint == "=":
            raise Exception("Cannot negate equality")
            # return '='; #!= is not directly supported in opb format, need to handle it separately
        else:
            raise Exception("Unknown operator " + constraint)

    def AssertNotEqualPB(self, clause, value, weights=None):
        Assert(Not(self.twoSidedPB(clause, value, "=", weights)))

    def twoSidedRangePB(self, clause, lowerBound, upperBound, weights=None):
        # Because opb doesn't support an inequality operator, we are instead going to make yet another choice, and then conditionally enforce that either > or < holds

        if lowerBound == upperBound:
            return self.twoSidedPB(clause, lowerBound, "=", weights)
        elif lowerBound is None:
            return self.twoSidedPB(clause, upperBound, "<=", weights)
        elif upperBound is None:
            return self.twoSidedPB(clause, lowerBound, ">=", weights)

        condition = Var()
        self.conditionalRangePB(clause, lowerBound, upperBound, weights, condition)

        # If the condition that we are within (inclusive) the range of [lowerBound,upperBound] is false, then we must be outside that range

        condition_g = (
            Var()
        )  # If variable condition_g is true, AND variable condition is false, then enforce the > condition.
        conditional_clause = []
        for v in clause:
            c = Var()
            conditional_clause.append(c)
            AssertImplies(And(Not(condition), condition_g), Equal(c, v))
        self.AssertPB(conditional_clause, upperBound, ">", weights)

        conditional_clause = []
        for v in clause:
            c = Var()
            conditional_clause.append(c)
            AssertImplies(And(Not(condition), Not(condition_g)), Equal(c, v))
        self.AssertPB(conditional_clause, lowerBound, "<", weights)

        return condition

    def twoSidedPB(self, clause, val, constraint, weights=None, condition=None):
        if constraint == "!=" or constraint == "<>":
            return Not(self.twoSidedPB(clause, val, "=", weights, condition))

        if condition is None:
            condition = Var()
        # elif not condition.isInput:
        #    v=Var()
        #    Assert(v==condition)
        #    condition=v

        nclause = []
        for v in clause:
            # if(not  v.isInput()):
            #    v2 = Var()
            #    AssertEq(v2,v)
            #    nclause.append(v2)
            # else:
            nclause.append(v)
        if weights is None:
            weights = []
        while len(weights) < len(clause):
            weights.append(1)

        if constraint != "=":
            # a two sided constraint is just two one-sided conditional constraints
            self.conditionalPB(nclause, val, constraint, weights, condition)
            self.conditionalPB(
                nclause, val, self._negate(constraint), weights, Not(condition)
            )
        else:
            # Is there a more efficient way to build a two sided equality constraint out of onesided constraints?
            self.conditionalPB(nclause, val, ">=", weights, condition)
            self.conditionalPB(nclause, val, "<=", weights, condition)
            self.conditionalPB(nclause, val, ">", weights, Not(condition))
            self.conditionalPB(nclause, val, "<", weights, Not(condition))

        return condition

    def conditionalPB_old(self, clause, val, constraint, weights=None, condition=None):
        if condition is None:
            condition = Var()

        conditional_clause = []
        for v in clause:
            c = Var()
            conditional_clause.append(c)
            Assert(Implies(condition, Equal(c, v)))
        self.AssertPB(conditional_clause, val, constraint, weights)
        return condition

    def conditionalPB(self, clause, val, constraint, weights=None, condition=None):
        if constraint == "!=" or constraint == "<>":
            v = Or(
                self.conditionalPB(clause, val, "<", weights),
                self.conditionalPB(clause, val, ">", weights),
            )
            if condition is not None:
                AssertEq(v, condition)
            return v

        if condition is None:
            condition = Var()
        else:
            v = Var()
            AssertEq(v, Not(condition))
            condition = v

        if weights is None:
            weights = []
        while len(weights) < len(clause):
            weights.append(1)
        nclause = []
        for v in clause:
            # For now, it is a limitation of the pb constraint solver that all variables it deals with must be input variables in the circuit.
            # So, we create them here if needed
            # if(not  v.isInput()):
            #    v2 = Var()
            #    AssertEq(v2,v)
            #    nclause.append(v2)
            # else:
            nclause.append(v)

        if (
                constraint == ">"
                or constraint == ">="
                or constraint == "="
                or constraint == "=="
        ):
            negWeightSum = 0
            for w in weights:
                if w < 0:
                    negWeightSum += abs(w)

            total = negWeightSum + val
            if constraint == ">":
                total += 1

            nclause.append(condition)
            weights.append(total)

            if constraint == ">" or constraint == ">=":
                self.AssertPB(nclause, val, constraint, weights)
            else:
                self.AssertPB(nclause, val, ">=", weights)
            nclause.pop()
            weights.pop()
        if (
                constraint == "<"
                or constraint == "<="
                or constraint == "="
                or constraint == "=="
        ):
            posWeightSum = 0
            for w in weights:
                if w > 0:
                    posWeightSum += abs(w)

            total = posWeightSum + val
            if constraint == "<":
                total += 1

            nclause.append(condition)
            weights.append(-total)
            if constraint == "<" or constraint == "<=":
                self.AssertPB(nclause, val, constraint, weights)
            else:
                self.AssertPB(nclause, val, "<=", weights)
            nclause.pop()
            weights.pop()
        return Not(
            condition
        )  # The property is enforced if the condition variable is false

    def conditionalRangePB(
            self, clause, lowerBound, upperBound, weights=None, condition=None
    ):

        if lowerBound == upperBound:
            return self.EqualPB(clause, lowerBound, weights, condition)
        elif lowerBound is None:
            return self.conditionalPB(clause, upperBound, "<=", weights, condition)
        elif upperBound is None:
            return self.conditionalPB(clause, lowerBound, ">=", weights, condition)

        if condition is None:
            condition = Var()
        conditional_clause = []
        for v in clause:
            c = Var()
            conditional_clause.append(c)
            AssertImplies(condition, Equal(c, v))
        self.AssertPB(conditional_clause, lowerBound, ">=", weights)
        self.AssertPB(conditional_clause, upperBound, "<=", weights)
        return condition

    def LessThanPB(self, clause, val, weights=None, condition=None):
        return self.conditionalPB(clause, val, "<", weights, condition)

    def GreaterThanPB(self, clause, val, weights=None, condition=None):
        # self.conditional_constraints.append((clause,val,condition,'>',weights))
        return self.conditionalPB(clause, val, ">", weights, condition)

    def LessEqPB(self, clause, val, weights=None, condition=None):
        # self.conditional_constraints.append((clause,val,condition,'<=',weights))
        return self.conditionalPB(clause, val, "<=", weights, condition)

    def GreaterEqPB(self, clause, val, weights=None, condition=None):
        # self.conditional_constraints.append((clause,val,condition,'>=',weights))
        return self.conditionalPB(clause, val, ">=", weights, condition)

    def EqualPB(self, clause, val, weights=None, condition=None):
        # self.conditional_constraints.append((clause,val,condition,'=',weights))
        return self.conditionalPB(clause, val, "=", weights, condition)

    def hasConstraints(self):
        return False

    # not required, will be called automatically inside monosat
    def write(self):
        self._monosat.flushPB()


class MonosatTheoryPB:

    def __init__(self):
        self.constraints = []

    def clear(self):
        self.constraints = []

    def hasConstraints(self):
        return len(self.constraints) > 0

    def AssertLessThanPB(self, clause, val, weights=None):
        self.AssertPB(clause, val, "<", weights)

    def AssertGreaterThanPB(self, clause, val, weights=None):
        self.AssertPB(clause, val, ">", weights)

    def AssertLessEqPB(self, clause, val, weights=None):
        self.AssertPB(clause, val, "<=", weights)

    def AssertGreaterEqPB(self, clause, val, weights=None):
        self.AssertPB(clause, val, ">=", weights)

    def AssertEqualPB(self, clause, val, weights=None):
        self.AssertPB(clause, val, "=", weights)

    def AssertRangePB(self, clause, lowerBound, upperBound, weights=None):

        if lowerBound == upperBound:
            self.AssertEqualPB(clause, lowerBound, weights)
        else:
            self.AssertGreaterEqPB(clause, lowerBound, weights)
            self.AssertLessEqPB(clause, upperBound, weights)

    def AssertPB(self, clause, val, constraint, weights=None):
        nclause = list(clause)
        self.constraints.append(
            (
                nclause,
                val,
                constraint,
                list(weights) if weights is not None else [],
                None,
                False,
            )
        )

    def twoSidedRangePB(
            self, clause, lowerBound, upperBound, weights=None, condition=None
    ):
        if condition is None:
            condition = Var()
        # elif not condition.isInput:
        #    v=Var()
        #    Assert(v==condition)
        #    condition=v

        if lowerBound == upperBound:
            return self.twoSided(clause, lowerBound, "=", weights, condition)

        c1 = self.twoSidedPB(clause, lowerBound, ">=", weights)
        c2 = self.twoSidedPB(clause, upperBound, "<=", weights)
        if condition is None:
            condition = And(c1, c2)
        else:
            AssertEq(condition, And(c1, c2))

        return condition

    def twoSidedPB(self, clause, val, constraint, weights=None, condition=None):
        if condition is None:
            condition = Var()
        # elif not condition.isInput:
        #    v=Var()
        #    Assert(v==condition)
        #    condition=v

        nclause = []
        for v in clause:
            # if(not  v.isInput()):
            #    v2 = Var()
            #    Assert(Equal(v2,v))
            #    nclause.append(v2)
            # else:
            nclause.append(v)
        self.constraints.append(
            (
                nclause,
                val,
                constraint,
                list(weights) if weights is not None else [],
                condition,
                False,
            )
        )

        return condition

    def conditionalPB(self, clause, val, constraint, weights=None, condition=None):
        if condition is None:
            condition = Var()
        # elif not condition.isInput:
        #    v=Var()
        #    Assert(v==condition)
        #    condition=v

        nclause = []
        for v in clause:
            # if(not  v.isInput()):
            #    v2 = Var()
            #    Assert(Equal(v2,v))
            #    nclause.append(v2)
            # else:
            nclause.append(v)
        self.constraints.append(
            (
                nclause,
                val,
                constraint,
                list(weights) if weights is not None else [],
                condition,
                True,
            )
        )

        return condition

    def conditionalRangePB(
            self, clause, lowerBound, upperBound, weights=None, condition=None
    ):
        if lowerBound == upperBound:
            return self.conditionalPB(clause, lowerBound, "=", weights, condition)

        c1 = self.twoSidedPB(clause, lowerBound, ">=", weights)
        c2 = self.twoSidedPB(clause, upperBound, "<=", weights)
        if condition is None:
            condition = And(c1, c2)
        else:
            Assert(Equal(condition, And(c1, c2)))

        return condition

    def LessThanPB(self, clause, val, weights=None, condition=None):
        return self.conditionalPB(clause, val, "<", weights, condition)

    def GreaterThanPB(self, clause, val, weights=None, condition=None):
        # self.conditional_constraints.append((clause,val,condition,'>',weights))
        return self.conditionalPB(clause, val, ">", weights, condition)

    def LessEqPB(self, clause, val, weights=None, condition=None):
        # self.conditional_constraints.append((clause,val,condition,'<=',weights))
        return self.conditionalPB(clause, val, "<=", weights, condition)

    def GreaterEqPB(self, clause, val, weights=None, condition=None):
        # self.conditional_constraints.append((clause,val,condition,'>=',weights))
        return self.conditionalPB(clause, val, ">=", weights, condition)

    def EqualPB(self, clause, val, weights=None, condition=None):
        # self.conditional_constraints.append((clause,val,condition,'=',weights))
        return self.conditionalPB(clause, val, "=", weights, condition)

    def write(self, filename):
        f = open(filename, "a")
        for (
                clause,
                val,
                constraint,
                weights,
                conditionVar,
                oneSided,
        ) in self.constraints:
            assert len(clause) == len(weights) or len(weights) == 0
            assert len(clause) > 0
            f.write("pb_lt " + str(len(clause)) + " ")
            for v in clause:
                f.write(str(v.getInputLiteral()) + " ")
            f.write(str(len(weights)) + " ")
            for w in weights:
                f.write(str(w) + " ")

            f.write(constraint + " " + str(val) + " ")
            if conditionVar is None:
                f.write("0")
            else:
                if oneSided:
                    f.write("1 ")
                else:
                    f.write("2 ")

                f.write(str(conditionVar.getInputLiteral()))
            f.write("\n")
            pass


class MinisatPlus:

    def __init__(self):
        self.constraints = []

    def clear(self):
        self.constraints = []

    def AssertLessThanPB(self, clause, val, weights=None):
        self.AssertPB(clause, val, "<", weights)
        # self.constraints.append((clause,val,'<',weights))

    def AssertGreaterThanPB(self, clause, val, weights=None):
        self.AssertPB(clause, val, ">", weights)
        # self.constraints.append((clause,val,'>',weights))

    def AssertLessEqPB(self, clause, val, weights=None):
        self.AssertPB(clause, val, "<=", weights)
        # self.constraints.append((clause,val,'<=',weights))

    def AssertGreaterEqPB(self, clause, val, weights=None):
        self.AssertPB(clause, val, ">=", weights)
        # self.constraints.append((clause,val,'>=',weights))

    def AssertEqualPB(self, clause, val, weights=None):
        # self.constraints.append((clause,val,'=',weights))
        self.AssertPB(clause, val, "=", weights)

    def AssertRangePB(self, clause, lowerBound, upperBound, weights=None):
        if lowerBound == upperBound:
            self.AssertEqualPB(clause, lowerBound, weights)
        else:
            self.AssertGreaterEqPB(clause, lowerBound, weights)
            self.AssertLessEqPB(clause, upperBound, weights)
            # self.constraints.append((clause,lowerBound,'>=',weights))
            # self.constraints.append((clause,upperBound,'<=',weights))

    def AssertPB(self, clause, val, constraint, weights=None):
        if constraint == "!=" or constraint == "<>":
            self.AssertNotEqualPB(clause, val, weights)
            return
        nclause = []
        nweights = []
        # need all the variables to be in positive polarity...
        for i, l in enumerate(clause):
            if weights is not None:
                w = weights[i]
            else:
                w = 1
            if w == 0:
                continue
            if l.isConstTrue():
                val -= w
                continue
            elif l.isConstFalse():
                continue

            if Monosat().isPositive(l.getLit()):
                nclause.append(l)
                nweights.append(w)
            else:
                nclause.append(Not(l))
                nweights.append(-w)
                val -= w

        self.constraints.append((nclause, val, constraint, nweights))

    def _negate(self, constraint):
        if constraint == "<":
            return ">="
        elif constraint == "<=":
            return ">"
        elif constraint == ">":
            return "<="
        elif constraint == ">=":
            return "<"
        elif constraint == "=":
            raise Exception("Cannot negate equality")
            # return '='; #!= is not directly supported in opb format, need to handle it separately
        else:
            raise Exception("Unknown operator " + constraint)

    def AssertNotEqualPB(self, clause, value, weights=None):
        Assert(Not(self.twoSidedPB(clause, value, "=", weights)))

    def twoSidedRangePB(self, clause, lowerBound, upperBound, weights=None):
        # Because opb doesn't support an inequality operator, we are instead going to make yet another choice, and then conditionally enforce that either > or < holds

        if lowerBound == upperBound:
            return self.twoSidedPB(clause, lowerBound, "=", weights)
        elif lowerBound is None:
            return self.twoSidedPB(clause, upperBound, "<=", weights)
        elif upperBound is None:
            return self.twoSidedPB(clause, lowerBound, ">=", weights)

        condition = Var()
        self.conditionalRangePB(clause, lowerBound, upperBound, weights, condition)

        # If the condition that we are within (inclusive) the range of [lowerBound,upperBound] is false, then we must be outside that range

        condition_g = (
            Var()
        )  # If variable condition_g is true, AND variable condition is false, then enforce the > condition.
        conditional_clause = []
        for v in clause:
            c = Var()
            conditional_clause.append(c)
            AssertImplies(And(Not(condition), condition_g), Equal(c, v))
        self.AssertPB(conditional_clause, upperBound, ">", weights)

        conditional_clause = []
        for v in clause:
            c = Var()
            conditional_clause.append(c)
            AssertImplies(And(Not(condition), Not(condition_g)), Equal(c, v))
        self.AssertPB(conditional_clause, lowerBound, "<", weights)

        return condition

    def twoSidedPB(self, clause, val, constraint, weights=None, condition=None):
        if constraint == "!=" or constraint == "<>":
            return Not(self.twoSidedPB(clause, val, "=", weights, condition))

        if condition is None:
            condition = Var()
        # elif not condition.isInput:
        #    v=Var()
        #    Assert(v==condition)
        #    condition=v

        nclause = []
        for v in clause:
            # if(not  v.isInput()):
            #    v2 = Var()
            #    AssertEq(v2,v)
            #    nclause.append(v2)
            # else:
            nclause.append(v)
        if weights is None:
            weights = []
        while len(weights) < len(clause):
            weights.append(1)

        """
        conditional_clause=[]
        for v in clause:
            c = Var()
            conditional_clause.append(c)
            Assert(Implies(condition,Equal(c,v)))
        self.AssertPB(conditional_clause,val,constraint,weights)
        if(constraint != '='):
            conditional_clause=[]
            for v in clause:
                c = Var()
                conditional_clause.append(c)
                Assert(Implies(Not(condition),Equal(c,v)))        
            self.AssertPB(conditional_clause,val,self._negate(constraint),weights)
            #This doesn't work. If constraint is (1 a 1 b <= 2), then the second side becomes
            #1 c 1 d > 2, which is not UNSAT.
        else:
            #Because opb doesn't support an inequality operator, we are instead going to make yet another choice, and then conditionally enforce that either > or < holds
            condition_g = Var() #If variable condition_g is true, AND variable condition is false, then enforce the > condition.
            conditional_clause=[]
            for v in clause:
                c = Var()
                conditional_clause.append(c)
                Assert(Implies(And(Not(condition),condition_g),Equal(c,v)))        
            self.AssertPB(conditional_clause,val,'>',weights)            
            
            conditional_clause=[]
            for v in clause:
                c = Var()
                conditional_clause.append(c)
                Assert(Implies(And(Not(condition),Not(condition_g)),Equal(c,v)))        
            self.AssertPB(conditional_clause,val,'<',weights)              
        """
        if constraint != "=":
            # a two sided constraint is just two one-sided conditional constraints
            self.conditionalPB(nclause, val, constraint, weights, condition)
            self.conditionalPB(
                nclause, val, self._negate(constraint), weights, Not(condition)
            )
        else:
            # Is there a more efficient way to build a two sided equality constraint out of onesided constraints?
            self.conditionalPB(nclause, val, ">=", weights, condition)
            self.conditionalPB(nclause, val, "<=", weights, condition)
            self.conditionalPB(nclause, val, ">", weights, Not(condition))
            self.conditionalPB(nclause, val, "<", weights, Not(condition))

        return condition

    def conditionalPB_old(self, clause, val, constraint, weights=None, condition=None):
        if condition is None:
            condition = Var()

        conditional_clause = []
        for v in clause:
            c = Var()
            conditional_clause.append(c)
            Assert(Implies(condition, Equal(c, v)))
        self.AssertPB(conditional_clause, val, constraint, weights)
        return condition

    def conditionalPB(self, clause, val, constraint, weights=None, condition=None):
        if constraint == "!=" or constraint == "<>":
            v = Or(
                self.conditionalPB(clause, val, "<", weights),
                self.conditionalPB(clause, val, ">", weights),
            )
            if condition is not None:
                AssertEq(v, condition)
            return v

        if condition is None:
            condition = Var()
        else:
            v = Var()
            AssertEq(v, Not(condition))
            condition = v

        """
        conditional_clause=[]
        for v in clause:
            c = Var()
            conditional_clause.append(c)
            Assert(Implies(condition,Equal(c,v)))
        self.AssertPB(conditional_clause,val,constraint,weights)
        """
        if weights is None:
            weights = []
        while len(weights) < len(clause):
            weights.append(1)
        nclause = []
        for v in clause:
            # For now, it is a limitation of the pb constraint solver that all variables it deals with must be input variables in the circuit.
            # So, we create them here if needed
            # if(not  v.isInput()):
            #    v2 = Var()
            #    AssertEq(v2,v)
            #    nclause.append(v2)
            # else:
            nclause.append(v)

        if (
                constraint == ">"
                or constraint == ">="
                or constraint == "="
                or constraint == "=="
        ):
            negWeightSum = 0
            for w in weights:
                if w < 0:
                    negWeightSum += abs(w)

            total = negWeightSum + val
            if constraint == ">":
                total += 1

            nclause.append(condition)
            weights.append(total)

            if constraint == ">" or constraint == ">=":
                self.AssertPB(nclause, val, constraint, weights)
            else:
                self.AssertPB(nclause, val, ">=", weights)
            nclause.pop()
            weights.pop()
        if (
                constraint == "<"
                or constraint == "<="
                or constraint == "="
                or constraint == "=="
        ):
            posWeightSum = 0
            for w in weights:
                if w > 0:
                    posWeightSum += abs(w)

            total = posWeightSum + val
            if constraint == "<":
                total += 1

            nclause.append(condition)
            weights.append(-total)
            if constraint == "<" or constraint == "<=":
                self.AssertPB(nclause, val, constraint, weights)
            else:
                self.AssertPB(nclause, val, "<=", weights)
            nclause.pop()
            weights.pop()
        return Not(
            condition
        )  # The property is enforced if the condition variable is false

    def conditionalRangePB(
            self, clause, lowerBound, upperBound, weights=None, condition=None
    ):

        if lowerBound == upperBound:
            return self.EqualPB(clause, lowerBound, weights, condition)
        elif lowerBound is None:
            return self.conditionalPB(clause, upperBound, "<=", weights, condition)
        elif upperBound is None:
            return self.conditionalPB(clause, lowerBound, ">=", weights, condition)

        if condition is None:
            condition = Var()
        conditional_clause = []
        for v in clause:
            c = Var()
            conditional_clause.append(c)
            AssertImplies(condition, Equal(c, v))
        self.AssertPB(conditional_clause, lowerBound, ">=", weights)
        self.AssertPB(conditional_clause, upperBound, "<=", weights)
        return condition

    def LessThanPB(self, clause, val, weights=None, condition=None):
        return self.conditionalPB(clause, val, "<", weights, condition)

    def GreaterThanPB(self, clause, val, weights=None, condition=None):
        # self.conditional_constraints.append((clause,val,condition,'>',weights))
        return self.conditionalPB(clause, val, ">", weights, condition)

    def LessEqPB(self, clause, val, weights=None, condition=None):
        # self.conditional_constraints.append((clause,val,condition,'<=',weights))
        return self.conditionalPB(clause, val, "<=", weights, condition)

    def GreaterEqPB(self, clause, val, weights=None, condition=None):
        # self.conditional_constraints.append((clause,val,condition,'>=',weights))
        return self.conditionalPB(clause, val, ">=", weights, condition)

    def EqualPB(self, clause, val, weights=None, condition=None):
        # self.conditional_constraints.append((clause,val,condition,'=',weights))
        return self.conditionalPB(clause, val, "=", weights, condition)

    def hasConstraints(self):
        return len(self.constraints) > 0

    def write(self):
        if len(self.constraints) == 0:
            return

        tmpfile = NamedTemporaryFile(delete=False, suffix=".opb")
        tmpopb = tmpfile.name
        tmpfile.close()

        tmpfile = NamedTemporaryFile(delete=False, suffix=".cnf")
        tmpcnf = tmpfile.name
        tmpfile.close()
        longest_constraint = 0

        try:
            minisat_plus_path = shutil.which("monosatpb")
        except:
            minisat_plus_path = which("monosatpb")
        if minisat_plus_path is None:
            raise RuntimeError(
                "In order to use PB constraints, monosatpb must be installed and on the path (see README).\n"
            )

        print(
            "Encoding pb constraints using %s using temporary files %s and %s "
            % (minisat_plus_path, tmpopb, tmpcnf)
        )
        fopb = open(tmpopb, "w")
        invarmap = dict()
        outvarmap = dict()
        varmap = dict()
        nvars = 0
        for (clause, val, op, weights) in self.constraints:
            if not isinstance(val, int):
                raise TypeError(
                    "PB constraints weights must compare to integers, but found "
                    + str(type(val))
                )
            if len(clause) == 0:
                clause = [false()]
            if len(clause) > longest_constraint:
                longest_constraint = len(clause)
            for v in clause:
                l = dimacs(v.getLit())
                if l not in invarmap:
                    nvars += 1
                    invarmap[l] = nvars
                    varmap[nvars] = v.getLit()

        n_pbs = 0
        nvars += 1
        fopb.write(
            "* #variable= "
            + str(nvars)
            + " #constraint= "
            + str(len(self.constraints))
            + "\n"
        )
        for (clause, val, op, weights) in self.constraints:
            if weights is None:
                weights = []
            if len(weights) > 0:
                pass
            if len(clause) == 0:
                clause = [false()]
            n_pbs += 1
            while len(weights) < len(clause):
                weights.append(1)  # Default weight

            for (v, w) in zip(clause, weights):
                l = dimacs(v.getLit())
                fopb.write(
                    ("+" if w >= 0 else "") + str(w) + " x" + str(invarmap[l]) + " "
                )
            fopb.write(op + " " + str(val) + " ;\n")

        fopb.close()

        os.system(
            "ulimit -s %d" % (32768 + longest_constraint)
        )  # need to give minisat+ more stack space, or it will crash when trying to create BDDs...
        opts = ""
        if longest_constraint >= 5000:
            opts = (
                "-ca"
            )  # Switch to simple adder encoding if there are very large pb constraints, or else minisat+ will be inordinately slow
            print(
                "Note: using minisat+ %s encoding because largest pb constraint has %d > 5000 arguments"
                % (opts, longest_constraint)
            )
        # This call is designed for the original Minisat+; it can easily be adapted for other pb constraint->cnf encoders.

        status = os.system(
            "%s %s -a -s -v0 -cnf=" % (minisat_plus_path, opts) + tmpcnf + " " + tmpopb
        )  # a turns off ansi-codes in output, -s turns off sat competition output, and -v0 sets verbosity to 0.
        if status != 0:
            raise RuntimeError(
                "In order to use PB constraints, minisat+ (1.0) must be installed and on the path.\nEither it wasn't found, or it was but still returned an exit code of %d\n"
                % (status)
            )

        print("Importing %d pseudoboolean constraints into Monosat..." % (n_pbs))
        t = time.process_time()

        Monosat().comment("pseudoboolean constraints")
        n_cls = 0
        convcnf = open(tmpcnf, "r")
        for line in convcnf.readlines():
            if line.startswith("p cnf"):  # Header
                header = line.split()
                # vars = int(header[2])
            elif line.startswith("c"):
                pass
            elif len(line.split()) == 0:
                pass
            else:
                n_cls += 1
                clause = list(map(int, line.split()))
                newclause = []
                assert len(clause) > 1
                assert clause[-1] == 0

                for i in range(len(clause)):
                    l = clause[i]
                    if l == 0:
                        continue
                    v = abs(l)
                    if v not in varmap:
                        varmap[v] = Monosat().newLit()

                    newl = varmap[v]
                    if l < 0:
                        newl = Monosat().Not(newl)
                    newclause.append(newl)
                    # fcnf.write(str(v) + " ")
                # fcnf.write("0\n")
                Monosat().addClause(newclause)
        Monosat().comment("end of pseudoboolean constraints")

        os.remove(tmpopb)
        os.remove(tmpcnf)
        PBManager().import_time += time.process_time() - t
        print("Imported pseudoboolean constraints into Monosat (%d clauses)" % (n_cls))


class PBSugar:

    def __init__(self):
        self.constraints = []

    def clear(self):
        self.constraints = []

    def hasConstraints(self):
        return len(self.constraints) > 0

    def AssertLessThanPB(self, clause, val, weights=None):
        self.AssertPB(clause, val, "<", weights)
        # self.constraints.append((clause,val,'<',weights))

    def AssertGreaterThanPB(self, clause, val, weights=None):
        self.AssertPB(clause, val, ">", weights)
        # self.constraints.append((clause,val,'>',weights))

    def AssertLessEqPB(self, clause, val, weights=None):
        self.AssertPB(clause, val, "<=", weights)
        # self.constraints.append((clause,val,'<=',weights))

    def AssertGreaterEqPB(self, clause, val, weights=None):
        self.AssertPB(clause, val, ">=", weights)
        # self.constraints.append((clause,val,'>=',weights))

    def AssertEqualPB(self, clause, val, weights=None):
        # self.constraints.append((clause,val,'=',weights))
        self.AssertPB(clause, val, "=", weights)

    def AssertRangePB(self, clause, lowerBound, upperBound, weights=None):
        if lowerBound == upperBound:
            self.AssertEqualPB(clause, lowerBound, weights)
        else:
            self.AssertGreaterEqPB(clause, lowerBound, weights)
            self.AssertLessEqPB(clause, upperBound, weights)
            # self.constraints.append((clause,lowerBound,'>=',weights))
            # self.constraints.append((clause,upperBound,'<=',weights))

    def AssertPB(self, clause, val, constraint, weights=None):
        nclause = list(clause)
        self.constraints.append(
            (nclause, val, constraint, list(weights) if weights is not None else [])
        )

    def _negate(self, constraint):
        if constraint == "<":
            return ">="
        elif constraint == "<=":
            return ">"
        elif constraint == ">":
            return "<="
        elif constraint == ">=":
            return "<"
        elif constraint == "=":
            raise Exception("Cannot negate equality")
            # return '='; #!= is not directly supported in opb format, need to handle it separately
        else:
            raise Exception("Unknown operator " + constraint)

    def twoSidedRangePB(self, clause, lowerBound, upperBound, weights=None):
        # Because opb doesn't support an inequality operator, we are instead going to make yet another choice, and then conditionally enforce that either > or < holds

        if lowerBound == upperBound:
            return self.twoSidedPB(clause, lowerBound, "=", weights)
        elif lowerBound is None:
            return self.twoSidedPB(clause, upperBound, "<=", weights)
        elif upperBound is None:
            return self.twoSidedPB(clause, lowerBound, ">=", weights)

        condition = Var()
        self.conditionalRangePB(clause, lowerBound, upperBound, weights, condition)

        # If the condition that we are within (inclusive) the range of [lowerBound,upperBound] is false, then we must be outside that range

        condition_g = (
            Var()
        )  # If variable condition_g is true, AND variable condition is false, then enforce the > condition.
        conditional_clause = []
        for v in clause:
            c = Var()
            conditional_clause.append(c)
            Assert(Implies(And(Not(condition), condition_g), Equal(c, v)))
        self.AssertPB(conditional_clause, upperBound, ">", weights)

        conditional_clause = []
        for v in clause:
            c = Var()
            conditional_clause.append(c)
            Assert(Implies(And(Not(condition), Not(condition_g)), Equal(c, v)))
        self.AssertPB(conditional_clause, lowerBound, "<", weights)

        return condition

    def twoSidedPB(self, clause, val, constraint, weights=None):
        condition = Var()
        conditional_clause = []
        for v in clause:
            c = Var()
            conditional_clause.append(c)
            Assert(Implies(condition, Equal(c, v)))
        self.AssertPB(conditional_clause, val, constraint, weights)
        if constraint != "=":
            conditional_clause = []
            for v in clause:
                c = Var()
                conditional_clause.append(c)
                Assert(Implies(Not(condition), Equal(c, v)))
            self.AssertPB(conditional_clause, val, self._negate(constraint), weights)
        else:
            # Because opb doesn't support an inequality operator, we are instead going to make yet another choice, and then conditionally enforce that either > or < holds
            condition_g = (
                Var()
            )  # If variable condition_g is true, AND variable condition is false, then enforce the > condition.
            conditional_clause = []
            for v in clause:
                c = Var()
                conditional_clause.append(c)
                Assert(Implies(And(Not(condition), condition_g), Equal(c, v)))
            self.AssertPB(conditional_clause, val, ">", weights)

            conditional_clause = []
            for v in clause:
                c = Var()
                conditional_clause.append(c)
                Assert(Implies(And(Not(condition), Not(condition_g)), Equal(c, v)))
            self.AssertPB(conditional_clause, val, "<", weights)

        return condition

    def conditionalPB_old(self, clause, val, constraint, weights=None, condition=None):
        if condition is None:
            condition = Var()

        conditional_clause = []
        for v in clause:
            c = Var()
            conditional_clause.append(c)
            Assert(Implies(condition, Equal(c, v)))
        self.AssertPB(conditional_clause, val, constraint, weights)
        return condition

    def conditionalPB(self, clause, val, constraint, weights=None, condition=None):
        if condition is None:
            condition = Var()
        else:
            v = Var()
            AssertEq(v, Not(condition))
            condition = v

        """
        conditional_clause=[]
        for v in clause:
            c = Var()
            conditional_clause.append(c)
            Assert(Implies(condition,Equal(c,v)))
        self.AssertPB(conditional_clause,val,constraint,weights)
        """
        if weights is None:
            weights = []
        while len(weights) < len(clause):
            weights.append(1)
        nclause = []
        for v in clause:
            # For now, it is a limitation of the pb constraint solver that all variables it deals with must be input variables in the circuit.
            # So, we create them here if needed
            # if(not  v.isInput()):
            #    v2 = Var()
            #    Assert(Equal(v2,v))
            #    nclause.append(v2)
            # else:
            nclause.append(v)

        if (
                constraint == ">"
                or constraint == ">="
                or constraint == "="
                or constraint == "=="
        ):
            negWeightSum = 0
            for w in weights:
                if w < 0:
                    negWeightSum += abs(w)

            total = negWeightSum + val
            if constraint == ">":
                total += 1

            nclause.append(condition)
            weights.append(total)

            if constraint == ">" or constraint == ">=":
                self.AssertPB(nclause, val, constraint, weights)
            else:
                self.AssertPB(nclause, val, ">=", weights)
            nclause.pop()
            weights.pop()
        if (
                constraint == "<"
                or constraint == "<="
                or constraint == "="
                or constraint == "=="
        ):
            posWeightSum = 0
            for w in weights:
                if w > 0:
                    posWeightSum += abs(w)

            total = posWeightSum + val
            if constraint == "<":
                total += 1

            nclause.append(condition)
            weights.append(-total)
            if constraint == "<" or constraint == "<=":
                self.AssertPB(nclause, val, constraint, weights)
            else:
                self.AssertPB(nclause, val, "<=", weights)
            nclause.pop()
            weights.pop()
        return Not(
            condition
        )  # The property is enforced if the condition variable is false

    def conditionalRangePB(
            self, clause, lowerBound, upperBound, weights=None, condition=None
    ):

        if lowerBound == upperBound:
            return self.EqualPB(clause, lowerBound, weights, condition)
        elif lowerBound is None:
            return self.conditionalPB(clause, upperBound, "<=", weights, condition)
        elif upperBound is None:
            return self.conditionalPB(clause, lowerBound, ">=", weights, condition)

        if condition is None:
            condition = Var()
        conditional_clause = []
        for v in clause:
            c = Var()
            conditional_clause.append(c)
            Assert(Implies(condition, Equal(c, v)))
        self.AssertPB(conditional_clause, lowerBound, ">=", weights)
        self.AssertPB(conditional_clause, upperBound, "<=", weights)
        return condition

    def LessThanPB(self, clause, val, weights=None, condition=None):
        return self.conditionalPB(clause, val, "<", weights, condition)

    def GreaterThanPB(self, clause, val, weights=None, condition=None):
        # self.conditional_constraints.append((clause,val,condition,'>',weights))
        return self.conditionalPB(clause, val, ">", weights, condition)

    def LessEqPB(self, clause, val, weights=None, condition=None):
        # self.conditional_constraints.append((clause,val,condition,'<=',weights))
        return self.conditionalPB(clause, val, "<=", weights, condition)

    def GreaterEqPB(self, clause, val, weights=None, condition=None):
        # self.conditional_constraints.append((clause,val,condition,'>=',weights))
        return self.conditionalPB(clause, val, ">=", weights, condition)

    def EqualPB(self, clause, val, weights=None, condition=None):
        # self.conditional_constraints.append((clause,val,condition,'=',weights))
        return self.conditionalPB(clause, val, "=", weights, condition)

    def write(self, cnf):
        if len(self.constraints) == 0:
            return

        print("Encoding pb constraints")
        fopb = open("convert.opb", "w")
        nvars = 0
        for (clause, val, op) in self.constraints:
            for v in clause:
                if v.getInputLiteral() > nvars:
                    nvars = v.getInputLiteral()
        nvars += 1
        fopb.write(
            "* #variable= "
            + str(nvars)
            + " #constraint= "
            + str(len(self.constraints))
            + "\n"
        )
        for (clause, val, op) in self.constraints:
            for v in clause:
                fopb.write("+1 x" + str(v.getInputLiteral()) + " ")
            fopb.write(op + " " + str(val) + " ;\n")

        fopb.close()
        os.system(
            "./pbsugar -map convert.map -n -sat convert.cnf -jar ./pbsugar-v1-1-1.jar "
            + "convert.opb"
        )
        fmap = open("convert.map", "r")
        count = int(fmap.readline())
        varmap = dict()
        maxvarmap = 0
        for line in fmap.readlines():
            line = line.split()
            invar = int(line[0][1:])
            outvar = int(line[1])
            varmap[invar] = outvar
            maxvarmap = max(maxvarmap, outvar)
        fmap.close()

        # now find the largest var in the cnf, so we can cleanly place the pb cnf in its own namespace
        fcnf = open(cnf, "r")
        vars = 0
        for line in fcnf.readlines():
            if line.startswith("p cnf"):  # Header
                header = line.split()
                vars = int(header[2])
        fcnf.close()
        fcnf = open(cnf, "a")
        fcnf.write("c pseudoboolean constraints:\n")
        vars = max(vars, maxvarmap)
        nextVar = vars + 1

        convcnf = open("convert.cnf", "r")
        for line in convcnf.readlines():
            if line.startswith("p cnf"):  # Header
                header = line.split()
                # vars = int(header[2])
            elif line.startswith("c"):
                pass
            elif len(line.split()) == 0:
                pass
            else:
                clause = list(map(int, line.split()))
                assert len(clause) > 1
                assert clause[-1] == 0

                for l in clause:
                    if l == 0:
                        continue
                    v = abs(l)
                    if v not in varmap:
                        varmap[v] = nextVar
                        nextVar += 1
                    v = varmap[v]
                    if l < 0:
                        v = -v

                    fcnf.write(str(v) + " ")
                fcnf.write("0\n")
        fcnf.write("c end of pseudoboolean constraints\n")
        fcnf.close()


def AssertLessThanPB(clause, val, weights=None):
    PBManager().AssertPB(clause, val, "<", weights)


def AssertGreaterThanPB(clause, val, weights=None):
    PBManager().AssertPB(clause, val, ">", weights)


def AssertLessEqPB(clause, val, weights=None):
    PBManager().AssertPB(clause, val, "<=", weights)


def AssertGreaterEqPB(clause, val, weights=None):
    PBManager().AssertPB(clause, val, ">=", weights)


def AssertEqualPB(clause, val, weights=None):
    PBManager().AssertPB(clause, val, "=", weights)


def AssertNotEqualPB(clause, val, weights=None):
    PBManager().AssertPB(clause, val, "!=", weights)


def AssertRangePB(clause, lowerBound, upperBound, weights=None):
    PBManager().pb.AssertRangePB(clause, lowerBound, upperBound, weights)


def AssertPB(clause, val, constraint, weights=None):
    PBManager().pb.AssertPB(clause, val, constraint, weights)


def twoSidedRangePB(clause, lowerBound, upperBound, weights=None):
    return PBManager().pb.twoSidedRangePB(clause, lowerBound, upperBound, weights)


def twoSidedPB(clause, val, constraint, weights=None, condition=None):
    return PBManager().pb.twoSidedPB(clause, val, constraint, weights, condition)


def conditionalPB(clause, val, constraint, weights=None, condition=None):
    return PBManager().pb.conditionalPB(clause, val, constraint, weights, condition)


def conditionalRangePB(clause, lowerBound, upperBound, weights=None, condition=None):
    return PBManager().pb.conditionalRangePB(
        clause, lowerBound, upperBound, weights, condition
    )


def LessThanPB(clause, val, weights=None, condition=None):
    return PBManager().conditionalPB(clause, val, "<", weights, condition)


def GreaterThanPB(clause, val, weights=None, condition=None):
    return PBManager().conditionalPB(clause, val, ">", weights, condition)


def LessEqPB(clause, val, weights=None, condition=None):
    return PBManager().conditionalPB(clause, val, "<=", weights, condition)


def GreaterEqPB(clause, val, weights=None, condition=None):
    return PBManager().conditionalPB(clause, val, ">=", weights, condition)


def EqualPB(clause, val, weights=None, condition=None):
    return PBManager().conditionalPB(clause, val, "=", weights, condition)


def AssertExactlyOne(clause):
    PBManager().AssertExactlyOne(clause)


def AssertAtMostOne(clause):
    PBManager().AssertAtMostOne(clause)


# shutil.which, backported from the python 3.3 sources, if python version < 3.3
# The 'which' function included here falls under Python's PSF
def which(cmd, mode=os.F_OK | os.X_OK, path=None):
    """Given a command, mode, and a PATH string, return the path which
    conforms to the given mode on the PATH, or None if there is no such
    file.

    `mode` defaults to os.F_OK | os.X_OK. `path` defaults to the result
    of os.environ.get("PATH"), or can be overridden with a custom search
    path.

    """

    # Check that a given file can be accessed with the correct mode.
    # Additionally check that `file` is not a directory, as on Windows
    # directories pass the os.access check.
    def _access_check(fn, mode):
        return os.path.exists(fn) and os.access(fn, mode) and not os.path.isdir(fn)

    # Short circuit. If we're given a full path which matches the mode
    # and it exists, we're done here.
    if _access_check(cmd, mode):
        return cmd

    path = (path or os.environ.get("PATH", os.defpath)).split(os.pathsep)

    if sys.platform == "win32":
        # The current directory takes precedence on Windows.
        if not os.curdir in path:
            path.insert(0, os.curdir)

        # PATHEXT is necessary to check on Windows.
        pathext = os.environ.get("PATHEXT", "").split(os.pathsep)
        # See if the given file matches any of the expected path extensions.
        # This will allow us to short circuit when given "python.exe".
        matches = [cmd for ext in pathext if cmd.lower().endswith(ext.lower())]
        # If it does match, only test that one, otherwise we have to try
        # others.
        files = [cmd] if matches else [cmd + ext.lower() for ext in pathext]
    else:
        # On other platforms you don't have things like PATHEXT to tell you
        # what file suffixes are executable, so just pass on cmd as-is.
        files = [cmd]

    seen = set()
    for dir in path:
        dir = os.path.normcase(dir)
        if not dir in seen:
            seen.add(dir)
            for thefile in files:
                name = os.path.join(dir, thefile)
                if _access_check(name, mode):
                    return name
    return None
