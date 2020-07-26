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
from monosat.logic import *
from monosat.manager import Manager

debug = False


# Collects a set of finite state machines and strings to encode together into a formula
class FSMManager(metaclass=Manager):

    def __init__(self):
        self.strings = []
        self.fsms = []

    # Each "string" must actually be a list of positive integers
    def newString(self, string):
        ints = []
        for s in string:
            if (isinstance(s, int)):
                ints.append(s + 1)
            else:
                ints.append(ord(s) - ord('a') + 1)

        return monosat.monosat_c.Monosat().newString(ints)


class FSM():
    def __init__(self, in_labels, out_labels=0, hasEpsilon=True):
        self._monosat = monosat.monosat_c.Monosat()
        self.pid = self._monosat.newFSM(in_labels, out_labels)
        self.hasEpsilon = hasEpsilon
        self.n_states = 0
        self._transitions = []

    def getID(self):
        return self.pid

    def addState(self):
        self.n_states += 1
        return self._monosat.newState(self.getID())

    def nStates(self):
        return self.n_states;

    def addTransition(self, u, v, in_label, out_label=0):
        assert (u < self.n_states)
        assert (u >= 0)
        assert (v < self.n_states)
        assert (v >= 0)
        var = Var(self._monosat.newTransition(self.getID(), u, v, in_label, out_label))
        self._transitions.append((u, v, in_label, out_label, var))
        return var

    def getTransitions(self):
        return self._transitions

    def accepts(self, starting_state, accepting_state, stringID):
        return Var(self._monosat.fsmAcceptsString(self.getID(), starting_state, accepting_state, stringID))

    def acceptsFSM(self, generator, generator_source, generator_dest, starting_state, accepting_state):
        return Var(self._monosat.fsmCompositionAccepts(generator.getID(), self.getID(), generator_source, generator_dest,
                                                       starting_state, accepting_state, -1))

    def draw(self):
        print("digraph{")
        for (u, v, label, output, var) in self._transitions:
            if var.value() is not False:
                print("n%d -> n%d [label=\"%d:%d\"]" % (u, v, label, output))

        print("}")
