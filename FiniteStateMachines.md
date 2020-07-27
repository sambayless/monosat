# Finite State Machines in MonoSAT
MonoSAT has experimental support for a theory of finite state machines.
This support is less well tested than MonoSAT's graph theory,
and some features are likely to have bugs, so use at your own risk!

MonoSAT's FSM theory is different from and more powerful than
string theory solvers (such as on [z3str][z3str]).
z3str supports constraints with strings variables, and several
common string operations (such as contains, startswith, etc),
and also supports a predicate that determines whether a string
is accepted by a regular expression. Regular expressions in z3str
must be constant (ground terms) in the formula - they cannot contain
variables.

MonoSAT doesn't directly support a notion of a regular expressions
or other string operations (though it might in the future).
Instead, MonoSAT supports constraints over strings and non-deterministic
finite state machines (which can encode regular expressions).
In MonoSAT's finite state machine theory, the transitions of the finite
state machines may (optionally) be variables, constrained by Boolean
constrains in the same way as edges in MonoSAT's graph theory.
These constraints can be used to synthesize (bounded-size) finite state
machines that accept or reject certain sets of strings.

For example, here is how to define a finite state machine with two symbolic
transitions (using GNF format; you can also create FSM constraitns
programatically - see the C/Python APIs):

```
c FSM with ID '1', that an input language of 26 characters:
fsm 1 26 0 0
c
c Transition for FSM '1', from state 2 to state 3, accepting character
c 22, bound to Boolean variable '5'
transition 1 2 3 22 0 5
c
c Transition for FSM '1', from state 0 to state 4, accepting character
c 23, bound to Boolean variable '6'
transition 1 0 4 23 0 6
c
c Constraint asserting that at least one of these transitions must be
c assigned false (eg, disabled)
-5 -6 0
```

The formula may define multiple finite state machines - each must be
assigned a unique non-negative integer ID (the above FSM has id '1', but
'0' is also a valid ID).
Finite state machines may be deterministic or non-deterministic and may
have loops;
the special character '0' is reserved and denotes a non-comsuming
'epsilon'-transition. As in MonoSAT's graph theory, each FSM has a
fixed, finite number of states (typically, < 100), and only transitions
that are defined in the formula and are assigned 'true' by the solver
may be used. The states in the FSM do not need to be explicitly defined;
they are inferred from the transitions that are defined.

You may also define string constants, and assert that a given string
must, or must not, be accepted by a given finite state machine (under
assignment to its transitions). It is also possible to define variable
strings; I'll describe that a bit lower down.

```
c The constant string "1234", with ID '1'.
c String characters must be positive integers.
c 0 terminates the string, but is not an element of the string
c (and is not a valid character)
str 1 1 2 3 4 0
c
c A second string, with ID '2', containing only character '23':
str 2 23 0
c
c An empty string, with ID '3':
str 3 0
c
c A string acceptance predicate, bound to variable '7', that evaluates
c to 'true' iff FSM 1, starting at state 0, accepts string '2',
c where state '4' of the finite state machine is treated as the accepting
c state.
accepts 1 0 4 2 7
c
c Constraint enforcing that the FSM must accept the string.
7 0
c Constraint enforcing that string 1 is not accepted:
accepts 1 0 4 1 8
-8 0
```

Each string constant must have a unique, non-negative ID, and must be a fixed
length. Characters in the string are positive integers (sot they cannot
be '0').

A valid solution to the above formula assigns transition '6' to
true and tranisiton '7' to false.

MonoSAT does not directly support a notion of string variables.
Instead, MonoSAT supports a notion of finite state machine generators,
that produce a string rather than consuming a string.
(MonoSAT also supports a more general notion of finite string
transducers, that both consume and generate a string, but that support
is currently only partially functional).
A linear, acyclic finite state generator may be used to encode string
variables
(specifically, the finite state generator must have no loops, no
epsilon transitions, and if there are any transitions from state 's',
they must all go to the same state.)


```
c A finite state machine generator, with ID 2.
c The output language of this FSM has 26 characters.
fsm 2 0 26 0
c Transition for FSM '2', from state 1 to state 5, producing character
c 23, bound to Boolean variable '9'
transition 1 1 5 0 23 9
c
c A second transition for FSM '2', also from state 1 to state 5,
c but producing character 17 (bound to variable '10')
transition 1 1 5 0 17 10
c If transition '9' is not enabled, then transition '10' also must not be enabled.
9 -10 0
c
c Constraint asserting that FSM 1 must accept the string produced by
c FSM 2, bound to variable 11. The must generator start at state 1 and
c end at state 5; the acceptor must start at state 0 and finish at
c state 4:
accepts_composition 2 1 1 5 0 4 0 11
c
c FSM 1 must accept the string produced by FSM 2:
11 0
```

### FSM Syntax Guide

```
fsm FSM_ID IN_LANGUAGE_SIZE OUT_LANGUAGE_SIZE 0 // last value is not used
transition FSM_ID FROM_STATE TO_STATE CONSUME_CHAR PRODUCE_CHAR VAR
str STR_ID [0 or more positive integers] 0
accepts FSM_ID START_STATE END_STATE STR_ID VAR
accepts_composition GEN_FSM_ID ACC_FSM_ID GEN_START GEN_END ACC_START ACC_END 0 VAR
```

Currently, the accepts_composition method is only implemented for linear
generator FSMs (which are powerful enough to encode string variables,
if combined with a pseudoBoolean constraint enforcing that at most one
transition between each state of the generator may be used - you may use
MonoSAT's built-in pseudoBoolean constraints to do so):

```
c At most 1 of the 2 literals in the set {9,10} may be true.
pb <= 1 2 9 10 0
c At least one of the 2 literals in the set {9, 10} must be true
9 10 0
```


# References
[z3str]: https://github.com/z3str/Z3-str, https://ece.uwaterloo.ca/~vganesh/Publications_files/vg2017-Z3str2-FMSD2017-Journal.pdf
