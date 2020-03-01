#!/usr/bin/env python3
# Authors Rémi Delmas & Christophe Garion
#License CC BY-NC-SA 3.0
#This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike
#3.0 Unported license (CC BY-NC-SA 3.0)
#You are free to Share (copy, distribute and transmite) and to Remix (adapt) this work under the following conditions:
#
#Attribution – You must attribute the work in the manner specified by the author or licensor (but not in
#any way that suggests that they endorse you or your use of the work).
#Noncommercial – You may not use this work for commercial purposes.
#Share Alike – If you alter, transform, or build upon this work, you may distribute the resulting work only
#under the same or similar license to this one.
#See http://creativecommons.org/licenses/by-nc-sa/3.0/ for more details.
import sys
import time
from monosat import *

if len(sys.argv) != 2:
    print("You should give the number of queens as parameter!")
    exit(1)

print("N-queens in MonoSAT")

VERBOSE = False

nb_queens = int(sys.argv[1])

if VERBOSE:
    print("Problem with " + str(nb_queens) + " queens")

def print_cell(row, col):
    if QMATRIX[row][col].value():
        print("1 ", end="")
    else:
        print("0 ", end="")


def print_matrix():
    for row in range(nb_queens):
        print()
        for col in range(nb_queens):
            print_cell(row, col)
    print()

time_start_gen = time.process_time()

# a matrix storing variables
QMATRIX = [[Var() for x in range(nb_queens)] for x in range(nb_queens)]

# row clauses
for row in range(nb_queens):
    # build clause for existence of a queen in a row
    c = [QMATRIX[row][col] for col in range(nb_queens)]
    AssertClause(c)
    # if queen is in a column, no other queen in the row
    for col in range(nb_queens):
        for other in range(nb_queens):
            if other != col:
                if VERBOSE:
                    print("~({0},{1}), ~({2},{3})".format(row, col, row, other))
                AssertClause([Not(QMATRIX[row][col]), Not(QMATRIX[row][other])])

# column clauses
for col in range(nb_queens):
    # build clause for existence of a queen in a column
    c = [QMATRIX[row][col] for row in range(nb_queens)]
    AssertClause(c)
    # if queen is in a row, no other queen in the column
    for row in range(nb_queens):
        for other in range(nb_queens):
            if other != row:
                if VERBOSE:
                    print("~({0},{1}), ~({2},{3})".format(row, col, other, col))
                AssertClause([Not(QMATRIX[row][col]), Not(QMATRIX[other][col])])

# diag clauses: setting a queen, compute all exclusion for diags
for row in range(nb_queens):
    for col in range(nb_queens):
        for x in range(1, min(nb_queens - row, nb_queens - col)):
            if VERBOSE:
                print('~({0}, {1}), ~({2}, {3})'.format(row, col, row+x, col+x))
            AssertClause([Not(QMATRIX[row][col]), Not(QMATRIX[row+x][col+x])])
        for x in range(1, min(row, col) + 1):
            if VERBOSE:
                print('~({0}, {1}), ~({2}, {3})'.format(row, col, row-x, col-x))
            AssertClause([Not(QMATRIX[row][col]), Not(QMATRIX[row-x][col-x])])
        for x in range(1, min(row, nb_queens - 1 - col) + 1):
            if VERBOSE:
                print('~({0}, {1}), ~({2}, {3})'.format(row, col, row-x, col+x))
            AssertClause([Not(QMATRIX[row][col]), Not(QMATRIX[row-x][col+x])])
        for x in range(1, min(nb_queens - 1 - row, col) + 1):
            if VERBOSE:
                print('~({0}, {1}), ~({2}, {3})'.format(row, col, row-x, col-x))
            AssertClause([Not(QMATRIX[row][col]), Not(QMATRIX[row+x][col-x])])

time_end_gen = time.process_time()

time_start_solve = time.process_time()

result = Solve()

time_end_solve = time.process_time()

if result:
    print_matrix()
else:
    print("UNSAT!")

print("time needed to generate clauses: " + str(time_end_gen - time_start_gen) + "s")
print("time needed to solve problem: " + str(time_end_solve - time_start_solve) + "s")
