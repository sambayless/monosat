#!/usr/bin/env python3
# @author Rémi Delmas remi.delmas.3000@gmail.com
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
import colorama
import sys
from monosat import *

colorama.init()

if len(sys.argv) != 5:
    print(
"""
usage: %s <grid_size> <entry_col> <exit_col> <min_path_length>
- grid_size:       INT
- entry_col:       INT in [0, grid_size-1]
- exit_col:        INT in [0, grid_size-1]
- min_path_lenght: INT in [0, grid_size**2]
""" % sys.argv[0])
    exit(1)

print("Colored Maze generation in MonoSAT with bitvectors & graph theories")

VERBOSE = 0 # verbosity level
GRID_SIZE = int(sys.argv[1])
ENTRY_COL = int(sys.argv[2])
EXIT_COL = int(sys.argv[3])
MIN_PATH_LENGTH = int(sys.argv[4])

assert( 0 <= ENTRY_COL and ENTRY_COL <= GRID_SIZE-1 )
assert( 0 <= EXIT_COL and EXIT_COL <= GRID_SIZE-1 )
assert( 0 <= MIN_PATH_LENGTH and MIN_PATH_LENGTH <= GRID_SIZE**2 ) 

# color wheel RED -> BLUE -> GREEN -> YELLOW
RED = 0 
GREEN = 1
BLUE = 2
YELLOW = 3

def verbose(s, level):
    """
    Prints s if current verbose level is >= to given level
    """
    if ( VERBOSE >= l ):
        print( s )


print( "GRID_SIZE: %d" % GRID_SIZE )

Monosat().newSolver("-no-decide-theories")

# the graph representing the terrain grid
graph = Graph()

# dictionary mapping integer coordinates to node objects (i,j): n
nodes = dict()

# dictionary mapping pairs of node objects  to edge variables (n1,n2): e
edges = dict()

# dictionary associating a color each node n:bv
color = dict() 

def print_color_map(entry_col):
    """
    Prints the current value of the bitvector variable c to stdout as an integer.
    Uses colors to indicate values of terrain material, base or treasure presence.
    :param c:
    :return:
    """
    from colorama import Fore, Back, Style

    for i in range( GRID_SIZE ):

        print( "\ni:%2d | " % i, end="" )

        for j in range( GRID_SIZE ):
            node = nodes[i,j]
            # background color
            color_idx = color[node].value()
            bcolor = Back.RESET
            if color_idx == 0:
                bcolor = Back.RED
            elif color_idx == 1:
                bcolor = Back.GREEN
            elif color_idx == 2:
                bcolor = Back.BLUE
            else : # color_idx == 3 
                bcolor = Back.YELLOW
            if i==0 and j == entry_col:
                print( Fore.BLACK + bcolor + "**" , end="" )
            else:
                print( Fore.BLACK + bcolor + "  " , end="" )

        print( Style.RESET_ALL, end="" )

# generate nodes and their associated decision variables and local constraints
for i in range( GRID_SIZE ):
    for j in range( GRID_SIZE ):
        n = graph.addNode()
        nodes[i, j] = n
        color[n] = BitVector(2) # 4 colors


# generate edges of the graph
for i in range( GRID_SIZE ):
    for j in range( GRID_SIZE ):
        n = nodes[i, j]
        if j < GRID_SIZE-1: # north edge
            other_n = nodes[i, j+1]

            # pointing north
            edge1 = graph.addEdge(n, other_n, 1)
            edges[n, other_n] = edge1

            # pointing south
            edge2 = graph.addEdge(other_n, n, 1)
            edges[other_n,n] = edge2

            # adjacent cells must not have the same color
            Assert( Not( Equal( color[n], color[other_n] ) ) )

            # colors must follow the correct wheel order
            AssertEq( 
                edge1,
                Or(
                    And( color[n] == RED, color[other_n] == GREEN ),
                    And( color[n] == GREEN, color[other_n] == BLUE ),
                    And( color[n] == BLUE, color[other_n] == YELLOW ),
                    And( color[n] == YELLOW, color[other_n] == RED )
                )
            )

            # colors must follow the correct wheel order
            AssertEq( 
                edge2,
                Or(
                    And( color[other_n] == RED, color[n] == GREEN ),
                    And( color[other_n] == GREEN, color[n] == BLUE ),
                    And( color[other_n] == BLUE, color[n] == YELLOW ),
                    And( color[other_n] == YELLOW, color[n] == RED )
                )
            )

        if i < GRID_SIZE-1: # add east edge
            other_n = nodes[i+1, j]

            # pointing east
            edge1 = graph.addEdge( n, other_n, 1 )
            edges[n, other_n] = edge1

            # pointing west
            edge2 = graph.addEdge( other_n, n, 1 )
            edges[other_n,n] = edge2
            
            # adjacent cells must not have the same color 
            Assert( Not( Equal( color[n], color[other_n] ) ) )

            # colors must follow the correct wheel order
            AssertEq( 
                edge1,
                Or(
                    And( color[n] == RED, color[other_n] == GREEN ),
                    And( color[n] == GREEN, color[other_n] == BLUE ),
                    And( color[n] == BLUE, color[other_n] == YELLOW ),
                    And( color[n] == YELLOW, color[other_n] == RED )
                )
            )

            # colors must follow the correct wheel order
            AssertEq( 
                edge2,
                Or(
                    And( color[other_n] == RED, color[n] == GREEN ),
                    And( color[other_n] == GREEN, color[n] == BLUE ),
                    And( color[other_n] == BLUE, color[n] == YELLOW ),
                    And( color[other_n] == YELLOW, color[n] == RED )
                )
            )

        # starting point in the maze
#ENTRY_COL = random.randrange( GRID_SIZE-1 )
print( "Maze entry in col %d on row 0" % ENTRY_COL )
start = nodes[0,ENTRY_COL]
Assert( color[start] == RED )

# exit point in the maze
#EXIT_COL = random.randrange( GRID_SIZE-1 )
end = nodes[GRID_SIZE-1, EXIT_COL]
print( "Maze exit in col %d on row %d" % (EXIT_COL, GRID_SIZE-1) )

# There is a path between selected cells
Assert( graph.reaches( start, end ) )

# It is at least MIN_PATH_LENGTH long
Assert( Not( graph.distance_leq( start, end, MIN_PATH_LENGTH ) ) )

# There is no path between other cells of the first and last row
for col1 in range(GRID_SIZE-1):
    for col2 in range(GRID_SIZE-1):
        if col1 != ENTRY_COL and col2 != EXIT_COL:
            n1 = nodes[0, col1]
            n2 = nodes[GRID_SIZE-1, col2]
            Assert( Not( graph.reaches( n1, n2 ) ) )



result = Solve()
if result:
    print( "\ncolor" )
    print_color_map(ENTRY_COL)
    print( "\n" )
    sys.stdout.flush()
else:
    print("UNSAT")
    sys.stdout.flush()
