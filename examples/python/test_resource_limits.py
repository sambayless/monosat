import os
import sys

from monosat import *

import math
import random


print("begin encode");
width=21
height=21
min_steps=221
max_steps=width*height
random.seed(1)
chromatic_graph=Graph();

class Block:
    def __init__(self,graph,x,y):
        self.color =Color(x,y)
        self.node = graph.addNode()

class Color:
    
    def __init__(self, x,y):
        
        self.red=Var("cell(red,"+str(x)+"," + str(y)+")");
        self.yellow=Var("cell(yellow,"+str(x)+"," + str(y)+")");
        self.green=Var("cell(green,"+str(x)+"," + str(y)+")");
        self.cyan=Var("cell(cyan,"+str(x)+"," + str(y)+")");
        self.blue=Var("cell(blue,"+str(x)+"," + str(y)+")");
        self.magenta=Var("cell(magenta,"+str(x)+"," + str(y)+")");
        
        self.colorbits=(self.red,self.yellow,self.green,self.cyan,self.blue,self.magenta)
        self.colornames = ("r","y","g","c","b","m")
        
        
startNode=chromatic_graph.addNode()
exitNode=chromatic_graph.addNode()

room=[]
for x in range(0,width):
    room.append([])
for y in range(0,height): 
    for x in range (0,width):
        room[x].append(Block(chromatic_graph, x,y))


#Exactly 1 color selected for each block
for y in range(0,height): 
    for x in range (0,width):
        AssertEqualPB(room[x][y].color.colorbits,1)
        


def passable(color1,color2):
    same = And(color1.red,color2.red)
    same = Or(same,And(color1.yellow,color2.yellow)) 
    same = Or(same,And(color1.green,color2.green)) 
    same = Or(same,And(color1.cyan,color2.cyan)) 
    same = Or(same,And(color1.blue,color2.blue)) 
    same = Or(same,And(color1.magenta,color2.magenta)) 
    
    nextColor = And(color1.red,color2.yellow)
    nextColor = Or(nextColor, And(color1.yellow,color2.green))
    nextColor = Or(nextColor, And(color1.green,color2.cyan))
    nextColor = Or(nextColor, And(color1.cyan,color2.blue))
    nextColor = Or(nextColor, And(color1.blue,color2.magenta))
    nextColor = Or(nextColor, And(color1.magenta,color2.red))
    
    nextColor = Or(nextColor, And(color1.yellow,color2.red))
    nextColor = Or(nextColor, And(color1.green,color2.yellow))
    nextColor = Or(nextColor, And(color1.cyan,color2.green))
    nextColor = Or(nextColor, And(color1.blue,color2.cyan))
    nextColor = Or(nextColor, And(color1.magenta,color2.blue))
    nextColor = Or(nextColor, And(color1.red,color2.magenta))
    
    return Or(same,nextColor)
    
#Add edges to the graph
for x in range(width):
    for y in range(height):
        if(x<width-1):
            e =chromatic_graph.addEdge( room[x][y].node,room[x+1][y].node)
            #this edge is enabled IFF you can transition between the colors in this room
            Assert(e==passable(room[x][y].color,room[x+1][y].color))
        if(y<height-1):
            e =chromatic_graph.addEdge( room[x][y].node,room[x][y+1].node)
            #this edge is enabled IFF you can transition between the colors in this room
            Assert(e==passable(room[x][y].color,room[x][y+1].color))
        if(x>0):
            e =chromatic_graph.addEdge( room[x][y].node,room[x-1][y].node)
            #this edge is enabled IFF you can transition between the colors in this room
            Assert(e==passable(room[x][y].color,room[x-1][y].color))
        if(y>0):
            e =chromatic_graph.addEdge( room[x][y].node,room[x][y-1].node)
            #this edge is enabled IFF you can transition between the colors in this room
            Assert(e==passable(room[x][y].color,room[x][y-1].color))

#connect start and exit to every node

entranceEdges =[]
exitEdges = []
for x in range(width):
    for y in range(height):
        entranceEdges.append( chromatic_graph.addEdge(startNode, room[x][y].node))        
        exitEdges.append(chromatic_graph.addEdge( room[x][y].node,exitNode))



#Add two, to account for the entrance and exit nodes
Assert(Not(chromatic_graph.distance_lt(startNode, exitNode, min_steps+2))) #The exit must not be reachable in less than min_steps (+2) steps

Assert( chromatic_graph.distance_leq(startNode, exitNode, max_steps+2))
#The second node in the reaches array is the exit node, because it is the second node that was created.

AssertEqualPB (entranceEdges,1)
AssertEqualPB(exitEdges,1)
try:
    r = Solve(time_limit_seconds=5)
    if r:
        for y in range(height):
            for x in range(width):
                count=0
                for color,name in zip( room[x][y].color.colorbits,room[x][y].color.colornames) :
                    if color.value():
                        print(name,end=" ")
                        count+=1
                assert(count==1)    
            print("")
    else:
        print("UNSAT")
except RuntimeError as e:
    print("Terminated early (as expected)")
print("Python Script Done!\n")