#!/usr/bin/env python3

#The MIT License (MIT)
#
#Copyright (c) 2014, Sam Bayless
#
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
#associated documentation files (the "Software"), to deal in the Software without restriction,
#including without limitation the rights to use, copy, modify, merge, publish, distribute,
#sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all copies or
#substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
#NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
#NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
#DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
#OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

if __name__ == "__main__":
    from os import path
    import os
    import sys
    #Setup PYTHONPATH... if anyone knows a better way to do this without a shell script, I'm all ears..
    sys.path.append(os.path.abspath(os.path.join( path.dirname(__file__),os.pardir)))
    
from monosat.logic import *
from monosat.gnf import writeGNF
from monosat.graphcircuit import Graph, GraphManager
from monosat.pbtheory import PBManager, MinisatPlus
import argparse
import tempfile


    
parser = argparse.ArgumentParser(description='Convert GNF format to plain CNF format')
parser.add_argument('input', type=str,
                   help='Input file, in GNF format')

parser.add_argument('output', type=str,nargs='?',default=None,
                   help='Output file, in CNF format')


args = parser.parse_args()

pbm = PBManager()
pbm.setPB(MinisatPlus())
mgr = GraphManager()

infile = open(args.input,'r')
output = open(args.output,'w') if args.input is not None else sys.stdout

graphs = []

#Read in the input file; echo the cnf right back to stdout - note that this will result in bad headers! Fix that later...

tmpcnf = tempfile.NamedTemporaryFile(mode='w',delete=False)
tmpcnfname = tmpcnf.name
#tmpcnf.close()

nclauses = 0

def getGraph(g):
    global graphs

    return graphs[g]

varmap = dict()
invarmap = dict()
    
maxVar = 0

varmap[truenode.getInputLiteral()] = 0
invarmap[0]=truenode


def getVar(v):
    global varmap
    global invarmap
    global maxVar
    assert(v>=0)    
    if v in invarmap:
        return invarmap[v]
    else:
        v2 = Var()
        varmap[v2.getInputLiteral()]=v
        invarmap[v]=v2
        if(v>maxVar):
            maxVar=v
        return v2

def getLit(l):
    global varmap
    global invarmap
    global maxVar
    v=abs(l)
    
    if v in invarmap:
        v2 = invarmap[v]
    else:
        v2 = Var()
        varmap[v2.getInputLiteral()]=v
        invarmap[v]=v2
        if(v>maxVar):
            maxVar=v
    if(l<0):
        return Not(v2)
    else:
        return v2

used=0

for line in infile:
    lstr = line.lstrip()
    ln = lstr.split()
    if lstr.startswith('c'):
        tmpcnf.write(line)        
    elif lstr.startswith('p graph') or  lstr.startswith('p cnf'):
        #dont do anything
        pass
    elif ln[0] == 'g':
        #f.write("g digraph " + str(self.graphs[gid].nodes) + " " + str(self.graphs[gid].numedges) + " " + str(gid) + "\n") # + " " + str(firstEdge)
        assert(len(ln)==5)
        graph = Graph(mgr)
        nodes = int(ln[2])
        edges = int(ln[3])
        g = int(ln[4])
        for n in range(nodes):
            graph.addNode()
        
        while(len(graphs)<=g):
            graphs.append(None)
        assert(graphs[g] is None)
        graphs[g]=graph
        
        
    elif ln[0] == 'e':
        #f.write("e " + str(gid) + " " + str(v) + " " + str(w) + " " +  str(s.getInputLiteral())  + "\n")
        assert(len(ln)==5)
        g = int(ln[1])
        v = int(ln[2])
        w = int(ln[3])
        var = getVar(int(ln[4]))
        
        getGraph(g).addEdge(v,w,None,var)
    elif ln[0] ==  'w':    
        assert(len(ln)==6)
        g = int(ln[1])
        v = int(ln[2])
        w = int(ln[3])       
        var = getVar(int(ln[4]))
        weight = int(ln[5])
        getGraph(g).addEdge(v,w,weight,var)
    elif ln[0] == 'r' :
        #f.write("r " + str(gid) + " " + str(start) + " " + str(reach) + " " + str(v.getInputLiteral()) + "\n" );
        assert(len(ln)==5)
        g = int(ln[1])
        start = int(ln[2])
        to = int(ln[3])       
        var = getVar(int(ln[4]))
        v = getGraph(g).reaches(start,to)
        Assert(var==v)
         
    elif ln[0]=='d':
        #f.write("r " + str(gid) + " " + str(start) + " " + str(reach) + " " + str(v.getInputLiteral()) + "\n" );
        assert(len(ln)==6)
        g = int(ln[1])
        start = int(ln[2])
        to = int(ln[3])       
        var = getVar(int(ln[4]))
        dist = int(ln[5])
        v = getGraph(g).reaches(start,to,dist)
        Assert(var==v)
    elif ln[0] == 'u':
        assert(len(ln)==5)
        g = int(ln[1])
        start = int(ln[2])
        to = int(ln[3])       
        var = getVar(int(ln[4]))
        v = getGraph(g).connects(start,to)
        Assert(var==v)

    elif ln[0] =='h':
        assert(len(ln)==6)
        g = int(ln[1])
        start = int(ln[2])
        to = int(ln[3])       
        var = getVar(int(ln[4]))
        dist = int(ln[5])
        v = getGraph(g).connects(start,to,dist)
        Assert(var==v)

    elif ln[0] == 'm':
        #f.write("m " + str(gid) + " " + str(minweight) + " " + str(v.getInputLiteral()) + "\n")
        assert(len(ln)==4)
        g = int(ln[1])
        minweight = int(ln[2])
        var = getVar(int(ln[3]))       

        v = getGraph(g).minimumSpanningTreeLessEq(minweight)
        Assert(v==var) 
    elif ln[0] == 'x':
        # f.write("x " + str(gid) + " " + str(edgeVar) + " " + str(var.getInputLiteral()) + "\n")   
        assert(len(ln)==4)
        g = int(ln[1])
        edgeVar =getVar( int(ln[2]))
        var = getVar(int(ln[3]))      
        v = getGraph(g).edgeInMinimumSpanningTree(edgeVar)
        Assert(v==var) 
    elif ln[0]== 'f':
        #f.write("f " + str(gid) + " " + str(s) + " " + str(t) + " " + str(flow)  + " " + str(var.getInputLiteral()) + "\n")
        assert(len(ln)==6)
        g = int(ln[1])
        s = int(ln[2])
        t = int(ln[3])       
        flow = int(ln[4])
        var =  getVar(int(ln[5]))              
        v = getGraph(g).maxFlowGreaterOrEqualToTwoSided(edgeVar)
        Assert(v==var) 
    elif ln[0]== 'n':
        #f.write("n " + str(gid) + " " + str(min_components) + " " + str(v.getInputLiteral()) + "\n")
        assert(len(ln)==4)
        g = int(ln[1])
        min_components = int(ln[2])
        var = getVar(int(ln[3]))       
            
        v = getGraph(g).connectedComponentsGreaterThan(edgeVar)
        Assert(v==var) 
     
    elif ln[0] == 'o':
        #pb.
        size = int(ln[1])
        clause = []
        weights=[]
        for i in range(size):
            clause.append(getLit(int(ln[2+i])))
        wsize = int(ln[2+size])
        for i in range(wsize):
            weights.append(getLit(int(ln[2+size+i +1])))
        pos = 2+size+wsize +1
        op = ln[pos]
        
        rhs = int(ln[pos+1])
        type = int(ln[pos+2])
        headlit=None
        if(type>0):
            headlit = getLit(int(ln[pos+3]))
            
        if(type==0):
            pbm.AssertPB(clause,rhs,op,weights)
   
        elif type==1:
            pbm.conditionalPB(clause,rhs,op,weights,headlit)
    
        elif type==2:

            pbm.twoSidedPB(clause,rhs,op,weights,headlit)
    else:
        nclauses+=1
        #This is a clause
        clause = ln[0:-1]
        assert(len(clause)==len(ln)-1)
        for lit in clause:
            if(abs(int(lit))>=maxVar):
                maxVar = abs(int(lit))
        tmpcnf.write(line)
        
tmpfile = tempfile.NamedTemporaryFile(mode='w',delete=False)
tmpfilename = tmpfile.name
tmpfile.close()


writeGNF( And(getAssertions()),mgr,tmpfilename)

nextVar = maxVar+1

#read back in the gnf, remap the variables as needed 
tmpfile = open(tmpfilename,'r')
for line in tmpfile.readlines():
    if line.startswith("p cnf"): #Header
        header =line.split()
        #vars = int(header[2]) 
    elif line.startswith('c'):
        pass
    elif len(line.split())==0:
        pass
    else:
        nclauses+=1
        clause =list(map(int, line.split()))
        assert(len(clause)>1)
        assert(clause[-1]==0)
        
        for l in clause:
            if l==0:
                
                continue
            v = abs(l)
            if v not in varmap:
                varmap[v]= nextVar
                nextVar+=1
            v = varmap[v]
            if l<0:
                v=-v
            
        #    output.write(str(v) + " ")
        #output.write("0\n")
tmpfile.close()


#ok, now that we know the number of vars and clauses, write everything to output
output.write("p cnf " + str(nextVar) + " " + str(nclauses) + "\n")

#Echo the cnf to the output
tmpcnf.close()
tmpcnf = open( tmpcnfname,'r')
for ln in tmpcnf:
    output.write(ln) 

output.write("c Converted gnf constraints\n")

tmpfile = open(tmpfilename,'r')
for line in tmpfile.readlines():
    if line.startswith("p cnf"): #Header
        header =line.split()
        #vars = int(header[2]) 
    elif line.startswith('c'):
        pass
    elif len(line.split())==0:
        pass
    else:
        nclauses+=1
        clause =list(map(int, line.split()))
        assert(len(clause)>1)
        assert(clause[-1]==0)
        
        for l in clause:
            if l==0:
                
                continue
            v = abs(l)
            assert(v in varmap)
            
            v = varmap[v]
            if l<0:
                v=-v
            output.write(str(v) + " ")
        output.write("0\n")
output.close()  
os.remove(tmpcnfname)
os.remove(tmpfilename)

