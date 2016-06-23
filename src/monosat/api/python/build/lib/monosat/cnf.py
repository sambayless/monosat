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

'''
Functions to write theory constraints to pure CNF, for use in SAT solvers 
'''
from math import *
from monosat.aiger import *
from monosat.graphcircuit import CNFGraph
from monosat.pbtheory import PBManager
from tempfile import NamedTemporaryFile, mktemp
import monosat.graphtheory
import os
import platform
import sys




def writeCNF(circuit, graphmgr, filestr,comments=None, use_unary_encoding=True):    
    print("Writing to " +  filestr + "...\n")
    fileout = open(filestr,'w')
    print("c Created with monosat_py", file=fileout)
    if comments is not None:
        print("c " + str(comments),file=fileout);
    fileout.close()
    print("Writing to " +  filestr + "...\n")
    
    if(VAR(circuit).isConstFalse()):
        print("Warning, constraints are trivially false");
    elif(VAR(circuit).isConstTrue()):
        print("Warning, constraints are trivially true");

    print("Building CNF graph constraints:")
    old_assertions = list(assertions)    
    clearAssertions()
    Assert(circuit)
    for gid in range(len(graphmgr.graphs)):
        
        g = CNFGraph(gid,use_unary_encoding)
        for _ in range(graphmgr.graphs[gid].numNodes()):
            g.addNode()
        
        for (v,w,var,weight) in graphmgr.graphs[gid].alledges:
            g.addEdge(v,w,weight,var)
                
        for (start,steps,reach,v) in graphmgr.graphs[gid].queries:
            Assert(v== g.reaches(start,reach,steps))
                
                
        for (start,steps,reach,v) in graphmgr.graphs[gid].undirected_queries:
            Assert(v== g.connects(start,reach,steps))
        
        for (minweight,v) in graphmgr.graphs[gid].mstQueries:
            if v is truenode:
                g.AssertMinimumSpanningTreeLessEq(minweight)
            else:
                Assert(v== g.minimumSpanningTreeLessEq(minweight))
                   
        for (edgeVar,v) in graphmgr.graphs[gid].mstEdgeQueries:     
           
            Assert(v== g.edgeInMinimumSpanningTree(edgeVar))  
            
        
        for (s,t, flow,v) in graphmgr.graphs[gid].flowQueries:
            #catch constants here 
            if v is truenode:
                g.AssertMaxFlowGreaterOrEqualToOneSided(s,t,flow)
            elif v is falsenode:
                g.AssertMinCutLessOrEqualToOneSided(s,t,flow-1)
            else:
                Assert(v== g.maxFlowGreaterOrEqualTo(s,t,flow))
            
        for (min_components,v) in graphmgr.graphs[gid].componentsQueries :
            if v is truenode:
                g.AssertConnectedComponentsGreaterThan(min_components)
            elif v is falsenode:
                g.AssertConnectedComponentsLessEq(min_components)
            else:
                Assert(v== g.connectedComponentsGreaterThan(min_components))                   
      
        for steiner in graphmgr.graphs[gid].steiners:
            #Not implemented.
            assert(False)
            
     
    
    
    varmap = WriteCNF(And(assertions),filestr)
    setLiteralMap(varmap)
    
    pbmgr = PBManager()
    pbmgr.write(filestr)    
    fileout = open(filestr,'a')
    all_symbols = list(getSymbols().keys())
    all_symbols.sort()
    for symbol in all_symbols:
        assert(symbol in getSymbols())
        var = getSymbols()[symbol]
        fileout.write("c var " + str(var.getInputLiteral()) + " " + symbol + "\n")
    fileout.close()
    print("Done writing to " +  filestr + "\n")   
    
    clearAssertions()
    setAssertions(old_assertions)
    clearLiteralMap()
     