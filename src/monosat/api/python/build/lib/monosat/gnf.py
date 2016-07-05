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
from __future__ import division, print_function
from monosat.logic import *
from monosat.bvtheory import BVManager
from monosat.fsmtheory import FSMManager
from monosat.geometrytheory import GeometryManager
from monosat.graphtheory import GraphManager
from monosat.pbtheory import PBManager
from tempfile import NamedTemporaryFile, mktemp
import os
import platform
import sys
import tempfile




def clear():
    clearAIG()
    GraphManager().clear()
    GeometryManager().clear()
    PBManager().clear()    
    FSMManager().clear()
    BVManager().clear()
    

def solve(managers=None,filename=None, ** kwargs):
    witfile=None
    keep_witfile=False
    comments=None
    seed=None
    args=" -decide-theories -no-decide-graph-rnd   -lazy-maxflow-decisions -conflict-min-cut -conflict-min-cut-maxflow -reach-underapprox-cnf"
    
    if kwargs is not None:
        if "comments" in kwargs:
            comments=kwargs["comments"]
        if "seed" in kwargs:
            seed=kwargs["seed"]            
        if "args" in kwargs:
            args=kwargs["args"]     
        if "witness" in kwargs:            
            witfile=kwargs["witness"]
            if witfile is not None:
                keep_witfile=True       


    keepFile = filename is not None
    if not keepFile:
        filename= tempfile.mktemp(".gnf")
    writeGNF(And(getAssertions()),managers,filename,comments)
    
    if witfile is None:
        witfile = tempfile.mktemp()

        
    print("writing witness to " + witfile)
    
    if(seed is None):
        seed=1
    else:
        print("Solver's random seed is %f"%(seed))
        
    ret = os.system("monosat " + args + "  -rnd-seed=%f -theory-witness-file=%s %s"%(seed, witfile,filename))

    #ret = os.system("~/workspaceC/modsat/Debug/monosat -dist=dijkstra -maxflow=edmondskarp-adj -rnd-seed=%f -decide-theories -no-decide-graph-rnd   -lazy-maxflow-decisions -conflict-min-cut -conflict-min-cut-maxflow -reach-underapprox-cnf -theory-witness-file=%s %s"%(seed, witfile,filename))
    #ret = os.system("monosat -no-allow-reach-decision -rnd-seed=2 -decide-theories -reach-underapprox-cnf -conflict-min-cut -witness-file=%s %s"%(witfile,filename))
    sys.stdout.flush()
    sys.stderr.flush()
    if(ret>>8 == 10):
        print("SAT")
        model = loadSolution(witfile)
        if not keepFile:
            os.remove(filename)
        if not keep_witfile:
            os.remove(witfile)
        return model
    if(ret>>8 == 20):
        if not keepFile:
            os.remove(filename)        

        print("UNSAT")
        return False
    print("Error: %d"%(ret) )
    return None

def writeGNF(circuit, managers, filestr,comments=None):
    if not circuit or isTrue(circuit):
        print("Warning: no constraints enforced or written to file (trivially satisfiable)\n")
        return 
    
    if managers is None:
        managers=[]
    if isinstance(managers,tuple):
        managers=list(managers)
    if not isinstance(managers,list):
        managers=[managers]
        
    for m in [GraphManager(),GeometryManager(),FSMManager(),BVManager()]:
        if m not in managers:
            managers.append(m)
    
    print("Writing to " +  filestr + "...\n")
    fileout = open(filestr,'w')
    print("c Created with monosat_py", file=fileout)
    if comments is not None:
        print("c " + str(comments),file=fileout);
    
    
    #aigfile = mktemp(suffix=".aig");
    #print("Writing to " +  aigfile + "...\n")
    if(VAR(circuit).isConstFalse()):
        print("Warning, constraints are trivially false");
    elif(VAR(circuit).isConstTrue()):
        print("Warning, constraints are trivially true");
    #print ("Writing aig to file " + aigfile)


    cnffile = mktemp(suffix=".cnf");
    varmap = WriteCNF(circuit,cnffile)
    setLiteralMap(varmap)    
    
    #append the cnf to the file we are building
    cnfin = open(cnffile,"r");    
    for line in cnfin:
        fileout.write(line);
    #fileout.write(cnfin.read())
    cnfin.close();
    os.remove(cnffile)
    fileout.flush();
    pbmgr = PBManager()
    pbmgr.write(filestr)
    
    print ("Writing graph...")
    graphfile = open(filestr,'a+');
    for manager in managers:
        if(manager is not None):
            manager.write(graphfile)
    all_symbols = list(getSymbols().keys())
    all_symbols.sort()
    for symbol in all_symbols:
        assert(symbol in getSymbols())
        var = getSymbols()[symbol]
        graphfile.write("c var " + str(var.getInputLiteral()) + " " + symbol + "\n")
    
    graphfile.close()
    print("Done writing to " +  filestr + "\n")
    
class Model(dict):
    def __init__(self,*arg,**kw):
        super(Model, self).__init__(*arg, **kw)    
        self.varmap = dict()
        self.inputvarmap = dict()
    
    def __contains__(self, item):
        if isinstance(item,Var):
            if(item.isInput()):
                inlit =  item.getInputLiteral()
                return inlit in self.inputvarmap
            else:
                return item.getLit() in self.varmap
        return super(Model, self).__contains__(item)
    
    def __getitem__(self,item):
        #Storing variables by their literal, because the equality function is over-ridden for variables and causes problems.
        if isinstance(item,Var):
            if(item.isInput()):
                return self.inputvarmap[item.getInputLiteral()]
            else:
                return self.varmap[item.getLit()]
            #
        return super(Model, self).__getitem__(item)
    
    def _setInputValue(self,input_lit, value):
        self.inputvarmap[input_lit]=value;
       
    def __setitem__(self, key, item):
        if isinstance(key,Var):            
            self.varmap[key.getLit()]=item
        else:
            super(Model, self).__setitem__(key,item)

def loadSolution(witfile):
    gmgr = None
    fmgr = None
    bvmgr = None
    managers = [GraphManager(),GeometryManager(),FSMManager(),BVManager()]
    for mgr in managers:
        if isinstance(mgr, GraphManager):
            gmgr = mgr;
        if isinstance(mgr,FSMManager):
            fmgr = mgr
        if isinstance(mgr,BVManager):
            bvmgr = mgr
    
   
    model = Model()    
    try:
        f = open(witfile,'r')            
        for ln in f:          
            if(ln.startswith("var ")):  
                symbol = ln[4:]
                symbol = symbol.strip()
                model[symbol]=True
            elif(ln.startswith("v ")):    
                wit=ln.split()        
                
                if(wit[0] != "v" or wit[-1] != "0"):
                    print("Bad witness!")
                    print(witfile)
                    return None
                
                litmap = getLiteralMap()
                varmap = dict((v,k) for k, v in litmap.items())
                model[truenode]=True
                #model[falsenode]=False                    
                for n in wit[1:-1]:
                    literal = int(n)
                    if literal==35:
                        pass
                   
                    #if (literal in varmap):
                    #    l = varmap[literal]
                        #v =Var(l)
                    model._setInputValue(literal,True)
                        #model[~v]=False
            elif ln.strip().startswith("bv"):
                if bvmgr is not None:
                    line = ln.strip();
                    tokens = line[2:].split()
                    bvID = int(tokens[0])
                    value = int(tokens[2])
                    if bvmgr.hasBitVector(bvID):
                        model[bvmgr.getBitVector(bvID)]=value
            else:
                               
                if gmgr is not None:
                    if  re.match("""Graph \d+ maxflow""",ln):                    
                        flowline= ln.split()
                        
                        graph = int(flowline[1])
                        source = int(flowline[3])
                        sink=int(flowline[5])
                        
                        if(flowline[6]=="is"):
                            #s-t maxflow
                            assert(len(flowline)==8);
                            flow=int(flowline[7])
                            model[(gmgr.getGraph(graph),'maxflow',(source,sink))]=flow              
                        else:          
                            #s-t edge flow assignment      
                            assert(len(flowline)==13);
                            u=int(flowline[8])
                            v = int(flowline[10])
                            flow = int(flowline[12])
                            model[(gmgr.getGraph(graph),'flow',(source,sink),(u,v))]=flow                   
                    elif  re.match("""Graph \d+ Path""",ln):                    
                        pathline= ln.split()
                        
                        graph = int(pathline[1])
                        source = int(pathline[4])
                        dest=int(pathline[6])
                        pathlist =  pathline[9].split(",")                   
                        path = list(map(int,filter(None, pathlist)))         
                        model[(gmgr.getGraph(graph),'reach',(source,dest))]=path              
 
                   
    #except Exception as e:
        #print(e)
        #print("Bad witness")
        #print(witfile)
        
        #return None
    finally:
        pass
    assert(model[truenode]==True)

    return model
