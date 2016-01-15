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
from __future__ import division
from __future__ import print_function
from ctypes import *
from os import path

from monosat.singleton import Singleton

import os
#module_path = os.path.abspath(path.dirname(__file__))
try:
    _monosat_c= cdll.LoadLibrary("libmonosat.so")
except Exception as e:
    print (e)

    _monosat_c= cdll.LoadLibrary('monosat.dll')

#Python interface to MonoSAT
    
#libc = CDLL('libc.so.6')
null_ptr = POINTER(c_int)()

c_uint_p = POINTER(c_int)
c_int_p = POINTER(c_int)

c_solver_p = c_void_p
c_graph_p = c_void_p
c_bv_p = c_void_p

c_literal = c_int
c_literal_p = c_int_p
c_bvID=c_int
c_bvID_p = c_int_p
c_var = c_int
c_var_p = c_int_p

def dimacs(l):
    assert(isinstance(l,int))
    if l & 1:
        return -(l//2 +1)
    else:
        return l//2+1
        
class Solver():
    def __init__(self,monosat_c,arguments=None):
        self.elapsed_time=0
        if arguments is None:
            arguments=[]    
        elif isinstance(arguments,str):
            #split this up by whitespace
            arguments=arguments.split()  
        self.arguments=arguments 
        arguments=["monosat"]+arguments
        arr = (c_char_p*len(arguments))()     
        for i,arg in enumerate(arguments):
            arr[i]=c_char_p(bytes(arg, 'utf-8')) 

        #In the future, allow multiple solvers to be instantiated...
        self.comments=[]
        
        self._ptr = monosat_c.newSolver_args(len(arguments),arr)
        self.bvtheory=  monosat_c.initBVTheory(self._ptr)
        self.output=None
        self.arguments=None
        self.symbolmap=dict()  
        self.graphs = []
        self.graph_ids=dict()
        self._true = None
      
    def delete(self):        
        Monosat().monosat_c.deleteSolver(self._ptr)
        self._ptr=None  
        if Monosat().solver is self:
            Monosat().solver=None

#very simple, low-level python to Monosat's C interface.
class Monosat(metaclass=Singleton):       
    def __init__(self):
        self._managers=dict()
        self.monosat_c=_monosat_c
        self.elapsed_time=0
        self._int_array = (c_int * (1024))()
        
        #Set the return types for each function
        self.monosat_c.newSolver.argtypes=[]
        self.monosat_c.newSolver.restype=c_solver_p
        
        self.monosat_c.newSolver_arg.argtypes=[c_char_p]
        self.monosat_c.newSolver_arg.restype=c_solver_p
        
        self.monosat_c.newSolver_args.argtypes=[c_int, POINTER(c_char_p)]
        self.monosat_c.newSolver_args.restype=c_solver_p
        
        self.monosat_c.deleteSolver.argtypes=[c_solver_p]
        
        self.monosat_c.readGNF.argtypes=[c_solver_p, c_char_p]
        
        self.monosat_c.solve.argtypes=[c_solver_p]
        self.monosat_c.solve.restype=c_bool
        
        self.monosat_c.solveAssumptions.argtypes=[c_solver_p,c_literal_p,c_int]
        self.monosat_c.solveAssumptions.restype=c_bool       
        
        self.monosat_c.solveAssumptions_MinBVs.argtypes=[c_solver_p,c_literal_p,c_int,c_int_p,c_int]
        self.monosat_c.solveAssumptions_MinBVs.restype=c_bool       

        self.monosat_c.solveLimited.argtypes=[c_solver_p]
        self.monosat_c.solveLimited.restype=c_int       

        self.monosat_c.solveAssumptionsLimited.argtypes=[c_solver_p,c_literal_p,c_int]
        self.monosat_c.solveAssumptionsLimited.restype=c_int      

        self.monosat_c.solveAssumptionsLimited_MinBVs.argtypes=[c_solver_p,c_literal_p,c_int,c_int_p,c_int]
        self.monosat_c.solveAssumptionsLimited_MinBVs.restype=c_int  
        
        
        self.monosat_c.setTimeLimit.argtypes=[c_solver_p,c_int]
        self.monosat_c.setMemoryLimit.argtypes=[c_solver_p,c_int]
        self.monosat_c.setConflictLimit.argtypes=[c_solver_p,c_int]
        self.monosat_c.setPropagationLimit.argtypes=[c_solver_p,c_int]    
        
        self.monosat_c.backtrack.argtypes=[c_solver_p]
        

        self.monosat_c.newVar.argtypes=[c_solver_p]
        self.monosat_c.newVar.restype=c_int

        self.monosat_c.setDecisionVar.argtypes=[c_solver_p,c_var,c_bool]
    
        self.monosat_c.isDecisionVar.argtypes=[c_solver_p,c_var]
        self.monosat_c.isDecisionVar.restype=c_bool
        
        self.monosat_c.setDecisionPriority.argtypes=[c_solver_p,c_var,c_int]
   
        self.monosat_c.getDecisionPriority.argtypes=[c_solver_p,c_var]
        self.monosat_c.getDecisionPriority.restype=c_int

        self.monosat_c.setDecisionPolarity.argtypes=[c_solver_p,c_var,c_bool]
   
        self.monosat_c.getDecisionPolarity.argtypes=[c_solver_p,c_var]
        self.monosat_c.getDecisionPolarity.restype=c_bool
 
        self.monosat_c.disallowLiteralSimplification.argtypes=[c_solver_p,c_literal]
        self.monosat_c.disallowLiteralSimplification.restype=c_bool
        
        self.monosat_c.addClause.argtypes=[c_solver_p,c_literal_p,c_int]
        self.monosat_c.addClause.restype=c_bool
 
        self.monosat_c.addUnitClause.argtypes=[c_solver_p,c_literal]
        self.monosat_c.addUnitClause.restype=c_bool

        self.monosat_c.addBinaryClause.argtypes=[c_solver_p,c_literal,c_literal]
        self.monosat_c.addBinaryClause.restype=c_bool

        self.monosat_c.addTertiaryClause.argtypes=[c_solver_p,c_literal,c_literal,c_literal]
        self.monosat_c.addTertiaryClause.restype=c_bool
        
        self.monosat_c.true_lit.argtypes=[c_solver_p]
        self.monosat_c.true_lit.restype=c_int 
        
        self.monosat_c.at_most_one.argtypes=[c_solver_p,c_var_p,c_int]

        
        self.monosat_c.initBVTheory.argtypes=[c_solver_p];
        self.monosat_c.initBVTheory.restype=c_bv_p;
        
        self.monosat_c.newBitvector_const.argtypes=[c_solver_p,c_bv_p, c_int, c_long]
        self.monosat_c.newBitvector_const.restype=c_bvID
        
        self.monosat_c.newBitvector.argtypes=[c_solver_p,c_bv_p, c_var_p, c_int]
        self.monosat_c.newBitvector.restype=c_bvID
        
        self.monosat_c.nBitvectors.argtypes=[c_solver_p,c_bv_p]
        self.monosat_c.nBitvectors.restype=c_int
        self.monosat_c.newBVComparison_const_lt.argtypes=[c_solver_p,c_bv_p,c_bvID, c_long]
        self.monosat_c.newBVComparison_const_lt.restype=c_literal
        
        self.monosat_c.newBVComparison_bv_lt.argtypes=[c_solver_p,c_bv_p,c_bvID, c_bvID]
        self.monosat_c.newBVComparison_bv_lt.restype=c_literal
        
        self.monosat_c.newBVComparison_bv_leq.argtypes=[c_solver_p,c_bv_p,c_bvID, c_bvID]
        self.monosat_c.newBVComparison_bv_leq.restype=c_literal
        
        self.monosat_c.newBVComparison_bv_gt.argtypes=[c_solver_p,c_bv_p,c_bvID, c_bvID]
        self.monosat_c.newBVComparison_bv_gt.restype=c_literal
        
        self.monosat_c.newBVComparison_bv_geq.argtypes=[c_solver_p,c_bv_p,c_bvID, c_bvID]
        self.monosat_c.newBVComparison_bv_geq.restype=c_literal
        
        self.monosat_c.bv_addition.argtypes=[c_solver_p,c_bv_p,c_bvID, c_bvID, c_bvID]
        self.monosat_c.bv_subtraction.argtypes=[c_solver_p,c_bv_p,c_bvID, c_bvID, c_bvID]
        
        self.monosat_c.bv_ite.argtypes=[c_solver_p,c_bv_p,c_literal, c_bvID,c_bvID,c_bvID]

        self.monosat_c.bv_min.argtypes=[c_solver_p,c_bv_p, c_bvID_p, c_int,c_bvID]
        self.monosat_c.bv_max.argtypes=[c_solver_p,c_bv_p, c_bvID_p, c_int, c_bvID]
        

        self.monosat_c.bv_not.argtypes=[c_solver_p,c_bv_p,c_bvID,c_bvID]
        self.monosat_c.bv_and.argtypes=[c_solver_p,c_bv_p,c_bvID,c_bvID,c_bvID]
        self.monosat_c.bv_nand.argtypes=[c_solver_p,c_bv_p,c_bvID,c_bvID,c_bvID]    
        self.monosat_c.bv_or.argtypes=[c_solver_p,c_bv_p,c_bvID,c_bvID,c_bvID]   
        self.monosat_c.bv_nor.argtypes=[c_solver_p,c_bv_p,c_bvID,c_bvID,c_bvID]        
        self.monosat_c.bv_xor.argtypes=[c_solver_p,c_bv_p,c_bvID,c_bvID,c_bvID]        
        self.monosat_c.bv_xnor.argtypes=[c_solver_p,c_bv_p,c_bvID,c_bvID,c_bvID]

        self.monosat_c.bv_concat.argtypes=[c_solver_p,c_bv_p,c_bvID,c_bvID,c_bvID]
        self.monosat_c.bv_popcount.argtypes=[c_solver_p,c_bv_p,c_literal_p,c_int ,c_bvID]
    
        
        self.monosat_c.bv_slice.argtypes=[c_solver_p,c_bv_p,c_bvID,c_int,c_int,c_bvID]

        self.monosat_c.newGraph.argtypes=[c_solver_p]
        self.monosat_c.newGraph.restype=c_graph_p

        self.monosat_c.newNode.argtypes=[c_solver_p,c_graph_p]
        self.monosat_c.newNode.restype=c_int
        
        self.monosat_c.newEdge.argtypes=[c_solver_p,c_graph_p, c_int, c_int, c_long]
        self.monosat_c.newEdge.restype=c_literal

        
        self.monosat_c.newEdge_double.argtypes=[c_solver_p,c_graph_p, c_int, c_int, c_double]
        self.monosat_c.newEdge_double.restype=c_literal
       
        self.monosat_c.newEdge_bv.argtypes=[c_solver_p,c_graph_p, c_int, c_int, c_bvID]
        self.monosat_c.newEdge_bv.restype=c_literal
       
        self.monosat_c.reaches.argtypes=[c_solver_p,c_graph_p, c_int, c_int]
        self.monosat_c.reaches.restype=c_literal

        self.monosat_c.shortestPathUnweighted_lt_const.argtypes=[c_solver_p,c_graph_p, c_int, c_int,c_int]
        self.monosat_c.shortestPathUnweighted_lt_const.restype=c_literal

        self.monosat_c.shortestPathUnweighted_leq_const.argtypes=[c_solver_p,c_graph_p, c_int, c_int,c_int]
        self.monosat_c.shortestPathUnweighted_leq_const.restype=c_literal


        self.monosat_c.shortestPath_lt_const.argtypes=[c_solver_p,c_graph_p, c_int, c_int,c_long]
        self.monosat_c.shortestPath_lt_const.restype=c_literal

        self.monosat_c.shortestPath_leq_const.argtypes=[c_solver_p,c_graph_p, c_int, c_int,c_long]
        self.monosat_c.shortestPath_leq_const.restype=c_literal
              

        self.monosat_c.shortestPath_lt_bv.argtypes=[c_solver_p,c_graph_p, c_int, c_int,c_bvID]
        self.monosat_c.shortestPath_lt_bv.restype=c_literal

        self.monosat_c.shortestPath_leq_bv.argtypes=[c_solver_p,c_graph_p, c_int, c_int,c_bvID]
        self.monosat_c.shortestPath_leq_bv.restype=c_literal        
        
        self.monosat_c.nVars.argtypes=[c_solver_p]
        self.monosat_c.nVars.restype=c_int      

        self.monosat_c.nClauses.argtypes=[c_solver_p]
        self.monosat_c.nClauses.restype=c_int   

        self.monosat_c.maximumFlow_geq.argtypes=[c_solver_p,c_graph_p, c_int, c_int,c_long]
        self.monosat_c.maximumFlow_geq.restype=c_literal               
        
        self.monosat_c.maximumFlow_gt.argtypes=[c_solver_p,c_graph_p, c_int, c_int,c_long]
        self.monosat_c.maximumFlow_gt.restype=c_literal      


        self.monosat_c.maximumFlow_geq_bv.argtypes=[c_solver_p,c_graph_p, c_int, c_int,c_bvID]
        self.monosat_c.maximumFlow_geq_bv.restype=c_literal   
        
        self.monosat_c.maximumFlow_gt_bv.argtypes=[c_solver_p,c_graph_p, c_int, c_int,c_bvID]
        self.monosat_c.maximumFlow_gt_bv.restype=c_literal           
        
        self.monosat_c.minimumSpanningTree_leq.argtypes=[c_solver_p,c_graph_p, c_long]
        self.monosat_c.minimumSpanningTree_leq.restype=c_literal    
             
        self.monosat_c.minimumSpanningTree_lt.argtypes=[c_solver_p,c_graph_p, c_long]
        self.monosat_c.minimumSpanningTree_lt.restype=c_literal            

        self.monosat_c.acyclic_undirected.argtypes=[c_solver_p,c_graph_p]
        self.monosat_c.acyclic_undirected.restype=c_literal      
        
        self.monosat_c.acyclic_directed.argtypes=[c_solver_p,c_graph_p]
        self.monosat_c.acyclic_directed.restype=c_literal            
        
        self.monosat_c.getModel_Literal.argtypes=[c_solver_p,c_literal]
        self.monosat_c.getModel_Literal.restype=c_int      

        self.monosat_c.getModel_BV.argtypes=[c_solver_p,c_bv_p, c_bvID, c_bool]
        self.monosat_c.getModel_BV.restype=c_long    

        self.monosat_c.getModel_MaxFlow.argtypes=[c_solver_p,c_graph_p, c_literal]
        self.monosat_c.getModel_MaxFlow.restype=c_long            

        self.monosat_c.getModel_EdgeFlow.argtypes=[c_solver_p,c_graph_p, c_literal, c_literal]
        self.monosat_c.getModel_EdgeFlow.restype=c_long        

        self.monosat_c.getModel_AcyclicEdgeFlow.argtypes=[c_solver_p,c_graph_p, c_literal, c_literal]
        self.monosat_c.getModel_AcyclicEdgeFlow.restype=c_long     

        self.monosat_c.getModel_MinimumSpanningTreeWeight.argtypes=[c_solver_p,c_graph_p, c_literal]
        self.monosat_c.getModel_MinimumSpanningTreeWeight.restype=c_long     
        self.newSolver()
        #For many (but not all) instances, the following settings may give good performance: 
        #self.init("-verb=0 -verb-time=0 -rnd-theory-freq=0.99 -no-decide-bv-intrinsic  -decide-bv-bitwise  -decide-graph-bv -decide-theories -no-decide-graph-rnd   -lazy-maxflow-decisions -conflict-min-cut -conflict-min-cut-maxflow -reach-underapprox-cnf -check-solution ")

    def getGID(self,graph):
        return self.solver.graph_ids[graph]
    
    
    def _getManagers(self):
        solver = self.getSolver()
        if solver not in self._managers:
            self._managers[solver]=dict()
        return self._managers[solver]
    
    def getSolver(self):
        return self.solver
    
    def setSolver(self,solver):
        self.solver=solver
    
    #Convenience method to re-initialize monosat with new solver
    def init(self,arguments=None):
        if self.solver is not None:
            self.solver.delete()        
        return self.newSolver(arguments)
    
    #Until a better system is created, this can be used to re-initialize Monosat with a new configuration.    
    def newSolver(self,arguments=None):
        self.solver = Solver(self.monosat_c,arguments)
        self.solver._true = self.getTrue()
        return self.solver


    def readGNF(self, filename):
        self.monosat_c.readGNF(self.solver._ptr,c_char_p(filename.encode('ascii')))
    
    def getOutputFile(self):
        return self.solver.output
    
    def setOutputFile(self, file):
        self.solver.output=file
        if file is not None:#Otherwise, this unit clause will have been skipped in the output file
            if self.solver.arguments is not None:
                self._echoOutput("c monosat " + " ".join(self.solver.arguments) +"\n")           
            self.addUnitClause(self.solver._true) 
    
    def _echoOutput(self,line):
        if self.solver.output:
            self.solver.output.write(line)
    
    def comment(self,c):
        self.solver.has_comments=True
        self.solver.comments.append(c)
        if self.solver.output:
            self._echoOutput("c " + c +"\n")           
    
    def getIntArray(self,nums):
        if len(nums)>len(self._int_array):
            self._int_array = (c_int * len(nums))()
        for i,n in enumerate(nums):            
            self._int_array[i]=c_int(n)
        return self._int_array
        

    def setTimeLimit(self,seconds):
        if seconds is None or seconds<0:
            self.monosat_c.setTimeLimit(self.solver._ptr,-1)
        else:
            self.monosat_c.setTimeLimit(self.solver._ptr,seconds)
            
    def setMemoryLimit(self,mb):
        if mb is None or mb<0:
            self.monosat_c.setMemoryLimit(self.solver._ptr,-1)
        else:
            self.monosat_c.setMemoryLimit(self.solver._ptr,mb)

    def setConflictLimit(self,conflicts):
        if conflicts is None or conflicts<0:
            self.monosat_c.setConflictLimit(self.solver._ptr,-1)
        else:
            self.monosat_c.setConflictLimit(self.solver._ptr,conflicts)

    def setPropagationLimit(self,propagations):
        if propagations is None or propagations<0:
            self.monosat_c.setPropagationLimit(self.solver._ptr,-1)
        else:
            self.monosat_c.setPropagationLimit(self.solver._ptr,propagations)



    def solve(self,assumptions=None,minimize_bvs=None):
        self.backtrack()
        if assumptions is None:
            assumptions=[]

        lp = self.getIntArray(assumptions)
        
        if minimize_bvs is not None and len(minimize_bvs)>0:
            if self.solver.output:
                for i,n in enumerate(minimize_bvs):  
                    self._echoOutput("minimize bv " + str(minimize_bvs[i])+"\n")
                self._echoOutput("solve" + " ".join((str(dimacs(c)) for c in assumptions))+"\n")
                self.solver.output.flush()
            lp2 = (c_int * len(minimize_bvs))()
            for i,n in enumerate(minimize_bvs):            
                lp2[i]=c_int(n)
            return self.monosat_c.solveAssumptions_MinBVs(self.solver._ptr,lp,len(assumptions),lp2,len(minimize_bvs))
        else:
            if self.solver.output:
                self._echoOutput("solve" + " ".join((str(dimacs(c)) for c in assumptions))+"\n")
                self.solver.output.flush()
            return self.monosat_c.solveAssumptions(self.solver._ptr,lp,len(assumptions))
        
    def solveLimited(self,assumptions=None,minimize_bvs=None):
        self.backtrack()        
        if assumptions is None:
            assumptions=[]
            
        lp = self.getIntArray(assumptions)
        
        if minimize_bvs is not None and len(minimize_bvs)>0:
            if self.solver.output:
                for i,n in enumerate(minimize_bvs):  
                    self._echoOutput("minimize bv " + str(minimize_bvs[i])+"\n")
                self._echoOutput("solve "+ " ".join((str(dimacs(c)) for c in assumptions))+"\n")
                self.solver.output.flush()
            lp2 = (c_int * len(minimize_bvs))()
            for i,n in enumerate(minimize_bvs):            
                lp2[i]=c_int(n)
            r= self.monosat_c.solveAssumptionsLimited_MinBVs(self.solver._ptr,lp,len(assumptions),lp2,len(minimize_bvs))
        else:
            if self.solver.output:
                self._echoOutput("solve " + " ".join((str(dimacs(c)) for c in assumptions))+"\n")
                self.solver.output.flush()
            r= self.monosat_c.solveAssumptionsLimited(self.solver._ptr,lp,len(assumptions))  
              
        if r==0:
            return True
        elif r==1:
            return False
        else:
            assert(r==2)
            return None  
    
    def backtrack(self):
        if self.solver.output:
            self.solver.output.flush()
        return self.monosat_c.backtrack(self.solver._ptr)        
    
    def addUnitClause(self,clause):
        self.backtrack()        
        if isinstance(clause, int ):
            if self.solver.output:
                self._echoOutput("%d 0\n"%(dimacs(clause)))            
            self.monosat_c.addUnitClause(self.solver._ptr,clause)
        elif len(clause)==1:
            if self.solver.output:     
                self._echoOutput("%d 0\n"%(dimacs(clause[0])))       
            self.monosat_c.addUnitClause(self.solver._ptr,clause[0])
        else:
            assert(False)

    def addBinaryClause(self,l0,l1):
        self.backtrack()        
        if self.solver.output:       
            self._echoOutput(" ".join((str(dimacs(c)) for c in (l0,l1)))+" 0\n")    
        self.monosat_c.addBinaryClause(self.solver._ptr,l0,l1)           


    def addTertiaryClause(self,l0,l1,l2):
        self.backtrack()
        if self.solver.output:   
            self._echoOutput(" ".join((str(dimacs(c)) for c in (l0,l1,l2)))+" 0\n")      
        self.monosat_c.addTertiaryClause(self.solver._ptr,l0,l1,l2)  

    
    def addClause(self,clause):
        if isinstance(clause, int ):            
            self.addUnitClause(clause)
        elif len(clause)==1:          
            self.addUnitClause(clause[0])
        elif len(clause)==2:            
            self.addBinaryClause(clause[0],clause[1])            
        elif len(clause)==3:            
            self.addTertiaryClause(clause[0],clause[1],clause[2])  
        else:         
            self.backtrack()
            if self.solver.output:
                self._echoOutput(" ".join((str(dimacs(c)) for c in clause))+" 0\n")
            lp = self.getIntArray(clause)
            self.monosat_c.addClause(self.solver._ptr,lp,len(clause))  
            
    #convenience code to and together to lits and return a new 
    def addAnd(self, lit1, lit2):
        out = self.newLit()
        self.addTertiaryClause(out, self.Not(lit1), self.Not(lit2))
        self.addBinaryClause(self.Not(out), lit1)
        self.addBinaryClause(self.Not(out), lit2)
        return out
                       
    def newLit(self, allow_simplification=False):
        varID = int(self.monosat_c.newVar(self.solver._ptr))
        if not allow_simplification:
            self.monosat_c.disallowLiteralSimplification(self.solver._ptr,varID*2)
        return varID*2
    
    def setDecisionVar(self,v, b):        
        self.monosat_c.setDecisionVar(self.solver._ptr,c_var(v),c_bool(b))
    
    def isDecisionVar(self,var):      
        return self.monosat_c.isDecisionVar(self.solver._ptr,c_var(v))

    def setDecisionPriority(self,v,priority):
        self.monosat_c.setDecisionPriority(self.solver._ptr,c_var(v),c_int(priority))     
        
    def getDecisionPriority(self,v):
        self.monosat_c.getDecisionPriority(self.solver._ptr,c_var(v))     

    def setDecisionPolarity(self,v, b):        
        self.monosat_c.setDecisionPolarity(self.solver._ptr,c_var(v),c_bool(b))
    
    def getDecisionPolarity(self,var):      
        return self.monosat_c.getDecisionPolarity(self.solver._ptr,c_var(v))

    def nVars(self):
        return self.monosat_c.nVars(self.solver._ptr)
    
    def nClauses(self):
        return self.monosat_c.nClauses(self.solver._ptr)
            
    def true(self):
        return self.solver._true
    
    def false(self):
        return self.Not(self.solver._true)
    
    def Not(self,lit):       
        return lit ^1
    
    def isPositive(self,lit):
        return (~lit)&1;
    
    def isNegative(self,lit):
        return lit &1
    
    def setSymbol(self,lit,symbol):
        pass


    def getSymbol(self,lit):
        pass
    
       
    def AssertAtMostOne(self,clause):
        self.backtrack()     
        newclause=[]
        for l in clause:
            if self.isPositive(l):
                newclause.append(l)
            else:
                #Create a new, positive literal and assert that it is equal to the old, negative one.
                l2 = self.newLit()
                self.addBinaryClause(self.Not(l), l2)
                self.addBinaryClause(self.Not(l2), l)
                newclause.append(l2)
                
        if self.solver.output:
            self._echoOutput("amo " + " ".join((str(dimacs(c)) for c in newclause))+" 0\n")
        lp = self.getIntArray(newclause)        
        for i in range(len(newclause)):
            l = lp[i]
            assert(self.isPositive(l))
            assert(l%2==0)#all of the literals must be positive
            lp[i]=l//2

        self.monosat_c.at_most_one(self.solver._ptr,lp,len(newclause))  
    
    #def preprocess(self,disable_future_preprocessing=False):
    #    self.monosat_c.preprocess(disable_future_preprocessing)
    
    #bv interface

    def newBitvector_const(self, width,val):
        self.backtrack()
        bvID = self.monosat_c.newBitvector_const(self.solver._ptr,self.solver.bvtheory, width, c_long(val))
        if self.solver.output:
            self._echoOutput("bv const %d %d "%(bvID, width) + str((val)) +"\n" )   
        return bvID
    
    
    def nBitvectors(self):
        return self.monosat_c.nBitvectors(self.solver._ptr,self.solver.bvtheory)
    
    def newBitvector(self, bits):
        self.backtrack()
        arr = self.getIntArray(bits)
        bvID = self.monosat_c.newBitvector(self.solver._ptr,self.solver.bvtheory,arr,c_int(len(bits)))
        self._num_bvs=bvID+1
        if self.solver.output:       
            self._echoOutput("bv %d %d "%(bvID, len(bits)) + " ".join((str(l+1) for l in bits)) +"\n" )   
        return bvID
    
    def getTrue(self):
        self.backtrack()
        if self.solver._true is None:
            l = self.monosat_c.true_lit(self.solver._ptr)
            self.addUnitClause(l)#This isn't technically required by the solver; only doing it to ensure the unit clause gets recorded in the GNF...
            self.solver._true =l
        return self.solver._true 
    
    def newBVComparison_const_lt(self, bvID, val):
        self.backtrack()
        l = self.monosat_c.newBVComparison_const_lt(self.solver._ptr, self.solver.bvtheory,c_bvID(bvID),c_long(val))
        if self.solver.output:
            self._echoOutput("bv const < %d %d %d\n"%(dimacs(l),bvID, val))
        return l
        
    def newBVComparison_bv_lt(self, bvID1, bvID2):
        self.backtrack()
        l= self.monosat_c.newBVComparison_bv_lt(self.solver._ptr, self.solver.bvtheory, c_bvID(bvID1), c_bvID(bvID2))
        if self.solver.output:
            self._echoOutput("bv < %d %d %d\n"%(dimacs(l),bvID1,bvID2))
        return l

    def newBVComparison_const_leq(self, bvID, val):
        self.backtrack()
        l= self.monosat_c.newBVComparison_const_leq(self.solver._ptr, self.solver.bvtheory,c_bvID(bvID),c_long(val))
        if self.solver.output:
            self._echoOutput("bv const <= %d %d %d\n"%(dimacs(l),bvID, val))
        return l
    
    def newBVComparison_bv_leq(self, bvID1, bvID2):
        self.backtrack()
        l= self.monosat_c.newBVComparison_bv_leq(self.solver._ptr, self.solver.bvtheory, c_bvID(bvID1), c_bvID(bvID2))
        if self.solver.output:
            self._echoOutput("bv <= %d %d %d\n"%(dimacs(l),bvID1,bvID2))
        return l
    
    def newBVComparison_const_gt(self, bvID, val):
        self.backtrack()
        l= self.monosat_c.newBVComparison_const_gt(self.solver._ptr, self.solver.bvtheory,c_bvID(bvID),c_long(val))
        if self.solver.output:
            self._echoOutput("bv const > %d %d %d\n"%(dimacs(l),bvID, val))
        return l
    
    def newBVComparison_bv_gt(self, bvID1, bvID2):
        self.backtrack()
        l= self.monosat_c.newBVComparison_bv_gt(self.solver._ptr, self.solver.bvtheory, c_bvID(bvID1), c_bvID(bvID2))
        if self.solver.output:
            self._echoOutput("bv > %d %d %d\n"%(dimacs(l),bvID1,bvID2))
        return l
    
    def newBVComparison_const_geq(self, bvID, val):
        self.backtrack()
        l= self.monosat_c.newBVComparison_const_geq(self.solver._ptr, self.solver.bvtheory,c_bvID(bvID),c_long(val))
        if self.solver.output:
            self._echoOutput("bv const >= %d %d %d\n"%(dimacs(l),bvID, val))
        return l
    
    def newBVComparison_bv_geq(self, bvID1, bvID2):
        self.backtrack()
        l= self.monosat_c.newBVComparison_bv_geq(self.solver._ptr, self.solver.bvtheory, c_bvID(bvID1), c_bvID(bvID2))
        if self.solver.output:
            self._echoOutput("bv >= %d %d %d\n"%(dimacs(l),bvID1,bvID2))
        return l
        
    def bv_addition(self, aID,bID, resultID):
        self.backtrack()
        self.monosat_c.bv_addition(self.solver._ptr, self.solver.bvtheory, c_bvID(aID), c_bvID(bID), c_bvID(resultID))
        if self.solver.output:
            self._echoOutput("bv + %d %d %d\n"%(resultID,aID,bID))

    def bv_subtraction(self, aID,bID, resultID):
        self.backtrack()
        self.monosat_c.bv_subtraction(self.solver._ptr, self.solver.bvtheory, c_bvID(aID), c_bvID(bID), c_bvID(resultID))
        if self.solver.output:
            self._echoOutput("bv - %d %d %d\n"%(resultID,aID,bID))

    def bv_ite(self, condition_lit, thnID,elsID, resultID):
        self.backtrack()
        self.monosat_c.bv_ite(self.solver._ptr, self.solver.bvtheory,condition_lit, c_bvID(thnID), c_bvID(elsID), c_bvID(resultID))
        if self.solver.output:
            self._echoOutput("bv ite %d %d %d %d\n"%(dimacs(condition_lit),thnID,elsID,resultID))
    
    def bv_min(self, args, resultID):
        self.backtrack()        
        lp = self.getIntArray(args)                
        self.monosat_c.bv_min(self.solver._ptr, self.solver.bvtheory,lp, len(args), c_bvID(resultID))
        if self.solver.output:            
            argstr = " ".join((str(a) for a in args))
            self._echoOutput("bv min %d %d %s\n"%(resultID,len(args),argstr))

    def bv_max(self, args, resultID):
        self.backtrack()        
        lp = self.getIntArray(args)                
        self.monosat_c.bv_max(self.solver._ptr, self.solver.bvtheory, lp,len(args), c_bvID(resultID))
        if self.solver.output:            
            argstr = " ".join((str(a) for a in args))
            self._echoOutput("bv max %d %d %s\n"%(resultID,len(args),argstr))

    
    def bv_not(self, aID,resultID):
        self.backtrack()
        self.monosat_c.bv_not(self.solver._ptr, self.solver.bvtheory, c_bvID(aID),c_bvID(resultID))
        if self.solver.output:
            self._echoOutput("bv not %d %d\n"%(aID,resultID))    
    
    def bv_and(self, aID,bID,resultID):
        self.backtrack()
        self.monosat_c.bv_and(self.solver._ptr, self.solver.bvtheory, c_bvID(aID), c_bvID(bID),c_bvID(resultID))
        if self.solver.output:
            self._echoOutput("bv and %d %d %d \n"%(aID,bID,resultID))
        
    
    def bv_nand(self, aID,bID,resultID):
        self.backtrack()
        self.monosat_c.bv_nand(self.solver._ptr, self.solver.bvtheory, c_bvID(aID), c_bvID(bID),c_bvID(resultID))
        if self.solver.output:
            self._echoOutput("bv nand %d %d %d \n"%(aID,bID,resultID))
     
    
    def bv_or(self, aID,bID,resultID):
        self.backtrack()
        self.monosat_c.bv_or(self.solver._ptr, self.solver.bvtheory, c_bvID(aID), c_bvID(bID),c_bvID(resultID))
        if self.solver.output:
            self._echoOutput("bv or %d %d %d \n"%(aID,bID,resultID))
 
    
    def bv_nor(self, aID,bID,resultID):
        self.backtrack()
        self.monosat_c.bv_nor(self.solver._ptr, self.solver.bvtheory, c_bvID(aID), c_bvID(bID),c_bvID(resultID))
        if self.solver.output:
            self._echoOutput("bv nor %d %d %d \n"%(aID,bID,resultID))

        
    def bv_xor(self, aID,bID,resultID):
        self.backtrack()
        self.monosat_c.bv_xor(self.solver._ptr, self.solver.bvtheory, c_bvID(aID), c_bvID(bID),c_bvID(resultID))
        if self.solver.output:
            self._echoOutput("bv xor %d %d %d \n"%(aID,bID,resultID))

    
    def bv_xnor(self, aID,bID,resultID):
        self.backtrack()
        self.monosat_c.bv_xnor(self.solver._ptr, self.solver.bvtheory, c_bvID(aID), c_bvID(bID),c_bvID(resultID))
        if self.solver.output:
            self._echoOutput("bv xnor %d %d %d \n"%(aID,bID,resultID))


    def bv_concat(self, aID,bID,resultID):
        self.backtrack()
        self.monosat_c.bv_concat(self.solver._ptr, self.solver.bvtheory, c_bvID(aID), c_bvID(bID),c_bvID(resultID))
        if self.solver.output:
            self._echoOutput("bv concat %d %d %d\n"%(aID,bID,resultID))


    def bv_slice(self, aID,lower, upper,resultID):
        self.backtrack()
        self.monosat_c.bv_slice(self.solver._ptr, self.solver.bvtheory, c_bvID(aID),lower,upper,c_bvID(resultID))
        if self.solver.output:
            self._echoOutput("bv slice %d %d %d %d\n"%(aID,lower,upper, resultID))
            
    def bv_popcount(self,args,resultID):
        self.backtrack()         
        newargs=[]
        for l in args:
            if self.isPositive(l):
                #newargs.append(l)
                l2 = self.newLit()
                self.addBinaryClause(self.Not(l), l2)
                self.addBinaryClause(self.Not(l2), l)
                newargs.append(l2)
            else:
                #Create a new, positive literal and assert that it is equal to the old, negative one.
                l2 = self.newLit()
                self.addBinaryClause(self.Not(l), l2)
                self.addBinaryClause(self.Not(l2), l)
                newargs.append(l2)
        
        lp = self.getIntArray(newargs)          
        self.monosat_c.bv_popcount(self.solver._ptr, self.solver.bvtheory,lp,len(newargs),c_bvID(resultID))
        if self.solver.output:
            self._echoOutput("bv popcount %d %d %s\n"%(resultID, len(newargs)," ".join((str(dimacs(l)) for l in newargs))))
    #Monosat graph interface
    
    def newGraph(self):
        self.backtrack()
        g = self.monosat_c.newGraph(self.solver._ptr)
        gid = len(self.solver.graphs)
        self.solver.graphs.append(g)
        self.solver.graph_ids[g]=gid
        if self.solver.output:
            self._echoOutput("digraph 0 0 %d\n"%(gid)) 
        
        return g
    
    def getGraph(self,id):
        return self.solver.graphs[id]
    
    def nGraphs(self):
        return len(self.solver.graphs)
    
    def hasGraph(self,gID):
        return gID in self.solver.graph_ids
    
    def newNode(self, graph):
        self.backtrack()
        return self.monosat_c.newNode(self.solver._ptr,graph)

    def newEdge(self, graph, u,v, weight):
        self.backtrack()
        l = self.monosat_c.newEdge(self.solver._ptr,graph,c_int(u),c_int(v),c_long(weight))
        if self.solver.output:
            self._echoOutput("edge " + str(self.getGID(graph)) + " " + str(u) + " " + str(v) + " " +  str(dimacs(l)) + " " + str((weight))  + "\n")
        return l
    
    def newEdge_double(self, graph, u,v, weight):
        self.backtrack()
        l = self.monosat_c.newEdge_double(self.solver._ptr,graph,c_int(u),c_int(v),c_double(weight))
        if self.solver.output:
            self._echoOutput("edge " + str(self.getGID(graph)) + " " + str(u) + " " + str(v) + " " +  str(dimacs(l)) + " " + str((weight))  + "\n")
        return l
        
    def newEdge_bv(self, graph, u,v, bvID):
        self.backtrack()
        l= self.monosat_c.newEdge_bv(self.solver._ptr,graph,c_int(u),c_int(v),c_bvID(bvID))
        if self.solver.output:
            self._echoOutput("edge_bv " + str(self.getGID(graph)) + " " + str(u) + " " + str(v) + " " +  str(dimacs(l)) + " " + str((bvID))  + "\n")
        return l

    
    def reaches(self, graph, u,v):
        self.backtrack()
        l= self.monosat_c.reaches(self.solver._ptr,graph,c_int(u),c_int(v))
        if self.solver.output:
            self._echoOutput("reach " + str(self.getGID(graph)) + " " + str(u) + " " + str(v) + " " + str(dimacs(l)) + "\n" );
        return l

    def shortestPathUnweighted_lt_const(self, graph, u,v,dist):
        self.backtrack()
        l= self.monosat_c.shortestPathUnweighted_lt_const(self.solver._ptr,graph,c_int(u),c_int(v),c_long(dist))
        if self.solver.output:
            self._echoOutput("distance_lt " + str(self.getGID(graph)) + " " + str(u) + " " + str(v) + " " + str(dimacs(l)) + " " + str((dist)) + "\n" );
        return l
    
    def shortestPathUnweighted_leq_const(self, graph, u,v,dist):
        self.backtrack()
        l= self.monosat_c.shortestPathUnweighted_leq_const(self.solver._ptr,graph,c_int(u),c_int(v),c_long(dist))
        if self.solver.output:
            self._echoOutput("distance_leq " + str(self.getGID(graph)) + " " + str(u) + " " + str(v) + " " + str(dimacs(l)) + " " + str((dist)) + "\n" );
        return l

    def shortestPath_lt_const(self, graph, u,v,dist):
        self.backtrack()
        l= self.monosat_c.shortestPath_lt_const(self.solver._ptr,graph,c_int(u),c_int(v),c_long(dist))
        if self.solver.output:
            self._echoOutput("weighted_distance_lt " + str(self.getGID(graph)) + " " + str(u) + " " + str(v) + " " + str(dimacs(l)) + " " + str((dist)) + "\n" );
        return l
    
    def shortestPath_leq_const(self, graph, u,v,dist):
        self.backtrack()
        l= self.monosat_c.shortestPath_leq_const(self.solver._ptr,graph,c_int(u),c_int(v),c_long(dist))
        if self.solver.output:
            self._echoOutput("weighted_distance_leq " + str(self.getGID(graph)) + " " + str(u) + " " + str(v) + " " + str(dimacs(l)) + " " + str((dist)) + "\n" );
        return l
    
    def shortestPath_lt_bv(self, graph, u,v,bvID):
        self.backtrack()
        l= self.monosat_c.shortestPath_lt_bv(self.solver._ptr,graph,c_int(u),c_int(v),c_bvID(bvID))
        if self.solver.output:
            self._echoOutput("weighted_distance_bv_lt " + str(self.getGID(graph)) + " " + str(u) + " " + str(v) + " " + str(dimacs(l)) + " " + str((bvID)) + "\n" );
        return l
                  
    def shortestPath_leq_bv(self, graph, u,v,bvID):
        self.backtrack()
        l= self.monosat_c.shortestPath_leq_bv(self.solver._ptr,graph,c_int(u),c_int(v),c_bvID(bvID))
        if self.solver.output:
            self._echoOutput("weighted_distance_bv_leq " + str(self.getGID(graph)) + " " + str(u) + " " + str(v) + " " + str(dimacs(l)) + " " + str((bvID)) + "\n" );
        return l
    
    def maximumFlow_geq(self, graph, s,t,flow):
        self.backtrack()
        l= self.monosat_c.maximumFlow_geq(self.solver._ptr,graph,c_int(s),c_int(t),c_long(flow))
        if self.solver.output:  
            self._echoOutput("maximum_flow_geq " + str(self.getGID(graph)) + " " + str(s) + " " + str(t) + " " + str(dimacs(l)) + " " + str((flow)) + "\n" );
        return l
    
    def maximumFlow_gt(self, graph, s,t,flow):
        self.backtrack()
        l= self.monosat_c.maximumFlow_gt(self.solver._ptr,graph,c_int(s),c_int(t),c_long(flow))
        if self.solver.output: 
            self._echoOutput("maximum_flow_gt " + str(self.getGID(graph)) + " " + str(s) + " " + str(t) + " " + str(dimacs(l)) + " " + str((flow)) + "\n" );
        return l
               
    def maximumFlow_geq_bv(self, graph, s,t,bvID):
        self.backtrack()
        l= self.monosat_c.maximumFlow_geq_bv(self.solver._ptr,graph,c_int(s),c_int(t),c_bvID(bvID))
        if self.solver.output:
            self._echoOutput("maximum_flow_bv_geq " + str(self.getGID(graph)) + " " + str(s) + " " + str(t) + " " + str(dimacs(l)) + " " + str((bvID)) + "\n" );
        return l
    
    def maximumFlow_gt_bv(self, graph, s,t,bvID):
        self.backtrack()
        l= self.monosat_c.maximumFlow_gt_bv(self.solver._ptr,graph,c_int(s),c_int(t),c_bvID(bvID))
        if self.solver.output:
            self._echoOutput("maximum_flow_bv_gt " + str(self.getGID(graph)) + " " + str(s) + " " + str(t) + " " + str(dimacs(l)) + " " + str((bvID)) + "\n" );
        return l
    
    def minimumSpanningTree_leq(self, graph,weight):
        self.backtrack()
        l= self.monosat_c.minimumSpanningTree_leq(self.solver._ptr,graph,c_long(weight))
        if self.solver.output:      
            self._echoOutput("mst_weight_leq " + str(self.getGID(graph)) + " " +  str(dimacs(l)) + " " + str((weight)) + "\n" );
        return l
            
    def minimumSpanningTree_lt(self, graph,weight):
        self.backtrack()
        l= self.monosat_c.minimumSpanningTree_lt(self.solver._ptr,graph,c_long(weight))
        if self.solver.output:      
            self._echoOutput("mst_weight_lt " + str(self.getGID(graph)) + " " +  str(dimacs(l)) + " " + str((weight)) + "\n" );
        return l    
    
    def acyclic_undirected(self, graph):
        self.backtrack()
        l= self.monosat_c.acyclic_undirected(self.solver._ptr,graph)
        if self.solver.output:  
            self._echoOutput("forest " + str(self.getGID(graph)) + " " +  str(dimacs(l)) + "\n" );
        return l  
                 
    def acyclic_directed(self, graph):
        self.backtrack()
        l= self.monosat_c.acyclic_directed(self.solver._ptr,graph)
        if self.solver.output:                      
            self._echoOutput("acyclic " + str(self.getGID(graph)) + " " +  str(dimacs(l)) + "\n" );
        return l  
  
    #0 = true, 1=false, 2=unassigned
    def getModel_Literal(self, lit):
        return self.monosat_c.getModel_Literal(self.solver._ptr, lit);

    def getModel_BV(self, bvID,getMaximumValue=False):
        return self.monosat_c.getModel_BV(self.solver._ptr, self.solver.bvtheory,c_bvID(bvID),c_bool(getMaximumValue));
        
    def getModel_MaxFlow(self, graph, flowlit):
        return self.monosat_c.getModel_MaxFlow(self.solver._ptr, graph,flowlit);
        
    def getModel_EdgeFlow(self, graph, flowlit,edgelit):
        return self.monosat_c.getModel_EdgeFlow(self.solver._ptr, graph,flowlit, edgelit);

    def getModel_AcyclicEdgeFlow(self, graph, flowlit,edgelit):
        return self.monosat_c.getModel_AcyclicEdgeFlow(self.solver._ptr, graph,flowlit, edgelit);


    def getModel_MinimumSpanningTreeWeight(self, graph, mstlit):
        return self.monosat_c.getModel_MinimumSpanningTreeWeight(self.solver._ptr, graph,mstlit); 



    
        