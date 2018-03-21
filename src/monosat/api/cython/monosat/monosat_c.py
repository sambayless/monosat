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

import os
import platform

from monosat.singleton import Singleton
from os import path
from enum import Enum
#module_path = os.path.abspath(path.dirname(__file__))


import pyximport
pyximport.install()
import monosat.monosat_p
print("Loading cython monosat library")

def dimacs(l):
    assert(isinstance(l,int))
    if l & 1:
        return -(l//2 +1)
    else:
        return l//2+1

class Ineq(Enum):
    LT=-2
    LEQ=-1
    EQ=0
    GEQ=1
    GT=2

class Solver():
    def __init__(self,monosat_c,arguments=None):
        self.elapsed_time=0
        if arguments is None:
            arguments=[]    
        elif isinstance(arguments,str):
            #split this up by whitespace
            arguments=arguments.split()  
        self.arguments=arguments 
        #arguments=["monosat"]+arguments
        arr = "monosat " + " ".join(arguments)
        #for i,arg in enumerate(arguments):
        #    arr[i]=(bytes(arg, 'utf-8'))

        #In the future, allow multiple solvers to be instantiated...
        self.comments=[]

        self._ptr = monosat_c.newSolver_arg(bytes(arr, 'utf-8'))


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
        self.monosat_c = monosat.monosat_p
        self.elapsed_time=0

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
    
    def getVersion(self):
        return self.monosat_c.getVersion()

    #Until a better system is created, this can be used to re-initialize Monosat with a new configuration.
    def newSolver(self,arguments=None):
        self.solver = Solver(self.monosat_c,arguments)
        self.solver._true = self.getTrue()
        return self.solver


    def readGNF(self, filename):
        self.monosat_c.readGNF(self.solver._ptr,filename.encode('ascii'))
    
    def getOutputFile(self):
        return self.solver.output
    
    def setOutputFile(self, filename):
        self.monosat_c.setOutputFile(self.solver._ptr,(filename.encode('ascii')))
        # self.solver.output=file
        # if file is not None:#Otherwise, this unit clause will have been skipped in the output file
        #     if self.solver.arguments is not None:
        #         self._echoOutput("c monosat " + " ".join(self.solver.arguments) +"\n")
        #     self.addUnitClause(self.solver._true)
    
    def _echoOutput(self,line):
        if self.solver.output:
            self.solver.output.write(line)
    
    def comment(self,c):
        self.solver.has_comments=True
        self.solver.comments.append(c)
        if self.solver.output:
            self._echoOutput("c " + c +"\n")           

    def getEmptyIntArray(self,length):                 
        if length>len(self._int_array):
            self._int_array = (c_int * length)()
        return self._int_array
    
    def getIntArray(self,nums):
        #if len(nums)>len(self._int_array):
        #    self._int_array = (c_int * len(nums))()
        #for i,n in enumerate(nums):
        #    self._int_array[i]=(n)
        return nums

    def getIntArray2(self,nums):
        #if len(nums)>len(self._int_array2):
        #    self._int_array2 = (c_int * len(nums))()
        #for i,n in enumerate(nums):
        #    self._int_array2[i]=(n)
        #return self._int_array2
        return nums

    def getLongArray(self,nums):
        #if len(nums)>len(self._long_array):
        #   self._long_array = ( * len(nums))()
        #for i,n in enumerate(nums):
        #    self._long_array[i]=(n)
        #return self._long_array
        return nums

    def intArrayToList(self, array_pointer,length):
        ret = []
        for i in range(length):
            ret.append(array_pointer[i])
        return ret        

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

    def lastSolutionWasOptimal(self):
        return self.monosat_c.lastSolutionWasOptimal(self.solver._ptr)

    def getConflictClause(self):
        conflict_size =  self.monosat_c.getConflictClause(self.solver._ptr,null_ptr,0)
        if conflict_size<0:
            return None #No conflict clause in the solver
        
 
        conflict_ptr = self.getEmptyIntArray(conflict_size)
        length = self.monosat_c.getConflictClause(self.solver._ptr,conflict_ptr,conflict_size)
        if length != conflict_size:
            raise RuntimeError("Error reading conflict clause")
        
        return self.intArrayToList(conflict_ptr,conflict_size)



    def solve(self,assumptions=None):
        self.backtrack()
        if assumptions is None:
            assumptions=[]

        lp = self.getIntArray(assumptions)

        if self.solver.output:
            self._echoOutput("solve" + " ".join((str(dimacs(c)) for c in assumptions))+"\n")
            self.solver.output.flush()
        return self.monosat_c.solveAssumptions(self.solver._ptr,lp,len(assumptions))
        
    def solveLimited(self,assumptions=None):
        self.backtrack()        
        if assumptions is None:
            assumptions=[]
            
        lp = self.getIntArray(assumptions)
        

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

        if self.solver.output:       
            self._echoOutput(" ".join((str(dimacs(c)) for c in (l0,l1)))+" 0\n")    
        self.monosat_c.addBinaryClause(self.solver._ptr,l0,l1)

    #clauses should be a list of pairs of literals, each of which will be added as a binary clause
    def addBinaryClauses(self,clauses):

        if len(clauses)==1:
            self.addBinaryClause(clauses[0][0],clauses[0][1])
        else:
            firsts, seconds = zip(*clauses)

            first_p = self.getIntArray(firsts)
            second_p = self.getIntArray2(seconds)
            self.monosat_c.addBinaryClauses(self.solver._ptr,first_p,second_p,len(clauses))


    def addTertiaryClause(self,l0,l1,l2):

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

            if self.solver.output:
                self._echoOutput(" ".join((str(dimacs(c)) for c in clause))+" 0\n")
            lp = self.getIntArray(clause)
            self.monosat_c.addClause(self.solver._ptr,lp,len(clause))

    def clearOptimizationObjectives(self):
        if self.solver.output:
            self._echoOutput("clear_opt\n")
        self.monosat_c.clearOptimizationObjectives(self.solver._ptr)

    def maximizeBV(self, bvID):
        if self.solver.output:
            self._echoOutput("maximize bv %d\n"%(bvID))

        self.monosat_c.maximizeBV(self.solver._ptr, self.solver.bvtheory, (bvID))

    def minimizeBV(self, bvID):
        if self.solver.output:
            self._echoOutput("minimize bv %d\n"%(bvID))
        self.monosat_c.minimizeBV(self.solver._ptr, self.solver.bvtheory, (bvID))

    def maximizeLits(self, lits):
        lp = self.getIntArray(lits)
        self.monosat_c.maximizeLits(self.solver._ptr, lp, len(lits))

    def minimizeLits(self, lits):
        lp = self.getIntArray(lits)
        self.monosat_c.minimizeLits(self.solver._ptr, lp, len(lits))

    def maximizeWeightedLits(self, lits,weights):
        lp = self.getIntArray(lits)
        lp2 = self.getIntArray2(weights)
        self.monosat_c.maximizeWeightedLits(self.solver._ptr, lp,lp2, len(lits))

    def minimizeWeightedLits(self, lits,weights):
        lp = self.getIntArray(lits)
        lp2 = self.getIntArray2(weights)
        self.monosat_c.minimizeWeightedLits(self.solver._ptr, lp,lp2, len(lits))

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
        self.monosat_c.setDecisionVar(self.solver._ptr,(v),(b))
    
    def isDecisionVar(self,var):      
        return self.monosat_c.isDecisionVar(self.solver._ptr,(v))

    def setDecisionPriority(self,v,priority):
        self.monosat_c.setDecisionPriority(self.solver._ptr,(v),(priority))     
        
    def getDecisionPriority(self,v):
        self.monosat_c.getDecisionPriority(self.solver._ptr,(v))     

    def setDecisionPolarity(self,v, b):        
        self.monosat_c.setDecisionPolarity(self.solver._ptr,(v),(b))
    
    def getDecisionPolarity(self,var):      
        return self.monosat_c.getDecisionPolarity(self.solver._ptr,(v))

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


    def pbOpToStr(self,op):
        if op == Ineq.LT:
            return "<"
        elif op==Ineq.LEQ:
            return "<="
        elif op==Ineq.EQ:
            return "=="
        elif op==Ineq.GEQ:
            return ">="
        elif op==Ineq.GT:
            return ">"

    def AssertPB(self, lits, coefs, op, rhs):
        self.backtrack()
        lp = self.getIntArray(lits)
        lp2 = self.getIntArray2(coefs)
        crhs = (rhs)
        assert(len(lits) == len(coefs))
        if self.solver.output:
            self._echoOutput("pb " + self.pbOpToStr(op) + " %d %d "%(rhs, len(lits))   + " ".join((str(dimacs(c)) for c in lits)) + " " + str(len(lits)) + " " +  " ".join((str(w) for w in coefs)))
        if op==Ineq.LT:
            self.monosat_c.assertPB_lt(self.solver._ptr,crhs,len(lits),lp,lp2)
        elif op==Ineq.LEQ:
            self.monosat_c.assertPB_leq(self.solver._ptr,crhs,len(lits),lp,lp2)
        elif op==Ineq.EQ:
            self.monosat_c.assertPB_eq(self.solver._ptr,crhs,len(lits),lp,lp2)
        elif op==Ineq.GEQ:
            self.monosat_c.assertPB_geq(self.solver._ptr,crhs,len(lits),lp,lp2)
        elif op==Ineq.GT:
            self.monosat_c.assertPB_gt(self.solver._ptr,crhs,len(lits),lp,lp2)

    #Convert any psuedo-boolean constraints into cnf in the solver (will be called automatically during solving)
    def flushPB(self):
        self.monosat_c.flushPB(self.solver._ptr)

    #def preprocess(self,disable_future_preprocessing=False):
    #    self.monosat_c.preprocess(disable_future_preprocessing)
    
    #bv interface
    def newBitvector_anon(self,width):
        self.backtrack()
        bvID = self.monosat_c.newBitvector_anon(self.solver._ptr,self.solver.bvtheory, width)
        if self.solver.output:
            self._echoOutput("bv anon %d %d"%(bvID, width)+"\n" )   
        return bvID        
    
    def newBitvector_const(self, width,val):
        self.backtrack()
        bvID = self.monosat_c.newBitvector_const(self.solver._ptr,self.solver.bvtheory, width, (val))
        if self.solver.output:
            self._echoOutput("bv const %d %d "%(bvID, width) + str((val)) +"\n" )   
        return bvID
    
    
    def nBitvectors(self):
        return self.monosat_c.nBitvectors(self.solver._ptr,self.solver.bvtheory)
    
    def newBitvector(self, bits):
        self.backtrack()
        arr = self.getIntArray(bits)
        bvID = self.monosat_c.newBitvector(self.solver._ptr,self.solver.bvtheory,arr,(len(bits)))
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
        l = self.monosat_c.newBVComparison_const_lt(self.solver._ptr, self.solver.bvtheory,(bvID),(val))
        if self.solver.output:
            self._echoOutput("bv const < %d %d %d\n"%(dimacs(l),bvID, val))
        return l
        
    def newBVComparison_bv_lt(self, bvID1, bvID2):
        self.backtrack()
        l= self.monosat_c.newBVComparison_bv_lt(self.solver._ptr, self.solver.bvtheory, (bvID1), (bvID2))
        if self.solver.output:
            self._echoOutput("bv < %d %d %d\n"%(dimacs(l),bvID1,bvID2))
        return l

    def newBVComparison_const_leq(self, bvID, val):
        self.backtrack()
        l= self.monosat_c.newBVComparison_const_leq(self.solver._ptr, self.solver.bvtheory,(bvID),(val))
        if self.solver.output:
            self._echoOutput("bv const <= %d %d %d\n"%(dimacs(l),bvID, val))
        return l
    
    def newBVComparison_bv_leq(self, bvID1, bvID2):
        self.backtrack()
        l= self.monosat_c.newBVComparison_bv_leq(self.solver._ptr, self.solver.bvtheory, (bvID1), (bvID2))
        if self.solver.output:
            self._echoOutput("bv <= %d %d %d\n"%(dimacs(l),bvID1,bvID2))
        return l
    
    def newBVComparison_const_gt(self, bvID, val):
        self.backtrack()
        l= self.monosat_c.newBVComparison_const_gt(self.solver._ptr, self.solver.bvtheory,(bvID),(val))
        if self.solver.output:
            self._echoOutput("bv const > %d %d %d\n"%(dimacs(l),bvID, val))
        return l
    
    def newBVComparison_bv_gt(self, bvID1, bvID2):
        self.backtrack()
        l= self.monosat_c.newBVComparison_bv_gt(self.solver._ptr, self.solver.bvtheory, (bvID1), (bvID2))
        if self.solver.output:
            self._echoOutput("bv > %d %d %d\n"%(dimacs(l),bvID1,bvID2))
        return l
    
    def newBVComparison_const_geq(self, bvID, val):
        self.backtrack()
        l= self.monosat_c.newBVComparison_const_geq(self.solver._ptr, self.solver.bvtheory,(bvID),(val))
        if self.solver.output:
            self._echoOutput("bv const >= %d %d %d\n"%(dimacs(l),bvID, val))
        return l
    
    def newBVComparison_bv_geq(self, bvID1, bvID2):
        self.backtrack()
        l= self.monosat_c.newBVComparison_bv_geq(self.solver._ptr, self.solver.bvtheory, (bvID1), (bvID2))
        if self.solver.output:
            self._echoOutput("bv >= %d %d %d\n"%(dimacs(l),bvID1,bvID2))
        return l
        
    def bv_addition(self, aID,bID, resultID):
        self.backtrack()
        self.monosat_c.bv_addition(self.solver._ptr, self.solver.bvtheory, (aID), (bID), (resultID))
        if self.solver.output:
            self._echoOutput("bv + %d %d %d\n"%(resultID,aID,bID))

    def bv_subtraction(self, aID,bID, resultID):
        self.backtrack()
        self.monosat_c.bv_subtraction(self.solver._ptr, self.solver.bvtheory, (aID), (bID), (resultID))
        if self.solver.output:
            self._echoOutput("bv - %d %d %d\n"%(resultID,aID,bID))

    def bv_multiply(self, aID,bID, resultID):
        self.backtrack()
        self.monosat_c.bv_multiply(self.solver._ptr, self.solver.bvtheory, (aID), (bID), (resultID))
        if self.solver.output:
            self._echoOutput("bv * %d %d %d\n"%(resultID,aID,bID))

    def bv_divide(self, aID,bID, resultID):
        self.backtrack()
        self.monosat_c.bv_divide(self.solver._ptr, self.solver.bvtheory, (aID), (bID), (resultID))
        if self.solver.output:
            self._echoOutput("bv / %d %d %d\n"%(resultID,aID,bID))

    def bv_ite(self, condition_lit, thnID,elsID, resultID):
        self.backtrack()
        self.monosat_c.bv_ite(self.solver._ptr, self.solver.bvtheory,condition_lit, (thnID), (elsID), (resultID))
        if self.solver.output:
            self._echoOutput("bv ite %d %d %d %d\n"%(dimacs(condition_lit),thnID,elsID,resultID))
    
    def bv_min(self, args, resultID):
        self.backtrack()        
        lp = self.getIntArray(args)                
        self.monosat_c.bv_min(self.solver._ptr, self.solver.bvtheory,lp, len(args), (resultID))
        if self.solver.output:            
            argstr = " ".join((str(a) for a in args))
            self._echoOutput("bv min %d %d %s\n"%(resultID,len(args),argstr))

    def bv_max(self, args, resultID):
        self.backtrack()        
        lp = self.getIntArray(args)                
        self.monosat_c.bv_max(self.solver._ptr, self.solver.bvtheory, lp,len(args), (resultID))
        if self.solver.output:            
            argstr = " ".join((str(a) for a in args))
            self._echoOutput("bv max %d %d %s\n"%(resultID,len(args),argstr))

    
    def bv_not(self, aID,resultID):
        self.backtrack()
        self.monosat_c.bv_not(self.solver._ptr, self.solver.bvtheory, (aID),(resultID))
        if self.solver.output:
            self._echoOutput("bv not %d %d\n"%(aID,resultID))    
    
    def bv_and(self, aID,bID,resultID):
        self.backtrack()
        self.monosat_c.bv_and(self.solver._ptr, self.solver.bvtheory, (aID), (bID),(resultID))
        if self.solver.output:
            self._echoOutput("bv and %d %d %d \n"%(aID,bID,resultID))
        
    
    def bv_nand(self, aID,bID,resultID):
        self.backtrack()
        self.monosat_c.bv_nand(self.solver._ptr, self.solver.bvtheory, (aID), (bID),(resultID))
        if self.solver.output:
            self._echoOutput("bv nand %d %d %d \n"%(aID,bID,resultID))
     
    
    def bv_or(self, aID,bID,resultID):
        self.backtrack()
        self.monosat_c.bv_or(self.solver._ptr, self.solver.bvtheory, (aID), (bID),(resultID))
        if self.solver.output:
            self._echoOutput("bv or %d %d %d \n"%(aID,bID,resultID))
 
    
    def bv_nor(self, aID,bID,resultID):
        self.backtrack()
        self.monosat_c.bv_nor(self.solver._ptr, self.solver.bvtheory, (aID), (bID),(resultID))
        if self.solver.output:
            self._echoOutput("bv nor %d %d %d \n"%(aID,bID,resultID))

        
    def bv_xor(self, aID,bID,resultID):
        self.backtrack()
        self.monosat_c.bv_xor(self.solver._ptr, self.solver.bvtheory, (aID), (bID),(resultID))
        if self.solver.output:
            self._echoOutput("bv xor %d %d %d \n"%(aID,bID,resultID))

    
    def bv_xnor(self, aID,bID,resultID):
        self.backtrack()
        self.monosat_c.bv_xnor(self.solver._ptr, self.solver.bvtheory, (aID), (bID),(resultID))
        if self.solver.output:
            self._echoOutput("bv xnor %d %d %d \n"%(aID,bID,resultID))


    def bv_concat(self, aID,bID,resultID):
        self.backtrack()
        self.monosat_c.bv_concat(self.solver._ptr, self.solver.bvtheory, (aID), (bID),(resultID))
        if self.solver.output:
            self._echoOutput("bv concat %d %d %d\n"%(aID,bID,resultID))

    def bv_bitblast(self,bvID):
        self.backtrack();
        self.monosat_c.bv_bitblast(self.solver._ptr, self.solver.bvtheory, (bvID))
        if self.solver.output:
            self._echoOutput("bv bitblast %d\n"%(bvID))

    def bv_slice(self, aID,lower, upper,resultID):
        self.backtrack()
        self.monosat_c.bv_slice(self.solver._ptr, self.solver.bvtheory, (aID),lower,upper,(resultID))
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
        self.monosat_c.bv_popcount(self.solver._ptr, self.solver.bvtheory,lp,len(newargs),(resultID))
        if self.solver.output:
            self._echoOutput("bv popcount %d %d %s\n"%(resultID, len(newargs)," ".join((str(dimacs(l)) for l in newargs))))


    def bv_unary(self,args,resultID):
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
        self.monosat_c.bv_unary(self.solver._ptr, self.solver.bvtheory,lp,len(newargs),(resultID))
        if self.solver.output:
            self._echoOutput("bv unary %d %d %s\n"%(resultID, len(newargs)," ".join((str(dimacs(l)) for l in newargs))))

    #Monosat fsm interface

    def newFSM(self, in_labels,out_labels):
        return self.monosat_c.newFSM.argtypes(self.solver._ptr, null_ptr,in_labels,out_labels)

    def newState(self, fsm_id):
        return self.monosat_c.newState(self.solver._ptr,null_ptr, fsm_id)

    def newTransition(self, fsm_id, from_state, to_state,in_label,out_label):
        return self.monosat_c.newTransition(self.solver._ptr,null_ptr, fsm_id)

    def newTransition(self, fsm_id, from_state, to_state,in_label,out_label):
        return self.monosat_c.newTransition(self.solver._ptr,null_ptr, fsm_id)

    #A string is an array of positive integers
    def newString(self, string_int_array):
        lp = self.getIntArray(string_int_array)
        return self.monosat_c.newString(self.solver._ptr,null_ptr, lp, len(string_int_array))

    def fsmAcceptsString(self, fsm_id, starting_state, accepting_state, strID):
        return self.monosat_c.fsmAcceptsString(self.solver._ptr,null_ptr, fsm_id,starting_state,accepting_state,strID)

    def fsmCompositionAccepts(self, fsm_generator_id, fsm_acceptor_id, gen_starting_state, gen_accepting_state, accept_starting_state, accept_accepting_state, strID):
        return self.monosat_c.fsmCompositionAccepts(self.solver._ptr,null_ptr, fsm_generator_id,fsm_acceptor_id, gen_starting_state, gen_accepting_state, accept_starting_state, accept_accepting_state, strID)


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

    def nNodes(self, graph):
        return self.monosat_c.nNodes(self.solver._ptr,graph)

    def nEdges(self, graph):
        return self.monosat_c.nEdges(self.solver._ptr,graph)

    def newEdge(self, graph, u,v, weight):
        #self.backtrack()
        l = self.monosat_c.newEdge(self.solver._ptr,graph,(u),(v),(weight))
        if self.solver.output:
            self._echoOutput("edge " + str(self.getGID(graph)) + " " + str(u) + " " + str(v) + " " +  str(dimacs(l)) + " " + str((weight))  + "\n")
        return l
    
    def newEdge_double(self, graph, u,v, weight):
        self.backtrack()
        l = self.monosat_c.newEdge_double(self.solver._ptr,graph,(u),(v),(weight))
        if self.solver.output:
            self._echoOutput("edge " + str(self.getGID(graph)) + " " + str(u) + " " + str(v) + " " +  str(dimacs(l)) + " " + str((weight))  + "\n")
        return l
        
    def newEdge_bv(self, graph, u,v, bvID):
        self.backtrack()
        l= self.monosat_c.newEdge_bv(self.solver._ptr,graph,(u),(v),(bvID))
        if self.solver.output:
            self._echoOutput("edge_bv " + str(self.getGID(graph)) + " " + str(u) + " " + str(v) + " " +  str(dimacs(l)) + " " + str((bvID))  + "\n")
        return l

    def newEdgeSet(self,graph,edges,enforceEdgeAssignments=True):
        self.backtrack()
        if self.solver.output:
            edgestr = "edge_set %d %d "%(self.getGID(graph), len(edges))
            self._echoOutput(edgestr + " ".join((str(dimacs(c)) for c in edges))+"\n")
        lp = self.getIntArray(edges)
        self.monosat_c.newEdgeSet(self.solver._ptr,graph,lp,len(edges), (enforceEdgeAssignments))


    def assignWeightsTo(self,graph,weight):
        self.monosat_c.graph_setAssignEdgesToWeight(self.solver._ptr,graph,(weight))


    def enforceRouting(self,graph,source,destination,nets,maxflowlit):
        r_ptr = self.monosat_c.createFlowRouting(self.solver._ptr,graph,(source),(destination),maxflowlit)
        for dest_edge_lits,net_reach_lits,disabled_edge_lit in nets:
            lp = self.getIntArray(dest_edge_lits)
            lp2 = self.getIntArray2(net_reach_lits)
            self.monosat_c.addRoutingNet(self.solver._ptr,graph,r_ptr,disabled_edge_lit, len(dest_edge_lits),lp,lp2)


    def reaches(self, graph, u,v):
        self.checkNode(graph,u);
        self.checkNode(graph,v);
        self.backtrack()
        l= self.monosat_c.reaches(self.solver._ptr,graph,(u),(v))
        if self.solver.output:
            self._echoOutput("reach " + str(self.getGID(graph)) + " " + str(u) + " " + str(v) + " " + str(dimacs(l)) + "\n" );
        return l

    def shortestPathUnweighted_lt_const(self, graph, u,v,dist):
        self.checkNode(graph,u);
        self.checkNode(graph,v);
        self.backtrack()
        l= self.monosat_c.shortestPathUnweighted_lt_const(self.solver._ptr,graph,(u),(v),(dist))
        if self.solver.output:
            self._echoOutput("distance_lt " + str(self.getGID(graph)) + " " + str(u) + " " + str(v) + " " + str(dimacs(l)) + " " + str((dist)) + "\n" );
        return l
    
    def shortestPathUnweighted_leq_const(self, graph, u,v,dist):
        self.checkNode(graph,u);
        self.checkNode(graph,v);
        self.backtrack()
        l= self.monosat_c.shortestPathUnweighted_leq_const(self.solver._ptr,graph,(u),(v),(dist))
        if self.solver.output:
            self._echoOutput("distance_leq " + str(self.getGID(graph)) + " " + str(u) + " " + str(v) + " " + str(dimacs(l)) + " " + str((dist)) + "\n" );
        return l

    def shortestPath_lt_const(self, graph, u,v,dist):
        self.checkNode(graph,u);
        self.checkNode(graph,v);
        self.backtrack()
        l= self.monosat_c.shortestPath_lt_const(self.solver._ptr,graph,(u),(v),(dist))
        if self.solver.output:
            self._echoOutput("weighted_distance_lt " + str(self.getGID(graph)) + " " + str(u) + " " + str(v) + " " + str(dimacs(l)) + " " + str((dist)) + "\n" );
        return l
    
    def shortestPath_leq_const(self, graph, u,v,dist):
        self.checkNode(graph,u);
        self.checkNode(graph,v);
        self.backtrack()
        l= self.monosat_c.shortestPath_leq_const(self.solver._ptr,graph,(u),(v),(dist))
        if self.solver.output:
            self._echoOutput("weighted_distance_leq " + str(self.getGID(graph)) + " " + str(u) + " " + str(v) + " " + str(dimacs(l)) + " " + str((dist)) + "\n" );
        return l
    
    def shortestPath_lt_bv(self, graph, u,v,bvID):
        self.checkNode(graph,u);
        self.checkNode(graph,v);
        self.checkBV(bvID);
        self.backtrack()
        l= self.monosat_c.shortestPath_lt_bv(self.solver._ptr,graph,(u),(v),(bvID))
        if self.solver.output:
            self._echoOutput("weighted_distance_bv_lt " + str(self.getGID(graph)) + " " + str(u) + " " + str(v) + " " + str(dimacs(l)) + " " + str((bvID)) + "\n" );
        return l
                  
    def shortestPath_leq_bv(self, graph, u,v,bvID):
        self.checkNode(graph,u);
        self.checkNode(graph,v);
        self.checkBV(bvID);
        self.backtrack()
        l= self.monosat_c.shortestPath_leq_bv(self.solver._ptr,graph,(u),(v),(bvID))
        if self.solver.output:
            self._echoOutput("weighted_distance_bv_leq " + str(self.getGID(graph)) + " " + str(u) + " " + str(v) + " " + str(dimacs(l)) + " " + str((bvID)) + "\n" );
        return l
    
    def maximumFlow_geq(self, graph, s,t,flow):
        self.checkNode(graph,s);
        self.checkNode(graph,t);
        self.backtrack()
        l= self.monosat_c.maximumFlow_geq(self.solver._ptr,graph,(s),(t),(flow))
        if self.solver.output:  
            self._echoOutput("maximum_flow_geq " + str(self.getGID(graph)) + " " + str(s) + " " + str(t) + " " + str(dimacs(l)) + " " + str((flow)) + "\n" );
        return l
    
    def maximumFlow_gt(self, graph, s,t,flow):
        self.checkNode(graph,s);
        self.checkNode(graph,t);
        self.backtrack()
        l= self.monosat_c.maximumFlow_gt(self.solver._ptr,graph,(s),(t),(flow))
        if self.solver.output: 
            self._echoOutput("maximum_flow_gt " + str(self.getGID(graph)) + " " + str(s) + " " + str(t) + " " + str(dimacs(l)) + " " + str((flow)) + "\n" );
        return l
               
    def maximumFlow_geq_bv(self, graph, s,t,bvID):
        self.checkNode(graph,s);
        self.checkNode(graph,t);
        self.checkBV(bvID)
        self.backtrack()
        l= self.monosat_c.maximumFlow_geq_bv(self.solver._ptr,graph,(s),(t),(bvID))
        if self.solver.output:
            self._echoOutput("maximum_flow_bv_geq " + str(self.getGID(graph)) + " " + str(s) + " " + str(t) + " " + str(dimacs(l)) + " " + str((bvID)) + "\n" );
        return l
    
    def maximumFlow_gt_bv(self, graph, s,t,bvID):
        self.checkNode(graph,s);
        self.checkNode(graph,t);
        self.checkBV(bvID)
        self.backtrack()
        l= self.monosat_c.maximumFlow_gt_bv(self.solver._ptr,graph,(s),(t),(bvID))
        if self.solver.output:
            self._echoOutput("maximum_flow_bv_gt " + str(self.getGID(graph)) + " " + str(s) + " " + str(t) + " " + str(dimacs(l)) + " " + str((bvID)) + "\n" );
        return l
    
    def minimumSpanningTree_leq(self, graph,weight):
        self.backtrack()
        l= self.monosat_c.minimumSpanningTree_leq(self.solver._ptr,graph,(weight))
        if self.solver.output:      
            self._echoOutput("mst_weight_leq " + str(self.getGID(graph)) + " " +  str(dimacs(l)) + " " + str((weight)) + "\n" );
        return l
            
    def minimumSpanningTree_lt(self, graph,weight):
        self.backtrack()
        l= self.monosat_c.minimumSpanningTree_lt(self.solver._ptr,graph,(weight))
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
        self.checkLit(lit);
        return self.monosat_c.getModel_Literal(self.solver._ptr, lit);

    def getModel_BV(self, bvID,getMaximumValue=False):
        self.checkBV(bvID)
        return self.monosat_c.getModel_BV(self.solver._ptr, self.solver.bvtheory,(bvID),(getMaximumValue));
        
    def getModel_MaxFlow(self, graph, flowlit):
        self.checkLit(flowlit);
        return self.monosat_c.getModel_MaxFlow(self.solver._ptr, graph,flowlit);
        
    def getModel_EdgeFlow(self, graph, flowlit,edgelit):
        self.checkLit(flowlit);
        self.checkLit(edgelit);
        return self.monosat_c.getModel_EdgeFlow(self.solver._ptr, graph,flowlit, edgelit);

    def getModel_AcyclicEdgeFlow(self, graph, flowlit,edgelit):
        self.checkLit(flowlit);
        self.checkLit(edgelit);
        return self.monosat_c.getModel_AcyclicEdgeFlow(self.solver._ptr, graph,flowlit, edgelit);


    def getModel_MinimumSpanningTreeWeight(self, graph, mstlit):
        self.checkLit(mstlit)
        return self.monosat_c.getModel_MinimumSpanningTreeWeight(self.solver._ptr, graph,mstlit); 

    def getModel_Path_Nodes(self, graph, reach_or_distance_lit):
        self.checkLit(reach_or_distance_lit)
        arg_length = self.monosat_c.getModel_Path_Nodes_Length(self.solver._ptr, graph,reach_or_distance_lit);
        if (arg_length <0):
            return None
        elif arg_length == 0:
            return []
        path = []
        #path = list(range(arg_length))
        #path_pointer = self.getIntArray(path)
        l = self.monosat_c.getModel_Path_Nodes(self.solver._ptr, graph,reach_or_distance_lit, arg_length,path);
        if l != arg_length:
            raise RuntimeError("Error reading path model")
        return path
        #return self.intArrayToList(path_pointer,arg_length)

    def getModel_Path_EdgeLits(self, graph, reach_or_distance_lit):
        self.checkLit(reach_or_distance_lit)
        arg_length = self.monosat_c.getModel_Path_EdgeLits_Length(self.solver._ptr, graph,reach_or_distance_lit);
        if (arg_length <0):
            return None
        elif arg_length == 0:
            return []
        path = []
        #path = list(range(arg_length))
        #path_pointer = self.getIntArray(path)
        l = self.monosat_c.getModel_Path_EdgeLits(self.solver._ptr, graph,reach_or_distance_lit, arg_length,path);
        if l != arg_length:
            raise RuntimeError("Error reading path model")
        return path
        #return self.intArrayToList(path_pointer,arg_length)

    def checkLit(self,l):
        if (l<0):
            raise RuntimeError("Bad literal %d"%(l))
        v = l//2
        if(v>= self.nVars()):
            raise RuntimeError("Bad variable %d (largest variable is %d)"%(v, self.nVars()-1))


    def checkBV(self,bvID):
        if(bvID<0 or bvID>=self.nBitvectors()):
            raise RuntimeError("Bitvector %d does not exist"%(bvID))

    def checkNode(self,g,s):

        if s<0 or s>=self.nNodes(g):
            raise RuntimeError("Node %d does not exist in this graph (the largest node is %d)"%(s,self.nNodes(g)-1))

        
