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

from enum import Enum
from monosat.logic import *
from monosat.manager import Manager
debug=False   
#Collects a set of graphs to encode together into a formula
class FSMManager(metaclass=Manager):
        
    def  __init__(self):
        self.strings = []
        self.fsms=[]
    
    def clear(self):
        self.strings = []
        self.fsms=[]        
    
    #Each "string" must actually be a list of positive integers
    def addString(self, string):
        ints = []
        for s in string:
            if(isinstance( s, int )):
                ints.append(s+1)
            else:
                ints.append(ord(s)-ord('a')+1)
            
        self.strings.append(ints)
        return len(self.strings)-1
    

    
    def newFSM(self, fsm):        
        self.fsms.append(fsm);
        return len(self.fsms)-1
    
    def write(self,f):
  
        
        for fsm in self.fsms:
            fsm.write(f)
        
        for id, string in enumerate(self.strings):
            f.write("str %d"%(id))
            for i in string:
                assert(i>0)                
                f.write(" %d"%(i))
            f.write("\n")       
                   
       
class FSM():    
    def __init__(self,mgr):
        self.pid = mgr.newFSM(self)
        self.hasEpsilon=True 
        self.n_labels=1
        self.n_out_labels=1
        self.n_states=0
        
        self._transitions=[]
        
        self._accepts =[]
        self._generates=[]
        self._transduces=[]
        self._accepts_generator=[]
        
    
    def getID(self):
        return self.pid

    def addInLabel(self):
        self.n_labels+=1
        return self.n_labels-1
        
    def addLabel(self):
        return self.addInLabel()

    def addOutLabel(self):
        self.n_out_labels+=1
        return self.n_out_labels-1
    
    def nInLabels(self):
        return self.n_labels
    
    def nOutLabels(self):
        return self.n_out_labels 
    
    def addState(self):
        self.n_states+=1
        return self.n_states-1
    
    def nStates(self):
        return self.n_states;
    def getTransitions(self):
        return self._transitions
    def addTransition(self, u,v, label, output=0):
        assert(u<self.n_states)
        assert(u>=0)
        assert(v<self.n_states)
        assert(v>=0)
        var = Var()
        self._transitions.append((u,v,label,output,var))
        return var
    
    def accepts(self,stringID,source,dest):
        v = Var()
        self._accepts.append((v,stringID,source,dest))
        return v
 
    def generates(self,stringID,source):
        v = Var()
        self._generates.append((v,stringID,source))
        return v
    
    def transduces(self,stringID1,stringID2,source,dest):
        v = Var()
        self._transduces.append((v,stringID1,stringID2,source,dest))
        return v  
    
    def accepts_generator(self,generator,generator_source, generator_dest,acceptor_source,acceptor_dest):
        v = Var()
        self._accepts_generator.append((v,generator,generator_source,generator_dest,acceptor_source,acceptor_dest,-1))
        return v           
    
    def write(self,f):
        f.write("fsm %d %d %d\n"%(self.pid, self.n_labels, 1 if self.hasEpsilon else 0))       
        
        for (u,v,label,output,var) in self._transitions:
            f.write("transition %d %d %d %d %d %d\n"%( self.pid,u,v,label,output, var.getInputLiteral()))
        
        for (v,stringID,source,dest) in self._accepts:
            f.write("accepts %d %d %d %d %d\n"%( self.pid,source,dest, stringID, v.getInputLiteral()))

        for (v,stringID,source) in self._generates:
            f.write("generates %d %d %d %d\n"%( self.pid,source, stringID, v.getInputLiteral()))
           
        for (v, stringID1,stringID2,source,dest) in self._transduces:
            f.write("transduces %d %d %d %d %d %d\n"%( self.pid,source,dest, stringID1,stringID2, v.getInputLiteral()))
            
        for (v,generator,generator_source,generator_dest,acceptor_source,acceptor_dest, stringID) in  self._accepts_generator:
            f.write("composition_accepts %d %d %d %d %d %d %d %d\n"%( generator.pid,self.pid,generator_source,generator_dest,acceptor_source,acceptor_dest,stringID, v.getInputLiteral()))
            
    def draw(self):
        
        print("digraph{")
        for (u,v,label,output,var) in self._transitions:
            print("n%d -> n%d [label=\"%d:%d\"]"%( u,v,label,output))

        
        print("}")
        
             