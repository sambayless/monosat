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
debug=False   
#Collects a set of graphs to encode together into a formula
class LSystemManager():
        
    def  __init__(self):
        self.strings = []
        self.systems=[]
    
    #Each "string" must actually be a list of positive integers
    def addString(self, string):
        ints = []
        for s in string:
            if(isinstance( s, int )):
                ints.append(s)
            else:
                ints.append(ord(s)-ord('a'))
            
        self.strings.append(ints)
        return len(self.strings)-1
    

    
    def newLSystem(self, fsm):        
        self.systems.append(fsm);
        return len(self.systems)-1
    
    def write(self,f):
  
        
        for fsm in self.systems:
            fsm.write(f)
        
        for id, string in enumerate(self.strings):
            f.write("lstr %d"%(id))
            for i in string:
                assert(i>0)                
                f.write(" %d"%(i))
            f.write("\n")       
                   
       
class LSystem():    
    def __init__(self,mgr):
        self.pid = mgr.newLSystem(self)
     
        self.n_characters=0
       
        self._rules=[]
        
        self._produces =[]

        
    
    def getID(self):
        return self.pid
    
    def addCharacter(self):
        self.n_characters+=1
    
    def addRule(self,c, rule):
        var = Var()
        self._rules.append((c,var,rule))
        return var;
 
    def getRules(self):
        return self._rules;

    def produces(self,atom,stringID):
        v = Var()
        self._produces.append((v,atom,stringID))
        return v
    
    
    
    def write(self,f):
        f.write("lsystem %d %d %d\n"%(self.pid, self.n_characters, len(self._rules)))       
        
        for (c,var,rule) in self._rules:
            f.write("production %d %d %d"%( self.pid,c, var.getInputLiteral()))
            for c in rule:
                f.write(" %d"%(c))
            f.write("\n")
        
        for (v,atom,stringID) in self._produces:
            f.write("produces %d %d %d %d\n"%( self.pid,atom, stringID, v.getInputLiteral()))

