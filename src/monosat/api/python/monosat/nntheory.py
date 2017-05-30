#The MIT License (MIT)
#
#Copyright (c) 2015, Sam Bayless
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


import monosat.monosat_c
from monosat.logic import *
from monosat.bvtheory import BitVector
from monosat.manager import Manager
import sys
debug=False   

#Collects a set of graphs to encode together into a formula
class NNManager(metaclass=Manager):
    
 
    
    def  __init__(self):
        self.nns=[]        
    
    def clear(self):
        self.nns=[]
    
    def addNN(self, g):
        self.nns.append(g)
    
    def getNN(self,gid):
        return self.nns[gid]
    


class NeuralNetwork():

    def __init__(self,prototxt, caffemodel):
        self._monosat = monosat.monosat_c.Monosat()
        manager = NNManager()
        manager.addNN(self)
        self.nn = self._monosat.newNeuralNetwork(prototxt,caffemodel)
        
        self.id = len(manager.nns)
        
    def getInput(self, index, bvwidth=8, min=0.0,max=1.0):  
        bv = BitVector(bvwidth,"anon")
        self._monosat.newNN_BV(self.nn,bv.getID(),True,index,min,max)        
        return bv
    
    def getOutput(self, index, bvwidth=8, min=0.0,max=1.0):  
        bv = BitVector(bvwidth,"anon")
        self._monosat.newNN_BV(self.nn,bv.getID(),False,index,min,max)        
        return bv
    
    def inputSize(self):
        return self._monosat.nn_layer_size(self.nn,True)
        
    def outputSize(self):
        return self._monosat.nn_layer_size(self.nn,False)
        