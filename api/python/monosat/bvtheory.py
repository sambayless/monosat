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

import monosat.monosat_c
from monosat.logic import *
from monosat.singleton import Singleton

import sys
debug=False  
#Collects a set of graphs to encode together into a formula
class BVManager(metaclass=Singleton):
        
    def  __init__(self):
        self.bvs = []
        self.aux_bvs=[]
        self.comparisons=[]
        self.consts=dict()
        
        self._monosat = monosat.monosat_c.Monosat()
        self._monosat.initBVTheory()
        self.bitblast_addition=False
        self.bitblast_addition_shadow=False
    
    def _setConstant(self,val,width,bv):
        self.consts[(val,width)]=bv
    
    def _hasConstant(self,val,width):
        return (val,width) in self.consts
    
    def _getConstant(self,val,width):
        return self.consts[(val,width)]
    
    def clear(self):
        self.bvs = []
        self.aux_bvs=[]
        self.comparisons=[]        
    
    #Each "string" must actually be a list of positive integers
    def bv(self,  width=8, const_value=None):

        return BitVector(self,width,const_value)
  
    def hasBitVector(self,id):
        return id>=0 and id< len(self.bvs) and self.bvs[id] is not None
    
    
    def getBitVector(self, id):
        return self.bvs[id]
    
    def getBitVectors(self):
        return self.bvs
    
    def write(self,f):
  
        
        for bv in self.bvs:
            bv.write(f)
        for i,bv in enumerate(self.aux_bvs):
            bv.pid = len(self.bvs)+i
            bv.write(f)
            
_bv_manager = BVManager() 
def bv(width, const_value=None):
    return _bv_manager.bv(width,const_value)                   
       
class BitVector():    
    def __init__(self,mgr,width,op=None,args=None):
        
        self.mgr = mgr
        self._bv=None
        self.symbol=None   

        self._width=width

        
        self.pid = None
              
        if args is None and isinstance(op, (int, float)):
            val = int(op)
            originalval=val

            if val<0:
                val=0
                print("Warning: negative bitvectors not yet supported, setting to 0", file=sys.stderr)
            
            if  val>= (1<<width):
                val = (1<<width)-1
                print("Warning: value %d is too large to represent with a width-%d bitvector, setting to %d"%(originalval,width, val), file=sys.stderr)

            self._constant=val


            #if(val>0):
            #    print("Warning: value %d is too large to represent with a width-%d bitvector, setting to %d"%(originalval,width, (1<<width)-1), file=sys.stderr)
            op=None  
            args=None   
            if not mgr._hasConstant(val,width):
                self.pid = mgr._monosat.newBitvector_const(width, val)
                mgr._setConstant(val,width,self)
            else:
                self.pid = mgr._getConstant(val,width).pid
            self._bv=[None]*width
            #fill bv with constants, for convenience elsewhere
            for i in range(width-1,-1,-1):
                v = 1<<i
                if val>=v:
                    val-=v
                    self._bv[i]=true 
                else:
                    self._bv[i]=false            
            
        else:
            self._bv=[] 
            for _ in range(width):
                self._bv.append(Var())   
            #arr = (c_int*width)()
            #for i,v in enumerate(self._bv):
            #    arr[i]=c_int(v.getLit()//2)
            bits = [v.getLit()//2 for v in self._bv]
                     
            self.pid = mgr._monosat.newBitvector(bits)

        if op == '+':
            assert(len(args)==2)
            if not mgr.bitblast_addition:
                mgr._monosat.bv_addition(args[0].getID(), args[1].getID(), self.getID())
            if mgr.bitblast_addition or mgr.bitblast_addition_shadow:
                carry = false
                for i, (a,b,out) in enumerate(zip(args[0],args[1],self)):
                    r,carry2=Add(a,b,carry)
                    Assert(out==r)
                    carry=carry2
                Assert(Not(carry))#disallow overflow.
    
    def checkValue(self,val):
        if val<0 or val>= 1<<self.width():
            print("Error: value %d is too large to represent with a width-%d bitvector"%(val,self.width()), file=sys.stderr)
            assert(False)
        
            
    def value(self):
        return self.mgr._monosat.getModel_BV(self.pid)
    
   
    
    def getSymbol(self):
        return self.symbol
    
    def setSymbol(self,s):        
        self.symbol = str(s)        

    """                
    def __eq__(self,other):
        if other and isinstance(other,self.__class__):
            return other.getID()==self.getID()
        return False"""
    
    def __hash__(self):
        return self.getID()
    
    def width(self):
        return self._width
    
    def getID(self):
        return self.pid
    
    def __getitem__(self,index):
        return self._bv[index]
    
    def __len__(self):
        return self._width
    
    def __add__(self,other):
        if not isinstance(other, BitVector):
            other = BitVector(self.mgr,self.width(),other)
        return BitVector(self.mgr,self.width(),'+',(self,other))
    
    __radd__ = __add__    
    
    def lt(self,compareTo):
        if  isinstance(compareTo, BitVector):
            return Var(self.mgr._monosat.newBVComparison_bv_lt(self.getID(),compareTo.getID()))
        else:
            self.checkValue(int(compareTo))
            return Var(self.mgr._monosat.newBVComparison_const_lt(self.getID(),int(compareTo)))

    def leq(self,compareTo):
        if  isinstance(compareTo, BitVector):
            return Var(self.mgr._monosat.newBVComparison_bv_leq(self.getID(),compareTo.getID()))
        else:
            self.checkValue(int(compareTo))
            return Var(self.mgr._monosat.newBVComparison_const_leq(self.getID(),int(compareTo)))

    
    def gt(self,compareTo):
        if  isinstance(compareTo, BitVector):
            return Var(self.mgr._monosat.newBVComparison_bv_gt(self.getID(),compareTo.getID()))
        else:
            self.checkValue(int(compareTo))
            return Var(self.mgr._monosat.newBVComparison_const_gt(self.getID(),int(compareTo)))
        #    compareTo = BitVector(self.mgr,self.width(),compareTo)
        
    
    def geq(self,compareTo):
        if  isinstance(compareTo, BitVector):
            return Var(self.mgr._monosat.newBVComparison_bv_geq(self.getID(),compareTo.getID()))
        else:
            self.checkValue(int(compareTo))
            return Var(self.mgr._monosat.newBVComparison_const_geq(self.getID(),int(compareTo)))
    
    def eq(self,compareTo):
        #if not isinstance(compareTo, BitVector):
        #    compareTo = BitVector(self.mgr,self.width(),compareTo)
        return And(self.leq(compareTo),self.geq(compareTo))
    
    def neq(self,compareTo):
        #if not isinstance(compareTo, BitVector):
        #    compareTo = BitVector(self.mgr,self.width(),compareTo)
        return Nand(self.leq(compareTo),self.geq(compareTo))
    
    def __lt__(self,other):
        return self.lt(other)

    def __le__(self,other):
        return self.leq(other)
    
    def __gt__(self,other):
        return self.gt(other)
    
    def __ge__(self,other):
        return self.geq(other)    
  
    def __eq__(self,other):
        return self.eq(other)     

    def __ne__(self,other):
        return self.neq(other)     
    
    def write(self,f):
        
        if self._constant is not None:
            f.write("bv const %d %d "%(self.pid, self.width()) + str(self._constant) +"\n" )   
             
        else:        
            f.write("bv %d %d"%(self.pid, self.width()))   
            for var in self._bv:
                f.write(" %d"%(var.getInputLiteral()))    
            f.write("\n")        
        if self.symbol is not None:
            f.write("c bv %d "%(self.pid) + str(self.symbol) +"\n" )
            
        if self.op == '+':
            f.write("bv + %d"%(self.pid))
            assert(len(self.args)==2)
            for arg in self.args:
                #if isinstance(arg,BitVector):
                f.write(" %d"%(arg.getID()))
                #else:
                #    f. write(" " +str(arg))
            f.write("\n")

        for (v,compareTo,op) in self._comparisons:
            #if (isinstance(compareTo,BitVector)):                       
            f.write("bv " + str(op) + " %d %d  "%(v.getInputLiteral(), self.pid) + " "+ str(compareTo.getID()) + "\n")         
            #else:
            #    f.write("bv " + str(op) + " %d bv%d  "%(v.getInputLiteral(), self.pid) + str(compareTo) + "\n") 