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
from monosat.singleton import Singleton
class GeometryManager(metaclass=Singleton):
        
    def  __init__(self):
        self.pointsets = []
        self.polygons=[]
        
    def clear(self):
        self.pointsets = []
        self.polygons=[]
    
    def addPointSet(self, g):
        self.pointsets.append(g)
    
    def newPolygon(self,polygon):
        self.polygons.append(polygon)
        return len(self.polygons)-1
    
    def write(self,f):
  
        
        for polygon in self.polygons:
            polygon.write(f)
        
        for pid in range(len(self.pointsets)):
            pointset = self.pointsets[pid]
            pointset.id = pid;        
        for pid in range(len(self.pointsets)):
            pointset = self.pointsets[pid]
            
            for point,v in pointset.points:         
                f.write("point %d %d %d "%(pid,v.getInputLiteral(), len(point)))
                for p in point:
                    f.write(str(p) + " ")
                f.write("\n")
                
        for pid in range(len(self.pointsets)):
            pointset = self.pointsets[pid]

            
            for (area,v) in pointset.constraints_convexHullAreas:
                f.write("c convex hull area < " + str(area) + "\n") 
                f.write("convex_hull_area_gt " + str(pid) + " " + str(v.getInputLiteral()) + " " + str(area)  + "\n")
            
            """   
            for (point,v) in pointset.constraints_convexHullContains:
                f.write("c convex hull contains point " + str(point) + "\n") 
                pstr = " "
                for p in point:
                    pstr += str(p) + " "                
                f.write("convex_hull_containment " + str(pid)  + " " + str(v.getInputLiteral()) + " " + str(len(point)) + " " +pstr + "\n")
            """
            for (points,v,inclusive) in pointset.constraints_convexHullIntersects:
                if not inclusive:
                    f.write("c convex hull intersects polygon " + str(points) + "\n") 
                    pstr = " "
                    dim = len(points[0])
                    
                    for p in points:
                        for d in p:
                            pstr += " " + str(d)               
                    f.write("convex_hull_intersects_polygon " + str(pid)  + " " + str(v.getInputLiteral()) + " " + str(len(points)) + " " + str(dim) +pstr + "\n")
                else:
                    f.write("c convex hull collides with polygon " + str(points) + "\n") 
                    pstr = " "
                    dim = len(points[0])
                    
                    for p in points:
                        for d in p:
                            pstr += " " + str(d)               
                    f.write("convex_hull_collides_polygon " + str(pid)  + " " + str(v.getInputLiteral()) + " " + str(len(points)) + " " + str(dim) +pstr + "\n")
                 
             
            for (point,pointVar,var) in pointset.constraints_pointOnHullOrDisabled:
                """
                pointVar = None
                for p,v in pointset.points:
                    if (p==point):
                        pointVar = v;
                        break
                """
                if pointVar is None:
                    f.write("c point is an element of convex hull constraint set to false because point is not in the pointset\n") 
                    f.write("-" + str(v.getInputLiteral()) + " 0\n")
                else:
                    f.write("c point is an element of convex hull " + str(pointVar.getInputLiteral()) + "\n") 
                    f.write("point_on_convex_hull " + str(pid) + " " + str(pointVar.getInputLiteral()) + " " + str(var.getInputLiteral()) + "\n")
                    
                
            for hull,v,inclusive in pointset.constraints_convexHullIntersectsHull:
                f.write("c two convex hulls intersecting each other\n")
                if inclusive:
                    f.write("convex_hulls_collide "+ str(pid) + " " + str(hull.id) + " "+ str(v.getInputLiteral()) + "\n")
                else:
                    f.write("convex_hulls_intersect "+ str(pid) + " " + str(hull.id) + " "+ str(v.getInputLiteral()) + "\n")
         

class SymbolicPolygon():    
    def __init__(self,mgr):
        self.pid = mgr.newPolygon(self) 
        self.area_lt_lits =[]
        self.overlap_lits=[]
        self.contains_lits = []
    
    def getID(self):
        return self.pid
    
    def areaLT(self,area):
        v = Var()
        self.area_lt_lits.append((v,area,False))
        return v
        
    def areaLEQ(self,area):
        v = Var()
        self.area_leq_lits.append((v,area,True))
        return v    
   
    def overlaps(self,polygon,inclusive=True):
        v = Var()
        self.overlap_lits.append((v,polygon,inclusive))
        return v    
    
    def contains(self,point,inclusive=True):
        v = Var();
        self.contains_lits.append((v,point,inclusive))
        return v


    def write(self,f):
        for (v,area,inclusive) in self.area_lt_lits:
            if inclusive:
                f.write("area_leq %d %d %d\n"%( self.pid, area, v.getInputLiteral()))
            else:
                f.write("area_lt %d%d %d\n"%( self.pid, area, v.getInputLiteral()))
            
        for (v,polygon,inclusive) in self.overlap_lits:
            if inclusive:
                f.write("overlaps %d %d %d\n"%(self.pid, polygon.getID(), v.getInputLiteral()))
            else:
                f.write("intersects %d %d %d\n"%(self.pid,  polygon.getID(), v.getInputLiteral()))         

        for (v,point,inclusive) in self.contains_lits:
            if inclusive:
                f.write("contains %d %d %d %d %d\n"%(self.pid, v.getInputLiteral(),len(point), point[0],point[1]))
            else:
                f.write("contains_exclusive %d %d %d %d\n"%(self.pid, v.getInputLiteral(),len(point), point[0],point[1]))    
           
"""    
class ConstantPolygon(SymbolicPolygon):
    def __init__(self,polygon):
        self.polygon=polygon
"""
       
class ConstantRect(SymbolicPolygon):
    def __init__(self,mgr,x=0,y=0,width=1,height=1):
        super().__init__(mgr)
        self.width=width
        self.height=height
        self.x = x
        self.y =y
        
    def write(self,f):
        f.write("rectangle %d %d %d %d %d %d\n"%(self.pid,2,self.x,self.y,self.width,self.height))
        super().write(f)        

class BinaryOperationPolygon(SymbolicPolygon):
    class BinOp(Enum):
        union=1
        intersection=2
        difference=3
        minkowskisum=4
    
    def __init__(self,mgr,r1,r2,operation):
        super().__init__(mgr)
        self.r1 = r1
        self.r2 = r2
        self.operation=operation
        
    def write(self,f):

        f.write("polygon_operation %d %d %d %s\n"%(self.pid,self.r1.getID(),self.r2.getID(),str(self.operation.name)))
        super().write(f)    

class ConditionalPolygon(SymbolicPolygon):
    def __init__(self,mgr,p1,p2,v):
        super().__init__(mgr)
        self.p1 = p1
        self.p2 = p2
        self.v =v
        
    def write(self,f):
        f.write("conditional_polygon %d %d %d %d\n"%(self.pid,self.v.getInputLiteral(),self.p1.getID() ,self.p2.getID()))
        super().write(f)              

class ConvexHull(SymbolicPolygon):
    def __init__(self,pointset):
        self.pointset = pointset

class PointSet():
    
    def __init__(self, manager):
        self.npoints=0
        self.numedges=0
        self.points=[]
        self.pointvars=[]    
        self.dimension=-1 
        self.constraints_convexHullAreas=[]
        #self.constraints_convexHullContains=[]
        self.constraints_convexHullIntersects=[]
        self.constraints_pointOnHullOrDisabled=[]
        self.constraints_convexHullIntersectsHull=[]
        self.convexHull=ConvexHull(self)
        
        manager.addPointSet(self)
    
    def addPoint(self, point):
        v = Var()
        if(self.dimension<0):
            self.dimension=len(point)
        if(len(point)!=self.dimension):
            assert(False)
        self.points.append((point,v))   
        self.pointvars.append(v)
        return v
    
    def getPoints(self):
        return self.pointvars
    
    def pointEnabled(self,index):
        return self.pointvars[index]
    
    def convexHullAreaGT(self,gt):
        v = Var()
        self.constraints_convexHullAreas.append((gt,v))
        return v

    def convexHullContains(self,point, inclusive=True):
        return self.convexHullIntersects([point],inclusive)
        #v = Var()
        #self.constraints_convexHullContains.append((point,v))
        #return v

    def convexHullIntersects(self,points, inclusive=True):
        v = Var()
        self.constraints_convexHullIntersects.append((points,v,inclusive))
        return v

    def convexHullIntersectsHull(self,hull, inclusive=True):
        v = Var()
        self.constraints_convexHullIntersectsHull.append((hull,v,inclusive))
        return v

    def getConvexHull(self):
        return self.convexHull

    def pointOnHull(self,point):
        pointVar = None
        for p,v in self.points:
            if (p==point):
                pointVar = v;
                break
        if pointVar is None:
            assert(False)
        v = Var()
                
        self.constraints_pointOnHullOrDisabled.append((point,pointVar,v))
        return And(pointVar, v)
   
