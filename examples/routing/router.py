#The MIT License (MIT)
#
#Copyright (c) 2017, Sam Bayless
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


from monosat import *
from time import time
import pcrt
import itertools

#There are many ways to perform circuit routing using MonoSAT.
#The approach uses just one graph, and uses a combination of MonoSAT's built-in reachability constraints to ensure
#nets are routed, while using propositional constraints over the edges in the graph to prevent nets from intersecting.
#This function also supports a slightly more complex router, which combines reachability and maximum flow constraints.
#The maximum flow constraint is not powerful enough to ensure correct routing, but is a (safe) overapproximative
#constraint, that allows the solver to prune large parts of the search space early.
#
#The variation described here supports on 2-terminal routing; use router_multi.py for multi-terminal routing.
#Finally, for the special case of Escape Routing (in which the destinations are interchangeable), see
#
#Bayless, Sam, Holger H. Hoos, and Alan J. Hu. "Scalable, high-quality, SAT-based multi-layer escape routing."
#Computer-Aided Design (ICCAD), 2016 IEEE/ACM International Conference on. IEEE, 2016.
def route(filename, monosat_args,use_maxflow=False, draw_solution=True, outputFile = None):
    (width, height),diagonals,nets,constraints,disabled = pcrt.read(filename)
    print(filename)
    print("Width = %d, Height = %d, %d nets, %d constraints, %d disabled"%(width,height,len(nets),len(constraints), len(disabled)))
    if diagonals:
        print("45-degree routing")
    else:
        print("90-degree routing")

    for net in nets:
        if(len(net)!=2):
            raise Exception("router.py  only supports routing nets with exactly 2 vertices. Use router_multi.py for routing nets with 2+ vertices.")

    if(len(monosat_args)>0):
        args = " ".join(monosat_args)
        print("MonoSAT args: " + args)
        Monosat().newSolver(args)

    if outputFile is not None:
        print("Writing output to " + outputFile)
        Monosat().setOutputFile(outputFile)

    g = Graph()

    print("Building grid")
    edges = dict()
    grid = dict()

    for x in range(width):
        for y in range(height):
            n = g.addNode("%d_%d"%(x,y))
            grid[(x,y)] = n
            edges[(x,y)]=[]


    disabled_nodes = set(disabled)
    undirected_edges = dict()
    #create undirected edges between neighbouring nodes
    if draw_solution:
        #the dicts here are only used to print the solution at the end. Disable for benchmarking.
        lefts = dict()
        rights = dict()
        ups = dict()
        downs = dict()

        down_lefts = dict()
        up_rights = dict()
        down_rights = dict()
        up_lefts = dict()



    start_nodes = set()
    end_nodes = set()
    net_nodes = set()
    for net in nets:
        assert(len(net)==2)
        start_nodes.add(net[0])
        end_nodes.add(net[1])
        for (x,y) in net:
            net_nodes.add((x,y))


    #create undirected edges between neighbouring nodes
    def makeEdge(n,r):
        e = None
        if n not in disabled_nodes and r not in disabled_nodes:
            if n in net_nodes or r in net_nodes:
                #only add directed edges for start/end nodes.
                allow_out=True
                allow_in = True
                if n in start_nodes or r in end_nodes:
                    allow_in=False
                if n in end_nodes or r in start_nodes:
                    allow_out = False
                assert(not (allow_in and allow_out))
                if allow_out:
                    e = g.addEdge(grid[n], grid[r]) #add a _directed_ edge from n to r
                    undirected_edges[e]=e
                elif allow_in:
                    e = g.addEdge(grid[r], grid[n])  # add a _directed_ edge from r to n
                    undirected_edges[e]=e
                else:
                    e=None
            else:
                #add undirected edges for other nodes in the grid
                e = g.addEdge(grid[n],grid[r]) #in monosat, undirected edges are specified by creating
                undirected_edges[e]=e
                e2 = g.addEdge(grid[r],grid[n]) #directed edges in both directions, and making them equal
                undirected_edges[e2]=e
                AssertEq(e,e2)
            if e is not None:
                edges[n].append(e)
                edges[r].append(e)
        return e

    for x in range(width):
        for y in range(height):
            n = (x,y)
            if n in disabled_nodes:
                continue

            if x<width-1:
                r = (x+1,y)
                e = makeEdge(n,r)
                if e is not None and draw_solution:
                    rights[n] = e
                    lefts[r] = e

            if y<height-1:
                r = (x,y+1)
                e = makeEdge(n,r)
                if e is not None and draw_solution:
                    downs[n] = e
                    ups[r] = e

            if diagonals:
                #if 45 degree routing is enabled
                diag_up=None
                diag_down=None
                if x<width-1 and y<height-1:
                    r = (x+1,y+1)
                    e = makeEdge(n,r)
                    diag_down = e
                    if e is not None and draw_solution:
                        down_rights[n] = e
                        up_lefts[r] = e

                if x>0 and y<height-1 and False:
                    r = (x-1,y+1)
                    e = makeEdge(n,r)
                    diag_up = e
                    if e is not None and draw_solution:
                        down_lefts[n] = e
                        up_rights[r] = e
                if diag_up and diag_down:
                    AssertNand(diag_up, diag_down) #cannot route both diagonals

    if len(constraints)>0:
        print("Enforcing constraints")
        vertex_used = dict()
        for x in range(width):
            for y in range(height):
                vertex_used[(x, y)] = Or(edges[(x, y)])  # A vertex is used exactly if one of its edges is enabled

        for constraint in constraints:

            if len(constraint)>1:
                #a constraint is a list of vertices of which at most one can be used
                vertex_used_list = []
                for node in constraint:
                    vertex_used_list.append(vertex_used[node])
                if(len(vertex_used_list)<20):
                    for a,b in itertools.combinations(vertex_used_list, 2):
                        AssertOr(~a,~b) #in every distinct pair of edges, at least one must be false
                else:
                    AssertAtMostOne(vertex_used_list) #use more expensive, but more scalable, built-in AMO theory

    print("Enforcing net separation")
    #enforce that at any node that is not a starting/ending node, either exactly 2, or exactly no, edges are enabled
    #this approach ensures acyclicity, and also ensures that nodes from one net do not reach any other net.
    #However, this approach cannot be used to route trees.
    for x in range(width):
        for y in range(height):
            n = (x,y)

            edge_list = edges[n]
            if n in net_nodes:
                #n is a start or end node in the net.
                #no constraint is required, as only directed edges are available at start/end nodes
                pass
            else:
                #n is not a start or end node; either exactly 2, or exactly 0, adjacent edges must be enabled
                #There are many ways to actually enforce that constraint, not clear what the best option is
                #This is slightly complicated by edge nodes/corner nodes having fewer neighbours.


                if len(edge_list)>2:
                    AssertNand(edge_list) #At least one edge must be disabled

                if len(edge_list)>=4:
                    #All but two edges must be disabled.
                    #This uses up to ~28 constraints. It might be better to implement this using PB constraints instead.
                    disabled_two = false()
                    for pair in itertools.combinations(edge_list, len(edge_list)-2):
                        disabled_two = Or(disabled_two, Nor(pair))
                    Assert(disabled_two)




    print("Enforcing net routing")
    #enforce reachability using MonoSAT's builtin reachability constraints
    #notice that unreachability between different nets does not need to be explicitly enforced here, because of the
    #edge constraints above.
    reach_constraints = []
    for net in nets:
        assert(len(net)==2)
        n1 = net[0]
        n2 = net[1]
        r = g.reaches(grid[n1],grid[n2])
        reach_constraints.append(r)
        Assert(r)
        r.setDecisionPriority(1); #decide reachability before considering regular variable decisions

    #optional: use a maximum flow constraint, as an overapprox of disjoint reachability.
    if use_maxflow:
        print("Enforcing maxflow constraint")
        #add a source and dest node, with 1 capacity from source to each net start vertex, and 1 capacity from each net end vertex to dest
        source = g.addNode()
        dest = g.addNode()
        for net in nets:
            Assert(g.addEdge(source,grid[net[0]],1)) #directed edges!
            Assert(g.addEdge(grid[net[1]],dest, 1))  #directed edges!
        m = g.maxFlowGreaterOrEqualTo(source,dest,len(nets))
        Assert(m)
        m.setDecisionPriority(-1);  # never make decisions on the maxflow predicate.

    if outputFile is not None:
        print("Wrote constraints to " + outputFile + ", exiting without solving")
        sys.exit(0)

    print("Solving...")
    if Solve():
        print("Solved!")
        if draw_solution:
            paths = set() #collect all the edges that make up the routing
            for r in reach_constraints:
                path = g.getPath(r,True)
                for e in path:
                    assert(e.value())
                    if e in undirected_edges.keys():
                        e = undirected_edges[e] #normalize the edge var
                        paths.add(e)

            #print the solution to the console
            for y in range(height):
                for x in range(width):
                    n = (x,y)

                    if (x,y) in net_nodes:
                        print("*",end="")
                    else:
                        print(" ",end="")
                    if x<width-1:
                        r = (x+1,y)
                        if n in rights:
                            e = rights[n]
                            if e in paths:
                                print("-",end="")
                            else:
                                print(" ",end="")
                        else:
                            print(" ", end="")
                print()
                for x in range(width):
                    n = (x, y)
                    if y<height-1:
                        r = (x,y+1)
                        if n in downs:
                            e = downs[n]
                            if e in paths:
                                print("|",end="")
                            else:
                                print(" ",end="")
                        else:
                            print(" ", end="")
                    drew_diag=False
                    if diagonals:
                        if y<height-1 and x<width-1:
                            r = (x+1,y+1)
                            if n in down_rights:
                                e = down_rights[n]
                                if e in paths:
                                    print("\\",end="")
                                    drew_diag=True
                        if y>0 and x<width-1:
                            n = (x, y+1)
                            r = (x+1,y)
                            if n in up_rights:
                                e = up_rights[n]
                                if e in paths:
                                    print("/",end="")
                                    assert(not drew_diag)
                                    drew_diag=True

                    if not drew_diag:
                        print(" ", end="")
                print()
        print("s SATISFIABLE")
        sys.exit(10)
    else:
        print("s UNSATISFIABLE")
        sys.exit(20)


if __name__ == '__main__':
    import sys
    monosat_args = ['-ruc'] #default argument for MonoSAT; enables the heuristics described in "Routing Under Constraints", FMCAD 2016, A. Nadel
    if len(sys.argv)<2:
        print("Usage: router.py [monosat arguments] filename.pcrt")
        sys.exit(1)

    use_maxflow=False
    for i,arg in enumerate(sys.argv):
        if sys.argv[i].startswith("--use-maxflow"):
            use_maxflow = True
            del(sys.argv[i])
            break

    outputFile = None
    for i,arg in enumerate(sys.argv):
        if sys.argv[i].startswith("--output"):
            outputFile = sys.argv[i+1]
            del(sys.argv[i])
            del(sys.argv[i])
            break


    if len(sys.argv) > 2:
        monosat_args = sys.argv[1:-1]
    route(sys.argv[-1], monosat_args,use_maxflow,True, outputFile)

