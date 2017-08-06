from monosat import *
import argparse
from time import time
from collections import defaultdict
import pcrt
import itertools

def route_multi(filename, monosat_args, maxflow_enforcement_level, flowgraph_separation_enforcement_style=0,graph_separation_enforcement_style=1,zeroEdgeWeights=False,draw_solution=True):
    (width, height),diagonals,nets,constraints,disabled = pcrt.read(filename)
    print(filename)
    print("Width = %d, Height = %d, %d nets, %d constraints" % (width, height, len(nets), len(constraints)))
    if diagonals:
        print("45-degree routing")
    else:
        print("90-degree routing")

    if (len(monosat_args) > 0):
        args = " ".join(monosat_args)
        print("Monosat args: " + args)
        Monosat().init(args)

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

    graphs = []
    all_graphs=[]
    for i in nets:
        graphs.append(Graph())
        if zeroEdgeWeights:
            graphs[-1].assignWeightsToZero() # this enables a heuristic on these graphs, from the RUC paper, which sets assigned edges to zero weight, to encourage edge-reuse in solutions
        all_graphs.append(graphs[-1])
    flow_graph=None
    flow_graph_edges = dict()
    flow_grap_edge_list = collections.defaultdict(list)
    if maxflow_enforcement_level>=1:
        flow_graph = Graph()
        all_graphs.append(flow_graph)

    out_grid = dict()
    print("Building grid")
    in_grid = dict()
    vertex_grid = dict()
    vertices=dict()
    fromGrid = dict()
    for g in all_graphs:
        out_grid[g] = dict()
        in_grid[g] = dict()
        vertex_grid[g] = dict()

    for x in range(width):
        for y in range(height):
            vertices[(x,y)]=[]
            for g in all_graphs:
                n = g.addNode("%d_%d" % (x, y))
                out_grid[g][(x, y)] = n

                in_grid[g][(x, y)] = g.addNode("in_%d_%d" % (x, y))
                fromGrid[out_grid[g][(x, y)]] = (x, y)
                fromGrid[in_grid[g][(x, y)]] = (x, y)
                e = g.addEdge(in_grid[g][(x, y)], out_grid[g][(x, y)], 1)  # add an edge with constant capacity of 1
                vertex_grid[g][(x,y)]=e

                if(g !=flow_graph):
                    vertices[(x,y)].append(e)

    print("Adding edges")
    disabled_nodes = set(disabled)
    undirected_edges = dict()
    start_nodes = set()
    end_nodes = set()
    net_nodes = set()
    for net in nets:

        start_nodes.add(net[0])
        for n in net[1:]:
            end_nodes.add(n)
        for (x, y) in net:
            net_nodes.add((x, y))

    def addEdge(n,r):
        e = None
        if n not in disabled_nodes and r not in disabled_nodes:
            if n in net_nodes or r in net_nodes:
                allow_out = True
                allow_in = True
                if n in start_nodes or r in end_nodes:
                    allow_in = False
                if n in end_nodes or r in start_nodes:
                    allow_out = False
                assert (not (allow_in and allow_out))
                if allow_out:

                    for g in graphs:
                        eg = (g.addEdge(out_grid[g][n], in_grid[g][r]))  # add a _directed_ edge from n to r
                        if e is None:
                            e = eg
                        undirected_edges[eg]=e
                        Assert(eg)
                    if flow_graph is not None:
                        ef = (flow_graph.addEdge(out_grid[flow_graph][n],
                                               in_grid[flow_graph][r])) # add a _directed_ edge from n to r

                        if flowgraph_separation_enforcement_style  > 0:
                            flow_graph_edges[(n, r)] = ef
                            flow_graph_edges[(r, n)] = ef
                            flow_grap_edge_list[n].append(ef)
                            flow_grap_edge_list[r].append(ef)
                        else:
                            Assert(ef)
                elif allow_in:
                    # e = g.addEdge(out_grid[r], in_grid[n])  # add a _directed_ edge from r to n

                    for g in graphs:

                        eg = (g.addEdge(out_grid[g][r], in_grid[g][n]))  # add a _directed_ edge from n to r
                        if e is None:
                            e = eg
                        undirected_edges[eg]=e
                        Assert(eg)

                    if flow_graph is not None:
                        ef = (flow_graph.addEdge(out_grid[flow_graph][r],
                                               in_grid[flow_graph][n])) # add a _directed_ edge from n to r

                        if flowgraph_separation_enforcement_style  > 0:
                            flow_graph_edges[(n, r)] = ef
                            flow_graph_edges[(r, n)] = ef
                            flow_grap_edge_list[n].append(ef)
                            flow_grap_edge_list[r].append(ef)
                        else:
                            Assert(ef)



                else:
                    e = None
            else:

                for g in graphs:
                    eg = (g.addEdge(out_grid[g][n], in_grid[g][r]))  # add a _directed_ edge from n to r
                    if e is None:
                        e = eg
                    undirected_edges[eg]=e
                    Assert(eg)
                    eg2 = (g.addEdge(out_grid[g][r], in_grid[g][n]))
                    Assert(eg2)
                    undirected_edges[eg2]=e #map e2 to e

                if flow_graph is not None:
                    ef = (flow_graph.addEdge(out_grid[flow_graph][r],
                                              in_grid[flow_graph][n]))  # add a _directed_ edge from n to r
                    ef2 = (flow_graph.addEdge(out_grid[flow_graph][n],
                                              in_grid[flow_graph][r]))  # add a _directed_ edge from r to n

                    if flowgraph_separation_enforcement_style > 0:
                        AssertEq(ef,ef2)
                        flow_grap_edge_list[n].append(ef)
                        flow_grap_edge_list[r].append(ef)
                        flow_graph_edges[(n,r)] = ef
                        flow_graph_edges[(r, n)] = ef
                    else:
                        Assert(ef)
                        Assert(ef2)

        return e


    for x in range(width):
        for y in range(height):
            n = (x, y)
            if n in disabled_nodes:
                continue
            if x < width - 1:
                r = (x + 1, y)
                e=addEdge(n,r)
                if e is not None and draw_solution:
                    rights[n] = e
                    lefts[r] = e

            if y < height - 1:
                r = (x, y + 1)
                e=addEdge(n, r)
                if e is not None and draw_solution:
                    downs[n] = e
                    ups[r] = e


            if diagonals:
                #if 45 degree routing is enabled
                diag_up=None
                diag_down=None
                if x<width-1 and y<height-1:
                    r = (x+1,y+1)
                    e = addEdge(n,r)
                    diag_down = e
                    if e is not None and draw_solution:
                        down_rights[n] = e
                        up_lefts[r] = e

                if x>0 and y<height-1 and False:
                    r = (x-1,y+1)
                    e = addEdge(n,r)
                    diag_up = e
                    if e is not None and draw_solution:
                        down_lefts[n] = e
                        up_rights[r] = e
                if diag_up and diag_down:
                    AssertNand(diag_up, diag_down) #cannot route both diagonals

    vertex_used=None

    if len(constraints) > 0:
        print("Enforcing constraints")
        vertex_used = dict()
        for x in range(width):
            for y in range(height):
                vertex_used[(x, y)] = Or(vertices[(x, y)])  # A vertex is used exactly if one of its edges is enabled

        for constraint in constraints:
            # a constraint is a list of vertices of which at most one can be used
            vertex_used_list = []
            for node in constraint:
                vertex_used_list.append(vertex_used[node])
            if (len(vertex_used_list) < 20):
                for a, b in itertools.combinations(vertex_used_list, 2):
                    AssertOr(~a, ~b)  # in every distinct pair of edges, at least one must be false
            else:
                AssertAtMostOne(vertex_used_list)  # use more expensive, but more scalable, built-in AMO theory

    uses_bv = (flow_graph and flowgraph_separation_enforcement_style >= 2) or (graph_separation_enforcement_style>=2)

    print("Enforcing separation")
    #force each vertex to be in at most one graph
    for x in range(width):
        for y in range(height):
            n = (x, y)
            if n not in net_nodes:
                if graph_separation_enforcement_style==1:
                    AMO(vertices[n])
                else:
                    #rely on the uniqueness bv encoding below to force at most one graph assign per vertex
                    assert(uses_bv)

                if flow_graph is not None or uses_bv:

                    if vertex_used is None: #only create this lazily, if needed
                        vertex_used = dict()
                        for x in range(width):
                            for y in range(height):
                                vertex_used[(x, y)] = Or(
                                    vertices[(x, y)])  # A vertex is used exactly if one of its edges is enabled
                    if flow_graph is not None:
                        #Assert that iff this vertex is in _any_ graph, it must be in the flow graph
                        AssertEq(vertex_used[(x, y)],vertex_grid[flow_graph][n])

    if uses_bv:

        vertices_bv = dict()
        bvwidth = math.ceil(math.log(len(nets)+1,2))
        print("Building BV (width = %d)" % (bvwidth))

        #bv 0 means unused
        for x in range(width):
            for y in range(height):
                #netbv = BitVector(bvwidth)
                netbv = [Var() for _ in range(bvwidth)]
                seen_bit = [False] * bvwidth # this is just for error checking this script
                for b in range(bvwidth):
                    #if the vertex is not used, set the bv to 0
                    AssertImplies(~vertex_used[(x, y)], ~netbv[b])

                for i in range(len(nets)):
                    net_n = i+1
                    seen_any_bits = False
                    for b in range(bvwidth):
                        bit = net_n & (1<<b)
                        if(bit):
                            AssertImplies(vertices[(x, y)][i], netbv[b] )
                            seen_bit[b]=True
                            seen_any_bits = True
                        else:
                            AssertImplies(vertices[(x, y)][i], ~netbv[b])
                    #AssertImplies(vertices[(x,y)][i],(netbv.eq(net_n)))
                    assert (seen_any_bits)
                if graph_separation_enforcement_style<3:
                    #rely on the above constraint AssertImplies(~vertex_used[(x, y)], ~netbv[b])
                    #to ensure that illegal values of netbv are disallowed
                    pass
                elif graph_separation_enforcement_style==3:
                    # directly rule out disallowed bit patterns
                    # len(nets)+1, because 1 is added to each net id for the above calculation (so that 0 can be reserved for no net)
                    for i in range(len(nets)+1,(1<<bvwidth)):
                        bits = []
                        for b in range(bvwidth):
                            bit = i & (1 << b)
                            if bit>0:
                                bits.append(netbv[b])
                            else:
                                bits.append(~netbv[b])
                        AssertNand(bits) #at least one of these bits cannot be assigned this way

                assert(all(seen_bit)) #all bits must have been set to 1 at some point, else the above constraints are buggy
                vertices_bv[(x, y)] = netbv



    if flow_graph and flowgraph_separation_enforcement_style==1:
        print("Enforcing (redundant) flow separation")
        #if two neighbouring nodes belong to different graphs, then
        #disable the edge between them.
        for x in range(width):
            for y in range(height):
                n = (x, y)

                if x < width - 1:
                    r = (x + 1, y)
                    if (n, r) in flow_graph_edges:
                        # if either end point is not is disabled, disable this edge... this is not technically required (since flow already cannot pass through unused vertices),
                        # but cheap to enforce and slightly reduces the search space.
                        AssertImplies(Or(Not(vertex_used[n]),Not(vertex_used[r])),Not(flow_graph_edges[(n, r)]))
                        any_same=false()
                        for i in range(len(vertices[n])):
                            # Enable this edge if both vertices belong to the same graph
                            same_graph = And(vertices[n][i], vertices[r][i])
                            AssertImplies(same_graph,flow_graph_edges[(n,r)])
                            any_same = Or(any_same,same_graph)
                            # Assert that if vertices[n] != vertices[r], then flow_graph_edges[(n,r)] = false
                            for j in range(i+1,len(vertices[r])):
                                # if the end points of this edge belong to different graphs, disable them.
                                AssertImplies(And(vertices[n][i], vertices[r][j]), Not(flow_graph_edges[(n, r)]))
                        AssertEq(flow_graph_edges[(n,r)],any_same)


                if y < height - 1:
                    r = (x, y + 1)
                    if (n, r) in flow_graph_edges:
                        # if either end point is not is disabled, disable this edge... this is not technically required (since flow already cannot pass through unused vertices),
                        # but cheap to enforce and slightly reduces the search space.
                        AssertImplies(Or(Not(vertex_used[n]), Not(vertex_used[r])), Not(flow_graph_edges[(n, r)]))
                        any_same = false()
                        for i in range(len(vertices[n])):
                            # Enable this edge if both vertices belong to the same graph
                            same_graph = And(vertices[n][i], vertices[r][i])
                            AssertImplies(same_graph, flow_graph_edges[(n, r)])
                            any_same = Or(any_same, same_graph)

                            # Assert that if vertices[n] != vertices[r], then flow_graph_edges[(n,r)] = false
                            for j in range(i + 1, len(vertices[r])):
                                # if the end points of this edge belong to different graphs, disable them.
                                AssertImplies(And(vertices[n][i], vertices[r][j]), Not(flow_graph_edges[(n, r)]))
                        AssertEq(flow_graph_edges[(n, r)], any_same)
    elif flow_graph and flowgraph_separation_enforcement_style == 2:
        print("Enforcing (redundant) BV flow separation")
        for x in range(width):
            for y in range(height):
                n = (x, y)

                if x < width - 1:
                    r = (x + 1, y)
                    if (n, r) in flow_graph_edges:
                        # if either end point is not is disabled, disable this edge... this is not technically required (since flow already cannot pass through unused vertices),
                        # but cheap to enforce and slightly reduces the search space.
                        # AssertImplies(Or(Not(vertex_used[n]),Not(vertex_used[r])),Not(flow_graph_edges[(n, r)]))

                        same_graph = BVEQ(vertices_bv[n], vertices_bv[r])  # And(vertices[n][i], vertices[r][i])
                        AssertEq(And(vertex_used[n], same_graph), flow_graph_edges[(n, r)])

                if y < height - 1:
                    r = (x, y + 1)
                    if (n, r) in flow_graph_edges:
                        # if either end point is not is disabled, disable this edge... this is not technically required (since flow already cannot pass through unused vertices),
                        # but cheap to enforce and slightly reduces the search space.
                        # AssertImplies(Or(Not(vertex_used[n]), Not(vertex_used[r])), Not(flow_graph_edges[(n, r)]))
                        same_graph = BVEQ(vertices_bv[n], vertices_bv[r])  # And(vertices[n][i], vertices[r][i])
                        AssertEq(And(vertex_used[n], same_graph), flow_graph_edges[(n, r)])

    for i,net in enumerate(nets):
        for n in net:
            Assert(vertices[n][i]) #terminal nodes must be assigned to this graph
            if(flow_graph):
                Assert(vertex_grid[flow_graph][n])


    print("Enforcing reachability")
    reachset = dict()
    # enforce reachability
    for i,net in enumerate(nets):
        reachset[i] = dict()
        n1 = net[0]
        for n2 in net[1:]:
            g = graphs[i]
            r = g.reaches(in_grid[g][n1], out_grid[g][n2])
            Assert(r)
            r.setDecisionPriority(1);  # decide reachability before considering regular variable decisions
            reachset[i][n2]=r


    if maxflow_enforcement_level >= 1:
        print("Enforcing flow")
        # add a source and dest node, with 1 capacity from source to each net start vertex, and 1 capacity from each net end vertex to dest
        g = flow_graph
        source = g.addNode()
        dest = g.addNode()
        for net in nets:
            Assert(g.addEdge(source, in_grid[g][net[0]], 1))  # directed edges!
            Assert(g.addEdge(out_grid[g][net[1]], dest, 1))  # directed edges!
        m = g.maxFlowGreaterOrEqualTo(source, dest, len(nets))
        Assert(m)
        if maxflow_enforcement_level == 3:
            m.setDecisionPriority(1);  # sometimes make decisions on the maxflow predicate.
        elif maxflow_enforcement_level == 4:
            m.setDecisionPriority(2);  # always make decisions on the maxflow predicate.
        else:
            m.setDecisionPriority(-1);  # never make decisions on the maxflow predicate.

    print("Solving...")

    if Solve():

        print("Solved!")
        if draw_solution:
            paths = set() #collect all the edges that make up the routing
            on_path = dict()
            drawn = set()
            for y in range(height):
                for x in range(width):
                    on_path[(x, y)] = -1

            for i, net in enumerate(nets):
                for n2 in net[1:]:
                    r = reachset[i][n2]
                    g = graphs[i]
                    path = g.getPath(r,True)
                    for e in path:
                        assert(e.value())
                        if e in undirected_edges.keys():
                            e = undirected_edges[e] #normalize the edge var
                            assert(e in up_rights.values() or e in  down_rights.values() or e in  downs.values() or e in rights.values())
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
                            e = rights[n];assert(e in undirected_edges.keys())
                            if e in paths:
                                print("-",end=""); drawn.add(e)
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
                            e = downs[n];assert(e in undirected_edges.keys())
                            if e in paths:
                                print("|",end=""); drawn.add(e)
                            else:
                                print(" ",end="")
                        else:
                            print(" ", end="")
                    drew_diag=False
                    if diagonals:
                        if y<height-1 and x<width-1:
                            r = (x+1,y+1)
                            if n in down_rights:
                                e = down_rights[n];assert(e in undirected_edges.keys())
                                if e in paths:
                                    print("\\",end=""); drawn.add(e)
                                    drew_diag=True
                        if y>0 and x<width-1:
                            n = (x, y+1)
                            r = (x+1,y)
                            if n in up_rights:
                                e = up_rights[n];assert(e in undirected_edges.keys())
                                if e in paths:
                                    print("/",end=""); drawn.add(e)
                                    assert(not drew_diag)
                                    drew_diag=True

                    if not drew_diag:
                        print(" ", end="")
                print()

            for e in paths:
                assert(e in drawn)
        print("s SATISFIABLE")
        sys.exit(10)
    else:
        print("s UNSATISFIABLE")
        sys.exit(20)

if __name__ == '__main__':
    import sys

    monosat_args = ['-ruc'] #default argument for MonoSAT; enables the heuristics described in "Routing Under Constraints", FMCAD 2016, A. Nadel

    parser = argparse.ArgumentParser(description='SAT-based, constrained multi-terminal routing')

    parser.add_argument('--use-maxflow', default = 0, type=int,  help='Set to >= 1 to enable redundant, over-approximative maximum flow constraints, which can help the solver prune bad solutions early. Options 2,3,4 control heuristic interactions between the flow constraints and the routing constraints in the solver.',choices=range(0,5))
    parser.add_argument('--separate-graphs',default=2, type=int, help='This controls the type of constraint that prevents nets from intersecting. All three are reasonable choices.',choices=range(1,4))
    parser.add_argument('--enforce-separate',default=0, type=int,  help='This controls the type of constraint used to prevent nets from intersecting with each other in the maximum flow constraint, IF maxflow constraints are used.',choices=range(0,4))
    parser.add_argument('--zero-edges',default=0, type=int,
                        help='This enables a heuristic which sets assigned edges to zero weight, to encourage edge-reuse in solutions in the solver.',choices=range(0,2))

    parser.add_argument('filename', type=str)

    args, unknown = parser.parse_known_args()


    if len(unknown)>0:
        print("Passing unrecognized arguments to monosat: " + str(unknown))
        monosat_args = unknown
    route_multi(args.filename, monosat_args, args.use_maxflow,args.enforce_separate,args.separate_graphs,args.zero_edges,True)
