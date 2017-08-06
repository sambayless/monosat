from monosat import *
from time import time
from collections import defaultdict
import pcrt
import itertools

def route_multi(filename, monosat_args, maxflow_enforcement_level, routing_separation_enforcement_style=0,graph_separation_enforcement_style=1,zeroEdgeWeights=False):
    (width, height), nets, constraints = pcrt.read(filename)
    print(filename)
    print("Width = %d, Height = %d, %d nets, %d constraints" % (width, height, len(nets), len(constraints)))

    if (len(monosat_args) > 0):
        args = " ".join(monosat_args)
        print("Monosat args: " + args)
        Monosat().init(args)



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
                #if maxflow_enforcement_level >= 2 and g==flow_graph:
                in_grid[g][(x, y)] = g.addNode("in_%d_%d" % (x, y))
                fromGrid[out_grid[g][(x, y)]] = (x, y)
                fromGrid[in_grid[g][(x, y)]] = (x, y)
                e = g.addEdge(in_grid[g][(x, y)], out_grid[g][(x, y)], 1)  # add an edge with constant capacity of 1
                vertex_grid[g][(x,y)]=e
                #else:
                #    in_grid[g][(x, y)] = n
                if(g !=flow_graph):
                    vertices[(x,y)].append(e)

    print("Adding edges")

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
                    Assert(g.addEdge(out_grid[g][n], in_grid[g][r]))  # add a _directed_ edge from n to r
                if flow_graph is not None:
                    e1 = (flow_graph.addEdge(out_grid[flow_graph][n],
                                           in_grid[flow_graph][r])) # add a _directed_ edge from n to r
                    if routing_separation_enforcement_style  > 0:
                        flow_graph_edges[(n, r)] = e1
                        flow_graph_edges[(r, n)] = e1
                        flow_grap_edge_list[n].append(e1)
                        flow_grap_edge_list[r].append(e1)
                    else:
                        Assert(e1)
            elif allow_in:
                # e = g.addEdge(out_grid[r], in_grid[n])  # add a _directed_ edge from r to n

                for g in graphs:
                    Assert(g.addEdge(out_grid[g][r], in_grid[g][n]))  # add a _directed_ edge from n to r

                if flow_graph is not None:
                    e1 = (flow_graph.addEdge(out_grid[flow_graph][r],
                                           in_grid[flow_graph][n])) # add a _directed_ edge from n to r
                    #AssertEq(Or(edge_set), e)  # there is an edge in the flow graph iff one of the other graphs has an edge

                    if routing_separation_enforcement_style  > 0:
                        flow_graph_edges[(n, r)] = e1
                        flow_graph_edges[(r, n)] = e1
                        flow_grap_edge_list[n].append(e1)
                        flow_grap_edge_list[r].append(e1)
                    else:
                        Assert(e1)



            else:
                e = None
        else:

            for g in graphs:
                Assert(g.addEdge(out_grid[g][n], in_grid[g][r]))  # add a _directed_ edge from n to r
                Assert(g.addEdge(out_grid[g][r], in_grid[g][n]))
            if flow_graph is not None:
                e1 = (flow_graph.addEdge(out_grid[flow_graph][r],
                                          in_grid[flow_graph][n]))  # add a _directed_ edge from n to r
                e2 = (flow_graph.addEdge(out_grid[flow_graph][n],
                                          in_grid[flow_graph][r]))  # add a _directed_ edge from r to n
                if routing_separation_enforcement_style > 0:
                    AssertEq(e1,e2)
                    flow_grap_edge_list[n].append(e1)
                    flow_grap_edge_list[r].append(e1)
                    flow_graph_edges[(n,r)] = e1
                    flow_graph_edges[(r, n)] = e1
                else:
                    Assert(e1)
                    Assert(e2)



    for x in range(width):
        for y in range(height):
            n = (x, y)
            if x < width - 1:
                r = (x + 1, y)
                addEdge(n,r)

            if y < height - 1:
                r = (x, y + 1)
                addEdge(n, r)

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

    uses_bv = (flow_graph and routing_separation_enforcement_style >= 2) or (graph_separation_enforcement_style>=2)

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



    if flow_graph and routing_separation_enforcement_style==1:
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
    elif flow_graph and routing_separation_enforcement_style == 2:
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

        on_path = dict()
        for y in range(height):
            for x in range(width):
                on_path[(x, y)] = -1

        for i, net in enumerate(nets):
            for n2 in net[1:]:
                r = reachset[i][n2]
                g = graphs[i]
                path = g.getPath(r,False)
                for n in path:
                    x,y = fromGrid[n]
                    on_path[(x,y)]=i

        for y in range(height):
            for x in range(width):
                n = (x, y)
                if (x, y) in net_nodes:
                    print("*", end="")
                else:
                    print(" ", end="")
                if x < width - 1:
                    r = (x + 1, y)

                    if on_path[n] >=0 and on_path[n]==on_path[r]:

                        print("-", end="")

                    else:
                        print(" ", end="")
            print()
            for x in range(width):
                n = (x, y)
                if y < height - 1:
                    d = (x, y + 1)
                    if  on_path[n] >=0 and on_path[n] == on_path[d]:

                        print("| ", end="")

                    else:
                        print("  ", end="")
            print()
        print("s SATISFIABLE")
        sys.exit(10)
    else:
        print("s UNSATISFIABLE")
        print("UNSAT")
        sys.exit(20)

if __name__ == '__main__':
    import sys
    print("Command: " + str(sys.argv))
    print("Command: " + str(sys.argv), file=sys.stderr)
    maxflow_enforcement_level = 0
    graph_separation_enforcement_style=0
    routing_separation_enforcement_style=0
    zeroEdgeWeights=True
    monosat_args = ['-ruc'] #default argument for MonoSAT; enables the heuristics described in "Routing Under Constraints", FMCAD 2016, A. Nadel
    found=True
    while(found):

        found=False
        if len(sys.argv[1])==0:
            found=True
            sys.argv = sys.argv[0:1] + sys.argv[2:]

        elif sys.argv[1].startswith("--use-maxflow="):
            found=True

            maxflow_enforcement_level = int(sys.argv[1][len("--use-maxflow="):])
            sys.argv = sys.argv[0:1] +  sys.argv[2:]
        elif sys.argv[1].startswith( "--separate-graphs="):
            found = True
            graph_separation_enforcement_style = int(sys.argv[1][len("--separate-graphs="):])

            sys.argv = sys.argv[0:1] + sys.argv[2:]
        elif sys.argv[1].startswith( "--enforce-separate="):
            found = True
            routing_separation_enforcement_style = int(sys.argv[1][len("--enforce-separate="):])

            sys.argv = sys.argv[0:1] + sys.argv[2:]
        elif sys.argv[1].startswith( "--zero-edges="):
            found = True
            zeroEdgeWeights = int(sys.argv[1][len("--zero-edges="):])>0
            sys.argv = sys.argv[0:1] + sys.argv[2:]

    if len(sys.argv) > 2:
        monosat_args = sys.argv[1:-1]
    route_multi(sys.argv[-1], monosat_args, maxflow_enforcement_level,routing_separation_enforcement_style,graph_separation_enforcement_style,zeroEdgeWeights)
