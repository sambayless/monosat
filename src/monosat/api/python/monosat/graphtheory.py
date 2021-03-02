# The MIT License (MIT)
#
# Copyright (c) 2014, Sam Bayless
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
# associated documentation files (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge, publish, distribute,
# sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or
# substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
# NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
# OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


import monosat.monosat_c
import sys
from monosat.bvtheory import BitVector
from monosat.logic import *
from monosat.manager import Manager

debug = False


# Collects a set of graphs to encode together into a formula
class GraphManager(metaclass=Manager):

    def __init__(self):
        self.graphs = []

    def clear(self):
        self.graphs = []

    def addGraph(self, g):
        self.graphs.append(g)


class Graph:
    class GraphType:
        int = 1
        float = 2
        rational = 3

    def __init__(self, graph_type=1):
        self._monosat = monosat.monosat_c.Monosat()
        manager = GraphManager()
        manager.addGraph(self)
        self.graph = self._monosat.newGraph()

        self.id = len(manager.graphs)
        self.has_any_bv_edges = False
        self.has_any_non_bv_edges = False
        self.nodes = 0
        self.numedges = 0
        self.names = dict()
        self.nodemap = dict()  # maps node names to node objects
        self.out_edges = []
        self.in_edges = []
        self.out_edge_map = []
        self.in_edge_map = []
        self.queries = []
        self.queryLookup = dict()
        self.alledges = []
        self.distance_rational_queries = []
        self.distance_float_queries = []

        self.dbg_reaches = []
        self.graph_type = graph_type
        self.edge_priorities = []
        # map from variables to edges
        self.all_undirectededges = []

        self.edgemap = dict()
        self.acyclic_querries = []

    def __getitem__(self, key):
        if isinstance(key, str):
            return self.nodemap[key]
        if isinstance(key, tuple) or isinstance(obj, list):
            if len(key) == 2:
                # interpret the key as an edge
                return self.getEdge(
                    self.nodemap[str(key[0])], self.nodemap[str(key[1])]
                )
        return self.nodemap[str(key)]

    def __contains__(self, item):
        return str(item) in self.nodemap.keys()

    def assignWeightsTo(self, w):
        self._monosat.assignWeightsTo(self.graph, w)

    def enforceRouting(self, source, destination, nets, maxflowlit):
        netlits = []

        for edge_lits, reach_lits, disabled_edge in nets:
            netlits.append(
                (
                    [x.getLit() for x in edge_lits],
                    [x.getLit() for x in reach_lits],
                    disabled_edge.getLit(),
                )
            )
        self._monosat.enforceRouting(
            self.graph, source, destination, netlits, maxflowlit.getLit()
        )

    def writeDot(self, out=sys.stdout, writeModel=True):
        print("digraph{", file=out)
        for n in range(self.nodes):
            print("n%d" % (n), file=out)

        for (v, w, var, weight) in self.getEdges():
            if not writeModel:
                if weight:
                    print(
                        'n%d -> n%d [label="v%s, w%s"]' % (v, w, str(var), str(weight)),
                        file=out,
                    )
                else:
                    print('n%d -> n%d [label="v%s"]' % (v, w, str(var)), file=out)
            else:
                edgeVal = var.value()
                if edgeVal is None:
                    edgecol = "black"
                elif edgeVal:
                    edgecol = "red"
                else:
                    edgecol = "blue"

                if weight:
                    if edgeVal is not None:
                        weightVal = weight.value()
                        print(
                            'n%d -> n%d [label="v%s, w%s=%s", color=%s]'
                            % (v, w, str(var), str(weight), str(weightVal), edgecol),
                            file=out,
                        )
                    else:
                        print(
                            'n%d -> n%d [label="v%s, w%s"]'
                            % (v, w, str(var), str(weight)),
                            file=out,
                        )
                else:
                    print(
                        'n%d -> n%d [label="v%s", color=%s]'
                        % (v, w, str(var), edgecol),
                        file=out,
                    )

        print("}", file=out)

    def addNode(self, name=None):
        n = self._monosat.newNode(self.graph)
        self.nodes = n + 1
        self.out_edges.append([])
        self.in_edges.append([])

        if name is None:
            name = str(n)
        else:
            name = str(name)

        if name is not None and name in self.nodemap:
            raise ValueError("Node %s already exists" % (str(name)))

        self.names[n] = name
        self.nodemap[name] = n
        self.out_edge_map.append(dict())
        self.in_edge_map.append(dict())

        return n

    def getSymbol(self, node):
        return self.names[node]

    def getMaxFlow(self, flowlit):
        return self._monosat.getModel_MaxFlow(self.graph, flowlit.getLit())

    def getEdgeFlow(self, flowlit, edgelit, force_acyclic_flow=False):
        if force_acyclic_flow:
            return self._monosat.getModel_AcyclicEdgeFlow(
                self.graph, flowlit.getLit(), edgelit.getLit()
            )
        else:
            return self._monosat.getModel_EdgeFlow(
                self.graph, flowlit.getLit(), edgelit.getLit()
            )

    """
    Get a path in the graph that satisfies the reachability or shortest path lit, if the shortest path lit is true in the model
    If 'return_edge_lits' is True, then return the path as a list of edge literals. Otherwise, returns the path as a list of nodes.
    Must not be caled before solve().
    """

    def getPath(self, reach_or_shortest_path_lit, return_edge_lits=False):
        if not return_edge_lits:
            return self._monosat.getModel_Path_Nodes(
                self.graph, reach_or_shortest_path_lit.getLit()
            )
        else:

            lits = self._monosat.getModel_Path_EdgeLits(
                self.graph, reach_or_shortest_path_lit.getLit()
            )
            if lits is None:
                return None
            lit_list = []
            for lit in lits:
                lit_list.append(Var(lit))

            return lit_list

    # Returns either a list containing all the edges from f to t, or, if there is only 1 edge from f to t, returns that edge.
    # throws a ValueError if there are no edges from f to t.
    def getEdge(self, f, t):
        if t in self.out_edge_map[f]:
            edges = self.out_edge_map[f][t]
            assert len(edges) > 0
            if len(edges) == 1:
                return edges[
                    0
                ]  # the common case is that the graph has no multiedges, so handle that specially
            else:
                return edges
        raise ValueError("No edge " + str(f) + "->" + str(t))

    def getAllEdges(self, undirected=False):
        if undirected:
            return self.all_undirectededges
        else:
            return self.alledges

    # returns the variable corresponding to the backward (directed) edge of the edge
    # corresponding to this variable, if such an edge exists. Returns None otherwise.
    def backEdgeVar(self, v):
        (v, w, var, weight) = self.getEdgeFromVar(v)
        if v in self.out_edge_map[w]:
            edges = self.out_edge_map[w][v]
            assert len(edges) > 0
            return edges[0][2]
        return None

    def hasEdge(self, f, t):
        return t in self.out_edge_map[f] and len(self.out_edge_map[f][t]) > 0
        # for (v,u,var,weight) in self.out_edges[f]:
        #   if(u==t):
        #       return True;

        # return False;

    def newEdgeSet(self, edges, enforceEdgeAssignments=True):
        for v in edges:
            assert v.getLit() in self.edgemap
        edgelits = [v.getLit() for v in edges]
        self._monosat.newEdgeSet(self.graph, edgelits, enforceEdgeAssignments)

    # add edge from v to w
    def addEdge(self, v, w, weight=1):
        while v >= self.numNodes() or w >= self.numNodes():
            self.addNode()

        if weight and isinstance(weight, float):
            assert (
                    self.graph_type == Graph.GraphType.float
                    or self.graph_type == Graph.GraphType.rational
            )
        elif weight and isinstance(weight, tuple):
            assert self.graph_type == Graph.GraphType.rational

        if weight and isinstance(weight, BitVector):
            assert self.graph_type == Graph.GraphType.int
            self.has_any_bv_edges = True
            assert not self.has_any_non_bv_edges
            var = Var(self._monosat.newEdge_bv(self.graph, v, w, weight.getID()))
        else:
            self.has_any_non_bv_edges = True
            assert not self.has_any_bv_edges
            if self.graph_type == Graph.GraphType.int:
                var = Var(self._monosat.newEdge(self.graph, v, w, weight))
            elif self.graph_type == Graph.GraphType.float:
                var = Var(self._monosat.newEdge_double(self.graph, v, w, weight))

        e = (v, w, var, weight)
        self.alledges.append(e)
        self.numedges = self.numedges + 1
        self.out_edges[v].append(e)
        self.in_edges[w].append(e)
        self.edgemap[e[2].getLit()] = e
        if w not in self.out_edge_map[v].keys():
            self.out_edge_map[v][w] = list()

        if v not in self.in_edge_map[w].keys():
            self.in_edge_map[w][v] = list()
        self.out_edge_map[v][w].append(e)
        self.in_edge_map[w][v].append(e)
        return e[2]

    def addUndirectedEdge(self, v, w, weight=1):
        while v >= self.numNodes() or w >= self.numNodes():
            self.addNode()
        if weight and isinstance(weight, float):
            assert (
                    self.graph_type == Graph.GraphType.float
                    or self.graph_type == Graph.GraphType.rational
            )
        elif weight and isinstance(weight, tuple):
            assert self.graph_type == Graph.GraphType.rational

        if weight and isinstance(weight, BitVector):
            assert self.graph_type == Graph.GraphType.int
            self.has_any_bv_edges = True
            assert not self.has_any_non_bv_edges
            v1 = Var(self._monosat.newEdge_bv(self.graph, v, w, weight.getID()))
            v2 = Var(self._monosat.newEdge_bv(self.graph, w, v, weight.getID()))
        else:
            self.has_any_non_bv_edges = True
            assert not self.has_any_bv_edges
            if self.graph_type == Graph.GraphType.int:
                v1 = Var(self._monosat.newEdge(self.graph, v, w, weight))
                v2 = Var(self._monosat.newEdge(self.graph, w, v, weight))
            elif self.graph_type == Graph.GraphType.float:
                v1 = Var(self._monosat.newEdge_double(self.graph, v, w, weight))
                v2 = Var(self._monosat.newEdge_double(self.graph, w, v, weight))

        e1 = (v, w, v1, weight)
        self.alledges.append(e1)
        self.numedges = self.numedges + 1
        self.out_edges[v].append(e1)
        self.in_edges[w].append(e1)
        e2 = (w, v, v2, weight)

        self.alledges.append(e2)
        self.numedges = self.numedges + 1
        self.out_edges[w].append(e2)
        self.in_edges[v].append(e2)
        self.edgemap[v1.getLit()] = e1
        self.edgemap[v2.getLit()] = e2
        AssertEq(v1, v2)
        self.all_undirectededges.append(e1)

        if w not in self.out_edge_map[v].keys():
            self.out_edge_map[v][w] = list()

        if v not in self.in_edge_map[w].keys():
            self.in_edge_map[w][v] = list()

        self.out_edge_map[v][w].append(e1)
        self.in_edge_map[w][v].append(e1)

        if v not in self.out_edge_map[w].keys():
            self.out_edge_map[w][v] = list()

        if w not in self.in_edge_map[v].keys():
            self.in_edge_map[v][w] = list()
        # add the edge twice, one for each direction.
        self.out_edge_map[w][v].append(e2)
        self.in_edge_map[v][w].append(e2)
        return v1

    def setEdgePriority(self, edgeVar, priority):
        self.edge_priorities.append((edgeVar, priority))

    def numNodes(self):
        return self.nodes

    def numEdges(self):
        return self.numedges

    def nNodes(self):
        return self.nodes

    def nEdges(self):
        return self.numedges

    def getNodes(self):
        return range(self.nodes)

    def getEdgeFromVar(self, var):
        return self.edgemap[var.getLit()]

    def getEdges(self, node=-1, undirected=False):
        if node >= 0:
            for edge in self.out_edges[node]:
                yield edge
            if undirected:
                for edge in self.in_edges[node]:
                    yield edge
        else:
            for node in self.out_edges:
                for edge in node:
                    yield edge

    def getOutgoingEdges(self, node=-1):
        if node >= 0:
            for edge in self.out_edges[node]:
                yield edge
        else:
            for node in self.out_edges:
                for edge in node:
                    yield edge

    def getIncomingEdges(self, node=-1):
        if node >= 0:
            for edge in self.in_edges[node]:
                yield edge
        else:
            for node in self.in_edges:
                for edge in node:
                    yield edge

    def getEdgeVars(self, node=-1):
        if node >= 0:
            for edge in self.out_edges[node]:
                yield edge[2]
        else:
            for node in self.out_edges:
                for edge in node:
                    yield edge[2]

    def getOutgoingEdgeVars(self, node=-1):
        if node >= 0:
            for edge in self.out_edges[node]:
                yield edge[2]
        else:
            for node in self.out_edges:
                for edge in node:
                    yield edge[2]

    def getIncomingEdgeVars(self, node=-1):
        if node >= 0:
            for edge in self.in_edges[node]:
                yield edge[2]
        else:
            for node in self.in_edges:
                for edge in node:
                    yield edge[2]

    # def createSteinerTree(self):
    # s = Graph.SteinerTree(len(self.steiners))
    # self.steiners.append(s)
    # return s

    # Shadow the graph theory with a (slow) implementation in the circuit, for debugging purposes
    """def _reachesAnyCircuit(self,start,n=None):
        assert(False)
        if(n is None):
            n=self.numNodes()
        n=int(n)
        #bellman-ford:
        if(n>self.numNodes()):
            n=self.numNodes();
        reaches=[Var(False)]*self.numNodes()
        reaches[start]=Var(True)
        for i in range(0,n):            
            #For each edge:
            for (v,w,var) in self.getEdges():
                reaches[w]=(var & reaches[v])| reaches[w]
        return reaches"""

    def distance_rational_leq(self, start, to, distance_fraction):
        v = (
            Var()
        )  # "distance_rational_leq(%d,%d,%d,%d,%d)"%(self.id, start,to,distance_numerator,distance_denominator))

        self.distance_rational_queries.append(("leq", start, distance_fraction, to, v))
        return v

    def distance_rational_lt(self, start, to, distance_fraction):
        v = (
            Var()
        )  # "distance_rational_leq(%d,%d,%d,%d,%d)"%(self.id, start,to,distance_numerator,distance_denominator))

        self.distance_rational_queries.append(("lt", start, distance_fraction, to, v))
        return v

    def distance_leq(self, start, to, distance):
        if isinstance(distance, tuple):
            return self.distance_rational_leq(start, to, distance)
        if isinstance(distance, BitVector):
            v = Var(
                self._monosat.shortestPath_leq_bv(
                    self.graph, start, to, distance.getID()
                )
            )
        else:
            v = Var(
                self._monosat.shortestPath_leq_const(self.graph, start, to, distance)
            )
        # v = Var() #"distance_float_leq(%d,%d,%d,%d)"%(self.id, start,to,distance))
        # self.distance_float_queries.append(('leq',start,distance,to,v))
        return v

    def distance_lt(self, start, to, distance):
        if isinstance(distance, tuple):
            return self.distance_rational_lt(start, to, distance)
        if isinstance(distance, BitVector):
            v = Var(
                self._monosat.shortestPath_lt_bv(
                    self.graph, start, to, distance.getID()
                )
            )
        else:
            v = Var(
                self._monosat.shortestPath_lt_const(self.graph, start, to, distance)
            )
        # v = Var() #"distance_float_leq(%d,%d,%d,%d)"%(self.id, start,to,distance))
        # self.distance_float_queries.append(('lt',start,distance,to,v))
        return v

    def reaches(self, start, to, withinSteps=None):
        if withinSteps is None or withinSteps < 0:
            withinSteps = None
        else:
            withinSteps = int(withinSteps)

        if withinSteps is not None and withinSteps >= self.numNodes():
            withinSteps = None

        # Check to see if we have already encoded this property
        if (start, to, withinSteps) in self.queryLookup:
            return self.queryLookup[(start, to, withinSteps)]
        if withinSteps is None:
            v = Var(self._monosat.reaches(self.graph, start, to))
        else:
            v = Var(
                self._monosat.shortestPathUnweighted_leq_const(
                    self.graph, start, to, withinSteps
                )
            )
            # v = Var("distance_leq(%d,%d,%d,%d)"%(self.id, start,to,withinSteps) if withinSteps is not None else "reaches(%d,%d,%d)"%(self.id, start,to))
        # self.queries.append((start,withinSteps,to,v))
        self.queryLookup[(start, to, withinSteps)] = v
        return v

    # Check each node for reachability from v in at most n steps
    def reachesAny(self, start, n=None):
        if n is None or n < 0:
            n = None
        else:
            n = int(n)
        if n is not None and n > self.numNodes():
            n = self.numNodes()
        reaches = []
        for i in range(self.numNodes()):
            if (start, i, n) in self.queryLookup:
                reaches.append(self.queryLookup[(start, i, n)])
            else:
                if n is None:
                    v = Var(self._monosat.reaches(self.graph, start, i))
                else:
                    v = Var(
                        self._monosat.shortestPathUnweighted_leq_const(
                            self.graph, start, i, n
                        )
                    )
                reaches.append(v)
                self.queryLookup[(start, i, n)] = v

        return reaches

    def reachesBackward(self, start, to):
        """
        A reaches query that traverses edges backwards (eg, u reaches v if the edge v->u is in the graph)
        """
        v = Var(self._monosat.reachesBackward(self.graph, start, to))
        return v

    def onPath(self, nodeOnPath, start, to):
        """
        True iff there exists a path from start to nodeOnPath, AND there exists a path from nodeOnPath to 'to'
        """
        v = Var(self._monosat.onPath(self.graph, nodeOnPath, start, to))
        return v

    def minimumSpanningTreeLessEq(self, minweight):
        v = Var(self._monosat.minimumSpanningTree_leq(self.graph, minweight))
        return v

    def acyclic(self, directed=True):
        if directed:
            v = Var(self._monosat.acyclic_directed(self.graph))
        else:
            v = Var(self._monosat.acyclic_undirected(self.graph))

        return v

    def AssertMinimumSpanningTreeLessEq(self, minweight):
        Assert(self.minimumSpanningTreeLessEq(minweight))

    def edgeInMinimumSpanningTree(self, edgeVar):
        varInTreeOrDisabled = Var()

        self.mstEdgeQueries.append((edgeVar, varInTreeOrDisabled))

        return And(varInTreeOrDisabled, edgeVar)

    def AssertEdgeInMinimumSpanningTree(self, minweight):
        Assert(self.minimumSpanningTree(minweight))

    def maxFlowGreaterOrEqualTo(self, s, t, flow):
        if isinstance(flow, BitVector):
            v = Var(self._monosat.maximumFlow_geq_bv(self.graph, s, t, flow.getID()))
        else:
            v = Var(self._monosat.maximumFlow_geq(self.graph, s, t, flow))
        return v

    def AssertMaxFlowGreaterOrEqualTo(self, s, t, flow):
        v = self.maxFlowGreaterOrEqualTo(s, t, flow)
        Assert(v)

    def AssertMaxFlowLessOrEqualTo(self, s, t, flow):
        v = self.maxFlowGreaterOrEqualTo(s, t, flow + 1)
        Assert(Not(v))


    def connectedComponentsGreaterOrEqualTo(self, components):
        v =  Var(self._monosat.connectedComponents_geq_const(self.graph, components))
        return v

    def connectedComponentsLessOrEqualTo(self, components):
        v = Not(self.connectedComponentsGreaterOrEqualTo(components+1))
        return v

    def AssertConnectedComponentsGreaterOrEqualTo(self, min_components):
        Assert(self.connectedComponentsGreaterOrEqualTo(min_components))

    def AssertConnectedComponentsEqualTo(self, min_components):
        Assert(self.connectedComponentsGreaterOrEqualTo(min_components))
        Assert(self.connectedComponentsLessOrEqualTo(min_components))

    def AssertConnectedComponentsLessOrEqualto(self, min_components):
        Assert(self.connectedComponentsLessOrEqualTo(min_components))

    def draw(self):

        print("digraph{")

        for n in range(self.nodes):
            print("n%d" % (n))

        for (v, w, var, weight) in self.getAllEdges():
            if weight is not None:
                print("""n%d->n%d [label="%d w=%d"]""" % (v, w, var.getVar(), weight))
            else:
                print("""n%d->n%d [label="%d"]""" % (v, w, var.getVar()))

        print("}")
