/***************************************************************************************[Solver.cc]
 The MIT License (MIT)

 Copyright (c) 2018, Sam Bayless

 Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
 associated documentation files (the "Software"), to deal in the Software without restriction,
 including without limitation the rights to use, copy, modify, merge, publish, distribute,
 sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all copies or
 substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
 NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
 OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 **************************************************************************************************/

package monosat;

import java.nio.IntBuffer;
import java.util.*;
import java.util.stream.Collectors;

public final class Graph {
    Solver solver;
    int bitwidth = -1;
    private long graphPtr;
    private Map<String, Integer> nodeMap = new HashMap<String, Integer>();
    private ArrayList<String> nodeNames = new ArrayList<String>();
    private Set<Integer> nodes = new HashSet<Integer>();//consider arraylist<Integer>
    private ArrayList<Map<Integer, LinkedList<Edge>>> edges = new ArrayList<>();
    private Map<Lit,Edge> edgeLitMap = new HashMap<>();
    private ArrayList<Edge> all_edges = new ArrayList<>();
    private ArrayList<Lit> all_edge_lits = new ArrayList<>();
    private ArrayList<ArrayList<Lit>> all_out_edge_lits = new ArrayList<>();
    private ArrayList<ArrayList<Lit>> all_in_edge_lits = new ArrayList<>();
    private ArrayList<ArrayList<Lit>> all_node_edge_lits = new ArrayList<>();

    public Graph(Solver solver) {
        this.solver = solver;
        graphPtr = MonosatJNI.newGraph(solver.solverPtr);
    }

    public Graph(Solver solver, int bitwidth) {
        this.solver = solver;
        graphPtr = MonosatJNI.newGraph(solver.solverPtr);
        assert (bitwidth >= 0);
        this.bitwidth = bitwidth;
    }

    public Solver getSolver(){
        return solver;
    }

    public int nNodes() {
        return MonosatJNI.nNodes(solver.solverPtr, graphPtr);
    }

    public int nEdges() {
        return MonosatJNI.nEdges(solver.solverPtr, graphPtr);
    }

    /**
     * Gets the name associated with a node
     * @return
     */
    public String getName(int node){
        validateNode(node);
        return nodeNames.get(node);
    }

    public int addNode() {
        return addNode(Integer.toString(nNodes()));
    }

    public int addNode(String name) {
        if (nodeMap.containsKey(name)) {
            throw new RuntimeException("Node names must be unique");
        }
        int n = MonosatJNI.newNode(solver.solverPtr, graphPtr);
        while (edges.size() <= n) {
            edges.add(null);
        }
        assert (edges.get(n) == null);
        edges.set(n, new TreeMap<Integer, LinkedList<Edge>>());//Consider hash map or arraylist here, depending on density of graph...
        nodes.add(n);
        nodeMap.put(name, n);
        nodeNames.add(name);
        while(all_in_edge_lits.size()<=n){
            all_in_edge_lits.add(new ArrayList<>());
            all_out_edge_lits.add(new ArrayList<>());
            all_node_edge_lits.add(new ArrayList<>());
        }
        return n;

    }

    public Set<Integer> nodes() {
        return Collections.unmodifiableSet(nodes);
    }

    public boolean hasNode(String name) {
        return nodeMap.containsKey(name);
    }

    public boolean hasNode(int n) {
        validateNode(n);
        return nodes.contains(n);
    }

    public int getNode(String name) {
        return nodeMap.get(name);
    }

    public boolean hasEdge(int from, int to) {
        assert (hasNode(from));
        assert (hasNode(to));
        LinkedList<Edge> edge_list = edges.get(from).get(to);
        if (edge_list == null) {
            return false;
        } else {
            return !edge_list.isEmpty();
        }
    }

    public boolean hasEdge(String from, String to) {
        return hasEdge(getNode(from), getNode(to));
    }

    public Edge getEdge(int from, int to) {
        if (!hasNode(from)) {
            throw new RuntimeException("No such node");
        }
        if (!hasNode(to)) {
            throw new RuntimeException("No such node");
        }
        LinkedList<Edge> edge_list = edges.get(from).get(to);
        if (edge_list == null || edge_list.isEmpty()) {
            throw new RuntimeException("No such edge");
        } else {
            return edge_list.getFirst();
        }
    }

    public Edge getEdge(String from, String to) {
        return getEdge(getNode(from), getNode(to));
    }

    public Edge getEdge(Lit edgeVar) {
        if(!edgeLitMap.containsKey(edgeVar)){
            throw new RuntimeException("Literal " + edgeVar + " does not correspond to an edge in the graph.");
        }
        return edgeLitMap.get(edgeVar);
    }

    public List<Edge> getAllEdges(int from, int to) {
        if (!hasNode(from)) {
            throw new RuntimeException("No such node");
        }
        if (!hasNode(to)) {
            throw new RuntimeException("No such node");
        }
        LinkedList<Edge> edge_list = edges.get(from).get(to);
        if (edge_list == null || edge_list.isEmpty()) {
            return Collections.<Edge>emptyList();
        } else {
            return Collections.<Edge>unmodifiableList(edge_list);
        }
    }

    public Collection<Edge> getAllEdges(String from, String to) {
        return getAllEdges(getNode(from), getNode(to));
    }
    public Collection<Edge> getAllEdges() {
        return Collections.unmodifiableList(all_edges);
    }
    public Collection<Edge> getAllEdges(int from) {
        if (!hasNode(from)) {
            throw new RuntimeException("No such node");
        }
        Map<Integer, LinkedList<Edge>> edge_map = edges.get(from);
        //is this efficient for large graphs?
        return edge_map.values().stream().flatMap(Collection::stream).collect(Collectors.toList());
    }

    public Collection<Edge> getAllEdges(String from) {
        return getAllEdges(getNode(from));
    }

    /**
     * Get all edge literals in the graph
     * @return
     */
    public Collection<Lit> getEdgeVars(){
        return Collections.unmodifiableList(all_edge_lits);
    }
    /**
     * If e is an edge f->t, and there exists an edge t->f in the graph, return the literal for such an edge.
     * If there is no such edge, throw a RuntimeException
     */
    public Lit getBackEdgeVar(Lit e){
        return getBackEdge(e).l;
    }
    /**
     * If e is an edge f->t, and there exists an edge t->f in the graph, return such an edge.
     * If there is no such edge, throw a RuntimeException
     */
    public Edge getBackEdge(Lit e){
        Edge edge = getEdge(e);
        return getEdge(edge.to,edge.from);
    }
    private void validateNode(int node){
        if (node<0 || node>=nNodes()){
            throw new RuntimeException("Node " + node + " does not exist");
        }
    }
    /**
     * Get all outgoing and all incoming edges from this node
     * @param node
     * @return
     */
    public Collection<Lit> getEdgeVars(int node){
        validateNode(node);

        return Collections.unmodifiableList(all_node_edge_lits.get(node));
    }

    public Collection<Lit> getOutgoingEdgeVars(int node){
        validateNode(node);
        return Collections.unmodifiableList(all_out_edge_lits.get(node));
    }
    public Collection<Lit> getIncomingEdgeVars(int node){
        validateNode(node);
        return Collections.unmodifiableList(all_in_edge_lits.get(node));
    }



    public int bitwidth() {
        return bitwidth;
    }

    public Lit addEdge(int from, int to) {
        return addEdge(from, to, 1);
    }

    public Lit addEdge(int from, int to, long constant_weight) {
        validateNode(from);
        validateNode(to);
        if (bitwidth >= 0) {
            return addEdge(from, to, solver.bv(bitwidth, constant_weight));
        } else {
            Lit l = solver.toLit(MonosatJNI.newEdge(solver.solverPtr, graphPtr, from, to, constant_weight));
            Map<Integer, LinkedList<Edge>> edge_map = edges.get(from);
            if (edge_map.get(to) == null) {
                edge_map.put(to, new LinkedList<>());
            }
            Edge e = new Edge(from, to, l, constant_weight);
            edge_map.get(to).add(e);
            edgeLitMap.put(l,e);
            all_edges.add(e);
            all_edge_lits.add(l);
            all_out_edge_lits.get(from).add(l);
            all_in_edge_lits.get(to).add(l);
            all_node_edge_lits.get(from).add(l);
            all_node_edge_lits.get(to).add(l);
            return l;
        }
    }

    public Lit addEdge(int from, int to, BitVector weight) {
        validateNode(from);
        validateNode(to);
        if (this.bitwidth < 0) {
            throw new RuntimeException("In order to use bitvector edge weights, the bitwidth must be passed to the graph constructor, eg:" +
                    "Graph g = new Graph(solver, 8); will accept edges with bitvectors of size 8. Otherwise, edge weights are assumed to be constant integers.");
        }
        Lit l = solver.toLit(MonosatJNI.newEdge_bv(solver.solverPtr, graphPtr, from, to, weight.id));
        Map<Integer, LinkedList<Edge>> edge_map = edges.get(from);
        if (edge_map.get(to) == null) {
            edge_map.put(to, new LinkedList<>());
        }
        Edge e = new Edge(from, to, l, weight);
        edge_map.get(to).add(e);
        all_edges.add(e);
        edgeLitMap.put(l,e);
        all_edge_lits.add(l);
        all_out_edge_lits.get(from).add(l);
        all_in_edge_lits.get(to).add(l);
        all_node_edge_lits.get(from).add(l);
        all_node_edge_lits.get(to).add(l);
        return l;
    }

    public Lit addUndirectedEdge(int from, int to) {
        return addUndirectedEdge(from, to, 1);
    }

    public Lit addUndirectedEdge(int from, int to, long constant_weight) {
        validateNode(to);
        validateNode(from);
        if (bitwidth >= 0) {
            return addUndirectedEdge(from, to, solver.bv(bitwidth, constant_weight));
        } else {
            Lit l = addEdge(from,to,constant_weight);
            Lit l2 = addEdge(to,from,constant_weight);
            solver.assertEqual(l,l2);
/*            Lit l = solver.toLit(MonosatJNI.newEdge(solver.solverPtr, graphPtr, from, to, constant_weight));
            Lit l2 = solver.toLit(MonosatJNI.newEdge(solver.solverPtr, graphPtr, to, from,constant_weight));
            solver.assertEqual(l,l2);
            Map<Integer, LinkedList<Edge>> edge_map = edges.get(from);
            if (edge_map.get(to) == null) {
                edge_map.put(to, new LinkedList<>());
            }
            Edge e1 = new Edge(from, to, l, constant_weight);
            Edge e2 = new Edge(to, from, l2, constant_weight);
            edge_map.get(to).add(e1);
            all_edges.add(e1);
            allEdgeMap.put(l,e);
            allEdgeMap.put(l2,e);
            all_edge_lits.add(l);
            all_out_edge_lits.get(from).add(l);
            all_in_edge_lits.get(to).add(l);
            all_out_edge_lits.get(to).add(l);
            all_in_edge_lits.get(from).add(l);
            all_node_edge_lits.get(from).add(l);
            all_node_edge_lits.get(to).add(l);*/
            return l;
        }
    }

    public Lit addUndirectedEdge(int from, int to, BitVector weight) {
        validateNode(to);
        validateNode(from);
        if (this.bitwidth < 0) {
            throw new RuntimeException("In order to use bitvector edge weights, the bitwidth must be passed to the graph constructor, eg:" +
                    "Graph g = new Graph(solver, 8); will accept edges with bitvectors of size 8. Otherwise, edge weights are assumed to be constant integers.");
        }
        Lit l = addEdge(from,to,weight);
        Lit l2 = addEdge(to,from,weight);
        solver.assertEqual(l,l2);
        /*Lit l = solver.toLit(MonosatJNI.newEdge_bv(solver.solverPtr, graphPtr, from, to, weight.id));
        Lit l2 = solver.toLit(MonosatJNI.newEdge_bv(solver.solverPtr, graphPtr, to, from, weight.id));
        solver.assertEqual(l,l2);
        Map<Integer, LinkedList<Edge>> edge_map = edges.get(from);
        if (edge_map.get(to) == null) {
            edge_map.put(to, new LinkedList<>());
        }
        Edge e = new Edge(from, to, l, weight);
        edge_map.get(to).add(e);
        all_edges.add(e);
        allEdgeMap.put(l,e);
        all_edge_lits.add(l);
        all_out_edge_lits.get(from).add(l);
        all_in_edge_lits.get(to).add(l);
        all_out_edge_lits.get(to).add(l);
        all_in_edge_lits.get(from).add(l);
        all_node_edge_lits.get(from).add(l);
        all_node_edge_lits.get(to).add(l);*/
        return l;
    }

    /**
     * True if the graph is acyclic, false otherwise.
     * If directed is false, all edges are interpreted as undirected for the purposes of cycle checking.
     *
     * @param directed
     * @return
     */
    public Lit acyclic(boolean directed) {
        if (directed) {
            return solver.toLit(MonosatJNI.acyclic_directed(solver.solverPtr, this.graphPtr));
        } else {
            return solver.toLit(MonosatJNI.acyclic_undirected(solver.solverPtr, this.graphPtr));
        }
    }

    /**
     * True if the graph is (directed) acyclic, false otherwise.
     *
     * @return
     */
    public Lit acyclic() {
        return acyclic(true);
    }

    public Lit reaches(int from, int to) {
        validateNode(from);
        validateNode(to);
        return solver.toLit(MonosatJNI.reaches(solver.solverPtr, graphPtr, from, to));
    }
    /**
    * True if there exists a backwards path from -> to.
    * For example, if there is an edge u->v enabled in the graph,
    * then there is a backwards path from v to u.
    */
    public Lit reachesBackward(int from, int to) {
        validateNode(from);
        validateNode(to);
        return solver.toLit(MonosatJNI.reachesBackward(solver.solverPtr, graphPtr, from, to));
    }
    /**
    * True if there exists a path from 'from' to nodeOnPath, and a path from 'nodeOnPath' to 'to'
    * This constraint is logically equivalent to And(g.reaches(from,nodeOnPath),g.reaches(nodeOnPath,to)),
    * but may be more efficient if lots of nodeOnPath checks are made between the same 'from' and 'to' nodes.
    */
    public Lit onPath(int nodeOnPath,int from, int to) {
        validateNode(from);
        validateNode(to);
        return solver.toLit(MonosatJNI.onPath(solver.solverPtr, graphPtr, nodeOnPath,from, to));
    }
    public Lit compareDistance(int from, int to,Comparison comparison, long compareTo) {
        validateNode(from);
        validateNode(to);
        switch (comparison) {
            case GEQ:
                return solver.toLit(MonosatJNI.shortestPath_lt_const(solver.solverPtr, graphPtr, from, to, compareTo)).not();
            case GT:
                return solver.toLit(MonosatJNI.shortestPath_leq_const(solver.solverPtr, graphPtr, from, to, compareTo)).not();
            case LEQ:
                return solver.toLit(MonosatJNI.shortestPath_leq_const(solver.solverPtr, graphPtr, from, to, compareTo));
            case LT:
                return solver.toLit(MonosatJNI.shortestPath_lt_const(solver.solverPtr, graphPtr, from, to, compareTo));
            case EQ: {
                Lit l1 = solver.toLit(MonosatJNI.shortestPath_leq_const(solver.solverPtr, graphPtr, from, to, compareTo));
                Lit l2 = solver.toLit(MonosatJNI.shortestPath_lt_const(solver.solverPtr, graphPtr, from, to, compareTo));
                return solver.and(l1, l2);
            }
            case NEQ: {
                Lit l1 = solver.toLit(MonosatJNI.shortestPath_leq_const(solver.solverPtr, graphPtr, from, to, compareTo));
                Lit l2 = solver.toLit(MonosatJNI.shortestPath_lt_const(solver.solverPtr, graphPtr, from, to, compareTo));
                return solver.nand(l1, l2);
            }
        }
        throw new RuntimeException("Unknown comparison");
    }

    public Lit compareDistance(int from, int to,Comparison comparison, BitVector compareTo) {
        validateNode(from);
        validateNode(to);
        switch (comparison) {
            case GEQ:
                return solver.toLit(MonosatJNI.shortestPath_lt_bv(solver.solverPtr, graphPtr, from, to, compareTo.id)).not();
            case GT:
                return solver.toLit(MonosatJNI.shortestPath_leq_bv(solver.solverPtr, graphPtr, from, to, compareTo.id)).not();
            case LEQ:
                return solver.toLit(MonosatJNI.shortestPath_leq_bv(solver.solverPtr, graphPtr, from, to, compareTo.id));
            case LT:
                return solver.toLit(MonosatJNI.shortestPath_lt_bv(solver.solverPtr, graphPtr, from, to, compareTo.id));
            case EQ: {
                Lit l1 = solver.toLit(MonosatJNI.shortestPath_leq_bv(solver.solverPtr, graphPtr, from, to, compareTo.id));
                Lit l2 = solver.toLit(MonosatJNI.shortestPath_lt_bv(solver.solverPtr, graphPtr, from, to, compareTo.id));
                return solver.and(l1, l2);
            }
            case NEQ: {
                Lit l1 = solver.toLit(MonosatJNI.shortestPath_leq_bv(solver.solverPtr, graphPtr, from, to, compareTo.id));
                Lit l2 = solver.toLit(MonosatJNI.shortestPath_lt_bv(solver.solverPtr, graphPtr, from, to, compareTo.id));
                return solver.nand(l1, l2);
            }
        }
        throw new RuntimeException("Unknown comparison");
    }

    public BitVector distance(int from, int to) {
        return distance(from, to, this.bitwidth);
    }

    public BitVector distance(int from, int to, int bitwidth) {
        validateNode(from);
        validateNode(to);
        if (bitwidth < 0) {
            //compute a 'large enough' bitwidth to hold the maximum possible distance
            bitwidth = 8;//temporary
        }

        BitVector result = new BitVector(solver, bitwidth);

        Lit l1 = solver.toLit(MonosatJNI.shortestPath_leq_bv(solver.solverPtr, graphPtr, from, to, result.id));
        Lit l2 = solver.toLit(MonosatJNI.shortestPath_lt_bv(solver.solverPtr, graphPtr, from, to, result.id));
        //result is geq to the shortest path, and is not greater than the shortest path, and so it is exactly equal to the shortestPath
        solver.assertTrue(l1);
        solver.assertFalse(l2);
        return result;
    }

    /**
     * Comparison a constant or bitvector to the maximum flow of a graph.
     * Note that if your goal is to assert or reason about a one-sided comparison between the maximum flow and a constant or bitvector,
     * the compareMaximumFlow form is more efficient then the direct maximumFlow() form.
     *
     * @param from
     * @param to
     * @return
     */
    public Lit compareMaximumFlow(int from, int to,Comparison comparison, long compareTo) {
        validateNode(from);
        validateNode(to);
        switch (comparison) {
            case GEQ:
                return solver.toLit(MonosatJNI.maximumFlow_geq(solver.solverPtr, graphPtr, from, to, compareTo));
            case GT:
                return solver.toLit(MonosatJNI.maximumFlow_gt(solver.solverPtr, graphPtr, from, to, compareTo));
            case LEQ:
                return solver.toLit(MonosatJNI.maximumFlow_gt(solver.solverPtr, graphPtr, from, to, compareTo)).not();
            case LT:
                return solver.toLit(MonosatJNI.maximumFlow_geq(solver.solverPtr, graphPtr, from, to, compareTo)).not();
            case EQ: {
                Lit l1 = solver.toLit(MonosatJNI.maximumFlow_geq(solver.solverPtr, graphPtr, from, to, compareTo));
                Lit l2 = solver.toLit(MonosatJNI.maximumFlow_gt(solver.solverPtr, graphPtr, from, to, compareTo));
                return solver.and(l1, l2);
            }
            case NEQ: {
                Lit l1 = solver.toLit(MonosatJNI.maximumFlow_geq(solver.solverPtr, graphPtr, from, to, compareTo));
                Lit l2 = solver.toLit(MonosatJNI.maximumFlow_gt(solver.solverPtr, graphPtr, from, to, compareTo));
                return solver.nand(l1, l2);
            }
        }
        throw new RuntimeException("Unknown comparison");
    }

    /**
     * Comparison a constant or bitvector to the maximum flow of a graph.
     * Note that if your goal is to assert or reason about a one-sided comparison between the maximum flow and a constant or bitvector,
     * the compareMaximumFlow form is more efficient then the direct maximumFlow() form.
     *
     * @param from
     * @param to
     * @return
     */
    public Lit compareMaximumFlow(int from, int to, Comparison comparison, BitVector compareTo) {
        validateNode(from);
        validateNode(to);
        switch (comparison) {
            case GEQ:
                return solver.toLit(MonosatJNI.maximumFlow_geq_bv(solver.solverPtr, graphPtr, from, to, compareTo.id));
            case GT:
                return solver.toLit(MonosatJNI.maximumFlow_gt_bv(solver.solverPtr, graphPtr, from, to, compareTo.id));
            case LEQ:
                return solver.toLit(MonosatJNI.maximumFlow_gt_bv(solver.solverPtr, graphPtr, from, to, compareTo.id)).not();
            case LT:
                return solver.toLit(MonosatJNI.maximumFlow_geq_bv(solver.solverPtr, graphPtr, from, to, compareTo.id)).not();
            case EQ: {
                Lit l1 = solver.toLit(MonosatJNI.maximumFlow_geq_bv(solver.solverPtr, graphPtr, from, to, compareTo.id));
                Lit l2 = solver.toLit(MonosatJNI.maximumFlow_gt_bv(solver.solverPtr, graphPtr, from, to, compareTo.id));
                return solver.and(l1, l2);
            }
            case NEQ: {
                Lit l1 = solver.toLit(MonosatJNI.maximumFlow_geq_bv(solver.solverPtr, graphPtr, from, to, compareTo.id));
                Lit l2 = solver.toLit(MonosatJNI.maximumFlow_gt_bv(solver.solverPtr, graphPtr, from, to, compareTo.id));
                return solver.nand(l1, l2);
            }
        }
        throw new RuntimeException("Unknown comparison");
    }

    public BitVector maximumFlow(int from, int to) {
        return maximumFlow(from, to, this.bitwidth);
    }

    /**
     * Create a new bitvector, equal in value to the maximum flow between from and to.
     * Note that if your goal is to  assert or reason about a one-sided comparison between the maximum flow and a constant or bitvector,
     * the compareMaximumFlow form is more efficient then the direct maximumFlow() form.
     *
     * @param from
     * @param to
     * @return
     */
    public BitVector maximumFlow(int from, int to, int bitwidth) {
        validateNode(from);
        validateNode(to);
        if (bitwidth < 0) {
            //compute a 'large enough' bitwidth to hold the maximum possible flow
            bitwidth = 8;//temporary
        }
        BitVector result = new BitVector(solver, bitwidth);
        Lit l1 = solver.toLit(MonosatJNI.maximumFlow_geq_bv(solver.solverPtr, graphPtr, from, to, result.id));
        Lit l2 = solver.toLit(MonosatJNI.maximumFlow_gt_bv(solver.solverPtr, graphPtr, from, to, result.id));
        //result is geq to the max flow, and is not greater than the max flow (and so it is exactly equal to the maxflow)
        solver.assertTrue(l1);
        solver.assertFalse(l2);
        return result;
    }

    public ArrayList<Integer> getPathNodes(Lit reach_or_distance_literal) {
        getSolver().validate(reach_or_distance_literal);
        ArrayList<Integer> store = new ArrayList<>();
        if (!MonosatJNI.hasModel(solver.solverPtr)) {
            throw new RuntimeException("Solver has no model (this may indicate either that the solve() has not yet been called, or that the most recent call to solve() returned a value other than true, or that a constraint was added into the solver after the last call to solve()).");
        }
        int length = MonosatJNI.getModel_Path_Nodes_Length(solver.solverPtr, graphPtr, reach_or_distance_literal.l);
        if (length < 0) {
            throw new RuntimeException("No path (this may indicate the solver has no model)");
        }
        IntBuffer buf = solver.getBuffer(0, length);
        int sz = MonosatJNI.getModel_Path_Nodes(solver.solverPtr, graphPtr, reach_or_distance_literal.l, length, buf);
        assert (sz == length);
        for (int i = 0; i < length; i++) {
            int n = buf.get(i);
            store.add(n);
        }
        return store;
    }

    public ArrayList<Lit> getPathEdges(Lit reach_or_distance_literal) {
        getSolver().validate(reach_or_distance_literal);
        if (!MonosatJNI.hasModel(solver.solverPtr)) {
            throw new RuntimeException("Solver has no model (this may indicate either that the solve() has not yet been called, or that the most recent call to solve() returned a value other than true, or that a constraint was added into the solver after the last call to solve()).");
        }
        int length = MonosatJNI.getModel_Path_EdgeLits_Length(solver.solverPtr, graphPtr, reach_or_distance_literal.l);
        if (length < 0) {
            throw new RuntimeException("No path (this may indicate the solver has no model)");
        }
        ArrayList<Lit> store = new ArrayList<>();
        IntBuffer buf = solver.getBuffer(0, length);
        int sz = MonosatJNI.getModel_Path_EdgeLits(solver.solverPtr, graphPtr, reach_or_distance_literal.l, length, buf);
        assert (sz == length);
        for (int i = 0; i < length; i++) {
            int l = buf.get(i);
            Lit lit = solver.toLit(l);
            store.add(lit);
        }
        return store;
    }



    public long getMaxFlow(Lit flowComparisonLiteral){
        return MonosatJNI.getModel_MaxFlow(solver.solverPtr,graphPtr, flowComparisonLiteral.l);
    }

    public long getEdgeFlow(Lit flowComparisonLiteral, Lit edgeLiteral){
        return getEdgeFlow(flowComparisonLiteral,edgeLiteral,false);
    }
    public long getEdgeFlow(Lit flowComparisonLiteral, Lit edgeLiteral, boolean forceAcyclicFlow){
        if(forceAcyclicFlow) {
            return MonosatJNI.getModel_AcyclicEdgeFlow(solver.solverPtr, graphPtr, flowComparisonLiteral.l,edgeLiteral.l);
        }else{
            return MonosatJNI.getModel_EdgeFlow(solver.solverPtr, graphPtr, flowComparisonLiteral.l,edgeLiteral.l);
        }
    }

    /**
     * Prints out a graphviz/.dot formatted representation of this graph
     * @return
     */
    public String draw(boolean showModel, boolean showConstants){
        if(showModel && !solver.hasModel()) {
            showModel= false;//there is no model, show all edges
        }
        StringBuilder writer = new StringBuilder();
        writer.append("digraph{\n");
        for(int n= 0;n<nNodes();n++){
            writer.append("n"+Integer.toString(n) + " [label=\""+getName(n) + "\"]\n");
        }


        for (Edge e:getAllEdges()){
            writer.append("n"+Integer.toString(e.from) + "->" + "n"+Integer.toString(e.to) + "[label=\""+e.l.toString() + "\"");
            if(showModel || showConstants){
                Optional<Boolean> possibleValue = solver.getPossibleValue(e.l);
                long weight = e.getWeight();
                if (possibleValue.isPresent() && !possibleValue.get()){
                    writer.append(",color=white");
                }else if (!possibleValue.isPresent() && showModel){
                    writer.append(",color=gray,style=dashed");
                }else{
                    writer.append(",color=black");
                    if(weight!=1){
                        //supress weight 1, as that is the default
                        writer.append(",label=\"" + Long.toString(weight) + "\"");
                    }
                }
            }else{

            }
            writer.append("]\n");
        }
        writer.append("}\n");
        return writer.toString();
    }

    public String draw(){
        return draw(true,true);
    }
    public String draw(boolean showModel){
        return draw(showModel,true);
    }

    public class Edge {
        public int from;
        public int to;
        public Lit l;
        private long weight = -1;
        private BitVector bv = null;

        protected Edge(int from, int to, Lit l, BitVector weight) {
            this.from = from;
            this.to = to;
            this.l = l;
            this.bv = weight;
            this.weight = -1;
        }

        protected Edge(int from, int to, Lit l, long weight) {
            this.from = from;
            this.to = to;
            this.l = l;
            this.bv = null;
            this.weight = weight;
        }

        public boolean hasBV() {
            return bv != null;
        }

        public BitVector getBV() {
            if (!hasBV()) {
                throw new RuntimeException("Attempt to access bitvector weight on edge with constant weight");
            }
            return bv;
        }
        public long getWeight() {
            if (hasBV()) {
                return getBV().value();
                //throw new RuntimeException("Attempt to access constant weight on edge with bitvector weight");
            }
            assert (weight >= 0);
            return weight;
        }
    }
}
