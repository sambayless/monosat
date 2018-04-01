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

public class Graph {
    Solver solver;
    int bitwidth = -1;
    private long graphPtr;
    private Map<String, Integer> nodeMap = new HashMap<String, Integer>();
    private Set<Integer> nodes = new HashSet<Integer>();//consider arraylist<Integer>
    private ArrayList<Map<Integer, LinkedList<Edge>>> edges = new ArrayList<>();

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

    public int nNodes() {
        return MonosatJNI.nNodes(solver.solverPtr, graphPtr);
    }

    public int nEdges() {
        return MonosatJNI.nEdges(solver.solverPtr, graphPtr);
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
        return n;

    }

    public Set<Integer> nodes() {
        return Collections.unmodifiableSet(nodes);
    }

    public boolean hasNode(String name) {
        return nodeMap.containsKey(name);
    }

    public boolean hasNode(Integer n) {
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

    public int bitwidth() {
        return bitwidth;
    }

    public Lit addEdge(int from, int to) {
        return addEdge(from, to, 1);
    }

    public Lit addEdge(int from, int to, long constant_weight) {
        if (bitwidth >= 0) {
            return addEdge(from, to, solver.bv(bitwidth, constant_weight));
        } else {
            Lit l = solver.toLit(MonosatJNI.newEdge(solver.solverPtr, graphPtr, from, to, constant_weight));
            Map<Integer, LinkedList<Edge>> edge_map = edges.get(from);
            if (edge_map.get(to) == null) {
                edge_map.put(to, new LinkedList<>());
            }
            edge_map.get(to).add(new Edge(from, to, l, constant_weight));
            return l;
        }
    }

    public Lit addEdge(int from, int to, BitVector weight) {
        if (this.bitwidth < 0) {
            throw new RuntimeException("In order to use bitvector edge weights, the bitwidth must be passed to the graph constructor, eg:" +
                    "Graph g = new Graph(solver, 8); will accept edges with bitvectors of size 8. Otherwise, edge weights are assumed to be constant integers.");
        }
        Lit l = solver.toLit(MonosatJNI.newEdge_bv(solver.solverPtr, graphPtr, from, to, weight.id));
        Map<Integer, LinkedList<Edge>> edge_map = edges.get(from);
        if (edge_map.get(to) == null) {
            edge_map.put(to, new LinkedList<>());
        }
        edge_map.get(to).add(new Edge(from, to, l, weight));
        return l;
    }

    public Lit addUndirectedEdge(int from, int to) {
        return addUndirectedEdge(from, to, 1);
    }

    public Lit addUndirectedEdge(int from, int to, long constant_weight) {
        if (bitwidth >= 0) {
            return addUndirectedEdge(from, to, solver.bv(bitwidth, constant_weight));
        } else {
            Lit l = solver.toLit(MonosatJNI.newEdge(solver.solverPtr, graphPtr, from, to, constant_weight));
            Lit l2 = solver.toLit(MonosatJNI.newEdge(solver.solverPtr, graphPtr, to, from,constant_weight));
            solver.assertEqual(l,l2);
            Map<Integer, LinkedList<Edge>> edge_map = edges.get(from);
            if (edge_map.get(to) == null) {
                edge_map.put(to, new LinkedList<>());
            }
            edge_map.get(to).add(new Edge(from, to, l, constant_weight));
            return l;
        }
    }

    public Lit addUndirectedEdge(int from, int to, BitVector weight) {
        if (this.bitwidth < 0) {
            throw new RuntimeException("In order to use bitvector edge weights, the bitwidth must be passed to the graph constructor, eg:" +
                    "Graph g = new Graph(solver, 8); will accept edges with bitvectors of size 8. Otherwise, edge weights are assumed to be constant integers.");
        }
        Lit l = solver.toLit(MonosatJNI.newEdge_bv(solver.solverPtr, graphPtr, from, to, weight.id));
        Lit l2 = solver.toLit(MonosatJNI.newEdge_bv(solver.solverPtr, graphPtr, to, from, weight.id));
        solver.assertEqual(l,l2);
        Map<Integer, LinkedList<Edge>> edge_map = edges.get(from);
        if (edge_map.get(to) == null) {
            edge_map.put(to, new LinkedList<>());
        }
        edge_map.get(to).add(new Edge(from, to, l, weight));
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
        return solver.toLit(MonosatJNI.reaches(solver.solverPtr, graphPtr, from, to));
    }

    public Lit compareDistance(int from, int to, long compareTo, Compare comparison) {
        switch (comparison) {
            case GEQ:
                return solver.toLit(MonosatJNI.shortestPath_lt_const(solver.solverPtr, graphPtr, from, to, compareTo)).negate();
            case GT:
                return solver.toLit(MonosatJNI.shortestPath_leq_const(solver.solverPtr, graphPtr, from, to, compareTo)).negate();
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

    public Lit compareDistance(int from, int to, BitVector compareTo, Compare comparison) {
        switch (comparison) {
            case GEQ:
                return solver.toLit(MonosatJNI.shortestPath_lt_bv(solver.solverPtr, graphPtr, from, to, compareTo.id)).negate();
            case GT:
                return solver.toLit(MonosatJNI.shortestPath_leq_bv(solver.solverPtr, graphPtr, from, to, compareTo.id)).negate();
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
     * Compare a constant or bitvector to the maximum flow of a graph.
     * Note that if your goal is to assert or reason about a one-sided comparison between the maximum flow and a constant or bitvector,
     * the compareMaximumFlow form is more efficient then the direct maximumFlow() form.
     *
     * @param from
     * @param to
     * @return
     */
    public Lit compareMaximumFlow(int from, int to, long compareTo, Compare comparison) {
        switch (comparison) {
            case GEQ:
                return solver.toLit(MonosatJNI.maximumFlow_geq(solver.solverPtr, graphPtr, from, to, compareTo));
            case GT:
                return solver.toLit(MonosatJNI.maximumFlow_gt(solver.solverPtr, graphPtr, from, to, compareTo));
            case LEQ:
                return solver.toLit(MonosatJNI.maximumFlow_gt(solver.solverPtr, graphPtr, from, to, compareTo)).negate();
            case LT:
                return solver.toLit(MonosatJNI.maximumFlow_geq(solver.solverPtr, graphPtr, from, to, compareTo)).negate();
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
     * Compare a constant or bitvector to the maximum flow of a graph.
     * Note that if your goal is to assert or reason about a one-sided comparison between the maximum flow and a constant or bitvector,
     * the compareMaximumFlow form is more efficient then the direct maximumFlow() form.
     *
     * @param from
     * @param to
     * @return
     */
    public Lit compareMaximumFlow(int from, int to, BitVector compareTo, Compare comparison) {
        switch (comparison) {
            case GEQ:
                return solver.toLit(MonosatJNI.maximumFlow_geq_bv(solver.solverPtr, graphPtr, from, to, compareTo.id));
            case GT:
                return solver.toLit(MonosatJNI.maximumFlow_gt_bv(solver.solverPtr, graphPtr, from, to, compareTo.id));
            case LEQ:
                return solver.toLit(MonosatJNI.maximumFlow_gt_bv(solver.solverPtr, graphPtr, from, to, compareTo.id)).negate();
            case LT:
                return solver.toLit(MonosatJNI.maximumFlow_geq_bv(solver.solverPtr, graphPtr, from, to, compareTo.id)).negate();
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
        reach_or_distance_literal.validate();
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
        reach_or_distance_literal.validate();
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
                throw new RuntimeException("Attempt to access constant weight on edge with bitvector weight");
            }
            assert (weight >= 0);
            return weight;
        }
    }
}
