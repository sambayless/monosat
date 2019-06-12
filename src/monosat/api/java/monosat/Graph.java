/*
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
*/

package monosat;

import java.nio.IntBuffer;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Each Solver may have one or more associated Graphs. Graphs may have directed or undirected edges,
 * with edges that may be unweighted, or weighted by constants (longs) or by BitVectors. Graphs may
 * have cycles, and may also have multiple edges between the same nodes (eg, the graph may be a
 * multigraph). <br>
 * Example usage:
 *
 * <blockquote>
 *
 * <pre>{@code
 * Solver solver = new Solver();
 * //Create a new graph in solver
 * Graph graph = new Graph(solver);
 *
 * //create 3 nodes in the graph
 * int n0 = graph.addNode();
 * int n1 = graph.addNode();
 * int n2 = graph.addNode();
 *
 * //Create 3 edges, which will be included in the graph
 * //iff their corresponding literals are assigned true by solver.
 * Lit e0 = graph.addEdge(n0, n1);
 * Lit e1 = graph.addEdge(n0, n2);
 * Lit e2 = graph.addEdge(n1, n2);
 *
 * //You can use arbitrary Boolean logic constraints to restrict which combinations of edges
 * //can be enabled in the solver.
 * solver.addClause(e0, e1.not());
 * solver.addClause(e0.not(), e2);
 *
 * //r is true iff there is a path from node '0' to node '3', in g.
 * Lit r = graph.reaches(0, 2);
 * //Assert that there must be a path from node 0 to node 2 under the assignment selected by the solver.
 * solver.addClause(r);
 * }</pre>
 *
 * </blockquote>
 *
 * Many common graph predicates are supported:
 *
 * <ul>
 *   <li>@see #Graph.reaches(int,int)
 *   <li>@see #Graph.distance()
 *   <li>@see #Graph.maxflow()
 *   <li>@see #Graph.onPath()
 *   <li>@see #Graph.acyclic()
 * </ul>
 */
public final class Graph {
  /** The solver instance this graph belongs to. */
  private final Solver solver;
  /** A handle to this graph in the native library. */
  private final long graphPtr;
  /**
   * If this graph uses bitvector edges, then those bitvectors must all have width = bitwidth. If
   * bitwidth is <0, then this graph uses constant longs as edge weights.
   */
  private final int bitwidth;
  /** Map from node name strings, to the integers representing nodes in the graph. */
  private final Map<String, Integer> nodeMap = new HashMap<String, Integer>();
  /** Iterable list of names for each node. */
  private final ArrayList<String> nodeNames = new ArrayList<String>();
  /** Iterable list of nodes in this graph. */
  private final ArrayList<Integer> nodes = new ArrayList<Integer>();
  /** An adjacency list representation fo the edges in this graph. */
  private final ArrayList<Map<Integer, LinkedList<Edge>>> adjacencyList = new ArrayList<>();
  /** A map from the literals controlling each edge, to their corresponding edge. */
  private final Map<Lit, Edge> edgeLitMap = new HashMap<>();
  /** An iterable list of the edges in this graph. */
  private final ArrayList<Edge> all_edges = new ArrayList<>();
  /** An iterable list of the edge literals in this graph. */
  private final ArrayList<Lit> all_edge_lits = new ArrayList<>();
  /** An adjacency list representation of just the outgoing edge literals in this graph. */
  private final ArrayList<ArrayList<Lit>> all_out_edge_lits = new ArrayList<>();
  /** An adjacency list representation of just the incoming edge literals in this graph. */
  private final ArrayList<ArrayList<Lit>> all_in_edge_lits = new ArrayList<>();
  /**
   * For each node in the graph, a list of all the incoming and outgoing edge literals connecting to
   * that node.
   */
  private final ArrayList<ArrayList<Lit>> all_node_edge_lits = new ArrayList<>();

  private final String _name;

  public String toString() {
    if (this._name.length() > 0) {
      return _name + "(Graph)";
    } else {
      return "Graph";
    }
  }


  /**
   * Instantiate a graph in solver. This graph can have unweighted edges, or constant (long)
   * weighted edges. To create a graph with variable (BitVector) weighted edges @see #Graph(Solver,
   * int).
   *
   * @param solver The solver instance in which to create this graph.
   */
  public Graph(Solver solver) {
    this(solver, "");
  }

  /**
   * Instantiate a graph in solver. This graph can have unweighted edges, or constant (long)
   * weighted edges. To create a graph with variable (BitVector) weighted edges @see #Graph(Solver,
   * int).
   *
   * @param solver The solver instance in which to create this graph.
   * @param name A name for the graph. If empty, the graph will be unnamed. Otherwise, the name must
   *     be unique.
   */
  public Graph(Solver solver, String name) {
    this(solver, -1, name, false);
  }

  /**
   * Instantiate a graph in solver, with BitVector weighted edges. Each edge in this graph must be
   * weighted by a BitVector of width 'bitwidth'. To create a graph with unweighted or constant
   * weighted edges @see #Graph(Solver).
   *
   * @param solver The solver instance in which to create this graph.
   */
  public Graph(Solver solver, int bitwidth) {
    this(solver, bitwidth, "");
  }

  /**
   * Instantiate a graph in solver, with BitVector weighted edges. Each edge in this graph must be
   * weighted by a BitVector of width 'bitwidth'. To create a graph with unweighted or constant
   * weighted edges @see #Graph(Solver).
   *
   * @param solver The solver instance in which to create this graph.
   * @param name A name for the graph. If empty, the graph will be unnamed. Otherwise, the name must
   *     be unique.
   */
  public Graph(Solver solver, int bitwidth, String name) {
    this(solver, bitwidth, name, true);
  }

  /**
   * For internal use only. Reconstruct a graph from an existing graph pointer in the solver.
   *
   * @param solver
   * @param graphPtr
   */
  protected Graph(Solver solver, long graphPtr) {
    this.solver = solver;
    this.graphPtr = graphPtr;
    this._name = MonosatJNI.getGraphName(solver.getSolverPtr(), graphPtr);
    solver.allGraphs.put(graphPtr, this);

    this.bitwidth = MonosatJNI.getGraphWidth(solver.getSolverPtr(), graphPtr);

    for (int n = 0; n < nNodes(); n++) {
      addNode(MonosatJNI.getNodeName(solver.getSolverPtr(), graphPtr, n), n);
    }
    //create edge objects for each edge in the solver.
    //it might be better to do this lazilly
    for(int n = 0;n<nEdges();n++){

      Lit l = solver.toLit(MonosatJNI.getEdgeLiteralN(solver.getSolverPtr(), graphPtr,n));
      int from = MonosatJNI.getEdge_from(solver.getSolverPtr(), graphPtr,l.toInt());
      int to = MonosatJNI.getEdge_to(solver.getSolverPtr(), graphPtr,l.toInt());
      boolean hasBV = MonosatJNI.edgeHasBVWeight(solver.getSolverPtr(), graphPtr,l.toInt());
      Edge e;
      if(hasBV){
        int bvID = MonosatJNI.getEdge_weight_bv(solver.getSolverPtr(), graphPtr,l.toInt());
        BitVector bv = solver.toBitVector(bvID);
        e = new Edge(from, to, l, bv);
      }else{
        long weight = MonosatJNI.getEdge_weight_const(solver.getSolverPtr(), graphPtr,l.toInt());
        e = new Edge(from, to, l, weight);
      }

      Map<Integer, LinkedList<Edge>> edge_map = adjacencyList.get(from);
      if (edge_map.get(to) == null) {
        edge_map.put(to, new LinkedList<>());
      }

      edge_map.get(to).add(e);
      edgeLitMap.put(l, e);
      all_edges.add(e);
      all_edge_lits.add(l);
      all_out_edge_lits.get(from).add(l);
      all_in_edge_lits.get(to).add(l);
      all_node_edge_lits.get(from).add(l);
      all_node_edge_lits.get(to).add(l);
    }
  }

  /**
   * Instantiate a graph in solver, with either weighted or unweighted edges. Each edge in this
   *
   * @param solver The solver instance in which to create this graph.
   * @param bitwidth The width of the bitvectors on the edges, or -1 if weighted edges are not used.
   * @param name A name for the graph. If empty, the graph will be unnamed. Otherwise, the name must
   *     be unique.
   * @weighted_edges Whether or not to use weighted edges.
   */
  private Graph(Solver solver, int bitwidth, String name, boolean weighted_edges) {
    this.solver = solver;

    this.bitwidth = bitwidth;

    if (weighted_edges) {
      if (bitwidth <= 0) {
        throw new IllegalArgumentException(
            "Graphs may either have unweighted edges, or have bit-vector weighted "
                + "edges with bit-width >= 0");
      } else if (bitwidth > 64) {
        throw new IllegalArgumentException(
            "Graphs may either have unweighted edges, or have bit-vector weighted "
                + "edges with bit-width <= 64");
      }
    } else {
      assert (bitwidth == -1);
    }

    if (name != null && !name.isEmpty()) {
      if (name == "True") {
        throw new IllegalArgumentException("Only the built-in True literal may be named \"True\"");
      } else if (name == "False") {
        throw new IllegalArgumentException(
            "Only the built-in False literal may be named \"False\"");
      }

      this._name = name;
      /*
      Check that the string contains only printable ascii characters
       */

      graphPtr = MonosatJNI.newGraph_Named(solver.getSolverPtr(), MonosatJNI.validID(name), bitwidth);

    } else {
      this._name = "";
      graphPtr = MonosatJNI.newGraph_Named(solver.getSolverPtr(), "", bitwidth);
    }
    solver.allGraphs.put(graphPtr, this);
  }

  /**
   * Get the solver that this BitVector belongs to.
   *
   * @return The solver that this BitVector belongs to.
   */
  public Solver getSolver() {
    return solver;
  }

  public String name() {
    return _name;
  }

  /**
   * Get the number of nodes in this graph.
   *
   * @return The number of nodes in this graph.
   */
  public int nNodes() {
    return MonosatJNI.nNodes(solver.getSolverPtr(), graphPtr);
  }

  /**
   * Get the number of edges in this graph.
   *
   * @return The number of edges in this graph.
   */
  public int nEdges() {
    return MonosatJNI.nEdges(solver.getSolverPtr(), graphPtr);
  }

  /**
   * Gets the name associated with a node.
   *
   * @return The name associated with a node.
   */
  public String getName(int node) {
    validateNode(node);
    if(node<nodeNames.size() && nodeNames.get(node).length()>0){
        return nodeNames.get(node);
    }else{
        return Integer.toString(node);
    }

  }

  /**
   * Create a new node. The name of this node will be set to <code>Integer.toString(nNodes())</code>
   *
   * @return A new node, represented by an integer.
   */
  public int addNode() {
    return addNode("");
  }

  /**
   * Create a new node, with the given name.
   *
   * @return A new node, represented by an integer.
   */
  public int addNode(String name) {
    return addNode(name, -1);
  }

  /**
   * Create a new node, with the given name. The name may either be empty, or it must be unique (within this graph).
   *
   *
   * @param nodeID If the node already exists in the solver, specify its ID here. Otherwise, specify
   *     -1.
   * @return A new node, represented by an integer.
   */
  private int addNode(String name, int nodeID) {

    if (nodeID < 0) {
        name = MonosatJNI.validID(name);

        if (name.length() > 0 && MonosatJNI.hasNamedNode(solver.getSolverPtr(), graphPtr,name)) {
            throw new RuntimeException("Node names must be unique");
        }
      if (name.length() > 0) {
        nodeID = MonosatJNI.newNode_Named(solver.getSolverPtr(), graphPtr,name);
      }else{
          nodeID = MonosatJNI.newNode(solver.getSolverPtr(), graphPtr);
      }
    }
    while (adjacencyList.size() <= nodeID) {
      adjacencyList.add(null);
    }
    assert (adjacencyList.get(nodeID) == null);
    // Consider hash map or ArrayList here, depending on
    adjacencyList.set(nodeID, new TreeMap<Integer, LinkedList<Edge>>());
    // density of graph...
    nodes.add(nodeID);
    if (name.length() > 0) {
      nodeMap.put(name, nodeID);
      nodeNames.add(name);
    }
    while (all_in_edge_lits.size() <= nodeID) {
      all_in_edge_lits.add(new ArrayList<>());
      all_out_edge_lits.add(new ArrayList<>());
      all_node_edge_lits.add(new ArrayList<>());
    }
    return nodeID;
  }

  /**
   * Get an (unmodifiable) view of the nodes in the graph. Note that nodes are guaranteed to be
   * represented by consecutive integers starting at 0, so it is safe and efficient to iterate
   * through the nodes in a graph using a for loop, eg:
   *
   * <blockquote>
   *
   * <pre>{@code
   * for (int n = 0;n<graph.nNodes();n++){
   *     System.out.println(graph.getName(n));
   * }
   * }</pre>
   *
   * </blockquote>
   *
   * @return An iterable list of nodes in this solver.
   */
  public List<Integer> nodes() {
    return Collections.unmodifiableList(nodes);
  }

  /**
   * Check if a node with the given string name exists in the graph.
   *
   * @param name The string to check for.
   * @return True if a node with this name exists in the graph.
   */
  public boolean hasNode(String name) {
    return nodeMap.containsKey(name);
  }

  /**
   * Check if the graph contains node 'n'.
   *
   * @param n The node to check for, represented by an integer.
   * @return True if node 'n' exists in this graph.
   */
  public boolean hasNode(int n) {
    validateNode(n);
    return n >= 0 && n < nNodes();
  }

  /**
   * Get a node with the given name, if one exists in the graph.
   *
   * @param name The string to check for.
   * @return The integer representing this node, if this node exists in the graph.
   * throws IllegalArgumentException if the node does not exist in the graph
   */
  public int getNode(String name) {
    if (nodeMap.containsKey(name)) {
      return nodeMap.get(name);
    }else{
        throw new IllegalArgumentException("No node " + name + " in graph " + toString());
    }
  }

  /**
   * Check if there is a directed edge from the node 'from' to the node 'to'.
   *
   * @param from The source node.
   * @param to The destination node.
   * @return True if there is a directed edge from the node 'from' to the node 'to'.
   */
  public boolean hasEdge(int from, int to) {
    assert (hasNode(from));
    assert (hasNode(to));
    LinkedList<Edge> edge_list = adjacencyList.get(from).get(to);
    if (edge_list == null) {
      return false;
    } else {
      return !edge_list.isEmpty();
    }
  }

  /**
   * Check if there is a directed edge from the node 'from' to the node 'to'.
   *
   * @param from The name of the source node.
   * @param to The name of the destination node.
   * @return True if there is a directed edge from the node 'from' to the node 'to'.
   */
  public boolean hasEdge(String from, String to) {
    return hasEdge(getNode(from), getNode(to));
  }

  /**
   * If there is a directed edge from the node 'from' to the node 'to', return that node (otherwise,
   * throw an exception). If there are multiple edges between 'from' and 'to', return an arbitrary
   * edge from that set.
   *
   * @param from The source node.
   * @param to The destination node.
   * @return An edge, if there is a directed edge from the node 'from' to the node 'to'.
   * @throws RuntimeException If there is no directed edge from 'from' to 'to'.
   */
  public Edge getEdge(int from, int to) {
    if (!hasNode(from)) {
      throw new RuntimeException("No such node");
    }
    if (!hasNode(to)) {
      throw new RuntimeException("No such node");
    }
    LinkedList<Edge> edge_list = adjacencyList.get(from).get(to);
    if (edge_list == null || edge_list.isEmpty()) {
      throw new RuntimeException("No such edge");
    } else {
      return edge_list.getFirst();
    }
  }

  /**
   * If there is a directed edge with the given name, return it.
   *
   * @param name The name of the edge to check for.
   * @return An edge, if there is a directed edge with the given name in this graph.
   * @throws RuntimeException If there is no directed edge with this name.
   */
  public Edge getEdge(String name) {
    Lit l = solver.getLiteral(name);
    return getEdge(l);
  }

  /**
   * If there is a directed edge from the node 'from' to the node 'to', return one such node
   * (otherwise, throw an exception). If there are multiple edges between 'from' and 'to', return an
   * arbitrary edge from that set.
   *
   * @param from The name of the source node.
   * @param to The name of the destination node.
   * @return True if there is a directed edge from the node 'from' to the node 'to'.
   * @throws RuntimeException If there is no directed edge from 'from' to 'to'.
   */
  public Edge getEdge(String from, String to) {
    return getEdge(getNode(from), getNode(to));
  }

  /**
   * Get the edge controlled by literal 'edgeVar'.
   *
   * @param edgeVar The literal associated with the edge to be returned.
   * @return The edge controlled by literal 'edgeVar'.
   */
  public Edge getEdge(Lit edgeVar) {
    if (!edgeLitMap.containsKey(edgeVar)) {
      throw new RuntimeException(
          "Literal " + edgeVar + " does not correspond to an edge in the graph.");
    }
    return edgeLitMap.get(edgeVar);
  }

  /**
   * A graph may have multiple directed edges between the same two nodes. This method returns a
   * (possibly empty) list of all edges from 'from' to 'to'.
   *
   * @param from The source node.
   * @param to The destination node.
   * @return A list of all the directed edges from 'from' to 'to'. If there are no such edges,
   *     return an empty list.
   */
  public List<Edge> getAllEdges(int from, int to) {
    if (!hasNode(from)) {
      throw new RuntimeException("No such node");
    }
    if (!hasNode(to)) {
      throw new RuntimeException("No such node");
    }
    LinkedList<Edge> edge_list = adjacencyList.get(from).get(to);
    if (edge_list == null || edge_list.isEmpty()) {
      return Collections.<Edge>emptyList();
    } else {
      return Collections.<Edge>unmodifiableList(edge_list);
    }
  }

  /**
   * A graph may have multiple directed edges between the same two nodes. This method returns a
   * (possibly empty) list of all edges from 'from' to 'to'.
   *
   * @param from The name of the source node.
   * @param to The name of the destination node.
   * @return A list of all the directed edges from 'from' to 'to'. If there are no such edges,
   *     return an empty list.
   */
  public Collection<Edge> getAllEdges(String from, String to) {
    return getAllEdges(getNode(from), getNode(to));
  }

  /**
   * Return an unmodifiable view of all the edges in this graph.
   *
   * @return An unmodifiable view of all the edges in this graph.
   */
  public Collection<Edge> getAllEdges() {
    return Collections.unmodifiableList(all_edges);
  }

  /**
   * Return an unmodifiable view of all the directed edges incident to the node 'from'.
   *
   * @param from The source node.
   * @return An unmodifiable view of all the directed edges incident to the node 'from'.
   */
  public Collection<Edge> getAllEdges(int from) {
    if (!hasNode(from)) {
      throw new RuntimeException("No such node");
    }
    Map<Integer, LinkedList<Edge>> edge_map = adjacencyList.get(from);
    // is this efficient for large graphs?
    return edge_map.values().stream().flatMap(Collection::stream).collect(Collectors.toList());
  }

  /**
   * Return an unmodifiable view of all the directed edges incident to the node 'from'.
   *
   * @param from The name of the source node.
   * @return An unmodifiable view of all the directed edges incident to the node 'from'.
   */
  public Collection<Edge> getAllEdges(String from) {
    return getAllEdges(getNode(from));
  }

  /**
   * Return an unmodifiable view of all the literals in this graph.
   *
   * @return An unmodifiable view of all the literals in this graph.
   */
  public Collection<Lit> getEdgeVars() {
    return Collections.unmodifiableList(all_edge_lits);
  }

  /**
   * If e is an edge f->t, and there exists an edge t->f in the graph, return the literal for such
   * an edge. (If there are multiple such edges, return the literal of an arbitrarily chosen edge
   * from among that set.) If there is no such edge, throw a RuntimeException.
   */
  public Lit getBackEdgeVar(Lit e) {
    return getBackEdge(e).l;
  }

  /**
   * If e is an edge f->t, and there exists an edge t->f in the graph, return such an edge. (If
   * there are multiple such edges, return the literal of an arbitrarily chosen edge from among that
   * set.) If there is no such edge, throw a RuntimeException
   */
  public Edge getBackEdge(Lit e) {
    Edge edge = getEdge(e);
    return getEdge(edge.to, edge.from);
  }

  /**
   * Check if this node is a node in the graph.
   *
   * @param node The node to test.
   */
  private void validateNode(int node) {
    if (node < 0 || node >= nNodes()) {
      throw new RuntimeException("Node " + node + " does not exist");
    }
  }
  /**
   * Get all outgoing and all incoming edges from this node.
   *
   * @param node The source node.
   * @return All outgoing and all incoming edges from this node.
   */
  public Collection<Lit> getEdgeVars(int node) {
    validateNode(node);
    return Collections.unmodifiableList(all_node_edge_lits.get(node));
  }

  /**
   * Get all outgoing edges from this node.
   *
   * @param node The source node.
   * @return All outgoing edges from this node.
   */
  public Collection<Lit> getOutgoingEdgeVars(int node) {
    validateNode(node);
    return Collections.unmodifiableList(all_out_edge_lits.get(node));
  }

  /**
   * Get all incoming edges to this node.
   *
   * @param node The source node.
   * @return All incoming edges to this node.
   */
  public Collection<Lit> getIncomingEdgeVars(int node) {
    validateNode(node);
    return Collections.unmodifiableList(all_in_edge_lits.get(node));
  }

  /**
   * If this edge has BitVector weighted edges, return the width of the bitvectors in this graph.
   * Else, return -1.
   *
   * @return The width of the bitvectors in this graph (or -1, if this graph has constant weights or
   *     is unweighted).
   */
  public int bitwidth() {
    return bitwidth;
  }

  /**
   * Add a new directed edge to the graph, from node 'from' to node 'to'. The edge will have weight
   * '1'.
   *
   * @param from The source node.
   * @param to The destination node.
   * @return The literal that controls whether this edge is included in the graph.
   */
  public Lit addEdge(int from, int to) {
    return addEdge(from, to, 1);
  }

  /**
   * Add a new directed edge to the graph, from node 'from' to node 'to'. The edge will have weight
   * '1'.
   *
   * @param from The source node.
   * @param to The destination node.
   * @param name An (optional) name for the edge literal. If empty, the literal will be unnamed.
   * @return The literal that controls whether this edge is included in the graph.
   */
  public Lit addEdge(int from, int to,String name) {
    return addEdge(from, to, 1,name);
  }
  /**
   * Add a new directed edge to the graph, from node 'from' to node 'to', with a constant weight.
   *
   * <p>Implementation note: Currently, reachability queries in MonoSAT are much faster if all of
   * the edges in a graph are strictly > 0.
   *
   * @param from The source node.
   * @param to The destination node.
   * @param constantWeight The constant weight of this edge. Must be >=0.
   * @return The literal that controls whether this edge is included in the graph.
   */
  public Lit addEdge(int from, int to, long constantWeight) {
    return addEdge(from,to,constantWeight,"");
  }
  /**
   * Add a new directed edge to the graph, from node 'from' to node 'to', with a constant weight.
   *
   * <p>Implementation note: Currently, reachability queries in MonoSAT are much faster if all of
   * the edges in a graph are strictly > 0.
   *
   * @param from The source node.
   * @param to The destination node.
   * @param constantWeight The constant weight of this edge. Must be >=0.
   * @param name An (optional) name for the edge literal. If empty, the literal will be unnamed.
   * @return The literal that controls whether this edge is included in the graph.
   */
  public Lit addEdge(int from, int to, long constantWeight, String name) {
    validateNode(from);
    validateNode(to);
    if (bitwidth >= 0) {
      return addEdge(from, to, solver.bv(bitwidth, constantWeight),name);
    } else {
      if(solver.hasLiteral(name)){
        throw new IllegalArgumentException("There is already a literal with name " + name + " in the solver");
      }
      Lit l =
          solver.toLit(MonosatJNI.newEdge(solver.getSolverPtr(), graphPtr, from, to, constantWeight));
      Map<Integer, LinkedList<Edge>> edge_map = adjacencyList.get(from);
      if (edge_map.get(to) == null) {
        edge_map.put(to, new LinkedList<>());
      }
      Edge e = new Edge(from, to, l, constantWeight);
      edge_map.get(to).add(e);
      edgeLitMap.put(l, e);
      all_edges.add(e);
      all_edge_lits.add(l);
      all_out_edge_lits.get(from).add(l);
      all_in_edge_lits.get(to).add(l);
      all_node_edge_lits.get(from).add(l);
      all_node_edge_lits.get(to).add(l);
      getSolver().addName(l,name);
      return l;
    }
  }

  /**
   * Add a new directed edge to the graph, from node 'from' to node 'to', with a BitVector edge
   * weight. The BitVector must have width = Graph.bitwidth(), and may be either a variable
   * BitVector, or a constant BitVector.
   *
   * @param from The source node.
   * @param to The destination node.
   * @param weight The BitVector weight of this edge.
   * @return The literal that controls whether this edge is included in the graph.
   */
  public Lit addEdge(int from, int to, BitVector weight) {
    return addEdge(from,to,weight,"");
  }
  /**
   * Add a new directed edge to the graph, from node 'from' to node 'to', with a BitVector edge
   * weight. The BitVector must have width = Graph.bitwidth(), and may be either a variable
   * BitVector, or a constant BitVector.
   *
   * @param from The source node.
   * @param to The destination node.
   * @param weight The BitVector weight of this edge.
   * @param name An (optional) name for the edge literal. If empty, the literal will be unnamed.
   * @return The literal that controls whether this edge is included in the graph.
   */
  public Lit addEdge(int from, int to, BitVector weight, String name) {
    validateNode(from);
    validateNode(to);
    if (this.bitwidth < 0) {
      throw new RuntimeException(
          "In order to use bitvector edge weights, the bitwidth must be passed to the graph constructor, eg:"
              + "Graph g = new Graph(solver, 8); will accept edges with bitvectors of size 8. Otherwise, edge weights are assumed to be constant integers.");
    }
    if(solver.hasLiteral(name)){
      throw new IllegalArgumentException("There is already a literal with name " + name + " in the solver");
    }
    Lit l = solver.toLit(MonosatJNI.newEdge_bv(solver.getSolverPtr(), graphPtr, from, to, weight.id));
    Map<Integer, LinkedList<Edge>> edge_map = adjacencyList.get(from);
    if (edge_map.get(to) == null) {
      edge_map.put(to, new LinkedList<>());
    }
    Edge e = new Edge(from, to, l, weight);
    edge_map.get(to).add(e);
    all_edges.add(e);
    edgeLitMap.put(l, e);
    all_edge_lits.add(l);
    all_out_edge_lits.get(from).add(l);
    all_in_edge_lits.get(to).add(l);
    all_node_edge_lits.get(from).add(l);
    all_node_edge_lits.get(to).add(l);
    getSolver().addName(l,name);
    return l;
  }

  /**
   * Add a new undirected edge to the graph, from node 'from' to node 'to'. The edge will have
   * weight '1'.
   *
   * @param from The source node.
   * @param to The destination node.
   * @return The literal that controls whether this edge is included in the graph.
   */
  public Lit addUndirectedEdge(int from, int to) {
    return addUndirectedEdge(from, to, 1);
  }

  /**
   * Add a new undirected edge to the graph, from node 'from' to node 'to'. The edge will have
   * weight '1'.
   *
   * @param from The source node.
   * @param to The destination node.
   * @param name An (optional) name for the edge literal. If empty, the literal will be unnamed.
   * @return The literal that controls whether this edge is included in the graph.
   */
  public Lit addUndirectedEdge(int from, int to, String name) {
    return addUndirectedEdge(from, to, 1,name);
  }

  /**
   * Add a new undirected edge to the graph, from node 'from' to node 'to', with a constant weight.
   *
   * <p>Implementation note: Currently, reachability queries in MonoSAT are much faster if all of
   * the edges in a graph are strictly > 0.
   *
   * @param from The source node.
   * @param to The destination node.
   * @param constantWeight The constant weight of this edge. Must be >=0.
   * @return The literal that controls whether this edge is included in the graph.
   */
  public Lit addUndirectedEdge(int from, int to, long constantWeight) {
    return addUndirectedEdge(from,to,constantWeight,"");
  }
  /**
   * Add a new undirected edge to the graph, from node 'from' to node 'to', with a constant weight.
   *
   * <p>Implementation note: Currently, reachability queries in MonoSAT are much faster if all of
   * the edges in a graph are strictly > 0.
   *
   * @param from The source node.
   * @param to The destination node.
   * @param constantWeight The constant weight of this edge. Must be >=0.
   * @param name An (optional) name for the edge literal. If empty, the literal will be unnamed.
   * @return The literal that controls whether this edge is included in the graph.
   */
  public Lit addUndirectedEdge(int from, int to, long constantWeight, String name) {
    validateNode(to);
    validateNode(from);
    if (bitwidth >= 0) {
      return addUndirectedEdge(from, to, solver.bv(bitwidth, constantWeight),name);
    } else {
      if(solver.hasLiteral(name)){
        throw new IllegalArgumentException("There is already a literal with name " + name + " in the solver");
      }
      Lit l = addEdge(from, to, constantWeight);
      Lit l2 = addEdge(to, from, constantWeight);
      solver.assertEqual(l, l2);
      getSolver().addName(l,name);
      getSolver().addName(l2,name+"_back");
      return l;
    }
  }

  /**
   * Add a new undirected edge to the graph, from node 'from' to node 'to', with a BitVector edge
   * weight. The BitVector must have width = Graph.bitwidth(), and may be either a variable
   * BitVector, or a constant BitVector.
   *
   * @param from The source node.
   * @param to The destination node.
   * @param weight The BitVector weight of this edge.
   * @return The literal that controls whether this edge is included in the graph.
   */
  public Lit addUndirectedEdge(int from, int to, BitVector weight) {
    return addUndirectedEdge(from,to,weight,"");
  }

  /**
   * Add a new undirected edge to the graph, from node 'from' to node 'to', with a BitVector edge
   * weight. The BitVector must have width = Graph.bitwidth(), and may be either a variable
   * BitVector, or a constant BitVector.
   *
   * @param from The source node.
   * @param to The destination node.
   * @param weight The BitVector weight of this edge.
   * @param name An (optional) name for the edge literal. If empty, the literal will be unnamed.
   * @return The literal that controls whether this edge is included in the graph.
   */
  public Lit addUndirectedEdge(int from, int to, BitVector weight, String name) {
    validateNode(to);
    validateNode(from);
    if (this.bitwidth < 0) {
      throw new RuntimeException(
          "In order to use bitvector edge weights, the bitwidth must be passed to the graph constructor, eg:"
              + "Graph g = new Graph(solver, 8); will accept edges with bitvectors of size 8. Otherwise, edge weights are assumed to be constant integers.");
    }
    if(solver.hasLiteral(name)){
      throw new IllegalArgumentException("There is already a literal with name " + name + " in the solver");
    }
    Lit l = addEdge(from, to, weight);
    Lit l2 = addEdge(to, from, weight);
    solver.assertEqual(l, l2);
    getSolver().addName(l,name);
    getSolver().addName(l2,name+"_back");
    return l;
  }

  /**
   * Returns a literal that evaluates to true if the graph is acyclic, in the assignment chosen by
   * the solver. If 'directed' is false, all edges are interpreted as undirected for the purposes of
   * cycle checking.
   *
   * @param directed If true, this predicate detects directed cycles in the graph. If false, then
   *     all edges in the graph are interpreted as undirected for the purposes of this predicate,
   *     and only undirected cycles are detected.
   * @return A literal that evaluates to true if the graph is acyclic, in the assignment chosen by
   *     the solver.
   */
  public Lit acyclic(boolean directed) {
    if (directed) {
      return solver.toLit(MonosatJNI.acyclic_directed(solver.getSolverPtr(), this.graphPtr));
    } else {
      return solver.toLit(MonosatJNI.acyclic_undirected(solver.getSolverPtr(), this.graphPtr));
    }
  }

  /**
   * Returns a literal that evaluates to true if the graph contains no directed cycles, in the
   * assignment chosen by the solver.
   *
   * @return A literal that evaluates to true if the graph is acyclic, in the assignment chosen by
   *     the solver.
   */
  public Lit acyclic() {
    return acyclic(true);
  }

  /**
   * Returns a literal that evaluates to true if there exists a (directed) path from node 'from' to
   * node 'to, in the assignment chosen by the solver.
   *
   * @param from The source node.
   * @param to The destination node.
   * @return A literal that evaluates to true if there exists a (directed) path from node 'from' to
   *     node 'to, in the assignment chosen by the solver.
   */
  public Lit reaches(int from, int to) {
    validateNode(from);
    validateNode(to);
    return solver.toLit(MonosatJNI.reaches(solver.getSolverPtr(), graphPtr, from, to));
  }
  /**
   * Returns a literal that evaluates to true if there exists a backwards path from -> to. For
   * example, if there is an edge u->v enabled in the graph, then there is a backwards path from v
   * to u.
   *
   * @param from The source node.
   * @param to The destination node.
   * @return A literal that evaluates to true if there exists a directed, backward path from node
   *     'from' to node 'to, in the assignment chosen by the solver.
   */
  public Lit reachesBackward(int from, int to) {
    validateNode(from);
    validateNode(to);
    return solver.toLit(MonosatJNI.reachesBackward(solver.getSolverPtr(), graphPtr, from, to));
  }
  /**
   * Returns a literal that evaluates to true if there exists a path from 'from' to nodeOnPath, and
   * a path from 'nodeOnPath' to 'to'. This constraint is logically equivalent to
   * And(g.reaches(from,nodeOnPath),g.reaches(nodeOnPath,to)), but may be more efficient if lots of
   * nodeOnPath checks are made between the same 'from' and 'to' nodes.
   *
   * @param nodeOnPath A node that should be on the path from 'from' to 'to'.
   * @param from The source node.
   * @param to The destination node.
   * @return A literal that evaluates to true if there exists a path from 'source' to 'destination'
   *     that passes through 'nodeOnPath'.
   */
  public Lit onPath(int nodeOnPath, int from, int to) {
    validateNode(from);
    validateNode(to);
    validateNode(nodeOnPath);
    return solver.toLit(MonosatJNI.onPath(solver.getSolverPtr(), graphPtr, nodeOnPath, from, to));
  }

  /**
   * Returns a literal that evaluates to true if the shortest path from 'from' to 'to' (in the
   * assignment chosen by the solver) has a (weighted) path length that satisfies 'comparison' with
   * respect to the constant 'compareTo'. *
   *
   * <p>Implementation note: Performing one sided comparisons (GT, GEQ, LT,LEQ) will be more
   * efficient than two sided comparisons (EQ,NEQ). One-sided comparisons are also more efficient
   * than directly instantiating a BitVector representation of the distance (@see
   * #Graph.distance(int,int)).
   *
   * @param from The source node.
   * @param to The destination node.
   * @param comparison The comparison to check (eg, GT, LT, EQ).
   * @param compareTo The constant to compare the shortest path length to.
   * @return A literal that evaluates to true if the shortest path from 'from' to 'to' (in the
   *     assignment chosen by the solver) has a (weighted) path length that satisfies 'comparison'
   *     with respect to the constant 'compareTo'.
   */
  public Lit compareDistance(int from, int to, Comparison comparison, long compareTo) {
    validateNode(from);
    validateNode(to);
    switch (comparison) {
      case GEQ:
        return solver
            .toLit(
                MonosatJNI.shortestPath_lt_const(solver.getSolverPtr(), graphPtr, from, to, compareTo))
            .not();
      case GT:
        return solver
            .toLit(
                MonosatJNI.shortestPath_leq_const(solver.getSolverPtr(), graphPtr, from, to, compareTo))
            .not();
      case LEQ:
        return solver.toLit(
            MonosatJNI.shortestPath_leq_const(solver.getSolverPtr(), graphPtr, from, to, compareTo));
      case LT:
        return solver.toLit(
            MonosatJNI.shortestPath_lt_const(solver.getSolverPtr(), graphPtr, from, to, compareTo));
      case EQ:
        {
          Lit l1 =
              solver.toLit(
                  MonosatJNI.shortestPath_leq_const(
                      solver.getSolverPtr(), graphPtr, from, to, compareTo));
          Lit l2 =
              solver.toLit(
                  MonosatJNI.shortestPath_lt_const(
                      solver.getSolverPtr(), graphPtr, from, to, compareTo));
          return solver.and(l1, l2);
        }
      case NEQ:
        {
          Lit l1 =
              solver.toLit(
                  MonosatJNI.shortestPath_leq_const(
                      solver.getSolverPtr(), graphPtr, from, to, compareTo));
          Lit l2 =
              solver.toLit(
                  MonosatJNI.shortestPath_lt_const(
                      solver.getSolverPtr(), graphPtr, from, to, compareTo));
          return solver.nand(l1, l2);
        }
    }
    throw new RuntimeException("Unknown comparison");
  }

  /**
   * Returns a literal that evaluates to true if the shortest path from 'from' to 'to' (in the
   * assignment chosen by the solver) has a (weighted) path length that satisfies 'comparison' with
   * respect to the BitVector 'compareTo'. *
   *
   * <p>Implementation note: Performing one sided comparisons (GT, GEQ, LT,LEQ) will be more
   * efficient than two sided comparisons (EQ,NEQ). One-sided comparisons are also more efficient
   * than directly instantiating a BitVector representation of the distance (@see
   * #Graph.distance(int,int)).
   *
   * @param from The source node.
   * @param to The destination node.
   * @param comparison The comparison to check (eg, GT, LT, EQ).
   * @param compareTo The constant to compare the shortest path length to.
   * @return A literal that evaluates to true if the shortest path from 'from' to 'to' (in the
   *     assignment chosen by the solver) has a (weighted) path length that satisfies 'comparison'
   *     with respect to the constant 'compareTo'.
   */
  public Lit compareDistance(int from, int to, Comparison comparison, BitVector compareTo) {
    validateNode(from);
    validateNode(to);
    switch (comparison) {
      case GEQ:
        return solver
            .toLit(
                MonosatJNI.shortestPath_lt_bv(solver.getSolverPtr(), graphPtr, from, to, compareTo.id))
            .not();
      case GT:
        return solver
            .toLit(
                MonosatJNI.shortestPath_leq_bv(solver.getSolverPtr(), graphPtr, from, to, compareTo.id))
            .not();
      case LEQ:
        return solver.toLit(
            MonosatJNI.shortestPath_leq_bv(solver.getSolverPtr(), graphPtr, from, to, compareTo.id));
      case LT:
        return solver.toLit(
            MonosatJNI.shortestPath_lt_bv(solver.getSolverPtr(), graphPtr, from, to, compareTo.id));
      case EQ:
        {
          Lit l1 =
              solver.toLit(
                  MonosatJNI.shortestPath_leq_bv(
                      solver.getSolverPtr(), graphPtr, from, to, compareTo.id));
          Lit l2 =
              solver.toLit(
                  MonosatJNI.shortestPath_lt_bv(
                      solver.getSolverPtr(), graphPtr, from, to, compareTo.id));
          return solver.and(l1, l2);
        }
      case NEQ:
        {
          Lit l1 =
              solver.toLit(
                  MonosatJNI.shortestPath_leq_bv(
                      solver.getSolverPtr(), graphPtr, from, to, compareTo.id));
          Lit l2 =
              solver.toLit(
                  MonosatJNI.shortestPath_lt_bv(
                      solver.getSolverPtr(), graphPtr, from, to, compareTo.id));
          return solver.nand(l1, l2);
        }
    }
    throw new RuntimeException("Unknown comparison");
  }

  /**
   * Returns a literal that evaluates to true if the shortest path from 'from' to 'to' (in the
   * assignment chosen by the solver) has a (weighted) path length that satisfies 'comparison' with
   * respect to the BitVector 'compareTo'. *
   *
   * <p>Implementation note: If your goal is to assert or reason about a one-sided comparison
   * between the shortest path length and a bitvector or constant, then this method is less
   * efficient than enforcing one-sided comparisons ((@see #Graph.distanceCompare(int,int,
   * Comparison,long))), as this method is equivalent to enforcing two one-sided comparisons between
   * this bitvector and the shortest path length in the graph.
   *
   * @param from The source node.
   * @param to The destination node.
   * @return A BitVector of width Graph.bitwidth(), which will be equal to the length of the
   *     shortest path from 'from' to 'to'.
   */
  public BitVector distance(int from, int to) {
    return distance(from, to, this.bitwidth);
  }

  /**
   * Returns a literal that evaluates to true if the shortest path from 'from' to 'to' (in the
   * assignment chosen by the solver) has a (weighted) path length that satisfies 'comparison' with
   * respect to the BitVector 'compareTo'. *
   *
   * <p>Implementation note: If your goal is to assert or reason about a one-sided comparison
   * between the shortest path length and a bitvector or constant, then this method is less
   * efficient than enforcing one-sided comparisons ((@see #Graph.distanceCompare(int,int,
   * Comparison,long))), as this method is equivalent to enforcing two one-sided comparisons between
   * this bitvector and the shortest path length in the graph.
   *
   * @param from The source node.
   * @param to The destination node.
   * @param bitwidth The bitwidth of the returned bitvector.
   * @return A BitVector of width bitwidth, which will be equal to the length of the shortest path
   *     from 'from' to 'to'.
   */
  public BitVector distance(int from, int to, int bitwidth) {
    validateNode(from);
    validateNode(to);
    if (bitwidth < 0) {
      // compute a 'large enough' bitwidth to hold the maximum possible distance
      bitwidth = 8; // temporary
    }

    BitVector result = new BitVector(solver, bitwidth);

    Lit l1 =
        solver.toLit(
            MonosatJNI.shortestPath_leq_bv(solver.getSolverPtr(), graphPtr, from, to, result.id));
    Lit l2 =
        solver.toLit(
            MonosatJNI.shortestPath_lt_bv(solver.getSolverPtr(), graphPtr, from, to, result.id));
    // result is geq to the shortest path, and is not greater than the shortest path, and so it is
    // exactly equal to the shortestPath
    solver.assertTrue(l1);
    solver.assertFalse(l2);
    return result;
  }

  /**
   * Returns a literal that compares the the maximum s-t- flow in this graph to a constant.
   *
   * <p>Implementation note: If your goal is to assert or reason about a one-sided comparison
   * between the maximum flow and a constant or bitvector, the compareMaximumFlow form is more
   * efficient then the direct maximumFlow() form.
   *
   * @param from The source node.
   * @param to The destination node.
   * @param comparison The comparison to check (eg, GT, LT, EQ).
   * @param compareTo The constant to compare the maximum flow to.
   * @return A literal that compares the the maximum s-t- flow in this graph to a constant.
   */
  public Lit compareMaximumFlow(int from, int to, Comparison comparison, long compareTo) {
    validateNode(from);
    validateNode(to);
    switch (comparison) {
      case GEQ:
        return solver.toLit(
            MonosatJNI.maximumFlow_geq(solver.getSolverPtr(), graphPtr, from, to, compareTo));
      case GT:
        return solver.toLit(
            MonosatJNI.maximumFlow_gt(solver.getSolverPtr(), graphPtr, from, to, compareTo));
      case LEQ:
        return solver
            .toLit(MonosatJNI.maximumFlow_gt(solver.getSolverPtr(), graphPtr, from, to, compareTo))
            .not();
      case LT:
        return solver
            .toLit(MonosatJNI.maximumFlow_geq(solver.getSolverPtr(), graphPtr, from, to, compareTo))
            .not();
      case EQ:
        {
          Lit l1 =
              solver.toLit(
                  MonosatJNI.maximumFlow_geq(solver.getSolverPtr(), graphPtr, from, to, compareTo));
          Lit l2 =
              solver.toLit(
                  MonosatJNI.maximumFlow_gt(solver.getSolverPtr(), graphPtr, from, to, compareTo));
          return solver.and(l1, l2);
        }
      case NEQ:
        {
          Lit l1 =
              solver.toLit(
                  MonosatJNI.maximumFlow_geq(solver.getSolverPtr(), graphPtr, from, to, compareTo));
          Lit l2 =
              solver.toLit(
                  MonosatJNI.maximumFlow_gt(solver.getSolverPtr(), graphPtr, from, to, compareTo));
          return solver.nand(l1, l2);
        }
    }
    throw new RuntimeException("Unknown comparison");
  }

  /**
   * Returns a literal that compares the the maximum s-t- flow in this graph to a BitVector.
   *
   * <p>Implementation note: If your goal is to assert or reason about a one-sided comparison
   * between the maximum flow and a constant or bitvector, the compareMaximumFlow form is more
   * efficient then the direct maximumFlow() form.
   *
   * @param from The source node.
   * @param to The destination node.
   * @param comparison The comparison to check (eg, GT, LT, EQ).
   * @param compareTo The BitVector to compare the maximum flow to.
   * @return A literal that compares the the maximum s-t- flow in this graph to a constant.
   */
  public Lit compareMaximumFlow(int from, int to, Comparison comparison, BitVector compareTo) {
    validateNode(from);
    validateNode(to);
    switch (comparison) {
      case GEQ:
        return solver.toLit(
            MonosatJNI.maximumFlow_geq_bv(solver.getSolverPtr(), graphPtr, from, to, compareTo.id));
      case GT:
        return solver.toLit(
            MonosatJNI.maximumFlow_gt_bv(solver.getSolverPtr(), graphPtr, from, to, compareTo.id));
      case LEQ:
        return solver
            .toLit(MonosatJNI.maximumFlow_gt_bv(solver.getSolverPtr(), graphPtr, from, to, compareTo.id))
            .not();
      case LT:
        return solver
            .toLit(
                MonosatJNI.maximumFlow_geq_bv(solver.getSolverPtr(), graphPtr, from, to, compareTo.id))
            .not();
      case EQ:
        {
          Lit l1 =
              solver.toLit(
                  MonosatJNI.maximumFlow_geq_bv(
                      solver.getSolverPtr(), graphPtr, from, to, compareTo.id));
          Lit l2 =
              solver.toLit(
                  MonosatJNI.maximumFlow_gt_bv(solver.getSolverPtr(), graphPtr, from, to, compareTo.id));
          return solver.and(l1, l2);
        }
      case NEQ:
        {
          Lit l1 =
              solver.toLit(
                  MonosatJNI.maximumFlow_geq_bv(
                      solver.getSolverPtr(), graphPtr, from, to, compareTo.id));
          Lit l2 =
              solver.toLit(
                  MonosatJNI.maximumFlow_gt_bv(solver.getSolverPtr(), graphPtr, from, to, compareTo.id));
          return solver.nand(l1, l2);
        }
    }
    throw new RuntimeException("Unknown comparison");
  }

  /**
   * Returns a BitVector that will be equal to the maximum s-t flow in this graph.
   *
   * <p>Implementation note: If your goal is to assert or reason about a one-sided comparison
   * between the maximum flow and a constant or bitvector, the compareMaximumFlow form is more
   * efficient then the direct maximumFlow() form.
   *
   * @param from The source node.
   * @param to The destination node.
   * @return A BitVector of width=Graph.bitwidth(), that will be equal to the maximum s-t flow in
   *     this graph.
   */
  public BitVector maximumFlow(int from, int to) {
    return maximumFlow(from, to, this.bitwidth);
  }

  /**
   * Returns a BitVector that will be equal to the maximum s-t flow in this graph, with the
   * specified bitwidth.
   *
   * <p>Implementation note: If your goal is to assert or reason about a one-sided comparison
   * between the maximum flow and a constant or bitvector, the compareMaximumFlow form is more
   * efficient then the direct maximumFlow() form.
   *
   * @param from The source node.
   * @param to The destination node.
   * @param bitwidth The width of the BitVector to create.
   * @return A BitVector of width=bitwidth, that will be equal to the maximum s-t flow in this
   *     graph.
   */
  public BitVector maximumFlow(int from, int to, int bitwidth) {
    validateNode(from);
    validateNode(to);
    if (bitwidth < 0) {
      // compute a 'large enough' bitwidth to hold the maximum possible flow
      bitwidth = 8; // temporary
    }
    BitVector result = new BitVector(solver, bitwidth);
    Lit l1 =
        solver.toLit(
            MonosatJNI.maximumFlow_geq_bv(solver.getSolverPtr(), graphPtr, from, to, result.id));
    Lit l2 =
        solver.toLit(MonosatJNI.maximumFlow_gt_bv(solver.getSolverPtr(), graphPtr, from, to, result.id));
    // result is geq to the max flow, and is not greater than the max flow (and so it is exactly
    // equal to the maxflow)
    solver.assertTrue(l1);
    solver.assertFalse(l2);
    return result;
  }

  /**
   * Retrieve a satisfying path (of nodes) from the SAT solver. This method may only be called after
   * a satisfying call to solver.solve().
   *
   * @param reach_or_distance_literal A reaches, distance, or onPath predicate belonging to this
   *     graph.
   * @return A new ArrayList, containing the nodes along a satisfying shortest path in the graph,
   *     corresponding to the given reaches, distance, or onPath literal.
   * @throws NoModelException If the solver does not have a satisfying assignment.
   */
  public ArrayList<Integer> getPathNodes(Lit reach_or_distance_literal) {
    getSolver().validate(reach_or_distance_literal);
    ArrayList<Integer> store = new ArrayList<>();
    if (!MonosatJNI.hasModel(solver.getSolverPtr())) {
      throw new NoModelException(
          "Solver has no model (this may indicate either that the solve() has not yet been called, or "
              + "that the most recent call to solve() returned a value other than true, or that a constraint was added into the solver after the last call to solve()).");
    }
    int length =
        MonosatJNI.getModel_Path_Nodes_Length(
            solver.getSolverPtr(), graphPtr, reach_or_distance_literal.l);
    if (length < 0) {
      throw new NoModelException("No path (this may indicate the solver has no model)");
    }
    IntBuffer buf = solver.getBuffer(0, length);
    int sz =
        MonosatJNI.getModel_Path_Nodes(
            solver.getSolverPtr(), graphPtr, reach_or_distance_literal.l, length, buf);
    assert (sz == length);
    for (int i = 0; i < length; i++) {
      int n = buf.get(i);
      store.add(n);
    }
    return store;
  }

  /**
   * Retrieve a satisfying path (of edge literals) from the SAT solver. This method may only be
   * called after a satisfying call to solver.solve().
   *
   * @param reach_or_distance_literal A reaches, distance, or onPath predicate belonging to this
   *     graph.
   * @return A new ArrayList, containing the literals corresponding to edges along a satisfying
   *     shortest path in the graph, corresponding to the given reaches, distance, or onPath
   *     literal.
   * @throws NoModelException If the solver does not have a satisfying assignment.
   */
  public ArrayList<Lit> getPathEdges(Lit reach_or_distance_literal) {
    getSolver().validate(reach_or_distance_literal);
    if (!MonosatJNI.hasModel(solver.getSolverPtr())) {
      throw new NoModelException(
          "Solver has no model (this may indicate either that the solve() has not yet been called, or that the most recent call to solve() returned a value other than true, or that a constraint was added into the solver after the last call to solve()).");
    }
    int length =
        MonosatJNI.getModel_Path_EdgeLits_Length(
            solver.getSolverPtr(), graphPtr, reach_or_distance_literal.l);
    if (length < 0) {
      throw new NoModelException("No path (this may indicate the solver has no model)");
    }
    ArrayList<Lit> store = new ArrayList<>();
    IntBuffer buf = solver.getBuffer(0, length);
    int sz =
        MonosatJNI.getModel_Path_EdgeLits(
            solver.getSolverPtr(), graphPtr, reach_or_distance_literal.l, length, buf);
    assert (sz == length);
    for (int i = 0; i < length; i++) {
      int l = buf.get(i);
      Lit lit = solver.toLit(l);
      store.add(lit);
    }
    return store;
  }

  /**
   * Retrieve the value of the maximum flow in this graph, from a satisfying solution in the solver,
   * corresponding to a maximum flow comparison literal.
   *
   * @param flowComparisonLiteral A maximum flow comparison literal.
   * @return The value of the maximum flow in this graph, from a satisfying solution in the solver.
   * @throws NoModelException If the solver does not have a satisfying assignment.
   */
  public long getMaxFlow(Lit flowComparisonLiteral) {
    return MonosatJNI.getModel_MaxFlow(solver.getSolverPtr(), graphPtr, flowComparisonLiteral.l);
  }

  /**
   * Retrieve the value of the flow along a given edge in the graph, corresponding to a specified
   * maximum flow comparison literal, from a satisfying solution in the solver.
   *
   * @param flowComparisonLiteral A maximum flow comparison literal.
   * @param edgeLiteral The edge literal to retrieve a flow for.
   * @return The value of the flow along this edge in this graph, from a satisfying solution in the
   *     solver.
   * @throws NoModelException If the solver does not have a satisfying assignment.
   */
  public long getEdgeFlow(Lit flowComparisonLiteral, Lit edgeLiteral) {
    return getEdgeFlow(flowComparisonLiteral, edgeLiteral, false);
  }

  /**
   * Retrieve the value of the flow along a given edge in the graph, corresponding to a specified
   * maximum flow comparison literal, from a satisfying solution in the solver.
   *
   * @param flowComparisonLiteral A maximum flow comparison literal.
   * @param edgeLiteral The edge literal to retrieve a flow for.
   * @param forceAcyclicFlow If true, then additional effort will be expended to ensure that the
   *     returned flow values are acyclic.
   * @return The value of the flow along this edge in this graph, from a satisfying solution in the
   *     solver.
   * @throws NoModelException If the solver does not have a satisfying assignment.
   */
  public long getEdgeFlow(Lit flowComparisonLiteral, Lit edgeLiteral, boolean forceAcyclicFlow) {
    if (forceAcyclicFlow) {
      return MonosatJNI.getModel_AcyclicEdgeFlow(
          solver.getSolverPtr(), graphPtr, flowComparisonLiteral.l, edgeLiteral.l);
    } else {
      return MonosatJNI.getModel_EdgeFlow(
          solver.getSolverPtr(), graphPtr, flowComparisonLiteral.l, edgeLiteral.l);
    }
  }

  /**
   * Create a string representation of this graph in Graphviz/.dot format.
   *
   * @param showModel If true, then only edges that are assigned non-false will be drawn. If false,
   *     all edges will be drawn.
   * @param showConstants If true, then edges will be labelled with weights.
   * @return A string representing this graph in Graphviz format.
   */
  public String draw(boolean showModel, boolean showConstants) {
    if (showModel && !solver.hasModel()) {
      showModel = false; // there is no model, show all edges
    }
    StringBuilder writer = new StringBuilder();
    writer.append("digraph{\n");
    for (int n = 0; n < nNodes(); n++) {
      writer.append("n" + Integer.toString(n) + " [label=\"" + getName(n) + "\"]\n");
    }

    for (Edge e : getAllEdges()) {
      writer.append(
          "n"
              + (e.from)
              + "->"
              + "n"
              +  (e.to)
              + "[label=\""
              + e.l.toString()
              + "\"");
      if (showModel) {
        Optional<Boolean> possibleValue = e.l.possibleValue();
        long weight = e.getWeight();
        if (possibleValue.isPresent() && !possibleValue.get()) {
          writer.append(",color=white");
        } else if (!possibleValue.isPresent() && showModel) {
          writer.append(",color=gray,style=dashed");
        } else {
          writer.append(",color=black");
          if (weight != 1) {
            // suppress weight 1, as that is the default
            writer.append(",label=\"" + Long.toString(weight) + "\"");
          }
        }
      } else {
        // do nothing
      }

      writer.append("]\n");
    }
    writer.append("}\n");
    return writer.toString();
  }

  /**
   * Create a string representation of this graph in Graphviz/.dot format.
   *
   * @return A string representing this graph in Graphviz format.
   */
  public String draw() {
    return draw(true, true);
  }

  /**
   * Create a string representation of this graph in Graphviz/.dot format.
   *
   * @param showModel If true, edges that are assigned false by the solver will be hidden. If false,
   *     all edges will be drawn.
   * @return A string representing this graph in Graphviz format.
   */
  public String draw(boolean showModel) {
    return draw(showModel, true);
  }

  /** Represents a directed edge in the graph. */
  public final class Edge {
    /** The source node of this directed edge. */
    public final int from;
    /** The destination node of this directed edge. */
    public final int to;
    /** The literal that controls this edge. */
    public final Lit l;
    /** The constant edge weight of this edge, or -1 if a BitVector edge weight is used. */
    private final long weight;
    /** The BitVector edge weight of this edge, or null if a constant edge weight is used. */
    private final BitVector bv;

    /**
     * For internal use only. To add an edge to the graph, @see #Graph.addEdge(int,int).
     *
     * @param from Source node.
     * @param to Destination node.
     * @param l Literal controlling this edge.
     * @param weight BitVector edge weight.
     */
    protected Edge(int from, int to, Lit l, BitVector weight) {
      this.from = from;
      this.to = to;
      this.l = l;
      this.bv = weight;
      this.weight = -1;
    }

    /**
     * For internal use only. To add an edge to the graph, @see #Graph.addEdge(int,int).
     *
     * @param from Source node.
     * @param to Destination node.
     * @param l Literal controlling this edge.
     * @param weight Constant edge weight.
     */
    protected Edge(int from, int to, Lit l, long weight) {
      this.from = from;
      this.to = to;
      this.l = l;
      this.bv = null;
      this.weight = weight;
    }

    /**
     * Check if this edge has a BitVector edge weight (as opposed to a constant edge weight).
     *
     * @return True if this edge has a BitVector edge weight (as opposed to a constant edge weight).
     */
    public boolean hasBV() {
      return bv != null;
    }

    /**
     * Return the BitVector weight of this edge, if this edge has a BitVector edge weight.
     *
     * @return The BitVector weight of this edge, if this edge has a BitVector edge weight.
     */
    public BitVector getBV() {
      if (!hasBV()) {
        throw new RuntimeException(
            "Attempt to access bitvector weight on edge with constant weight");
      }
      return bv;
    }

    /**
     * Get the weight of this edge. If the edge has a bitvector weight, return the value of that
     * bitvector.
     *
     * @return The weight of this edge.
     * @throws NoModelException If this edge has a BitVector weight, and the solver has no
     *     satisfying assignment for that Bitvector.
     */
    public long getWeight() {
      if (hasBV()) {
        return getBV().value();
        // throw new RuntimeException("Attempt to access constant weight on edge with bitvector
        // weight");
      }
      assert (weight >= 0);
      return weight;
    }
  }
}
