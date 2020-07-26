# MonoSAT File Format

MonoSAT is a SAT Modulo Theory solver for Boolean Monotonic Theories, featuring support for a wide set of graph properties, as well as a number of other theories, including finite state machine synthesis (still experimental).

This file documents MonoSAT's __*.GNF__ file format, which is a superset of the common [DIMACS CNF format](http://www.cs.ubc.ca/~hoos/SATLIB/Benchmarks/SAT/satformat.ps) (specifically, a superset of DIMACS as parsed by Minisat).

Just like a DIMACS CNF, a GNF instance begins with a `p cnf <#Vars> <#Clauses>` header line, followed by a number of clauses in 0-terminated format, one per line:

```
p cnf 4 3
1 3 -4 0
4 0
2 -3 0
```

In addition to clauses, a GNF instance may declare any number of graphs. A graph is declared with a line `digraph <weight type> <# nodes> <# edges> <GraphID>`, followed by at most "# edges" edge declarations. GraphID is a non-negative integer that must be unique for each graph; <weight type> may be one of 'int', 'float', 'rational', or may be ommited, in which case 'int' is assumed. After the graph is declared, edges may be declared for that graph on subsequent lines, as `edge <GraphID> <from> <to> <CNF Variable> [weight]`, where 'weight' is optional. If weight is specified, its type must match the graph's declared weight type. Rational values are read using the GMP [mpq_set_str] function. If weight is ommited, the edge has unit weight.

```
p cnf 4 3
1 3 -4 0
4 0
2 -3 0
digraph int 3 4 0
edge 0 0 1 1
edge 0 1 0 2
edge 0 1 2 3 
edge 0 0 2 4 4
```

This creates an integer-weighted graph (with ID '0'), with 3 nodes (numbered 0..2) and 4 directed edges, three of them unit weight, and the last with edge weight 4. Each edge is paired with a unique varaible from the CNF; any variables in the CNF may be used (that is, any variable <= #Vars in the CNF header). The first edge connects node 0 to node 1, IFF variable 1 is assigned true. The second edge connects node 1 to node 0, IFF varible 2 is assigned true. The third edge connects node 1 to node 2, IFF variable 3 is assigned true. The fourth edge connects node 0 to node 2, IFF variable 4 is assigned true. 

Two important restrictions to note are that 
  - edge variables are variables, *not* literals (that is, they must be positive integers)
  - every edge must be assigned a variable 
  - no two edges (or any other theory element from *any* graph in the instance) may share the same variable.

These restrictions reflect implementation choices in MonoSAT, and may be lifted in the future. 

You can specify constraints on the graph using any of the graph properties that MonoSAT supports.  
* Reachability Properties (*e.g.,* node `a` must reach node `b`, IFF `variable` is true): 
    *  `reach <GraphID> <a> <b> <variable>`
* Shortest Path Properties (*e.g.,* shortest path from node `a` to node `b` must have length <= `distance`, IFF `variable` is true):
    * Unweighted:  `distance_leq <GraphID> <a> <b> <variable> <distance>`, with `distance` a non-negative integer. All edges are treated as having unit weight, regardless of their specified weight.
    * Weighted:  `weighted_distance_leq <GraphID> <a> <b> <variable> <distance>`, where `distance `is parsed according to the graph's weight type.
* Maximum Flow Properties (*e.g.,* the maximum flow from node `s` to node `t` must be >= `flow`, IFF `variable` is true):
    *  `maximum_flow_geq <GraphID> <s> <t> <variable> <flow>`
* Minimum Weight Spanning Tree Properties (*e.g.,* the minimum weight spanning tree of graph `GraphID` must be <= `mstweight`, IFF `variable` is true)
    *  `mst_weight_leq <GraphID> <variable> <mstweight>`
* Acyclicity Properties (*e.g.,* `variable` is true IFF graph `GraphID` is a directed acyclic graph (or, for undirected acyclicity, is a forest)
    *  `acyclic <GraphID> <variable>`
    *  `forest <GraphID> <variable>`
    
Variants of these properties  with strict comparison operators are also supported, *e.g.,* `distance_lt`,  `maximum_flow_gt`, `mst_weight_lt`, .

Each of these graph properties is associated with a (unique) Boolean variable in the CNF, which must be true if and only if that graph property holds. This means that by asserting those variables to be true or false in the CNF, one can force any of these graph properties to be either true or false; but one can also trigger complex conditions in the CNF based on the truth-value of the graph property.

Graph reachability properties are specified as ```reach <GraphID> <from> <to> <variable>```. The variable (as with edge variables, a unique positive integer) must be assigned true if from reaches to, and must be false otherwise. For example, to specify that the there must exist a path from node 0 to node 2, we can create a reach property attached to variable 5, and assert that variable 5 must be true in the CNF:

```
p cnf 5 3
1 3 -4 0
4 0
2 -3 0
5 0
digraph int 3 4 0
edge 0 0 1 1
edge 0 1 0 2
edge 0 1 2 3 
edge 0 0 2 4 4
reach 0 0 2 5
```

By instead adding the unit clause '-5 0', you could have forced the graph to *not* have a path from node 0 to node 2. (This would be UNSAT in this example).

Finally, you can combine multiple graph properties, and combine them together with arbitrary Boolean logic. For example, you could assert that *either* node 2 must not be reachable from node 0, *or* the weighted shortest path from node 0 to node 2 must be <= 3:

```
p cnf 5 3
5 6 0
-5 -6 0
digraph int 3 4 0
edge 0 0 1 1
edge 0 1 0 2
edge 0 1 2 3 
edge 0 0 2 4 4
reach 0 0 2 5
weighted_distance_leq 0 0 2 6 3 
```

These are the graph properties that are currently well-supported by MonoSAT; many other useful graph properties are Boolean monotonic with respect to the edges in a graph, and could be supported in the future. Interesting possibilities include planarity detection, connected components, global minimum cuts, and many variatons of network flow properties. 

MonoSAT has experimental support for symbolic finite state machines;
see 'FiniteStateMachines.md' for an explanation of the file format.

[mpq_set_str]:https://gmplib.org/manual/Initializing-Rationals.html#Initializing-Rationals
