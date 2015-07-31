# MonoSAT
MonoSAT is a SAT Modulo Theory solver for *[monotonic theories]*, over Booleans and bitvectors. It supports a wide set of graph predicates (including reachability, shortest paths, maximum *s-t* flow, minimum spanning tree, and acyclicity constraints). MonoSAT also has limited support for geometric constraints involving convex hulls of point sets, and experimental support for constraints on finite state machines. 

MonoSAT now comes with a simplified, Z3-inspired Python 3 interface (see api/python). See installation instructions below; see also the [tutorial].

To see further examples of use cases for MonoSAT, and details on the (very simple) input file format that MonoSAT accepts, see  [FORMAT].

###Building
From the root directory, build MonoSAT with:

```
$cd Release
$make
```

or

```
$cd Static
$make
```

or (to build the shared library, required by the Python interface):

```
$cd SharedLibrary
$make
```

MonoSAT requires C++11 support, zlib, and GMP >= 5.1.3. Tested on Ubuntu 14.04, with G++ 4.8.2 and G++ 4.9. The python library requires Python 3.3+.

If you get errors along the lines of "error: could not convert ‘x’ from ‘__gmp_expr<__mpq_struct [1], __mpq_struct [1]>’ to ‘bool’", then you likely need to install a more recent version of GMP.

If you build MonoSAT without using the provided makefiles, it is critically important to compile with `NDEBUG` set (*i.e.,* `-DNDEBUG`), as otherwise many very expensive debugging assertions will be enabled. 

###Install the Python Library

To install the Python library (system-wide), first install the shared library, and then use Python's setuptools to install the Python library.

On Ubuntu (14.04):

$cd SharedLibrary
$make
$sudo cp libmonosat.so /usr/local/lib/
$cd ../api/python
$sudo python3 setup.py install


###Usage
MonoSAT is based on [MiniSat 2][Minisat], and supports many of the same calling conventions:

```
$monosat [-witness|-witness-file=filename] input_file.gnf
```

Where input_file.gnf is a file in [GNF format][FORMAT] (a very simple extension of DIMACS CNF format to support graph and geometry predicates). Use `-witness` to print the solution (if one exists) to stdout, or `-witness-file` to save it to file.

MonoSAT includes a very large set of configuration options - most of which you should stay away from unless you know what you are doing or want to explore the internals of MonoSAT (also, some of those configuration options might lead to buggy behaviour). Two options that often have a large impact on performance are `-decide-theories` and `-conflict-min-cut`:

```
$monosat -decide-theories -conflict-min-cut input_file.gnf
```

The `-decide-theories` option will cause the solver to make heuristic decisions that try to satisfy the various SMT predicates, which will often lead to improved performance, but can be pathologically bad in some common cases, and so is disabled by default. `-conflict-min-cut` will cause the solver to use a much slower, but more aggressive, clause learning strategy for reachability predicates; it may be useful for small, dificult instances.

###Source Overview
MonoSAT is written in C++. Core SAT solver functionality is in the `core/` and `simp/` directories; in particular, note `core/Config.cpp`, which is a central listing of all the configuration options available to MonoSAT. 

The graph and geometry theory solvers can be found in `geometry/` and `graph/`. Many of the graph algorithsms used by MonoSAT are collected in  `dgl/` (for 'Dynamic Graph Library'). 

`dgl/` incldudes C++ implementations of several dynamic graph algorithms (as well as some more common graph algorithms), and is well-optimized for medium sized (<20,000 nodes, < 100,000 edges), sparse graphs. The algorithms in dgl are designed for the case where the set of *possible* edges (and nodes) is fixed and known in advance (or only changes infrequently), and from that fixed set of possible edges many subsets of edges will subsequently be selected to be included in or excluded from the graph. 'dgl' supports templated edge weights and edge capacities, and has been tested successfully with integers, floats, and GMP arbitrary precision rationals.

######Algirthms implemented in 'dgl/' include:
* Reachability/Shortest Path
    * [Ramalingam-Reps](#ramalingam1996incremental) dynamic single-source shortest paths algorithm (with improvements described in [Buriol et al. 2008](#buriol2008speeding))
    * [Thorup](#thorup2000near)'s dynamic connectivity algorithm (for undirected graphs)
    * [Dijkstra](#dijkstra1959note)'s algorithm
    * DFS/BFS
* Minimum Spanning Tree
    * [Spira and Pan](#spira1975finding)'s dynamic minimum spanning tree algorithm
    * [Kruskal](#kruskal)'s algorithm
    * [Prim](#prims)'s algorithm
* Maximum *s-t* Flow
    * [Kohli and Torr](#kohlitorr)'s dynamic maximum flow algorithm
    * [Edmonds-Karp](#edmondskarp) algorithm (including a [dynamic variant](#dynamic_edmonds_karp))
    * [Dinitz](#dinitz)'s algorithm (including the [dynamic tree](#dynamic_tree) variant)
* And many supporting algorithms and data structures, including:
    * [Euler Tree](http://en.wikipedia.org/wiki/Euler_tour_technique)
    * [Link Cut Tree](http://en.wikipedia.org/?title=Link/cut_tree)
    * Splay Tree
    * Disjoint Set



###Licensing
The majority of MonoSAT is released under the [MIT license] (as documented in individual source files). However, some of the code, including some important libraries, fall under the GPL, and as a result, 
MonoSAT as a whole is currently released under the [GPL] (version 3 or later).

###Acknowledgements

MonoSAT was made possible by the use of several open-source projects, including the afore-mentioned [MiniSat], as well as a high-performance [dynamic maximum-flow algorithm] by Pushmeet Kohli and Philip Torr, Emil Stefanov's implementation of [Disjoint Sets], a [Link-Cut Tree] implementation by Daniel Sleator, and a [computational geometry library] by Chelton Evans.

[monotonic theories]: http://www.cs.ubc.ca/labs/isd/Projects/monosat/smmt.pdf
[FORMAT]: FORMAT.md
[tutorial]: TUTORIAL.md
[MiniSat]:http://minisat.se/

[MIT license]: http://opensource.org/licenses/MIT
[GPL]: http://www.gnu.org/licenses/gpl.html
[dynamic maximum-flow algorithm]:http://research.microsoft.com/en-us/um/people/pkohli/code/rrr.txt
[Link-Cut Tree]: http://codeforces.com/contest/117/submission/860934
[computational geometry library]:http://www.fluxionsdividebyzero.com/p1/math/geometry/geom.html
[Disjoint Sets]: http://web.rememberingemil.org/Projects/DisjointSets.aspx.html


###References


<a name="buriol2008speeding">
[Buriol, Luciana S., Mauricio GC Resende, and Mikkel Thorup. "Speeding up dynamic shortest-path algorithms." INFORMS Journal on Computing 20.2 (2008): 191-204.](http://dx.doi.org/10.1287/ijoc.1070.0231)</a>


<a name="dijkstra1959note">
[Dijkstra, Edsger W. "A note on two problems in connexion with graphs." Numerische mathematik 1.1 (1959): 269-271](http://dx.doi.org/10.1007%2FBF01386390)</a>

<a name="dinitz">
[Dinitz, Y. "Algorithm for solution of a problem of maximum flow in a network with power estimation". Doklady Akademii nauk SSSR 11: 1277–1280  (1970)](http://www.cs.bgu.ac.il/~dinitz/D70.pdf) </a>

<a name="edmondskarp">[Edmonds, Jack, and Richard M. Karp. "Theoretical improvements in algorithmic efficiency for network flow problems." Journal of the ACM (JACM) 19.2 (1972): 248-264](http://dx.doi.org/10.1145%2F321694.321699)</a>

<a name="kohli2005efficiently">
[Kohli, Pushmeet, and Philip HS Torr. "Efficiently solving dynamic markov random fields using graph cuts." Computer Vision, 2005. ICCV 2005. Tenth IEEE International Conference on. Vol. 2. IEEE, 2005](http://dx.doi.org/10.1109/ICCV.2005.81)</a>

<a name="dynamic_edmonds_karp">
[Korduban, D. "Incremental Maximum Flow in Dynamic graphs." Theoretical Computer Science Stack Exchange. http://cstheory.stackexchange.com/q/10186, 2012](http://cstheory.stackexchange.com/a/10186)</a>

<a name="kruskal">
[Kruskal, Joseph B. "On the shortest spanning subtree of a graph and the traveling salesman problem." Proceedings of the American Mathematical society 7.1 (1956): 48-50](http://dx.doi.org/10.1090%2FS0002-9939-1956-0078686-7)</a>

<a name="prim">
Prim, Robert Clay. "Shortest connection networks and some generalizations." Bell system technical journal 36.6 (1957): 1389-140](http://dx.doi.org/10.1002/j.1538-7305.1957.tb01515.x)</a>

<a name="ramalingam1996incremental"> [Ramalingam, Ganesan, and Thomas Reps. "An incremental algorithm for a generalization of the shortest-path problem." Journal of Algorithms 21.2 (1996): 267-305.](http://dx.doi.org/10.1006/jagm.1996.0046)</a>


<a name="dynamic_tree">
[Sleator, Daniel D., and Robert Endre Tarjan. "A data structure for dynamic trees." Proceedings of the thirteenth annual ACM symposium on Theory of computing. ACM, 1981.](http://dx.doi.org/10.1145/800076.802464)</a>


<a name="spira1975finding">
[Spira, Philip M., and A. Pan. "On finding and updating spanning trees and shortest paths." SIAM Journal on Computing 4.3 (1975): 375-380.](http://dx.doi.org/10.1137/0204032)</a>

<a name="thorup2000near">
[Thorup, Mikkel. "Near-optimal fully-dynamic graph connectivity." Proceedings of the thirty-second annual ACM symposium on Theory of computing. ACM, 2000.](http://dx.doi.org/10.1145/335305.335345)</a>
